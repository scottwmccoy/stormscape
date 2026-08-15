"""Near-real-time burn severity (dNBR) from CIMSS **BRISK**.

Where :mod:`stormscape.mrms` gives the rain that falls *on* a burn scar, this
module gives the scar. **BRISK** ("Burned area Rapid Interim Severity risK",
UW-Madison CIMSS/SSEC) maps every large (>~5,000 acre) US wildfire **daily**,
which is what makes it usable while a fire is still burning: the authoritative
products land far too late for that -- BAER soil burn severity arrives days to
weeks after containment and only for fires that get an assessment, MTBS a year
or more later. BRISK closes the gap that matters for post-fire debris-flow
hazard, where the rain can arrive before the fire is out.

It is a Google Earth Engine data-fusion dNBR composite over nine satellites
(GOES-East/West ABI, S-NPP / NOAA-20 / NOAA-21 VIIRS, Landsat 8/9, Sentinel-2a/b),
so a scar gets a coarse-but-immediate GOES look that sharpens as Landsat and
Sentinel-2 overpasses accumulate.

Access
------
The portal (``cimss.ssec.wisc.edu/brisk``) is a RealEarth viewer, and RealEarth's
WMTS tiles are **rendered PNG** -- pretty, but not data. The raw field is in the
open archive behind it, one GeoTIFF per fire per day::

    https://bin.ssec.wisc.edu/pub/realearth/brisk/<year>/<Fire>-<ST>-dNBR_<YYYYMMDD>_235959.tif

so this module reads the Apache directory index, parses those names into a
catalog, and downloads only the scenes that intersect the AOI. Finding which
fires those are is cheap because GDAL can range-read a GeoTIFF header over HTTP
(``/vsicurl/``): ~0.02 s per scene across a thread pool, so a whole day's ~60
fires are screened in under two seconds without downloading a pixel.

The companion BAER **soil** burn severity archive (``baer-data/``, ``product=
"sbs"``) is read the same way -- see the warning about it below.

Gotchas
-------
* **dNBR is a *vegetation* index, not soil burn severity.** The USGS post-fire
  debris-flow models are calibrated on **soil** burn severity (BAER SBS), which
  is dNBR *adjusted by field crews* for soil hydrophobicity, ground cover and
  duff consumption. BRISK is explicitly an **interim** product: use it to act
  early, then supersede it with SBS when the BAER assessment lands.
* **The scenes are EPSG:3857 at "60 m", which is not 60 m of ground.** Web
  Mercator metres shrink with latitude, so a 60 m cell is ~46.5 m at 39 deg N and
  ~42 m at 45 deg N. Fields are kept on their native 3857 grid (resampling the
  science data to make a rounder number would only lose fidelity), and every
  scene tested lands on an exact 60 m multiple, so mosaicking neighbouring fires
  is a paste rather than a warp.
* **NaN marks the area outside the burn, but the files do not tag a nodata
  value.** Reading with ``masked=True`` therefore masks nothing; test
  ``np.isfinite`` instead. Typically only ~10-15% of a scene's pixels are valid.
* **dNBR is unscaled here** (about -0.3 to 1.0), not the x1000 integer form the
  MTBS/USGS severity thresholds are usually quoted in. :data:`SEVERITY_SCHEMES`
  holds the same breaks divided by 1000.
* A fire's footprint **grows day to day**, so the same fire's scenes differ in
  extent; :func:`find_scenes` takes the latest per fire by default.

Outputs mirror :mod:`stormscape.mrms`: a result ``dict`` (``fields / transform /
crs / profile / meta``) ready for :func:`stormscape.mrms.save_fields` and
:func:`stormscape.plot.drape_i15`. No new dependencies.
"""

from __future__ import annotations

import datetime as dt
import json
import os
import re
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from typing import Optional, Sequence

import numpy as np
import pandas as pd

from .aoi import bbox_polygon, load_aoi

#: per-fire daily dNBR composites (CIMSS BRISK)
BRISK_BASE = "https://bin.ssec.wisc.edu/pub/realearth/brisk"
#: BAER soil burn severity, when an assessment exists (sparse)
BAER_BASE = "https://bin.ssec.wisc.edu/pub/realearth/baer-data"

PRODUCTS = {
    "dnbr": dict(base=BRISK_BASE, token="dNBR", kind="continuous",
                 label="dNBR", scale=1.0, first_year=2025),
    "sbs": dict(base=BAER_BASE, token="sbs", kind="class",
                label="soil burn severity class", scale=1.0, first_year=2025),
    # The BAER teams' own dNBR, published beside their soil-burn-severity
    # rasters. int16 **x1000** (the BARC convention -- NOAA's own
    # ``BARC256 = dNBR*5 - 275`` identity confirms the scaling), so it is
    # divided back to a plain index on read. 2025 only so far.
    "baer_dnbr": dict(base=BAER_BASE, token="dNBR", kind="continuous",
                      label="BAER dNBR", scale=1000.0, first_year=2025),
}

WORKERS = 12                 # parallel header reads / downloads
CACHE_DIR = "brisk_cache"    # sits at the event root, like nexrad_cache/
INDEX_TTL_H = 6.0            # re-list the *current* year this often
TIMEOUT = 60.0

#: dNBR severity breaks. The USGS/MTBS thresholds are published on the x1000
#: integer scale (100 / 270 / 440 / 660); BRISK's own portal colour scale steps
#: at 0.10 / 0.40 / 0.70. Class 0 is always "unburned or below the low break".
SEVERITY_SCHEMES = {
    "usgs": dict(
        breaks=(0.10, 0.27, 0.44, 0.66),
        labels=("unburned", "low", "moderate-low", "moderate-high", "high")),
    "brisk": dict(
        breaks=(0.10, 0.40, 0.70),
        labels=("unburned", "low", "moderate", "high")),
}

#: Days a BRISK composite needs before its *magnitude* is trustworthy.
#: Measured against the BAER teams' own dNBR on the 39 fires that carry both
#: (``examples/brisk_vs_baer.py``): BRISK's pattern is right immediately, but a
#: composite scored on the BAER assessment date reads **low** until it has
#: ingested a clear post-fire Landsat/Sentinel overpass. Every poor performer
#: recovered to r = 0.80-0.96 given a scene 5-21 days later, and a fixed +14 d
#: rule took fires at r >= 0.90 from 14/39 to 24/39. Advisory by default; a hard
#: filter with ``--min-age``.
MATURITY_DAYS = 14

#: BAER soil burn severity is delivered already classified.
SBS_LABELS = {1: "unburned/very low", 2: "low", 3: "moderate", 4: "high"}
#: Valid SBS class range. The rasters carry an embedded palette that colours 1-4
#: with BRISK's own four severity colours and paints **0 and 5+ the same black**,
#: i.e. the product itself treats anything outside 1-4 as not-a-severity (a
#: water/inholding/unmapped mask, ~3% of the one NV scene checked). So values
#: outside this range are read as missing rather than charted as a class.
SBS_VALID = (1, 4)

#: The **BAER burn-severity class colours**, verified from the products
#: themselves: all 77 of the 2025 soil-burn-severity rasters -- written by many
#: different BAER teams, in ERDAS Imagine -- carry this identical embedded
#: palette for classes 1-4 (only the class-5 mask colour varies). BRISK's own
#: ``qgis_BRISK_dNBR_colorscale_v2.txt`` uses exactly these four colours, so the
#: portal and the BAER deliverables already share one scheme; this is it.
BAER_CLASS_COLORS = ((0, 128, 128),      # 1 unburned / very low   teal
                     (82, 204, 204),     # 2 low                   cyan
                     (255, 232, 32),     # 3 moderate              yellow
                     (168, 0, 0))        # 4 high                  dark red

#: The same palette as a continuous dNBR ramp, as (dNBR, RGB) anchors -- the
#: interpolated form BRISK publishes, with its hard step at the 0.10 unburned
#: break. Registered as the ``"baer"`` colormap over [0, 1].
BAER_ANCHORS = ((0.00, (0, 128, 128)), (0.10, (0, 128, 128)),
                (0.10, (82, 204, 204)), (0.40, (255, 232, 32)),
                (0.70, (168, 0, 0)), (1.00, (114, 0, 0)))
BRISK_ANCHORS = BAER_ANCHORS             # back-compat alias; one scheme

_HREF = re.compile(r'href="([^"?/][^"]*\.tif)"')


# --------------------------------------------------------------------------- #
# catalog: parse the archive's directory index
# --------------------------------------------------------------------------- #
def _product(product: str) -> dict:
    try:
        return PRODUCTS[product]
    except KeyError:
        raise ValueError(f"product must be one of {sorted(PRODUCTS)}, "
                         f"got {product!r}") from None


def parse_name(filename: str, token: str = "dNBR") -> Optional[dict]:
    """Split ``Hidden-Valley-NV-dNBR_20260814_235959.tif`` into its parts.

    Returns ``{fire, state, date, filename}`` or ``None`` if the name does not
    match. The state is the trailing two-letter code when present -- a few
    scenes are named by incident number alone (``0231-OR``) or carry a
    non-state code, so it is informational, never a filter.
    """
    m = re.match(rf"^(?P<fire>.+?)-{re.escape(token)}_"
                 r"(?P<date>\d{8})_(?P<time>\d{6})\.tif$", filename)
    if not m:
        return None
    # BAER names carry a "-prelim-" infix (Alder-Springs-OR-prelim-dNBR_...).
    # Strip it before reading the state code, or the state is lost and the fire
    # name stops matching the same fire's BRISK scenes.
    fire = re.sub(r"-prelim$", "", m.group("fire"), flags=re.I)
    state = None
    sm = re.match(r"^(?P<name>.+)-(?P<st>[A-Za-z]{2})$", fire)
    if sm:
        fire, state = sm.group("name"), sm.group("st").upper()
    return dict(fire=fire, state=state,
                date=dt.datetime.strptime(m.group("date"), "%Y%m%d").date(),
                filename=filename)


def _list_url(url: str, timeout: float = TIMEOUT) -> list:
    """Filenames linked from an Apache directory index."""
    with urllib.request.urlopen(url, timeout=timeout) as resp:   # noqa: S310
        html = resp.read().decode("utf-8", "replace")
    return _HREF.findall(html)


def list_year(year: int, product: str = "dnbr", cache_dir: str = CACHE_DIR,
              ttl_h: float = INDEX_TTL_H, verbose: bool = True) -> pd.DataFrame:
    """Catalog one archive year as a DataFrame ``[fire, state, date, filename,
    url]``.

    Cached to ``<cache_dir>/index_<product>_<year>.csv``. Past years never
    change, so their cache never expires; the **current** year is re-listed
    once ``ttl_h`` hours have passed, which is what keeps "near real-time"
    honest without a request per call.
    """
    cfg = _product(product)
    os.makedirs(cache_dir, exist_ok=True)
    path = os.path.join(cache_dir, f"index_{product}_{year}.csv")
    stale = True
    if os.path.exists(path):
        age_h = (time.time() - os.path.getmtime(path)) / 3600.0
        stale = (year >= dt.date.today().year) and age_h > ttl_h
    if not stale and os.path.exists(path):
        df = pd.read_csv(path, parse_dates=["date"])
        df["date"] = df["date"].dt.date
        return df

    url = f"{cfg['base']}/{year}/"
    if verbose:
        print(f"listing {url}")
    try:
        names = _list_url(url)
    except Exception as exc:                                     # noqa: BLE001
        if os.path.exists(path):        # serve a stale cache rather than fail
            if verbose:
                print(f"note: {exc}; using cached index {path}")
            df = pd.read_csv(path, parse_dates=["date"])
            df["date"] = df["date"].dt.date
            return df
        raise RuntimeError(f"failed to list {url}: {exc}") from exc

    rows = []
    for n in names:
        rec = parse_name(n, cfg["token"])
        if rec:
            rec["url"] = f"{cfg['base']}/{year}/{n}"
            rows.append(rec)
    df = pd.DataFrame(rows, columns=["fire", "state", "date", "filename", "url"])
    df.to_csv(path, index=False)
    if verbose:
        print(f"  {len(df)} scenes, {df.fire.nunique()} fires")
    return df


def catalog(product: str = "dnbr", years: Optional[Sequence[int]] = None,
            cache_dir: str = CACHE_DIR, ttl_h: float = INDEX_TTL_H,
            verbose: bool = True) -> pd.DataFrame:
    """Catalog several archive years at once (default: this year and last)."""
    if years is None:
        y = dt.date.today().year
        years = [y - 1, y]
    frames = []
    for yr in years:
        try:
            frames.append(list_year(yr, product, cache_dir, ttl_h, verbose))
        except RuntimeError as exc:            # a year that does not exist yet
            if verbose:
                print(f"note: {exc}")
    if not frames:
        return pd.DataFrame(columns=["fire", "state", "date", "filename", "url"])
    return pd.concat(frames, ignore_index=True)


# --------------------------------------------------------------------------- #
# which scenes touch the AOI  (header range-reads, no pixels)
# --------------------------------------------------------------------------- #
def _bounds_cache_path(cache_dir: str, product: str) -> str:
    return os.path.join(cache_dir, f"bounds_{product}.json")


def _read_bounds(url: str):
    """``(west, south, east, north)`` in EPSG:4326 from a remote GeoTIFF header."""
    import rasterio
    from rasterio.warp import transform_bounds
    with rasterio.open("/vsicurl/" + url) as ds:
        return tuple(float(v) for v in
                     transform_bounds(ds.crs, "EPSG:4326", *ds.bounds))


def scene_bounds(urls: Sequence[str], cache_dir: str = CACHE_DIR,
                 product: str = "dnbr", workers: int = WORKERS,
                 verbose: bool = True) -> dict:
    """Map each scene URL to its lon/lat bounds, reading headers in parallel.

    Results are memoised in ``<cache_dir>/bounds_<product>.json``: a scene is
    immutable once published (its name carries the date), so a hit is always
    valid and the AOI screen gets free after the first run.
    """
    import rasterio
    os.makedirs(cache_dir, exist_ok=True)
    path = _bounds_cache_path(cache_dir, product)
    known = {}
    if os.path.exists(path):
        try:
            with open(path) as fh:
                known = json.load(fh)
        except (OSError, ValueError):
            known = {}
    todo = [u for u in urls if u not in known]
    if todo:
        if verbose:
            print(f"reading {len(todo)} scene headers "
                  f"({len(urls) - len(todo)} cached)")
        env = dict(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
                   CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",
                   GDAL_HTTP_MAX_RETRY="2", GDAL_HTTP_RETRY_DELAY="1")

        def one(u):
            try:
                return u, _read_bounds(u)
            except Exception:                                    # noqa: BLE001
                return u, None                    # unreadable scene -> skip it

        with rasterio.Env(**env), ThreadPoolExecutor(workers) as ex:
            for u, b in ex.map(one, todo):
                if b is not None:
                    known[u] = list(b)
        with open(path, "w") as fh:
            json.dump(known, fh)
    return {u: tuple(known[u]) for u in urls if u in known}


def _empty_scenes():
    """An empty scene table with the real columns.

    An AOI that has not burned is the ordinary answer, not an error, so every
    exit path returns the same shape -- callers can test ``len(...)`` and still
    read ``.fire`` / ``.geometry`` without a special case.
    """
    import geopandas as gpd
    return gpd.GeoDataFrame(
        {c: [] for c in ("fire", "state", "date", "url", "age_days",
                         "overlap")},
        geometry=[], crs=4326)


def find_scenes(aoi, product: str = "dnbr", date=None, since=None,
                latest: bool = True, years: Optional[Sequence[int]] = None,
                min_age_days: Optional[float] = None,
                cache_dir: str = CACHE_DIR, workers: int = WORKERS,
                verbose: bool = True):
    """Scenes whose footprint intersects the AOI, as a GeoDataFrame.

    Parameters
    ----------
    aoi
        Anything :func:`stormscape.aoi.load_aoi` accepts.
    date
        Consider only scenes published on or before this date (default: all).
        Pass a storm date to ask "what did the scar look like *then*" rather
        than today -- the archive keeps every day.
    since
        Drop scenes older than this date. Useful to exclude last season's fires.
    latest
        Keep only the most recent scene per fire (the default, and almost always
        what you want -- a fire has one file per day). ``False`` returns the
        full time series, e.g. to watch a scar grow.
    years
        Archive years to search (default: this year and last).

    Returns
    -------
    geopandas.GeoDataFrame
        Columns ``fire, state, date, url, geometry`` (footprint boxes,
        EPSG:4326), sorted largest-overlap first. Empty if nothing matches.
    """
    import geopandas as gpd

    bounds, geom = load_aoi(aoi)
    target = geom if geom is not None else bbox_polygon(bounds)

    df = catalog(product, years=years, cache_dir=cache_dir, verbose=verbose)
    if not len(df):
        return _empty_scenes()
    # A scene's "age" is days since the fire first appeared in the archive --
    # a proxy for how much post-fire imagery the composite has had a chance to
    # absorb. Computed on the FULL catalog, before any date/since filtering,
    # or trimming the early scenes would make an old fire look brand new.
    # (A fire that started before the earliest searched year reads younger than
    # it is; widen `years` if that matters.)
    first_seen = df.groupby("fire")["date"].min()
    if date is not None:
        d = date if isinstance(date, dt.date) else _as_date(date)
        df = df[df.date <= d]
    if since is not None:
        s = since if isinstance(since, dt.date) else _as_date(since)
        df = df[df.date >= s]
    if latest and len(df):
        df = df.sort_values("date").groupby("fire", as_index=False).last()
    if not len(df):
        return _empty_scenes()

    df = df.assign(age_days=[(d - first_seen[f]).days
                             for f, d in zip(df.fire, df.date)])

    bmap = scene_bounds(list(df.url), cache_dir, product, workers, verbose)
    rows = []
    for r in df.itertuples(index=False):
        b = bmap.get(r.url)
        if b is None:
            continue
        poly = bbox_polygon(b)
        if not poly.intersects(target):
            continue
        rows.append(dict(fire=r.fire, state=r.state, date=r.date, url=r.url,
                         age_days=int(r.age_days),
                         overlap=poly.intersection(target).area,
                         geometry=poly))
    if not rows:                 # scenes existed, none of them reach this AOI
        return _empty_scenes()
    out = gpd.GeoDataFrame(rows, geometry="geometry", crs=4326)
    out = out.sort_values("overlap", ascending=False).reset_index(drop=True)

    # Maturity is screened *after* the AOI intersection, so the message names
    # only fires the caller actually asked about -- filtering the whole catalog
    # first would report hundreds of irrelevant fires nationwide.
    if min_age_days is not None:
        young = out[out.age_days < min_age_days]
        out = out[out.age_days >= min_age_days].reset_index(drop=True)
        if len(young) and verbose:
            # Dropping a fire from a hazard map is loud on purpose: the caller
            # asked for mature composites, and quietly substituting an immature
            # one -- which under-reads severity -- would defeat the request.
            named = ", ".join(f"{r.fire} ({r.age_days} d)"
                              for r in young.head(8).itertuples(index=False))
            more = f", +{len(young) - 8} more" if len(young) > 8 else ""
            print(f"--min-age {min_age_days:g} d drops {len(young)} fire(s) in "
                  f"this AOI whose composite is too young: {named}{more}")
        if not len(out):
            return _empty_scenes()
    if verbose:
        print(f"{len(out)} scene(s) intersect the AOI"
              + (f": {', '.join(out.fire.astype(str))}" if len(out) else ""))
    return out


def _as_date(value) -> dt.date:
    """Accept ``date``/``datetime``/``'YYYYMMDD'``/``'YYYY-MM-DD'``."""
    if isinstance(value, dt.datetime):
        return value.date()
    if isinstance(value, dt.date):
        return value
    s = str(value).strip().replace("-", "")
    return dt.datetime.strptime(s, "%Y%m%d").date()


# --------------------------------------------------------------------------- #
# fetch + cache the scenes, mosaic onto one grid
# --------------------------------------------------------------------------- #
def fetch_scene(url: str, cache_dir: str = CACHE_DIR, retries: int = 2,
                verbose: bool = True) -> str:
    """Download one scene into ``cache_dir`` (skipped if already there).

    Scenes are immutable -- the filename carries the date -- so a cached file is
    never stale and re-running a map costs no network at all.
    """
    os.makedirs(cache_dir, exist_ok=True)
    path = os.path.join(cache_dir, os.path.basename(url))
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path
    if verbose:
        print(f"downloading {os.path.basename(url)}")
    tmp = path + ".part"
    for attempt in range(retries + 1):
        try:
            urllib.request.urlretrieve(url, tmp)                  # noqa: S310
            os.replace(tmp, path)
            return path
        except Exception as exc:                                 # noqa: BLE001
            if os.path.exists(tmp):
                os.remove(tmp)
            if attempt >= retries:
                raise RuntimeError(f"failed to download {url}: {exc}") from exc
            time.sleep(2 * (attempt + 1))
    return path


def _target_grid(paths, bounds, dst_crs, resolution=None):
    """Output grid: the AOI in ``dst_crs``, snapped outward to the source cell
    size so a same-CRS scene lands on it exactly (a paste, not a resample)."""
    import rasterio
    from rasterio.warp import transform_bounds
    with rasterio.open(paths[0]) as ds:
        src_crs = ds.crs
        res = resolution or abs(ds.res[0])
    dst_crs = dst_crs or src_crs
    w, s, e, n = transform_bounds("EPSG:4326", dst_crs, *bounds)
    left = np.floor(w / res) * res
    bottom = np.floor(s / res) * res
    right = np.ceil(e / res) * res
    top = np.ceil(n / res) * res
    width = max(1, int(round((right - left) / res)))
    height = max(1, int(round((top - bottom) / res)))
    transform = rasterio.transform.from_origin(left, top, res, res)
    return dst_crs, transform, width, height, res


def mosaic(paths: Sequence[str], bounds, dst_crs=None, resolution=None,
           categorical: bool = False):
    """Merge scenes onto one grid clipped to ``bounds`` (lon/lat).

    Each scene is reprojected into a shared destination grid and combined with a
    NaN-aware **maximum**, so where two fires overlap the more severe value
    wins and a scene's NaN surround never erases a neighbour's data. For BRISK
    the reprojection is an identity copy (every scene is already EPSG:3857 on
    the same 60 m grid); the general path is what lets the UTM-tiled BAER
    scenes mosaic too.
    """
    import rasterio
    from rasterio.enums import Resampling
    from rasterio.warp import reproject

    dst_crs, transform, width, height, res = _target_grid(
        paths, bounds, dst_crs, resolution)
    out = np.full((height, width), np.nan, dtype="float32")
    how = Resampling.nearest if categorical else Resampling.bilinear
    for p in paths:
        with rasterio.open(p) as ds:
            src = ds.read(1).astype("float32")
            if ds.nodata is not None and np.isfinite(ds.nodata):
                src = np.where(src == ds.nodata, np.nan, src)
            if categorical:      # see SBS_VALID: 0 and 5+ are mask, not severity
                lo, hi = SBS_VALID
                src = np.where((src < lo) | (src > hi), np.nan, src)
            tmp = np.full((height, width), np.nan, dtype="float32")
            reproject(source=src, destination=tmp,
                      src_transform=ds.transform, src_crs=ds.crs,
                      dst_transform=transform, dst_crs=dst_crs,
                      src_nodata=np.nan, dst_nodata=np.nan, resampling=how)
        out = np.fmax(out, tmp)          # fmax ignores NaN on either side
    return out, transform, dst_crs


def classify(dnbr, scheme: str = "usgs"):
    """Severity class index for a dNBR array (0 = unburned, rising with break).

    NaN in, NaN out -- the unburned *surround* stays missing rather than being
    called class 0, which is a real distinction: "we looked and it did not burn"
    versus "we have no observation here".
    """
    try:
        breaks = SEVERITY_SCHEMES[scheme]["breaks"]
    except KeyError:
        raise ValueError(f"scheme must be one of {sorted(SEVERITY_SCHEMES)}, "
                         f"got {scheme!r}") from None
    a = np.asarray(dnbr, dtype="float32")
    cls = np.zeros(a.shape, dtype="float32")
    for b in breaks:
        cls += (a >= b).astype("float32")
    return np.where(np.isfinite(a), cls, np.nan).astype("float32")


def burn_severity(aoi, date=None, product: str = "dnbr", scheme: str = "usgs",
                  since=None, fires: Optional[Sequence[str]] = None,
                  years: Optional[Sequence[int]] = None,
                  min_age_days: Optional[float] = None,
                  cache_dir: str = CACHE_DIR, pad_deg: float = 0.02,
                  workers: int = WORKERS, verbose: bool = True):
    """Burn severity over an AOI as an ``mrms``-style result dict.

    Finds every fire intersecting the AOI, downloads and caches those scenes,
    mosaics them, and returns ``fields`` with:

    ``dnbr``
        the composite index on its native EPSG:3857 grid (``sbs`` returns the
        delivered classes here instead);
    ``severity``
        the class index from :func:`classify` (dNBR only).

    Returns ``None`` when no fire intersects the AOI -- an ordinary answer, not
    an error, and the common one for an AOI that has not burned.
    """
    cfg = _product(product)
    scenes = find_scenes(aoi, product=product, date=date, since=since,
                         years=years, min_age_days=min_age_days,
                         cache_dir=cache_dir, workers=workers, verbose=verbose)
    if fires:
        want = {f.lower() for f in fires}
        scenes = scenes[scenes.fire.str.lower().isin(want)]
    if not len(scenes):
        if verbose:
            print("no burn-severity scenes intersect this AOI")
        return None

    paths = [fetch_scene(u, cache_dir, verbose=verbose) for u in scenes.url]
    bounds, _ = load_aoi(aoi, pad_deg=pad_deg)
    categorical = cfg["kind"] == "class"
    arr, transform, crs = mosaic(paths, bounds, dst_crs="EPSG:3857",
                                 categorical=categorical)
    scale = float(cfg.get("scale", 1.0))
    if scale != 1.0:                     # BAER dNBR ships x1000; see PRODUCTS
        arr = arr / scale

    # Advisory, not a filter: an immature composite reads LOW, which for
    # post-fire hazard work is the dangerous direction to be wrong in.
    young = scenes[scenes.age_days < MATURITY_DAYS] if "age_days" in scenes else []
    if verbose and len(young) and product == "dnbr":
        print(f"note: {len(young)} composite(s) younger than {MATURITY_DAYS} d "
              f"({', '.join(f'{r.fire} {r.age_days} d' for r in young.itertuples(index=False))})"
              f"; pattern is reliable but the magnitude may under-read "
              f"-- see --min-age")

    fields = {"dnbr": arr}
    if not categorical:
        fields["severity"] = classify(arr, scheme)
    else:
        fields = {"severity": arr}

    h, w = arr.shape
    profile = dict(driver="GTiff", height=h, width=w, count=1, dtype="float32",
                   crs=crs, transform=transform, nodata=np.nan, compress="LZW")
    finite = arr[np.isfinite(arr)]
    meta = dict(source="CIMSS BRISK" if product == "dnbr" else "BAER SBS",
                product=product, scheme=scheme if not categorical else "baer",
                fires=list(scenes.fire), states=list(scenes.state),
                scene_dates=[str(d) for d in scenes.date],
                scenes=[os.path.basename(p) for p in paths],
                burned_px=int(finite.size),
                dnbr_max=float(finite.max()) if finite.size else float("nan"),
                dnbr_p98=float(np.percentile(finite, 98))
                if finite.size else float("nan"))
    if not categorical and finite.size:
        breaks = SEVERITY_SCHEMES[scheme]["breaks"]
        labels = SEVERITY_SCHEMES[scheme]["labels"]
        cls = classify(arr, scheme)
        cf = cls[np.isfinite(cls)]
        meta["class_fraction"] = {labels[i]: float((cf == i).mean())
                                  for i in range(len(breaks) + 1)}
    return dict(fields=fields, transform=transform, crs=crs, profile=profile,
                meta=meta)


# --------------------------------------------------------------------------- #
# the portal's own colour ramp
# --------------------------------------------------------------------------- #
def severity_colors(scheme: str = "usgs"):
    """``(ListedColormap, BoundaryNorm, tick_positions, labels)`` for a **classed**
    dNBR map drawn in the BAER colours.

    This is what makes a stormscape burn map read like a BAER deliverable: the
    field is banded at the scheme's severity breaks rather than shaded
    continuously, and the colour bar is labelled with the class names instead of
    dNBR numbers.

    A scheme with three breaks (four classes) maps exactly onto the four
    official :data:`BAER_CLASS_COLORS`. The five-class ``usgs`` scheme has no
    official fifth colour, so its colours are **sampled from the same ramp** at
    class midpoints -- still the BAER scheme, but interpolated rather than the
    literal published table.
    """
    import numpy as _np
    from matplotlib.colors import BoundaryNorm, ListedColormap

    spec = SEVERITY_SCHEMES.get(scheme)
    if spec is None:
        raise ValueError(f"scheme must be one of {sorted(SEVERITY_SCHEMES)}, "
                         f"got {scheme!r}")
    breaks, labels = list(spec["breaks"]), list(spec["labels"])
    n = len(labels)
    if n == len(BAER_CLASS_COLORS):
        colors = [tuple(c / 255.0 for c in rgb) for rgb in BAER_CLASS_COLORS]
    else:
        ramp = register_baer_cmap()
        edges = [0.0] + breaks + [1.0]
        colors = [ramp((edges[i] + edges[i + 1]) / 2.0)[:3] for i in range(n)]
    # bound the outer classes for display: dNBR runs about -0.3 to 1.5
    edges = [-1.0] + breaks + [2.0]
    norm = BoundaryNorm(edges, ncolors=n)
    ticks = [(edges[i] + edges[i + 1]) / 2.0 for i in range(n)]
    # the open-ended end classes would otherwise centre their tick off-scale
    ticks[0] = (breaks[0] + max(-0.3, edges[0])) / 2.0
    ticks[-1] = (breaks[-1] + 1.0) / 2.0
    return ListedColormap(colors), norm, ticks, labels


def register_baer_cmap(name: str = "baer"):
    """Register the BAER / BRISK dNBR ramp as a matplotlib colormap.

    Anchored over dNBR **0 to 1**, the domain
    :func:`stormscape.plot.drape_i15` draws in (``vmin=0``), so a continuous map
    made with ``cmap="baer", vmax=1.0`` matches both the CIMSS portal and the
    BAER products' own colours. For the *classed* rendering BAER actually
    publishes, use :func:`severity_colors`.
    """
    import matplotlib
    from matplotlib.colors import LinearSegmentedColormap
    try:
        return matplotlib.colormaps[name]
    except KeyError:
        pass
    pts, seen = [], set()
    for v, rgb in BAER_ANCHORS:
        # the ramp steps discontinuously at 0.10; nudge the duplicate stop so
        # from_list keeps the hard edge instead of rejecting the repeat
        while v in seen:
            v += 1e-6
        seen.add(v)
        pts.append((v, tuple(c / 255.0 for c in rgb)))
    # 1024 levels rather than the default 256: the ramp steps discontinuously at
    # the unburned break, and a coarse lookup table smears that edge across a
    # visible band of dNBR instead of holding it at 0.10
    cmap = LinearSegmentedColormap.from_list(name, pts, N=1024)
    matplotlib.colormaps.register(cmap, name=name)
    return cmap


def register_brisk_cmap(name: str = "brisk"):
    """Back-compat alias for :func:`register_baer_cmap`.

    BRISK's published ramp and the BAER products' embedded palette are the same
    four colours, so ``"brisk"`` and ``"baer"`` name one scheme.
    """
    return register_baer_cmap(name)
