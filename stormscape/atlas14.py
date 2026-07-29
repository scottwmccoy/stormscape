"""NOAA **Atlas 14** rainfall climatology as raster fields for an AOI.

Where :mod:`stormscape.mrms` / :mod:`stormscape.nexrad` give the *observed* peak
15/30/60-min rainfall-intensity fields of a storm, this module gives the
**climatological** intensity for the same durations at a chosen recurrence
interval (default the **1-year** ARI, the reference established by the group's
prior work) -- so a storm can be put in context (observed vs climatology) and an
**anomaly** (observed / climatology) computed cell-wise.

Access
------
The ``pfdf.data.noaa.atlas14`` helper is **point-only** (one CSV per lat/lon), so
it cannot return fields. Instead this module fetches NOAA's authoritative
**gridded** Atlas 14 product -- the zipped ESRI-ASCII grids served from the PFDS
GIS tree::

    https://hdsc.nws.noaa.gov/pub/hdsc/data/<region>/<region><ari>yr<dur><stat>.zip

* ``region`` -- an Atlas 14 volume/region code (``sw`` = semiarid Southwest +
  California, covers NV/CA/AZ/UT/NM; full table in ``data/atlas14_regions.csv``);
* ``ari``    -- recurrence interval in years (``1``, ``2``, ``5``, ...);
* ``dur``    -- duration code (``15m`` ``30m`` ``60m`` ``06h`` ``24h`` ...);
* ``stat``   -- ``a`` mean (default), ``al`` / ``au`` lower / upper 90% CI;
  an ``_ams`` suffix selects the annual-maximum series (default is partial
  duration).

Each grid is geographic **EPSG:4269** (NAD83), ~800 m cells, NODATA ``-999``,
**units = 1000ths of an inch (depth)**. We convert depth to an intensity that
matches the observed mm/h fields::

    i_mmph = mils / 1000 * 25.4 * (60 / duration_min)

Outputs mirror :mod:`stormscape.mrms`: a result ``dict``
(``fields / transform / crs / profile / meta``) ready for
:func:`stormscape.mrms.save_fields` and :func:`stormscape.plot.drape_i15`.
No extra dependencies -- rasterio's AAIGrid driver + rioxarray do the work.
"""

from __future__ import annotations

import os
import urllib.request
import zipfile
from importlib import resources
from typing import Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
import rioxarray  # noqa: F401  (registers the .rio accessor)
import xarray as xr

from .aoi import load_aoi

BASE_URL = "https://hdsc.nws.noaa.gov/pub/hdsc/data/"

# duration in minutes -> NOAA grid duration code
DURATIONS = {
    5: "05m", 10: "10m", 15: "15m", 30: "30m", 60: "60m",
    120: "02h", 180: "03h", 360: "06h", 720: "12h", 1440: "24h",
    2880: "48h", 4320: "03d",
}
# statistic name -> NOAA grid suffix (mean estimate / 90% CI bounds)
STAT = {"mean": "a", "lower": "al", "upper": "au"}

MM_PER_INCH = 25.4

_REGIONS: Optional[pd.DataFrame] = None


# --------------------------------------------------------------------------- #
# region table  (data/atlas14_regions.csv, built from the grid headers)
# --------------------------------------------------------------------------- #
def _regions() -> pd.DataFrame:
    """Bundled Atlas 14 region table (region, description, W, S, E, N)."""
    global _REGIONS
    if _REGIONS is None:
        with resources.files("stormscape").joinpath(
                "data/atlas14_regions.csv").open() as fh:
            _REGIONS = pd.read_csv(fh)
    return _REGIONS


def region_for_bounds(bounds: Sequence[float]) -> str:
    """Pick the Atlas 14 region code covering an AOI ``(W, S, E, N)``.

    Prefers a region whose extent fully contains the AOI; if several do (the
    volumes overlap at their seams) the one whose centre is nearest the AOI
    centre wins. Falls back to regions merely containing the AOI centre, and
    raises with the available codes if none match (pass ``region=`` then).
    """
    w, s, e, n = bounds
    clon, clat = (w + e) / 2.0, (s + n) / 2.0
    r = _regions()
    contains = r[(r.west <= w) & (r.east >= e) & (r.south <= s) & (r.north >= n)]
    pool = contains if len(contains) else r[
        (r.west <= clon) & (r.east >= clon)
        & (r.south <= clat) & (r.north >= clat)]
    if not len(pool):
        raise ValueError(
            f"no NOAA Atlas 14 region covers bounds {tuple(bounds)}; pass "
            f"region= explicitly. available: {sorted(r.region)}")
    cx = (pool.west + pool.east) / 2.0
    cy = (pool.south + pool.north) / 2.0
    dist = (cx - clon) ** 2 + (cy - clat) ** 2
    return str(pool.loc[dist.idxmin(), "region"])


# --------------------------------------------------------------------------- #
# fetch + cache the gridded ASCII
# --------------------------------------------------------------------------- #
def _dur_code(duration: int) -> str:
    if duration not in DURATIONS:
        raise ValueError(
            f"unsupported duration {duration} min; choose from "
            f"{sorted(DURATIONS)}")
    return DURATIONS[duration]


def grid_url(region: str, ari: int, duration: int,
             stat: str = "a", series: str = "pds") -> str:
    """URL of the zipped ASCII grid for one region/ARI/duration/statistic."""
    suffix = "" if series == "pds" else "_ams"
    return (f"{BASE_URL}{region}/"
            f"{region}{ari}yr{_dur_code(duration)}{stat}{suffix}.zip")


def fetch_grid(region: str, ari: int, duration: int, stat: str = "a",
               series: str = "pds", cache_dir: str = "atlas14_cache",
               retries: int = 2, verbose: bool = True) -> str:
    """Download (if needed), unzip, and return the path to the ``.asc`` grid.

    The zip ships the ``.asc`` plus a ``.prj`` (EPSG:4269) and ``.xml``; all are
    extracted next to each other so rasterio reads the CRS from the sidecar.
    """
    suffix = "" if series == "pds" else "_ams"
    stem = f"{region}{ari}yr{_dur_code(duration)}{stat}{suffix}"
    os.makedirs(cache_dir, exist_ok=True)
    zpath = os.path.join(cache_dir, stem + ".zip")
    apath = os.path.join(cache_dir, stem + ".asc")
    if not os.path.exists(apath):
        if not os.path.exists(zpath):
            url = grid_url(region, ari, duration, stat=stat, series=series)
            if verbose:
                print(f"downloading {url}")
            last = None
            for attempt in range(retries + 1):
                try:
                    urllib.request.urlretrieve(url, zpath)
                    break
                except Exception as exc:                       # noqa: BLE001
                    last = exc
                    if os.path.exists(zpath):
                        os.remove(zpath)
                    if attempt >= retries:
                        raise RuntimeError(
                            f"failed to download {url}: {exc}") from last
        with zipfile.ZipFile(zpath) as zf:
            zf.extractall(cache_dir)
    return apath


# --------------------------------------------------------------------------- #
# climatology field(s) over an AOI
# --------------------------------------------------------------------------- #
def _load_da(src: Union[str, xr.DataArray], clip_box=None) -> xr.DataArray:
    """Load a raster path (or pass a DataArray through) into memory and close
    the file promptly -- avoids the benign ``sys.excepthook`` GC error rioxarray
    throws when a dataset is left open until interpreter exit. ``clip_box`` is an
    optional ``(minx, miny, maxx, maxy)`` to crop before loading."""
    if not isinstance(src, str):
        return src
    with rioxarray.open_rasterio(src, masked=True) as da:
        if da.rio.crs is None:
            da = da.rio.write_crs(4269)
        if clip_box is not None:
            da = da.rio.clip_box(*clip_box)
        return da.load()


def climatology_field(aoi, durations: Sequence[int] = (15, 30, 60),
                      ari: int = 1, region: Optional[str] = None,
                      stat: str = "mean", series: str = "pds",
                      cache_dir: str = "atlas14_cache", pad_deg: float = 0.05,
                      verbose: bool = True) -> dict:
    """NOAA Atlas 14 intensity climatology over an AOI, as a result dict.

    For each duration (minutes) the matching gridded grid is fetched, cropped to
    the AOI bbox (+ ``pad_deg``), and converted from depth (mils) to intensity
    (mm/h). All durations are put on one grid (the first duration's), so the
    fields share a single transform/profile -- ``mrms.save_fields(result,
    out_dir, f"{key}_clim")`` then writes ``<key>_clim_i15.tif`` etc.
    """
    bounds, _ = load_aoi(aoi)
    region = region or region_for_bounds(bounds)
    w, s, e, n = bounds
    bb = (w - pad_deg, s - pad_deg, e + pad_deg, n + pad_deg)
    if stat not in STAT:
        raise ValueError(f"stat must be one of {sorted(STAT)}, got {stat!r}")
    scode = STAT[stat]
    if verbose:
        print(f"NOAA Atlas 14 region '{region}', {ari}-yr {stat}, "
              f"durations {tuple(durations)} min")

    ref = None
    fields = {}
    for d in durations:
        asc = fetch_grid(region, ari, d, stat=scode, series=series,
                         cache_dir=cache_dir, verbose=verbose)
        da = _load_da(asc, clip_box=bb)
        if ref is None:
            ref = da
        else:                              # keep every duration on the same grid
            da = da.rio.reproject_match(ref)
        arr = np.asarray(da.squeeze().values, dtype="float32")
        fields[f"i{d}"] = (arr / 1000.0) * MM_PER_INCH * (60.0 / d)

    transform = ref.rio.transform()
    h, w_ = fields[f"i{durations[0]}"].shape
    profile = dict(driver="GTiff", height=h, width=w_, count=1, dtype="float32",
                   crs="EPSG:4269", transform=transform, nodata=np.nan,
                   compress="LZW")
    meta = dict(source="NOAA Atlas 14", region=region, ari_yr=int(ari),
                stat=stat, series=series, durations_min=list(durations),
                **{f"i{d}_max_mmph": float(np.nanmax(fields[f"i{d}"]))
                   for d in durations})
    return dict(fields=fields, transform=transform, crs="EPSG:4269",
                profile=profile, meta=meta)


# --------------------------------------------------------------------------- #
# anomaly  =  observed / climatology
# --------------------------------------------------------------------------- #
def anomaly(observed: Union[str, xr.DataArray],
            clim: Union[str, xr.DataArray], eps: float = 0.1) -> xr.DataArray:
    """Cell-wise anomaly ratio ``observed / climatology`` as a DataArray.

    ``clim`` is reprojected/resampled onto the ``observed`` grid first (so the
    NAD83 climatology aligns with the EPSG:4326 radar field). Cells where the
    climatology is <= ``eps`` mm/h or the observation is missing become NaN.
    Write it with ``ratio.rio.to_raster(path)``.
    """
    obs = _load_da(observed)
    clim_on = _load_da(clim).rio.reproject_match(obs)
    ratio = obs / clim_on
    ratio = ratio.where(np.isfinite(obs) & (clim_on > eps))
    ratio = ratio.rio.write_nodata(np.nan)
    return ratio


# --------------------------------------------------------------------------- #
# point precipitation-frequency curve + recurrence interval
# --------------------------------------------------------------------------- #
PFDS_URL = "https://hdsc.nws.noaa.gov/cgi-bin/hdsc/new/fe_text_{stat}.csv"
# Atlas 14 partial-duration recurrence intervals (years); fallback if the
# response header is missing.
_DEFAULT_ARIS = (1, 2, 5, 10, 25, 50, 100, 200, 500, 1000)


def pf_point(lat: float, lon: float, durations: Sequence[int] = (15, 30, 60),
             stat: str = "mean", series: str = "pds",
             timeout: float = 30.0) -> pd.DataFrame:
    """NOAA **PFDS point** precipitation-frequency curve at one location.

    Queries NOAA's Hydrometeorological Design Studies Center point service --
    the authoritative source the precipitation-frequency *maps* are tiled from
    (the same data ``pfdf.data.noaa.atlas14`` wraps) -- and returns a DataFrame
    of **intensity (mm/h)** indexed by recurrence interval (ARI, years) with one
    column ``i{d}`` per requested duration. ``stat`` is ``mean`` (default) /
    ``lower`` / ``upper`` (90% CI bounds); ``series`` is ``pds`` (partial
    duration -- matches the gridded :func:`climatology_field`) or ``ams``
    (annual maximum).

    Unlike the gridded product (one ARI per download), this returns the **whole
    ARI curve** in one request -- needed to invert an observed intensity into a
    recurrence interval (:func:`recurrence_interval`). The point estimate is more
    accurate than sampling the ~800 m grid, and needs no region lookup (the
    service locates the volume from ``lat``/``lon``).
    """
    if stat not in STAT:
        raise ValueError(f"stat must be one of {list(STAT)}, got {stat!r}")
    url = (PFDS_URL.format(stat=stat) +
           f"?lat={float(lat):.4f}&lon={float(lon):.4f}"
           f"&data=intensity&units=metric&series={series}")
    with urllib.request.urlopen(url, timeout=timeout) as resp:   # noqa: S310
        text = resp.read().decode("utf-8", "replace")
    aris, rows = None, {}
    want = {f"{int(d)}-min": int(d) for d in durations}
    for ln in text.splitlines():
        head = ln.split(":")[0].strip()
        if aris is None and "ARI" in ln and "year" in ln.lower():
            aris = [float(x) for x in ln.split(":", 1)[1].split(",") if x.strip()]
        if head in want:
            rows[want[head]] = [float(x) for x in
                                ln.split(":", 1)[1].split(",") if x.strip()]
    if not rows:
        raise ValueError(f"no PF estimates parsed from PFDS for ({lat}, {lon})")
    aris = aris or list(_DEFAULT_ARIS)
    return pd.DataFrame({f"i{d}": rows[d] for d in durations if d in rows},
                        index=pd.Index(aris, name="ari_yr"))


def recurrence_interval(intensity: float, aris: Sequence[float],
                        curve: Sequence[float]) -> float:
    """Recurrence interval (years) of an observed ``intensity`` against a PF
    ``curve`` (intensity at each ARI in ``aris``, monotonically increasing).

    Atlas 14 publishes the fitted quantile *curve*, not a closed-form inverse, so
    the standard practice is to interpolate the observed value against it. This
    uses **log-log** interpolation (log intensity vs log ARI). **Below the 1-yr
    quantile it extrapolates** the lowest log-log segment -- the partial-duration
    series admits sub-annual recurrence (more than one event per year), so values
    near and below 1 yr are valid (e.g. a 0.5-yr "twice-a-year" intensity). Returns
    ``inf`` above the largest tabulated ARI, and ``nan`` only for non-positive /
    missing input or a degenerate curve.
    """
    inten = np.asarray(curve, float)
    a = np.asarray(aris, float)
    if not np.isfinite(intensity) or intensity <= 0 or inten.size < 2:
        return float("nan")
    if intensity > inten[-1]:
        return float("inf")
    li, la = np.log(inten), np.log(a)
    if intensity < inten[0]:                      # extrapolate below the 1-yr point
        denom = li[1] - li[0]
        if denom == 0:
            return float("nan")
        slope = (la[1] - la[0]) / denom
        return float(np.exp(la[0] + (np.log(intensity) - li[0]) * slope))
    return float(np.exp(np.interp(np.log(intensity), li, la)))
