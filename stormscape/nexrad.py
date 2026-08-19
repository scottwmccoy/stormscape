"""Single-radar NEXRAD **Level II** access for an AOI (true elevation tilts).

Where :mod:`stormscape.mrms` pulls the gridded NOAA **MRMS** *mosaic* (already
quality-controlled and blended across radars into a 2-min QPE), this module
reaches the **raw single-radar volumes** -- the WSR-88D Level II archive -- so
you can look at reflectivity / velocity at each elevation tilt for the radar
nearest an AOI. It is the in-suite complement to MRMS: MRMS answers *"what was
the blended rain rate?"*, Level II answers *"what did the nearest radar actually
see at the lowest tilt?"* -- which is what's needed to check whether MRMS's
convective/hail Z-R typing (the driver of the radar-vs-gauge i15 over-read) is
justified by the underlying reflectivity.

Transport is the **nexradaws** package (v2.0+, which points at the current
``unidata-nexrad-level2`` S3 bucket -- the older ``noaa-nexrad-level2`` Big-Data
bucket we first tried by hand was deprecated). Reading the volumes uses
**Py-ART**. Both are optional dependencies (``pip install stormscape[nexrad]`` /
``conda install -c conda-forge arm_pyart`` + ``pip install nexradaws``); they are
imported lazily so the rest of the toolkit works without them.

Outputs mirror :mod:`stormscape.mrms`: gridded fields come back as a result
``dict`` (``fields / transform / crs / profile / meta``) in EPSG:4326, ready for
:func:`stormscape.mrms.save_fields` and :func:`stormscape.plot.drape_i15`.

Archive note: the Level II record on AWS runs from roughly **1991 to present**
(coverage is spotty before the early 2000s); the 2-min MRMS cadence this toolkit
otherwise relies on only starts 2020-10-14, so Level II reaches *further back*
than MRMS for older events.
"""

from __future__ import annotations

import os
import re
import warnings
from datetime import datetime, timezone
from importlib import resources
from typing import List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
import rioxarray  # noqa: F401  (registers the .rio accessor used in _to_result)
import xarray as xr

from .aoi import load_aoi

# pyart field name -> short tag used in output filenames (<key>_<tag>.tif)
_SHORT = {
    "reflectivity": "refl",
    "velocity": "vel",
    "spectrum_width": "sw",
    "differential_reflectivity": "zdr",
    "differential_phase": "phidp",
    "cross_correlation_ratio": "rhohv",
}

# Z = a R**b reflectivity-rate relations (Z in mm^6 m^-3, R in mm/h). Invert to
# R = (z/a)**(1/b) with z = 10**(dBZ/10). The WSR-88D operational default is the
# convective (300, 1.4); the others are offered for sensitivity tests.
RATE_RELATIONS = {
    "convective": (300.0, 1.4),    # WSR-88D operational default (Z = 300 R^1.4)
    "stratiform": (200.0, 1.6),    # Marshall-Palmer
    "tropical": (250.0, 1.2),      # Rosenfeld tropical
    "marshall_palmer": (200.0, 1.6),
}

# Volume count above which download_scans() warns that the window looks like a
# whole storm day rather than a storm. Level II is ~10 volumes/hour, so this is
# roughly "more than half a day of scans".
BULK_SCAN_WARN = 120


# --------------------------------------------------------------------------- #
# NEXRAD site table  (data/nexrad_sites.csv, from NCEI HOMR)
# --------------------------------------------------------------------------- #
_SITES: Optional[pd.DataFrame] = None


def _sites() -> pd.DataFrame:
    """The bundled WSR-88D site table (id, lat, lon, elev_m), cached."""
    global _SITES
    if _SITES is None:
        with resources.files("stormscape").joinpath(
                "data/nexrad_sites.csv").open() as fh:
            _SITES = pd.read_csv(fh)
    return _SITES


def radar_location(radar_id: str) -> Tuple[float, float, float]:
    """``(lat, lon, elev_m)`` for a four-letter radar id (e.g. ``"KRGX"``)."""
    s = _sites()
    row = s[s["id"] == radar_id.upper()]
    if not len(row):
        raise KeyError(f"unknown NEXRAD radar id {radar_id!r}")
    return (float(row["lat"].iloc[0]), float(row["lon"].iloc[0]),
            float(row["elev_m"].iloc[0]))


def _haversine_km(lat1, lon1, lat2, lon2) -> np.ndarray:
    """Great-circle distance (km); scalars or arrays for the second point."""
    r = 6371.0
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dphi = np.radians(np.asarray(lat2) - lat1)
    dlmb = np.radians(np.asarray(lon2) - lon1)
    a = (np.sin(dphi / 2) ** 2
         + np.cos(p1) * np.cos(p2) * np.sin(dlmb / 2) ** 2)
    return 2 * r * np.arcsin(np.sqrt(a))


def nearest_radar(aoi) -> Tuple[str, float, float, float]:
    """Closest WSR-88D to the AOI centroid -> ``(id, dist_km, lat, lon)``."""
    bounds, _ = load_aoi(aoi)
    clat = (bounds[1] + bounds[3]) / 2.0
    clon = (bounds[0] + bounds[2]) / 2.0
    s = _sites()
    d = _haversine_km(clat, clon, s["lat"].to_numpy(), s["lon"].to_numpy())
    i = int(np.argmin(d))
    return (str(s["id"].iloc[i]), float(d[i]),
            float(s["lat"].iloc[i]), float(s["lon"].iloc[i]))


def _radar_id(radar) -> str:
    """Resolve a radar argument to an id: a 4-letter id passes through; an AOI
    (tuple / geometry / vector path) resolves to its :func:`nearest_radar`."""
    if isinstance(radar, str) and re.fullmatch(r"[A-Za-z]{4}", radar):
        rid = radar.upper()
        if rid in set(_sites()["id"]):
            return rid
    return nearest_radar(radar)[0]


# --------------------------------------------------------------------------- #
# transport (nexradaws -> unidata-nexrad-level2)
# --------------------------------------------------------------------------- #
def _conn():
    try:
        import nexradaws
    except ImportError as e:                                   # optional dep
        raise ImportError(
            "nexradaws is required for NEXRAD access -- `pip install nexradaws` "
            "(or install the 'nexrad' extra).") from e
    return nexradaws.NexradAwsInterface()


def _utc(t: datetime) -> datetime:
    """Coerce a datetime to tz-aware UTC (nexradaws/pyart compare tz-aware)."""
    return t.replace(tzinfo=timezone.utc) if t.tzinfo is None \
        else t.astimezone(timezone.utc)


def _iso(t: datetime) -> str:
    return _utc(t).replace(microsecond=0).isoformat()


def available_scans(radar, start: datetime, end: datetime) -> list:
    """List archived volume scans for ``radar`` over ``[start, end]`` (UTC).

    ``radar`` is a four-letter id or any AOI (resolved to the nearest radar).
    Returns nexradaws ``AwsNexradFile`` objects (``.scan_time``, ``.filename``,
    ``.key`` ...) with the ``_MDM`` metadata-only entries filtered out.
    """
    rid = _radar_id(radar)
    scans = _conn().get_avail_scans_in_range(_utc(start), _utc(end), rid)
    return [s for s in scans if not str(s.filename).endswith("_MDM")]


def download_scans(scans: Sequence, cache_dir: str = "nexrad_cache") -> List[str]:
    """Download ``scans`` to ``cache_dir`` (flat); return local file paths.

    Already-present files are not re-fetched, so repeated runs over the same
    storm window reuse the cache.
    """
    os.makedirs(cache_dir, exist_ok=True)
    have, todo = [], []
    for s in scans:
        p = os.path.join(cache_dir, s.filename)
        (have if os.path.exists(p) else todo).append(p if os.path.exists(p) else s)
    if len(todo) > BULK_SCAN_WARN:
        # Level II runs ~10 volumes/hour at ~10 MB each, so a storm-DAY window is
        # ~2 GB and tens of minutes of Py-ART time -- nearly all of it dry. Bound
        # the window first with stormscape.mrms.wet_window (hourly QPE, a few KB
        # per hour) or gauges.storm_window. Warn rather than refuse: some callers
        # legitimately want a long span.
        warnings.warn(
            f"downloading {len(todo)} NEXRAD volumes (~{len(todo) * 10 / 1024:.1f} GB); "
            "if this is a whole storm day, bound the window with "
            "stormscape.mrms.wet_window() first -- the CLI does this automatically",
            stacklevel=2)
    if todo:
        res = _conn().download(todo, cache_dir)
        have.extend(lf.filepath for lf in res.iter_success())
        if res.failed_count:
            warnings.warn(f"{res.failed_count} NEXRAD download(s) failed")
    return sorted(have)


def nearest_scan(radar, when: datetime, search_min: float = 20.0):
    """The single ``AwsNexradFile`` whose scan time is closest to ``when``."""
    import datetime as _dt
    when = _utc(when)
    scans = available_scans(radar, when - _dt.timedelta(minutes=search_min),
                            when + _dt.timedelta(minutes=search_min))
    if not scans:
        raise ValueError(f"no scans within {search_min:g} min of {when:%Y-%m-%d %H:%MZ}")
    return min(scans, key=lambda s: abs((_utc(s.scan_time) - when).total_seconds()))


def read_sweep(scan: Union[str, object]):
    """Open a downloaded Level II volume -> a Py-ART ``Radar`` object.

    Accepts a local path, a nexradaws ``LocalNexradFile``, or anything with a
    ``.filepath``. (Pass an undownloaded ``AwsNexradFile`` through
    :func:`download_scans` first.)
    """
    try:
        import pyart
    except ImportError as e:                                   # optional dep
        raise ImportError(
            "Py-ART is required to read Level II -- "
            "`conda install -c conda-forge arm_pyart`.") from e
    if hasattr(scan, "open_pyart"):
        return scan.open_pyart()
    path = scan if isinstance(scan, str) else getattr(scan, "filepath", None)
    if not path:
        raise TypeError(f"cannot read {scan!r}; pass a path or LocalNexradFile")
    return pyart.io.read_nexrad_archive(path)


# --------------------------------------------------------------------------- #
# polar sweep -> AOI raster (azimuthal-equidistant about the radar -> 4326)
# --------------------------------------------------------------------------- #
def _aeqd(rlat: float, rlon: float) -> str:
    return (f"+proj=aeqd +lat_0={rlat} +lon_0={rlon} +x_0=0 +y_0=0 "
            "+datum=WGS84 +units=m +no_defs")


def _extent_for(rlat, rlon, bounds, res_m, margin_m=3000.0):
    """Regular x/y grid (m, AEQD about the radar) spanning the AOI bbox."""
    from pyproj import Transformer
    crs = _aeqd(rlat, rlon)
    tf = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
    w, s, e, n = bounds
    xc, yc = tf.transform([w, e, w, e], [s, s, n, n])
    xs = np.arange(min(xc) - margin_m, max(xc) + margin_m + res_m, res_m)
    ys = np.arange(min(yc) - margin_m, max(yc) + margin_m + res_m, res_m)
    return xs, ys, crs


def _grid_sweep(radar, field, sweep, xs, ys):
    """Rasterise one sweep onto the (xs, ys) AEQD grid by nearest gate.

    Gates are placed at their ground (x, y) about the radar; each grid cell takes
    the value of the nearest gate, but only if that gate is within the local
    beam spacing (~1 deg azimuthally) so the field is masked beyond real data
    instead of smeared across the whole disk.
    """
    from scipy.spatial import cKDTree
    if field not in radar.fields:
        raise ValueError(f"field {field!r} not in scan; have {sorted(radar.fields)}")
    one = radar.extract_sweeps([sweep])
    gx, gy, _ = one.get_gate_x_y_z(0)
    data = np.ma.filled(one.fields[field]["data"].astype("float64"), np.nan)
    ok = np.isfinite(gx) & np.isfinite(gy)
    tree = cKDTree(np.column_stack([gx[ok].ravel(), gy[ok].ravel()]))
    vals = data[ok].ravel()
    xx, yy = np.meshgrid(xs, ys)
    q = np.column_stack([xx.ravel(), yy.ravel()])
    dist, idx = tree.query(q, k=1)
    out = vals[idx]
    res = float(xs[1] - xs[0])
    rng = np.hypot(q[:, 0], q[:, 1])
    out[dist > np.maximum(1.5 * res, rng * np.radians(1.0))] = np.nan
    return out.reshape(yy.shape)


def _to_result(arr2d, xs, ys, aeqd_crs, bounds, meta, short):
    """AEQD array + grid -> EPSG:4326 result dict (clipped to the AOI bbox)."""
    da = xr.DataArray(arr2d.astype("float32"),
                      coords={"y": ys, "x": xs}, dims=("y", "x"))
    da = da.sortby("y", ascending=False).rio.write_crs(aeqd_crs)
    da = da.rio.write_nodata(np.nan).rio.reproject("EPSG:4326")
    try:
        da = da.rio.clip_box(*bounds)
    except Exception:                                          # noqa: BLE001
        pass
    arr = np.asarray(da.values, dtype="float32")
    if arr.ndim == 3:
        arr = arr[0]
    tr = da.rio.transform()
    profile = dict(driver="GTiff", height=arr.shape[0], width=arr.shape[1],
                   count=1, dtype="float32", crs="EPSG:4326", transform=tr,
                   nodata=float("nan"))
    return dict(fields={short: arr}, transform=tr, crs="EPSG:4326",
                profile=profile, meta=meta)


def _radar_meta_id(radar) -> str:
    name = getattr(radar, "metadata", {}).get("instrument_name", "")
    if isinstance(name, bytes):
        name = name.decode(errors="ignore")
    return str(name).strip() or "?"


def lowest_tilt_grid(radar, aoi, field: str = "reflectivity", sweep: int = 0,
                     res_m: float = 500.0) -> dict:
    """Grid one elevation tilt of a Py-ART ``Radar`` over the AOI.

    ``sweep=0`` is the lowest tilt (~0.5 deg) -- the radar's closest look at the
    surface, the relevant level for QPE. Returns a result ``dict`` in EPSG:4326
    (``fields={'refl': ...}`` for reflectivity) compatible with
    :func:`stormscape.mrms.save_fields` / :func:`stormscape.plot.drape_i15`.
    """
    bounds, _ = load_aoi(aoi)
    rlat = float(radar.latitude["data"][0])
    rlon = float(radar.longitude["data"][0])
    xs, ys, crs = _extent_for(rlat, rlon, bounds, res_m)
    arr = _grid_sweep(radar, field, sweep, xs, ys)
    elev = float(radar.extract_sweeps([sweep]).fixed_angle["data"][0])
    meta = {"source": "NEXRAD Level II (unidata-nexrad-level2)",
            "radar": _radar_meta_id(radar), "field": field, "sweep": int(sweep),
            "elevation_deg": round(elev, 2), "res_m": float(res_m),
            "scan_time": radar.time["units"].split("since")[-1].strip()}
    return _to_result(arr, xs, ys, crs, bounds, meta, _SHORT.get(field, field))


# --------------------------------------------------------------------------- #
# convenience: AOI + time -> gridded field
# --------------------------------------------------------------------------- #
def reflectivity_field(aoi, when: datetime, radar: Optional[str] = None,
                       field: str = "reflectivity", sweep: int = 0,
                       res_m: float = 500.0, cache_dir: str = "nexrad_cache",
                       search_min: float = 20.0) -> dict:
    """One scan nearest ``when``: nearest radar -> download -> grid lowest tilt.

    The single-scan analogue of :func:`stormscape.mrms.i15_storm_day`.
    """
    rid = radar.upper() if isinstance(radar, str) else nearest_radar(aoi)[0]
    scan = nearest_scan(rid, when, search_min=search_min)
    path = download_scans([scan], cache_dir)[0]
    res = lowest_tilt_grid(read_sweep(path), aoi, field=field, sweep=sweep,
                           res_m=res_m)
    res["meta"].update(radar=rid, scan_time=_iso(scan.scan_time),
                       requested=_iso(when))
    return res


def reflectivity_composite(aoi, start: datetime, end: datetime,
                           radar: Optional[str] = None,
                           field: str = "reflectivity", sweep: int = 0,
                           res_m: float = 500.0, cache_dir: str = "nexrad_cache",
                           max_scans: int = 40) -> dict:
    """Per-cell **maximum** of the lowest tilt over every scan in the window.

    The "storm-peak reflectivity" map -- the single-radar analogue of the i15
    running max, useful for spotting convective/hail cores. Scans are capped at
    ``max_scans`` (evenly subsampled) to bound the work.
    """
    rid = radar.upper() if isinstance(radar, str) else nearest_radar(aoi)[0]
    scans = available_scans(rid, start, end)
    if not scans:
        raise ValueError(f"no {rid} scans in [{_iso(start)}, {_iso(end)}]")
    if max_scans and len(scans) > max_scans:
        scans = scans[::int(np.ceil(len(scans) / max_scans))]
    paths = download_scans(scans, cache_dir)
    bounds, _ = load_aoi(aoi)
    rlat, rlon, _ = radar_location(rid)
    xs, ys, crs = _extent_for(rlat, rlon, bounds, res_m)
    acc, n = None, 0
    for p in paths:
        try:
            arr = _grid_sweep(read_sweep(p), field, sweep, xs, ys)
        except Exception as exc:                              # noqa: BLE001
            warnings.warn(f"skip {os.path.basename(p)}: {exc}")
            continue
        acc = arr if acc is None else np.fmax(acc, arr)
        n += 1
    if acc is None:
        raise ValueError("no usable scans in window")
    meta = {"source": "NEXRAD Level II (unidata-nexrad-level2)", "radar": rid,
            "field": field, "sweep": int(sweep), "res_m": float(res_m),
            "composite": "max", "n_scans": n,
            "start": _iso(start), "end": _iso(end)}
    return _to_result(acc, xs, ys, crs, bounds, meta,
                      _SHORT.get(field, field) + "max")


# --------------------------------------------------------------------------- #
# point sampling + Z-R  (radar-vs-gauge diagnostic)
# --------------------------------------------------------------------------- #
def sample_radar_at_points(radar, points, field: str = "reflectivity",
                           sweep: int = 0) -> np.ndarray:
    """Sample one sweep at point locations by nearest gate -> ndarray.

    ``points`` is a GeoDataFrame (any CRS) or an ``(N, 2)`` lon/lat array. Reads
    the value of the closest radar gate to each point on the chosen tilt;
    no-echo gates come back as NaN. Pairs with :func:`z_to_rate` to convert
    sampled reflectivity to an implied rain rate at each gauge.
    """
    from scipy.spatial import cKDTree
    if field not in radar.fields:
        raise ValueError(f"field {field!r} not in scan; have {sorted(radar.fields)}")
    one = radar.extract_sweeps([sweep])
    glat, glon, _ = one.get_gate_lat_lon_alt(0)
    vals = np.ma.filled(one.fields[field]["data"].astype("float64"), np.nan).ravel()
    if hasattr(points, "geometry"):
        pts = points.to_crs(4326)
        lons = pts.geometry.x.to_numpy(float)
        lats = pts.geometry.y.to_numpy(float)
    else:
        arr = np.asarray(points, dtype=float)
        lons, lats = arr[:, 0], arr[:, 1]
    cos = np.cos(np.radians(float(np.nanmean(lats)) if len(lats) else 0.0))
    tree = cKDTree(np.column_stack([glon.ravel() * cos, glat.ravel()]))
    _, idx = tree.query(np.column_stack([lons * cos, lats]), k=1)
    return vals[idx]


def z_to_rate(dbz, a: float = 300.0, b: float = 1.4,
              dbz_cap: Optional[float] = None) -> np.ndarray:
    """Convert reflectivity (dBZ) to rain rate (mm/h) via ``Z = a R**b``.

    Defaults are the WSR-88D operational convective relation (300, 1.4); pick a
    pair from :data:`RATE_RELATIONS` for other regimes. ``dbz_cap`` (e.g. 53)
    applies the operational hail cap before conversion -- the difference between
    capped and uncapped rate at a high-dBZ cell is exactly the hail
    over-estimation that inflates radar-only QPE.
    """
    dbz = np.asarray(dbz, dtype=float)
    if dbz_cap is not None:
        dbz = np.minimum(dbz, dbz_cap)
    z = np.power(10.0, dbz / 10.0)
    return np.power(z / a, 1.0 / b)


# --------------------------------------------------------------------------- #
# v1 intensity stack: Level II -> rate-per-scan -> i15/i30/i60 (MRMS-style)
# --------------------------------------------------------------------------- #
def _scan_epoch(radar) -> float:
    """Volume-start time of a Py-ART radar as POSIX seconds (UTC)."""
    iso = radar.time["units"].split("since")[-1].strip().replace("Z", "+00:00")
    return datetime.fromisoformat(iso).timestamp()


def _low_sweeps(radar, tol: float = 0.25, sweep: Optional[int] = None):
    """[(sweep_index, time_offset_s)] for the lowest-tilt reflectivity cut(s).

    With one cut per volume this is just the base scan; under SAILS the lowest
    elevation is revisited mid-volume, so this returns *each* low cut (deduped
    to >90 s apart, which drops the split-cut Doppler twin) -- recovering the
    ~2.5-3 min effective low-level cadence. ``sweep`` forces a fixed index.
    """
    if sweep is not None:
        sl = radar.get_slice(sweep)
        return [(sweep, float(radar.time["data"][sl.start]))]
    fa = np.asarray(radar.fixed_angle["data"], dtype=float)
    lo = float(np.nanmin(fa))
    cand = []
    for s in range(radar.nsweeps):
        if fa[s] > lo + tol:
            continue
        sl = radar.get_slice(s)
        d = radar.fields["reflectivity"]["data"][sl]
        if np.isfinite(np.ma.filled(d, np.nan)).any():        # has real echo room
            cand.append((s, float(radar.time["data"][sl.start])))
    cand.sort(key=lambda x: x[1])
    out = []
    for s, off in cand:
        if out and off - out[-1][1] < 90.0:                   # drop split-cut twin
            continue
        out.append((s, off))
    return out


def _stack_to_intensities(times: np.ndarray, rates: np.ndarray,
                          durations=(15, 30, 60)) -> dict:
    """Per-scan rate grids -> running-max i15/i30/i60 + total + peak (per cell).

    Mirrors :func:`stormscape.mrms.compute_i15` / the gauge reducer: the rate
    series is interpolated to a regular 1-minute grid (handles the irregular
    ~3-5 min single-radar cadence the same way the gauge side handles irregular
    tips), then trailing-window accumulations give i15 = mean(i16, i14) and plain
    i30/i60, maximised over the storm. No-echo cells are treated as 0 mm/h
    (assumes the AOI is within radar coverage -- a v2 beam-blockage mask will
    relax this).
    """
    from scipy.interpolate import interp1d
    r0 = np.nan_to_num(rates, nan=0.0)
    tg = np.arange(times[0], times[-1] + 1.0, 60.0)
    rm = interp1d(times, r0, axis=0, bounds_error=False, fill_value=0.0)(tg)
    m = rm.shape[0]
    am = rm / 60.0                                            # mm per minute
    s = np.concatenate([np.zeros((1,) + rm.shape[1:]), np.cumsum(am, axis=0)])
    nan2d = np.full(rm.shape[1:], np.nan, dtype="float32")
    out = {}
    for d in durations:
        if d == 15:
            if m > 16:
                i16 = (s[16:] - s[:-16]) * 60.0 / 16.0
                i14 = (s[16:] - s[2:m - 13]) * 60.0 / 14.0
                out["i15max"] = np.nanmax((i16 + i14) / 2.0, axis=0)
            else:
                out["i15max"] = nan2d
        elif m > d:
            out[f"i{d}max"] = np.nanmax((s[d:] - s[:-d]) * 60.0 / d, axis=0)
        else:
            out[f"i{d}max"] = nan2d
    out["total_mm"] = s[-1]
    out["peakrate_mmph"] = np.nanmax(r0, axis=0)
    return out


def _multi_to_result(fields: dict, xs, ys, aeqd_crs, bounds, meta) -> dict:
    """Several AEQD field grids -> one EPSG:4326 result dict (shared grid)."""
    names = list(fields)
    arr = np.stack([np.asarray(fields[n], dtype="float32") for n in names])
    da = xr.DataArray(arr, coords={"band": names, "y": ys, "x": xs},
                      dims=("band", "y", "x"))
    da = da.sortby("y", ascending=False).rio.write_crs(aeqd_crs)
    da = da.rio.write_nodata(np.nan).rio.reproject("EPSG:4326")
    try:
        da = da.rio.clip_box(*bounds)
    except Exception:                                          # noqa: BLE001
        pass
    tr = da.rio.transform()
    h, w = da.shape[1], da.shape[2]
    profile = dict(driver="GTiff", height=h, width=w, count=1, dtype="float32",
                   crs="EPSG:4326", transform=tr, nodata=float("nan"))
    out_fields = {n: np.asarray(da.sel(band=n).values, dtype="float32")
                  for n in names}
    return dict(fields=out_fields, transform=tr, crs="EPSG:4326",
                profile=profile, meta=meta)


# S-band R(Kdp) = alpha * Kdp**beta (Bringi & Chandrasekar 2001). Kdp is derived
# from PhiDP with Py-ART's variational kdp_maesaka, which does its own PhiDP
# regularisation -- NEXRAD PhiDP is raw (system offset + folding) and the Z-PHI
# specific-attenuation path can't take it directly (ΔΦDP inflates, A blows up).
# Kdp is clipped to a physical S-band range to suppress derivative noise spikes.
_RKDP_S = (44.0, 0.822)                                         # alpha, beta
_KDP_CLIP = (0.0, 7.0)                                          # deg/km

# WSR-88D operational dual-pol R(Z, ZDR) for rain (Giangrande & Ryzhkov 2008):
# R = 0.0142 Z^0.770 Zdr_lin^-1.67  (Z linear mm6/m3, Zdr as a linear ratio).
# The negative ZDR exponent is the point: at a given Z, larger drops (higher
# ZDR) mean fewer drops and LESS rain -- the correction a fixed Z-R cannot make.
_RZZDR = (0.0142, 0.770, -1.67)
_ZDR_RAIN_DB = (0.0, 4.0)       # rain-plausible ZDR; outside = noise/cal issues


def zzdr_to_rate(dbz, zdr_db, dbz_cap: Optional[float] = None) -> np.ndarray:
    """Rain rate (mm/h) from reflectivity + differential reflectivity.

    The operational WSR-88D dual-pol relation (:data:`_RZZDR`). ZDR is clipped
    to the rain-plausible window (:data:`_ZDR_RAIN_DB`); where ZDR is not
    finite the caller should fall back to Z-R (this function returns NaN
    there). ``dbz_cap`` applies the same hail cap as :func:`z_to_rate`.
    """
    dbz = np.asarray(dbz, dtype=float)
    zdr = np.asarray(zdr_db, dtype=float)
    if dbz_cap is not None:
        dbz = np.minimum(dbz, dbz_cap)
    z = np.power(10.0, dbz / 10.0)
    zdr_lin = np.power(10.0, np.clip(zdr, *_ZDR_RAIN_DB) / 10.0)
    c, az, bz = _RZZDR
    out = c * np.power(z, az) * np.power(zdr_lin, bz)
    return np.where(np.isfinite(zdr), out, np.nan)


def _rates_zzdr(radar, low, xs, ys, a, b, dbz_cap, z_blend):
    """v3 rate grids: R(Z, ZDR) where moderate/heavy, capped Z-R elsewhere.

    Same blend policy as :func:`_rates_kdp` -- dual-pol relations are noisy in
    light rain, so R(Z,ZDR) applies only at/above ``z_blend`` dBZ. ZDR is a
    directly recorded field, so unlike Kdp no retrieval step is needed;
    single-pol volumes fall back to pure capped Z-R.
    """
    sub = radar.extract_sweeps([s for s, _ in low])
    out = []
    for k in range(sub.nsweeps):
        dbz = _grid_sweep(sub, "reflectivity", k, xs, ys)
        rz = z_to_rate(dbz, a=a, b=b, dbz_cap=dbz_cap)
        if "differential_reflectivity" not in sub.fields:
            out.append(rz)
            continue
        zdr = _grid_sweep(sub, "differential_reflectivity", k, xs, ys)
        rzz = zzdr_to_rate(dbz, zdr, dbz_cap=dbz_cap)
        use = np.isfinite(rzz) & (dbz >= z_blend)
        out.append(np.where(use, rzz, rz))
    return out


# Hydro blend: per-gate relation selection (CSU-HIDRO / WSR-88D dual-pol QPE
# style decision tree). Class codes written to the `relmode` product.
HYDRO_DRY, HYDRO_ZR, HYDRO_ZZDR, HYDRO_KDP_HAIL, HYDRO_CENSORED = 0, 1, 2, 3, 4
_HAIL_DBZ = 45.0     # hail suspicion: strong Z ...
_HAIL_ZDR = 0.8      # ... with small ZDR (dB) -- tumbling ice looks isotropic
_RHOHV_MIN = 0.85    # below this the echo is not meteorological
_KDP_GUARD = 1.5     # cap R(Kdp) at this x capped Z-R (backscatter-phase guard)


def hydro_select(dbz, zdr, rhohv, r_z, r_kdp, r_zzdr, z_blend: float = 35.0):
    """Per-cell relation selection -> (rate, lo, hi, class).

    The decision tree, each branch validated on the Aug 2026 Stallion storms:

    - ``rhohv < 0.85``            -> censor to 0 (non-meteorological echo)
    - ``dbz < z_blend``           -> capped Z-R (dual-pol too noisy in light rain)
    - ``dbz >= 45 & zdr <= 0.8``  -> **hail**: R(Kdp) (ice-blind), capped at
      1.5x the capped Z-R -- the backscatter-phase guard, added because a
      marginal-Z R(Kdp) blow-up (4-6x both Z-based estimators) was observed
      where no hail signature existed
    - otherwise (rain, moderate+) -> R(Z,ZDR) (DSD-robust)

    ``lo``/``hi`` bound each cell by the relations *defensibly applicable*
    there: {Z-R, R(Z,ZDR)} in rain, {R(Kdp), Z-R} in hail (R(Z,ZDR) blows up
    in hail and is excluded), the single Z-R in light rain. A reconnaissance
    field that carries its own spread cannot quietly overstate confidence.
    """
    dbz = np.asarray(dbz, dtype="float64")
    fin = np.isfinite(dbz)
    rate = np.where(fin, r_z, np.nan)
    lo = rate.copy()
    hi = rate.copy()
    cls = np.where(fin, HYDRO_ZR, HYDRO_DRY).astype("uint8")

    zdr_ok = zdr is not None and np.isfinite(np.asarray(zdr, dtype="float64")).any()
    if zdr_ok:
        zdr = np.asarray(zdr, dtype="float64")
        m_rain = fin & (dbz >= z_blend) & np.isfinite(r_zzdr)
        rate = np.where(m_rain, r_zzdr, rate)
        lo = np.where(m_rain, np.fmin(r_z, r_zzdr), lo)
        hi = np.where(m_rain, np.fmax(r_z, r_zzdr), hi)
        cls = np.where(m_rain, HYDRO_ZZDR, cls).astype("uint8")

        m_hail = fin & (dbz >= _HAIL_DBZ) & np.isfinite(zdr) & (zdr <= _HAIL_ZDR)
        rk = np.where(np.isfinite(r_kdp) & (r_kdp > 0),
                      np.fmin(r_kdp, _KDP_GUARD * r_z), r_z)
        rate = np.where(m_hail, rk, rate)
        lo = np.where(m_hail, np.fmin(rk, r_z), lo)
        hi = np.where(m_hail, np.fmax(rk, r_z), hi)
        cls = np.where(m_hail, HYDRO_KDP_HAIL, cls).astype("uint8")

    if rhohv is not None:
        cc = np.asarray(rhohv, dtype="float64")
        m_cen = fin & np.isfinite(cc) & (cc < _RHOHV_MIN)
        rate = np.where(m_cen, 0.0, rate)
        lo = np.where(m_cen, 0.0, lo)
        hi = np.where(m_cen, 0.0, hi)
        cls = np.where(m_cen, HYDRO_CENSORED, cls).astype("uint8")
    return rate, lo, hi, cls


def _rates_hydro(radar, low, xs, ys, a, b, dbz_cap, z_blend):
    """Per-cut (rate, lo, hi, class) grids for the hydro blend.

    Single-pol volumes degrade to pure capped Z-R with a point envelope, so
    the same call works across eras.
    """
    sub = radar.extract_sweeps([s for s, _ in low])
    dual = "differential_reflectivity" in sub.fields
    kdp_ok = False
    if dual and "differential_phase" in sub.fields:
        try:
            import pyart
            kd = pyart.retrieve.kdp_maesaka(sub, psidp_field="differential_phase")[0]
            kd["data"] = np.clip(np.ma.filled(kd["data"], 0.0), *_KDP_CLIP)
            sub.add_field("kdp", kd, replace_existing=True)
            rr = pyart.retrieve.est_rain_rate_kdp(sub, alpha=_RKDP_S[0],
                                                  beta=_RKDP_S[1], kdp_field="kdp")
            sub.add_field("rate_kdp", rr, replace_existing=True)
            kdp_ok = True
        except Exception as exc:                               # noqa: BLE001
            warnings.warn(f"R(Kdp) failed ({exc}); hail branch falls to Z-R")
    out = []
    for k in range(sub.nsweeps):
        dbz = _grid_sweep(sub, "reflectivity", k, xs, ys)
        r_z = z_to_rate(dbz, a=a, b=b, dbz_cap=dbz_cap)
        if not dual:
            cls = np.where(np.isfinite(dbz), HYDRO_ZR, HYDRO_DRY).astype("uint8")
            out.append((r_z, r_z.copy(), r_z.copy(), cls))
            continue
        zdr = _grid_sweep(sub, "differential_reflectivity", k, xs, ys)
        r_zzdr = zzdr_to_rate(dbz, zdr, dbz_cap=dbz_cap)
        r_kdp = (_grid_sweep(sub, "rate_kdp", k, xs, ys) if kdp_ok
                 else np.full_like(dbz, np.nan))
        cc = (_grid_sweep(sub, "cross_correlation_ratio", k, xs, ys)
              if "cross_correlation_ratio" in sub.fields else None)
        out.append(hydro_select(dbz, zdr, cc, r_z, r_kdp, r_zzdr,
                                z_blend=z_blend))
    return out


def _rates_kdp(radar, low, xs, ys, a, b, dbz_cap, z_blend):
    """v2 dual-pol rate grids for the given low cuts (one per cut).

    Rain rate from **specific differential phase** R(Kdp) -- hail-robust, since
    Kdp responds to liquid not ice -- used where the cell is moderate/heavy
    (``Z >= z_blend``) and Kdp is reliably positive; capped convective Z-R fills
    the rest. Kdp comes from ``pyart.retrieve.kdp_maesaka`` (variational, takes
    raw NEXRAD PhiDP). Volumes without dual-pol (pre-~2012 ``PhiDP``) fall back to
    pure capped Z-R, so the same call works across eras.
    """
    sub = radar.extract_sweeps([s for s, _ in low])
    zr = [z_to_rate(_grid_sweep(sub, "reflectivity", k, xs, ys), a=a, b=b,
                    dbz_cap=dbz_cap) for k in range(sub.nsweeps)]
    if "differential_phase" not in sub.fields:                 # single-pol era
        return zr
    try:
        import pyart
        kdp = pyart.retrieve.kdp_maesaka(sub, psidp_field="differential_phase")[0]
        kdp["data"] = np.clip(np.ma.filled(kdp["data"], 0.0), *_KDP_CLIP)
        sub.add_field("kdp", kdp, replace_existing=True)
        rr = pyart.retrieve.est_rain_rate_kdp(sub, alpha=_RKDP_S[0],
                                              beta=_RKDP_S[1], kdp_field="kdp")
        sub.add_field("rate_kdp", rr, replace_existing=True)
    except Exception as exc:                                   # noqa: BLE001
        warnings.warn(f"R(Kdp) failed ({exc}); Z-R fallback for this volume")
        return zr
    out = []
    for k in range(sub.nsweeps):
        rk = _grid_sweep(sub, "rate_kdp", k, xs, ys)
        dbz = _grid_sweep(sub, "reflectivity", k, xs, ys)
        use = np.isfinite(rk) & (rk > 0) & (dbz >= z_blend)
        out.append(np.where(use, rk, zr[k]))
    return out


def intensity_stack(aoi, start: datetime, end: datetime,
                    radar: Optional[str] = None, a: float = 300.0,
                    b: float = 1.4, dbz_cap: Optional[float] = 53.0,
                    sweep: Optional[int] = None, res_m: float = 500.0,
                    durations=(15, 30, 60), cache_dir: str = "nexrad_cache",
                    elev_tol: float = 0.25, method: str = "za",
                    z_blend: float = 35.0, rate_cap: Optional[float] = None,
                    blockage_dem=None, cbb_max: float = 0.5) -> dict:
    """Single-radar **i15/i30/i60** stack from Level II -- the MRMS analogue.

    For every volume in ``[start, end]`` the lowest-tilt reflectivity (all SAILS
    low cuts) is gridded over the AOI and turned into a rain rate by a capped
    convective Z-R (:func:`z_to_rate`, ``Z = a R**b`` with the ``dbz_cap`` hail
    cap). The per-scan rate grids are stacked, interpolated to 1-minute, and
    reduced to peak ``i15`` (= mean of trailing 16/14-min), ``i30``, ``i60``,
    storm ``total_mm`` and ``peakrate_mmph`` -- the **same fields and result
    dict** as :func:`stormscape.mrms.i15_storm_day`, so it drops straight into
    ``save_fields`` / ``drape_i15`` / ``compare``.

    ``method`` selects the rate retrieval: ``"za"`` (v1) is the fixed capped
    convective Z-R -- one relation for every era; ``"kdp"`` (v2) uses dual-pol
    **R(Kdp)** where ``Z >= z_blend`` (hail-robust), falling back to capped Z-R in
    light rain and for pre-dual-pol volumes. ``rate_cap`` (mm/h) clips every
    method's per-scan rate to a physical maximum before stacking -- an operational
    hail-cap analogue that bounds the i15 peak the same way for every method.
    Pass ``blockage_dem`` (a DEM path) to
    mask cells whose cumulative DEM beam-blockage exceeds ``cbb_max`` and emit a
    ``cbb`` quality field (the single-radar analogue of MRMS RQI).
    """
    rid = _radar_id(radar) if radar is not None else nearest_radar(aoi)[0]
    scans = available_scans(rid, start, end)
    if not scans:
        raise ValueError(f"no {rid} scans in [{_iso(start)}, {_iso(end)}]")
    paths = download_scans(scans, cache_dir)
    bounds, _ = load_aoi(aoi)
    rlat, rlon, _ = radar_location(rid)
    xs, ys, crs = _extent_for(rlat, rlon, bounds, res_m)

    times, grids, nvol = [], [], 0
    glo, ghi = [], []
    cls_counts = None                       # (5, H, W) usage counts for relmode
    for p in paths:
        try:
            rad = read_sweep(p)
            if "reflectivity" not in rad.fields:
                continue
            base = _scan_epoch(rad)
            low = _low_sweeps(rad, elev_tol, sweep)
            if not low:
                continue
            if method == "hydro":
                quads = _rates_hydro(rad, low, xs, ys, a, b, dbz_cap, z_blend)
                if cls_counts is None:
                    cls_counts = np.zeros((5,) + quads[0][0].shape, "uint16")
                for (si, off), (rate, lo_g, hi_g, cl) in zip(low, quads):
                    grids.append(rate)
                    glo.append(lo_g)
                    ghi.append(hi_g)
                    times.append(base + off)
                    wetc = np.nan_to_num(rate, nan=0.0) > 1.0
                    for c in (HYDRO_ZR, HYDRO_ZZDR, HYDRO_KDP_HAIL):
                        cls_counts[c] += ((cl == c) & wetc)
                    cls_counts[HYDRO_CENSORED] += (cl == HYDRO_CENSORED)
                nvol += 1
                continue
            if method == "kdp":
                rates = _rates_kdp(rad, low, xs, ys, a, b, dbz_cap, z_blend)
            elif method == "zzdr":
                rates = _rates_zzdr(rad, low, xs, ys, a, b, dbz_cap, z_blend)
            else:
                rates = [z_to_rate(_grid_sweep(rad, "reflectivity", si, xs, ys),
                                   a=a, b=b, dbz_cap=dbz_cap) for si, _ in low]
            for (si, off), rate in zip(low, rates):
                grids.append(rate)
                times.append(base + off)
            nvol += 1
        except Exception as exc:                               # noqa: BLE001
            warnings.warn(f"skip {os.path.basename(p)}: {exc}")
    if len(times) < 2:
        raise ValueError("need >= 2 usable scans for an intensity stack")

    times = np.asarray(times, dtype=float)
    order = np.argsort(times, kind="stable")
    streams = {"": np.stack(grids)[order]}
    if method == "hydro":
        streams["lo"] = np.stack(glo)[order]
        streams["hi"] = np.stack(ghi)[order]
    times = times[order]
    uniq = np.unique(times)
    if len(uniq) < len(times):                                 # merge equal times
        streams = {k: np.stack([v[times == t].mean(axis=0) for t in uniq])
                   for k, v in streams.items()}
        times = uniq
    if rate_cap is not None:                                  # uniform max-rate cap
        streams = {k: np.minimum(v, rate_cap) for k, v in streams.items()}
    fields = _stack_to_intensities(times, streams[""], durations)
    if method == "hydro":
        for tag in ("lo", "hi"):
            f = _stack_to_intensities(times, streams[tag], durations)
            fields[f"i15max_{tag}"] = f["i15max"]
            fields[f"total_{tag}"] = f["total_mm"]
        if cls_counts is not None:
            used = cls_counts[1:]
            relmode = np.where(used.sum(axis=0) > 0,
                               used.argmax(axis=0) + 1, 0)
            fields["relmode"] = relmode.astype("float32")
    blk = {}
    if blockage_dem is not None:                               # DEM beam-blockage
        try:
            rad0 = read_sweep(paths[0])
            s0 = _low_sweeps(rad0, elev_tol, sweep)[0][0]
            cbb = _blockage_cbb(rad0, str(blockage_dem), s0, xs, ys)
            blocked = np.isfinite(cbb) & (cbb > cbb_max)
            for k in list(fields):
                fields[k] = np.where(blocked, np.nan, fields[k])
            fields["cbb"] = cbb.astype("float32")
            blk = {"blockage": "DEM beam-blockage (wradlib)", "cbb_max": cbb_max,
                   "dem": os.path.basename(str(blockage_dem)),
                   "blocked_cells": int(blocked.sum())}
        except Exception as exc:                               # noqa: BLE001
            warnings.warn(f"beam-blockage skipped: {exc}")
    cad = float(np.median(np.diff(times)) / 60.0) if len(times) > 1 else float("nan")
    src = {"kdp": "NEXRAD L2 i15 stack: R(Kdp) blended with capped Z-R",
           "zzdr": "NEXRAD L2 i15 stack: R(Z,ZDR) blended with capped Z-R",
           "hydro": "NEXRAD L2 i15 stack: per-gate hydro blend "
                    "(Z-R / R(Z,ZDR) / R(Kdp)-in-hail, CSU-HIDRO style)",
           }.get(method, "NEXRAD L2 i15 stack: capped convective Z-R")
    meta = {"source": src, "radar": rid, "method": method,
            "zr_a": a, "zr_b": b, "dbz_cap": dbz_cap, "rate_cap": rate_cap,
            "res_m": float(res_m), "n_volumes": nvol, "n_low_cuts": int(len(times)),
            "cadence_min": round(cad, 2), "start": _iso(start), "end": _iso(end),
            "aoi_i15max": round(float(np.nanmax(fields["i15max"])), 1)}
    if method == "kdp":
        meta.update(kdp_alpha=_RKDP_S[0], kdp_beta=_RKDP_S[1], z_blend_dbz=z_blend)
    elif method == "zzdr":
        meta.update(zzdr_c=_RZZDR[0], zzdr_zexp=_RZZDR[1], zzdr_zdrexp=_RZZDR[2],
                    z_blend_dbz=z_blend)
    elif method == "hydro":
        meta.update(z_blend_dbz=z_blend, hail_dbz=_HAIL_DBZ, hail_zdr=_HAIL_ZDR,
                    rhohv_min=_RHOHV_MIN, kdp_guard=_KDP_GUARD)
        if cls_counts is not None:
            tot = int(cls_counts[1:4].sum())
            if tot:
                meta["relation_usage_pct"] = {
                    "zr_light": round(100 * int(cls_counts[HYDRO_ZR].sum()) / tot, 1),
                    "zzdr_rain": round(100 * int(cls_counts[HYDRO_ZZDR].sum()) / tot, 1),
                    "kdp_hail": round(100 * int(cls_counts[HYDRO_KDP_HAIL].sum()) / tot, 1)}
    meta.update(blk)
    return _multi_to_result(fields, xs, ys, crs, bounds, meta)


# --------------------------------------------------------------------------- #
# v2 step 2: DEM beam-blockage (single-radar quality; the MRMS-RQI analogue)
# --------------------------------------------------------------------------- #
def _blockage_cbb(radar, dem_path, sweep, xs, ys, beamwidth=0.925):
    """Cumulative beam-blockage fraction (0-1) for one tilt, gridded to (xs, ys).

    Terrain is sampled from ``dem_path`` at each gate; partial blockage
    (wradlib) uses the half-power beam radius and the beam-centre height, then
    accumulates along each ray.
    """
    import rasterio
    import rasterio.warp
    import wradlib
    from pyproj import Transformer
    sub = radar.extract_sweeps([sweep])
    lat, lon, alt = sub.get_gate_lat_lon_alt(0)               # beam-centre MSL [m]
    rng = np.broadcast_to(np.asarray(sub.range["data"], float)[None, :], lat.shape)
    a = rng * np.radians(beamwidth / 2.0)                     # half-power radius [m]
    th = np.full(lat.size, np.nan)
    with rasterio.open(dem_path) as ds:
        w, s, e, n = rasterio.warp.transform_bounds(ds.crs, "EPSG:4326", *ds.bounds)
        flon, flat = lon.ravel(), lat.ravel()
        inb = (flon >= w) & (flon <= e) & (flat >= s) & (flat <= n)
        if inb.any():
            tf = Transformer.from_crs("EPSG:4326", ds.crs, always_xy=True)
            gx, gy = tf.transform(flon[inb], flat[inb])
            v = np.array([t[0] for t in ds.sample(np.column_stack([gx, gy]))],
                         dtype=float)
            if ds.nodata is not None:
                v = np.where(v == ds.nodata, np.nan, v)
            th[inb] = v
    th = th.reshape(lat.shape)
    pbb = np.nan_to_num(
        np.asarray(wradlib.qual.beam_block_frac(th, alt, a), dtype=float), nan=0.0)
    cbb = wradlib.qual.cum_beam_block_frac(pbb)
    sub.add_field("cbb", {"data": np.ma.masked_invalid(cbb)}, replace_existing=True)
    return _grid_sweep(sub, "cbb", 0, xs, ys)


def beam_blockage(radar, aoi, dem, sweep: int = 0, res_m: float = 500.0,
                  beamwidth: float = 0.925) -> dict:
    """Cumulative beam-blockage (0-1) over the AOI for ``radar`` at ``sweep``.

    Geometry-only (time-independent): a DEM-derived quality field flagging where
    terrain blocks the beam -- the single-radar analogue of MRMS RQI. Returns a
    result dict (``fields={'cbb': ...}``) in EPSG:4326.
    """
    bounds, _ = load_aoi(aoi)
    rlat = float(radar.latitude["data"][0])
    rlon = float(radar.longitude["data"][0])
    xs, ys, crs = _extent_for(rlat, rlon, bounds, res_m)
    cbb = _blockage_cbb(radar, str(dem), sweep, xs, ys, beamwidth)
    meta = {"source": "NEXRAD DEM beam-blockage (wradlib)",
            "radar": _radar_meta_id(radar), "sweep": int(sweep),
            "dem": os.path.basename(str(dem)), "res_m": float(res_m)}
    return _multi_to_result({"cbb": cbb}, xs, ys, crs, bounds, meta)


# --------------------------------------------------------------------------- #
# virtual gauges from single-radar Level II (the pre-2020 / MRMS-gap fallback)
# --------------------------------------------------------------------------- #
def virtual_gauge_timeseries(points, start: datetime, end: datetime,
                             radar: Optional[str] = None, method: str = "kdp",
                             durations=(5, 15, 30, 60), a: float = 300.0,
                             b: float = 1.4, dbz_cap: Optional[float] = 53.0,
                             z_blend: float = 35.0,
                             rate_cap: Optional[float] = None,
                             cache_dir: str = "nexrad_cache") -> dict:
    """Virtual-gauge rainfall time series from single-radar Level II.

    The NEXRAD fallback for :func:`stormscape.mrms.virtual_gauge_timeseries`
    (pre-2020 events / MRMS gaps): samples the lowest-tilt rate (``method`` ``za``
    capped Z-R or ``kdp`` R(Kdp), the same recipes as :func:`intensity_stack`) at
    each point per scan, interpolates to 1-minute, and reduces to ``i{d}``.
    Returns ``{name: DataFrame}`` with the same columns as the MRMS version minus
    the QPE overlay (a single radar has no gauge-corrected product).
    """
    import geopandas as gpd
    import pandas as pd
    from shapely.geometry import Point

    from .mrms import _as_points
    pts = _as_points(points)
    lons = [lo for _, lo, _ in pts]
    lats = [la for _, _, la in pts]
    rid = (_radar_id(radar) if radar is not None
           else nearest_radar((min(lons), min(lats), max(lons), max(lats)))[0])
    paths = download_scans(available_scans(rid, start, end), cache_dir)
    g = gpd.GeoDataFrame({"name": [n for n, _, _ in pts]},
                         geometry=[Point(lo, la) for _, lo, la in pts], crs=4326)

    rec = {n: [] for n, _, _ in pts}
    for p in paths:
        try:
            rad = read_sweep(p)
        except Exception:                                      # noqa: BLE001
            continue
        if "reflectivity" not in rad.fields:
            continue
        low = _low_sweeps(rad)
        if not low:
            continue
        sub = rad.extract_sweeps([s for s, _ in low])
        base = _scan_epoch(rad)
        kdp_ok = False
        if method in ("kdp", "hydro") and "differential_phase" in sub.fields:
            try:
                import pyart
                kd = pyart.retrieve.kdp_maesaka(sub, psidp_field="differential_phase")[0]
                kd["data"] = np.clip(np.ma.filled(kd["data"], 0.0), *_KDP_CLIP)
                sub.add_field("kdp", kd, replace_existing=True)
                kdp_ok = True
            except Exception:                                  # noqa: BLE001
                kdp_ok = False
        for k, (si, off) in enumerate(low):
            dbz = sample_radar_at_points(sub, g, field="reflectivity", sweep=k)
            rz = z_to_rate(dbz, a=a, b=b, dbz_cap=dbz_cap)
            if kdp_ok:
                kv = sample_radar_at_points(sub, g, field="kdp", sweep=k)
                rk = _RKDP_S[0] * np.where(np.isfinite(kv) & (kv > 0), kv, 0.0
                                           ) ** _RKDP_S[1]
                rate = np.where(np.isfinite(rk) & (rk > 0) & (dbz >= z_blend),
                                rk, rz)
            elif (method == "zzdr"
                  and "differential_reflectivity" in sub.fields):
                zv = sample_radar_at_points(
                    sub, g, field="differential_reflectivity", sweep=k)
                rzz = zzdr_to_rate(dbz, zv, dbz_cap=dbz_cap)
                rate = np.where(np.isfinite(rzz) & (dbz >= z_blend), rzz, rz)
            elif (method == "hydro"
                  and "differential_reflectivity" in sub.fields):
                zv = sample_radar_at_points(
                    sub, g, field="differential_reflectivity", sweep=k)
                rzz = zzdr_to_rate(dbz, zv, dbz_cap=dbz_cap)
                if kdp_ok:
                    kv = sample_radar_at_points(sub, g, field="kdp", sweep=k)
                    rk = _RKDP_S[0] * np.where(np.isfinite(kv) & (kv > 0),
                                               kv, 0.0) ** _RKDP_S[1]
                else:
                    rk = np.full_like(np.atleast_1d(dbz), np.nan, dtype=float)
                cc = (sample_radar_at_points(
                          sub, g, field="cross_correlation_ratio", sweep=k)
                      if "cross_correlation_ratio" in sub.fields else None)
                rate, _, _, _ = hydro_select(dbz, zv, cc, rz, rk, rzz,
                                             z_blend=z_blend)
            else:
                rate = rz
            if rate_cap is not None:
                rate = np.minimum(rate, rate_cap)
            for (name, _, _), val in zip(pts, np.atleast_1d(rate)):
                rec[name].append((base + off, float(val)))

    t1 = pd.date_range(start, end, freq="1min")
    out = {}
    for name, _, _ in pts:
        if len(rec[name]) < 2:
            continue
        s = pd.Series({pd.Timestamp(e, unit="s"): v
                       for e, v in sorted(rec[name])})
        s = s[~s.index.duplicated()].sort_index()
        rate1 = (s.reindex(s.index.union(t1)).interpolate("time")
                 .reindex(t1).fillna(0.0))
        total = (rate1 / 60.0).cumsum()
        df = pd.DataFrame({"rate_mmph": rate1, "total_mm": total})
        for d in durations:
            if d == 15:
                df["i15_mmph"] = (total.diff(16) * 60.0 / 16.0
                                  + total.diff(14) * 60.0 / 14.0) / 2.0
            else:
                df[f"i{d}_mmph"] = total.diff(d) * 60.0 / d
        out[name] = df
    return out
