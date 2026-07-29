"""Spatial smoothing of radar-derived rainfall-intensity fields.

The observed MRMS / single-radar NEXRAD ``i15/i30/i60`` fields are peaky (~1 km
native, sharp convective cores) next to the smooth NOAA Atlas 14 climatology
(~800 m PRISM interpolation). Smoothing the radar field is useful for two
distinct, legitimate reasons:

1. **Display / context** -- a lightly smoothed field reads against the smooth
   climatology without the 1-km speckle dominating the eye.
2. **Radar-gauge agreement** -- mild spatial averaging reduces the three errors
   that inflate radar-vs-gauge scatter: point-vs-pixel *representativeness*
   error, wind-drift / parallax *registration* error, and pixel-scale *noise*.
   The literature finds agreement improves up to an optimum (typically ~3-6 km;
   the classic approach is an N x N averaging window centred on the gauge) and
   degrades beyond it -- :func:`gauge_skill_sweep` finds that optimum.

**Caveat for this project.** ``i15`` is a *peak-intensity* metric (debris-flow
triggering). Smoothing mechanically lowers peaks and shrinks the radar's known
positive bias. So when judging "agreement", trust the **correlation** (up) and
**RMSE** (down); the bias *ratio* falls monotonically through 1.0 as a pure
side-effect of peak reduction, not as a skill gain.

Methods
-------
Four isotropic low-pass kernels, all **NaN-aware** so masked / no-coverage cells
never bleed zeros into their neighbours:

- ``gaussian`` -- Gaussian kernel (recommended default: smooth, fast, no blocky
  or rank artefacts);
- ``uniform``  -- boxcar mean over a square window (the literature's N x N
  averaging window);
- ``median``   -- edge-preserving median (de-speckles, keeps the convective edge
  sharper than a mean);
- ``idw``      -- inverse-distance weighting, kernel ``1/d**power`` inside a
  search radius (the moving-window analogue of point IDW interpolation).

The linear filters (gaussian / uniform / idw) use *normalised convolution* --
``smooth(data) / smooth(valid_mask)`` with zero-fill -- which simultaneously
prevents zero-bleed across NaN cells, renormalises at coverage edges, and (for
the Gaussian) corrects the truncated-kernel weight sum for free. The median is a
rank filter, so it instead uses :func:`scipy.ndimage.generic_filter` with a
NaN-ignoring callback over a disk footprint.

Smoothing scale
---------------
One physical knob, ``radius_km`` -- the *nominal smoothing scale* (about the
Gaussian sigma). It maps to a pixel extent per method so that a given
``radius_km`` is roughly visually comparable across methods:

- gaussian: ``sigma_px = radius_km / cell_km``;
- uniform / median: half-width ``h = round(sqrt(3) * radius_km / cell_km)`` then
  an **odd** window ``2h+1`` (the sqrt(3) makes the box's standard deviation equal
  ``radius_km``; the median shares the boxcar footprint -- comparable support, not
  comparable variance);
- idw: search radius ``r_px = round(2 * radius_km / cell_km)`` (~2-sigma support).

``radius_km = 0`` is the identity (no smoothing). Cell size is read from the
raster transform at the AOI-centre latitude (a 0.01 deg MRMS cell ~ 0.99 km).
"""

from __future__ import annotations

import math
import os
import warnings
from typing import Dict, Optional, Sequence

import numpy as np
import pandas as pd
import rasterio
from scipy import ndimage as ndi

# rainfall intensity / depth fields -> smoothable
DEFAULT_SMOOTH_FIELDS = ("i15max", "i30max", "i60max", "i2max", "total",
                         "peakrate_mmph", "mstotal")
# categorical / quality fields (angles, indices, beam heights) -> never smooth
NEVER_SMOOTH = ("tpki15", "rqi", "shsr", "cbb")

# method key -> display label (used by the figures + CLI choices)
METHODS: Dict[str, str] = {
    "gaussian": "Gaussian",
    "uniform": "Boxcar mean",
    "median": "Median",
    "idw": "Inverse-distance",
}

_KM_PER_DEG = 111.32
_SQRT3 = math.sqrt(3.0)


# --------------------------------------------------------------------------- #
# cell size
# --------------------------------------------------------------------------- #
def cell_size_km(transform, lat: float) -> float:
    """Mean ground cell size (km) of a **geographic** (degree) ``transform``.

    ``transform`` is a rasterio :class:`~affine.Affine`; ``lat`` is the latitude
    at which to evaluate the east-west cell width (cells shrink as ``cos(lat)``).
    Returns the mean of the east-west and north-south cell sizes.
    """
    dx = abs(transform.a) * _KM_PER_DEG * math.cos(math.radians(lat))
    dy = abs(transform.e) * _KM_PER_DEG
    return 0.5 * (dx + dy)


def _cell_km_from_dataset(ds) -> float:
    """Cell size (km) for an open rasterio dataset -- geographic or projected."""
    t = ds.transform
    if ds.crs is not None and ds.crs.is_geographic:
        lat = t.f + t.e * (ds.height / 2.0)          # AOI-centre latitude
        return cell_size_km(t, lat)
    return 0.5 * (abs(t.a) + abs(t.e)) / 1000.0       # projected: assume metres


# --------------------------------------------------------------------------- #
# NaN-aware kernels
# --------------------------------------------------------------------------- #
def _normalized(data: np.ndarray, lin_filter) -> np.ndarray:
    """NaN-aware linear smoothing via normalised convolution.

    ``lin_filter`` applies a linear filter to a zero-filled array (it must use
    ``mode="constant", cval=0`` so NaN / off-array cells contribute nothing).
    Returns ``lin_filter(data) / lin_filter(valid)``, NaN where no valid cell is
    in range.
    """
    valid = np.isfinite(data).astype(np.float64)
    filled = np.where(valid > 0, data, 0.0).astype(np.float64)
    num = lin_filter(filled)
    den = lin_filter(valid)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(den > 1e-12, num / den, np.nan)
    return out


def _gaussian(data: np.ndarray, sigma_px: float, truncate: float = 2.0) -> np.ndarray:
    return _normalized(data, lambda a: ndi.gaussian_filter(
        a, sigma=sigma_px, mode="constant", cval=0.0, truncate=truncate))


def _uniform(data: np.ndarray, half_px: int) -> np.ndarray:
    size = 2 * int(half_px) + 1                        # odd: no half-cell shift
    return _normalized(data, lambda a: ndi.uniform_filter(
        a, size=size, mode="constant", cval=0.0))


def _idw_kernel(r_px: float, power: float) -> np.ndarray:
    """``1/d**power`` weights inside radius ``r_px``; centre weight = 1.0."""
    r = max(1, int(np.ceil(r_px)))
    y, x = np.ogrid[-r:r + 1, -r:r + 1]
    d = np.sqrt(x * x + y * y)
    k = np.zeros_like(d, dtype=np.float64)
    inside = d <= r_px
    with np.errstate(divide="ignore"):
        k[inside] = 1.0 / np.power(d[inside], power)
    k[r, r] = 1.0                                       # centre: d=0 -> avoid inf
    return k


def _idw(data: np.ndarray, r_px: float, power: float = 2.0) -> np.ndarray:
    kernel = _idw_kernel(r_px, power)
    return _normalized(data, lambda a: ndi.correlate(
        a, kernel, mode="constant", cval=0.0))


def _disk(r: int) -> np.ndarray:
    y, x = np.ogrid[-r:r + 1, -r:r + 1]
    return (x * x + y * y) <= r * r + 1e-9             # odd-sided disk footprint


def _median(data: np.ndarray, half_px: int) -> np.ndarray:
    r = int(half_px)
    if r < 1:
        return data.astype(np.float64, copy=True)
    foot = _disk(r)

    def _nanmed(x):                                     # NaN-ignoring, all-NaN safe
        return np.nanmedian(x) if np.isfinite(x).any() else np.nan

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return ndi.generic_filter(data.astype(np.float64), _nanmed,
                                  footprint=foot, mode="constant", cval=np.nan)


# --------------------------------------------------------------------------- #
# dispatch
# --------------------------------------------------------------------------- #
def smooth_array(arr: np.ndarray, cell_km: float, method: str,
                 radius_km: float, power: float = 2.0) -> np.ndarray:
    """Smooth a 2-D field (NaN = missing) at ``radius_km`` with ``method``.

    ``cell_km`` is the raster cell size in km (see :func:`cell_size_km`).
    ``radius_km <= 0`` returns an unchanged copy. The kernel is clamped so it
    cannot exceed the array (heavily-averaged but still valid on tiny AOIs).
    """
    arr = np.asarray(arr, dtype=np.float64)
    if radius_km is None or radius_km <= 0:
        return arr.copy()
    if method not in METHODS:
        raise ValueError(f"unknown method '{method}'; choose from {list(METHODS)}")
    if not np.isfinite(arr).any():                      # all-NaN field -> passthrough
        return arr.copy()
    cell_km = float(cell_km)
    maxhalf = max(1, (min(arr.shape) - 1) // 2)
    if method == "gaussian":
        sigma_px = min(radius_km / cell_km, maxhalf / 2.0)
        return _gaussian(arr, sigma_px)
    if method == "uniform":
        half = max(1, min(int(round(_SQRT3 * radius_km / cell_km)), maxhalf))
        return _uniform(arr, half)
    if method == "median":
        half = max(1, min(int(round(_SQRT3 * radius_km / cell_km)), maxhalf))
        return _median(arr, half)
    # idw
    r_px = min(2.0 * radius_km / cell_km, float(maxhalf))
    return _idw(arr, r_px, power=power)


# --------------------------------------------------------------------------- #
# raster I/O
# --------------------------------------------------------------------------- #
def _read_field(path: str):
    """Read a single-band raster -> (float64 array with NaN, transform, crs,
    cell_km, profile)."""
    with rasterio.open(path) as ds:
        arr = ds.read(1).astype(np.float64)
        nodata = ds.nodata
        transform, crs = ds.transform, ds.crs
        cell_km = _cell_km_from_dataset(ds)
        profile = ds.profile.copy()
    if nodata is not None and np.isfinite(nodata):
        arr = np.where(arr == nodata, np.nan, arr)
    return arr, transform, crs, cell_km, profile


def smooth_tif(in_path: str, out_path: str, method: str, radius_km: float,
               power: float = 2.0) -> str:
    """Smooth one field GeoTIFF -> ``out_path`` (same grid/CRS, float32, NaN nodata)."""
    arr, _, _, cell_km, profile = _read_field(in_path)
    out = smooth_array(arr, cell_km, method, radius_km, power=power).astype(np.float32)
    profile.update(dtype="float32", nodata=np.nan, count=1)
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    with rasterio.open(out_path, "w", **profile) as d:
        d.write(out, 1)
    return out_path


def smooth_event_fields(in_dir: str, key: str, out_dir: str, out_key: str,
                        method: str, radius_km: float,
                        fields: Sequence[str] = DEFAULT_SMOOTH_FIELDS,
                        power: float = 2.0):
    """Smooth every intensity field GeoTIFF of an event; skip missing ones.

    Reads ``<in_dir>/<key>_<field>.tif`` for each in ``fields`` and writes
    ``<out_dir>/<out_key>_<field>.tif``. Categorical / quality fields
    (:data:`NEVER_SMOOTH`) are never touched. Returns the written paths.
    """
    os.makedirs(out_dir, exist_ok=True)
    written = []
    for f in fields:
        if f in NEVER_SMOOTH:
            continue
        ip = os.path.join(in_dir, f"{key}_{f}.tif")
        if not os.path.exists(ip):
            continue
        op = os.path.join(out_dir, f"{out_key}_{f}.tif")
        smooth_tif(ip, op, method, radius_km, power=power)
        written.append(op)
    return written


def smooth_dataarray(src, method: str, radius_km: float, power: float = 2.0):
    """Smooth a single-band raster (path or rioxarray ``DataArray``) in memory.

    Returns a 2-D ``DataArray`` smoothed with :func:`smooth_array` (NaN-aware),
    keeping the source coords / CRS / transform -- so it can flow straight into
    a reproject/plot or into :func:`stormscape.atlas14.anomaly` (which accepts a
    ``DataArray``) with no intermediate GeoTIFF. The cell size is read from the
    geotransform (geographic grids use the centre-latitude east-west size). A
    falsy ``method`` or ``radius_km <= 0`` returns the (loaded) field unchanged.
    """
    import rioxarray  # noqa: F401  -- registers the .rio accessor; lazy import
    if isinstance(src, str):
        with rioxarray.open_rasterio(src, masked=True) as d:   # close promptly
            da = d.load()
    else:
        da = src
    if not method or radius_km is None or radius_km <= 0:
        return da
    # squeeze only to read the 2-D grid + transform; keep ``da``'s own dim order
    # (e.g. band,y,x from a path) so the result writes back / divides cleanly.
    arr2d = da.squeeze()
    vals = np.asarray(arr2d.values, dtype=np.float64)
    t = arr2d.rio.transform()
    crs = arr2d.rio.crs
    if crs is not None and crs.is_geographic:
        lat = t.f + t.e * (vals.shape[0] / 2.0)               # AOI-centre lat
        cell_km = cell_size_km(t, lat)
    else:
        cell_km = 0.5 * (abs(t.a) + abs(t.e)) / 1000.0        # projected: metres
    sm = smooth_array(vals, cell_km, method, radius_km, power=power)
    return da.copy(data=sm.reshape(da.shape).astype(da.dtype))


# --------------------------------------------------------------------------- #
# gauge-skill sweep
# --------------------------------------------------------------------------- #
def _sample_array_at_points(arr: np.ndarray, transform, crs,
                            pts_gdf) -> np.ndarray:
    """Nearest-cell sample of ``arr`` at point geometries.

    Matches :func:`stormscape.compare.sample_raster_at_points` semantics: points
    are reprojected into ``crs`` and indexed with ``rasterio.transform.rowcol``
    (default ``op=floor`` -- the same indexing ``ds.sample`` uses; ``round`` would
    be off by half a cell). Off-grid / NaN cells return NaN.
    """
    pts = pts_gdf.to_crs(crs) if pts_gdf.crs is not None else pts_gdf
    xs = np.asarray([g.x for g in pts.geometry], dtype=np.float64)
    ys = np.asarray([g.y for g in pts.geometry], dtype=np.float64)
    rows, cols = rasterio.transform.rowcol(transform, xs, ys)   # op=floor
    rows, cols = np.asarray(rows).astype(int), np.asarray(cols).astype(int)
    out = np.full(rows.shape, np.nan, dtype=np.float64)
    ok = (rows >= 0) & (rows < arr.shape[0]) & (cols >= 0) & (cols < arr.shape[1])
    out[ok] = arr[rows[ok], cols[ok]]
    return out


def gauge_skill_sweep(gauges, in_dir: str, key: str,
                      methods: Sequence[str] = tuple(METHODS),
                      radii_km: Sequence[float] = (0, 0.5, 1, 2, 3, 4, 6, 8),
                      durations: Sequence[int] = (15, 30, 60),
                      rqi_min: Optional[float] = None,
                      max_report_min: Optional[float] = None,
                      power: float = 2.0) -> pd.DataFrame:
    """Sweep smoothing radius vs radar-gauge skill, per method and duration.

    For each duration ``d`` the field ``<in_dir>/<key>_i{d}max.tif`` is loaded
    once and the gauge cell indices precomputed once; then for every
    ``(method, radius)`` the field is smoothed in memory, sampled at the gauges,
    and scored with :func:`stormscape.compare.comparison_stats` (so the RQI and
    sub-hourly cadence screens are identical to production ``compare``). The
    ``radius_km = 0`` row is the raw (unsmoothed) baseline.

    Returns a tidy DataFrame with columns ``method, radius_km, duration, metric,
    n, gauge_mean, radar_mean, bias, rmse, mae, corr, ratio``.
    """
    from .compare import comparison_stats, sample_raster_at_points

    gcol_for = {d: f"i{d}_mmph" for d in durations}
    have_report = "report_min" in getattr(gauges, "columns", [])

    rqi_vals = None
    if rqi_min is not None:
        rqi_path = os.path.join(in_dir, f"{key}_rqi.tif")
        if os.path.exists(rqi_path):
            rqi_vals = sample_raster_at_points(gauges, rqi_path)

    out_rows = []
    for d in durations:
        gcol = gcol_for[d]
        if gcol not in getattr(gauges, "columns", []):
            continue
        path = os.path.join(in_dir, f"{key}_i{d}max.tif")
        if not os.path.exists(path):
            continue
        arr, transform, crs, cell_km, _ = _read_field(path)
        pts = gauges.to_crs(crs) if gauges.crs is not None else gauges
        xs = np.asarray([g.x for g in pts.geometry], dtype=np.float64)
        ys = np.asarray([g.y for g in pts.geometry], dtype=np.float64)
        rows, cols = rasterio.transform.rowcol(transform, xs, ys)   # op=floor
        rows, cols = np.asarray(rows).astype(int), np.asarray(cols).astype(int)
        ok = (rows >= 0) & (rows < arr.shape[0]) & (cols >= 0) & (cols < arr.shape[1])
        gauge_vals = pd.to_numeric(gauges[gcol], errors="coerce").to_numpy()

        for method in methods:
            for r in radii_km:
                sm = smooth_array(arr, cell_km, method, r, power=power)
                radar = np.full(rows.shape, np.nan, dtype=np.float64)
                radar[ok] = sm[rows[ok], cols[ok]]
                tbl = pd.DataFrame({
                    gcol: gauge_vals,
                    f"radar_i{d}max": radar,
                    f"resid_i{d}max": radar - gauge_vals,
                })
                if rqi_vals is not None:
                    tbl["rqi"] = rqi_vals
                if have_report:
                    tbl["report_min"] = pd.to_numeric(
                        gauges["report_min"], errors="coerce").to_numpy()
                st = comparison_stats(tbl, pairs={gcol: f"i{d}max"},
                                      rqi_min=rqi_min,
                                      max_report_min=max_report_min)
                row = st.iloc[0].to_dict()
                row.update(method=method, radius_km=float(r), duration=int(d))
                out_rows.append(row)

    cols_order = ["method", "radius_km", "duration", "metric", "n",
                  "gauge_mean", "radar_mean", "bias", "rmse", "mae", "corr",
                  "ratio"]
    df = pd.DataFrame(out_rows)
    if len(df):
        df = df[[c for c in cols_order if c in df.columns]]
    return df


def best_radius(sweep: pd.DataFrame, by: str = "rmse",
                duration: int = 15) -> Dict[str, tuple]:
    """Best (radius_km, score) per method for one duration, excluding raw (r=0).

    ``by`` is minimised for ``rmse``/``mae``/``bias`` and maximised for ``corr``.
    """
    sub = sweep[(sweep["duration"] == duration) & (sweep["radius_km"] > 0)]
    minimise = by in ("rmse", "mae")
    out: Dict[str, tuple] = {}
    for m, g in sub.groupby("method"):
        g = g.dropna(subset=[by])
        if not len(g):
            out[m] = (float("nan"), float("nan"))
            continue
        idx = g[by].idxmin() if minimise else g[by].idxmax()
        out[m] = (float(g.loc[idx, "radius_km"]), float(g.loc[idx, by]))
    return out
