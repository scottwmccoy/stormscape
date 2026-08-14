"""Radar QPE vs ground-gauge comparison.

Samples the MRMS radar fields (storm ``total`` and the peak ``i15/i30/i60``
intensities) at rain-gauge locations and computes, per gauge, the residual
``radar - gauge`` plus aggregate skill statistics (bias, RMSE, MAE,
correlation, mass ratio). Gauges far from the radar or behind terrain are
unreliable -- the beam overshoots low rain -- so samples can be screened by the
radar quality index (RQI) at the gauge cell.

Pairs the per-gauge metrics from :func:`stormscape.gauges.gauge_fields` with
the radar GeoTIFFs written by :func:`stormscape.mrms.save_fields`. All sampling
reprojects gauges into each raster's CRS, so the native-EPSG:4326 MRMS rasters
and EPSG:4326 gauge points line up without manual reprojection.
"""

from __future__ import annotations

import os
from typing import Dict, Optional

import geopandas as gpd
import numpy as np
import pandas as pd

from .layout import find
import rasterio

# gauge metric column -> radar field (GeoTIFF stem ``<key>_<field>.tif``)
DEFAULT_PAIRS: Dict[str, str] = {
    "total_mm": "total",
    "i15_mmph": "i15max",
    "i30_mmph": "i30max",
    "i60_mmph": "i60max",
}


def sample_raster_at_points(points: gpd.GeoDataFrame,
                            raster_path: str) -> np.ndarray:
    """Sample a single-band raster at point geometries -> float ndarray.

    Points are reprojected into the raster's CRS first; off-raster cells and
    the raster's nodata value come back as NaN.
    """
    with rasterio.open(raster_path) as ds:
        pts = points.to_crs(ds.crs)
        coords = [(geom.x, geom.y) for geom in pts.geometry]
        vals = np.array([v[0] for v in ds.sample(coords)], dtype="float64")
        nodata = ds.nodata
    if nodata is not None and np.isfinite(nodata):
        vals = np.where(vals == nodata, np.nan, vals)
    return vals


def radar_vs_gauge(gauges: gpd.GeoDataFrame, rasters: Dict[str, str],
                   pairs: Dict[str, str] = DEFAULT_PAIRS,
                   rqi_path: Optional[str] = None,
                   rqi_min: Optional[float] = None) -> gpd.GeoDataFrame:
    """Per-gauge radar-vs-gauge table.

    Parameters
    ----------
    gauges
        Gauge GeoDataFrame with the metric columns named in ``pairs`` keys
        (e.g. from :func:`stormscape.gauges.gauge_fields`).
    rasters
        ``{radar_field: path}`` for the radar GeoTIFFs to sample.
    pairs
        ``{gauge_column: radar_field}`` mapping (default :data:`DEFAULT_PAIRS`).
    rqi_path, rqi_min
        Optional RQI raster to sample at each gauge; ``rqi_min`` adds a boolean
        ``rqi_ok`` column (it does not drop rows -- screening happens in
        :func:`comparison_stats`).

    Returns
    -------
    GeoDataFrame
        ``gauges`` plus, for each pair, ``radar_<field>`` and
        ``resid_<field>`` (= radar - gauge), and optionally ``rqi``/``rqi_ok``.
    """
    out = gauges.copy()
    for gcol, rfield in pairs.items():
        path = rasters.get(rfield)
        if not path or gcol not in out.columns:
            continue
        out[f"radar_{rfield}"] = sample_raster_at_points(out, path)
        out[f"resid_{rfield}"] = out[f"radar_{rfield}"] - pd.to_numeric(
            out[gcol], errors="coerce")
    if rqi_path:
        out["rqi"] = sample_raster_at_points(out, rqi_path)
        if rqi_min is not None:
            out["rqi_ok"] = out["rqi"] >= rqi_min
    return out


def comparison_stats(table: gpd.GeoDataFrame,
                     pairs: Dict[str, str] = DEFAULT_PAIRS,
                     rqi_min: Optional[float] = None,
                     max_report_min: Optional[float] = None) -> pd.DataFrame:
    """Aggregate skill statistics per metric from a :func:`radar_vs_gauge` table.

    One row per metric with ``n`` (valid pairs), gauge/radar means, ``bias``
    (mean radar - gauge), ``rmse``, ``mae``, Pearson ``corr``, and ``ratio``
    (Σradar / Σgauge). When ``rqi_min`` is given and an ``rqi`` column exists,
    only gauges with ``rqi >= rqi_min`` are counted. When ``max_report_min`` is
    given and a ``report_min`` column exists, the **sub-hourly** metrics
    (everything except ``total``) additionally drop gauges whose native
    reporting interval exceeds it -- coarse reporters smear short-duration
    peaks, so this keeps the i15/i30/i60 comparison fair while leaving the
    cadence-insensitive storm total on the full sample.
    """
    rqi_mask = None
    if rqi_min is not None and "rqi" in table.columns:
        rqi_mask = table["rqi"] >= rqi_min
    cadence = (pd.to_numeric(table["report_min"], errors="coerce")
               if max_report_min is not None and "report_min" in table.columns
               else None)
    rows = []
    for gcol, rfield in pairs.items():
        rcol, dcol = f"radar_{rfield}", f"resid_{rfield}"
        if dcol not in table.columns:
            continue
        g = pd.to_numeric(table[gcol], errors="coerce")
        r = pd.to_numeric(table[rcol], errors="coerce")
        m = g.notna() & r.notna()
        if rqi_mask is not None:
            m &= rqi_mask
        if cadence is not None and rfield != "total":   # cadence screen
            m &= cadence <= max_report_min
        gg, rr = g[m].to_numpy(), r[m].to_numpy()
        n = int(m.sum())
        if n == 0:
            rows.append(dict(metric=rfield, n=0, gauge_mean=np.nan,
                             radar_mean=np.nan, bias=np.nan, rmse=np.nan,
                             mae=np.nan, corr=np.nan, ratio=np.nan))
            continue
        resid = rr - gg
        corr = (float(np.corrcoef(gg, rr)[0, 1])
                if n > 1 and gg.std() > 0 and rr.std() > 0 else np.nan)
        rows.append(dict(
            metric=rfield, n=n,
            gauge_mean=float(gg.mean()), radar_mean=float(rr.mean()),
            bias=float(resid.mean()),
            rmse=float(np.sqrt((resid ** 2).mean())),
            mae=float(np.abs(resid).mean()), corr=corr,
            ratio=float(rr.sum() / gg.sum()) if gg.sum() else np.nan))
    return pd.DataFrame(rows)


def compare_storm(gauges: gpd.GeoDataFrame, out_dir: str, key: str,
                  pairs: Dict[str, str] = DEFAULT_PAIRS,
                  rqi_min: Optional[float] = None,
                  max_report_min: Optional[float] = None,
                  multisensor: bool = False):
    """Compare gauges against the radar GeoTIFFs saved under ``out_dir/<key>_*``.

    Convenience wrapper that resolves the ``<key>_<field>.tif`` paths written
    by :func:`stormscape.mrms.save_fields` (plus ``<key>_rqi.tif``) and runs
    :func:`radar_vs_gauge` + :func:`comparison_stats`.

    With ``multisensor=True`` and a ``<key>_mstotal.tif`` present (gauge-corrected
    MRMS total from :func:`stormscape.mrms.multisensor_total`), an extra
    ``mstotal`` row compares the gauge storm total against it -- the
    gauge-corrected counterpart to the radar-only ``total`` row.

    Returns ``(table, stats)``.
    """
    rasters = {}
    for rfield in set(pairs.values()):
        p = find(out_dir, f"{key}_{rfield}.tif")
        if os.path.exists(p):
            rasters[rfield] = p
    rqi_p = find(out_dir, f"{key}_rqi.tif")
    rqi_path = rqi_p if os.path.exists(rqi_p) else None
    table = radar_vs_gauge(gauges, rasters, pairs=pairs, rqi_path=rqi_path,
                           rqi_min=rqi_min)
    ms_path = find(out_dir, f"{key}_mstotal.tif")
    if multisensor and os.path.exists(ms_path) and "total_mm" in table.columns:
        table["radar_mstotal"] = sample_raster_at_points(table, ms_path)
        table["resid_mstotal"] = table["radar_mstotal"] - pd.to_numeric(
            table["total_mm"], errors="coerce")
    stats = comparison_stats(table, pairs=pairs, rqi_min=rqi_min,
                             max_report_min=max_report_min)
    if "resid_mstotal" in table.columns:                 # gauge-corrected total
        ms = comparison_stats(table, pairs={"total_mm": "mstotal"},
                              rqi_min=rqi_min)
        stats = pd.concat([stats, ms], ignore_index=True)
    return table, stats


def gauge_recurrence_table(gauges: gpd.GeoDataFrame,
                           durations=(15, 30, 60), stat: str = "mean",
                           series: str = "pds",
                           peak_times: Optional[Dict[str, object]] = None,
                           aoi_bounds=None, name_col: str = "name") -> pd.DataFrame:
    """Per-(wet)-gauge recurrence table vs the NOAA Atlas 14 climatology.

    For every **wet** gauge (peak ``i15_mmph`` > 0) in ``gauges`` (an EPSG:4326
    point GeoDataFrame with ``i{d}_mmph`` peak columns, e.g. from
    :func:`stormscape.gauges.gauge_fields`), records the observed peak ``I{d}``,
    the **anomaly** (observed / 1-yr climatology) and the **recurrence interval**
    of each observed peak. The climatology comes from the NOAA PFDS **point**
    service per gauge (:func:`stormscape.atlas14.pf_point`) -- the full ARI curve,
    so the anomaly's 1-yr reference and the recurrence interval share one
    point-accurate source. ``peak_times`` is an optional ``{gauge_name:
    timestamp}`` map (e.g. the time of the I15 peak from the gauge series) added
    as an ``i15_peak_time`` column. Returns a tidy DataFrame sorted by I15.

    Columns: ``gauge, station_id, lon, lat, report_min`` (native gauge reporting
    cadence, min -- coarse reporters smear i15/i30 low, so their short-duration RIs
    read low) ``[, i15_peak_time]`` then per duration ``i{d}`` (observed mm/h),
    ``clim1yr_i{d}``, ``anom_i{d}``, ``RI_i{d}`` (years; ``nan`` = below the 1-yr
    quantile, ``inf`` = above the top tabulated ARI).
    """
    from . import atlas14
    g = gauges.to_crs(4326) if getattr(gauges, "crs", None) is not None else gauges
    i15 = pd.to_numeric(g.get("i15_mmph"), errors="coerce")
    wet = g[i15.fillna(0) > 0].copy()
    peak_times = peak_times or {}
    has_geo_tpk = "i15_peak_time" in wet.columns

    def _dist_aoi_km(lon, lat):                          # 0 inside, else km to bbox
        if not aoi_bounds:
            return None
        w, s, e, n = aoi_bounds
        dx = max(w - lon, 0.0, lon - e) * 111.32 * np.cos(np.radians(lat))
        dy = max(s - lat, 0.0, lat - n) * 111.32
        return float((dx * dx + dy * dy) ** 0.5)

    out = []
    for _, r in wet.iterrows():
        name = r.get(name_col)
        lon, lat = float(r.geometry.x), float(r.geometry.y)
        try:                                             # native gauge tempo (min)
            rm = float(r.get("report_min"))
            if not np.isfinite(rm):
                rm = None
        except (TypeError, ValueError):
            rm = None
        rec = {"gauge": name, "station_id": r.get("station_id"),
               "lon": lon, "lat": lat, "report_min": rm}
        # I15 time-of-peak: prefer the geojson column (from fetch_gauge_event),
        # else the supplied peak_times map (e.g. read from the series CSVs)
        tpk = r.get("i15_peak_time") if has_geo_tpk else peak_times.get(name)
        if tpk is not None or peak_times or has_geo_tpk:
            rec["i15_peak_time"] = tpk
        if aoi_bounds:
            d = _dist_aoi_km(lon, lat)
            rec["dist_to_aoi_km"] = round(d, 2)
            rec["in_aoi"] = d <= 0.0
        try:
            pf = atlas14.pf_point(rec["lat"], rec["lon"], durations=durations,
                                  stat=stat, series=series)
        except Exception:                                # noqa: BLE001 (network)
            pf = None
        for d in durations:
            obs = float(pd.to_numeric(r.get(f"i{d}_mmph"), errors="coerce"))
            rec[f"i{d}"] = obs
            if pf is not None and f"i{d}" in pf.columns:
                col = pf[f"i{d}"]
                c1 = float(col.loc[1.0] if 1.0 in col.index else col.iloc[0])
                rec[f"clim1yr_i{d}"] = c1
                rec[f"anom_i{d}"] = obs / c1 if c1 > 0 else float("nan")
                rec[f"RI_i{d}"] = atlas14.recurrence_interval(
                    obs, pf.index.values, col.values)
        out.append(rec)
    df = pd.DataFrame(out)
    if "i15" in df.columns:
        df = df.sort_values("i15", ascending=False).reset_index(drop=True)
    return df
