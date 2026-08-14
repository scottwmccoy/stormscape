"""Ground rain-gauge rainfall (Synoptic / MesoWest) for an AOI + storm window.

Point-gauge counterpart to :mod:`stormscape.mrms`. Where that engine stacks
radar returns into an i15 *field*, this pulls rain-gauge observations from the
**Synoptic Data / MesoWest** API over the same AOI and storm window, and
reduces each gauge to the quantities used for the radar-vs-gauge comparison:
storm **total** and peak **15 / 30 / 60-minute** intensity (mm/h).

The Synoptic transport (:func:`get_stations`, :func:`get_rainfall`) and the
irregular-observation -> 1-minute interpolation
(:func:`precipitation_per_minute`) are adapted from the USGS **FlowAlert**
package (King, Rengers, Wedell & Fee, 2024; CC0 / U.S. public domain --
https://code.usgs.gov/ghsc/lhp/FlowAlert). The peak-intensity reducer re-uses
stormscape's own ``(i16 + i14) / 2`` estimator (see :mod:`stormscape.mrms`) so
the gauge i15 is defined identically to the radar ``i15max``; the 30/60-min
peaks use plain trailing windows.

A Synoptic API token is required (free for academic / research use): pass it
explicitly as ``token=`` or set ``$SYNOPTIC_TOKEN``.
"""

from __future__ import annotations

import os
import warnings
from datetime import datetime, timezone
from typing import Dict, List, Optional, Sequence, Tuple

import geopandas as gpd
import numpy as np
import pandas as pd
import requests

from .aoi import load_aoi
from .layout import out_path

warnings.filterwarnings("ignore")

SYNOPTIC_BASE = "https://api.synopticdata.com/v2"
META_URL = f"{SYNOPTIC_BASE}/stations/metadata"
TS_URL = f"{SYNOPTIC_BASE}/stations/timeseries"

DEFAULT_DURATIONS = (15, 30, 60)

# Synoptic derived-precip variables requested with precip=1 (units mm).
ACC_VAR = "precip_accumulated_set_1d"     # cumulative precip over the request
INT_VAR = "precip_intervals_set_1d"       # precip since the previous observation


class SynopticError(Exception):
    """Synoptic API returned an error (bad token, account access, etc.)."""


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _token(token: Optional[str]) -> str:
    tok = (token or os.environ.get("SYNOPTIC_TOKEN")
           or os.environ.get("MESOWEST_TOKEN"))
    if not tok:
        raise SynopticError(
            "No Synoptic token. Pass token=... or set $SYNOPTIC_TOKEN "
            "(free academic/research access at synopticdata.com).")
    return tok


def _utc(t: datetime) -> datetime:
    """Coerce a datetime to UTC (naive datetimes are assumed already UTC)."""
    return t.replace(tzinfo=timezone.utc) if t.tzinfo is None \
        else t.astimezone(timezone.utc)


def _raise_for_synoptic(r: requests.Response) -> None:
    """Surface Synoptic HTTP errors with their human-readable message.

    A revoked token or an account without an active data plan returns HTTP 403
    ``"Unauthorized. Review your account access settings in the customer
    console."`` -- propagate that text rather than a downstream KeyError.
    """
    if r.status_code == 200:
        return
    try:
        body = r.json()
        msg = body.get("message") or body.get("ERROR") or str(body)
    except Exception:                                  # noqa: BLE001
        msg = r.text[:200]
    raise SynopticError(f"HTTP {r.status_code}: {msg}")


def _stations_json_to_gdf(payload: dict, status: Optional[str]) -> gpd.GeoDataFrame:
    """Parse a Synoptic stations/metadata payload -> GeoDataFrame(EPSG:4326)."""
    cols = ["station_id", "name", "state", "elevation", "network", "provider",
            "timezone", "status", "geometry"]
    code = payload.get("SUMMARY", {}).get("RESPONSE_CODE")
    if code == 2 or "STATION" not in payload:          # 2 = zero results
        return gpd.GeoDataFrame(columns=cols, geometry="geometry", crs=4326)
    df = pd.DataFrame(payload["STATION"])
    lon = pd.to_numeric(df["LONGITUDE"], errors="coerce")
    lat = pd.to_numeric(df["LATITUDE"], errors="coerce")
    out = gpd.GeoDataFrame(
        {
            "station_id": df["STID"].astype(str),
            "name": df.get("NAME"),
            "state": df.get("STATE"),
            "elevation": pd.to_numeric(df.get("ELEVATION"), errors="coerce"),
            "network": pd.to_numeric(df.get("MNET_ID"), errors="coerce"),
            "provider": df.get("SHORTNAME"),
            "timezone": df.get("TIMEZONE"),
            "status": df.get("STATUS"),
        },
        geometry=gpd.points_from_xy(lon, lat), crs=4326)
    if status:
        out = out[out["status"].fillna("").str.upper() == status.upper()]
    return out.reset_index(drop=True)


def _timeseries_json_to_frames(payload: dict) -> Dict[str, pd.DataFrame]:
    """Parse a Synoptic stations/timeseries payload -> {stid: DataFrame}."""
    out: Dict[str, pd.DataFrame] = {}
    for st in payload.get("STATION", []):
        d = pd.DataFrame(st.get("OBSERVATIONS", {}))
        if "date_time" in d:
            d["date_time"] = pd.to_datetime(d["date_time"])
        out[str(st["STID"])] = d
    return out


# --------------------------------------------------------------------------- #
# Synoptic transport (adapted from FlowAlert, CC0)
# --------------------------------------------------------------------------- #
def get_stations(aoi, token: Optional[str] = None, status: Optional[str] = "ACTIVE",
                 pad_deg: float = 0.0) -> gpd.GeoDataFrame:
    """Synoptic stations within the AOI bbox -> GeoDataFrame(EPSG:4326).

    One row per station with ``station_id, name, state, elevation`` (m),
    ``network`` (Synoptic MNET_ID), ``provider, timezone, status`` and a Point
    ``geometry`` (lon/lat). ``status='ACTIVE'`` drops decommissioned sites;
    pass ``status=None`` to keep all. Adapted from FlowAlert's
    ``mesowest.get_stations``; returns a GeoDataFrame (not a plain DataFrame)
    to match :mod:`stormscape.refdata` conventions.
    """
    bounds, _ = load_aoi(aoi, pad_deg=pad_deg)
    bbox = ",".join(f"{v:.6f}" for v in bounds)
    r = requests.get(META_URL, timeout=60,
                     params={"token": _token(token), "bbox": bbox, "complete": "1"})
    _raise_for_synoptic(r)
    return _stations_json_to_gdf(r.json(), status)


def get_rainfall(station_ids: Sequence[str], starttime: datetime,
                 endtime: datetime, token: Optional[str] = None,
                 url: str = TS_URL) -> Dict[str, pd.DataFrame]:
    """Per-station precip time series from Synoptic over ``[start, end]`` (UTC).

    Returns ``{station_id: DataFrame}`` whose columns include ``date_time``
    (UTC), :data:`ACC_VAR` (cumulative mm) and :data:`INT_VAR` (mm since the
    previous observation). Over-large requests (Synoptic "Querying too many
    station hours") are split in half and recombined. Adapted from FlowAlert's
    ``mesowest.get_rainfall``.
    """
    tok = _token(token)
    fmt = "%Y%m%d%H%M"
    params = {
        "start": _utc(starttime).strftime(fmt),
        "end": _utc(endtime).strftime(fmt),
        "stid": ",".join(station_ids),
        "obtimezone": "UTC", "output": "json", "precip": "1",
        "qc": "on", "qc_flags": "on", "qc_remove_data": "on",
        "units": "metric,precip|mm", "token": tok,
    }
    r = requests.get(url, params=params, timeout=120)
    _raise_for_synoptic(r)
    payload = r.json()
    if "STATION" in payload:
        return _timeseries_json_to_frames(payload)

    msg = payload.get("SUMMARY", {}).get("RESPONSE_MESSAGE", "")
    if "Querying too many station hours" in msg:       # split and recurse
        mid = starttime + (endtime - starttime) / 2
        first = get_rainfall(station_ids, starttime, mid, token=tok, url=url)
        second = get_rainfall(station_ids, mid, endtime, token=tok, url=url)
        for sid, d in second.items():
            if sid in first and len(first[sid]):
                first[sid] = (pd.concat([first[sid], d])
                              .drop_duplicates(subset="date_time")
                              .sort_values("date_time").reset_index(drop=True))
            else:
                first[sid] = d
        return first
    raise SynopticError(msg or "empty Synoptic response")


# --------------------------------------------------------------------------- #
# gauge -> per-minute -> peak intensity  (intensity reducer)
# --------------------------------------------------------------------------- #
def precipitation_per_minute(
        data: pd.DataFrame, precip_var: str = INT_VAR,
        time_var: str = "date_time") -> Tuple[list, np.ndarray, np.ndarray]:
    """Irregular gauge obs -> 1-minute cumulative + per-minute precip (mm).

    Interpolates the cumulative gauge record onto a regular 1-minute grid, so
    tipping-bucket records reported at uneven intervals become a clean
    minute-resolution series the rolling-window estimators can run on. Adapted
    from FlowAlert's ``rainfall.precipitation_per_minute``.

    Returns ``(times, cumulative, per_minute)``; the first ``per_minute`` value
    is 0 by construction. Empty inputs (< 2 valid obs) return empties.
    """
    if precip_var not in data.columns or time_var not in data.columns:
        return [], np.array([]), np.array([])
    data = data.dropna(subset=[precip_var])
    if len(data) < 2:
        return [], np.array([]), np.array([])
    precip = data[precip_var].clip(lower=0).to_numpy(dtype="float64")
    times = np.array([pd.Timestamp(t).timestamp() for t in data[time_var]])
    cumulative = np.cumsum(precip)
    grid = np.arange(times[0], times[-1] + 60, 60)
    interp = np.interp(grid, times, cumulative)
    per_minute = np.concatenate(([0.0], np.diff(interp)))
    out_times = [datetime.fromtimestamp(t, tz=timezone.utc) for t in grid]
    return out_times, interp, per_minute


def _peak_window(per_minute: np.ndarray, minutes: int) -> float:
    """Peak trailing ``minutes``-window intensity (mm/h) over a 1-min series."""
    if len(per_minute) < minutes:
        return float("nan")
    sums = np.convolve(per_minute, np.ones(minutes), mode="valid")
    return float(np.nanmax(sums) * 60.0 / minutes)


def _peak_i15(per_minute: np.ndarray) -> float:
    """Peak 15-min intensity via stormscape's ``(i16 + i14) / 2`` estimator.

    Mirrors :func:`stormscape.mrms.compute_i15` on a 1-minute series: the mean
    of the trailing 16-minute and trailing 14-minute intensities, end-aligned
    to the same final minute, maximised over the storm.
    """
    n = len(per_minute)
    if n < 16:
        return _peak_window(per_minute, 15)
    s16 = np.convolve(per_minute, np.ones(16), mode="valid")   # len n-15
    s14 = np.convolve(per_minute, np.ones(14), mode="valid")   # len n-13
    # a 16-min sum at index j spans minutes [j, j+15]; the 14-min sum ending on
    # the same minute (j+15) is s14[j+2], spanning [j+2, j+15].
    i16 = s16 * 60.0 / 16.0
    i14 = s14[2:2 + len(s16)] * 60.0 / 14.0
    return float(np.nanmax((i16 + i14) / 2.0))


def _report_min(obs) -> float:
    """Native precip **cadence** (minutes): the median interval between
    *precip-bearing* observations; ``NaN`` for <2 such obs / no precip variable.

    Measured on the rows that actually carry ``precip_intervals`` -- so a
    multi-variable station (an ASOS reporting 5-min temperature but precip only at
    hourly METAR times) is timed by its hourly *precip* interval, not its 5-min obs
    cadence. Coarse (hourly/daily) reporters have their bursts smeared by the
    1-minute interpolation, biasing the short-duration peaks (i15/i30) low; the
    storm total is cadence-insensitive. Surfaced alongside the metrics (the
    ``report_min`` gauge field) so a comparison can screen on it
    (``--max-report-min``) and the figures can report each gauge's tempo.
    """
    if (obs is None or not hasattr(obs, "columns") or INT_VAR not in obs.columns
            or "date_time" not in obs.columns):
        return float("nan")
    vt = pd.to_datetime(obs.dropna(subset=[INT_VAR])["date_time"]).sort_values()
    if len(vt) < 2:
        return float("nan")
    return float(vt.diff().dropna().dt.total_seconds().median() / 60.0)


def gauge_intensities(obs: pd.DataFrame,
                      durations: Sequence[int] = DEFAULT_DURATIONS) -> dict:
    """Reduce one gauge's observations to total + peak intensities.

    Returns ``{'total_mm', 'n_obs', 'n_min', 'report_min', 'i15_mmph',
    'i30_mmph', 'i60_mmph', ...}`` where ``n_obs``/``report_min`` count and time
    the **precip-bearing** observations (the native precip cadence; for a
    multi-variable ASOS this is the hourly METAR precip interval, not the 5-min
    obs cadence). The 15-minute peak uses the radar-matching ``(i16 + i14) / 2``
    estimator; the others use plain trailing windows. Stations reporting no
    precip variable yield all-NaN metrics.
    """
    nan_out = {"total_mm": float("nan"), "n_obs": 0, "n_min": 0,
               "report_min": float("nan"),
               **{f"i{d}_mmph": float("nan") for d in durations}}
    if obs is None or not hasattr(obs, "columns") or INT_VAR not in obs.columns:
        return nan_out                       # station reported no precip variable
    # Use only precip-bearing obs: with precip=1 a multi-variable station (an
    # ASOS, say) reports many sub-hourly rows but precip_intervals only at METAR
    # times, so cadence and total must come from the rows that carry precip --
    # the same rows precipitation_per_minute interpolates over.
    valid = obs.dropna(subset=[INT_VAR])
    total = float(pd.to_numeric(valid[INT_VAR], errors="coerce")
                  .clip(lower=0).sum())
    report_min = _report_min(obs)        # native precip cadence (see _report_min)
    _, _, per_minute = precipitation_per_minute(obs)
    out = {"total_mm": total, "n_obs": int(len(valid)),
           "n_min": int(len(per_minute)), "report_min": report_min}
    for d in durations:
        out[f"i{d}_mmph"] = (_peak_i15(per_minute) if d == 15
                             else _peak_window(per_minute, d))
    return out


def gauge_fields(aoi, starttime: datetime, endtime: datetime,
                 token: Optional[str] = None,
                 durations: Sequence[int] = DEFAULT_DURATIONS,
                 status: Optional[str] = "ACTIVE", batch: int = 50,
                 pad_deg: float = 0.05) -> gpd.GeoDataFrame:
    """Gauges in the AOI with their storm total + peak intensities.

    Discovers stations in the AOI, downloads their precip over
    ``[starttime, endtime]`` (UTC; use the radar storm window from
    :func:`stormscape.mrms.find_wet_hours`), and joins the per-gauge metrics
    from :func:`gauge_intensities`. One row per gauge (Point geometry,
    EPSG:4326) with ``total_mm`` and ``i{d}_mmph`` columns -- ready to overlay
    on the radar figure or feed the radar-vs-gauge comparison.
    """
    stations = get_stations(aoi, token=token, status=status, pad_deg=pad_deg)
    if not len(stations):
        return stations
    ids: List[str] = stations["station_id"].tolist()
    rain: Dict[str, pd.DataFrame] = {}
    for i in range(0, len(ids), batch):                # chunk to limit hours
        rain.update(get_rainfall(ids[i:i + batch], starttime, endtime,
                                 token=token))
    recs = [gauge_intensities(rain.get(sid, pd.DataFrame()), durations)
            for sid in ids]
    metrics = pd.DataFrame(recs, index=stations.index)
    return stations.join(metrics)


# --------------------------------------------------------------------------- #
# per-gauge rainfall TIME SERIES (counterpart to mrms.virtual_gauge_timeseries)
# --------------------------------------------------------------------------- #
def _series_from_per_minute(per_minute, times, starttime, endtime, durations):
    """1-min per-minute precip (mm) -> df (rate/total/i{d}) on the full window.

    Outside the gauge's own record the rate is 0 (no precip reported). The
    15-minute intensity uses the ``(i16+i14)/2`` estimator; others are plain
    trailing windows -- identical to the radar/virtual-gauge reducers so the two
    are directly comparable.
    """
    idx = pd.to_datetime(times)
    if idx.tz is not None:
        idx = idx.tz_localize(None)
    grid = pd.date_range(starttime, endtime, freq="1min")
    if grid.tz is not None:                      # match the (naive) obs index, so
        grid = grid.tz_localize(None)            # tz-aware start/end can't zero out
    pm = pd.Series(np.asarray(per_minute, float), index=idx).reindex(grid).fillna(0.0)
    total = pm.cumsum()
    df = pd.DataFrame({"rate_mmph": pm * 60.0, "total_mm": total}, index=grid)
    for d in durations:
        if d == 15:
            df["i15_mmph"] = (total.diff(16) * 60.0 / 16.0
                              + total.diff(14) * 60.0 / 14.0) / 2.0
        else:
            df[f"i{d}_mmph"] = total.diff(d) * 60.0 / d
    return df


def gauge_timeseries(aoi, starttime, endtime, token: Optional[str] = None,
                     durations: Sequence[int] = (5, 15, 30, 60),
                     status: Optional[str] = "ACTIVE", batch: int = 50,
                     pad_deg: float = 0.05) -> Tuple[Dict[str, "pd.DataFrame"],
                                                     gpd.GeoDataFrame]:
    """Per-gauge rainfall **time series** for every Synoptic station in the AOI.

    The real-gauge counterpart to :func:`stormscape.mrms.virtual_gauge_timeseries`
    -- returns ``({name: DataFrame}, stations)`` where each DataFrame is indexed
    at 1-minute UTC with ``rate_mmph``, ``total_mm`` and ``i{d}_mmph`` (matching
    the radar reducers), and ``df.attrs`` carries ``station_id``/``lon``/``lat``
    for matching against virtual gauges. Stations reporting no precip are skipped.
    """
    stations = get_stations(aoi, token=token, status=status, pad_deg=pad_deg)
    series: Dict[str, pd.DataFrame] = {}
    if not len(stations):
        return series, stations
    ids = stations["station_id"].tolist()
    rain: Dict[str, pd.DataFrame] = {}
    for i in range(0, len(ids), batch):
        rain.update(get_rainfall(ids[i:i + batch], starttime, endtime, token=token))
    for _, row in stations.iterrows():
        obs = rain.get(row["station_id"])
        times, _cum, per_min = (precipitation_per_minute(obs)
                                if obs is not None else ([], None, []))
        if len(times) < 2:
            continue
        df = _series_from_per_minute(per_min, times, starttime, endtime, durations)
        df.attrs.update(station_id=str(row["station_id"]),
                        lon=float(row.geometry.x), lat=float(row.geometry.y),
                        report_min=_report_min(obs))
        name = str(row["name"]) if row.get("name") else str(row["station_id"])
        series[name] = df
    return series, stations


def _safe_name(name) -> str:
    """Filesystem-safe gauge stem (matches ``stormscape.cli._safe_name``)."""
    return "".join(ch if ch.isalnum() else "_" for ch in str(name))


def unique_safe_names(names) -> Dict[str, str]:
    """Map gauge names -> filesystem stems that stay unique even on a
    case-INSENSITIVE filesystem (macOS APFS/HFS+): two stations whose names differ
    only in case -- e.g. ``'Virginia City'`` vs ``'VIRGINIA CITY'`` -- have
    :func:`_safe_name` stems differing only in case, which collapse to one file and
    silently overwrite (one gauge then vanishes from the per-gauge CSVs / detail
    figures). Suffix the 2nd+ of each case-insensitive group with ``_2``, ``_3`` ...
    Deterministic (sorted input) so the writer and reader agree."""
    seen: Dict[str, int] = {}
    out: Dict[str, str] = {}
    for n in sorted(names, key=str):
        base = _safe_name(n)
        k = seen.get(base.lower(), 0)
        out[n] = base if k == 0 else f"{base}_{k + 1}"
        seen[base.lower()] = k + 1
    return out


def fetch_gauge_event(aoi, starttime, endtime, out_dir: str, key: str,
                      token: Optional[str] = None,
                      durations: Sequence[int] = DEFAULT_DURATIONS,
                      series_durations: Sequence[int] = (5, 15, 30, 60),
                      status: Optional[str] = "ACTIVE", batch: int = 50,
                      pad_deg: float = 0.05, write_series: bool = True,
                      layout=None
                      ) -> Tuple[gpd.GeoDataFrame, Dict[str, "pd.DataFrame"]]:
    """**One** Synoptic fetch -> the canonical gauge store for an event.

    A single ``get_stations`` + ``get_rainfall`` whose shared raw obs produce
    *both* the per-gauge peak metrics and the 1-minute series, written to one
    consistent place that every other step reads (no second fetch, no drift):

    * ``<out_dir>/<key>_gauges.geojson`` -- one row per station: coords + peak
      ``total_mm`` / ``i{d}_mmph`` / native ``report_min`` + the **I15
      time-of-peak** (``i15_peak_time``, ISO-UTC).
    * ``<out_dir>/RainGaugeData/<key>_gauge_<name>.csv`` -- the 1-min series
      (``rate_mmph``/``total_mm``/``i{d}_mmph``) with ``lon``/``lat``/
      ``station_id`` columns embedded so each file is self-describing.

    Returns ``(stations_gdf_with_metrics, {name: series_df})``. This supersedes
    calling :func:`gauge_fields` and :func:`gauge_timeseries` separately (which
    each hit Synoptic independently and could return different station sets).
    """
    stations = get_stations(aoi, token=token, status=status, pad_deg=pad_deg)
    if not len(stations):
        return stations, {}
    ids = stations["station_id"].tolist()
    rain: Dict[str, pd.DataFrame] = {}
    for i in range(0, len(ids), batch):
        rain.update(get_rainfall(ids[i:i + batch], starttime, endtime, token=token))

    metrics, series = [], {}
    for _, row in stations.iterrows():
        obs = rain.get(row["station_id"])
        m = gauge_intensities(obs if obs is not None else pd.DataFrame(), durations)
        tpk = None
        if obs is not None:
            times, _cum, per_min = precipitation_per_minute(obs)
            if len(times) >= 2:
                df = _series_from_per_minute(per_min, times, starttime, endtime,
                                             series_durations)
                df.attrs.update(station_id=str(row["station_id"]),
                                lon=float(row.geometry.x),
                                lat=float(row.geometry.y),
                                report_min=m["report_min"])
                series[str(row["name"]) if row.get("name")
                       else str(row["station_id"])] = df
                if "i15_mmph" in df and df["i15_mmph"].notna().any():
                    tpk = df["i15_mmph"].idxmax()
        m["i15_peak_time"] = tpk.isoformat() if tpk is not None else None
        metrics.append(m)

    gdf = stations.join(pd.DataFrame(metrics, index=stations.index))
    os.makedirs(out_dir, exist_ok=True)
    gpath = out_path(out_dir, f"{key}_gauges.geojson", layout=layout)
    gdf.to_file(gpath, driver="GeoJSON")
    if write_series:
        rgd = os.path.join(out_dir, "RainGaugeData")
        os.makedirs(rgd, exist_ok=True)
        stems = unique_safe_names(series.keys())   # case-insensitive-unique stems
        for name, df in series.items():
            out = df.copy()
            out["station_id"] = df.attrs.get("station_id")
            out["lon"] = df.attrs.get("lon")
            out["lat"] = df.attrs.get("lat")
            out["report_min"] = df.attrs.get("report_min")   # native gauge tempo
            out.to_csv(os.path.join(rgd, f"{key}_gauge_{stems[name]}.csv"))
    return gdf, series


def load_event_series(rgd_dir: str, key: str, gauges: gpd.GeoDataFrame
                      ) -> Dict[str, "pd.DataFrame"]:
    """Reload an event's per-gauge series from the saved RainGaugeData CSVs.

    The temporal counterpart to reusing saved rasters -- no Synoptic re-fetch.
    Reads ``<rgd_dir>/<key>_gauge_<name>.csv`` for each station in ``gauges`` and
    re-attaches ``lon``/``lat``/``station_id`` to ``df.attrs`` (from the gauges
    GeoDataFrame) so the series match against virtual gauges exactly as a fresh
    :func:`gauge_timeseries` would. Returns ``{name: DataFrame}``.
    """
    series: Dict[str, pd.DataFrame] = {}
    if not os.path.isdir(rgd_dir):
        return series
    g = gauges.to_crs(4326) if getattr(gauges, "crs", None) is not None else gauges
    names = [str(r.get("name") or r.get("station_id")) for _, r in g.iterrows()]
    stems = unique_safe_names(names)               # must match fetch_gauge_event
    for _, r in g.iterrows():
        name = str(r.get("name") or r.get("station_id"))
        p = os.path.join(rgd_dir, f"{key}_gauge_{stems[name]}.csv")
        if not os.path.exists(p):
            continue
        df = pd.read_csv(p, index_col=0, parse_dates=True)
        # native gauge tempo: prefer the authoritative geojson field, fall back to
        # the CSV column (older stores embedded neither -> NaN, handled downstream).
        try:
            rm = float(r.get("report_min"))
        except (TypeError, ValueError):
            rm = float("nan")
        if not np.isfinite(rm) and "report_min" in df.columns:
            vals = pd.to_numeric(df["report_min"], errors="coerce").dropna()
            rm = float(vals.iloc[0]) if len(vals) else float("nan")
        df.attrs.update(station_id=str(r.get("station_id")),
                        lon=float(r.geometry.x), lat=float(r.geometry.y),
                        report_min=rm)
        series[name] = df
    return series


def storm_window(series: Dict[str, "pd.DataFrame"], pad_min: int = 30,
                 rate_col: str = "rate_mmph", cover: float = 0.99) -> Optional[Tuple]:
    """Tight ``(start, end)`` window holding the central ``cover`` fraction of the
    total rainfall *mass* across a ``{name: 1-min DataFrame}`` set (from
    :func:`load_event_series` / :func:`fetch_gauge_event`), padded by ``pad_min``
    minutes; ``None`` if dry.

    Mass-weighted (not any-nonzero) so isolated single-gauge tips far from the
    storm don't stretch the window -- it captures the storm *duration*, letting a
    full storm-DAY record be trimmed for analyses like the virtual-gauge
    comparison (the temporal analog of a spatial zoom).
    """
    agg = None
    for df in series.values():
        if rate_col not in getattr(df, "columns", []):
            continue
        inc = pd.to_numeric(df[rate_col], errors="coerce").fillna(0.0) / 60.0  # mm/min
        agg = inc.copy() if agg is None else agg.add(inc, fill_value=0.0)
    if agg is None or float(agg.sum()) <= 0:
        return None
    cum = agg.cumsum() / float(agg.sum())
    tail = (1.0 - float(cover)) / 2.0
    lo = agg.index[(cum >= tail).to_numpy()][0]
    hi = agg.index[(cum <= 1.0 - tail).to_numpy()][-1]
    pad = pd.Timedelta(minutes=int(pad_min))
    return (lo - pad, hi + pad)
