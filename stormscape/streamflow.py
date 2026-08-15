"""USGS stream gauges -- discharge and stage -- for an AOI + storm window.

The channel-side counterpart to :mod:`stormscape.gauges`. Where that pulls what
fell out of the sky at a point, this pulls what came down the channel: the
**discharge** and **gage height** records of the USGS stream gauges inside the
AOI, over the same window, so a hydrograph can be read against the rainfall
that produced it.

For post-fire work the pairing is the whole point. A burned catchment converts
rain to runoff at a rate an unburned one does not, so the interesting quantity
is not the peak i15 or the peak discharge alone but the relationship between
them -- and whether the response arrived minutes or hours after the cell.

Two sources, one interface
--------------------------
USGS is mid-migration, so :data:`SOURCES` carries both and the choice is a
keyword rather than a rewrite:

``nwis`` (default)
    The legacy **Water Services** endpoints (``waterservices.usgs.gov``).
    No API key, no rate limit worth worrying about, and the format the whole
    hydrology world already parses. USGS has announced these will be
    decommissioned, but the timeline is explicitly uncertain.

``ogc``
    The modernized **Water Data OGC API** (``api.waterdata.usgs.gov``) that
    replaces it. Works anonymously at **100 requests/hour per IP**, which is
    ample here (one site query plus a handful of series per AOI); a free key
    raises that to 1,000/hour and is read from ``$STORMSCAPE_USGS_API_KEY``.
    It is still served under a ``/v0/`` path, which is why it is not yet the
    default -- but it is implemented and tested, not a stub, so switching is
    ``source="ogc"`` on the day the legacy service goes away.

Units
-----
Every series carries **both** unit systems -- ``discharge_cms``/``discharge_cfs``
and ``stage_m``/``stage_ft`` -- because the gauges are published in cubic feet
per second and feet, while the rest of stormscape is metric (mm/h, km²).
Summaries and figures default to **SI**, matching the rainfall side; pass
``units="cfs"`` (CLI ``--units cfs``) to read them the way the USGS gauge page
does. Nothing is lost either way: the conversion is exact and both columns are
always written.

Gotchas worth knowing
---------------------
* NWIS returns timestamps in the **gauge's local time with an offset**
  (``2026-08-12T17:00:00.000-07:00``), not UTC. Everything here is normalised to
  UTC on arrival so it aligns with MRMS, NEXRAD and the rain gauges, all of
  which are UTC.
* Missing values are the sentinel **-999999**, not null.
* The site service returns **one row per available time series**, so a gauge
  with discharge, stage and temperature appears three times -- dedupe on
  ``site_no`` or the count is meaningless.
* Cadence varies by gauge (5-minute on the Truckee, 15-minute on the Carson);
  ``report_min`` records the measured interval, the same way it does for rain
  gauges, because a coarse gauge smooths a flashy peak.
"""

from __future__ import annotations

import io
import os
import warnings
from datetime import datetime, timezone
from typing import Dict, Optional, Sequence

import geopandas as gpd
import numpy as np
import pandas as pd
import requests

from .aoi import load_aoi
from .layout import find_subdir, out_path, subdir

warnings.filterwarnings("ignore")

NWIS_BASE = "https://waterservices.usgs.gov/nwis"
OGC_BASE = "https://api.waterdata.usgs.gov/ogcapi/v0"

#: USGS parameter codes we read. Discharge and gage height are the two a
#: storm-response question actually needs; everything else is noise here.
DISCHARGE = "00060"
STAGE = "00065"
PARAMETERS = {DISCHARGE: "discharge", STAGE: "stage"}

CFS_TO_CMS = 0.028316846592      # exact
FT_TO_M = 0.3048                 # exact
MI2_TO_KM2 = 2.589988110336      # exact
NWIS_NODATA = -999999.0

#: series column -> (SI column, native column, SI unit, native unit)
UNIT_COLUMNS = {
    "discharge": ("discharge_cms", "discharge_cfs", "m$^3$ s$^{-1}$", "ft$^3$ s$^{-1}$"),
    "stage": ("stage_m", "stage_ft", "m", "ft"),
}

SOURCES = {
    "nwis": dict(label="USGS NWIS Water Services (legacy)", key_env=None,
                 note="no API key; announced for decommission, timeline uncertain"),
    "ogc": dict(label="USGS Water Data OGC API v0", key_env="STORMSCAPE_USGS_API_KEY",
                note="100 req/hour anonymous, 1000 with a free key; v0 path"),
}

STREAM_SUBDIR = "StreamGaugeData"

#: The legacy NWIS *site* service is slow over a bounding box -- 30-45 s is
#: normal, against sub-second for the OGC equivalent -- so it gets its own,
#: longer default. The value-series endpoints are quick.
SITE_TIMEOUT = 180

#: A gauge counts as active if its discharge record reaches within this many
#: days of now. Used to make the OGC source's site set match what NWIS returns
#: for ``siteStatus=active``.
ACTIVE_WITHIN_DAYS = 365


def _utc(t) -> datetime:
    """Coerce to a timezone-aware UTC datetime."""
    ts = pd.Timestamp(t)
    ts = ts.tz_localize("UTC") if ts.tz is None else ts.tz_convert("UTC")
    return ts.to_pydatetime()


def _api_key(source: str, api_key: Optional[str]) -> Optional[str]:
    """Explicit key wins, else the source's env var. Never logged."""
    if api_key:
        return api_key
    env = SOURCES.get(source, {}).get("key_env")
    return os.environ.get(env) if env else None


def _get(url, params, timeout=60, what="usgs"):
    """GET returning parsed JSON (or text for RDB); warns + returns None on failure."""
    try:
        r = requests.get(url, params=params, timeout=timeout)
    except Exception as e:                                     # noqa: BLE001
        _fail(what, repr(e)[:120])
        return None
    if r.status_code != 200:
        # 404 from NWIS means "no sites match", which is a legitimate answer
        if r.status_code == 404:
            return ""
        _fail(what, f"HTTP {r.status_code} {r.text[:120]}")
        return None
    return r


def _fail(what, detail):
    """Report loudly: this module's warnings would otherwise be swallowed by the
    module-level ``filterwarnings('ignore')`` that the geo stack pulls in."""
    import sys
    msg = f"USGS query failed ({what}): {detail}"
    warnings.warn(msg)
    print(f"warning: {msg}", file=sys.stderr)


# --------------------------------------------------------------------------- #
# sites
# --------------------------------------------------------------------------- #
def _rdb_to_frame(text: str) -> pd.DataFrame:
    """Parse USGS RDB: '#' comments, a header row, then a format row to drop."""
    lines = [ln for ln in (text or "").splitlines() if not ln.startswith("#")]
    if len(lines) < 2:
        return pd.DataFrame()
    body = "\n".join([lines[0]] + lines[2:])          # drop the '5s/15s' row
    return pd.read_csv(io.StringIO(body), sep="\t", dtype=str)


def _num(s):
    return pd.to_numeric(s, errors="coerce") if s is not None else np.nan


def _sites_nwis(bounds, timeout, active_only) -> pd.DataFrame:
    # `outputDataTypeCd` and `siteOutput=expanded` are mutually exclusive -- the
    # pair is an HTTP 400. `hasDataTypeCd=iv` already restricts to gauges with
    # instantaneous values, so the former is redundant anyway.
    params = dict(format="rdb",
                  bBox=",".join(f"{b:.6f}" for b in bounds),
                  siteType="ST", hasDataTypeCd="iv",
                  parameterCd=DISCHARGE,
                  siteOutput="expanded",
                  siteStatus="active" if active_only else "all")
    r = _get(f"{NWIS_BASE}/site/", params, timeout, "nwis sites")
    if r is None:
        return pd.DataFrame()
    df = _rdb_to_frame(r if isinstance(r, str) else r.text)
    if not len(df):
        return df
    # one row per available time series -> a gauge appears once per parameter
    df = df.drop_duplicates("site_no")
    return pd.DataFrame(dict(
        site_no=df.site_no.astype(str).str.strip(),
        name=df.station_nm.astype(str).str.strip(),
        lat=_num(df.dec_lat_va), lon=_num(df.dec_long_va),
        drain_area_km2=_num(df.get("drain_area_va")) * MI2_TO_KM2,
        alt_ft=_num(df.get("alt_va")),
        huc=df.get("huc_cd", pd.Series(index=df.index, dtype=str)),
        site_type=df.get("site_tp_cd", pd.Series(index=df.index, dtype=str)),
    ))


def _ogc_pages(url, params, timeout, what, api_key, max_pages=40):
    """Yield features across the OGC API's `next`-link pagination."""
    p = dict(params)
    if api_key:
        p["api_key"] = api_key
    for _ in range(max_pages):
        r = _get(url, p, timeout, what)
        if r is None or r == "":
            return
        try:
            payload = r.json()
        except Exception as e:                                 # noqa: BLE001
            _fail(what, f"bad JSON: {repr(e)[:80]}")
            return
        feats = payload.get("features") or []
        yield from feats
        nxt = next((l.get("href") for l in payload.get("links", [])
                    if l.get("rel") == "next"), None)
        if not nxt or not feats:
            return
        url, p = nxt, {}          # the next href already carries the query


def _ogc_discharge_sites(bounds, timeout, active_only, api_key) -> pd.DataFrame:
    """Gauges with an *instantaneous* discharge series, + period of record.

    ``monitoring-locations`` alone answers "is this a stream site", not "does it
    measure discharge every 15 minutes", so filtering on it gives a set six
    times larger than the NWIS equivalent -- mostly gauges that never recorded
    discharge, plus records that ended decades ago. The time-series metadata
    collection carries both facts, so the two sources return the same gauges.
    ``computation_period_identifier="Points"`` is the OGC name for what NWIS
    calls instantaneous values.
    """
    params = dict(bbox=",".join(f"{b:.6f}" for b in bounds),
                  parameter_code=DISCHARGE,
                  computation_period_identifier="Points", limit=1000, f="json")
    span: Dict[str, list] = {}
    for f in _ogc_pages(f"{OGC_BASE}/collections/time-series-metadata/items",
                        params, timeout, "ogc time-series metadata", api_key):
        p = f.get("properties", {})
        sid = str(p.get("monitoring_location_id") or "").replace("USGS-", "").strip()
        if not sid:
            continue
        beg = pd.to_datetime(p.get("begin") or p.get("begin_utc"), utc=True,
                             errors="coerce")
        end = pd.to_datetime(p.get("end") or p.get("end_utc"), utc=True,
                             errors="coerce")
        prev = span.get(sid)
        span[sid] = [min(beg, prev[0]) if prev and pd.notna(prev[0]) else beg,
                     max(end, prev[1]) if prev and pd.notna(prev[1]) else end]
    if not span:
        return pd.DataFrame()
    df = pd.DataFrame([dict(site_no=s, begin=v[0], end=v[1])
                       for s, v in span.items()])
    if active_only:
        cutoff = pd.Timestamp.now(tz="UTC") - pd.Timedelta(days=ACTIVE_WITHIN_DAYS)
        df = df[df.end.notna() & (df.end >= cutoff)]
    return df.reset_index(drop=True)


def _sites_ogc(bounds, timeout, active_only, api_key) -> pd.DataFrame:
    want = _ogc_discharge_sites(bounds, timeout, active_only, api_key)
    if not len(want):
        return pd.DataFrame()
    keep = set(want.site_no)
    params = dict(bbox=",".join(f"{b:.6f}" for b in bounds),
                  site_type_code="ST", limit=1000, f="json")
    rows = []
    for f in _ogc_pages(f"{OGC_BASE}/collections/monitoring-locations/items",
                        params, timeout, "ogc sites", api_key):
        p = f.get("properties", {})
        sid = str(p.get("monitoring_location_number") or "").strip()
        if sid not in keep:
            continue
        geom = f.get("geometry") or {}
        crd = geom.get("coordinates") or [np.nan, np.nan]
        rows.append(dict(
            site_no=sid,
            name=str(p.get("monitoring_location_name") or "").strip(),
            lon=crd[0], lat=crd[1],
            drain_area_km2=_num(pd.Series([p.get("drainage_area")])).iloc[0]
            * MI2_TO_KM2,
            alt_ft=_num(pd.Series([p.get("altitude")])).iloc[0],
            huc=p.get("hydrologic_unit_code"),
            site_type=p.get("site_type_code")))
    df = pd.DataFrame(rows)
    if not len(df):
        return df
    return df.drop_duplicates("site_no").merge(want, on="site_no", how="left")


def stream_sites(aoi, source="nwis", active_only=True, api_key=None,
                 pad_deg=0.0, timeout=SITE_TIMEOUT) -> gpd.GeoDataFrame:
    """USGS stream gauges inside ``aoi`` -> GeoDataFrame (EPSG:4326).

    Columns: ``site_no``, ``name``, ``drain_area_km2``, ``alt_ft``, ``huc``,
    ``site_type``. ``active_only`` keeps gauges currently reporting; pass
    ``False`` to include discontinued records (useful for historical events).
    """
    if source not in SOURCES:
        raise ValueError(f"unknown source {source!r}; have {sorted(SOURCES)}")
    bounds, _ = load_aoi(aoi, pad_deg=pad_deg)
    key = _api_key(source, api_key)
    df = (_sites_nwis(bounds, timeout, active_only) if source == "nwis"
          else _sites_ogc(bounds, timeout, active_only, key))
    if not len(df):
        return _empty_sites()
    df = df[df.lat.notna() & df.lon.notna()]
    # the OGC bbox filter is generous at the edges; clip to the AOI either way
    w, s, e, n = bounds
    df = df[(df.lon >= w) & (df.lon <= e) & (df.lat >= s) & (df.lat <= n)]
    g = gpd.GeoDataFrame(df.reset_index(drop=True),
                         geometry=gpd.points_from_xy(df.lon, df.lat), crs=4326)
    g["source"] = source
    return g.reset_index(drop=True)


def _empty_sites() -> gpd.GeoDataFrame:
    cols = ("site_no", "name", "lat", "lon", "drain_area_km2", "alt_ft", "huc",
            "site_type", "source")
    return gpd.GeoDataFrame({c: pd.Series(dtype="object") for c in cols},
                            geometry=gpd.GeoSeries([], crs=4326), crs=4326)


# --------------------------------------------------------------------------- #
# time series
# --------------------------------------------------------------------------- #
def _empty_series() -> pd.DataFrame:
    return pd.DataFrame(
        {c: pd.Series(dtype="float64") for c in
         ("discharge_cms", "discharge_cfs", "stage_m", "stage_ft")},
        index=pd.DatetimeIndex([], tz="UTC", name="time"))


def _frame_from_pairs(pairs: Dict[str, list]) -> pd.DataFrame:
    """``{param_code: [(utc_time, value), ...]}`` -> one aligned UTC frame."""
    out = None
    for code, rows in pairs.items():
        if not rows:
            continue
        name = PARAMETERS.get(code, code)
        idx = pd.DatetimeIndex([r[0] for r in rows], tz="UTC", name="time")
        s = pd.Series([r[1] for r in rows], index=idx, dtype="float64")
        s = s[~s.index.duplicated(keep="last")].sort_index()
        col = pd.DataFrame({name: s})
        out = col if out is None else out.join(col, how="outer")
    if out is None:
        return _empty_series()
    # native units come off the wire; carry both so nothing has to be re-derived
    if "discharge" in out:
        out["discharge_cfs"] = out.pop("discharge")
        out["discharge_cms"] = out.discharge_cfs * CFS_TO_CMS
    if "stage" in out:
        out["stage_ft"] = out.pop("stage")
        out["stage_m"] = out.stage_ft * FT_TO_M
    cols = [c for c in ("discharge_cms", "discharge_cfs", "stage_m", "stage_ft")
            if c in out]
    return out[cols].sort_index()


def _series_nwis(site_ids, start, end, timeout) -> Dict[str, pd.DataFrame]:
    params = dict(format="json", sites=",".join(site_ids),
                  parameterCd=f"{DISCHARGE},{STAGE}",
                  startDT=_utc(start).strftime("%Y-%m-%dT%H:%MZ"),
                  endDT=_utc(end).strftime("%Y-%m-%dT%H:%MZ"))
    r = _get(f"{NWIS_BASE}/iv/", params, timeout, "nwis iv")
    if r is None or r == "":
        return {}
    try:
        payload = r.json()
    except Exception as e:                                     # noqa: BLE001
        _fail("nwis iv", f"bad JSON: {repr(e)[:80]}")
        return {}
    per_site: Dict[str, Dict[str, list]] = {}
    for ts in payload.get("value", {}).get("timeSeries", []):
        site = ts["sourceInfo"]["siteCode"][0]["value"]
        code = ts["variable"]["variableCode"][0]["value"]
        rows = []
        for v in ts.get("values", [{}])[0].get("value", []):
            val = pd.to_numeric(v.get("value"), errors="coerce")
            if not np.isfinite(val) or val == NWIS_NODATA:
                continue
            # NWIS stamps local time WITH an offset; normalise to UTC
            rows.append((pd.Timestamp(v["dateTime"]).tz_convert("UTC"), float(val)))
        per_site.setdefault(site, {})[code] = rows
    return {s: _frame_from_pairs(p) for s, p in per_site.items()}


def _series_ogc(site_ids, start, end, timeout, api_key) -> Dict[str, pd.DataFrame]:
    span = (f"{_utc(start).strftime('%Y-%m-%dT%H:%M:%SZ')}/"
            f"{_utc(end).strftime('%Y-%m-%dT%H:%M:%SZ')}")
    out: Dict[str, pd.DataFrame] = {}
    for site in site_ids:
        pairs: Dict[str, list] = {}
        for code in (DISCHARGE, STAGE):
            params = dict(monitoring_location_id=f"USGS-{site}",
                          parameter_code=code, datetime=span,
                          limit=10000, f="json")
            rows = []
            for f in _ogc_pages(f"{OGC_BASE}/collections/continuous/items",
                                params, timeout, f"ogc {site} {code}", api_key):
                p = f.get("properties", {})
                val = pd.to_numeric(p.get("value"), errors="coerce")
                if not np.isfinite(val):
                    continue
                rows.append((pd.Timestamp(p["time"]).tz_convert("UTC"), float(val)))
            pairs[code] = rows
        frame = _frame_from_pairs(pairs)
        if len(frame):
            out[site] = frame
    return out


def stream_series(site_ids: Sequence[str], start, end, source="nwis",
                  api_key=None, timeout=60) -> Dict[str, pd.DataFrame]:
    """Discharge + stage series per gauge -> ``{site_no: DataFrame}``.

    Frames are indexed by **UTC** time and carry ``discharge_cms``,
    ``discharge_cfs``, ``stage_m``, ``stage_ft`` (whichever the gauge reports).
    """
    if source not in SOURCES:
        raise ValueError(f"unknown source {source!r}; have {sorted(SOURCES)}")
    ids = [str(s).strip() for s in site_ids if str(s).strip()]
    if not ids:
        return {}
    key = _api_key(source, api_key)
    return (_series_nwis(ids, start, end, timeout) if source == "nwis"
            else _series_ogc(ids, start, end, timeout, key))


# --------------------------------------------------------------------------- #
# summaries
# --------------------------------------------------------------------------- #
def _report_min(idx) -> float:
    """Median minutes between observations -- the gauge's measured cadence.

    Measured rather than taken from metadata, for the same reason the rain-gauge
    side measures it: a coarse gauge smooths a flashy peak, and the reader needs
    to know which peaks to trust.
    """
    if idx is None or len(idx) < 2:
        return float("nan")
    d = pd.Series(idx).diff().dropna().dt.total_seconds() / 60.0
    return float(d.median()) if len(d) else float("nan")


def flow_summary(sites: gpd.GeoDataFrame, series: Dict[str, pd.DataFrame],
                 units="si") -> gpd.GeoDataFrame:
    """Per-gauge storm summary joined onto the site geometry.

    Adds peak discharge and its time, peak stage, the rise above the first
    observation in the window (the storm response, as distinct from the
    absolute peak on a big river), unit discharge per km² where a drainage area
    is published, and the measured reporting cadence.
    """
    si = str(units).lower() != "cfs"
    qcol = "discharge_cms" if si else "discharge_cfs"
    hcol = "stage_m" if si else "stage_ft"
    rows = []
    for site, df in (series or {}).items():
        q = df[qcol] if qcol in df else pd.Series(dtype="float64")
        h = df[hcol] if hcol in df else pd.Series(dtype="float64")
        q = q.dropna()
        rec = dict(site_no=str(site),
                   peak_discharge=float(q.max()) if len(q) else np.nan,
                   peak_time=(q.idxmax().to_pydatetime() if len(q) else None),
                   start_discharge=float(q.iloc[0]) if len(q) else np.nan,
                   peak_stage=float(h.dropna().max()) if len(h.dropna()) else np.nan,
                   n_obs=int(len(df)), report_min=_report_min(df.index))
        rec["rise_discharge"] = (rec["peak_discharge"] - rec["start_discharge"]
                                 if len(q) else np.nan)
        rec["rise_ratio"] = (rec["peak_discharge"] / rec["start_discharge"]
                             if len(q) and rec["start_discharge"] > 0 else np.nan)
        # A peak sitting on the last observation is not a peak -- the window cut
        # the hydrograph while it was still rising, and every quantity derived
        # from it reads low. Nothing downstream can tell, so say so here.
        rec["peak_at_edge"] = bool(
            len(q) and rec["peak_time"] is not None
            and (df.index[-1] - q.idxmax()) <= pd.Timedelta(
                minutes=max(rec["report_min"] or 0, 1) * 1.5))
        rows.append(rec)
    summ = pd.DataFrame(rows)
    if sites is None or not len(sites):
        return gpd.GeoDataFrame(summ, geometry=gpd.GeoSeries([], crs=4326),
                                crs=4326) if not len(summ) else summ
    g = sites.copy()
    g["site_no"] = g.site_no.astype(str)
    if len(summ):
        g = g.merge(summ, on="site_no", how="left")
    else:
        for c in ("peak_discharge", "peak_time", "start_discharge", "peak_stage",
                  "n_obs", "report_min", "rise_discharge", "rise_ratio",
                  "peak_at_edge"):
            g[c] = np.nan
    with np.errstate(invalid="ignore", divide="ignore"):
        g["unit_discharge"] = g.peak_discharge / g.drain_area_km2
    g.attrs["units"] = "si" if si else "cfs"
    return gpd.GeoDataFrame(g, geometry=g.geometry, crs=4326)


# --------------------------------------------------------------------------- #
# canonical event store (mirrors gauges.fetch_gauge_event)
# --------------------------------------------------------------------------- #
def fetch_stream_event(aoi, start, end, out_dir: str, key: str, source="nwis",
                       active_only=True, api_key=None, units="si",
                       pad_deg=0.0, layout=None, timeout=SITE_TIMEOUT):
    """One fetch, one store: sites + series + summary, written once and reused.

    Writes ``<key>_streamgauges.geojson`` (locations + peak summary) and one
    self-describing CSV per gauge under ``StreamGaugeData/``. Returns
    ``(summary_gdf, {site_no: DataFrame})``.

    Deliberately mirrors :func:`stormscape.gauges.fetch_gauge_event`: the store
    is the whole requested window, and each analysis clips to what it needs.
    Two independent fetches of the same event are how the rain-gauge coordinates
    drifted apart once already.
    """
    sites = stream_sites(aoi, source=source, active_only=active_only,
                         api_key=api_key, pad_deg=pad_deg, timeout=timeout)
    if not len(sites):
        return _empty_sites(), {}
    series = stream_series(sites.site_no.tolist(), start, end, source=source,
                           api_key=api_key, timeout=timeout)
    summary = flow_summary(sites, series, units=units)

    os.makedirs(out_dir, exist_ok=True)
    gpath = out_path(out_dir, f"{key}_streamgauges.geojson", layout=layout)
    out = summary.copy()
    if "peak_time" in out:
        out["peak_time"] = out.peak_time.astype(str)
    out.to_file(gpath, driver="GeoJSON")

    sdir = subdir(out_dir, STREAM_SUBDIR, layout=layout)
    names = dict(zip(sites.site_no.astype(str), sites.name.astype(str)))
    for site, df in series.items():
        d = df.copy()
        d["site_no"] = site
        d["name"] = names.get(site, site)
        d.to_csv(os.path.join(sdir, f"{key}_stream_{site}.csv"))
    return summary, series


def load_event_series(out_dir: str, key: str) -> Dict[str, pd.DataFrame]:
    """Reload the per-gauge CSVs written by :func:`fetch_stream_event`."""
    sdir = find_subdir(out_dir, STREAM_SUBDIR)
    if not sdir or not os.path.isdir(sdir):
        return {}
    out = {}
    prefix = f"{key}_stream_"
    for fn in sorted(os.listdir(sdir)):
        if not (fn.startswith(prefix) and fn.endswith(".csv")):
            continue
        df = pd.read_csv(os.path.join(sdir, fn), index_col=0, parse_dates=True)
        if df.index.tz is None:
            df.index = df.index.tz_localize("UTC")
        site = fn[len(prefix):-4]
        if "name" in df:
            df.attrs["name"] = str(df.name.iloc[0]) if len(df) else site
        out[site] = df
    return out


def describe_sources():
    """Print the source registry -- which is default, which needs a key."""
    for name, spec in SOURCES.items():
        env = spec.get("key_env")
        state = "no key needed" if not env else (
            "key set" if os.environ.get(env) else f"optional ${env} (not set)")
        print(f"{name:<6} {spec['label']}")
        print(f"{'':<6}  {spec['note']}")
        print(f"{'':<6}  {state}")
