"""MRMS radar rainfall -> peak 15-minute intensity (i15) fields.

Builds, for one storm-day over an AOI, a map of the peak 15-minute rainfall
intensity ``i15_max`` (mm/h) plus companion fields, by stacking the 2-minute
MRMS ``PrecipRate`` returns through a rolling 15-minute estimator. Adapted
from D. Cavagna's ``MRMS_stack.py``: the i15 estimator and the NOAA-S3 URL
scheme are reused; here it is AOI-windowed, parallelised, and driven by a
date rather than a fixed file list.

i15 estimator
-------------
MRMS ``PrecipRate`` is a 2-min instantaneous rate (mm/h); ``a2 = rate * 2/60``
is the 2-min accumulation (mm). Over a trailing 16-min window (8 steps):
``i16 = sum(8) * 60/16`` and ``i14 = sum(last 7) * 60/14``; ``i15 = mean(i16,
i14)``. The running maximum over the storm gives ``i15_max``. The 30- and
60-minute peaks (``i30max``/``i60max``) use plain trailing windows (15 / 30
steps) scaled to mm/h, kept as per-cell running maxima alongside ``i15max``.

Storm-window detection
----------------------
Hourly ``RadarOnly_QPE_01H`` is scanned over the UTC window covering the
local calendar day (or an explicit ``window``); each hour's areal maximum over
the AOI is taken, the wettest hours (> ``qpe_thresh``, capped at
``max_wet_hours``) are kept, and 2-min ``PrecipRate`` is stacked over each
contiguous run. The cap ranks by *intensity*, so it discards a long storm's
weakest hours and shortens the stacked span; :func:`find_wet_hours` warns
whenever it binds.

**Hourly QPE is stamped at the END of the hour it describes.** ``QPE(HH)`` is
the rain that fell in ``[HH-1, HH]``. Verified against the 2-min ``PrecipRate``
accumulation on independent hours: the ratio over ``[HH-1, HH]`` is 1.000,
while ``[HH, HH+1]`` gives 0.13-0.33. So a run of wet stamps ``[h0..hn]``
describes rain over ``[h0-1h, hn]``, and that -- plus a 14-minute lead so the
rolling i15 is defined from the run's first wet minute -- is the span stacked.
Getting this backwards silently truncates the front of a storm.

Outputs are AOI-clipped GeoTIFFs in EPSG:4326 (the native MRMS grid):
``i15max, i30max, i60max, i2max, total, tpki15`` (UTC hour of peak i15),
``rqi`` (radar quality index), ``shsr`` (seamless hybrid-scan beam height, km).
"""

from __future__ import annotations

import datetime as dt
import gzip
import io
import os
import time
import warnings
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
import rasterio
import requests
from rasterio.windows import Window

from .aoi import load_aoi
from .layout import out_path

warnings.filterwarnings("ignore")

S3 = "https://noaa-mrms-pds.s3.amazonaws.com/CONUS"

# CONUS MRMS grid (from grib metadata): UL (-130, 55), 0.01 deg, 7000 x 3500.
G_W, G_N, G_RES, G_NX, G_NY = -130.0, 55.0, 0.01, 7000, 3500

# (product subdirectory, file prefix); S3 key = <dir>/<YYYYMMDD>/<prefix>_<dt>
PRODUCTS = {
    "PrecipRate": ("PrecipRate_00.00", "MRMS_PrecipRate_00.00"),
    "RadarOnly":  ("RadarOnly_QPE_01H_00.00", "MRMS_RadarOnly_QPE_01H_00.00"),
    "RQI":        ("RadarQualityIndex_00.00", "MRMS_RadarQualityIndex_00.00"),
    "SHSR":       ("SeamlessHSRHeight_00.00", "MRMS_SeamlessHSRHeight_00.00"),
    # gauge-corrected hourly QPE (Pass-2 = final, more gauges; Pass-1 = early)
    "MultiSensor":  ("MultiSensor_QPE_01H_Pass2_00.00",
                     "MRMS_MultiSensor_QPE_01H_Pass2_00.00"),
    "MultiSensor1": ("MultiSensor_QPE_01H_Pass1_00.00",
                     "MRMS_MultiSensor_QPE_01H_Pass1_00.00"),
    # precip-typing / hail products (QPE-bias diagnosis): rain category, the
    # Z-R/rate relationship picked per cell, and max estimated hail size
    "PrecipFlag":   ("PrecipFlag_00.00", "MRMS_PrecipFlag_00.00"),
    "PrecipRateID": ("SyntheticPrecipRateID_00.00",
                     "MRMS_SyntheticPrecipRateID_00.00"),
    "MESH":         ("MESH_00.50", "MRMS_MESH_00.50"),
    "RadarOnly15M": ("RadarOnly_QPE_15M_00.00", "MRMS_RadarOnly_QPE_15M_00.00"),
}

# Defaults (overridable per call).
QPE_THRESH = 2.5            # mm; hourly areal-max above this = a "wet hour"
MAX_WET_HRS = 8             # cap processed wet hours (cost + i15-peak capture)
SCAN_PAD_H = (4, 10)        # UTC scan = [day 04:00, next-day 10:00] ~ local day
WORKERS = 12                # parallel MRMS downloads


class Missing(Exception):
    """File genuinely absent on the server (HTTP 404) -- do not retry."""


def parse_date(date) -> dt.date:
    """Accept a date/datetime, 'YYYYMMDD', or 'YYYY-MM-DD' -> datetime.date."""
    if isinstance(date, dt.datetime):
        return date.date()
    if isinstance(date, dt.date):
        return date
    s = str(date).strip().replace("-", "")
    return dt.datetime.strptime(s, "%Y%m%d").date()


def aoi_window(bounds):
    """rasterio Window + transform for a (W,S,E,N) AOI on the CONUS grid."""
    w, s, e, n = bounds
    c0 = max(int(np.floor((w - G_W) / G_RES)), 0)
    c1 = min(int(np.ceil((e - G_W) / G_RES)), G_NX)
    r0 = max(int(np.floor((G_N - n) / G_RES)), 0)
    r1 = min(int(np.ceil((G_N - s) / G_RES)), G_NY)
    win = Window(c0, r0, c1 - c0, r1 - r0)
    tr = rasterio.transform.from_origin(G_W + c0 * G_RES, G_N - r0 * G_RES,
                                        G_RES, G_RES)
    return win, tr


def fetch(product, t, win):
    """Download one MRMS grib2 -> AOI-windowed array (values < 0 -> NaN).

    Retries only transient failures (timeouts, dropped connections, 5xx). A
    404 means the timestep does not exist; we raise :class:`Missing` at once
    so absent files never burn backoff time.
    """
    date, hms = t.strftime("%Y%m%d"), t.strftime("%H%M%S")
    pdir, prefix = PRODUCTS[product]
    url = f"{S3}/{pdir}/{date}/{prefix}_{date}-{hms}.grib2.gz"
    err = "?"
    for k in range(3):
        try:
            r = requests.get(url, timeout=60)
            if r.status_code == 200:
                with gzip.GzipFile(fileobj=io.BytesIO(r.content)) as g:
                    raw = g.read()
                with rasterio.open(io.BytesIO(raw)) as ds:
                    a = ds.read(1, window=win).astype("float32")
                a[a < 0] = np.nan          # MRMS no-coverage / missing flags
                return a
            if r.status_code == 404:
                raise Missing(f"{product} {date}-{hms}")
            err = f"HTTP {r.status_code}"
        except Missing:
            raise
        except Exception as e:             # noqa: BLE001
            err = repr(e)[:120]
        time.sleep(1.5 * (k + 1))
    raise RuntimeError(f"{product} {date}-{hms}: {err}")


def fetch_many(product, times, win, workers=WORKERS):
    """Parallel fetch; returns {t: array} (absent/failed timesteps omitted)."""
    out = {}

    def one(t):
        try:
            return t, fetch(product, t, win)
        except Exception:                  # noqa: BLE001 - Missing or network
            return t, None

    with ThreadPoolExecutor(max_workers=workers) as ex:
        for t, a in ex.map(one, list(times)):
            if a is not None:
                out[t] = a
    return out


def compute_i15(stack):
    """i15 (mm/h) from a trailing list of >=8 2-minute accumulations (mm)."""
    s = np.dstack(stack[-8:])
    i16 = np.nansum(s, axis=2) * 60 / 16
    i14 = np.nansum(s[:, :, 1:], axis=2) * 60 / 14
    return (i16 + i14) / 2


def window_hours(date0=None, scan_pad_h=SCAN_PAD_H, window=None):
    """The hourly stamps to scan: an explicit ``window`` or a storm-day span.

    ``window`` is ``(start, end)`` UTC datetimes and wins when given; otherwise
    the span is ``[date0 scan_pad_h[0]Z, next-day scan_pad_h[1]Z]``, which
    covers the local calendar day.
    """
    if window is not None:
        start, end = window
        start = start.replace(minute=0, second=0, microsecond=0)
        n = int((end - start).total_seconds() // 3600) + 1
        return [start + dt.timedelta(hours=h) for h in range(max(n, 1))]
    if date0 is None:
        raise ValueError("need a date or an explicit window")
    start = dt.datetime(date0.year, date0.month, date0.day, scan_pad_h[0])
    return [start + dt.timedelta(hours=h)
            for h in range(24 + scan_pad_h[1] - scan_pad_h[0])]


def find_wet_hours(date0, win, qpe_thresh=QPE_THRESH, max_wet_hours=MAX_WET_HRS,
                   scan_pad_h=SCAN_PAD_H, workers=WORKERS, window=None):
    """Scan hourly QPE; return (wet-hour DataFrame, full scan DataFrame).

    Warns when more than ``max_wet_hours`` hours are wet: the cap keeps the most
    *intense* ones, so the storm's weak opening/closing hours are dropped and
    the stacked span shortens with nothing downstream to show for it.
    """
    hours = window_hours(date0, scan_pad_h, window)
    arrs = fetch_many("RadarOnly", hours, win, workers=workers)
    rec = pd.DataFrame(
        [(t, float(np.nanmax(arrs[t])) if t in arrs
          and np.isfinite(arrs[t]).any() else 0.0) for t in hours],
        columns=["t", "qmax"])
    wet = rec[rec.qmax > qpe_thresh].sort_values("qmax", ascending=False)
    if wet.empty:                          # fall back to the single best hour
        wet = rec.sort_values("qmax", ascending=False).head(1)
    if len(wet) > max_wet_hours:
        # The cap ranks by INTENSITY, so what it discards are the storm's
        # weakest hours -- typically its opening and closing tails. Those hours
        # also bound the stacked span: dropping a trailing wet stamp shortens
        # the contiguous run, so `total` loses that hour's rain outright and the
        # rolling i15/i30/i60 never see it. Nothing downstream can tell this
        # happened -- `n_wet_hr` reports what was KEPT -- so warn rather than
        # truncate quietly. Warn, don't refuse: the cap is a real cost control
        # and a caller may genuinely want the most intense hours only.
        dropped = wet.iloc[max_wet_hours:].sort_values("t")
        warnings.warn(
            f"{len(wet)} wet hours found but only the {max_wet_hours} most "
            f"intense are kept: dropping "
            f"{', '.join(t.strftime('%m-%d %HZ') for t in dropped.t)} "
            f"(weakest kept {wet.qmax.iloc[max_wet_hours - 1]:.1f} mm, "
            f"strongest dropped {dropped.qmax.max():.1f} mm). "
            "Storm totals and the stacked span will be short -- raise "
            "max_wet_hours (CLI --max-wet-hours) to keep the whole storm. "
            "Multi-storm windows essentially always need it raised.",
            stacklevel=2)
    wet = wet.head(max_wet_hours).sort_values("t")
    return wet, rec


def contiguous_runs(hour_list):
    """Group sorted wet hours into contiguous runs (a gap > 1 h splits)."""
    runs, cur = [], [hour_list[0]]
    for t in hour_list[1:]:
        if (t - cur[-1]) <= dt.timedelta(hours=1):
            cur.append(t)
        else:
            runs.append(cur)
            cur = [t]
    runs.append(cur)
    return runs


def wet_window(aoi, start, end, pad_deg=0.05, qpe_thresh=QPE_THRESH,
               pad_min=30, workers=WORKERS):
    """Tight ``(start, end)`` spanning the wet hours in ``[start, end]``; ``None``
    if the whole span is dry.

    The radar-side counterpart of :func:`stormscape.gauges.storm_window`, for
    bounding an expensive fetch *before* making it. One hourly ``RadarOnly`` QPE
    grid per hour over the AOI window is a few KB, so probing a 30-hour span
    costs far less than a single NEXRAD Level II volume -- which is the point: a
    storm-DAY window handed straight to
    :func:`stormscape.nexrad.virtual_gauge_timeseries` pulls ~10 volumes an hour
    for the whole day, nearly all of them dry.

    Spans the first to the last wet hour (plus ``pad_min``), so a day with two
    separate cells returns one window covering both rather than dropping the gap
    between them.
    """
    bounds, _ = load_aoi(aoi, pad_deg=pad_deg)
    win, _ = aoi_window(bounds)
    hours, t = [], start.replace(minute=0, second=0, microsecond=0)
    while t <= end:
        hours.append(t)
        t += dt.timedelta(hours=1)
    arrs = fetch_many("RadarOnly", hours, win, workers=workers)
    wet = [h for h in hours
           if h in arrs and np.isfinite(arrs[h]).any()
           and float(np.nanmax(arrs[h])) > qpe_thresh]
    if not wet:
        return None
    pad = dt.timedelta(minutes=int(pad_min))
    # the hourly QPE stamped at HH covers the hour ENDING at HH, so the wet
    # period itself starts an hour before the first wet stamp
    return (max(min(wet) - dt.timedelta(hours=1) - pad, start),
            min(max(wet) + pad, end))


def i15_storm_day(aoi, date=None, pad_deg=0.05, qpe_thresh=QPE_THRESH,
                  max_wet_hours=MAX_WET_HRS, scan_pad_h=SCAN_PAD_H,
                  workers=WORKERS, verbose=True, window=None):
    """Build the peak-i15 field (and companions) for one storm-day over an AOI.

    Parameters
    ----------
    aoi
        Anything :func:`stormscape.aoi.load_aoi` accepts.
    date
        Storm-day (local calendar day): date/datetime, 'YYYYMMDD', or
        'YYYY-MM-DD'. Optional when ``window`` is given.
    window
        Explicit ``(start, end)`` UTC datetimes to scan and stack, overriding
        ``date``/``scan_pad_h``. Use it when the storm does not line up with a
        local day -- back-to-back evening storms otherwise share the ~30 h
        storm-day span, and the previous evening's tail gets stacked into
        today's peak-intensity maps.
    pad_deg
        Degrees to pad the AOI bounds before windowing the MRMS grid.

    Returns
    -------
    dict
        ``fields`` ({name: ndarray} for i15max, i30max, i60max, i2max, total,
        tpki15, rqi, shsr), ``transform`` (rasterio Affine), ``crs``
        ('EPSG:4326'),
        ``profile`` (rasterio profile ready for ``rasterio.open(... , **p)``),
        and ``meta`` (scalar summary: peak UTC time, areal-max QPE, wet-hour
        count, AOI-max i15, median RQI).
    """
    bounds, _ = load_aoi(aoi, pad_deg=pad_deg)
    if date is None and window is None:
        raise ValueError("i15_storm_day needs a date or an explicit window")
    date0 = parse_date(date) if date is not None else window[0].date()
    win, tr = aoi_window(bounds)
    shape = (int(win.height), int(win.width))
    if shape[0] <= 0 or shape[1] <= 0:
        raise ValueError(f"AOI {bounds} is empty or outside the CONUS grid.")

    wet, scout = find_wet_hours(date0, win, qpe_thresh, max_wet_hours,
                                scan_pad_h, workers, window=window)
    peak_t = scout.loc[scout.qmax.idxmax(), "t"]
    if verbose:
        print(f"  {date0}: peak {peak_t:%m-%d %H}Z "
              f"qmax={scout.qmax.max():.1f} mm, {len(wet)} wet hr", flush=True)

    i15_max = np.zeros(shape, np.float32)
    i30_max = np.zeros(shape, np.float32)
    i60_max = np.zeros(shape, np.float32)
    i2_max = np.zeros(shape, np.float32)
    total = np.zeros(shape, np.float32)
    tpki15 = np.full(shape, np.nan, np.float32)

    for run in contiguous_runs(list(wet.t)):
        # Hourly QPE stamped HH covers the hour ENDING at HH (see the module
        # docstring), so a run of wet stamps [h0..hn] describes rain over
        # [h0-1h, hn] -- NOT [h0, hn+1h]. Reading the latter skipped up to
        # 46 min at the start of the first wet hour and spent a fetch on a
        # usually-dry trailing hour. The extra 14 min is the lead-in that makes
        # the rolling i15 valid from the run's first wet minute.
        t0 = run[0] - dt.timedelta(hours=1, minutes=14)
        t1 = run[-1]
        steps = list(pd.date_range(t0, t1, freq="2min"))
        arrs = fetch_many("PrecipRate", steps, win, workers=workers)
        stack = []
        for t in steps:
            if t not in arrs:              # missing timestep -> reset stack
                stack = []
                continue
            a = arrs[t]
            a2 = np.nan_to_num(np.clip(a, 0, None)) * 2 / 60
            i2_max = np.fmax(i2_max, a)
            total = total + a2
            stack.append(a2)
            stack = stack[-30:]            # keep the longest window (60 min)
            if len(stack) >= 8:            # i15 = mean(i16, i14), trailing 16 min
                i15 = compute_i15(stack)
                newmax = i15 > i15_max
                i15_max = np.where(newmax, i15, i15_max)
                tpki15 = np.where(newmax, t.hour + t.minute / 60.0, tpki15)
            if len(stack) >= 15:           # i30: trailing 30 min (15 steps)
                i30 = np.nansum(np.dstack(stack[-15:]), axis=2) * 60 / 30
                i30_max = np.fmax(i30_max, i30)
            if len(stack) >= 30:           # i60: trailing 60 min (30 steps)
                i60 = np.nansum(np.dstack(stack[-30:]), axis=2) * 60 / 60
                i60_max = np.fmax(i60_max, i60)

    hour_t = peak_t.replace(minute=0, second=0)
    try:
        rqi = fetch("RQI", hour_t, win)
    except Exception:                      # noqa: BLE001
        rqi = np.full(shape, np.nan, np.float32)
    try:
        shsr = fetch("SHSR", hour_t, win)
    except Exception:                      # noqa: BLE001
        shsr = np.full(shape, np.nan, np.float32)

    fields = {"i15max": i15_max, "i30max": i30_max, "i60max": i60_max,
              "i2max": i2_max, "total": total, "tpki15": tpki15,
              "rqi": rqi, "shsr": shsr}
    profile = dict(driver="GTiff", height=shape[0], width=shape[1], count=1,
                   dtype="float32", crs="EPSG:4326", transform=tr,
                   nodata=np.nan, compress="LZW")
    meta = dict(date=str(date0), peak_utc=f"{peak_t:%Y%m%d-%H%M}",
                qmax_mm=float(scout.qmax.max()), n_wet_hr=int(len(wet)),
                i15max_aoi=float(np.nanmax(i15_max)),
                i30max_aoi=float(np.nanmax(i30_max)),
                i60max_aoi=float(np.nanmax(i60_max)),
                rqi_med=float(np.nanmedian(rqi)))
    return dict(fields=fields, transform=tr, crs="EPSG:4326",
                profile=profile, meta=meta)


def multisensor_total(aoi, date=None, pad_deg=0.05, scan_pad_h=SCAN_PAD_H,
                      pass2=True, workers=WORKERS, verbose=True, window=None):
    """Gauge-corrected storm total from MRMS MultiSensor QPE (hourly).

    Sums MRMS ``MultiSensor_QPE_01H`` (Pass-2 by default -- the final,
    gauge-corrected pass) over the local-day UTC window, AOI-clipped, as a
    counterpart to the radar-only ``total`` from :func:`i15_storm_day`.
    Comparing gauges against both isolates radar-only QPE bias from the
    gauge-corrected product. MultiSensor is hourly, so there is no i15/i30/i60
    analogue -- only the accumulation.

    Returns the same dict shape as :func:`i15_storm_day` (``fields`` has the
    single key ``mstotal``), so :func:`save_fields` writes ``<key>_mstotal.tif``.
    """
    bounds, _ = load_aoi(aoi, pad_deg=pad_deg)
    if date is None and window is None:
        raise ValueError("multisensor_total needs a date or an explicit window")
    date0 = parse_date(date) if date is not None else window[0].date()
    win, tr = aoi_window(bounds)
    shape = (int(win.height), int(win.width))
    if shape[0] <= 0 or shape[1] <= 0:
        raise ValueError(f"AOI {bounds} is empty or outside the CONUS grid.")
    hours = window_hours(date0, scan_pad_h, window)
    product = "MultiSensor" if pass2 else "MultiSensor1"
    arrs = fetch_many(product, hours, win, workers=workers)
    total = np.zeros(shape, np.float32)
    n = 0
    for t in hours:                            # each hourly QPE = preceding hour
        a = arrs.get(t)
        if a is None:
            continue
        total = total + np.nan_to_num(np.clip(a, 0, None))
        n += 1
    if verbose:
        print(f"  MultiSensor {date0}: {n}/{len(hours)} hourly QPE, AOI max "
              f"{float(np.nanmax(total)):.1f} mm", flush=True)
    profile = dict(driver="GTiff", height=shape[0], width=shape[1], count=1,
                   dtype="float32", crs="EPSG:4326", transform=tr,
                   nodata=np.nan, compress="LZW")
    meta = dict(date=str(date0), n_hours=int(n), product=PRODUCTS[product][0],
                mstotal_aoi=float(np.nanmax(total)))
    return dict(fields={"mstotal": total}, transform=tr, crs="EPSG:4326",
                profile=profile, meta=meta)


def save_fields(result, out_dir, key, which=None, layout=None):
    """Write the storm-day fields to ``out_dir/<key>_<field>.tif``.

    ``which`` optionally restricts the fields written (default: all).
    Returns the list of written paths.
    """
    os.makedirs(out_dir, exist_ok=True)
    profile = {k: v for k, v in result["profile"].items()}
    paths = []
    for label, arr in result["fields"].items():
        if which and label not in which:
            continue
        path = out_path(out_dir, f"{key}_{label}.tif", layout=layout)
        with rasterio.open(path, "w", **profile) as d:
            d.write(arr.astype("float32"), 1)
        paths.append(path)
    return paths


# --------------------------------------------------------------------------- #
# virtual gauges: pull a rainfall time series from the radar at point(s)
# --------------------------------------------------------------------------- #
def _as_points(points):
    """Normalise points -> list of (name, lon, lat). Accepts a GeoDataFrame, a
    ``{name: (lon, lat)}`` dict, or a list of ``(name, lon, lat)`` / ``(lon, lat)``."""
    try:
        import geopandas as gpd
        if isinstance(points, gpd.GeoDataFrame):
            g = points.to_crs(4326)
            col = next((c for c in ("name", "station_id", "id") if c in g), None)
            names = g[col].astype(str).tolist() if col else \
                [f"VG{i + 1}" for i in range(len(g))]
            return [(n, float(p.x), float(p.y)) for n, p in zip(names, g.geometry)]
    except ImportError:                                        # pragma: no cover
        pass
    if isinstance(points, dict):
        return [(str(k), float(v[0]), float(v[1])) for k, v in points.items()]
    out = []
    for i, p in enumerate(points):
        out.append((str(p[0]), float(p[1]), float(p[2])) if len(p) == 3
                   else (f"VG{i + 1}", float(p[0]), float(p[1])))
    return out


def virtual_gauge_timeseries(points, start, end, durations=(5, 15, 30, 60),
                             multisensor=True, pad_deg=0.05, workers=WORKERS):
    """Sample MRMS at point 'virtual gauges' -> per-gauge rainfall time series.

    For each point, pulls the 2-min ``PrecipRate`` over ``[start, end]`` (UTC),
    interpolates to 1-minute, and reduces it to a cumulative total and trailing
    intensities ``i{d}`` (the 15-min one via stormscape's ``(i16+i14)/2``
    estimator, the rest plain windows). With ``multisensor`` it also samples the
    hourly gauge-corrected ``MultiSensor_QPE_01H`` (Pass-2 -> Pass-1 -> RadarOnly
    fallback) as a reference. A port of the VirtualGage timeseries in
    D. Cavagna's ``MRMS_stack``.

    Returns ``{name: DataFrame}`` indexed at 1-minute UTC with ``rate_mmph``,
    ``total_mm``, ``i{d}_mmph`` and (if multisensor) ``i60_qpe_mmph`` /
    ``total_qpe_mm`` (placed on the hour rows, NaN between).
    """
    pts = _as_points(points)
    lons = [lon for _, lon, _ in pts]
    lats = [lat for _, _, lat in pts]
    bounds = (min(lons) - pad_deg, min(lats) - pad_deg,
              max(lons) + pad_deg, max(lats) + pad_deg)
    win, tr = aoi_window(bounds)
    rc = [rasterio.transform.rowcol(tr, lon, lat) for _, lon, lat in pts]

    t2 = pd.date_range(start, end, freq="2min")
    grids = fetch_many("PrecipRate", t2, win, workers=workers)

    hrs, qhr = pd.DatetimeIndex([]), {}
    if multisensor:
        hrs = pd.date_range(pd.Timestamp(start).floor("h"), end, freq="1h")
        for prod in ("MultiSensor", "MultiSensor1", "RadarOnly"):
            todo = [h for h in hrs if h not in qhr]
            if not todo:
                break
            qhr.update(fetch_many(prod, todo, win, workers=workers))

    def _samp(a, r, c):
        return (float(a[r, c]) if a is not None and 0 <= r < a.shape[0]
                and 0 <= c < a.shape[1] else np.nan)

    t1 = pd.date_range(start, end, freq="1min")
    out = {}
    for (name, _, _), (r, c) in zip(pts, rc):
        rate1 = (pd.Series([_samp(grids.get(t), r, c) for t in t2], index=t2)
                 .reindex(t1).interpolate("time"))
        total = (rate1 / 60.0).cumsum()
        df = pd.DataFrame({"rate_mmph": rate1, "total_mm": total})
        for d in durations:
            if d == 15:
                df["i15_mmph"] = (total.diff(16) * 60.0 / 16.0
                                  + total.diff(14) * 60.0 / 14.0) / 2.0
            else:
                df[f"i{d}_mmph"] = total.diff(d) * 60.0 / d
        if multisensor and len(qhr):
            hs = pd.Series({h: _samp(qhr.get(h), r, c) for h in hrs}).sort_index()
            df["i60_qpe_mmph"] = hs.reindex(t1)
            df["total_qpe_mm"] = (hs.cumsum() - hs.iloc[0]).reindex(t1)
        out[name] = df
    return out
