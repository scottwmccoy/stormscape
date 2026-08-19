"""Sub-beam evaporation: how much of the radar-measured rain reaches the ground?

Every radar rate estimate -- Z-R, R(Kdp), R(Z,ZDR) alike -- measures rain *at
the beam*, hundreds of metres to kilometres above ground. In a dry sub-cloud
layer a substantial fraction evaporates on the way down (the limit case is
virga: everything does). Over the Great Basin in summer this is a first-order
correction, and no beam-level estimator can see it.

Model
-----
Bulk drop-evaporation integrated over an exponential DSP gives a rate loss of
the form ``dR/dz = -a (1 - RH) R^b`` with ``b ~ 0.65`` (evaporation scales
with total drop surface area, which grows sublinearly with rate -- light rain
loses proportionally more than heavy rain). Integrating analytically from the
beam down through depth ``d`` (km):

    R_ground = [ R_beam^(1-b) - (1-b) a (1-RH) d ]^(1/(1-b)),  floored at 0

Default ``a = 1.0`` is anchored so ~10 mm/h falling through 1 km of RH=20% air
loses ~35% -- the magnitude range reported for desert convection (Rosenfeld &
Mintz 1988). This is a FIRST-ORDER screen: coefficients are exposed, and the
right reading of the output is "how sensitive is this cell to the dry layer",
not a gauge-grade correction.

Moisture input
--------------
The NWS Reno radiosonde (REV / WMO 72489) launches at 00Z -- ~5 pm PDT, almost
exactly the peak hour of Great Basin monsoon convection -- and gives the RH
profile through the sub-beam layer directly. :func:`fetch_sounding` pulls the
University of Wyoming archive (network); :func:`parse_wyoming` handles a saved
listing; :func:`mean_rh` averages RH over the ground-to-beam layer. A scalar
``rh`` can also be passed straight through (e.g. from Synoptic surface obs).
"""
from __future__ import annotations

import datetime as dt
import os
from typing import Optional

import numpy as np

from .layout import out_path

#: dR/dz = -a (1-RH) R^b   (R mm/h, z km)
EVAP_A = 1.0
EVAP_B = 0.65

#: station id -> WMO number (used for saved Wyoming listings / reference)
WMO = {"REV": 72489, "OAK": 72493, "SLC": 72572, "LKN": 72582, "VEF": 72388}

# IEM RAOB JSON -- structured archive of the same launches, no HTML scraping.
# (The old weather.uwyo.edu cgi-bin endpoint was retired and 404s; a saved
# Wyoming TEXT:LIST page still parses via parse_wyoming.)
_IEM = "https://mesonet.agron.iastate.edu/json/raob.py?ts={ts}&station={stn}"


def rh_from_t_td(t_c, td_c):
    """Relative humidity (0..1) from temperature/dewpoint via Magnus."""
    t = np.asarray(t_c, dtype="float64")
    td = np.asarray(td_c, dtype="float64")
    es = lambda x: np.exp(17.625 * x / (243.04 + x))
    return np.clip(es(td) / es(t), 0.0, 1.0)


def evap_factor(rate_mmph, depth_km, rh, a: float = EVAP_A,
                b: float = EVAP_B) -> np.ndarray:
    """Fraction of the beam-level rate that survives to the ground (0..1).

    Vectorised over any mix of array/scalar inputs. RH is 0..1 over the
    sub-beam layer; ``depth_km`` the beam height AGL.
    """
    r = np.asarray(rate_mmph, dtype="float64")
    d = np.asarray(depth_km, dtype="float64")
    q = np.clip(1.0 - np.asarray(rh, dtype="float64"), 0.0, 1.0)
    with np.errstate(invalid="ignore"):
        core = np.power(np.maximum(r, 0.0), 1.0 - b) - (1.0 - b) * a * q * d
        rg = np.power(np.maximum(core, 0.0), 1.0 / (1.0 - b))
        f = np.where(r > 0, rg / np.where(r > 0, r, 1.0), 1.0)
    return np.clip(f, 0.0, 1.0)


def parse_wyoming(text: str):
    """Parse a University-of-Wyoming TEXT:LIST sounding into a DataFrame.

    Returns columns ``p_hpa, z_m (MSL), t_c, td_c, rh`` (rh 0..1). Handles the
    listing's fixed-width table between the header rule and the station block.
    """
    import pandas as pd
    rows = []
    for line in text.splitlines():
        parts = line.split()
        if len(parts) >= 5:
            try:
                p, z, t, td, rh = (float(parts[0]), float(parts[1]),
                                   float(parts[2]), float(parts[3]),
                                   float(parts[4]))
            except ValueError:
                continue
            if 100.0 <= p <= 1100.0 and -100.0 <= t <= 60.0 and 0 <= rh <= 100:
                rows.append((p, z, t, td, rh / 100.0))
    if not rows:
        raise ValueError("no sounding rows parsed; is this a TEXT:LIST page?")
    return pd.DataFrame(rows, columns=["p_hpa", "z_m", "t_c", "td_c", "rh"])


def fetch_sounding(station: str = "REV", when: Optional[dt.datetime] = None):
    """Fetch the launch nearest ``when`` (rounded to 00Z/12Z) -- network (IEM).

    ``station`` is the ICAO-less sounding id ("REV"); the IEM archive keys it
    as "K" + id. Wind-only levels (no T/Td) are dropped; RH is derived from
    T/Td via :func:`rh_from_t_td`.
    """
    import pandas as pd
    import requests
    when = when or dt.datetime.utcnow()
    day = when.date() + dt.timedelta(days=1 if when.hour >= 18 else 0)
    h = 0 if (when.hour >= 18 or when.hour < 6) else 12
    ts = f"{day:%Y%m%d}{h:02d}00"
    stn = station.upper() if station.upper().startswith("K")         else "K" + station.upper()
    r = requests.get(_IEM.format(ts=ts, stn=stn), timeout=60)
    r.raise_for_status()
    profs = r.json().get("profiles", [])
    if not profs or not profs[0].get("profile"):
        raise ValueError(f"no {stn} sounding at {ts}Z in the IEM archive")
    rows = [(lv["pres"], lv["hght"], lv["tmpc"], lv["dwpc"])
            for lv in profs[0]["profile"]
            if lv.get("tmpc") is not None and lv.get("dwpc") is not None
            and lv.get("hght") is not None]
    df = pd.DataFrame(rows, columns=["p_hpa", "z_m", "t_c", "td_c"])
    df["rh"] = rh_from_t_td(df.t_c, df.td_c)
    df.attrs["station"] = stn
    df.attrs["valid"] = profs[0].get("valid", f"{day} {h:02d}Z")
    return df


def mean_rh(profile, z0_m: float, z1_m: float) -> float:
    """Layer-mean RH (0..1) between two MSL heights, from a sounding profile."""
    z0, z1 = sorted((float(z0_m), float(z1_m)))
    d = profile[(profile.z_m >= z0 - 200) & (profile.z_m <= z1 + 200)]
    if len(d) < 2:
        d = profile.iloc[(profile.z_m - (z0 + z1) / 2).abs().argsort()[:3]]
    return float(np.interp([(z0 + z1) / 2], [0], [d.rh.mean()])[0]) \
        if False else float(d.rh.mean())


def subbeam_correct(rate_tif: str, dem_tif: str, radar_lonlat, radar_elev_m,
                    out_dir: str, key: str, rh: float,
                    tilt_deg: float = 0.5, a: float = EVAP_A,
                    b: float = EVAP_B, layout=None) -> dict:
    """Apply the sub-beam evaporation model to a rate-like raster.

    Beam height AGL per cell = 4/3-earth beam altitude at the cell's range
    (``radar_elev_m`` + r sin(tilt) + r^2/(2 (4/3) Re)) minus the DEM. Writes
    ``<key>_subbeam.tif`` (corrected field) and ``<key>_evapfrac.tif``
    (fraction LOST, 0..1). Returns a summary dict.
    """
    import pyproj
    import rasterio
    from rasterio.warp import reproject, Resampling

    with rasterio.open(rate_tif) as ds:
        rate = ds.read(1).astype("float64")
        prof, tr, crs = ds.profile, ds.transform, ds.crs
    with rasterio.open(dem_tif) as dd:
        dem = np.full(rate.shape, np.nan)
        reproject(dd.read(1), dem, src_transform=dd.transform, src_crs=dd.crs,
                  dst_transform=tr, dst_crs=crs, resampling=Resampling.average,
                  src_nodata=dd.nodata, dst_nodata=np.nan)

    rows, cols = np.indices(rate.shape)
    xs, ys = rasterio.transform.xy(tr, rows.ravel(), cols.ravel())
    t = pyproj.Transformer.from_crs(crs, 4326, always_xy=True)
    lons, lats = t.transform(np.asarray(xs), np.asarray(ys))
    rlon, rlat = radar_lonlat
    g = pyproj.Geod(ellps="WGS84")
    _, _, rng = g.inv(np.full_like(lons, rlon), np.full_like(lats, rlat),
                      lons, lats)
    rng = rng.reshape(rate.shape)
    beam_msl = (radar_elev_m + rng * np.sin(np.radians(tilt_deg))
                + rng ** 2 / (2.0 * (4.0 / 3.0) * 6371000.0))
    depth_km = np.clip((beam_msl - dem) / 1000.0, 0.0, None)

    f = evap_factor(rate, depth_km, rh, a=a, b=b)
    corrected = rate * f

    prof.update(count=1, dtype="float32", compress="LZW", nodata=np.nan)
    p_out = out_path(out_dir, f"{key}_subbeam.tif", layout)
    with rasterio.open(p_out, "w", **prof) as ds:
        ds.write(corrected.astype("float32"), 1)
        ds.update_tags(MODEL=f"dR/dz=-a(1-RH)R^b a={a} b={b}", RH=str(rh),
                       TILT_DEG=str(tilt_deg), SOURCE=os.path.basename(rate_tif))
    p_frac = out_path(out_dir, f"{key}_evapfrac.tif", layout)
    with rasterio.open(p_frac, "w", **prof) as ds:
        ds.write((1.0 - f).astype("float32"), 1)

    wet = np.isfinite(rate) & (rate > 1.0)
    return {"rh": rh, "tilt_deg": tilt_deg,
            "beam_agl_km": {"min": round(float(np.nanmin(depth_km)), 2),
                            "med": round(float(np.nanmedian(depth_km)), 2),
                            "max": round(float(np.nanmax(depth_km)), 2)},
            "loss_pct_wet_cells": {
                "median": round(float(np.nanmedian(100 * (1 - f[wet]))), 1)
                if wet.any() else None,
                "max": round(float(np.nanmax(100 * (1 - f[wet]))), 1)
                if wet.any() else None},
            "subbeam_tif": p_out, "evapfrac_tif": p_frac}
