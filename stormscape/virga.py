"""Virga-risk mask: does the mosaic's intensity have low-level support?

Motivation (Stallion Fire, Aug 2026). A field crew found little evidence of
rain over most of a burn where MRMS mapped 30+ mm/h 15-minute intensities.
The radar (KRGX) sits *inside* the perimeter; sampling its own lowest tilt
(0.0 deg, 0.1-1.4 km AGL over the burn) with two independent rate retrievals
gave ~1 mm/h at the same cells. The mosaic's intensity came from **elevated
scans** -- precipitation aloft that evaporated before reaching the ground
(virga), a routine failure mode in dry-boundary-layer convection over the
Great Basin. The reverse also occurred: shallow cells the hybrid scan
under-weighted read ~7x low against the local low tilt.

This module compares a mosaic intensity field (e.g. MRMS ``i15max``) against
a single-radar **lowest-tilt** field over the same event (e.g. the output of
``stormscape nexrad --intensity --method kdp``) and classifies each cell:

====== ============ ====================================================
value  class        meaning
====== ============ ====================================================
0      SUPPORTED    mosaic and low tilt agree within ``ratio``
1      VIRGA_RISK   mosaic >= ``ratio`` x low tilt: intensity exists
                    only aloft; suspect evaporation below the beam
2      UNDERREAD    low tilt >= ``ratio`` x mosaic: the local base scan
                    saw rain the mosaic discounted
255    NODATA       either field missing, or both below ``min_mmph``
====== ============ ====================================================

Cells where *both* fields are below ``min_mmph`` are NODATA, not SUPPORTED:
agreement about drizzle is not evidence, and flagging it would dilute the
mask. The mask is a *screen*, not a verdict -- a VIRGA_RISK cell means "check
before believing", with gauge or field evidence the arbiter.

Within ~2 km of the radar the gridded lowest tilt itself is unreliable (few
usable gates, clutter filtering), so cells inside ``exclude_km`` of the radar
are NODATA. Pass the radar location whenever it falls inside the AOI.
"""
from __future__ import annotations

import os

import numpy as np

from .layout import out_path

SUPPORTED, VIRGA_RISK, UNDERREAD, NODATA = 0, 1, 2, 255

CLASS_NAMES = {SUPPORTED: "supported", VIRGA_RISK: "virga_risk",
               UNDERREAD: "underread", NODATA: "nodata"}


def classify(mosaic, support, min_mmph: float = 10.0,
             ratio: float = 3.0) -> np.ndarray:
    """Classify aligned mosaic vs lowest-tilt intensity arrays (see module doc).

    ``ratio`` is the disagreement factor that flags a cell (default 3: the
    observed artefacts ran 4-30x while honest retrieval scatter stayed within
    ~2x). ``min_mmph`` keeps drizzle out of the mask entirely.
    """
    m = np.asarray(mosaic, dtype="float64")
    s = np.asarray(support, dtype="float64")
    if m.shape != s.shape:
        raise ValueError(f"shape mismatch {m.shape} vs {s.shape}; "
                         "regrid first (virga_mask does this for rasters)")
    out = np.full(m.shape, NODATA, dtype="uint8")
    ok = np.isfinite(m) & np.isfinite(s)
    big = ok & ((m >= min_mmph) | (s >= min_mmph))
    # eps floor so a hard zero on one side still yields a finite ratio
    eps = 0.1
    r = np.where(big, m / np.maximum(s, eps), np.nan)
    out[big & (r >= ratio)] = VIRGA_RISK
    out[big & (np.maximum(m, eps) <= s / ratio)] = UNDERREAD
    out[big & (out == NODATA)] = SUPPORTED
    return out


def virga_mask(mosaic_tif: str, support_tif: str, out_dir: str, key: str,
               min_mmph: float = 10.0, ratio: float = 3.0,
               radar_lonlat=None, exclude_km: float = 2.0,
               layout=None) -> dict:
    """Raster front-end: regrid ``mosaic_tif`` onto ``support_tif``'s grid,
    classify, and write ``<key>_virgarisk.tif`` (uint8 classes) plus
    ``<key>_supportratio.tif`` (mosaic / low-tilt, float).

    The *support* raster (single-radar lowest tilt) defines the output grid --
    it is the finer, local product. Returns a summary dict (cell counts, %,
    paths).
    """
    import rasterio
    from rasterio.warp import reproject, Resampling

    with rasterio.open(support_tif) as ds:
        sup = ds.read(1).astype("float64")
        prof = ds.profile
        tr, crs = ds.transform, ds.crs
    with rasterio.open(mosaic_tif) as dm:
        mos = np.full(sup.shape, np.nan, dtype="float64")
        reproject(dm.read(1), mos, src_transform=dm.transform, src_crs=dm.crs,
                  dst_transform=tr, dst_crs=crs,
                  resampling=Resampling.bilinear,
                  src_nodata=dm.nodata, dst_nodata=np.nan)

    cls = classify(mos, sup, min_mmph=min_mmph, ratio=ratio)

    if radar_lonlat is not None and exclude_km > 0:
        import pyproj
        rows, cols = np.indices(sup.shape)
        xs, ys = rasterio.transform.xy(tr, rows.ravel(), cols.ravel())
        t = pyproj.Transformer.from_crs(crs, 4326, always_xy=True)
        lons, lats = t.transform(np.asarray(xs), np.asarray(ys))
        g = pyproj.Geod(ellps="WGS84")
        rlon, rlat = radar_lonlat
        _, _, rng = g.inv(np.full_like(lons, rlon), np.full_like(lats, rlat),
                          lons, lats)
        cls[(rng.reshape(sup.shape) < exclude_km * 1000.0)] = NODATA

    with np.errstate(divide="ignore", invalid="ignore"):
        rat = (mos / np.maximum(sup, 0.1)).astype("float32")

    prof.update(count=1, compress="LZW")
    p_cls = out_path(out_dir, f"{key}_virgarisk.tif", layout)
    with rasterio.open(p_cls, "w", **{**prof, "dtype": "uint8",
                                      "nodata": NODATA}) as ds:
        ds.write(cls, 1)
        ds.update_tags(CLASSES="0=supported,1=virga_risk,2=underread,255=nodata",
                       MIN_MMPH=str(min_mmph), RATIO=str(ratio),
                       MOSAIC=os.path.basename(mosaic_tif),
                       SUPPORT=os.path.basename(support_tif))
    p_rat = out_path(out_dir, f"{key}_supportratio.tif", layout)
    with rasterio.open(p_rat, "w", **{**prof, "dtype": "float32",
                                      "nodata": np.nan}) as ds:
        ds.write(rat, 1)

    n = {c: int((cls == v).sum()) for v, c in CLASS_NAMES.items()}
    assessed = n["supported"] + n["virga_risk"] + n["underread"]
    pct = {c: (100.0 * n[c] / assessed if assessed else 0.0)
           for c in ("supported", "virga_risk", "underread")}
    return {"counts": n, "percent": {k: round(v, 1) for k, v in pct.items()},
            "assessed": assessed, "virgarisk_tif": p_cls,
            "supportratio_tif": p_rat}
