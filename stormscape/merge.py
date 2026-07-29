"""Gauge-radar merging: adjust a radar field toward gauge point values.

Two methods, plus leave-one-out validation:

* :func:`mean_field_bias` -- the classic single multiplicative factor
  (Sigma gauge / Sigma radar-at-gauges). Removes the *average* bias but adds no
  spatial structure (the baseline you've used before).
* :func:`conditional_merge` -- conditional merging after Sinclair & Pegram
  (2005), done in log space for the multiplicative rainfall bias: keep the
  radar's small-scale texture but swap its smooth component for the
  gauge-interpolated field,
  ``merged = expm1( log1p(R) - IDW(log1p(R@gauges)) + IDW(log1p(gauges)) )``.

Interpolation is inverse-distance weighting in a local-metre frame (no extra
dependencies). :func:`loo_cross_validate` scores raw radar vs MFB vs CM by
predicting each held-out gauge from the others.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from rasterio.transform import rowcol, xy


def _to_m(lon, lat, lon0, lat0):
    """Local equirectangular metres about (lon0, lat0)."""
    lon, lat = np.asarray(lon, float), np.asarray(lat, float)
    return ((lon - lon0) * 111320.0 * np.cos(np.radians(lat0)),
            (lat - lat0) * 111320.0)


def sample_field(field, transform, lons, lats):
    """Value of ``field`` at each (lon, lat); off-grid -> NaN."""
    H, W = field.shape
    out = []
    for lo, la in zip(np.atleast_1d(lons), np.atleast_1d(lats)):
        r, c = rowcol(transform, lo, la)
        r, c = int(r), int(c)
        out.append(field[r, c] if 0 <= r < H and 0 <= c < W else np.nan)
    return np.asarray(out, float)


def _idw(sx, sy, sv, dx, dy, power=2.0, eps=1.0):
    """Inverse-distance interpolation of (sx,sy,sv) onto (dx,dy) -> 1-D array."""
    m = np.isfinite(sv)
    sx, sy, sv = sx[m], sy[m], sv[m]
    if len(sv) == 0:
        return np.full(np.size(dx), np.nan)
    d = np.sqrt((dx[:, None] - sx[None, :]) ** 2 + (dy[:, None] - sy[None, :]) ** 2)
    w = 1.0 / np.maximum(d, eps) ** power
    return (w * sv[None, :]).sum(1) / w.sum(1)


def _grid_m(field, transform, lon0, lat0):
    H, W = field.shape
    rr, cc = np.meshgrid(np.arange(H), np.arange(W), indexing="ij")
    xs, ys = xy(transform, rr.ravel(), cc.ravel())   # cell centres (lon, lat)
    return _to_m(np.asarray(xs), np.asarray(ys), lon0, lat0)


def mean_field_bias(field, transform, lons, lats, vals):
    """Multiplicative mean-field-bias correction. Returns (field*, factor)."""
    r = sample_field(field, transform, lons, lats)
    v = np.asarray(vals, float)
    m = np.isfinite(r) & np.isfinite(v) & (r > 0)
    factor = float(v[m].sum() / r[m].sum()) if m.any() and r[m].sum() > 0 else 1.0
    return np.clip(field * factor, 0, None), factor


def local_bias(field, transform, lons, lats, vals, power=2.0, clip_log=2.5):
    """Spatially-variable multiplicative bias: ``radar * exp(IDW(log gauge/radar))``.

    Keeps the radar's spatial detail (it only rescales it, like MFB) but lets the
    correction factor vary in space. Uses only wet gauges (radar>0, gauge>0);
    the log-factor is clipped to ``+/-clip_log`` to tame outliers.
    """
    r = sample_field(field, transform, lons, lats)
    v = np.asarray(vals, float)
    m = np.isfinite(r) & np.isfinite(v) & (r > 0) & (v > 0)
    if m.sum() < 3:
        return mean_field_bias(field, transform, lons, lats, vals)[0]
    lon0, lat0 = float(np.nanmean(lons)), float(np.nanmean(lats))
    logf = np.clip(np.log(v[m] / r[m]), -clip_log, clip_log)
    sxm, sym = _to_m(np.asarray(lons)[m], np.asarray(lats)[m], lon0, lat0)
    gxm, gym = _grid_m(field, transform, lon0, lat0)
    lf = _idw(sxm, sym, logf, gxm, gym, power).reshape(field.shape)
    merged = np.clip(field, 0, None) * np.exp(lf)
    return np.where(np.isfinite(field), np.clip(merged, 0, None), np.nan)


def conditional_merge(field, transform, lons, lats, vals, power=2.0):
    """Conditional merging (log space) of gauges into a radar field."""
    lons, lats, vals = np.asarray(lons, float), np.asarray(lats, float), np.asarray(vals, float)
    lon0, lat0 = float(np.nanmean(lons)), float(np.nanmean(lats))
    gxm, gym = _grid_m(field, transform, lon0, lat0)
    sxm, sym = _to_m(lons, lats, lon0, lat0)
    r_at_g = sample_field(field, transform, lons, lats)
    lg = np.log1p(np.clip(vals, 0, None))
    lr = np.log1p(np.clip(r_at_g, 0, None))
    G = _idw(sxm, sym, lg, gxm, gym, power).reshape(field.shape)
    Rg = _idw(sxm, sym, lr, gxm, gym, power).reshape(field.shape)
    merged = np.expm1(np.log1p(np.clip(field, 0, None)) - Rg + G)
    return np.where(np.isfinite(field), np.clip(merged, 0, None), np.nan)


def loo_cross_validate(field, transform, lons, lats, vals, power=2.0):
    """Leave-one-out: predict each gauge from the others. Returns (table, stats).

    ``stats`` has RMSE / MAE / bias / mass-ratio for ``raw`` (radar), ``mfb``
    and ``cm`` against the held-out gauge values.
    """
    lons, lats, vals = np.asarray(lons, float), np.asarray(lats, float), np.asarray(vals, float)
    r_all = sample_field(field, transform, lons, lats)
    idx = np.where(np.isfinite(r_all) & np.isfinite(vals))[0]
    lon0, lat0 = float(np.nanmean(lons[idx])), float(np.nanmean(lats[idx]))
    rows = []
    for i in idx:
        o = idx[idx != i]
        ro, vo = r_all[o], vals[o]
        f = float(vo.sum() / ro.sum()) if ro.sum() > 0 else 1.0
        sxm, sym = _to_m(lons[o], lats[o], lon0, lat0)
        xi, yi = _to_m(lons[i:i + 1], lats[i:i + 1], lon0, lat0)
        Gi = _idw(sxm, sym, np.log1p(np.clip(vo, 0, None)), xi, yi, power)[0]
        Rgi = _idw(sxm, sym, np.log1p(np.clip(ro, 0, None)), xi, yi, power)[0]
        cm = np.expm1(np.log1p(max(r_all[i], 0.0)) - Rgi + Gi)
        mw = (ro > 0) & (vo > 0)
        if mw.sum() >= 3:
            lf = np.clip(np.log(vo[mw] / ro[mw]), -2.5, 2.5)
            sxw, syw = _to_m(lons[o][mw], lats[o][mw], lon0, lat0)
            lbc = max(r_all[i], 0.0) * np.exp(_idw(sxw, syw, lf, xi, yi, power)[0])
        else:
            lbc = r_all[i] * f
        rows.append(dict(gauge=vals[i], raw=r_all[i], mfb=r_all[i] * f,
                         lbc=lbc, cm=max(cm, 0.0)))
    d = pd.DataFrame(rows)

    def stat(col):
        e = d[col] - d.gauge
        return dict(rmse=float(np.sqrt((e ** 2).mean())), mae=float(e.abs().mean()),
                    bias=float(e.mean()),
                    ratio=float(d[col].sum() / d.gauge.sum()) if d.gauge.sum() else np.nan)

    return d, {k: stat(k) for k in ("raw", "mfb", "lbc", "cm")}
