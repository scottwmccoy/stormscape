"""Validation study: does BRISK compute dNBR the same way the BAER teams do?

BRISK gives a burn scar within a day; the BAER teams' own dNBR arrives days to
weeks later but is the product post-fire hazard work is built on. Both are
published in the same SSEC archive, so for every fire that has both we can put
them cell-for-cell on one grid and check.

Run it::

    python examples/brisk_vs_baer.py                    # all matched fires
    python examples/brisk_vs_baer.py --years 2025 --out-dir ./baer_study

It downloads nothing permanently -- both products are read straight over HTTP
with GDAL range requests -- and writes a per-fire CSV plus a two-panel figure.

Method
------
Four choices matter, and each is a place a careless comparison goes wrong:

1. **Scaling.** BAER dNBR is int16 **x1000** (the BARC convention; NOAA's own
   ``BARC256 = dNBR*5 - 275`` identity confirms it). BRISK is float32 unscaled.
   Divide BAER by 1000 or the correlation is fine and every ratio is 1000x off.

2. **Support.** BAER is 20 m (Sentinel-2 source) or 30 m (Landsat) in Albers;
   BRISK is 60 m Web-Mercator (~46 m of ground). BAER is **averaged down onto
   the BRISK grid**, aggregating the finer product to the coarser support --
   upsampling BRISK instead would invent detail in the thing under test.

3. **Which cells.** Scoring every valid cell inflates r by ~0.03-0.10, because
   the large unburned surround agrees near zero and does the work. Only
   **burned cells (BAER dNBR >= 0.1)** answer the question.

4. **Which date.** BRISK is a *running composite*, so it is scored twice: on the
   BAER assessment date, and at a fixed **+14 days**. The offset is set in
   advance rather than tuned per fire -- picking each fire's best-correlating
   scene would fit noise and inflate the answer.

What it found (39 fires, 2025)
------------------------------
BRISK computes the same quantity: median **slope 1.013, bias +0.003, ratio
1.011** over burned cells, and **r = 0.938 across 1.24 M cells** pooled over
well-matched fires. The poor performers are **compositing latency, not a
different algorithm** -- every one read *low* on the BAER date and recovered to
r = 0.80-0.96 given a scene 5-21 days later. The +14 d rule takes fires at
r >= 0.90 from 14/39 to 24/39 and fires below 0.60 from 5 to 1.

So: trust BRISK's *pattern* immediately, but give the composite about two weeks
before trusting its *magnitude* -- which is what ``burn --min-age`` enforces.
"""

from __future__ import annotations

import argparse
import datetime as dt
import os
import warnings
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import reproject

from stormscape.burn import catalog

warnings.filterwarnings("ignore")

MATURED_OFFSET_DAYS = 14        # fixed in advance; see the module docstring
BURNED_MIN_DNBR = 0.10          # BAER's own unburned break
MIN_CELLS = 100
GDAL_ENV = dict(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
                CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",
                GDAL_HTTP_MAX_RETRY="2", GDAL_HTTP_RETRY_DELAY="1")


def _matched_fires(years):
    """``{fire: (baer_rows, brisk_rows)}`` for fires carrying both products."""
    baer = catalog("baer_dnbr", years=years, verbose=False)
    brisk = catalog("dnbr", years=years, verbose=False)
    if not len(baer) or not len(brisk):
        return {}
    b = {f: g for f, g in baer.groupby("fire")}
    k = {f: g for f, g in brisk.groupby("fire")}
    return {f: (b[f], k[f]) for f in sorted(set(b) & set(k))}


def _read_pair(baer_url, brisk_url):
    """BAER dNBR resampled onto the BRISK grid; returns ``(brisk, baer)`` arrays."""
    with rasterio.Env(**GDAL_ENV):
        with rasterio.open(brisk_url) as kds:
            k = kds.read(1).astype("float64")
            dst_t, dst_crs = kds.transform, kds.crs
        with rasterio.open(baer_url) as bds:
            src = bds.read(1).astype("float64")
            if bds.nodata is not None and np.isfinite(bds.nodata):
                src = np.where(src == bds.nodata, np.nan, src)
            src = src / 1000.0                       # x1000 -> plain dNBR
            b = np.full(k.shape, np.nan)
            reproject(source=src, destination=b,
                      src_transform=bds.transform, src_crs=bds.crs,
                      dst_transform=dst_t, dst_crs=dst_crs,
                      src_nodata=np.nan, dst_nodata=np.nan,
                      resampling=Resampling.average)  # fine -> coarse support
            res_m = abs(bds.res[0])
    return k, b, res_m


def _stats(brisk, baer):
    m = np.isfinite(brisk) & np.isfinite(baer) & (baer >= BURNED_MIN_DNBR)
    if m.sum() < MIN_CELLS:
        return None
    x, y = brisk[m], baer[m]
    if np.std(x) == 0 or np.std(y) == 0:
        return None
    slope, intercept = np.polyfit(y, x, 1)           # BRISK = f(BAER)
    d = x - y
    return dict(n=int(m.sum()), r=float(np.corrcoef(x, y)[0, 1]),
                slope=float(slope), intercept=float(intercept),
                bias=float(d.mean()), rmse=float(np.sqrt((d ** 2).mean())),
                ratio=float(np.median(x / y)))


def compare_fire(fire, baer_rows, brisk_rows):
    """Score one fire on the BAER date and at +14 days."""
    br = baer_rows.sort_values("date").iloc[-1]      # latest BAER assessment
    cand = brisk_rows.sort_values("date")
    out = dict(fire=fire, baer_date=str(br.date))
    for tag, off in (("same", 0), ("d14", MATURED_OFFSET_DAYS)):
        want = br.date + dt.timedelta(days=off)
        pick = cand.iloc[int(np.argmin([abs((d - want).days) for d in cand.date]))]
        try:
            k, b, res_m = _read_pair(br.url, pick.url)
        except Exception as exc:                                  # noqa: BLE001
            print(f"  !! {fire} ({tag}): {exc}", flush=True)
            continue
        out[f"gap_{tag}"] = (pick.date - br.date).days
        out.setdefault("baer_res_m", res_m)
        st = _stats(k, b)
        if st:
            out.update({f"{key}_{tag}": v for key, v in st.items()})
    return out


def summarise(df):
    for tag, name in (("same", "BRISK on the BAER date"),
                      ("d14", f"BRISK +{MATURED_OFFSET_DAYS} days (matured)")):
        r = df.get(f"r_{tag}", pd.Series(dtype=float)).dropna()
        if not len(r):
            continue
        print(f"\n{name}:  n={len(r)} fires")
        print(f"   median r      {r.median():.3f}   "
              f"IQR {r.quantile(.25):.3f}..{r.quantile(.75):.3f}")
        print(f"   median slope  {df[f'slope_{tag}'].median():.3f}   "
              f"bias {df[f'bias_{tag}'].median():+.4f}   "
              f"ratio {df[f'ratio_{tag}'].median():.3f}")
        print(f"   r>=0.90: {(r >= .9).sum()}/{len(r)}   "
              f"r>=0.80: {(r >= .8).sum()}/{len(r)}   r<0.60: {(r < .6).sum()}/{len(r)}")
    if "baer_res_m" in df:
        print("\nby BAER source resolution (20 m = Sentinel-2, 30 m = Landsat):")
        for res, g in df.groupby("baer_res_m"):
            rr = g["r_same"].dropna()
            if len(rr):
                print(f"   {int(res)} m  n={len(g):2d}  median r {rr.median():.3f}")


def plot(df, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    d = df.dropna(subset=["r_same"]).sort_values("r_same")
    fig, ax = plt.subplots(figsize=(7.5, 0.22 * len(d) + 2.2))
    y = np.arange(len(d))
    ax.hlines(y, d.r_same, d.r_d14, color="0.75", lw=1.2, zorder=1)
    ax.scatter(d.r_same, y, s=26, color="#d62728",
               label="BRISK on the BAER date", zorder=2)
    ax.scatter(d.r_d14, y, s=26, color="#2ca02c",
               label=f"BRISK +{MATURED_OFFSET_DAYS} days", zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels(d.fire, fontsize=6.5)
    ax.axvline(0.9, color="0.5", lw=0.7, ls=":")
    ax.set_xlabel("correlation with BAER dNBR  (burned cells)")
    ax.set_title("BRISK vs BAER dNBR per fire:\nthe failures are compositing "
                 "latency, not a different algorithm", fontsize=11)
    ax.legend(loc="upper left", fontsize=9, framealpha=0.95)
    ax.set_xlim(-0.05, 1.02)
    ax.set_ylim(-1, len(d))
    ax.grid(axis="x", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_png, dpi=170, bbox_inches="tight")
    print(f"\nwrote {out_png}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--years", type=int, nargs="+", default=[2025],
                    help="archive years to search (default 2025 -- the only "
                         "year with BAER dNBR at the time of writing)")
    ap.add_argument("--out-dir", default="baer_study")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args(argv)

    os.makedirs(args.out_dir, exist_ok=True)
    matched = _matched_fires(args.years)
    print(f"{len(matched)} fires carry both a BRISK composite and a BAER dNBR")
    if not matched:
        return None

    with ThreadPoolExecutor(args.workers) as ex:
        rows = list(ex.map(lambda kv: compare_fire(kv[0], *kv[1]),
                           matched.items()))
    df = pd.DataFrame([r for r in rows if r])
    csv = os.path.join(args.out_dir, "brisk_vs_baer.csv")
    df.to_csv(csv, index=False)
    print(f"wrote {csv}")
    summarise(df)
    if "r_d14" in df:
        plot(df, os.path.join(args.out_dir, "brisk_vs_baer.png"))
    return df


if __name__ == "__main__":
    main()
