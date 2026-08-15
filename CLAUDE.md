# Project: stormscape — terrain + radar-rainfall mapping

A small, reusable Python package that turns an **area of interest** (a bbox, a
vector file, or a shapely geometry) plus a **storm-day date** into:

1. a USGS **3DEP / The National Map** DEM + hillshade (`stormscape/dem.py`);
2. **peak 15/30/60-minute rainfall-intensity (`i15`/`i30`/`i60`) fields** stacked
   from 2-minute NOAA **MRMS** `PrecipRate` returns (+ a gauge-corrected
   MultiSensor QPE total), with companion fields — storm total, peak 2-min rate,
   time-of-peak, radar quality (RQI), beam height (`stormscape/mrms.py`);
3. **ground rain-gauge rainfall** from the **Synoptic / MesoWest** API, reduced
   to the same metrics — storm total + peak 15/30/60-min intensity
   (`stormscape/gauges.py`);
4. a **figure** draping the radar field over the hillshade with optional
   labelled vector overlays (NHD streams, TIGER roads, GNIS places) and gauge
   points, plus a **radar-vs-gauge comparison** — per-gauge residuals + skill
   stats (`stormscape/plot.py`, `stormscape/compare.py`, `stormscape/refdata.py`);
5. **near-real-time burn severity** — CIMSS **BRISK** daily multi-satellite
   dNBR composites (+ the BAER soil-burn-severity archive) over the AOI, cached,
   classified and mapped (`stormscape/burn.py`);
6. **abandoned mine features** — USGS **USMIN** topo mine symbols (dumps,
   tailings, adits, shafts) over the AOI, grouped, filterable and mappable as
   points or a per-km2 density surface (`stormscape/mines.py`);
7. raw **single-radar NEXRAD Level II** tilts (reflectivity/velocity) for the
   radar nearest the AOI — gridded over the AOI or sampled at the gauges, with a
   Z–R diagnostic (`stormscape/nexrad.py`); the underlying radar behind the MRMS
   mosaic, and the way to reach pre-2020 events.

Driven entirely by AOI + date; no project-specific inputs. It was lifted out of
a post-fire debris-flow study and generalized — nothing here depends on that
project's fire perimeters or event inventory.

## Python environment — IMPORTANT

> **Note for anyone who cloned this repo:** `GISMan` below is the *maintainer's*
> local conda env name. Substitute your own env — or build a fresh one with
> `conda env create -f environment.yml` (creates `stormscape`), which installs
> every dependency named here. The rules that matter generally are: conda-forge
> with strict channel priority, prefer conda over pip, never install into `base`.

- **On the maintainer's machine, run all Python with the GISMan conda env:**
  `/opt/anaconda3/envs/GISMan/bin/python`. It already satisfies every
  dependency — py3dep, rioxarray, **rasterio built with the GRIB driver**,
  geopandas, shapely, contextily, matplotlib.
- GISMan is a **conda-forge (strict)** env. Add packages with
  `conda install -n GISMan -c conda-forge --override-channels <pkg>`. Only fall
  back to GISMan's pip if a package is unavailable on conda-forge — and say so.
- **Never install into the base env.** Keep base clean. Prefer conda over pip.
- **NEXRAD Level II** (`stormscape/nexrad.py`) adds `arm_pyart` + `boto3` +
  `pytz` + `scipy` (conda-forge) and **`nexradaws` via GISMan's pip** — nexradaws
  is *not on conda-forge*, so pip is the only route. All installed in GISMan
  (2026-06-25). The geo stack (`pyproj`/`rioxarray`/`scipy`) does the gridding.
- **Portable / for sharing:** `conda env create -f environment.yml` builds a
  standalone `stormscape` env from conda-forge; then `pip install -e .` installs the
  package + the `stormscape` CLI into whichever env is active.
- MRMS needs rasterio's **GRIB** driver (`libgdal-grib` on conda-forge). Verify:
  `python -c "import rasterio; assert 'GRIB' in rasterio.drivers.raster_driver_extensions().values()"`
- **GeoPDF export** (`stormscape/export.py`) needs GDAL's **PDF** driver, a
  separate conda-forge plugin (`libgdal-pdf`, installed in GISMan 2026-06-29;
  also pulled `poppler` for read-back). conda-forge's `libgdal-core` ships
  *without* PDF — verify `from osgeo import gdal; gdal.GetDriverByName("PDF")` is
  not None (`DCAP_CREATECOPY=YES`). Imported lazily, so `import stormscape` + the
  3857 GeoTIFF export work without it; only the GeoPDF path requires it.

## Running it

```bash
# whole pipeline: DEM -> hillshade -> i15 -> figure (+ gauges + comparison)
python -m stormscape run --aoi aoi.kmz --date 20260619 --resolution 10 \
    --out-dir ./out --key event --perimeters aoi.kmz --reference --clip \
    --gauges --compare --multisensor       # gauge steps need $SYNOPTIC_TOKEN
```

```bash
# raw single-radar NEXRAD Level II: storm-peak reflectivity over the nearest radar
python -m stormscape nexrad --aoi aoi.kmz --composite \
    --start 202606192000 --end 202606200200 \
    --hillshade ./out/event_hillshade.tif --out-dir ./out --key event
```

`run` chains `dem`, `i15`, and `map` (each usable on its own); `gauges`,
`compare`, and `nexrad` are also standalone subcommands. The AOI is `--aoi <vector file>`
(KMZ/GeoJSON/SHP/GPKG) or `--bbox W S E N` (lon/lat degrees). `--reference`
overlays labelled vectors; `--clip` crops to the AOI; `--basemap` underlays USGS
topo tiles (needs contextily). `--gauges` fetches + overlays Synoptic gauges;
`--compare` adds radar-vs-gauge residuals + skill stats + a residual map;
`--multisensor` adds the gauge-corrected MultiSensor total. Full flag table is
in `README.md`.

## Layout

- `stormscape/aoi.py` — AOI parsing (bbox / vector / geometry) + overlay loading
- `stormscape/dem.py` — 3DEP DEM download, 1 m availability check, hillshade
- `stormscape/mrms.py` — MRMS fetch/stack → i15/i30/i60 fields + MultiSensor total
- `stormscape/gauges.py` — Synoptic/MesoWest gauges → total + peak intensities
- `stormscape/compare.py` — sample radar at gauges → residuals + skill stats
- `stormscape/nexrad.py` — single-radar NEXRAD Level II tilts (nexradaws + Py-ART)
- `stormscape/atlas14.py` — NOAA Atlas 14 gridded climatology → intensity fields + anomaly
- `stormscape/smoothing.py` — NaN-aware field smoothing (gaussian/uniform/median/idw) + radar-gauge skill sweep
- `stormscape/export.py` — georeferenced exports for GIS/CalTopo: EPSG:3857 GeoTIFFs (raw float + colorized RGBA) + GeoPDF figures
- `stormscape/burn.py` — near-real-time burn severity (BRISK dNBR / BAER SBS)
- `stormscape/mines.py` — abandoned mine features (USGS USMIN) + density
- `stormscape/refdata.py` — AOI-scoped NHD streams / TIGER roads / GNIS places
- `stormscape/plot.py` — drape i15 over hillshade + basemap/vector/gauge overlays
- `stormscape/data/` — bundled tables (`nexrad_sites.csv` NCEI HOMR, `atlas14_regions.csv`)
- `stormscape/cli.py` — `stormscape {dem,i15,map,run,gauges,compare,nexrad,panels,vgauge,zoom,pick,climate,smooth,recurrence,export,burn,mines}`
- `environment.yml`, `pyproject.toml`, `examples/`

## Conventions & gotchas (learned during development)

- **i15 estimator** (after D. Cavagna's `MRMS_stack.py`): MRMS `PrecipRate` is a
  2-min instantaneous rate (mm/h); per step `a2 = rate·2/60` (mm). Over a
  trailing 16-min window (8 steps) `i16 = Σ(8)·60/16`, `i14 = Σ(last 7)·60/14`,
  and `i15 = mean(i16, i14)`; keep the per-cell running max over the storm.
- **i30/i60:** plain trailing windows (15 / 30 two-min steps `·60/N`) over the
  same stack, per-cell running max — so `i15max ≥ i30max ≥ i60max` cell-wise.
  The gauge-side i15 reuses the radar `(i16+i14)/2` estimator on the 1-min
  interpolated series so the two sides are directly comparable.
- **MRMS source:** NOAA S3 `noaa-mrms-pds.s3.amazonaws.com/CONUS/<PROD>/...`;
  CONUS grid UL(−130, 55), 0.01°, 7000×3500. A 404 means "missing" — reset the
  stack, never fabricate continuity. The 2-min cadence is reliable from ~Oct 2020.
- **Always inspect the `rqi` field.** Far from radar / behind terrain the beam
  overshoots low rain and i15 is unreliable; filter by RQI for quantitative work.
- **CRS:** DEMs are stored in **EPSG:5070** (CONUS Albers, equal-area m); MRMS
  GeoTIFFs are native **EPSG:4326**. **Figures display in an auto-selected UTM
  zone** (near north-up, with lat/long axis ticks via `scale_ticks`) — Albers is
  too convergent this far west (~15° at −119°) for straight lat/long ticks or an
  up-pointing north arrow. Set `drape_i15(work_crs="UTM")` (CLI figures already do).
- **Reference streams:** `add_reference` draws all watercourses but splits them
  — *named* mainstems at normal weight + labels, *unnamed* headwaters as finer,
  lighter lines (no labels) — so the dense NHD HR network adds drainage texture
  without burying the i15 field. River-name labels use the place-name font size.
- **Colormap:** the i15 field defaults to **`YlGnBu`** (colourblind-safe
  sequential); the residual map colours gauges by radar−gauge on diverging
  `RdBu_r`. `--cmap` overrides the field (incl. `cmc.*` Crameri maps, registered
  lazily). Never use jet/turbo — perceptually non-uniform (per scientific-plots).
- **GNIS layer 1 (incorporated places) returns polygons,** not points — use
  `geometry.representative_point()` for placement, not `.geometry.x/.y`.
- **contextily `add_basemap`** takes imshow kwargs (e.g. `zorder=0`) directly,
  **not** an `extra_imshow_args` dict (that errors on contextily ≥ 1.7).
- **Synoptic gauges:** the **Time Series API** with `precip=1` (the *Basic
  Precipitation Service*) → `precip_intervals_set_1d` (incremental mm), which is
  `gauges.INT_VAR`. Needs a free token (`$SYNOPTIC_TOKEN`; academic open-access
  is free). The dedicated `/stations/precipitation` endpoint is Enterprise-gated
  (403). Gauge transport + 1-min interpolation are vendored from USGS
  **FlowAlert** (CC0). Stations reporting no precip variable → NaN metrics
  (handle, don't crash).
- **Gauge cadence (`report_min`) must be measured on *precip-bearing* obs.** A
  multi-variable ASOS (e.g. KRNO) reports temp/wind every 5 min but precip only
  at hourly METAR times; counting all obs would mislabel it 5-min and wrongly
  pass the `--max-report-min` screen. Coarse reporters smear i15/i30 low under
  the 1-min interpolation; the storm total is cadence-insensitive.
- **MultiSensor QPE** (`MultiSensor_QPE_01H_Pass2`, hourly, gauge-corrected) is
  the comparison baseline that isolates radar-only bias. On Hidden Valley
  20260619 radar-only ran ~2.6× the gauges on storm total; MultiSensor ~2.0× —
  gauge correction roughly halves it, but point-vs-pixel leaves a residual gap,
  and cadence was *not* the driver of the large i15 gap (gauges already ≤15-min).
- **NEXRAD Level II** (`nexrad.py`): transport is **nexradaws ≥2.0**, which points
  at the current **`unidata-nexrad-level2`** S3 bucket — the older
  `noaa-nexrad-level2` Big-Data bucket was deprecated (that's why hand-rolled S3
  access 404'd earlier). Volumes read with **Py-ART**; both are *optional* deps
  (`pip install -e ".[nexrad]"`), imported lazily so `import stormscape` works
  without them. `nearest_radar` uses the bundled NCEI HOMR table
  (`data/nexrad_sites.csv`, 163 WSR-88D). The lowest sweep is gridded by
  nearest-gate ground projection (AEQD about the radar → reproject 4326), masked
  beyond ~1° beam spacing; output is an `mrms`-style result dict (so `save_fields`
  / `drape_i15` just work). Sweep indices are pyart's **split-cut order** (pairs
  per elevation: reflectivity cut, then Doppler cut), so `--field velocity
  --sweep 0` may be empty — use the Doppler sweep. **KRGX is a mountaintop radar
  (2559 m): its lowest tilt is ~0.0°**, not 0.5° (`fixed_angle [0,0,0.48,...]`) —
  an unusually clean low look, reinforcing that HV's i15 over-read is *not* beam
  overshoot. Archive reaches back to the ~1990s — *further than* MRMS's 2020 2-min
  cadence.
- **Z–R diagnostic** (`nexrad.z_to_rate`): `R=(z/a)^(1/b)`, default WSR-88D
  convective (a=300, b=1.4); `dbz_cap=53` applies the operational hail cap. The
  gap between capped and uncapped rate at high-dBZ cells *is* the hail
  over-estimation. **HV 20260619 live (KRGX) — point #7 closed:** lowest-tilt peak
  reflectivity hit **70.5 dBZ** in the AOI and **≥55 dBZ (hail) at 7/28 wet
  gauges** (peak 60–65 dBZ at the wettest), so MRMS's convective+hail `PrecipFlag`
  typing is *justified by the raw single-radar Z*. Uncapped convective Z–R gives
  350–750 mm/h at those cells vs ~104 capped → the i15 over-read is genuine radar
  QPE physics, not artefact. (Deliverables in the Box event folder:
  `HiddenValley_reflmax.tif`, `_nexrad.png`; `nexrad_cache/` holds the volumes —
  keep out of git.)
- **Single-radar i15 stack** (`nexrad.intensity_stack`, **v1**, CLI `nexrad
  --intensity`): Level II → per-low-tilt-sweep **capped convective Z–R** rate
  (harvests *all* SAILS low cuts, deduping the split-cut twin → ~3-min effective
  cadence; HV got 2.93 min from 62 volumes / 124 cuts) → 1-min interp → the same
  i15=(i16+i14)/2 / i30 / i60 running-max estimators as `mrms` → an mrms-style
  result dict (`<key>_i15max.tif`, `_i30max`, `_i60max`, `_total_mm`,
  `_peakrate_mmph`), so it flows straight into `compare`/`plot`/`merge`. No-echo
  cells = 0 mm/h (assumes AOI in radar coverage). **HV 20260619 validation vs
  MRMS:** L2 i15max 103.8 (MRMS 102.2); cell-wise L2/MRMS ≈ **1.56** (r 0.76),
  i30/i60 ≈ 1.6; at gauges median MRMS/gauge 2.30, **L2/gauge 3.97**, L2/MRMS
  1.50, corr(L2,gauge) 0.80 vs (MRMS,gauge) 0.89. So v1 tracks MRMS spatially but
  over-reads ~1.5× — expected for one fixed convective Z–R vs MRMS's
  dual-pol/VPR/multi-radar QC; **cadence is NOT the cause** (~3 min, near MRMS's
  2). Fig `HiddenValley_mrms_vs_nexrad.png`.
- **v2 dual-pol `--method kdp` (done):** R(Kdp) where Z≥`z_blend` (35 dBZ), Kdp
  from `pyart.retrieve.kdp_maesaka` (variational; clip 0–7°/km, use Kdp>0 only),
  S-band R=44·Kdp^0.822, blended to capped Z–R in light rain + pre-dual-pol →
  cross-era. **Closes the gap:** HV i15 **1.56→0.80× MRMS** (i30 0.83, i60 0.85),
  **2.0× gauges** (MRMS 2.30×, v1 3.97×), v2/MRMS 0.86 at gauges; less biased but
  noisier (at-gauge corr 0.68 vs MRMS 0.89 — Kdp texture; raise `z_blend` to
  localise). Fig `HiddenValley_v1v2_vs_mrms.png`. **R(A) tried first + abandoned:**
  raw NEXRAD PhiDP (offset+folding) makes `pyart.correct.calculate_attenuation_zphi`
  blow up (A ~10× high) — needs PhiDP preprocessing; `kdp_maesaka` self-regularises
  (and that pyart fn errors unless you pass `temp_ref="fixed_fzl", fzl=`).
- **v2 step 2 done — DEM beam-blockage (wradlib):** `intensity_stack(blockage_dem=,
  cbb_max=0.5)` + standalone `beam_blockage(radar, aoi, dem)`. Terrain sampled at
  gate lon/lat → `wradlib.qual.beam_block_frac`(th, beam-height from pyart
  `get_gate_lat_lon_alt`, half-power radius r·sin(bw/2)) → `cum_beam_block_frac`;
  masks cells CBB>`cbb_max` and emits a `cbb` field (MRMS-RQI analogue). HV ≈ 0
  (KRGX mountaintop 0° tilt → unobstructed). Needs `wradlib` (conda-forge).
  Full multi-tilt hybrid-scan (pick lowest *unblocked* tilt) + VPR are the
  remaining extensions — only matter where the low tilt is blocked (not HV).
  **Use reflectivity-only Z–R (`za`) for strict cross-era consistency** (dual-pol
  only post-~2012).
  **Don't use Level III QPE for multi-era work** — its archive mixes algorithms
  (legacy single-pol DPA, 4 km HRAP) vs (dual-pol DAA/DPR, 0.25 km); process
  Level II with one fixed algorithm instead.
- **Rate-method benchmark (HV, 2026-06-25 — `HiddenValley_rate_method_benchmark.md`):**
  compared single-radar i15 methods vs gauges (truth) + MRMS. **Maesaka R(Kdp)
  (current v2) wins — nothing beat it.** CSU `calc_blended_rain` decision tree was
  *worst* on skill (gauge corr 0.29 — raw uncalibrated NEXRAD Zdr + Kdp noise make
  its R(Z,Zdr)/R(Kdp,Zdr) spiky); Vulpiani Kdp smoother (corr 0.85) but over-reads
  (3.1× gauge); Beard-Chuang coeff (50.7, 0.85) matches MRMS bias (1.00×) at the
  same Maesaka noise (cheap optional swap). **Operational R(A) built + tested and
  it FAILS for i15:** the raw-PhiDP blow-up is fixed by feeding the clean
  reconstructed PhiDP `kdp_maesaka(...)[2]` (NOT `[1]`, which is Kdp stdev) to
  `calculate_attenuation_zphi`, and ρHV<0.9 censors hail — but R(A) still
  over-reads i15 **20× gauge / peak 1313** because Z-PHI attenuation spikes
  dominate the running-MAX (R(Kdp) is clipped 0–7; R(A) is unbounded). R(A) is for
  gauge-corrected *accumulation* in pure rain, not peak-i15 in hail — **don't
  re-attempt R(A) for i15.** **Uniform per-scan `rate_cap`** (added to
  `intensity_stack` / CLI `--rate-cap`; operational hail-cap analogue) at
  120 mm/h *helps R(Kdp)* — gauge corr 0.68→0.74 at no bias cost (clips Kdp
  spikes), so **recommend `--rate-cap 120`** — and bounds R(A) (peak 1313→120)
  but only halves its over-read (20→9× gauge), confirming R(A) over-reads
  *broadly* in hail, not just via spikes. Libs in GISMan: `wradlib` (blockage,
  package dep), `csu_radartools` (benchmark-only, not imported by the package).
- **Diagnostic panels + virtual gauges (ported from D. Cavagna's `MRMS_stack`,
  2026-06-26):** `plot.diagnostic_panels(dir, key, which=("tpki15","total","rqi",
  "shsr"))` tiles saved field tifs (CB-safe cmaps — NOT the original's jet/turbo;
  time-of-peak unwrapped past midnight); CLI `panels`. `mrms.virtual_gauge_timeseries(
  points, start, end)` samples PrecipRate (2-min → 1-min interp) + hourly
  MultiSensor QPE (Pass2→Pass1→RadarOnly fallback) at point(s) via
  `rasterio.transform.rowcol` on the AOI-window transform (MRMS lons are negative
  here — no +360 like the original) → per-gauge DataFrame (`total_mm`,
  `i5/i15/i30/i60_mmph` — 15 via the `(i16+i14)/2` estimator, others trailing
  windows — plus `i60_qpe_mmph`/`total_qpe_mm` on the hour rows);
  `plot.plot_virtual_gauge(df, name)` = 2-panel (intensities + cumulative
  PrecipRate vs QPE), matches `Figures/TruckeeRiver_VG2`. CLI `vgauge --point
  LON,LAT[,NAME]` (repeatable) or `--points-file`. **User chose 5/15/30/60-min
  intervals** (original used I2). **Validated:** VG i15 reproduces the MRMS
  i15max *field* (HV 68.0 vs 68.5; Six Mile 50.4 vs 50.8) — confirms the
  point-sampling + estimator. Deliverables: `HiddenValley_panels.png`,
  `HiddenValley_VG_*.png`.
- **Panels + VG enhancements (2026-06-26):** `diagnostic_panels` now drapes the
  **full main-map context** (hillshade, NHD/TIGER/GNIS reference, perimeter,
  north arrow, lat/long ticks, clip) per panel — same look as `drape_i15`. Its
  figure width is **auto-sized from the map aspect** (`panel_h·aspect + 1.5` per
  column, `panel_h=5.6`) so equal-aspect panels don't float in over-wide cells —
  that float was the dead white space between columns for tall/narrow AOIs (HV
  ~0.4 w:h). Pass `figsize=` to override.
  `vgauge --gauges` drops a VG at *every* real Synoptic station + overlays the
  real series (`gauges.gauge_timeseries` → real-gauge time-series counterpart,
  same df shape; needs `$SYNOPTIC_TOKEN`); `--atlas` → `plot.virtual_gauge_atlas`
  (one small panel/gauge, radar VG solid + matching real gauge dashed). Per-gauge
  CSVs (one file each, Cavagna-style) go to `<out-dir>/RainGaugeData/`
  (`*_vgauge_*.csv` + `*_gauge_*.csv`). **NEXRAD fallback:**
  `nexrad.virtual_gauge_timeseries` (CLI `vgauge --source nexrad --method
  za|kdp`) samples the L2 lowest-tilt rate at points (no QPE col) for pre-2020 /
  MRMS-gap events. **Decision — kept VG extraction SEPARATE from the raster
  stack** (`i15_storm_day`): the stack processes only *wet hours* over the full
  AOI whereas VGs want the *full window* at points, the NEXRAD fallback needs its
  own path anyway, and the re-fetch is a bounded one-pass cost over a small
  bbox — decoupling > the modest duplicate-fetch saving. Real-gauge series math +
  atlas overlay validated offline (no token); bulk MRMS atlas + NEXRAD fallback
  validated on HV. Deliverables: `HiddenValley_vg_atlas.png`, `RainGaugeData/`.
- **Multi-source 3-way atlas (2026-06-26):** `virtual_gauge_atlas` now takes
  `sources = {label: {name: df}}` (any number of radar sources, each a solid
  coloured line — MRMS blue, NEXRAD red) + `real_series=` (matching real gauge
  overlaid **dashed black**); a bare `{name: df}` is still accepted (back-compat).
  CLI **`vgauge --nexrad`** is an *additive* flag (distinct from `--source`):
  in the post-2020 both-available epoch it runs the single-radar NEXRAD R(Kdp)
  VGs *alongside* MRMS so the atlas shows MRMS+NEXRAD+gauge together. Per-source
  CSVs are now source-labelled — `<key>_vgauge_mrms_<name>.csv`,
  `<key>_vgauge_nexrad_<name>.csv`, plus real `<key>_gauge_<name>.csv` — and the
  atlas writes to **`<key>_vg_atlas.png`** (was `_atlas.png`; renamed for
  consistency with the function + prior deliverable). Token: pass via
  `$SYNOPTIC_TOKEN` or `--token`, used transiently — **never stored in repo/memory**.
  Cosmetic: `tight_layout(rect=[0,0,1,0.97])` keeps the suptitle off the top-row
  panel titles. **Validated on HV 20260619 (KRGX):** 51 stations, all 3 series
  per panel, 51/51 MRMS↔real name-matched; NEXRAD R(Kdp) generally rides just
  above MRMS, both above the gauge (consistent with the bias study). The atlas
  re-renders instantly from the `RainGaugeData/` CSVs (no re-fetch) — see
  `scratchpad/regen_atlas.py`.
- **Per-gauge detail figures (2026-06-26):** `plot.virtual_gauge_detail(sources,
  name, real_series=)` is the *big-figure* counterpart of one atlas panel — a
  4×1 vertical stack (**cumulative rainfall, then I60 / I30 / I15**) for one
  gauge, same `{label:{name:df}}`+`real_series` data model and line styles as the
  atlas (MRMS blue, NEXRAD red, gauge dashed black). Both it and the atlas share
  the `_overlay_vg_series` helper (atlas lw 1.0/1.1, detail 1.8/2.0). CLI
  **`vgauge --detail`** writes one per gauge → `<out-dir>/VirtualGaugeFigures/`
  `<key>_vgdetail_<name>.png`. Kept the existing `plot_virtual_gauge` (single
  MRMS source, 2-panel + QPE, for explicit `--point`s) *unchanged* — `_detail` is
  the multi-source comparison view, `plot_virtual_gauge` the QPE deep-dive.
  Validated on HV (51 figs from the CSVs, no re-fetch — `scratchpad/gen_detail.py`):
  e.g. Six Mile Canyon cumulative MRMS 19.5 / gauge 16.5 / NEXRAD 13.3 mm; Hidden
  Valley 48 / 20 / 28 mm (MRMS high, gauge low, NEXRAD between).
- **Gauge tempo on detail figures (2026-06-29):** each `virtual_gauge_detail`
  title now annotates the **ground gauge's native reporting cadence** —
  `(ground gauge: ~5-min reporting)` / `~hourly (60-min)` / `~daily (1440-min)`.
  Source is the per-gauge **`report_min`** field (median minutes between
  *precip-bearing* obs — **measured, not metadata**, so a multi-variable ASOS like
  KRNO reads ~hourly, not its 5-min temp cadence; Synoptic exposes no reliable
  per-gauge *precip* interval), already written into `<key>_gauges.geojson`. It is
  carried to the figure via the real series' **`df.attrs['report_min']`**: set by
  `gauges.fetch_gauge_event` + `gauge_timeseries`, persisted as a `report_min`
  column in the per-gauge `RainGaugeData/*_gauge_*.csv`, and restored by
  `load_event_series` (**preferring the geojson field**, so OLD stores whose CSVs
  lack the column still get it — HV: 51/52 restored). One shared cadence helper
  `gauges._report_min` (now also used by `gauge_intensities`); `plot._tempo_label`
  formats it (`None`→omitted, finite handles 5/10/15/60/1440). `virtual_gauge_detail`
  takes an explicit `report_min=` (CLI passes it from the series attrs) and falls
  back to the series attrs for direct library callers. **Why it matters:** coarse
  reporters get bursts smeared by the 1-min interpolation → i15/i30 read **low**
  (the daily RENO WFO draws one step vs MRMS/NEXRAD sub-hourly structure; total is
  cadence-insensitive), so the tempo flags how far to trust the peaks. **Detail
  figs only** — atlas panel titles + the `plot_virtual_gauge` 2-panel are
  unchanged. Validated offline on the HV Box store (no token / no re-fetch): real
  MRMS+NEXRAD+gauge detail figs across 5/10/15/60/1440-min in
  `scratchpad/TEMPO_montage.png`. **Cadence filter for the pipeline (2026-06-29):**
  `gauges`/`vgauge` gained **`--max-report-min MIN`** (same flag name as
  `run`/`compare`, here a *station-set* filter in `_filter_stations`): keep only
  gauges whose `report_min <= MIN` — **`--max-report-min 60` = hourly or finer**,
  dropping daily/coarse reporters from the atlas + detail set before the clean
  refresh. Reads `report_min` straight from the series `df.attrs` (so it composes
  with `--wet-only`/`--max-dist-km`); a NaN cadence (unmeasurable: no/sparse precip)
  is **left to `--wet-only`**, not dropped on cadence. HV: `60` drops just the daily
  RENO WFO (52→51), `15` keeps 5/10/15-min only (52→48), `--all-gauges`-off default
  +`60` → 20 wet hourly-or-finer.
- **`zoom` subcommand (2026-06-26) — re-view, don't re-run.** To zoom into a
  sub-region of a processed event, `zoom --from-dir <dir> --from-key <key>` +
  a smaller `--bbox/--aoi` re-renders the main map + diagnostic panels **clipped
  to the sub-AOI straight from the existing rasters — no MRMS re-download**
  (`drape_i15`/`diagnostic_panels` `clip=` is a *view* clip: sets ax x/y-lims +
  re-fetches reference at the tighter extent). Rationale: **MRMS is native ~1 km
  with no finer tier**, so a fresh run just re-downloads identical radar; the
  *only* product that benefits from zooming is terrain → `--refine-dem
  --resolution N` re-pulls just a finer 3DEP DEM/hillshade for the zoom extent.
  `--crop-rasters` writes cropped (CRS-aware via `_crop_rasters`, reprojecting the
  4326 zoom box per raster) copies → a self-contained zoom folder. `--gauges`
  overlays the source `<from-key>_gauges.geojson`. **Gotcha:** `rioxarray.open_
  rasterio` left open until GC throws a benign `Error in sys.excepthook` at
  interpreter exit — `_crop_rasters` uses `with open_rasterio(...) as da: ....load()`
  to close promptly (the lazy `_load` in plot.py is the same pattern but hasn't
  surfaced it). Validated on HV (sub-bbox of the big AOI): zoomed map + panels +
  11 cropped rasters, panels auto-size to the zoom aspect. `zoom` passes one
  `zoom_alpha` (= `--alpha`, else the project default) to **both** `drape_i15`
  (map) and `diagnostic_panels` (panels).
  **Climate maps in zoom (2026-06-28):** `zoom` now also produces the NOAA Atlas 14
  **climatology comparison + anomaly** maps for the sub-AOI, **on by default** (like
  map/panels; `--no-climate` to skip). Implemented by extracting `_run_climate(args,
  *, src, from_key, out_dir, key, bounds, clip_gdf, hs_path, gauges_arg, title=)`
  from `_cmd_climate` — both commands now call it; the climate CLI knobs live in a
  shared `_add_climate_opts(p)` parser helper used by **both** the `climate` and
  `zoom` parsers (`--ari/--durations/--region/--stat/--anomaly-cmap/--obs-smooth/
  --obs-smooth-radius/--eps/--shared-row-scale/--no-comparison/--no-anomaly`). The
  zoom clim is **fetched fresh for the tighter extent** (its own `atlas14_cache`)
  while the observed field is reused from the source event (MRMS has no finer res),
  Gaussian-1km-smoothed by default. The whole climate step is **guarded** in
  `_cmd_zoom` (try/except): a failed Atlas 14 fetch warns and still leaves the
  map+panels (same resilience pattern as `--refine-dem`). Outputs:
  `<key>_climate_compare.png`, `<key>_clim_i{d}.tif`, `<key>_anom_i{d}.png`/`.tif`.
  Validated on HV (sub-bbox): clim tifs are smaller than the full-AOI ones (clipped),
  the compare figure axes span only the zoom box; `--no-climate` writes 0 files / no
  fetch; the `climate` command itself regression-passes after the refactor. Since
  `pick`'s suggested `zoom` command is unchanged, the picker-driven workflow now
  yields climate maps automatically.
  **Smooth radar for ALL zoom rainfall maps (2026-06-28):** every rainfall *display*
  in the zoom folder is Gaussian-1km-smoothed by default, not just the climate ones —
  the i15 map (`drape_i15`) and the rainfall panels (`diagnostic_panels`), driven by
  the **same `--obs-smooth`/`--obs-smooth-radius` knob** (`osm`/`osr` resolved in
  `_cmd_zoom` and passed as `field_smooth=`/`field_smooth_radius_km=`). Added
  `field_smooth`/`field_smooth_radius_km`/`smooth_power` params to **`drape_i15`** (a
  `_load_i15()` helper smooths the field via `smooth_dataarray` before `_to_crs`/
  reproject) and **`diagnostic_panels`** (smooths only fields in **`_MASK_DRY`** — the
  intensities + depths incl. `total_mm`; the categorical `tpki15`/`rqi`/`shsr`/`cbb`
  stay raw — and appends "(Gaussian 1 km)" to those panel titles). Both default OFF
  at the library level (so standalone `map`/`run`/`panels` are unchanged); only
  `zoom` wires the on-by-default. The map title gains a `(zoom, gaussian 1 km)` note.
  **Display only — `--crop-rasters` still writes RAW cropped tifs** (data integrity;
  use `smooth --write` for smoothed field tifs). Validated on HV sub-bbox: i15 map +
  `total` panel visibly smoothed (peak ~80→65), `tpki15`/`rqi`/`shsr` panels raw.
  **Op tip — re-render a 1 m-refined zoom WITHOUT re-fetching the DEM:** point
  `zoom --from-dir`/`--from-key` at the **zoom folder itself** (e.g. `--from-dir
  <…>/HiddenValley_zoom --from-key HiddenValley_zoom`, same `--out-dir`/`--key`, no
  `--refine-dem`/`--crop-rasters`). It reuses that folder's existing 1 m hillshade +
  cropped observed tifs and just re-renders the figures (map/panels/climate) with
  current defaults — the slow/flaky 1 m 3DEP fetch is skipped, DEM/hillshade/cropped
  rasters untouched. Used 2026-06-28 to update the Box `HiddenValley_zoom` (bbox
  `-119.7534 39.4597 -119.6476 39.5585`, recovered from the original command) to a
  smoothed map+panels + new climate set on the kept 1 m terrain (~30 min, all 1 m
  renders). The exact zoom bbox is NOT recoverable from the cropped outputs alone
  (crop pad + 5070→4326 bulge) — get it from the original command / `pick`.
  **Hillshade render optimization (2026-06-28) — ~30 min → ~1 min (≈29×).** The
  zoom was slow because the 1 m hillshade is ~167 M cells (236 M after the UTM
  warp) — ~25× more than a 200–300 dpi figure can resolve — and it was re-loaded +
  re-reprojected + imshow'd **per figure** (6 figs / ~14 panels), with the rainfall
  field upsampled to that same 236 M-cell grid, and each ~1.9 GB array thrashing
  memory. Profiled: 1-panel imshow @1 m = 7.8 s vs 0.3 s decimated (~26×). Fix:
  **`plot._prepare_hillshade(hillshade, work_crs, max_px)`** — reproject to UTM +
  downsample to `max_px` on the long side in **one GDAL warp** (`Resampling.average`;
  never materializes full-res; extent preserved exactly), returning
  `(hillshade_da, resolved_crs)`. `_cmd_zoom` + `_cmd_climate` prepare it **once**
  and pass the same small array + explicit `work_crs` to every figure (`drape_i15`,
  `diagnostic_panels`, `climatology_comparison`, `anomaly_map`) — each one's
  `_to_crs(hs, wc)` is then a no-op, so they reuse the one prepared grid and the
  field upsamples to it (5.9 M not 236 M cells). Cap scales with `--dpi`
  (`cli._render_px = max(2500, dpi·12)`) so it never under-samples even at high
  dpi; **on-disk 1 m DEM/hillshade untouched** (render-only). HV zoom full run
  (map+panels+climate, incl. Atlas 14 fetch + 6 reference fetches) **62 s**, output
  pixel-indistinguishable from the 1 m render. Standalone `map`/`run` unchanged
  (still pass a path + `work_crs="UTM"`; they use small 10 m hillshades anyway).
- **Drape alpha — one project-wide default (2026-06-28):** the rainfall/field
  layer draped over the hillshade defaults to **`alpha=0.32`** *everywhere* (so
  terrain/basemap reads through), set in **one place** —
  `plot.DEFAULT_FIELD_ALPHA`. Every map function (`drape_i15`, `_draw_field`,
  `diagnostic_panels`, `climatology_comparison`, `anomaly_map`, `bbox_picker`)
  takes `alpha=None` and the two leaf drawers (`drape_i15`, `_draw_field`)
  resolve `None → DEFAULT_FIELD_ALPHA`; the wrappers just forward. The CLI passes
  `alpha=args.alpha` (None unless `--alpha` given) on **every** map command —
  `map`/`run`/`nexrad`/`zoom`/`climate`/`panels`/`compare` — so **`--alpha`
  overrides uniformly** (added `--alpha` to `panels` + `compare`, which lacked
  it; the `run --compare` residual map previously ignored `--alpha` — fixed).
  This replaced the old scattered per-figure defaults (0.45 main maps, 0.58/0.5
  drape, 0.6 nexrad/climate-compare, 0.72 anomaly). To change the global default,
  edit `DEFAULT_FIELD_ALPHA`. Note: hillshade opacity is separate
  (`--hillshade-alpha`); the anomaly map's integer contour lines stay crisp
  regardless of field alpha, so 0.32 keeps the recurrence-multiple breaks legible.
- **`--refine-dem` resilience (2026-06-26):** the **3DEP 1 m dynamic WMS is
  slow/flaky** — it times out at py3dep's hard-coded **120 s** request limit
  (not env-tunable, not exposed), and on a bad stretch even a small ~3 km box
  fails while 10 m fetches in ~6 s. So `dem.get_dem` / `fetch_dem_and_hillshade`
  now **retry** (default `retries=2`, linear `retry_wait` backoff) on the WMS
  timeout, and `zoom --refine-dem` **catches a persistent failure and falls back
  to the source hillshade** (renders the zoom at coarse terrain) instead of
  aborting with no figures; `_crop_rasters` then crops the source DEM/hillshade
  too (the `skip=` is gated on whether the refine actually *succeeded*, tracked by
  `refined`). **HV coverage:** 1 m = 100% (lidar) but `3 m = 0%` (only 1 m + 10 m
  + 30 m source footprints), so **1 m is the sole finer-than-10 m option there** —
  and `py3dep` has no fast 1 m path (`static_3dep_dem` is the 10 m seamless; the
  dynamic WMS is the only 1 m route). For radar-scale work 10 m already far
  out-resolves the ~1 km i15, so 1 m terrain is mostly cosmetic; retry 1 m later
  when 3DEP is responsive. **Key sizing gotcha:** the zoom's default
  `--pad-deg 0.05` (~5.5 km) *quintuples* the 1 m area (a ~9×11 km box → ~385 Mpx)
  and blows the 120 s timeout even on a healthy server — so `zoom --refine-dem`
  caps the DEM pad at `min(pad, 0.003)` (the view is clipped to the box anyway, so
  the hillshade only needs a sliver of edge). With the cap + a responsive server,
  HV's full zoom box fetched 1 m in ~14 min (the hillshade over ~100 Mpx is the
  slow part; a quick `get_dem(box, resolution=1, retries=0)` ~70 Mpx probe in ~60 s
  is a good health check before committing). 1 m DEM/hillshade tifs are big — keep
  them out of git.
- **NEVER let rioxarray pick the resampling on a DEM (2026-07-31) — `dem.get_dem`
  no longer calls `py3dep.get_dem` at 10/30/60 m.** Found from finished figures:
  every coarse hillshade carried *two* nearest-neighbour artefacts, while 1 m
  hillshades were clean.
  **(A) a ~45° "corduroy" hatch on smooth slopes** — `py3dep.get_dem` internally
  does `static_3dep_dem(...).rio.reproject(5070)` with rioxarray's **default
  (nearest)** resampling and **no target resolution**. The seamless 3DEP VRT is
  **EPSG:4269 at ⅓ arc-second**, so that warp is a rotation plus a non-integer
  scale and nearest aliases along diagonals. It also lands on an arbitrary
  **~9.38 m** grid, not the 10 m requested.
  **(B) an axis-aligned ~15-px grid** — warping that 9.38 m to 10.0 m, nearest
  again: 10/9.3817 = 1.0659, so one row and column is dropped every ~15 cells.
  **Why 1 m was clean:** resolutions **outside {10, 30, 60}** route to the
  *dynamic* 3DEP image service, which resamples server-side and returns 5070
  directly — no client-side nearest warp at all.
  **Fix:** read the VRT via `py3dep.static_3dep_dem` on its **native 4269 grid**
  (+20 cells so the kernel never reaches past the data), warp **exactly once**
  with `resampling=` (new arg, default `DEFAULT_RESAMPLING = bilinear`), then trim.
  `_needs_warp` skips the warp entirely when the DEM already arrives on the target
  CRS *and* within 1% of the target spacing — the dynamic path usually does.
  **Bilinear over cubic:** cubic overshoots at cliffs and the hillshade turns that
  into bright rims, for a 2% difference in elevation RMS.
  **Scored** against 3DEP 2.5 m averaged onto the same 10 m grid (an independent
  measurement of the same ground, not another resampling of the same numbers):
  elevation RMS **2.81 m → 0.91 m**, hillshade roughness **0.064 → 0.047**
  (reference 0.041), short grid period gone. The control that pins the blame on
  *nearest* rather than on the number of warps: VRT → **nearest** 10 m (one warp)
  still hatches at roughness 0.060.
  **CLI (2026-08-14):** `--resampling` on `dem` / `run` / `zoom --refine-dem`
  (`nearest|bilinear|cubic|cubic_spline|lanczos|average` — the ones that mean
  something on a continuous surface; `mode`/`sum`/`q1` are valid rasterio
  resamplers but nonsense for a DEM, so they are not offered). It defaults to
  **`None`, not a name**, so `dem.DEFAULT_RESAMPLING` stays the single source
  of truth — `_resolve_resampling(None)` returns it. `nearest` is deliberately
  still reachable: reproducing the hatch is the control experiment above.
  **Downstream warning:** hillshading differentiates, so the artefact is far
  louder in shaded relief than in elevations — but it is *in* the elevations.
  Anything derived from a 10/30/60 m `get_dem` before this date is affected, worst
  at high derivative order. Measured on 12 basins in the debris-flow study:
  relief unchanged, slope biased **high by ~1.8° (6%)** but r = 0.997 (a
  near-constant offset), aspect fractions move ~0.012, and **plan curvature
  halved with r = 0.79**. Re-derive second-derivative metrics.
  Tests: `tests/test_dem.py` (16, offline, `py3dep` replaced by a recording double).
- **`pick` subcommand (2026-06-26) — cross-platform browser bbox picker.**
  `plot.bbox_picker(i15, hillshade=, reference=, perimeters=, gauges=, cmd_prefix=,
  cmd_suffix=)` renders the event's i15 map **through `drape_i15`** so it carries
  the **full main-map context** (labelled NHD/TIGER/GNIS reference, perimeter,
  gauges, north arrow, lat/long ticks, colorbar) in the same **UTM** projection as
  the production maps, embeds it as a base64 PNG in a **self-contained HTML** file
  with a vanilla-JS canvas: drag a rectangle → `--bbox W S E N` + a ready-to-run
  `zoom` command (copy buttons). CLI `pick --from-dir --from-key` (reference ON by
  default; `--no-reference` skips the fetch; `--gauges`, `--perimeters`,
  `--local-roads`) writes `<from-key>_pick.html` and `webbrowser.open()`s it.
  **Coords stay accurate without a 4326 render:** capture the map *axes* rectangle
  (`ax.get_position()` after `fig.canvas.draw()`, as image fractions) + its four
  UTM corners transformed to lon/lat; the JS maps a click via **bilinear over the
  4 corners** (exact for an AOI this size). **Chosen for cross-platform robustness**
  (Windows/Linux users expected): a browser is the one GUI every desktop has — no
  matplotlib GUI backend / display / X11 dependency (rendering is headless
  Agg/PIL). **Gotcha:** without a clip, `drape_i15` frames to the hillshade
  footprint and the Albers(5070)→UTM rotation shows as gray wedges — so
  `bbox_picker` **clips tightly to the i15 extent** (like the production maps) →
  clean north-up background. `max_px=1700`/`dpi=160` keep the HTML lean (~0.7 MB
  no-ref). Validated on HV: tight corners (~-119.84/39.26/-119.37/39.95, small UTM
  tilt per-corner), tokens filled, background matches the main map, bilinear math
  checked. Matplotlib `RectangleSelector` (macosx/tk) was the alternative but is
  display-dependent → rejected for portability.
- **NOAA Atlas 14 rainfall climatology + anomaly (`atlas14.py`, CLI `climate`,
  2026-06-28):** puts a storm in climatological context — observed I15/I30/I60 vs
  the **1-yr** (default ARI) precipitation-frequency climatology, plus per-duration
  **anomaly = observed/climatology** maps. **`pfdf.data.noaa.atlas14` is
  point-only** (`download(lat,lon)` → one CSV; no raster path) so it was
  *bypassed*; instead fetch NOAA's authoritative **gridded ASCII** product:
  `hdsc.nws.noaa.gov/pub/hdsc/data/<region>/<region><ari>yr<dur><stat>.zip`
  (`stat` `a`=mean, `al`/`au`=90% CI; `_ams` suffix=annual-max vs default PDS).
  Each zip ≈1.5 MB → `.asc`+`.prj`+`.xml`, **EPSG:4269**, ~800 m, NODATA −999,
  **units 1000ths-inch depth** → convert `i_mmph = mils/1000·25.4·(60/dur_min)`.
  **No new deps** (rasterio AAIGrid + rioxarray). Region auto-picked from the AOI
  by bbox-containment then nearest-centroid (`data/atlas14_regions.csv`, 24
  regions built from the grid headers; `sw` = Vol 1 Southwest covers NV/CA/AZ/UT/NM;
  `--region` overrides). `climate` reuses an event's `<from-key>_i{d}max.tif` (MRMS
  *or* NEXRAD — source-agnostic), writes `<key>_clim_i{d}.tif`, a 3×N
  `climate_compare.png` (rows=durations, col0 clim col1 observed, **independent
  per-panel sequential scale** by default — clim is smooth/small vs the peaky
  observed, a shared scale washes it out; `--shared-row-scale` to share), and
  `<key>_anom_i{d}.tif`+`.png` (**diverging `cmc.vik` `TwoSlopeNorm(vcenter=1)`,
  integer contours** 1×/2×/3×… per the scientific-plots default). Always clips to
  the AOI (else Albers→UTM gray wedges). **Rotated-gray-box / white-wedge fix
  (2026-06-29) — `plot._fill_hillshade_nan`, used by ALL maps:** the hillshade
  (5070) reprojects to a **rotated** rectangle in UTM with NaN corners. Feeding the
  2-D masked array to `imshow(cmap="gray")` renders those corners either as a flat
  gray box (matplotlib's **data-stage resampling fills the mask** on the
  colorbar/`tight_layout` axes resize) OR as stark **white wedges** — and *which*
  you get is **layout/dpi-dependent + fragile** (the production `drape_i15` maps
  happened to get the gray-fill, which blends into terrain → looked clean; the
  climate maps got white wedges + a visible rotated rectangle → looked wrong). Tried
  + rejected: clipping the hillshade to the AOI (floats the map in white space);
  explicit RGBA alpha=0 (gives stark white wedges); `interpolation_stage="rgba"` /
  `"nearest"` (no effect — auto-downsample resampling ignores them). **Real fix:**
  `_fill_hillshade_nan(hsv)` replaces the NaN corners with the **mean terrain tone**
  (`np.nanmean`) before `imshow`, so the corners **blend into the hillshade →
  terrain fills the view frame** reliably (not layout-dependent). Applied in BOTH
  `_draw_field` (anomaly + `climatology_comparison` + `diagnostic_panels`) AND
  `drape_i15` (i15 / nexrad / compare maps) so EVERY map renders the hillshade
  identically — the user's requirement that all maps match `HiddenValley.png`.
  The hillshade is left UNclipped (terrain fills the frame). Validated: anomaly +
  climate-compare + panels + i15-with-gauges all fill the frame, no gray box, no
  white wedges, matching `HiddenValley.png`. **Extent auto-match (2026-06-29):** the
  observed (`i15max`) footprint = AOI + the ~0.05° MRMS fetch pad, so clipping
  `climate` to it framed *wider* than the `map`/`run` i15 maps (which clip to the
  `--aoi`/`--perimeters` polygon). Fix — **`climate` now auto-matches the event
  AOI**: with no `--aoi`/`--bbox` it calls `cli._find_event_aoi(from_dir, from_key)`
  (looks for `<from_key>_aoi.geojson`, then user-placed `<from_key>_AOI.{kmz,geojson,
  gpkg,shp}`), clips to it, and draws it as the perimeter — IDENTICAL extent to the
  i15 maps (measured: same UTM xlim/ylim, same `clip_margin` 0.04); falls back to the
  i15max footprint (with a note) only if none found. To make the AOI reliably
  present, **`dem`/`i15`/`run`/`nexrad` now save `<key>_aoi.geojson`** via
  `cli._save_event_aoi` (the resolved `--aoi`/`--bbox` as a GeoJSON polygon). HV's
  `HiddenValley_AOI.kmz` is found by the auto-match; its climate deliverables match
  `HiddenValley.png` exactly, AOI perimeter drawn.
  `plot._draw_field` is the shared
  per-panel context helper extracted from `diagnostic_panels` (now used by panels
  + both climate figures; regression-checked). **Anomaly() reproject_matches clim
  (4269) onto the observed (4326) grid then divides**, masks clim≤eps. **HV
  20260619 validated:** gridded 1-yr depth (15/30/60 = 4.65/6.27/7.75 mm) matches
  NOAA's **point** PFDS estimate (5/6/8 mm, integer-rounded — the service pfdf
  wraps), confirming units; anomaly median ≈1× (i15 0.99, i30 0.89, i60 0.79 —
  storm near the 1-yr overall) with the convective core ~5×. Keep grids/cache +
  clim/anom tifs **out of git** (cache → `<out-dir>/atlas14_cache/`).
  Deliverables: `HiddenValley_climate_compare.png`, `_anom_i{15,30,60}.png`.
  **Observed-field smoothing now ON by default (2026-06-28):** `climate` smooths
  the observed radar field **Gaussian, 1 km** before *both* the comparison figure
  (col1, labelled `Observed (Gaussian, 1 km)`) and the anomaly (the divide + the
  written `<key>_anom_i{d}.tif`), so the peaky ~1 km radar reads against the
  smooth ~800 m climatology and single-pixel spikes don't dominate the
  recurrence-multiple anomaly. CLI `--obs-smooth {gaussian,uniform,median,idw,none}`
  (default gaussian; `none`=raw) + `--obs-smooth-radius KM` (default 1.0). Wired
  via **`smoothing.smooth_dataarray(src, method, radius_km)`** (new) — accepts a
  path *or* a DataArray, NaN-aware on the native grid, **returns the input's own
  dim order** (don't squeeze to 2-D: `atlas14.anomaly`→`ratio.rio.to_raster`
  needs band,y,x; a squeezed 2-D obs raised `InvalidDimensionOrder` only at write,
  not in-memory). `climatology_comparison` gained `obs_smooth/obs_smooth_radius_km`
  (library default OFF; the CLI sets the 1 km default). **HV check:** i15 native
  peak 102.2→79.1 (mean 24.28→24.30 preserved, no NaN bleed); anomaly i15 core
  5.3×→4.0× (95th 3.27→2.90, median ~1×). NB this smooths *display/anomaly* only —
  the stored `<key>_i{d}max.tif` observed fields are untouched; for smoothed
  *field* tifs use `smooth --write`.
- **Radar-field smoothing (`smoothing.py`, CLI `smooth`, 2026-06-28):** smooths
  a peaky radar intensity field (for display next to the smooth Atlas 14
  climatology, and to test radar-gauge representativeness). **Four NaN-aware
  isotropic low-pass methods** in `METHODS`: `gaussian` (recommended default),
  `uniform` (boxcar mean = the literature's N×N window), `median`
  (edge-preserving rank), `idw` (`1/d^power`, center weight 1.0 — the
  moving-window analogue of point IDW). **NaN handling:** linear filters
  (gaussian/uniform/idw) use **normalised convolution** `smooth(data)/smooth(valid)`
  (zero-fill, `mode=constant`) — no zero-bleed across masked cells, edge +
  Gaussian-truncation renormalised for free; **median is a rank filter** so it
  uses `scipy.ndimage.generic_filter` + guarded `nanmedian` over an odd disk
  (normalised conv does NOT apply). One knob **`radius_km` = nominal scale ≈
  Gaussian σ**, mapped per method (gaussian σ=r/cell; box/median half-width
  `round(√3·r/cell)`, **odd** window `2h+1` — even windows shift the field; idw
  radius `2·r/cell`); `radius_km=0` = identity. `cell_size_km` from the transform
  at AOI-centre lat (0.01° MRMS ≈ 0.99 km). **Sampler gotcha:** the in-memory
  point sampler uses `rasterio.transform.rowcol` default **`op=floor`** to match
  `ds.sample`/`compare.sample_raster_at_points` (`round` is off by half a cell);
  parity-verified: `gauge_skill_sweep` `radius_km=0` == `compare.compare_storm`
  exactly. `gauge_skill_sweep` reuses `compare.comparison_stats` (same RQI +
  cadence screens). **Figures:** `plot.smoothing_comparison` (methods×radii grid,
  col0=raw, **single shared scale** so peak-flattening shows — unlike
  `climatology_comparison`'s per-panel) + `plot.smoothing_skill_plot` (corr/RMSE/
  ratio vs radius, optimum starred). **`--write METHOD --write-radius KM`** emits
  smoothed `<key>_<field>.tif` that flow into `compare`/`map`/`climate`. **SciPy
  only — no new deps.** **Interpretation (critical):** i15 is a *peak* metric, so
  smoothing mechanically lowers the radar's positive bias → judge agreement by
  **corr (↑) / RMSE (↓)**, NOT the bias *ratio* (falls toward 1× as a side
  effect; RMSE also partly confounded by bias). **HV 20260619 finding (MRMS):** corr
  peaks at raw/~0.5 km then *declines* for all methods (raw i15 r≈0.58), min-RMSE
  optimum is large (3–8 km, bias-driven), and the bias **ratio stays ≫1×
  (~3–3.5×) at every radius** — so smoothing does NOT meaningfully improve
  radar-gauge agreement at HV; the over-read is hail-driven QPE physics, not a
  pixel-scale representativeness artefact (consistent with the single-radar Z
  analysis). Deliverables: `HiddenValley_smoothing_compare.png`, `_smoothing_skill.png`+`.csv`.
  **HV finding (single-radar NEXRAD, 2026-06-28) — the benefit is product-dependent
  and scales with the field's intrinsic noise.** Ran the *same* skill sweep on the
  Level II intensity stack (`--from-key HiddenValley_nexkdp` / `_nex`, same 28-gauge
  `HiddenValley_gauges.geojson`): for **R(Kdp) (`nexkdp`, the current product)
  smoothing genuinely IMPROVES correlation** — raw i15 r 0.44 → **0.54 at ~1 km
  (+0.10, ~24% rel; all four methods)** then declines (textbook interior optimum),
  because Kdp is a noisy phase derivative ("Kdp texture") with real pixel-scale
  speckle that ~1 km (≈2 NEXRAD cells; grid is **0.475 km**, finer than MRMS 0.99)
  averaging removes. The reflectivity-only **`za` field gains only +0.01–0.03**
  (Z is smoother), and **MRMS gains nothing / loses** (−0.02 to −0.06; already
  multi-radar QC'd). So: smoothing helps where the retrieval is noisy (single-radar
  R(Kdp)), not where it's pre-smoothed (MRMS). **Bias is never fixed by smoothing**
  — ratio stays ~3.0× (R(Kdp)) / ~5.4× (za) / ~3.5× (MRMS) at every radius (hail
  QPE physics), reaffirming the MRMS conclusion for magnitude while flipping it for
  *pattern* agreement. min-RMSE optima stay large (3–8 km, bias-confounded — don't
  trust them). Deliverables: `HiddenValley_nexkdp_smoothing_compare.png`,
  `_nexkdp_smoothing_skill.png`+`.csv`. **CLI note:** NEXRAD field tifs are keyed
  `<event>_nexkdp_*` but the **hillshade + gauges are event-keyed** (`<event>_*`),
  so point `smooth` at a NEXRAD dir with `--hillshade <event>_hillshade.tif` +
  `--gauges-file <event>_gauges.geojson` (`--gauges-file` now also drives the
  figure's gauge overlay, not just the sweep).
- **Wet-gauge recurrence table (`compare.gauge_recurrence_table`, CLI `recurrence`,
  2026-06-28):** per wet gauge (peak i15>0) — observed peak I15/I30/I60, the **time
  of the I15 peak**, the **anomaly** (observed / 1-yr Atlas 14), and the
  **recurrence interval** of each peak. **Climatology source decision — NOAA PFDS
  *point* query, NOT pfdf and NOT raster sampling.** pfdf isn't installed (heavy);
  sampling our downloaded grids can't give the RI (we only have the 1-yr grids — RI
  needs the whole curve, 10 ARIs×3 dur = 30 grids). `atlas14.pf_point(lat,lon,
  stat="mean",series="pds")` hits NOAA's PFDS point service
  (`cgi-bin/hdsc/new/fe_text_{mean,upper,lower}.csv?...&data=intensity&units=metric&series=pds`)
  → DataFrame intensity(mm/h) indexed by ARI [1,2,5,10,25,50,100,200,500,1000],
  cols i{d}; **same Atlas 14 source as the maps** (validated: 15-min 1-yr 18 mm/h ≡
  gridded 4.65 mm depth), point-accurate, no region lookup (service locates by
  lat/lon), **no new deps** (urllib, like `fetch_grid`). **Recurrence interval =
  `atlas14.recurrence_interval(obs, aris, curve)` — log-log interp (NOAA publishes
  the quantile *curve*, not a closed-form inverse); **extrapolates the lowest
  log-log segment BELOW the 1-yr quantile → numeric sub-annual RIs** (PDS admits
  >1 event/yr; 18 mm/h→1.0 yr, 15→0.60, 10→0.19 cont. across 1 yr), inf above the
  top ARI (`>1000`), nan only for non-positive/missing input. (Earlier it returned
  nan below 1-yr → CLI showed a bare `<1`; the user flagged that sub-1-yr RIs are
  valid + wanted them filled, 2026-06-28.)** **Time-of-peak is offline:**
  read from the saved RainGaugeData per-gauge CSVs (`<key>_gauge_<safe>.csv`,
  `i15_mmph` column → `idxmax`) — **no Synoptic token needed**; CSV peak matches the
  geojson peak exactly for all 14 HV gauges. (NB the event-root `RainGaugeData/` is
  empty; the real one is `Raster_Data/RainGaugeData/`, 51 `_gauge_` + 102 `_vgauge_`
  CSVs from the 2026-06-26 `vgauge --gauges` run.) **"Wet" = i15>0 (14 at HV)**, not
  the 28 non-null (14 of those measured exactly 0). **HV finding:** only the 3
  convective-core gauges (Six Mile, Virginia Highlands, Hidden Valley) reached a
  real recurrence — i15 **~7–11 yr** (anom ~2×, peaks 21:45–22:15 UTC); every other
  wet gauge was **sub-1-year**. Gauge-grounded RI (~decadal at the core) is the
  ground truth vs the radar anomaly maps' ~5× (radar's known ~2× over-read). The RI
  column adds nonlinearity the anomaly hides (2× the 1-yr ≈ an 11-yr event).
  Deliverables: `HiddenValley_gauge_recurrence.csv` + `.md`. **Reusable for #4
  multi-storm.** **Gauge-tempo column (2026-06-29):** the table now includes
  **`report_min`** (CSV col 5; `tempo min` in the `.md`, after Gauge) so the reader
  can see which RIs to trust — i15-based RIs from a coarse reporter are smeared low.
  Pulled from the geojson `report_min` in `gauge_recurrence_table`; the md writer
  adds the column when present (`with_tempo`). **Refreshed against the 0.12-pad
  store (2026-06-29):** 22 wet, 7 ≥1-yr; **Dayton (5-min) I15 88.8 → RI 206.6 yr**,
  GW5980/EW1922 Stagecoach (5/15-min) 52/30 yr, core Six Mile/Virginia Highlands/
  Hidden Valley 7–11 yr; the tempo column flags the junk rows — **RENO WFO** (1440,
  RI ~2e-9) and **DESERT SPRINGS** (60-min, tell-tale flat I15=I30=I60=3.0). Run:
  `recurrence --from-dir <RD> --from-key HiddenValley --out-dir <RD> --key
  HiddenValley --aoi <RD>/HiddenValley_AOI.kmz` (PFDS point, no token).
- **Canonical gauge store + full-day/clip-to-rain-window model (2026-06-28).**
  Root problem found: two independent Synoptic fetches (`gauge_fields` for the
  `gauges`/`run`/`compare` geojson vs `gauge_timeseries` for the `vgauge` CSVs)
  drifted apart, and `vgauge` **discarded the stations GDF** (`real, _ =`) while
  `to_csv` dropped the `df.attrs` lon/lat — so the wider 51-station set's
  coordinates were lost (the 45-geojson missed the 7 wettest, incl. Dayton/
  Stagecoach). **Fix — one fetch, one store, reused everywhere:**
  **`gauges.fetch_gauge_event(aoi,start,end,out_dir,key,...)`** does a single
  `get_stations`+`get_rainfall` → writes `<key>_gauges.geojson` (coords + peak
  metrics + **`i15_peak_time`**) AND self-describing `RainGaugeData/<key>_gauge_
  <name>.csv` (series + `lon/lat/station_id` columns); the `gauges` CLI now uses
  it. **Design (user's, key):** the store is the **full storm-DAY record** (spatial
  analog: the whole AOI); each analysis **clips to the tight rain window** (temporal
  analog: a zoom) — so NO standardized window. **`gauges.storm_window(series,
  pad_min=30, cover=0.99)`** finds that window by **rain *mass* coverage** (central
  99% of Σ incremental depth — NOT any-nonzero, which stray single-gauge tips
  stretched to 27 h; mass-coverage gives HV 20:18–00:58, both cells).
  **`gauges.load_event_series(rgd,key,gauges_gdf)`** reloads the store's series +
  re-attaches coords (no re-fetch). `vgauge --from-dir/--from-key` now **reuses the
  store** (loads real series, auto-trims to `storm_window`, generates the virtual
  gauges over that window; `--refetch` to force live; `--pad-min`), and does NOT
  rewrite the full-day `_gauge_` CSVs (leaves the store intact). **Atlas filters
  (2026-06-28):** `vgauge --max-dist-km KM` (drop stations >KM from the AOI —
  `--bbox`/`--aoi`, else the `<from-key>_i15max` footprint; `_dist_to_aoi_km`/
  `_filter_stations`) + `--wet-only` (drop stations whose peak atlas-metric
  intensity < `--wet-min`, default 0.5 mm/h — a floor so traces like RENO WFO
  i15≈0.01 read as dry); filters the station set before the virtual gen + atlas
  (explicit `--point`s untouched). HV `--max-dist-km 3 --wet-only` → 68→21 gauges
  (vs 66+). **tz fix:**
  `_series_from_per_minute` now strips the grid tz so a tz-aware window can't zero
  the series. **HV re-fetch (token, transient from /tmp file, deleted after):
  45→116 stations, 14→22 wet** — recovered Dayton (i15 88.8, **RI ~207 yr**, 1 km S
  of AOI), Stagecoach (63.7/53.1, RI 30–52 yr, 0.3–2.2 km E, a later eastern cell
  peaking ~00:15), CW1177 Reno (22.1, 1.5 yr, in/near AOI). `recurrence` gained
  `aoi_bounds` → **`dist_to_aoi_km`/`in_aoi`** columns (+ `--max-dist-km`) and reads
  `i15_peak_time` from the geojson; **kept all 22, distance-flagged** (per user —
  near-AOI regionals are relevant). Window caveat: the `--date` default (~30 h) is
  fine for the store (clip handles it), but for a *clean* store prefer the actual
  storm day; the existing `_vgauge_` CSVs are still the old 6 h window until vgauge
  is re-run reusing the store. **Token is NEVER in memory** (the memory only
  *discusses* it; value-search empty) — pass `$SYNOPTIC_TOKEN`/`--token`, or a
  transient `/tmp` file read via `SYNOPTIC_TOKEN="$(cat …)"` then deleted.
- **Unified `gauges` pipeline (2026-06-29) — one command does it all.** Per the
  user, drawing gauges from Synoptic should by default *complete every step*:
  `gauges --aoi --date` now (1) writes the canonical store (as before), then (2)
  drops virtual gauges at every wet near-AOI station from **MRMS + NEXRAD (3-way by
  default)** and renders the comparison **atlas**, and (3) the per-gauge **detail
  figures** (`VirtualGaugeFigures/`). The store keeps the full storm DAY; steps 2–3
  clip to `storm_window` + filter to **wet** (default; `--all-gauges` to keep dry,
  `--wet-min` floor) near-AOI (`--max-dist-km`) stations. The atlas and the detail
  figures consume **one filtered station set** → their gauges are uniform by
  construction (the user's explicit requirement). Opt-outs: `--store-only` (old
  behaviour), `--no-atlas`, `--no-detail`, `--no-nexrad` (→ MRMS+real 2-way).
  Implemented by extracting **`cli._virtual_gauge_products(...)`** (the VG-source
  sampling + CSVs + atlas + detail body) from `_cmd_vgauge` and calling it from
  BOTH `_cmd_gauges` (default, NEXRAD/wet-only on) and `_cmd_vgauge` (unchanged
  behaviour: atlas with `--atlas`/`--gauges`, detail with `--detail`, NEXRAD opt-in)
  — `vgauge` stays the **replay/explicit-point** tool (reuse a store, no re-fetch).
  Parser opts via `_add_gauges_pipeline_opts(pg)` (`set_defaults(nexrad=True,
  wet_only=True, source="mrms", durations=[5,15,30,60])`). `gauges` is now a HEAVY
  command (Synoptic + MRMS + NEXRAD-volume fetch) by design.
- **Case-insensitive filename collision fix (2026-06-29).** Two distinct stations
  whose names differ only in case — HV's **`'Virginia City'` vs `'VIRGINIA CITY'`**
  — have `_safe_name` stems differing only in case, which **collapse to one file on
  macOS** (case-insensitive APFS/HFS+): one gauge's per-gauge CSV + detail figure
  silently overwrote the other (20 files for 21 atlas gauges → non-uniform). Fix:
  **`gauges.unique_safe_names(names)`** — deterministic (sorted) map suffixing the
  2nd+ of each case-insensitive group (`Virginia_City` → `Virginia_City_2`); used by
  `fetch_gauge_event` (store CSV write), `load_event_series` (read — same map so they
  agree), and `_virtual_gauge_products` (VG CSVs + detail + user-point figs). Only
  colliding names get a suffix, so non-colliding stems (and old stores) are
  unchanged. Verified: `unique_safe_names` case-insensitively unique + order-
  independent; the HV reuse atlas/detail now emit 21 files for 21 gauges.
- **Georeferenced export for GIS/CalTopo (`export.py`, CLI `export`, 2026-06-29):**
  re-packages a processed event — no radar re-run, modelled on `climate`/`zoom`
  (`--from-dir`/`--from-key` + AOI auto-match). Two products: **(1) EPSG:3857
  GeoTIFFs** (CalTopo's native Web-Mercator) of `--layers` (default
  **`anom_i15 i15max`**) — each as a **raw single-band float**
  `<key>_<field>_3857.tif` (`reproject_geotiff`, bilinear; nearest for categorical
  `tpki15`/`rqi`/`shsr`/`cbb`) **and** a **colorized RGBA** `_3857_rgb.tif`
  (`_colormap_rgba` → 4-band Byte, alpha tagged via `ColorInterp`; **NaN +
  dry<`wet_min` cells transparent** so it overlays a CalTopo basemap looking like
  the figure; i15→YlGnBu, anom→`cmc.vik` `TwoSlopeNorm(1)`, matching the maps).
  **(2) GeoPDFs** of the two primary figures — i15 map `<key>.pdf` + anomaly
  `<key>_anom_i{d}.pdf` (anom needs `<from-key>_anom_i{d}.tif` from `climate`).
  **`figure_to_geopdf(fig, ax, out, crs)`** rasterizes the full styled figure
  (savefig→PIL, **no `bbox_inches='tight'`** so the axes rect maps deterministically
  to page pixels), builds a page→CRS **affine from `ax.get_position()` +
  `get_xlim`/`get_ylim`** (spine convention, inversion-safe), and writes via GDAL's
  PDF driver with a **`NEATLINE`** = the map rectangle so the colorbar/title margins
  are *excluded* from georeferencing (not mis-placed). **Renders the PDF in UTM by
  default** (`--pdf-crs UTM` → identical look to the PNG deliverables; georef exact
  in whatever projected CRS) — `--pdf-crs EPSG:3857` to match the layers. **Needs
  GDAL's PDF driver** (`libgdal-pdf`, see env note); `geopdf_supported()` gates it
  and the GeoTIFF export + a note run without it. **HV 20260619 validated:** read-back
  GeoPDF embeds EPSG:32611, neatline→4326 spans the AOI (+0.04 clip margin); 3857
  tifs CRS-correct (i15max 0–93.7, anom 0–3.95); RGBA i15 79.7% painted (dry
  transparent), anom 100% (finite everywhere). Deliverables in the Box RD:
  `HiddenValley.pdf`, `HiddenValley_anom_i15.pdf`, `HiddenValley_{i15max,anom_i15}_3857{,_rgb}.tif`.
  **(3) Stream-network vector export (`export_streams`, CLI `export --streams`,
  2026-06-29):** the **full-resolution NHD network** for the AOI as a vector layer —
  reuses `refdata.streams` (NHDPlus HR flowlines, named creeks + unnamed headwaters,
  `watercourse_only`), **clips to the AOI polygon by default** (`gpd.clip`;
  `--streams-bbox` keeps whole flowlines over the bbox), writes EPSG:4326 (vectors
  reproject on import; GeoJSON is WGS84 by spec) via OGR — `--streams-format
  {geojson,gpkg,shp,kml}` (all four drivers present). `--streams-named-only` drops
  headwaters. Opt-in (network fetch, flaky endpoint degrades to empty). **HV
  validated:** clip → 1495 flowlines (206 named / 1289 unnamed), 100% within the AOI
  polygon; bbox → 1801 (whole flowlines). Deliverable `HiddenValley_streams.geojson`
  (3.5 MB) in the Box RD.
- **Wet-window pre-flight before NEXRAD (2026-08-14) — `mrms.wet_window()`.**
  A `--date` window spans ~30 h (`SCAN_PAD_H`) so the local day is covered, but
  handing that to the NEXRAD Level II fetch pulls **~300 volumes / ~2 GB**, of
  which a few hours are wet. `mrms.wet_window(aoi, start, end)` is the radar-side
  counterpart of `gauges.storm_window()`: it scans hourly `RadarOnly` QPE (a few
  KB per hour — cheaper than one Level II volume) and returns the tight
  `(start, end)` spanning the wet hours, `None` if dry. It spans **first→last**
  wet hour so a day with two cells returns one window covering both, and it
  opens an hour before the first wet stamp because hourly QPE at `HH` covers the
  hour *ending* at `HH`. Bug/Stallion 20260813: **30 h → 10 h**.
  **The bug it fixes:** `_cmd_vgauge` narrowed to the storm only `if start is
  None`, but `--date` sets `start` — so *every* `vgauge --date` run skipped the
  trim and fetched the whole day. `_cmd_gauges` narrowed correctly via
  `storm_window` but fell back to the **full span** when the gauges read dry,
  which is the same trap whenever gauges miss a cell the radar saw. Both now go
  through `cli._narrow_to_storm()`: explicit `--start`/`--end` always wins, then
  gauge rain-mass, then MRMS hourly QPE, then (only if all dry) the original
  window — never silently truncated to nothing. `nexrad.download_scans()` also
  warns above `BULK_SCAN_WARN` (120 volumes) so library callers see the cost
  coming. Tests in `tests/test_mrms.py` (bracketing, dry→None, two-cell span,
  clamping to the requested span, threshold).
- **Sorted output layout is the default (2026-08-14) — `stormscape/layout.py`.**
  A finished event is a few dozen files of four kinds, and one flat directory
  buries the figures. Products now sort into `figures/` (+ nested
  `figures/VirtualGaugeFigures/`), `rasters/`, `tables/`, `vectors/`;
  `RainGaugeData/` and the `nexrad_cache/` / `atlas14_cache/` stay at the event
  root (stores and inputs, not products — and the NEXRAD cache should stay
  obvious to delete). Writers go through `layout.out_path()` / `layout.subdir()`,
  readers through `layout.find()` / `layout.find_subdir()`.
  **Backward compatibility is the whole trick:** `find()` checks the sorted
  subdirectory then falls back to the flat path, so every pre-existing event
  folder keeps working through `--from-dir` / `--radar-dir` with **no
  migration** — verified by running `compare` against the same event laid out
  both ways and diffing the CSV (byte-identical). When a file is absent `find()`
  returns the *sorted* path so the error names where a fresh run would put it.
  Opt out with **`--flat`** on any writing subcommand, `layout="flat"` in the
  library, or `$STORMSCAPE_LAYOUT=flat`; `main()` exports the env var once
  rather than threading `layout=` through ~40 call sites, so `--flat` reaches
  library calls (`mrms.save_fields`, `gauges.fetch_gauge_event`) unchanged.
  `README.md` is pinned to the event root (it describes the folder; filing it
  under `tables/` with the CSVs buries it). Unmapped extensions (`.html` from
  `pick`) stay at the root. Tests in `tests/test_layout.py` (25) cover sorting,
  both read layouts, precedence, the env var, and a write→find round trip;
  `tests/test_export.py` now asserts on `export_geotiffs`' **returned** paths
  rather than hand-built ones.
- **Explicit analysis window for the stack (2026-08-15) — `--start`/`--end` on
  `i15`/`run`, `mrms.window_hours`.** `--date` scans a fixed ~30 h span
  (`SCAN_PAD_H = (4, 10)`, [day 04Z, next-day 10Z]) so the local day is covered.
  Back-to-back evening storms therefore **share** a storm-day window: on
  2026-08-14 the 04Z hour was the *tail of the 13 Aug storm* (6.0 mm areal max)
  and was stacked into "today's" peak-intensity maps. Verified against the real
  AOI — `--date` selects `[04Z, 23Z, 00Z, 01Z]`, the window selects
  `[23Z, 00Z, 01Z]`.
  `mrms.window_hours(date0, scan_pad_h, window)` centralises the choice;
  `i15_storm_day` and `multisensor_total` take `window=(start, end)` and make
  `date` optional. **`--start`/`--end` now scope the RADAR STACK, not just the
  gauges** — they previously existed on `run` (via `_add_gauge_opts`) and scoped
  only the gauge fetch while the radar quietly stacked the whole day, which is a
  trap; one pair of flags now means "the analysis window" everywhere. `i15` had
  no window flags at all and gets them via the new `_add_window_opts`. `--date`
  is no longer `required` on `i15`/`run` — give a date, a window, or both (the
  window wins); giving neither is rejected in `_stack_window`. `_event_label`
  supplies the `YYYYMMDD` for keys and titles from whichever was given.
- **`find()` climbs out of a layout subdirectory (2026-08-15).** Pointing
  `--from-dir` at `<event>/rasters` is the natural mistake — it is where the
  GeoTIFFs are — but the event AOI is in `<event>/vectors` and `RainGaugeData/`
  at the root, so the lookup missed them and `climate` silently degraded to
  framing on the i15 footprint instead of the event AOI. `layout._candidates`
  now appends the parent's sorted+flat pair when `basename(in_dir)` is in
  `RESERVED`, so `--from-dir <event>/rasters` behaves exactly like
  `--from-dir <event>`. Climbing is **only ever a fallback**: a local hit still
  wins, and an ordinary directory that merely happens to be named `rasters`
  gains nothing it would not otherwise find.
- **Hourly QPE is stamped at the END of its hour — the stack span was off by
  ~1 h (fixed 2026-08-15).** `QPE(HH)` is the rain in `[HH-1, HH]`. Measured,
  not assumed: the 2-min `PrecipRate` accumulation over `[HH-1, HH]` matches
  `QPE(HH)` at a **median ratio of 1.000**, while `[HH, HH+1]` gives 0.13–0.33
  (checked on two independent hours, 20260815-00Z and 20260813-21Z). So a run of
  wet stamps `[h0..hn]` describes rain over `[h0-1h, hn]`. `i15_storm_day` was
  stacking `[h0-14min, hn+1h]` — shifted a full hour late, so it **skipped up to
  46 min at the front of the first wet hour** and spent a fetch on a
  usually-dry trailing hour. The unread slice held **68.7%** of the first wet
  hour's rain on 20260813 (14.9 mm in a cell), 43.2% on 20260814, 14.9% on
  20260812. Now `t0 = h0 - 1h - 14min`, `t1 = hn` (the 14 min is the lead-in
  that makes the rolling i15 valid from the run's first wet minute).
  **No published number changed** — re-stacking 20260813 gives identical
  AOI values (i15max 86.60, i30max 62.41, i60max 41.67) and identical per-fire
  maxima; the 222 cells that do read higher all fall outside the perimeters.
  The bug predates the window work; making `wet_window`'s convention explicit is
  what surfaced the contradiction. Tests pin the span, and the module docstring
  now states the convention, because getting it backwards silently truncates the
  front of a storm rather than failing.
  **Gotcha for anyone testing this:** RQI/SHSR are read with `fetch` (singular),
  not `fetch_many` — mock both or the "offline" test hits S3 and CI's
  `-m "not network"` run breaks.
- **`max_wet_hours` warns when it truncates (2026-08-15) — `find_wet_hours`.**
  `MAX_WET_HRS = 8` caps the processed wet hours, but it ranks them by
  **intensity**, so what it discards is a long storm's weakest hours — normally
  its opening and closing tails. Those hours also *bound the stacked span*:
  dropping a trailing wet stamp shortens the contiguous run, so `total` loses
  that rain outright and the rolling i15/i30/i60 never see it. Nothing
  downstream could tell — `n_wet_hr` reports what was **kept**, so a truncated
  run looks like a complete one.
  **Found by consumers, not by tests.** Stacking the three Bug/Stallion storms
  (12–14 Aug 2026) into one composite made `max(3 individual storms) ==
  composite` checkable for the first time, and it failed: the composite was
  higher at 684 cells. Cause — the **published 13 Aug analysis has 9 wet hours,
  one over the cap, and silently lost its 04Z tail** (6.0 mm). Impact there was
  small (all per-fire maxima unchanged; Stallion i15 areal mean 11.58 → 11.84
  mm/h, total 5.73 → 5.89 mm) but it was invisible. On the 3-day window the cap
  drops **12 of 20** wet hours. Single-storm runs can never surface this — only
  a multi-storm stack cross-checked against its parts.
  Now warns naming the dropped stamps, the weakest kept vs strongest dropped
  qmax, and `--max-wet-hours`. **Warn, don't refuse** (same call as
  `nexrad.BULK_SCAN_WARN`): the cap is a real cost control and some callers do
  want the most intense hours only. Tests cover the truncating case, the dropped
  hour being the *weakest* not the last, the shortened span, and — importantly —
  that the sub-cap and all-dry-fallback paths stay silent.
- **Near-real-time burn severity (`burn.py`, CLI `burn`, 2026-08-14) — CIMSS
  BRISK.** Puts a scar under the rain *while the fire is still burning*, which the
  authoritative products cannot: BAER soil burn severity lands days-to-weeks after
  containment and only for assessed fires, MTBS a year+. BRISK is a Google Earth
  Engine **dNBR data-fusion composite over nine satellites** (GOES-E/W ABI, SNPP /
  NOAA-20 / NOAA-21 VIIRS, Landsat 8/9, Sentinel-2a/b), mapping every large
  (>~5,000 acre) US fire **daily**.
  **Access — the portal is a decoy.** `cimss.ssec.wisc.edu/brisk` is a RealEarth
  viewer whose WMTS (`re-brisk.ssec.wisc.edu/wmts/BurnScars-dNBR.xml`) serves
  **rendered PNG tiles, capped at zoom 7** — pretty, not data. The raw field is in
  the open Apache-indexed archive behind it, **one GeoTIFF per fire per day**:
  `bin.ssec.wisc.edu/pub/realearth/brisk/<year>/<Fire>-<ST>-dNBR_<YYYYMMDD>_235959.tif`
  (2025→present; 8k+ scenes / ~590 fires in 2026 alone, current through *today*).
  So `burn.py` scrapes the directory index into a cached catalog and downloads
  only the scenes that intersect the AOI.
  **AOI screening is header-only.** GDAL range-reads a remote GeoTIFF header via
  `/vsicurl/` in ~0.02 s across a 12-thread pool — a whole day's ~60 fires screen
  in **1.4 s without downloading a pixel**; footprints are memoised in
  `brisk_cache/bounds_dnbr.json`, so the first full-archive screen (~650 scenes,
  ~25 s) makes every later one **0.04 s**. Index cache TTL is **6 h for the
  current year and infinite for past years** (closed years never gain scenes),
  and a dead server falls back to the stale index instead of failing.
  **Gotchas, all measured:**
  (a) scenes are **EPSG:3857 at "60 m", which is not 60 m of ground** — Web
  Mercator metres shrink with latitude, so the cell is **46.5 m at 39°N**; areas
  computed on the raw transform run ~66% high there, hence the `cos(lat)`
  correction in `cli._burn_class_table`;
  (b) every scene tested lands on an **exact 60 m multiple**, so mosaicking
  neighbouring fires is a paste, not a warp;
  (c) **NaN marks outside-the-burn but the files tag no nodata**, so
  `masked=True` masks nothing — test `np.isfinite` (typically only 10-19% of a
  scene is valid);
  (d) dNBR is **unscaled** (~-0.3 to 1.0), not the ×1000 integer form the
  MTBS/USGS thresholds are quoted in — `SEVERITY_SCHEMES["usgs"]` is
  0.10/0.27/0.44/0.66, `["brisk"]` the portal's own 0.10/0.40/0.70;
  (e) a fire's footprint **grows day to day**, so `find_scenes` takes the latest
  scene per fire **on or before `--date`** — a storm-day map must not include
  severity mapped after the storm.
  **Mosaic rule: NaN-aware `np.fmax`.** Scenes are mostly NaN outside their own
  fire, so a plain "first/last wins" merge blanks whichever fire is written second
  wherever footprints overlap; `fmax` combines, is order-independent, and takes
  the more severe value where two fires genuinely overlap. Each scene is
  `reproject`d into one destination grid, so the same code path handles the
  BAER scenes (which are **per-fire UTM at 20 m**, not 3857).
  **BAER SBS (`--product sbs`, `baer-data/`)** is the authoritative soil product
  but **sparse** (8 fires in 2026 vs BRISK's ~590) and classified uint8. Its
  embedded palette colours 1-4 with BRISK's own four severity colours and paints
  **0 and 5+ the same black** — i.e. the product itself treats anything outside
  1-4 as a mask (water/inholding/unmapped, ~3% of the NV scene checked), so
  `SBS_VALID = (1, 4)` reads the rest as missing rather than charting it as a
  class. Categorical mosaics use **nearest** (averaging class 1 and 3 into 2
  would be a fabricated severity).
  **Display — the BAER palette is the default, and it IS the BRISK palette
  (2026-08-14).** Measured, not assumed: **all 77** of the 2025 BAER
  soil-burn-severity rasters — written by many different BAER teams, in ERDAS
  Imagine — embed an *identical* class palette, `1 (0,128,128) teal / 2
  (82,204,204) cyan / 3 (255,232,32) yellow / 4 (168,0,0) dark red` (only the
  class-5 mask colour varies: black/gray/white). Those are exactly the four
  colours in BRISK's `qgis_BRISK_dNBR_colorscale_v2.txt`, so the portal and the
  BAER deliverables already share one scheme — `BAER_CLASS_COLORS` /
  `BAER_ANCHORS` (`BRISK_ANCHORS` is an alias; `register_brisk_cmap` →
  `register_baer_cmap`). Maps are **classed by default**, the way BAER
  publishes: `burn.severity_colors(scheme)` returns
  `(ListedColormap, BoundaryNorm, ticks, labels)` so the field is banded at the
  severity breaks and the **colour bar is labelled with class names, not dNBR
  numbers**. A 4-class scheme uses the official table exactly; the 5-class
  `usgs` scheme has no official 5th colour so its colours are *sampled from the
  same ramp* at class midpoints (documented, not smuggled). `--continuous` gives
  the smooth ramp (same colours), an explicit `--cmap` opts out entirely.
  This required `plot.drape_i15` to gain **`norm=` / `cbar_ticks=` /
  `cbar_ticklabels=`** (a norm replaces the linear `vmin=0,vmax=` scale) — small
  and general, usable by any future classed map. `_burn_display_defaults`
  resolves the cut/scale **per product** — dNBR 0.10/1.0, SBS 1.5/4.0 — because
  one shared default would clip the class map to its lowest class or paint
  unburned ground. The drape reuses `drape_i15` (so `--reference`/`--clip`/north
  arrow/ticks all just work); `--alpha` is worth raising from the project 0.32
  when the scar is the subject.
  **BAER also publishes its own dNBR, and BRISK matches it (validated
  2026-08-14).** `baer-data/<year>/` holds `<Fire>-<ST>-prelim-dNBR_<date>_*.tif`
  (**118 for 2025**, none yet for 2026) alongside the `-sbs_` rasters — note the
  **`-prelim-` infix**, which `parse_name` strips. They are **int16 ×1000**
  (the BARC convention; NOAA's own `BARC256 = dNBR×5 − 275` identity confirms the
  scaling), **ESRI:102039** Albers, **20 m where the source was Sentinel-2 and
  30 m where it was Landsat**. **39 fires have both products.** Method: divide
  BAER by 1000 and `Resampling.average` it **down onto the BRISK grid** —
  aggregating the finer product to the coarser support rather than upsampling
  the thing under test — then score only **burned cells (BAER dNBR ≥ 0.1)**,
  because the huge unburned surround agreeing near zero inflates r by ~0.03–0.1.
  **Result: BRISK computes the same dNBR.** Pooled over 8 well-matched fires,
  **r = 0.938 on 1.24 M cells, slope 1.089**; across all 39, median **slope
  1.013, bias +0.003, ratio 1.011** — unity slope, no bias. **The disagreements
  are compositing latency, not a different algorithm**: every poor performer
  (middle-mesa-nm r 0.03, turkeyfeather-nm 0.12, blind-az 0.22, derby-co 0.27,
  dillon-ca 0.58, laguna-nm 0.60, island-creek-id 0.61) recovers to **r 0.80–0.96
  when given a BRISK scene 5–21 days later**, and every one of them read *low* on
  the BAER date — the composite had not yet ingested a clear post-fire overpass.
  A pre-set (not per-fire-tuned) **+14 d** rule lifts the distribution: median r
  0.874→0.913, r≥0.90 **14→24** of 39, r<0.60 **5→1**, IQR 0.765–0.934 → 0.866–0.932.
  Note it rescues the *tail* rather than improving the median fire (21 better,
  18 slightly worse) — a mature scar keeps darkening away from the BAER snapshot,
  which is why the +14 d slope rises to ~1.05–1.09. **BAER 20 m (Sentinel-2)
  agrees better than 30 m (Landsat)** — median r 0.923 vs 0.814, Mann-Whitney
  **p = 0.0005** — most likely because BRISK's composite and the Sentinel-based
  BAER share the same underlying acquisition. **Operational rule: trust BRISK's
  *pattern* immediately, but give the composite ~2 weeks before trusting its
  *magnitude*** — early scenes under-read. Study script + CSVs live in the
  repo as **`examples/brisk_vs_baer.py`** (reruns the whole study from the two
  catalogs; writes CSV + figure to `--out-dir`, default `baer_study/`, which is
  gitignored). The BAER dNBR is reachable as **`--product baer_dnbr`**
  (`PRODUCTS[...]["scale"] = 1000.0` divides it back to a plain index on read).
  **Acted on in the tool:** `burn.MATURITY_DAYS = 14`; `find_scenes` adds an
  **`age_days`** column (days since the fire's *first* appearance in the
  catalog, computed **before** any `--date`/`--since` trim, or trimming early
  scenes would make an old fire look new), `--list` prints it and stars
  immature composites, `burn_severity` prints an advisory, and **`--min-age
  DAYS`** hard-filters. The maturity screen runs **after** the AOI
  intersection -- filtering the national catalog first named 593 irrelevant
  fires in the drop message. Off by default: an immature scar still beats no
  scar, provided the run says so.
  **`dnbr`/`severity` are deliberately NOT in `plot._MASK_DRY`** — that mask cuts
  below 0.5, which on a 0-1 dNBR field would erase everything short of high
  severity, and it also drives smoothing, which a severity field does not want.
  **THE CAVEAT THAT MATTERS FOR THIS GROUP'S SCIENCE: dNBR is a *vegetation*
  index, not soil burn severity.** The USGS post-fire debris-flow models are
  calibrated on **soil** burn severity — dNBR adjusted by field crews for
  hydrophobicity, ground cover and duff consumption. BRISK is explicitly
  **interim**: act on it early, then supersede with `--product sbs`.
  Outputs are `mrms`-style result dicts (`save_fields` / `drape_i15` just work):
  `<key>_dnbr.tif`, `<key>_severity.tif`, `<key>_burn_classes.csv` (pixels + true
  km² + fraction per class), `<key>_burn_scenes.geojson` (provenance: which fire,
  which date), `<key>_burn.png`. Cache is `brisk_cache/` (added to
  `layout.RESERVED`, matches `nexrad_cache/`; keep out of git). **No new deps.**
  Validated live on Ward NV 20260814 (91.2 km² burned, 6.3% high severity),
  a 15-fire central-Oregon mosaic, and Cottonwood-Peak NV SBS.
  Tests: `tests/test_burn.py` (48, offline — synthetic listings + local GeoTIFFs).

- **Abandoned mine features (`mines.py`, CLI `mines`, 2026-08-15) — USGS USMIN.**
  Puts historic mining under the rain: in the Great Basin the steep catchments
  that produce debris flows are full of it, and **mine dumps and tailings** are
  loose, often contaminated material sitting on the channel network.
  **The headline constraint, which drove the whole design: there is NO public
  point-level abandoned-mine HAZARD database, by policy.** USGS **FS 2025-3003**
  (2025) says the national abandoned-mine-feature database being built under
  USMIN "will not publish specific location information of any abandoned mine
  workings, and the detailed national abandoned mine feature database will not be
  publicly available" — the locations could be used to enter hazardous workings or
  vandalise historic structures — and only aggregated derivatives (per county,
  per watershed) are planned. **Nevada matches:** NDOM's operational AML layers
  are real and findable (org `CXYUMoYknZtf5Qr3`; the internal *AML Field Map*
  web map references `NVPoints`, `NVSites`, `InternFieldDataCaptureLayers`
  (hazards / revisit-securing / non-hazards), `AMLCracForms`) but every one
  answers **`499 Token Required`** anonymously — their public AML items are a
  hazard-reporting form and a "Stay Out and Stay Alive" game. Their open-data
  portal's "AML" entry is a **Hub Page**, not a layer. So NDOM needs a
  credentialed request, not an endpoint. Public Nevada fallback: **NBMG OFR
  2001-03** (100k+ sites, shapefile, UTM 11 **NAD27**, ~115 MB, $30, compiled with
  BLM NV + NDOM, 2001, self-described preliminary).
  **What we use instead: USMIN "Prospect- and Mine-Related Features from USGS 7.5-
  and 15-minute topographic quadrangle maps"** — public precisely *because* it is
  digitised from already-published topo sheets, so it reveals nothing the printed
  maps did not. `energy.usgs.gov/arcgis/rest/services/Hosted/USMin_Prospect_and_
  mine_related_map_features/FeatureServer`, **layer 17 = points, 18 = polygons**,
  EPSG:4326. All 11 western states covered (NV 121,193 · CA 61,347 · AZ 43,741 ·
  CO 43,297 · MT 22,335 · NM 18,709 · UT 18,471 · ID 16,635 · OR 14,775 ·
  WA 11,016 · WY 7,690). Read it as a **historical map compilation, not a hazard
  inventory** — `topo_date` runs 1950-1994, no hazard ranking, no securing status.
  **THE GOTCHA THAT WOULD HAVE SHIPPED A BROKEN DEFAULT: the waste is in the
  POLYGON layer.** USMIN splits features by how the symbol was drawn, and dumps
  and tailings were nearly always drawn as an extent — **14,815 polygons vs 413
  points nationally; 778 vs 4 in Nevada**. So `geometry="points"` returns ~nothing
  for a `waste` query and looks exactly like an AOI with no mining.
  **`geometry="both"` is the default** for that reason (`_as_points` takes
  `representative_point()` for polygons when counting/labelling — guaranteed
  inside, unlike a crescent dump's centroid).
  **Paging bug found + fixed in `refdata` (latent, would have bitten this
  module):** ArcGIS reports truncation as `exceededTransferLimit` **at the top
  level on a MapServer** (NHD, TIGER, GNIS) but **only under `properties` on a
  hosted FeatureServer** (USMIN). `_query` checked the top level only, so USMIN
  returned exactly **2000 of 3134** features over the HV AOI with no error.
  `refdata._more_pages` now checks both; `_query` was renamed **`arcgis_query`**
  (public, `_query` kept as an alias) and gained `token=`, `paginates=` and
  `what=`. The existing layers were *unaffected* (all MapServers — verified by
  paging 26,037 NHD flowlines correctly), so nothing published changed.
  **Second silent-failure fix:** an ArcGIS error is **HTTP 200 with an
  `{"error": ...}` body**, which read as "no features here". It now warns *and*
  prints to stderr — printing because `refdata` installs a module-level
  `warnings.filterwarnings("ignore")` that would swallow the warning. (That
  global filter is a pre-existing wart: it silences warnings process-wide for
  anything imported after it.) This is how the **NBMG USMIN mirror**
  (`gisweb.unr.edu/nbmg/.../USMIN/MapServer`, v4.0, NV only) went unnoticed —
  it answers **"Pagination is not supported"** to any `resultOffset`, hence the
  `paginates=False` flag in its source spec. Prefer the national service.
  **Also `ftr_name` is `''`, not null** (2,990 of 3,134 over HV), so `.notna()`
  calls every feature named — `_blank_to_na` fixes it.
  **Grouping:** 55 national feature types collapse to six groups (`waste`,
  `openings`, `surface`, `aggregate`, `prospect`, `other`) via an exact map plus
  **prefix families** (`tailings*`→waste, `quarry*`→surface) so a subtype USMIN
  adds later lands correctly instead of falling to `other`. **Default
  `DEFAULT_KINDS = ("waste", "openings")`** — two-thirds of features over a
  Nevada AOI are prospect pits, which bury the rainfall field. The kind filter is
  **pushed into the service `where`** when expressible (`other` and "everything"
  cannot be enumerated → `None`, filter locally) and **always re-applied locally**
  so a source that ignores `where` still honours `kinds`.
  **Density:** `density_grid` (vector, → cell centres with `count`/`per_km2`) and
  `density_raster` (→ mrms-style dict, so `save_fields`/`drape_i15` just work) both
  bin in **EPSG:5070 (equal-area)** on a grid **anchored to the projection origin**,
  so a cell is the same km² everywhere and the same ground bins identically run to
  run. Both paths are cross-checked in tests to give identical totals and maxima.
  `plot.add_mines` draws points (polygons as real footprints) or graduated symbols
  with **marker area ∝ count** plus a size key — a graduated symbol is unreadable
  without one. `MINE_DENSITY_SWITCH = 400` drives `mode="auto"`.
  **Pluggable sources by design:** `SOURCES` + `register_source()`; `ndom` is
  already registered with its URLs and a `$STORMSCAPE_NDOM_TOKEN` hook, and its
  field mapping uses **tolerant candidate lists** (`_pick` returns an empty column
  rather than silently mislabelling one) because those column names are unverified
  guesses until access lands. Tokens are read from arg or env, **never logged**.
  **CLI:** `mines` writes `<key>_mines.geojson`, `_mine_classes.csv`,
  `_mine_density.tif`, `_mines.png`; `--mines` is an overlay flag on **map, run,
  nexrad, zoom, burn, export** (the commands that render through `drape_i15` —
  `climate`/`smooth` use `_draw_field`, and the `compare` residual maps are left
  clean deliberately). `_mine_kwargs(args)` uses `getattr` so it is safe to splat
  everywhere. Defaults `cmap="YlOrBr"`, `wet_min=0.5` (counts, not rainfall).
  `drape_i15` gained overlay params `mines`/`mines_mode`/`mines_kinds`/
  `mines_cell_km`/`mines_groups`/`mines_labels`; `mines=True` auto-fetches like
  `reference=True`. **`drape_i15` now also accepts `i15=None`** → terrain +
  overlays, no field and no colour bar (it raises only if hillshade *and* field
  are both None). That exists because the first `mines` figures draped the
  density raster *and* drew the features: blocky 1 km cells under graduated
  symbols of the same counts, i.e. one quantity encoded twice. The figure now
  carries the features alone; the raster stays a GIS product, and
  **`--density-map`** drapes it on request. Useful for any future
  vectors-over-terrain map.
  **Bug caught only by rendering, not by the unit tests:** `_cmd_mines` read
  `args.vmax`, which only `burn`'s parser defines — so the command fetched, wrote
  three files, then died at the last step. The CLI tests patch the command out, and
  the real path needs network, so neither saw it. Fixed by adding `--vmax` *and* by
  `test_cmd_mines_only_reads_args_the_parser_defines`, which **AST-walks
  `_cmd_mines` for `args.<name>` and checks each against the parsed namespace**
  (39 attributes). Worth copying for other commands.
  **No new deps.** Tests: `tests/test_mines.py` (75, offline — `requests.get`
  replaced by a replaying double).

- **Testing (`tests/`, pytest, added 2026-07-29) — run it before every push.**
  `pytest` (or `/opt/anaconda3/envs/GISMan/bin/python -m pytest`) — **480 tests,
  ~2 s, entirely offline** (no MRMS/NEXRAD/Synoptic/3DEP/Atlas 14/USMIN request,
  no token), so it is cheap enough to run constantly. Config lives in
  `pyproject.toml` (`[tool.pytest.ini_options]`, `--strict-markers`); install the
  deps with `pip install -e ".[test]"`. Tests encode the **documented invariants**
  from this file, not implementation details: the i15 estimator round-trip, the
  smoothing guarantees (radius-0 identity, no NaN zero-bleed, monotone peak
  reduction, δ-spike symmetry), sub-1-yr recurrence continuity, the precip-bearing
  cadence rule, case-only filename collisions, mass-weighted storm windows, the
  RQI/cadence screens (and that cadence must *not* filter `total`), the
  3857/RGBA-transparency export contract, and that all 15 subcommands still parse.
  Markers: `optional_deps` (Py-ART/wradlib/GDAL-PDF — self-skipping) and `network`
  (deselected in CI via `-m "not network"`). **Fixtures over live data:**
  `conftest.py` factories write synthetic GeoTIFFs/gauge frames to `tmp_path`; don't
  reach into the Box event folder from a test.
  **CI** (`.github/workflows/ci.yml`) runs the suite on **Python 3.10 + 3.13** on
  every push/PR via **plain pip** (`pip install -e ".[test]"`), which also guards
  two contracts: that `pip install stormscape` actually resolves, and that
  `import stormscape` + every CLI parser works with **none** of the optional radar
  stack installed. Two packaging bugs were found precisely this way — see below.
- **Packaging gotchas (found by testing `pip install` in a clean venv, 2026-07-29):**
  (1) in `pyproject.toml` the `dependencies` key must come **before** any
  sub-table like `[project.urls]`, or TOML parses it as `project.urls.dependencies`
  and the build dies; (2) under **PEP 639** a `license = "MIT"` expression and a
  `"License :: OSI Approved :: ..."` classifier **cannot coexist** — modern
  setuptools errors out, so the classifier is gone. Also `scipy` was missing from
  the core deps even though `smoothing.py` imports it at module level (so a
  no-extras `pip install` broke `import stormscape`); it is now required.
  **Always validate a packaging change with a throwaway venv**, not just the
  already-populated GISMan env, which masks missing dependencies.
- **`py3dep` is imported lazily** (`dem._py3dep()`, 2026-07-29) even though it is
  a required dependency — it drags in ~35 HyRiver packages, and deferring it keeps
  `import stormscape` fast and turns a broken HyRiver install into a clear
  "install py3dep" error instead of an unusable package. Same pattern as
  `nexrad`'s `arm_pyart`/`nexradaws` and `export`'s `osgeo`.
- **Verification is also visual** — the pytest suite covers the math and the
  plumbing, but figures still need looking at. Render them and inspect (FlowAlert
  fixtures for parsing parity; a live Hidden Valley run for end-to-end). Keep
  outputs (`*.tif`, large rasters, caches) out of git.

## Provenance

Earlier development notes (i15-estimator validation, radar-era storm analysis)
live in the originating post-fire debris-flow project; **this repo is the clean,
generic, canonical home** — it was renamed from `i15toolkit` to `stormscape`
(2026-06-29) once the scope grew past 15-minute intensities to cover single-radar
NEXRAD, gauges, climatology, smoothing, and GIS export. A stale local working
copy was deleted (2026-06-25) to avoid version confusion; the original shared
MRMS code (D. Cavagna's `MRMS_stack.py`) lives outside this repo.

Note the `i15`/`i30`/`i60` names throughout the code are the **rainfall metrics**
(peak 15/30/60-minute intensity) and are deliberately unchanged — only the
package name changed.
