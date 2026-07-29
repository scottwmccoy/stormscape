# stormscape

Terrain + radar-rainfall mapping from public US datasets, in one small package.

Given an **area of interest** (a bounding box, a shapefile/GeoJSON, or a
shapely geometry) and, for rainfall, a **storm-day date**, it will:

1. **Download a DEM** from USGS **3DEP / The National Map** at 1 m / 10 m /
   30 m and render a **hillshade** (`stormscape.dem`);
2. **Stack 2-minute NOAA MRMS `PrecipRate` radar returns** into maps of the
   **peak 15/30/60-minute rainfall intensity** (`i15`/`i30`/`i60`, mm/h) over
   the storm, plus companion fields — total accumulation, peak 2-min rate, time
   of peak, radar quality (RQI), beam height (`stormscape.mrms`);
3. **Pull ground rain-gauge rainfall** for the same AOI + storm from the
   **Synoptic / MesoWest** API and reduce each gauge to the same metrics —
   storm total and peak 15/30/60-min intensity (`stormscape.gauges`);
4. **Drape the radar field over the hillshade** with optional vector overlays,
   overlay the gauges, and **compare radar QPE against the gauges** — per-gauge
   residuals + skill stats, optionally screened by radar quality and gauge
   reporting cadence (`stormscape.plot`, `stormscape.compare`);
5. **Reach the raw single-radar NEXRAD Level II volumes** for the radar nearest
   the AOI — reflectivity / velocity at each elevation tilt, gridded over the AOI
   or sampled at the gauges (the underlying radar behind the MRMS mosaic;
   `stormscape.nexrad`).

It was lifted out of a post-fire debris-flow study and generalized: nothing
here is tied to that project's fire perimeters or event inventory — you supply
the AOI and date.

---

## Install

The geospatial stack reads best from **conda-forge** (notably `libgdal-grib`,
which gives rasterio's GDAL the GRIB driver MRMS needs):

```bash
conda env create -f environment.yml      # creates env "stormscape"
conda activate stormscape
pip install -e .                          # installs the package + CLI
```

Already have a working geospatial env (py3dep + rioxarray + rasterio-with-GRIB
+ geopandas)? Just install into it — no need for a new environment:

```bash
conda activate <your-env>
pip install -e .
```

Quick check that GRIB is available (needed for MRMS):

```python
import rasterio
assert "GRIB" in rasterio.drivers.raster_driver_extensions().values()
```

The [`export`](#georeferenced-export-for-gis--caltopo) command's GeoPDF output
additionally needs GDAL's PDF driver — `conda install -c conda-forge
libgdal-pdf` (optional; the EPSG:3857 GeoTIFF export and everything else work
without it).

---

## Command line

```bash
# everything at once: DEM -> hillshade -> i15 -> draped figure
python -m stormscape run \
    --bbox -105.55 40.55 -105.25 40.80 \
    --date 2021-07-20 --resolution 10 \
    --out-dir ./out --key cameron_peak

# or step by step
python -m stormscape dem --aoi burn_perimeter.geojson --resolution 10 --out-dir ./out --key fire
python -m stormscape i15 --aoi burn_perimeter.geojson --date 20210720      --out-dir ./out --key fire
python -m stormscape map --hillshade out/fire_hillshade.tif --i15 out/fire_i15max.tif \
    --perimeters burn_perimeter.geojson --out out/fire.png
```

`--bbox W S E N` is lon/lat degrees; `--aoi` is any vector file GeoPandas can
read. `run`/`map` accept overlays: `--perimeters` (bordered), `--basins` (thin
outlines), `--highlight` (bold cyan), `--points` (triangles).

**Full flag reference for every subcommand:** [`docs/cli.md`](docs/cli.md) (or run
`python -m stormscape <command> --help`).

### Basemap for orientation (named creeks/rivers, roads, place names)

Add `--basemap` to `run` or `map` to underlay an open-source basemap
(downloaded as tiles via [contextily](https://github.com/geopandas/contextily))
so you can tell *where* the storm fell:

```bash
python -m stormscape map --hillshade out/event_hillshade.tif \
    --i15 out/event_i15max.tif --perimeters aoi.kmz --basemap \
    --out out/event_basemap.png
```

The default provider is **`USGS.USTopo`** — the USGS National Map
topographic basemap (public-domain USGS tiles with **named creeks/rivers**,
roads, contours, and place names; the same National Map service the USGS
`pfdf` package draws from). Because those tiles are already labelled, the
hillshade is turned off so the labels read through the rain (the i15 drape uses
the project-wide default `--alpha 0.32`). Tune with `--basemap-provider`
(`USGS.USImageryTopo` for imagery + labels, `OpenStreetMap.Mapnik` for dense
roads/waterways) and `--basemap-zoom` (default auto). For a label-free base
(e.g. `CartoDB.VoyagerNoLabels`), add `--basemap-labels
CartoDB.PositronOnlyLabels` to float names on top instead.

Needs contextily: `conda install -c conda-forge contextily` (or
`pip install -e ".[basemap]"`).

### Vector reference overlays — streams / roads / place names (recommended)

For crisp, labelled features drawn directly on the hillshade (no tiles), add
`--reference`:

```bash
python -m stormscape map --hillshade out/event_hillshade.tif \
    --i15 out/event_i15max.tif --perimeters aoi.kmz --reference \
    --out out/event_reference.png
```

This fetches AOI-scoped vectors from public ArcGIS REST services and overlays
them with labels:

| layer | source | field |
|---|---|---|
| streams (named) | USGS **NHDPlus HR** flowlines | `gnis_name` |
| roads (primary + secondary) | US Census **TIGER/Line** | `NAME` |
| place names | USGS **GNIS** | `gaz_name` |

These are the same National Map / NHD data the USGS `pfdf` package uses, but
queried per bounding box (small, fast) rather than downloading whole
hydrologic-unit bundles. Named mainstems are drawn at normal weight and labelled; the smaller *unnamed*
headwater streams are drawn as finer, lighter lines so the dense NHD HR network
adds drainage texture without burying the i15 field. Add `--local-roads` for
residential streets, `--no-reference-labels` for lines/points without text.
Add `--clip` to crop the figure tightly to the AOI / `--perimeters` extent
(`--clip-margin` sets the padding fraction) and `--dpi` to set the export
resolution (default 200). The hillshade is `hillshade_vmin`/`hillshade_vmax`
in the API. (Note: the DEM-masking flag on `dem`/`run` is `--clip-dem`, kept
distinct from the figure-cropping `--clip`.)
**No extra dependency** — uses `requests` + `geopandas`. Fetch the layers
yourself with `stormscape.refdata.streams/roads/places(aoi)` (each returns a
GeoDataFrame with a tidy `name` column) and pass them to `drape_i15` via the
`streams=`/`roads=`/`places=` arguments.

---

### Rain gauges & radar-vs-gauge comparison

Pull ground rain-gauge data for the AOI + storm and compare it against the radar
field. Gauges come from the **Synoptic / MesoWest** API, which needs a free
token (academic / research "open access" is free — request it at
[synopticdata.com](https://synopticdata.com)). Set it once:

```bash
export SYNOPTIC_TOKEN=xxxxxxxx          # or pass --token on any gauge command
```

```bash
# 1) full gauge pipeline in ONE Synoptic draw:
#    (a) canonical store (reused by compare / recurrence / vgauge):
#        <key>_gauges.geojson (coords + total_mm, i15/i30/i60_mmph, report_min =
#        native precip cadence, i15_peak_time) + RainGaugeData/ series CSVs;
#    (b) virtual-gauge rainfall comparison atlas (<key>_vg_atlas.png), MRMS +
#        NEXRAD + real (3-way) at every wet near-AOI station; and
#    (c) per-gauge detail figures -> VirtualGaugeFigures/.
#    The store keeps the whole storm day; (b)/(c) clip to the rain window + use the
#    same wet-gauge set (so the atlas and detail figures are uniform).
python -m stormscape gauges --aoi aoi.kmz --date 20260619 \
    --out-dir ./out --key event              # --store-only / --no-detail / --no-nexrad to trim

# 2) sample the radar rasters at each gauge -> residuals + skill stats + a CSV
python -m stormscape compare --gauges out/event_gauges.geojson \
    --radar-dir ./out --key event --out out/event_compare.csv \
    --map out/event_resid.png             # optional radar-minus-gauge map

# or do it all in one shot (DEM -> i15 -> gauges -> figure -> comparison)
python -m stormscape run --aoi aoi.kmz --date 20260619 --resolution 10 \
    --out-dir ./out --key event --gauges --compare --multisensor
```

`compare` reports, per metric (storm total + i15/i30/i60), the **bias, RMSE,
MAE, Pearson correlation, and mass ratio** (Σradar / Σgauge). Two screens keep
the comparison honest:

- `--rqi-min 0.8` drops gauges where the radar beam is unreliable (low RQI at
  the gauge cell);
- `--max-report-min 15` drops gauges whose **native precip reporting interval**
  exceeds 15 min — for the *sub-hourly* metrics only. A coarse (hourly) reporter
  smears bursts under the 1-min interpolation and reads an artificially low i15;
  storm total is cadence-insensitive and always uses every gauge.

`--multisensor` adds a gauge-corrected **MRMS MultiSensor QPE** total row, to
separate radar-only QPE bias from the gauge-corrected product. (Gauge data uses
Synoptic's *Basic Precipitation Service* — the `precip=1` flag on the Time
Series API; the dedicated Precipitation Service is an Enterprise add-on and is
not used.)

On the figure the gauges are filled circles — coloured on the **i15 scale** (so
a gauge reads like the radar beneath it) on the main `run` map, or on a
diverging **radar − gauge** scale on the `compare --map` residual map.

---

### Single-radar NEXRAD Level II (raw tilts)

`run` / `i15` use the gridded MRMS *mosaic* — already QC'd and blended across
radars. To reach the **raw single-radar** archive — the WSR-88D **Level II**
volumes, with every elevation tilt — use `nexrad`. It finds the radar nearest
the AOI, pulls the volume scans for your time / window from AWS, grids the
lowest tilt over the AOI, and drapes it:

```bash
# storm-peak reflectivity (per-cell max over the window), nearest radar
python -m stormscape nexrad --aoi aoi.kmz --composite \
    --start 202606192000 --end 202606200200 \
    --hillshade out/event_hillshade.tif --out-dir ./out --key event

# a single scan nearest a time
python -m stormscape nexrad --aoi aoi.kmz --date 20260619 --time 2230 \
    --hillshade out/event_hillshade.tif --out-dir ./out --key event
```

Writes `<key>_refl.tif` (or `<key>_reflmax.tif` for `--composite`) in EPSG:4326
plus a draped `<key>_nexrad.png`. Choose the radar explicitly with `--radar
KRGX`, another tilt with `--sweep 1`, another moment with `--field velocity`, or
the grid step with `--res-m`. Volumes are cached under `<out-dir>/nexrad_cache/`.

The figure takes the **same context flags as `run`/`map`** — `--reference`
(labelled NHD streams / TIGER roads / GNIS places), `--perimeters`, `--clip`,
`--basemap` — plus `--gauges` (live Synoptic fetch; needs `$SYNOPTIC_TOKEN`) or
`--gauges-file gauges.geojson` (a precomputed file, no token) to overlay gauges
**coloured by the radar value sampled at each** on the field's dBZ scale, so a
reflectivity map reads like the rainfall maps.

Transport is [**nexradaws**](https://github.com/aarande/nexradaws) (v2+, the
current `unidata-nexrad-level2` S3 bucket — the older `noaa-nexrad-level2`
Big-Data bucket it replaced was deprecated); reading the volumes uses
[**Py-ART**](https://arm-doe.github.io/pyart/). Both are optional extras:

```bash
conda install -c conda-forge arm_pyart      # reader (binary stack)
pip install -e ".[nexrad]"                  # + nexradaws (pip-only transport)
```

The Level II archive on AWS reaches back to the 1990s — *further than* MRMS's
2020-on 2-min cadence — so `nexrad` is also how you look at **older** events.

**Pre-2020 i15/i30/i60 (the MRMS analogue).** Add `--intensity` to turn the
Level II volumes into the *same* peak-intensity fields the MRMS engine produces —
`<key>_i15max.tif`, `_i30max`, `_i60max`, `_total_mm`, `_peakrate_mmph` — so
single-radar and MRMS results are directly comparable and you can reach storms
that pre-date MRMS:

```bash
python -m stormscape nexrad --aoi aoi.kmz --intensity \
    --start 202606192000 --end 202606200200 \
    --hillshade out/event_hillshade.tif --out-dir ./out --key event_l2
```

Each volume's lowest tilt (all SAILS low cuts → ~3-min effective cadence) becomes
a rain rate via a **capped convective Z–R** (`Z = a R^b`; `--zr-a`/`--zr-b`, hail
cap `--dbz-cap 53`, or `--no-hail-cap`), stacked, interpolated to 1-minute, and
reduced with the same trailing-window estimators as `mrms`. This is **v1** — one
fixed reflectivity Z–R for cross-era consistency (dual-pol fields only exist
post-~2012), lowest tilt only. On Hidden Valley it tracks MRMS spatially
(r ≈ 0.76–0.8) but runs ~1.5× MRMS / ~4× gauges — the expected high bias of a
simple single-radar Z–R before dual-pol corrections.
Because the field names match MRMS, compare the two directly with
`stormscape.compare`. (Level II is the right basis for a *consistent* multi-era
record; the Level III QPE archive mixes algorithms/grids across eras — legacy
single-pol DPA on a 4 km grid vs dual-pol DAA/DPR at 0.25 km.)

**v2 — dual-pol R(Kdp)** (`--method kdp`, 2012-on data). Where `Z ≥ --z-blend`
(default 35 dBZ) the rate comes from **specific differential phase R(Kdp)** —
hail-robust, since Kdp tracks liquid not ice — computed with Py-ART's variational
`kdp_maesaka`, blending back to capped Z–R in light rain and for pre-dual-pol
volumes (so it stays cross-era):

```bash
python -m stormscape nexrad --aoi aoi.kmz --intensity --method kdp \
    --start 202606192000 --end 202606200200 --out-dir ./out --key event_l2kdp
```

On Hidden Valley this pulls the single-radar i15 from **1.56× → 0.80× MRMS** and
from ~4× to **2.0× the gauges** — closer to the gauges than MRMS itself (2.3×).
R(Kdp) is less biased but noisier than Z–R (Kdp-derivative texture); raise
`--z-blend` to apply it only in the heaviest cells. A uniform **`--rate-cap`**
(e.g. `120`, mm/h) clips per-scan rate spikes for *any* method before stacking —
an operational hail-cap analogue that also lifts R(Kdp)'s gauge correlation
(0.68→0.74 on Hidden Valley); recommended. (Note: NEXRAD `differential
phase` is raw — the R(A) specific-attenuation path needs PhiDP preprocessing and
isn't used; `kdp_maesaka` self-regularizes.)

**Beam-blockage masking** (`--blockage-dem dem.tif`). Pass a DEM and the stack
flags cells whose cumulative beam blockage exceeds `--cbb-max` (default 0.5),
masking them and emitting a `cbb` quality field (0–1) — the single-radar analogue
of MRMS's RQI, for terrain-shadowed AOIs. Needs **wradlib** (`conda install -c
conda-forge wradlib`). On Hidden Valley it's ~0 (KRGX is a mountaintop radar
seeing the valley unobstructed). Inspect blockage directly with
`stormscape.beam_blockage(radar, aoi, dem)`.

In Python: `nearest_radar(aoi)`, `reflectivity_field(aoi, when)` /
`reflectivity_composite(aoi, start, end)` → a result dict (feed it to
`save_fields` / `drape_i15`), and for the radar-vs-gauge diagnostic
`sample_radar_at_points(radar, gauges)` + `z_to_rate(dbz)` — the WSR-88D
convective Z–R, with an optional hail cap to expose hail over-estimation.

---

### Diagnostic panels & virtual gauges

Two overview tools ported from the original `MRMS_stack`.

**Multi-panel diagnostic map** — tile the stacked companion fields (time of peak
i15, QPE storm total, RQI, beam-height SHSR) to judge *where and when* the peak
fell and whether the radar could see it cleanly. It shares the **same context as
the main maps** — add `--hillshade`, `--reference`, `--perimeters`, `--clip`:

```bash
python -m stormscape panels --radar-dir ./out --key event \
    --hillshade out/event_hillshade.tif --reference --perimeters aoi.kmz --clip \
    --gauges out/event_gauges.geojson --out out/event_panels.png
```

`--fields` picks which `<key>_<field>.tif` to panel (default `tpki15 total rqi
shsr`; any saved field works — `i15max i30max i60max i2max cbb`). Colourblind-safe
colormaps; time-of-peak is unwrapped across midnight.

**Virtual gauges** — drop point(s) into the radar grid and pull a rainfall time
series there: intensities over 5/15/30/60-min windows and cumulative total,
radar-only PrecipRate vs the gauge-corrected MultiSensor QPE:

```bash
python -m stormscape vgauge --date 20260619 \
    --point -119.709,39.485,HiddenValley --point -119.557,39.309,SixMile \
    --out-dir ./out --key event          # or --points-file points.geojson
```

Add `--gauges` to also drop a virtual gauge at **every real Synoptic station** in
the AOI and overlay the real gauge series for comparison. This **reuses a saved
gauge store** when you pass `--from-dir`/`--from-key` (the `gauges`-built
`<key>_gauges.geojson` + `RainGaugeData/`) — no re-fetch, no token — auto-trimming
the full storm-day record to the storm's rain window; otherwise it fetches live
(needs `$SYNOPTIC_TOKEN` or `--token`). `--atlas` then renders an **atlas** subplot
of all gauges, written to
`<key>_vg_atlas.png`. Per-gauge time-series **CSVs** go to `<out-dir>/RainGaugeData/`
— one file per gauge per source (`*_vgauge_mrms_*.csv`, `*_vgauge_nexrad_*.csv`,
`*_gauge_*.csv`) — and explicit `--point`s also get the 2-panel figure (top:
I5/I15/I30/I60 + I60 from QPE; bottom: cumulative PrecipRate vs QPE).

Add `--detail` to also write a **big 4-row figure per gauge** — cumulative rainfall,
then I60 / I30 / I15 — into `<out-dir>/VirtualGaugeFigures/`
(`<key>_vgdetail_<name>.png`), reusing the atlas line styles (MRMS blue, NEXRAD red,
real gauge dashed black). It's the full-size, per-gauge counterpart of one atlas
panel for analysing each gauge in detail.

Sources combine: `--source nexrad` uses the **single-radar NEXRAD Level II** stack
instead of MRMS (pre-2020 / MRMS-gap fallback; `--method za|kdp`), while `--nexrad`
*adds* the NEXRAD series **alongside** MRMS for ≥2020 events — so the atlas overlays
**MRMS VG, NEXRAD VG and the real gauge** on each panel and the CSVs include both
radar sources. In Python: `mrms.virtual_gauge_timeseries` /
`nexrad.virtual_gauge_timeseries` (both → `{name: DataFrame}`),
`gauges.gauge_timeseries` (real gauges, same shape), `plot.virtual_gauge_atlas`,
which takes `{source_label: {name: df}}` plus a `real_series` overlay, and
`plot.virtual_gauge_detail(sources, name, real_series=)` for one gauge's 4-row figure.

### Zoom into a sub-region of a processed event

Once an event is processed, **don't re-run to zoom in** — MRMS has no finer
resolution than its native ~1 km grid, so a fresh run just re-downloads identical
radar data. Instead, `zoom` re-renders the figures clipped to a sub-AOI straight
from the existing rasters:

```bash
python -m stormscape zoom \
    --from-dir ./out --from-key event \      # the already-processed event
    --bbox -119.65 39.34 -119.48 39.52 \     # the zoom window (or --aoi sub.geojson)
    --out-dir ./out/zoom_south --key south \
    --reference --gauges --crop-rasters       # denser local labels at the zoom scale
```

This writes `south.png` (zoomed map), `south_panels.png`, **and the NOAA Atlas 14
climatology set for the sub-AOI** (`south_climate_compare.png`, `south_clim_i{d}.tif`,
`south_anom_i{d}.png`/`.tif`) — reusing the source DEM/MRMS/gauge data with **no
re-download** (the climatology is fetched fresh for the tighter extent; the observed
field is reused). **Every rainfall map in the zoom folder is Gaussian-1km-smoothed
for display by default** — the i15 map, the rainfall panels (the categorical
`tpki15`/`rqi`/`shsr` panels stay raw), and the climatology — all via one
`--obs-smooth`/`--obs-smooth-radius` knob (`--obs-smooth none` for the raw fields;
display only — `--crop-rasters` GeoTIFFs stay raw). Reference vectors
(streams/roads/place names) are re-fetched at the tighter extent so you gain local
detail. `--crop-rasters` also writes cropped GeoTIFFs for a self-contained zoom
folder. The climate maps are on by default — pass `--no-climate` to skip them (a
failed Atlas 14 fetch is non-fatal and still leaves the map + panels). The **only**
product worth re-fetching is terrain — add `--refine-dem --resolution 3` to pull a
finer 3DEP DEM+hillshade for the zoom extent (3DEP has 3 m / 1 m tiers; MRMS does not).

Don't want to eyeball the bbox numbers? `pick` opens an **interactive browser
picker** — drag a rectangle on the event's map and it hands you the `--bbox` and a
ready-to-run `zoom` command:

```bash
python -m stormscape pick --from-dir ./out --from-key event
```

It writes a self-contained `event_pick.html` — the event's i15 map with the **full
context of the production maps** (labelled streams/roads/place names, AOI perimeter,
gauges, north arrow, lat/long ticks) plus a tiny JS canvas — and opens it in your
default browser. No server, no GUI toolkit, no internet, so it works the same on
Windows, macOS, and Linux. (Reference labels are on by default; `--no-reference`
skips the network fetch.)

### Rainfall climatology — NOAA Atlas 14 (observed vs the 1-year storm)

Put a storm in context: how did its observed I15/I30/I60 compare to the
**climatological** intensity for the same durations? `climate` fetches NOAA
**Atlas 14** gridded precipitation-frequency data (default the **1-year**
recurrence interval — the group's reference) and pairs it with an
already-processed event's observed fields (MRMS or NEXRAD), no radar re-run:

```bash
python -m stormscape climate \
    --from-dir ./out --from-key event \      # the already-processed event
    --out-dir ./out --key event \
    --ari 1 --durations 15 30 60 \
    --reference --gauges
```

This writes the climatology rasters (`event_clim_i{15,30,60}.tif`), a 3×2
**comparison** figure (`event_climate_compare.png` — rows = durations, left =
Atlas 14 climatology, right = observed), and per-duration **anomaly** maps
(`event_anom_i{15,30,60}.png` + `.tif`) of *observed ÷ climatology* on a
diverging colormap centred at 1× with integer contours (so a "3×" contour marks
where the storm hit three times the 1-year intensity). The AOI defaults to the
observed footprint; the Atlas 14 region is picked automatically (override with
`--region`). Atlas 14's authoritative gridded ASCII grids are used directly
(`pfdf.data.noaa` is point-only), depth converted to intensity to match the
observed mm/h fields — no extra dependencies.

The observed radar field is **Gaussian-smoothed at a 1 km radius by default**
(both the comparison figure and the anomaly) so the peaky ~1 km field is visually
comparable to the smooth ~800 m climatology and single-pixel spikes don't
dominate the anomaly; the observed panel is labelled with the method/radius. Use
`--obs-smooth none` for the raw field, or `--obs-smooth`/`--obs-smooth-radius` to
tune it (same four methods as `smooth`).

---

### Smoothing a radar field — methods, comparison, gauge-skill

The observed radar fields are peaky (~1 km, sharp convective cores) next to the
smooth Atlas 14 climatology. `smooth` evaluates *how much* and *which way* to
smooth, reusing an already-processed event's rasters (MRMS or NEXRAD), no
re-run:

```bash
python -m stormscape smooth \
    --from-dir ./out --from-key event \
    --out-dir ./out --key event \
    --field i15max --methods gaussian uniform median idw --radii 0 1 2 4 \
    --gauge-analysis --gauges --clip          # skill sweep needs the gauges geojson
```

It writes a **comparison figure** (`event_smoothing_compare.png`) — a methods ×
radii grid (column 0 = raw) on one shared colour scale so you can see the
peak-flattening — and, with `--gauge-analysis`, a **radar–gauge skill** figure +
CSV (`event_smoothing_skill.{png,csv}`) of correlation / RMSE / bias-ratio vs
smoothing radius, per method and duration, with the optimum starred. Four
NaN-aware methods are available: `gaussian` (recommended default), `uniform`
(boxcar mean), `median` (edge-preserving), and `idw` (inverse-distance, the
moving-window analogue of point IDW). `radius_km` is the nominal smoothing scale
(≈ Gaussian σ), mapped to a comparable extent per method. `--write <method>
--write-radius <km>` emits smoothed `<key>_<field>.tif`s that flow straight into
`compare`/`map`/`climate`. Uses SciPy only — no new dependencies.

Because **I15 is a peak metric**, smoothing mechanically lowers the radar's
positive bias, so trust the **correlation** (up) and **RMSE** (down) as the skill
signal — not the bias *ratio* (which falls toward 1× as a side effect). On Hidden
Valley the correlation peaks at a small radius (~0.5–1 km) then declines, and the
bias ratio stays ≫ 1× at every radius — i.e. smoothing can't fix that event's
hail-driven QPE over-read (consistent with the single-radar Z analysis).

---

### Gauge recurrence — anomaly + return period vs NOAA Atlas 14

`recurrence` builds a per-gauge table for every **wet** gauge (peak I15 > 0): the
observed peak **I15/I30/I60**, the **time of the I15 peak**, the **anomaly**
(observed ÷ 1-yr Atlas 14), and the **recurrence interval** of each peak.

```bash
python -m stormscape recurrence \
    --from-dir ./out --from-key event --out-dir ./out
```

The climatology comes from the NOAA **PFDS point** service per gauge
(`atlas14.pf_point`), which returns the full duration × ARI curve — so the 1-yr
anomaly reference and the recurrence interval share one point-accurate source (the
same Atlas 14 data the maps are tiled from). The recurrence interval is obtained
by log-log interpolating the observed value against that curve (Atlas 14 publishes
the quantile curve, not a closed-form inverse); `<1` means below the 1-yr quantile,
`>1000` above the top tabulated ARI. The I15 time-of-peak is read offline from the
event's saved `RainGaugeData/` per-gauge series — **no Synoptic token needed**.
Writes `<key>_gauge_recurrence.csv` + `.md`. (This is a point query, not raster
sampling: the downloaded 1-yr grids can't yield a recurrence interval, which needs
the whole ARI curve.)

---

### Georeferenced export for GIS / CalTopo

`export` re-packages a processed event for a GIS or **CalTopo** — no radar re-run.
It writes the rainfall fields as **EPSG:3857 GeoTIFFs** (Web-Mercator, CalTopo's
native projection) and the two primary maps as **georeferenced PDFs** (GeoPDFs).

```bash
python -m stormscape export \
    --from-dir ./out --from-key event --out-dir ./out --key event \
    --aoi event_AOI.kmz --gauges --reference
```

For each `--layers` field (default `anom_i15 i15max`) it writes both a **raw
single-band float** GeoTIFF (`event_i15max_3857.tif`, the data values) and a
**colorized RGBA** GeoTIFF (`event_i15max_3857_rgb.tif`, the project colormap with
transparent dry / no-data cells, so it drops onto a CalTopo basemap looking like
the figure). It also writes GeoPDFs of the i15 map (`event.pdf`) and the anomaly
map (`event_anom_i15.pdf`) — the full styled figure with the **map frame**
registered to its coordinates (a neatline excludes the colorbar / margins). The
PDFs render in UTM by default (identical to the PNG deliverables) and frame to the
event AOI like the `map`/`run` figures.

The GeoPDF half needs GDAL's **PDF driver** (`conda install -n GISMan -c
conda-forge libgdal-pdf`); without it the GeoTIFF export still runs and the PDFs
are skipped with a note. The anomaly GeoPDF reuses `event_anom_i15.tif` from a
prior `climate` run.

Add `--streams` to also export the **full-resolution NHD stream network** for the
AOI as a vector layer (`event_streams.geojson` — the same dense NHDPlus HR
flowlines, named creeks and unnamed headwaters, clipped to the AOI polygon, EPSG:4326),
ready to import into CalTopo or a GIS. `--streams-format {geojson,gpkg,shp,kml}`
picks the format; `--streams-bbox` keeps whole flowlines over the bounding box;
`--streams-named-only` drops the unnamed headwaters.

---

## Python API

```python
from stormscape import (fetch_dem_and_hillshade, i15_storm_day,
                        save_fields, drape_i15)

aoi = (-105.55, 40.55, -105.25, 40.80)            # or "perimeter.geojson"

dem, hs = fetch_dem_and_hillshade(aoi, resolution=10,
                                  dem_path="dem.tif", hillshade_path="hs.tif")

res = i15_storm_day(aoi, "2021-07-20")            # dict of fields + metadata
save_fields(res, ".", key="event")                # writes event_i15max.tif, ...
print(res["meta"])                                # peak time, AOI-max i15, RQI

drape_i15("hs.tif", "event_i15max.tif", out_path="event.png",
          perimeters="perimeter.geojson")
```

Gauges and the radar-vs-gauge comparison (needs `$SYNOPTIC_TOKEN`):

```python
import datetime as dt
from stormscape import gauge_fields, compare_storm

gauges = gauge_fields(aoi, dt.datetime(2021, 7, 20, 4),
                      dt.datetime(2021, 7, 21, 10))   # one row per gauge
table, stats = compare_storm(gauges, ".", key="event",
                             rqi_min=0.8, max_report_min=15)
print(stats)                       # bias / rmse / mae / corr / ratio per metric
```

Check 1 m lidar availability before requesting it:

```python
from stormscape import coverage_fraction
coverage_fraction(aoi, res="1m")    # fraction of the AOI with 1 m 3DEP source
```

---

## The i15 estimator

MRMS `PrecipRate` is a 2-minute instantaneous rate (mm/h). For each 2-min
step, `a2 = rate · 2/60` is the accumulation (mm). Over a trailing 16-min
window (8 steps), `i16 = Σ(8)·60/16` and `i14 = Σ(last 7)·60/14`, and
`i15 = mean(i16, i14)`. The per-cell running maximum over the storm is
`i15max`. (Estimator after D. Cavagna's `MRMS_stack.py`.) The 30- and 60-minute
peaks (`i30max`, `i60max`) use plain trailing windows (15 / 30 two-minute
steps), scaled to mm/h and kept as running maxima alongside `i15max`. The gauge
side computes the matching `i15` with the same `(i16+i14)/2` estimator on the
1-minute-interpolated gauge series, so the two are directly comparable.

The storm window is found automatically: hourly `RadarOnly_QPE_01H` is scanned
over the UTC window covering the local calendar day, the wettest hours are
kept, and 2-min `PrecipRate` is stacked over each contiguous wet run.

---

## Caveats

- **Coverage in time.** MRMS `PrecipRate` is available from ~**2014**; the
  2-min cadence and these product paths are reliable from ~**Oct 2020** on the
  NOAA S3 archive. Earlier dates may have gaps (missing steps reset the i15
  stack rather than fabricate continuity).
- **Radar quality.** Always inspect the `rqi` field. Far from radar / behind
  terrain, the beam overshoots low rain (`shsr` = beam height, km) and i15 is
  unreliable. Filter by RQI for quantitative work.
- **Local-day → UTC.** The scan window (`scan_pad_h=(4,10)`) assumes a CONUS
  local day. For other longitudes, widen it.
- **1 m DEMs** exist only where lidar has been flown; check `coverage_fraction`
  first and avoid mixing 1 m with resampled fill (seam artifacts in
  slope/curvature). Large AOIs at 1 m are very large downloads.
- **CRS.** DEMs are stored in EPSG:5070 (CONUS Albers, equal-area metres); MRMS
  GeoTIFFs in native EPSG:4326. **Figures are drawn in an auto-selected UTM
  zone** (near north-up, with latitude/longitude axis ticks); the hillshade and
  i15 are reprojected to it on the fly.
- **Radar vs gauge.** Radar QPE is a ~1 km² grid value; a gauge is a point, so
  exact agreement is not expected (point-vs-pixel). For the sub-hourly peaks,
  only gauges reporting at sub-15-min cadence resolve a true i15 — screen with
  `--max-report-min`. Radar-**only** MRMS (`PrecipRate` / `RadarOnly_QPE`) tends
  to **overestimate** intense convective rain; the gauge-corrected MultiSensor
  total (`--multisensor`) usually agrees better.
- **Colour.** The i15 field uses a colourblind-safe sequential colormap
  (`YlGnBu`) by default; the residual map's gauges use a diverging `RdBu`.
  Change the field map with `--cmap` (e.g. `inferno`, `cmc.lajolla`,
  `cmc.oslo`); avoid `jet`/`turbo` (perceptually non-uniform).

## Data sources

- USGS 3DEP elevation via [`py3dep`](https://github.com/hyriver/py3dep)
  (HyRiver) — the layer behind The National Map.
- NOAA MRMS on the public S3 archive
  (`noaa-mrms-pds.s3.amazonaws.com/CONUS/...`) — `PrecipRate`, `RadarOnly_QPE`,
  and gauge-corrected `MultiSensor_QPE_01H` (Pass-2).
- Single-radar NEXRAD **Level II** volumes on AWS (`unidata-nexrad-level2`) via
  [`nexradaws`](https://github.com/aarande/nexradaws), read with
  [Py-ART](https://arm-doe.github.io/pyart/). WSR-88D site coordinates from NCEI
  HOMR (`data/nexrad_sites.csv`).
- Ground rain gauges via the **Synoptic / MesoWest** Time Series API
  (`precip=1` Basic Precipitation Service); needs a free `$SYNOPTIC_TOKEN`. The
  gauge transport + 1-minute interpolation adapt the USGS **FlowAlert** package
  (King, Rengers, Wedell & Fee, 2024; CC0).
- **NOAA Atlas 14** precipitation-frequency climatology via the gridded ASCII
  grids on the PFDS GIS server (`hdsc.nws.noaa.gov/pub/hdsc/data/...`); region
  extents bundled in `data/atlas14_regions.csv`.

## Layout

```
stormscape/
├── stormscape/
│   ├── aoi.py       AOI parsing (bbox / vector / geometry) + overlay loading
│   ├── dem.py       3DEP DEM download, 1 m availability, hillshade
│   ├── mrms.py      MRMS fetch/stack -> i15/i30/i60 fields + MultiSensor (engine)
│   ├── nexrad.py    single-radar NEXRAD Level II tilts + intensity stacks
│   ├── gauges.py    Synoptic/MesoWest gauges -> total + peak 15/30/60 intensities
│   ├── compare.py   sample radar at gauges -> residuals, skill stats, recurrence
│   ├── merge.py     radar-gauge bias correction + conditional merge
│   ├── atlas14.py   NOAA Atlas 14 climatology grids -> intensity fields + anomaly
│   ├── smoothing.py NaN-aware field smoothing + radar-gauge skill sweep
│   ├── export.py    EPSG:3857 GeoTIFFs, GeoPDF figures, NHD stream vectors
│   ├── refdata.py   AOI-scoped NHD streams / TIGER roads / GNIS places
│   ├── plot.py      drape fields over hillshade + basemap/vector/gauge overlays
│   ├── data/        bundled tables (nexrad_sites.csv, atlas14_regions.csv)
│   └── cli.py       15 subcommands (see docs/cli.md)
├── examples/            worked example + batch event templates
├── docs/cli.md          full CLI reference
├── environment.yml      conda-forge environment ("stormscape")
└── pyproject.toml
```

---

## License & credits

Released under the **[MIT License](LICENSE)** — use, modify, and redistribute
freely with attribution.

Third-party code and the terms of every data source are recorded in
**[NOTICE.md](NOTICE.md)**. In brief: the gauge transport is adapted from USGS
**FlowAlert** (CC0), the i15 estimator follows D. Cavagna's MRMS stacking
approach, and all data is public-domain US government data — *except* the
Synoptic / MesoWest gauge API, which needs your own free token.

If this contributes to published work, please cite the repository and the
underlying data products (see [NOTICE.md](NOTICE.md#citing-this-software)).
