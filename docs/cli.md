# `stormscape` — command-line reference

Invoke as `python -m stormscape <command>` (or `stormscape <command>` once the
package is installed). This page mirrors `python -m stormscape <command> --help`;
run that for the authoritative, always-current text.

```
stormscape {dem, i15, map, run, gauges, compare, nexrad, panels, vgauge, zoom,
            pick, climate, smooth, recurrence, export, burn}
```

| Command | Purpose |
|---|---|
| `dem` | Download a 3DEP DEM + hillshade for an AOI |
| `i15` | Stack MRMS → peak I15/I30/I60 (+ companion) fields for an AOI + storm-day |
| `map` | Drape an existing i15 GeoTIFF over an existing hillshade GeoTIFF |
| `run` | `dem` → `i15` → figure in one shot (+ optional gauges/compare) |
| `gauges` | Fetch Synoptic gauge storm-total + peak 15/30/60-min intensities |
| `compare` | Sample radar rasters at gauges → residuals + skill stats (+ residual map) |
| `nexrad` | Single-radar NEXRAD Level II: reflectivity field **or** an i15 intensity stack |
| `panels` | Multi-panel diagnostic map (time-of-peak i15, QPE total, RQI, SHSR) |
| `vgauge` | Virtual rain gauges: radar rainfall time series at point(s) / all stations |
| `zoom` | Re-render an existing event's figures clipped to a sub-AOI (reuses rasters; no MRMS re-download) |
| `pick` | Interactive browser bbox picker — drag a zoom box on an event's map (cross-platform, no GUI toolkit) |
| `climate` | NOAA Atlas 14 rainfall climatology vs observed I15/I30/I60 — comparison figure + anomaly maps (reuses rasters; no re-run) |
| `smooth` | Smooth a radar field: methods×radii comparison figure + optional radar–gauge skill sweep (reuses rasters; no re-run) |
| `recurrence` | Wet-gauge table: peak I15/I30/I60 + time-of-peak, anomaly + recurrence interval vs NOAA Atlas 14 (PFDS point query) |
| `export` | Georeferenced EPSG:3857 GeoTIFFs + GeoPDF figures + NHD stream vectors for GIS / CalTopo (reuses rasters; no re-run) |
| `burn` | Near-real-time burn severity over an AOI from CIMSS BRISK daily dNBR (or the BAER soil-burn-severity archive) |

The AOI is **always** given as either `--bbox W S E N` (lon/lat degrees) or
`--aoi <vector file>` (GeoJSON/SHP/GPKG/KMZ). Gauge steps need a free Synoptic
token via `$SYNOPTIC_TOKEN` or `--token`.

---

## Shared option groups

These repeat across commands; described once here.

### A. AOI & I/O — `dem, i15, run, gauges, compare, nexrad, vgauge, zoom, climate, smooth, recurrence, export` (and `--out-dir/--key` on `map`)
| Arg | Type / default | Meaning |
|---|---|---|
| `--bbox W S E N` | 4 floats | AOI bounding box, lon/lat degrees |
| `--aoi AOI` | path | AOI vector file (alternative to `--bbox`) |
| `--pad-deg` | float, `0.05` | degrees to pad the AOI bounds |
| `--dst-crs` | str, `EPSG:5070` | working CRS for DEM/figure |
| `--out-dir` | path | output directory |
| `--key` | str | filename stem for outputs |

### B. Gauge / time window — `gauges, compare, run, nexrad, vgauge`
| Arg | Type / default | Meaning |
|---|---|---|
| `--date` | `YYYYMMDD` | storm day; sets the default UTC window |
| `--start` | `YYYYMMDDHHMM` | explicit UTC window start (with `--end`, overrides `--date`) |
| `--end` | `YYYYMMDDHHMM` | explicit UTC window end |
| `--token` | str | Synoptic API token (else `$SYNOPTIC_TOKEN`) |
| `--durations MIN…` | ints, `15 30 60`† | peak-intensity window durations, minutes |

† `vgauge` defaults to `5 15 30 60`.

### C. Figure styling — `map, run, nexrad, zoom, climate, smooth, panels, compare, export`
| Arg | Type / default | Meaning |
|---|---|---|
| `--perimeters` | path | vector overlay outlined as a border |
| `--basins` | path | vector overlay as thin outlines |
| `--highlight` | path | vector overlay drawn bold (cyan) |
| `--points` | path | point vector overlay (triangles) |
| `--title` | str | figure title |
| `--cmap` | str, `YlGnBu` | i15 field colormap (colorblind-safe; e.g. `inferno`, `cmc.lajolla`) |
| `--wet-min` | float, `5.0` | mask i15 below this mm/h in the drape |
| `--basemap` | flag | underlay USGS topo tiles via contextily |
| `--basemap-provider` | str, `USGS.USTopo` | contextily provider key |
| `--basemap-labels` | str | labels-only provider on top (`none` to disable) |
| `--basemap-zoom` | int, auto | basemap tile zoom level |
| `--hillshade-alpha` | float, `1.0` (`0.0` w/ basemap) | hillshade opacity |
| `--alpha` | float, `0.32` | field drape opacity over terrain/basemap (project-wide default; same on `map`/`run`/`nexrad`/`zoom`/`climate`/`smooth`/`panels`/`compare`) |
| `--reference` | flag | auto-fetch + overlay NHD streams / TIGER roads / GNIS places |
| `--local-roads` | flag | include residential/local roads in `--reference` |
| `--no-reference-labels` | flag | reference lines/points without text labels |
| `--clip` | flag | clip figure tightly to AOI / `--perimeters` |
| `--clip-margin` | float, `0.04` | fractional margin around the clip extent |
| `--dpi` | int, `200` | output figure resolution |

---

## Command-specific arguments

### `dem` — DEM + hillshade
Group **A**, plus:
| Arg | Type / default | Meaning |
|---|---|---|
| `--resolution` | int, `10` | DEM resolution, metres |
| `--clip-dem` | flag | mask DEM outside the AOI polygon |

### `i15` — MRMS i15/i30/i60 stack  *(`--date` required)*
Group **A**, plus:
| Arg | Type / default | Meaning |
|---|---|---|
| `--date` | `YYYYMMDD` | **required** storm day |
| `--qpe-thresh` | float, `2.5` | hourly-QPE threshold (mm) for "wet hour" selection |
| `--max-wet-hours` | int, `8` | cap on wet hours stacked |
| `--workers` | int, `12` | parallel download workers |
| `--multisensor` | flag | also fetch gauge-corrected MultiSensor QPE total → `<key>_mstotal.tif` |

### `map` — drape i15 over hillshade  *(`--hillshade` + `--i15` required)*
Group **C**, plus `--out-dir`, `--key`, `--dst-crs`, and:
| Arg | Type / default | Meaning |
|---|---|---|
| `--hillshade` | path | **required** hillshade GeoTIFF |
| `--i15` | path | **required** i15 GeoTIFF |
| `--out` | path | output PNG path |

### `run` — full pipeline  *(`--date` required)*
Groups **A + B + C**, plus the `i15` options (`--resolution`, `--clip-dem`,
`--qpe-thresh`, `--max-wet-hours`, `--workers`, `--multisensor`) and:
| Arg | Type / default | Meaning |
|---|---|---|
| `--gauges` | flag | fetch Synoptic gauges + overlay on the figure |
| `--compare` | flag | also compute radar-vs-gauge stats + residual map (implies `--gauges`) |
| `--rqi-min` | float | ignore gauges whose radar RQI is below this (compare) |
| `--max-report-min` | float | for i15/i30/i60: ignore gauges coarser than this native reporting interval (min) |

### `gauges` — full gauge pipeline (store + atlas + detail figures)
Groups **A + B**, plus the pipeline options below. One Synoptic draw does the whole
gauge workflow:

1. **Canonical store** (reused by `compare` / `recurrence` / `vgauge`, no re-fetch):
   - `<key>_gauges.geojson` — one row per gauge: coords, `total_mm`,
     `i15/i30/i60_mmph`, `report_min` (native precip cadence), and the **I15
     time-of-peak** (`i15_peak_time`).
   - `RainGaugeData/<key>_gauge_<name>.csv` — the self-describing 1-minute series
     (`rate_mmph`/`total_mm`/`i{d}_mmph` + `lon`/`lat`/`station_id`).
2. **Virtual gauges + rainfall comparison atlas** (`<key>_vg_atlas.png`) — a VG at
   every gauge from MRMS **and NEXRAD** (3-way; blue/red/dashed-black) by default.
3. **Per-gauge detail figures** → `VirtualGaugeFigures/<key>_vgdetail_<name>.png`.
   Each title annotates the ground gauge's native **reporting tempo** — e.g.
   `(ground gauge: ~5-min reporting)` / `~hourly (60-min)` / `~daily (1440-min)` —
   from the `report_min` field (median interval between precip-bearing obs). Coarse
   reporters have their bursts smeared by the 1-min interpolation, so their i15/i30
   read low; the tempo flags how far to trust the short-duration peaks.

The store keeps the **whole storm day**; steps 2–3 clip to the storm's rain window
(`gauges.storm_window`) and run on the **wet, near-AOI** stations — the SAME station
set feeds the atlas and the detail figures, so their gauges are uniform. Use a
generous `--pad-deg` to include regional gauges (the store keeps them all; distance
flag via `recurrence`).

| arg | type | meaning |
|---|---|---|
| `--store-only` | flag | write only the store; skip the atlas + detail (the old behaviour) |
| `--no-atlas` / `--no-detail` | flag | skip just the atlas / just the detail figures |
| `--no-nexrad` | flag | MRMS + real only (default adds single-radar NEXRAD → 3-way) |
| `--all-gauges` | flag | include dry stations (default: wet only) |
| `--wet-min` | float | peak intensity (mm/h) below which a station is "dry" (default 0.5) |
| `--max-dist-km` | float | also drop stations farther than this from the AOI |
| `--max-report-min` | float | keep only gauges reporting at this cadence or finer, min (e.g. `60` = hourly or finer; drops daily/coarse reporters whose i15/i30 are smeared by the 1-min interpolation) |
| `--atlas-metric` | int | intensity duration (min) shown in the atlas panels (default 15) |
| `--no-multisensor` | flag | skip the MRMS MultiSensor QPE overlay |
| `--radar` / `--method` | str | NEXRAD radar id (nearest if unset) / rate recipe (`kdp` default) |
| `--cache-dir` | path | NEXRAD volume cache (default `<out-dir>/nexrad_cache`) |
| `--pad-min` | int | padding (min) around the auto-detected rain window (default 30) |
| `--no-series` | flag | store the geojson only, skip the per-gauge series CSVs |

Needs `$SYNOPTIC_TOKEN` (or `--token`). To regenerate the atlas/detail from a saved
store **without** re-fetching gauges, use `vgauge --from-dir/--from-key`.

### `compare` — radar vs gauges → CSV (+ optional map)
Groups **A + B**, plus:
| Arg | Type / default | Meaning |
|---|---|---|
| `--gauges` | path | precomputed gauges GeoJSON (skips live fetch) |
| `--radar-dir` | path, `--out-dir` | dir holding `<key>_*.tif` radar rasters |
| `--rqi-min` | float | drop gauges below this radar RQI |
| `--max-report-min` | float | i15/i30/i60 only: drop gauges coarser than this cadence (min) |
| `--multisensor` | flag | add a gauge-corrected MultiSensor total row (needs `<key>_mstotal.tif`) |
| `--out` | path | output CSV for the per-gauge table |
| `--map` | path | also render a residual-map PNG here |
| `--map-metric` | str, `i15max` | radar field used for `--map` residuals |
| `--hillshade` | path | hillshade background for `--map` |
| `--perimeters` | path | AOI/perimeter overlay for `--map` |
| `--reference` | flag | labelled vector overlay on `--map` |
| `--clip` | flag | clip `--map` to `--perimeters` |
| `--title` | str | title for `--map` |
| `--cmap` | str, `YlGnBu` | i15 background colormap |

### `nexrad` — single-radar Level II
Groups **A + C**, plus `--date/--token/--durations` from **B**, plus:
| Arg | Type / default | Meaning |
|---|---|---|
| `--time` | `HHMM` | UTC scan time for a single scan (default: storm-window midpoint) |
| `--start` / `--end` | `YYYYMMDDHHMM` | window for `--composite` / `--intensity` |
| `--radar` | 4-letter id | radar (default: nearest to AOI) |
| `--field` | str, `reflectivity` | pyart field (`velocity`, `cross_correlation_ratio`, …) |
| `--sweep` | int, lowest | elevation-tilt index |
| `--composite` | flag | per-cell max over all scans (storm-peak reflectivity) |
| `--intensity` | flag | build the i15/i30/i60 rainfall stack instead of a reflectivity field |
| `--zr-a` / `--zr-b` | float, `300` / `1.4` | Z–R coefficients in `Z=a·R^b` |
| `--dbz-cap` | float, `53` | hail cap (dBZ) before Z–R |
| `--no-hail-cap` | flag | disable the dBZ hail cap |
| `--method {za,kdp}` | `za` | intensity rate: `za`=capped convective Z–R (all eras); `kdp`=R(Kdp) blended (2012+) |
| `--z-blend` | float, `35` | `kdp`: use R(Kdp) at/above this dBZ, Z–R below |
| `--rate-cap` | float | clip every method's per-scan rate to this max mm/h (operational hail-cap analogue) |
| `--blockage-dem` | path | DEM for beam-blockage masking (adds a `cbb` field) |
| `--cbb-max` | float, `0.5` | cumulative beam-blockage fraction above which to mask |
| `--res-m` | int, `500` | output grid resolution, metres |
| `--max-scans` | int, `40` | cap on scans gridded for `--composite` |
| `--cache-dir` | path, `<out-dir>/nexrad_cache` | downloaded-volume cache |
| `--hillshade` | path | hillshade to drape over |
| `--gauges` | flag | fetch + overlay Synoptic gauges coloured by sampled radar value |
| `--gauges-file` | path | precomputed gauges GeoJSON to overlay |

### `panels` — multi-panel diagnostic map  *(`--key` required)*
| Arg | Type / default | Meaning |
|---|---|---|
| `--radar-dir` | path, `.` | dir with `<key>_<field>.tif` |
| `--key` | str | **required** field stem |
| `--fields FIELDS…` | strs, `tpki15 total rqi shsr` | which fields to panel |
| `--perimeters` | path | AOI/perimeter outline |
| `--gauges` | path | gauges GeoJSON overlaid as points |
| `--hillshade` | path | hillshade to drape the panels over |
| `--reference` / `--local-roads` / `--no-reference-labels` | flags | labelled NHD/TIGER/GNIS overlay controls |
| `--clip` / `--clip-margin` | flag / `0.04` | clip panels to `--perimeters` |
| `--out` | path, `<radar-dir>/<key>_panels.png` | output PNG |
| `--title` / `--dpi` | str / `200` | title, resolution |

### `vgauge` — virtual rain gauges
Groups **A + B** (durations default `5 15 30 60`), plus:
| Arg | Type / default | Meaning |
|---|---|---|
| `--point LON,LAT[,NAME]` | repeatable | virtual-gauge location (e.g. `-119.7,39.5,Site`) |
| `--points-file` | path | vector file of gauge points |
| `--no-multisensor` | flag | skip the hourly MultiSensor QPE overlay |
| `--gauges` | flag | also drop a VG at every real Synoptic station + overlay the real series |
| `--from-dir` / `--from-key` | path / str | **reuse a saved gauge store** (`<from-key>_gauges.geojson` + `RainGaugeData/`) instead of re-fetching; the full storm-day record is auto-trimmed to the rain window (no `$SYNOPTIC_TOKEN` needed) |
| `--refetch` | flag | ignore the saved store and fetch live from Synoptic |
| `--pad-min` | int, `30` | padding (min) around the auto-detected rain window when reusing |
| `--max-dist-km` | float | with `--gauges`: drop stations farther than this from the AOI (`--bbox`/`--aoi`, else the `<from-key>_i15max` footprint) |
| `--wet-only` | flag | with `--gauges`: drop stations whose peak atlas-metric intensity is below `--wet-min` (declutters the atlas) |
| `--wet-min` | float, `0.5` | `--wet-only` floor (mm/h) below which a station counts as dry (drops traces) |
| `--max-report-min` | float | with `--gauges`: keep only gauges reporting at this cadence or finer, min (e.g. `60` = hourly or finer; drops coarse reporters) |
| `--atlas` | flag | render the all-gauge atlas → `<key>_vg_atlas.png` |
| `--atlas-metric` | int, `15` | duration (min) shown in atlas panels |
| `--detail` | flag | write a big 4-row figure per gauge (cumulative + I60/I30/I15) → `VirtualGaugeFigures/` |
| `--source {mrms,nexrad}` | `mrms` | primary rainfall source |
| `--nexrad` | flag | *also* include NEXRAD L2 alongside MRMS on the atlas/CSVs (≥2020 events) |
| `--radar` | 4-letter id | radar for NEXRAD source (nearest if unset) |
| `--method {za,kdp}` | `kdp` | NEXRAD rate recipe |
| `--cache-dir` | path | NEXRAD volume cache dir |
| `--dpi` | int, `200` | figure resolution |

### `zoom` — re-render an existing event at a sub-AOI
Re-draws the main i15 map, diagnostic panels, **and the NOAA Atlas 14 climatology
comparison + anomaly maps** — all **clipped to a sub-AOI from the already-processed
rasters** (no MRMS re-download; MRMS has no finer resolution). The climatology is
fetched fresh for the tighter extent and the observed field is reused from the
source event. **Every rainfall map in the zoom folder is Gaussian-1km-smoothed for
display by default** — the i15 map, the rainfall panels (intensities + storm
total; the categorical `tpki15`/`rqi`/`shsr` panels stay raw), and the climatology
observed/anomaly (like [`climate`](#climate--noaa-atlas-14-rainfall-climatology-vs-observed)) — all driven by the single
`--obs-smooth`/`--obs-smooth-radius` knob (`--obs-smooth none` for raw). Display
only; the cropped GeoTIFFs (`--crop-rasters`) stay raw. The `--bbox`/`--aoi` here
is the *zoom window*; `--out-dir`/`--key` are the *new* output. Groups **A + C**
(it always clips to the sub-AOI, so `--clip` is implied), plus the
[`climate`](#climate--noaa-atlas-14-rainfall-climatology-vs-observed) knobs
(`--ari`, `--durations`, `--region`, `--stat`, `--anomaly-cmap`, `--obs-smooth`,
`--obs-smooth-radius`, `--eps`, `--shared-row-scale`), and:
| Arg | Type / default | Meaning |
|---|---|---|
| `--from-dir` | path, **required** | directory of the processed event (its `--out-dir`) |
| `--from-key` | str, **required** | key/stem of the processed event |
| `--fields FIELDS…` | strs, `tpki15 total rqi shsr` | panel field stems to re-render |
| `--gauges` | flag | overlay the source `<from-key>_gauges.geojson` if present |
| `--refine-dem` | flag | re-fetch a finer DEM+hillshade for the zoom extent (the one product that benefits; retries the slow 3DEP 1 m WMS, then falls back to the existing hillshade so the zoom still renders) |
| `--resolution` | int, `10` | DEM resolution for `--refine-dem`, metres |
| `--clip-dem` | flag | with `--refine-dem`, mask the DEM outside the zoom polygon |
| `--crop-rasters` | flag | also write cropped copies of the source GeoTIFFs to the new folder (self-contained zoom) |
| `--no-map` / `--no-panels` / `--no-climate` | flags | skip the main map / the panels / the climatology + anomaly maps |

The climate step is **on by default** (like the map and panels); pass
`--no-climate` to skip it (e.g. a quick terrain-only zoom, or when offline). A
failed Atlas 14 fetch is non-fatal — it warns and still writes the map + panels.

For speed, the hillshade is reprojected + **downsampled to render resolution once**
(scaled to `--dpi`) and reused across all the figures, rather than re-processed
per figure — a fine (e.g. 1 m) zoom hillshade would otherwise dominate the runtime
(it's far more detail than a figure can show). This is invisible in the output and
the on-disk DEM/hillshade GeoTIFFs are left at full resolution.

### `pick` — interactive browser bbox picker
Writes a **self-contained HTML** file and opens it in your browser; drag a
rectangle to read off `--bbox W S E N` and a ready-to-run `zoom` command (copy
buttons). The background is the event's i15 map rendered with the **full main-map
context** (labelled NHD/TIGER/GNIS reference, AOI perimeter, gauges, north arrow,
lat/long ticks, colorbar — same look as the production maps); the JS maps clicks
to accurate lon/lat via the map axes' corner coordinates. No server, no GUI
toolkit, no internet — works on Windows/macOS/Linux. Not a group-A command (it
points at an existing event, not an AOI):
| Arg | Type / default | Meaning |
|---|---|---|
| `--from-dir` | path, **required** | directory of the processed event |
| `--from-key` | str, **required** | key/stem of the processed event |
| `--i15` | path | i15 GeoTIFF (default `<from-dir>/<from-key>_i15max.tif`) |
| `--hillshade` | path | hillshade GeoTIFF (default `<from-dir>/<from-key>_hillshade.tif`) |
| `--out` | path | output HTML (default `<from-dir>/<from-key>_pick.html`) |
| `--cmap` | str, `YlGnBu` | colormap for the i15 background |
| `--wet-min` | float, `5.0` | mask i15 below this mm/h in the background |
| `--no-reference` | flag | skip the NHD/TIGER/GNIS overlay (**on by default**; disabling skips a network fetch) |
| `--local-roads` | flag | include residential/local roads in the reference overlay |
| `--no-reference-labels` | flag | reference lines/points without text labels |
| `--perimeters` | path | AOI/perimeter outline to draw on the map |
| `--gauges` | flag | overlay the source `<from-key>_gauges.geojson` if present |
| `--no-open` | flag | write the HTML but don't open a browser |

### `climate` — NOAA Atlas 14 climatology vs observed

Puts a storm in climatological context. Reuses an already-processed event's
**observed** fields (`<from-key>_i{15,30,60}max.tif`, from MRMS *or* NEXRAD) — no
radar re-run — fetches the matching NOAA **Atlas 14** gridded precipitation-
frequency climatology (default the **1-year** ARI), and writes (1) the
climatology rasters, (2) a 3×2 **comparison** figure (rows = durations; left =
Atlas 14 climatology, right = observed), and (3) per-duration **anomaly** maps
(observed ÷ climatology) on a diverging colormap centred at 1× with integer
contours. The climatology is the authoritative gridded ASCII product (EPSG:4269,
~800 m), depth converted to intensity (mm/h) to match the observed fields; the
region is chosen automatically from the AOI (override with `--region`).

**Extent / matching the i15 maps — automatic.** With no `--aoi`/`--bbox`, `climate`
now **auto-matches the event's AOI** so the figures frame *identically* to the
`map`/`run` i15 maps (and draws the AOI outline). It looks in `--from-dir` for
`<from-key>_aoi.geojson` (written automatically by `dem`/`i15`/`run`/`nexrad`) then
a user-placed `<from-key>_AOI.{kmz,geojson,gpkg,shp}`. Only if **no** AOI is found
does it fall back to the wider **observed (`i15max`) raster footprint** (= the AOI
**plus the ~0.05° MRMS fetch pad**). Override anytime with `--aoi`/`--bbox`
(+ optional `--perimeters`); the same `--clip-margin` (0.04) is used on both paths.
Groups **A + C**, plus:
| Arg | Type / default | Meaning |
|---|---|---|
| `--from-dir` | path, **required** | directory of the processed event (its `--out-dir`) |
| `--from-key` | str, **required** | key/stem of the processed event |
| `--ari` | int, `1` | recurrence interval, years (1, 2, 5, …) |
| `--durations MIN…` | ints, `15 30 60` | intensity durations, minutes |
| `--region` | str | Atlas 14 region code override (else auto; e.g. `sw`, `tx`, `se`, `mw`, `ne`, `orb`, `inw`) |
| `--stat {mean,lower,upper}` | `mean` | grid statistic (mean estimate or 90% CI bound) |
| `--hillshade` | path | hillshade (default `<from-dir>/<from-key>_hillshade.tif`) |
| `--gauges` | flag | overlay the source `<from-key>_gauges.geojson` if present |
| `--anomaly-cmap` | str, `cmc.vik` | diverging colormap for the anomaly maps |
| `--shared-row-scale` | flag | share the colour scale within each comparison row (default: each panel independent) |
| `--obs-smooth {gaussian,uniform,median,idw,none}` | `gaussian` | smooth the observed radar field before the comparison **and** anomaly (so the peaky ~1 km field reads against the smooth ~800 m climatology); `none` keeps the raw field |
| `--obs-smooth-radius KM` | float, `1.0` | observed-field smoothing radius (nominal scale ≈ Gaussian σ) |
| `--eps` | float, `0.1` | climatology floor (mm/h) below which the anomaly is masked |
| `--no-comparison` / `--no-anomaly` | flags | skip the comparison figure / the anomaly maps |
| `--cmap` | str, `YlGnBu` | sequential colormap for the comparison panels |

The AOI defaults to the observed `i15max` footprint; pass `--bbox`/`--aoi` to
restrict it (e.g. to a zoom region). Climatology grids are cached under
`<out-dir>/atlas14_cache/`.

By default the observed radar field is **Gaussian-smoothed at a 1 km radius**
before both the comparison figure and the anomaly (and the written
`<key>_anom_i{d}.tif`), so the peaky radar field is comparable to the smooth
climatology and single-pixel spikes don't dominate the recurrence-multiple
anomaly. The observed column is labelled with the method/radius. Pass
`--obs-smooth none` for the raw field, or tune `--obs-smooth`/`--obs-smooth-radius`
(see [`smooth`](#smooth--spatial-smoothing-of-a-radar-field) for the methods).

---

### `smooth` — spatial smoothing of a radar field

Evaluates **how much** and **which way** to smooth a radar intensity field, by
reusing an already-processed event's `<from-key>_<field>.tif` (MRMS *or* NEXRAD) —
no radar re-run. Two products: (1) a **comparison figure** — a methods × radii
grid (column 0 = raw) on a single shared colour scale so the peak-flattening is
visible; (2) with `--gauge-analysis`, a **radar–gauge skill sweep** (a CSV +
figure of correlation / RMSE / bias-ratio vs smoothing radius, per method and
duration, with the optimum starred). Four NaN-aware methods: `gaussian`
(recommended default), `uniform` (boxcar mean), `median` (edge-preserving), `idw`
(inverse-distance, the moving-window analogue of point IDW). `radius_km` is the
nominal smoothing scale (≈ Gaussian σ), mapped to a comparable pixel extent per
method. Optionally `--write` smoothed field tifs that flow straight into
`compare`/`map`/`climate`. Groups **A + C**, plus:
| Arg | Type / default | Meaning |
|---|---|---|
| `--from-dir` | path, **required** | directory of the processed event (its `--out-dir`) |
| `--from-key` | str, **required** | key/stem of the processed event |
| `--field` | str, `i15max` | field to smooth/compare (e.g. `i30max`, `i60max`, `total`) |
| `--methods …` | choices, all four | methods for the comparison grid (`gaussian uniform median idw`) |
| `--radii KM…` | floats, `0 1 2 4` | radii (km) for the comparison columns (`0` = raw) |
| `--power` | float, `2` | IDW distance power |
| `--no-shared-scale` | flag | scale each comparison panel independently (default: one shared scale) |
| `--gauge-analysis` | flag | sweep radius vs radar–gauge skill → `<key>_smoothing_skill.{png,csv}` |
| `--sweep KM…` | floats, `0 0.5 1 2 3 4 6 8` | radii for the skill sweep |
| `--gauges` | flag | overlay the source `<from-key>_gauges.geojson` on the maps |
| `--gauges-file` | path | gauges GeoJSON for `--gauge-analysis` (overrides the source geojson) |
| `--rqi-min` | float | ignore gauges whose radar RQI is below this (skill sweep) |
| `--max-report-min` | float | for I15/I30/I60: ignore gauges whose reporting interval (min) exceeds this |
| `--write {gaussian,uniform,median,idw}` | — | write smoothed field tifs at this method (needs `--write-radius`) |
| `--write-radius KM` | float | radius for `--write` |
| `--no-comparison` | flag | skip the comparison figure |
| `--durations MIN…` | ints, `15 30 60` | durations for the skill sweep |
| `--hillshade` | path | hillshade (default `<from-dir>/<from-key>_hillshade.tif`) |

Because I15 is a **peak** metric, smoothing mechanically lowers the radar's
positive bias — so read the **correlation** (up) and **RMSE** (down) as the skill
signal, not the bias *ratio* (which falls toward 1× as a side effect; the ratio
panel marks 1× but is not the criterion).

---

### `recurrence` — wet-gauge anomaly + recurrence interval
Builds a per-gauge table for every **wet** gauge (peak I15 > 0) in an event's
`<from-key>_gauges.geojson`: the observed peak **I15/I30/I60**, the **time of the
I15 peak**, the **anomaly** (observed ÷ 1-yr Atlas 14), and the **recurrence
interval** of each peak. The climatology is the NOAA **PFDS point** service per
gauge (`atlas14.pf_point`) — the full duration × ARI curve, so the 1-yr anomaly
reference and the recurrence interval share one point-accurate source (same Atlas
14 data as the maps). The recurrence interval is log-log interpolation of the
observed value against that curve (`<1` below the 1-yr quantile, `>1000` above the
top ARI). Time-of-peak comes from the geojson's `i15_peak_time` (or, for older
stores, the saved `RainGaugeData/` series) — no Synoptic token needed. Each gauge
is flagged by distance to the AOI (`dist_to_aoi_km`/`in_aoi`) so near-AOI regional
gauges are included but identifiable (`--max-dist-km` to drop far ones). Writes
`<key>_gauge_recurrence.csv` + `.md`.
| Arg | Type / default | Meaning |
|---|---|---|
| `--from-dir` | path, **required** | directory of the processed event (its `--out-dir`) |
| `--from-key` | str, **required** | key/stem of the processed event |
| `--max-dist-km` | float | drop gauges farther than this from the AOI (default: keep all, flagged) |
| `--gauges-file` | path | gauges geojson (default `<from-dir>/<from-key>_gauges.geojson`) |
| `--durations MIN…` | ints, `15 30 60` | intensity durations |
| `--series {pds,ams}` | `pds` | PF series: partial-duration (default; matches the maps) or annual-maximum |
| `--stat {mean,lower,upper}` | `mean` | PF estimate (mean, or 90% CI bound) |
| `--raingauge-dir` | path | per-gauge series CSVs for time-of-peak (default `<from-dir>/RainGaugeData`) |
| `--no-peak-time` | flag | skip the I15 time-of-peak column |

---

### `export` — georeferenced GeoTIFFs + GeoPDFs for GIS / CalTopo

Re-exports an already-processed event into formats a GIS or **CalTopo** can use —
no radar re-run. Two products:

1. **EPSG:3857 GeoTIFFs** (Web-Mercator, CalTopo's native projection) of the
   rainfall fields. For each `--layers` field (default `anom_i15 i15max`) it writes
   both a **raw single-band float** GeoTIFF (`<key>_<field>_3857.tif`, the data
   values for analysis / restyling) and a **colorized RGBA** GeoTIFF
   (`<key>_<field>_3857_rgb.tif`, the project colormap with **transparent dry /
   no-data cells**, so it drops onto a CalTopo basemap looking like the figure).

2. **Georeferenced PDFs** (GeoPDFs) of the two primary map figures — the i15 map
   (`<key>.pdf`) and the anomaly map (`<key>_anom_i{d}.pdf`). The full styled
   figure (hillshade, draped field, reference labels, gauges, north arrow, ticks,
   colorbar) is embedded and the **map frame** is registered to its projected
   coordinates with a **neatline** bounding the georeferenced area (so the colorbar
   / margins are excluded, not mis-placed). Readable as a georeferenced layer by
   CalTopo / Avenza / QGIS / ArcGIS.

3. With `--streams`, the **full-resolution NHD stream network** for the AOI as a
   vector layer (`<key>_streams.geojson` by default) — the same dense NHDPlus HR
   flowlines the figures' reference overlay draws (named creeks **and** unnamed
   headwaters), clipped to the AOI polygon, in EPSG:4326. Import straight into
   CalTopo / a GIS.

By default the GeoPDFs render in **UTM** (identical look to the PNG deliverables;
the georeferencing is exact in whatever projected CRS they render in) and frame to
the event AOI exactly like the `map`/`run` i15 maps (auto-matched like `climate`).
The anomaly GeoPDF needs `<from-key>_anom_i{d}.tif` from a prior [`climate`](#climate--noaa-atlas-14-climatology-vs-observed)
run. **The GeoPDF half needs GDAL's PDF driver** (`conda install -c conda-forge
libgdal-pdf`); without it the GeoTIFF export still runs and the PDFs
are skipped with a note. Groups **A + C**, plus:
| Arg | Type / default | Meaning |
|---|---|---|
| `--from-dir` | path, **required** | directory of the processed event (its `--out-dir`) |
| `--from-key` | str, **required** | key/stem of the processed event |
| `--layers FIELD…` | strs, `anom_i15 i15max` | field suffixes to export as GeoTIFF (e.g. add `i30max i60max total peakrate_mmph`) |
| `--crs` | str, `EPSG:3857` | GeoTIFF export CRS (CalTopo's native Web-Mercator) |
| `--pdf-crs` | str, `UTM` | projected CRS to render the GeoPDFs in (default UTM = the PNG look; pass `EPSG:3857` to match the layers) |
| `--anom-duration MIN` | int, `15` | anomaly duration for the GeoPDF (needs `<from-key>_anom_i<MIN>.tif`) |
| `--ari` | int, `1` | climatology ARI (years) labelling the anomaly |
| `--anomaly-cmap` | str, `cmc.vik` | diverging colormap for the anomaly |
| `--hillshade` | path | hillshade (default `<from-dir>/<from-key>_hillshade.tif`) |
| `--gauges` | flag | overlay the source `<from-key>_gauges.geojson` on the GeoPDFs |
| `--streams` | flag | also export the full-resolution NHD stream network for the AOI as a vector layer |
| `--streams-format {geojson,gpkg,shp,kml}` | `geojson` | vector format for `--streams` |
| `--streams-bbox` | flag | keep whole flowlines over the AOI bounding box (default: clip to the AOI polygon) |
| `--streams-named-only` | flag | keep only GNIS-named creeks/rivers (default: the full network incl. unnamed headwaters) |
| `--no-geotiff` | flag | skip the EPSG:3857 GeoTIFF layers |
| `--no-rgb` / `--no-raw` | flags | of each layer, skip the colorized RGBA / the raw float GeoTIFF |
| `--no-figures` | flag | skip the GeoPDFs (export only the GeoTIFF layers) |
| `--no-i15` / `--no-anom` | flags | skip the i15-map / anomaly-map GeoPDF |

```bash
# CalTopo-ready 3857 layers + georeferenced PDFs + the NHD stream network
python -m stormscape export --from-dir ./out --from-key HiddenValley \
    --out-dir ./out --key HiddenValley \
    --aoi HiddenValley_AOI.kmz --gauges --reference --streams
```

---

### `burn` — near-real-time burn severity (CIMSS BRISK dNBR)

Maps the **burn scar** over an AOI from **CIMSS BRISK**: a daily, nine-satellite
dNBR composite (GOES ABI, VIIRS, Landsat 8/9, Sentinel-2) covering every large US
fire, so severity is available *while the fire is still burning*. Scenes are
per-fire, per-day GeoTIFFs in the open SSEC archive; this command screens the
archive for fires intersecting the AOI (reading only GeoTIFF *headers* over HTTP
range requests — it downloads nothing until it knows what it needs), caches the
scenes it does need, mosaics them, classifies severity and draws the map.

> **dNBR is a *vegetation* index, not soil burn severity.** The USGS post-fire
> debris-flow models are calibrated on **soil** burn severity — dNBR adjusted by
> field crews for hydrophobicity, ground cover and duff consumption. BRISK is
> explicitly **interim**: act on it early, then supersede it with `--product sbs`
> when the BAER assessment lands.

Groups **A + C**, plus:
| Arg | Type / default | Meaning |
|---|---|---|
| `--date YYYYMMDD` | str | the scar **as of** this date — each fire is shown at its latest scene on or before it (default: latest available) |
| `--since YYYYMMDD` | str | ignore scenes older than this (e.g. to drop last season's fires) |
| `--product {dnbr,sbs,baer_dnbr}` | `dnbr` | `dnbr` = BRISK daily composite; `sbs` = BAER **soil** burn severity, authoritative but only for assessed fires; `baer_dnbr` = the BAER teams' own dNBR (2025 only so far) |
| `--min-age DAYS` | float | require the composite to be at least DAYS old (since the fire entered the archive) and skip fires that are not, naming them. Off by default — an immature scar still beats none, and the run says so |
| `--scheme {usgs,brisk}` | `usgs` | dNBR severity breaks: USGS/MTBS `0.10/0.27/0.44/0.66`, or the portal's `0.10/0.40/0.70` |
| `--fire NAME…` | strs | restrict to these fires (names as printed by `--list`) |
| `--years Y…` | ints | archive years to search (default: this year and last) |
| `--list` | flag | list the fires intersecting the AOI and stop — **downloads nothing** |
| `--all-dates` | flag | with `--list`, show every scene per fire instead of only the latest |
| `--cache-dir` | path | scene + index cache (default `<out-dir>/brisk_cache`) |
| `--hillshade` | path | hillshade for the map backdrop (default `<out-dir>/<key>_hillshade.tif`) |
| `--dem` | flag | fetch a DEM + hillshade for the AOI if none is found |
| `--resolution` | int, `10` | DEM resolution (m) for `--dem` |
| `--vmax` | float | colour-scale max (default `1.0` for dNBR, `4.0` for the `sbs` class field) |
| `--no-map` | flag | write the rasters + table only, no figure |
| `--workers` | int, `12` | parallel header reads / downloads |

| `--continuous` | flag | shade dNBR as a smooth ramp instead of banding it into severity classes (same colours) |

Maps default to the **BAER class palette**, banded into the severity classes BAER
publishes, with the colour bar labelled by class name rather than dNBR value — so a
stormscape map can be laid beside a BAER product. `--continuous` gives the smooth
ramp in the same colours, an explicit `--cmap` (e.g. `YlOrRd`) opts out of the class
colours entirely, and `--alpha` makes the scar read more strongly than the
project-wide 0.32 drape.

```bash
# what has burned here? (headers only, no downloads)
python -m stormscape burn --aoi event_AOI.kmz --list

# fetch, cache, classify and map it
python -m stormscape burn --aoi event_AOI.kmz --date 20260814 \
    --out-dir ./out --key event --dem --reference --clip
```

---

## Outputs at a glance

| Command | Writes |
|---|---|
| `dem` | `<key>_dem.tif`, `<key>_hillshade.tif` |
| `i15` | `<key>_i15max.tif`, `_i30max`, `_i60max`, `_i2max`, `_total`, `_tpki15`, `_rqi`, `_shsr` (+ `_mstotal.tif` with `--multisensor`) |
| `map` / `run` | draped figure PNG (`run` also writes the `dem`+`i15` rasters; `--compare` adds `<key>_compare.csv` + residual map) |
| `gauges` | `<key>_gauges.geojson` (coords + peaks + I15 time-of-peak) + `RainGaugeData/<key>_gauge_<name>.csv` series (unless `--no-series`) |
| `compare` | per-gauge CSV (`--out`) + optional residual map (`--map`) |
| `nexrad` | `<key>_refl.tif` + `<key>_nexrad.png` (or `<key>_i15max.tif` … with `--intensity`); volumes cached under `--cache-dir` |
| `panels` | `<key>_panels.png` |
| `vgauge` | per-gauge CSVs in `RainGaugeData/`; `<key>_vg_atlas.png` (`--atlas`); per-gauge figures in `VirtualGaugeFigures/` (`--detail`) |
| `burn` | `<key>_dnbr.tif`, `<key>_severity.tif`, `<key>_burn_classes.csv` (pixels + true km² + fraction per class), `<key>_burn_scenes.geojson` (which fire, which date), `<key>_burn.png`; scenes cached under `brisk_cache/` |
| `zoom` | `<key>.png` (zoomed map) + `<key>_panels.png` + (unless `--no-climate`) `<key>_climate_compare.png`, `<key>_clim_i{d}.tif`, `<key>_anom_i{d}.png`/`.tif`; `--crop-rasters` adds cropped `<key>_<field>.tif`; `--refine-dem` adds `<key>_dem/_hillshade.tif` |
| `pick` | `<from-key>_pick.html` (self-contained browser bbox picker) |
| `climate` | `<key>_clim_i{15,30,60}.tif`, `<key>_climate_compare.png`, `<key>_anom_i{15,30,60}.tif` + `.png`; grids cached under `atlas14_cache/` |
| `smooth` | `<key>_smoothing_compare.png`; with `--gauge-analysis` also `<key>_smoothing_skill.png` + `.csv`; with `--write` smoothed `<key>_<field>.tif` |
| `recurrence` | `<key>_gauge_recurrence.csv` + `.md` (per wet gauge: peak I15/30/60, time-of-peak, anomaly, recurrence interval) |
| `export` | `<key>_<field>_3857.tif` (raw) + `_3857_rgb.tif` (colorized) per `--layers`; `<key>.pdf` + `<key>_anom_i{d}.pdf` GeoPDFs (needs `libgdal-pdf`); with `--streams`, `<key>_streams.geojson` (NHD network) |
