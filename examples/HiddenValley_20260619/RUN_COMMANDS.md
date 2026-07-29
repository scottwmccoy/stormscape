# Worked example — Hidden Valley, NV (storm of 2026-06-19)

A real end-to-end run, reduced to commands you can copy. It produces a 10 m DEM
+ hillshade, the MRMS peak-intensity fields, and the **production figure**
(`HiddenValley.png`, in this folder) with labelled vector reference overlays
(NHD streams, TIGER roads, GNIS place names), clipped to the AOI.

**Inputs you supply:** an AOI vector file (KMZ/GeoJSON/SHP/GPKG) or a `--bbox`,
and the storm date. This example's AOI covers
W −119.7406, S 39.3304, E −119.4536, N 39.8777 (≈ 24.5 × 60.7 km) — so you can
reproduce it with `--bbox` alone if you have no vector file.

---

## 1. The whole pipeline in one command

```bash
stormscape run \
  --aoi ./HiddenValley_AOI.kmz \
  --date 20260619 \
  --resolution 10 \
  --out-dir ./out --key HiddenValley \
  --perimeters ./HiddenValley_AOI.kmz \
  --reference --clip \
  --title "Hidden Valley, NV  -  peak i15  (2026-06-19)"
```

No AOI file? Use the bounding box instead:

```bash
stormscape run \
  --bbox -119.7406 39.3304 -119.4536 39.8777 \
  --date 20260619 --resolution 10 \
  --out-dir ./out --key HiddenValley --reference --clip
```

If you have not installed the package, you can run it from a source checkout
without installing — from the repository root:

```bash
python -m stormscape run --bbox -119.7406 39.3304 -119.4536 39.8777 \
  --date 20260619 --out-dir ./out --key HiddenValley
```

## 2. Add rain gauges and validate the radar (needs a Synoptic token)

```bash
export SYNOPTIC_TOKEN=your_token_here     # free for academic use
stormscape gauges --aoi ./HiddenValley_AOI.kmz --date 20260619 \
  --out-dir ./out --key HiddenValley --pad-deg 0.12 --max-report-min 60
```

`--pad-deg 0.12` is deliberate here: this AOI is long and narrow, and the wettest
gauges of this storm sat just *outside* it (one of them recorded a ~200-year
15-minute intensity). The default 0.05° pad would silently miss them.
`--max-report-min 60` drops daily-reporting gauges, whose short-duration peaks
are smeared low by the 1-minute interpolation.

## 3. Climatological context, then export for GIS / CalTopo

```bash
stormscape climate  --from-dir ./out --from-key HiddenValley \
                    --out-dir ./out --key HiddenValley
stormscape export   --from-dir ./out --from-key HiddenValley \
                    --out-dir ./out --key HiddenValley --streams
```

`climate` writes the NOAA Atlas 14 comparison + anomaly maps (observed ÷ the
1-year storm). `export` writes EPSG:3857 GeoTIFFs, GeoPDFs of the two primary
figures, and the full-resolution NHD stream network as a vector layer.

---

## Flags used above

| flag | meaning |
|---|---|
| `run` | DEM → hillshade → MRMS i15 → figure (whole pipeline) |
| `--aoi` / `--bbox` | area of interest: vector file, or `W S E N` in lon/lat degrees |
| `--date` | storm day, `YYYYMMDD` |
| `--resolution` | DEM resolution in metres (10 here) |
| `--perimeters` | vector outline drawn on the figure |
| `--reference` | overlay labelled NHD streams + TIGER roads + GNIS places |
| `--clip` | crop the figure tightly to the AOI (`--clip-margin` sets padding) |
| `--key` / `--out-dir` | output filename stem / directory |
| `--pad-deg` | degrees to pad the AOI before fetching (gauges, radar) |

Other figure flags: `--local-roads`, `--no-reference-labels`, `--dpi`
(default 200), `--basemap` (USGS topo tiles instead of vectors). DEM masking is
`--clip-dem`, kept distinct from the figure-cropping `--clip`. Full reference:
[`docs/cli.md`](../../docs/cli.md).

## Outputs

`HiddenValley.png` (the figure in this folder) plus GeoTIFFs — `_i15max` (peak
15-min intensity, mm/h), `_i30max`, `_i60max`, `_i2max`, `_total`, `_tpki15`
(time of peak), `_rqi` (radar quality), `_shsr` (EPSG:4326, native MRMS grid),
and `_dem` / `_hillshade` (10 m, EPSG:5070). Roughly 200 MB in total, which is
why the `.gitignore` keeps `*.tif` out of the repository.

## What this storm looked like

Radar quality was high across the AOI (**RQI median 1.00**, min 0.90). Inside the
AOI polygon the peak **i15 reached ≈ 102 mm/h** (mean 34, p90 63), the peak 2-min
rate ≈ 125 mm/h, and the storm total up to 70 mm, peaking ≈ 22:48 UTC
(≈ 3:48 PM PDT) — a short, intense convective cell.

Worth knowing before you trust radar intensities at face value: for this storm
the radar read roughly **2–4× higher than the ground gauges**, because hail
inflates reflectivity-based rainfall estimates. That discrepancy is real physics,
not a bug — see the radar-vs-gauge and `smooth` sections of the README, and use
`compare` to quantify it for your own event.
