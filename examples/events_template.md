# Event inventory — CSV template

One **row = one storm event**. The CSV defines *what / where / when* for each
event; the analysis options that are the same across the whole batch (e.g.
`--reference`, `--clip`, `--gauges`, `--compare`, `--multisensor`, `--basemap`,
`--dpi`) are applied **once** by the batch runner rather than repeated as columns.
The handful of options that genuinely vary per event get their own columns below
(`source`, `radar`, `method`, `resolution`, `out_dir`).

Blank cell = use the default. Quote any field that contains a comma (e.g.
`name`, `notes`). The file is plain CSV — edit it in Excel/Sheets or any editor,
and `pandas.read_csv` reads it directly.

`examples/events_template.csv` ships with two rows: a **real** Hidden Valley row
(MRMS era, AOI from a KMZ) and an **EXAMPLE** pre-2020 row (single-radar NEXRAD,
`--bbox` form) — keep, edit, or delete them.

## Columns

| Column | Required? | Maps to | Meaning |
|---|---|---|---|
| `key` | **yes** | `--key` | short event id / output filename stem (also the per-event output subfolder) |
| `name` | optional | `--title` | human-readable label used as the figure title |
| `date` | **yes**\* | `--date` | storm day `YYYYMMDD` (drives MRMS wet-hour selection + the default gauge window) |
| `start` | optional | `--start` | explicit UTC window start `YYYYMMDDHHMM` — **required for `source=nexrad`** and to narrow windows |
| `end` | optional | `--end` | explicit UTC window end `YYYYMMDDHHMM` |
| `aoi` | AOI\*\* | `--aoi` | path to an AOI vector file (KMZ/GeoJSON/SHP/GPKG) |
| `bbox_w` `bbox_s` `bbox_e` `bbox_n` | AOI\*\* | `--bbox W S E N` | AOI bounding box in lon/lat degrees (use these **or** `aoi`) |
| `out_dir` | optional | `--out-dir` | output directory; blank → batch runner uses `<out-root>/<key>` |
| `perimeters` | optional | `--perimeters` | outline drawn on the figure (often the same file as `aoi`) |
| `resolution` | optional | `--resolution` | DEM resolution, metres (default `10`) |
| `source` | optional | — | `mrms` (default; ≥ Oct 2020) or `nexrad` (single-radar Level II i15 stack, for pre-2020 / MRMS-gap events) |
| `radar` | optional | `--radar` | NEXRAD 4-letter id for `source=nexrad` (blank → nearest to the AOI) |
| `method` | optional | `--method` | NEXRAD rate recipe: `kdp` (default, dual-pol, 2012+) or `za` (capped convective Z–R, all eras) |
| `notes` | optional | — | free text (event description, perimeter provenance, caveats) |

\* `date` may be omitted only if both `start` and `end` are given.
\*\* Each row needs **either** `aoi` **or** all four `bbox_*` values.

## What a row becomes

The Hidden Valley row is equivalent to:

```bash
python -m stormscape run \
  --aoi ./aoi/HiddenValley_AOI.kmz \
  --date 20260619 --resolution 10 \
  --perimeters ./aoi/HiddenValley_AOI.kmz \
  --out-dir ".../2026_06_19_HiddenValley_Vista/Raster_Data" --key HiddenValley \
  --title "Hidden Valley, NV — peak i15 (2026-06-19)" \
  --reference --clip            # <- batch-wide analysis flags, applied to every row
```

The pre-2020 example row, because `source=nexrad`, would instead build the
single-radar NEXRAD Level II i15 stack over its `--start`/`--end` window
(`--method kdp`) and drape that.

## Running a batch

A `batch` subcommand that consumes this CSV is **not built yet** — once you've
populated a handful of events, the intended call is roughly:

```bash
python -m stormscape batch --inventory examples/events_template.csv \
  --out-root ./events --reference --clip --gauges --compare
```

i.e. the inventory supplies the per-event definition, and the trailing flags are
the analysis profile applied to all of them. Tell me when the CSV has a few real
rows and I'll build the runner against this exact column contract.
