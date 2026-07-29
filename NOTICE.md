# Attributions and data sources

`stormscape` is released under the [MIT License](LICENSE). This file records
third-party code incorporated in the project and the terms of the public data
it accesses.

## Third-party code

| Where | Origin | Terms |
|---|---|---|
| `stormscape/gauges.py` — Synoptic/MesoWest gauge transport and the 1-minute precipitation interpolation | USGS **FlowAlert** | CC0 1.0 Universal (public domain dedication) |
| `stormscape/mrms.py` — the peak-i15 stacking estimator | D. Cavagna, `MRMS_stack.py` | Used with permission |

## Data sources

All of the following are **public-domain US government data**, accessed at run
time over public endpoints. No credentials are required and none are bundled.

| Data | Agency / product |
|---|---|
| Digital elevation models | USGS **3DEP** / The National Map |
| Mosaic radar QPE (`PrecipRate`, MultiSensor QPE, RQI) | NOAA **MRMS** |
| Single-radar volumes | NOAA **NEXRAD Level II** (via the `unidata-nexrad-level2` AWS bucket) |
| Precipitation-frequency climatology | NOAA **Atlas 14** (HDSC gridded + PFDS point) |
| Stream network | USGS **NHD / NHDPlus HR** |
| Roads | US Census **TIGER/Line** |
| Place names | USGS **GNIS** |

### One exception — rain gauges

Ground rain-gauge observations come from the **Synoptic Data / MesoWest** API,
which is *not* public domain. It requires **your own API token** (free for
academic and research "open access" use — request one at
[synopticdata.com](https://synopticdata.com)) and its use is governed by
Synoptic Data's terms of service.

Pass your token via the `SYNOPTIC_TOKEN` environment variable or the `--token`
flag. **Never commit a token to this or any repository.**

## Citing this software

If `stormscape` contributes to published work, please cite the repository:

> McCoy, S. W. (2026). *stormscape: storm rainfall from public US radar, gauges,
> and climatology, analyzed and mapped over terrain.*
> https://github.com/scottwmccoy/stormscape

Please also cite the underlying data products you relied on (MRMS, NEXRAD,
Atlas 14, 3DEP, Synoptic) as those agencies request.
