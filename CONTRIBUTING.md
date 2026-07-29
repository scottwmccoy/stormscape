# Contributing to stormscape

Welcome — this guide is written for **students and researchers** joining the
project, whether you want to fix a bug, add a feature, or just get it running on
your own storm event. No prior experience with this codebase is assumed.

If you only want to *use* stormscape, start with the [README](README.md) and the
[CLI reference](docs/cli.md). Come back here when you want to change something.

---

## Three rules that matter more than the rest

**1. Never commit an API token.** The Synoptic / MesoWest gauge API needs a
personal token. It belongs in an environment variable, never in a file you commit:

```bash
export SYNOPTIC_TOKEN=your_token_here
```

Every gauge command reads `$SYNOPTIC_TOKEN` (or takes `--token`). If you ever
paste a token into a script, a notebook, or a commit message, **assume it is
compromised** — git history is permanent, and deleting it in a later commit does
not remove it. Tell the maintainer so the token can be revoked and reissued.

**2. Never commit data.** This tool downloads gigabytes. A single event's rasters
run ~200 MB, and the request cache reached 4.6 GB during development. The
[`.gitignore`](.gitignore) already excludes `*.tif`, `cache/`, `*_cache/`,
`RainGaugeData/`, `scratchpad/`, and friends — please don't override it with
`git add -f`. Before committing, check what you're about to send:

```bash
git status
git diff --cached --stat
```

If you see a `.tif`, a cache directory, or anything measured in MB, stop and ask.
Small example figures (a PNG under a megabyte or so) are fine and welcome.

**3. Run the tests before you push.** They take about two seconds.

```bash
pytest
```

---

## Getting set up

The geospatial stack (GDAL, rasterio, geopandas, pyproj) installs far more
reliably from **conda-forge** than from pip. Use the provided environment:

```bash
git clone https://github.com/scottwmccoy/stormscape.git
cd stormscape
conda env create -f environment.yml     # creates an env named "stormscape"
conda activate stormscape
pip install -e ".[test]"                # editable install + test dependencies
```

`-e` (editable) means your edits take effect immediately — no reinstalling.

### Check it worked

```bash
pytest                                  # expect ~232 passing in a couple of seconds
python -m stormscape --help             # should list 15 subcommands
```

MRMS radar files are GRIB2, which needs GDAL's GRIB driver. The conda environment
includes it; verify with:

```python
import rasterio
assert "GRIB" in rasterio.drivers.raster_driver_extensions().values()
```

If that assertion fails, your GDAL lacks the GRIB driver — `conda install -c
conda-forge libgdal-grib` inside the environment. (GeoPDF export separately needs
`libgdal-pdf`; it's optional and skips itself when absent.)

### Already have a working geospatial environment?

Just `pip install -e ".[test]"` into it. You need Python **3.10 or newer**.

---

## Making a change

**If you don't have push access** (the default), fork the repository on GitHub,
then:

```bash
git clone https://github.com/YOUR-USERNAME/stormscape.git
cd stormscape
git remote add upstream https://github.com/scottwmccoy/stormscape.git
git checkout -b my-change
# ... edit, test ...
git push origin my-change
```

Then open a Pull Request against `scottwmccoy/stormscape` on GitHub.

**If you do have push access**, branch directly — but please don't commit to
`main`:

```bash
git checkout -b my-change
# ... edit, test ...
git push origin my-change
```

Keep your branch current with `git fetch upstream && git rebase upstream/main`
(or `git pull --rebase origin main` if you're working in the main repo).

### What happens next

Continuous integration runs the test suite on Python 3.10 and 3.13 for every pull
request. A green check means your change didn't break anything covered by tests;
a red X will tell you which test failed and why. Fix it and push again — the PR
updates automatically.

---

## What makes a good contribution here

This is scientific software, so **correctness beats cleverness**. A few
expectations:

- **Add a test for anything with a right answer.** If you can state the expected
  output — "a constant 60 mm/h rate must give i15 = 60" — encode it. See the
  worked example below.
- **Say why in comments, not what.** The code shows what it does; the interesting
  part is the physical or practical reason. Look at
  [`stormscape/mrms.py`](stormscape/mrms.py) for the tone.
- **Render the figure and look at it.** The tests cover the math and the file
  plumbing, but they cannot tell you a map looks wrong. If you touch anything in
  `plot.py`, generate a figure and inspect it before pushing.
- **Keep the docs in step.** A new flag needs an entry in
  [`docs/cli.md`](docs/cli.md); a notable new capability needs a README mention.
- **Update `CLAUDE.md`** if you learn something non-obvious — a gotcha, a failed
  approach and why it failed, a validated result. That file is the project's
  accumulated hard-won knowledge, and a dead end you document saves the next
  person days.

Small, focused pull requests get reviewed faster than large ones. If you're
planning something substantial, open an issue first so we can talk about the
approach before you invest the effort.

---

## Project conventions you need to know

These aren't style preferences — breaking them produces wrong or misleading
science.

**Colormaps must be perceptually uniform and colorblind-safe.** Never `jet` or
`turbo`: they create false edges and rank badly for the ~8% of readers with
colour-vision deficiency. Rainfall fields default to `YlGnBu`, anomalies to a
diverging map centred at 1× (`cmc.vik`), residuals to `RdBu_r`.

**`i15`, `i30`, `i60` are rainfall metrics, not the old package name.** They mean
peak 15/30/60-minute rainfall intensity. `drape_i15`, `i15max`, `i15_storm_day`
are all correct as written — don't "finish" the `i15toolkit` → `stormscape` rename
by touching them.

**Figures are drawn in UTM.** CONUS Albers (EPSG:5070) is too convergent in the
western US for straight lat/long ticks or an upward-pointing north arrow. DEMs are
*stored* in 5070 and MRMS rasters in 4326, but maps display in an auto-selected
UTM zone.

**Heavy optional dependencies are imported lazily**, inside the function that
needs them — see `dem._py3dep()` or the `arm_pyart` imports in `nexrad.py`. This
keeps `import stormscape` fast and makes a missing package produce a clear
"install this" message instead of breaking the whole package. If you add a
dependency that isn't needed at import time, follow that pattern.

**Always inspect the RQI field for quantitative work.** Far from the radar or
behind terrain, the beam overshoots low rain and the intensity estimate is
unreliable. `compare --rqi-min` exists for this reason.

**Radar is not truth.** On the Hidden Valley test event the radar read roughly
2–4× the ground gauges, because hail inflates reflectivity-based rainfall
estimates. That's real physics, not a bug. Validate against gauges before
drawing conclusions from a radar field.

**Test packaging changes in a clean virtual environment.** Two `pip install`-
breaking bugs once hid behind a fully-populated conda env. If you touch
`pyproject.toml`:

```bash
python -m venv /tmp/checkenv
/tmp/checkenv/bin/pip install .        # no extras — must succeed
/tmp/checkenv/bin/python -c "import stormscape"
```

---

## Adding a test — worked example

Tests live in `tests/`, one file per module, and run offline: **no network, no
API token.** Synthetic inputs with a known answer, via the fixtures in
`tests/conftest.py` (`field_tif`, `points_gdf`, `gauge_obs`, `minute_series`).

Say you want to pin down that smoothing can't invent rainfall:

```python
# tests/test_smoothing.py
def test_smoothing_never_raises_the_peak():
    """Low-pass filters cannot amplify — this underpins the caveat that smoothing
    lowers the radar's positive i15 bias mechanically, not by improving skill."""
    rng = np.random.RandomState(0)
    field = rng.rand(41, 41) * 100.0
    out = smoothing.smooth_array(field, cell_km=1.0, method="gaussian",
                                 radius_km=2.0)
    assert np.nanmax(out) <= np.nanmax(field) + 1e-9
```

Note the docstring explains *why the invariant matters*, not what the code does.
That's the house style: a failing test should tell you what broke scientifically.

Useful commands while developing:

```bash
pytest tests/test_smoothing.py            # one file
pytest -k "recurrence"                     # tests matching a name
pytest -x                                  # stop at the first failure
pytest -q                                  # quiet summary
```

If your test needs an optional dependency, mark it and let it skip cleanly rather
than fail on machines that lack it:

```python
@pytest.mark.optional_deps
def test_something_with_radar_volumes():
    pyart = pytest.importorskip("pyart")     # skips if Py-ART isn't installed
    ...
```

Anything requiring the network gets `@pytest.mark.network`; CI deselects those
with `-m "not network"`. Both markers are declared in `pyproject.toml`, and
`--strict-markers` means a typo in a marker name fails loudly instead of being
silently ignored.

---

## Where things live

```
stormscape/
├── aoi.py        area-of-interest parsing (bbox / vector / geometry)
├── dem.py        USGS 3DEP DEM download + hillshade
├── mrms.py       MRMS mosaic radar -> peak i15/i30/i60 fields
├── nexrad.py     single-radar NEXRAD Level II (incl. pre-2020 events)
├── gauges.py     Synoptic rain gauges -> totals + peak intensities
├── compare.py    radar vs gauges: residuals, skill stats, recurrence
├── merge.py      bias correction and radar-gauge merging
├── atlas14.py    NOAA Atlas 14 climatology + recurrence intervals
├── smoothing.py  NaN-aware field smoothing + gauge-skill sweep
├── export.py     EPSG:3857 GeoTIFFs, GeoPDFs, NHD stream vectors
├── refdata.py    NHD streams / TIGER roads / GNIS place names
├── plot.py       terrain-draped maps, panels, virtual-gauge figures
└── cli.py        the 15 subcommands
tests/            offline test suite (start here to understand behaviour)
docs/cli.md       full flag reference for every subcommand
examples/         a worked Hidden Valley example + batch templates
CLAUDE.md         accumulated conventions, gotchas, and validated findings
```

**Reading the tests is the fastest way to learn what a function actually
guarantees** — that's largely why they exist.

---

## Reporting a bug or requesting a feature

Open an issue at
[github.com/scottwmccoy/stormscape/issues](https://github.com/scottwmccoy/stormscape/issues).
A useful bug report includes:

- the **exact command** you ran (with your token redacted);
- the **AOI and date**, so it can be reproduced;
- the full **error message and traceback**;
- your OS and the output of `conda list | grep -E "gdal|rasterio|geopandas"`.

Remote data services are a common culprit — 3DEP's 1 m endpoint is genuinely slow
and flaky, MRMS returns 404 for missing intervals (which means "no data", never
"retry forever"), and NHD's endpoint sometimes returns nothing. If a failure
mentions a timeout or an empty result, try again before assuming the code is
broken.

---

## Getting help

Ask. Opening an issue with a half-formed question is completely fine, and so is
asking Scott directly. If something in this guide is unclear or wrong, that's a
documentation bug worth reporting too.

By contributing, you agree your work is licensed under the
[MIT License](LICENSE), the same terms as the rest of the project.
