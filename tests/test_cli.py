"""CLI wiring: every subcommand parses, and the flags that guard behaviour exist.

These are cheap but catch the most common regression when adding a command — a
parser that references a missing attribute, or an option quietly dropped from one
command while the docs still advertise it.
"""
from __future__ import annotations

import pytest

from stormscape import cli

SUBCOMMANDS = ["dem", "i15", "map", "run", "gauges", "compare", "nexrad",
               "panels", "vgauge", "zoom", "pick", "climate", "smooth",
               "recurrence", "export"]


@pytest.mark.parametrize("cmd", SUBCOMMANDS)
def test_every_subcommand_builds_its_help(cmd, capsys):
    """--help exercises the whole parser for that command."""
    with pytest.raises(SystemExit) as e:
        cli.main([cmd, "--help"])
    assert e.value.code == 0
    assert cmd in capsys.readouterr().out


def test_no_subcommand_is_an_error():
    with pytest.raises(SystemExit) as e:
        cli.main([])
    assert e.value.code != 0


def test_unknown_subcommand_is_an_error():
    with pytest.raises(SystemExit) as e:
        cli.main(["definitely-not-a-command"])
    assert e.value.code != 0


def _opts(cmd):
    """The set of long options a subcommand accepts."""
    with pytest.raises(SystemExit):
        cli.main([cmd, "--help"])
    return cmd


@pytest.mark.parametrize("cmd, flag", [
    # the project-wide drape opacity must be overridable on every map command
    ("map", "--alpha"), ("run", "--alpha"), ("nexrad", "--alpha"),
    ("zoom", "--alpha"), ("climate", "--alpha"), ("panels", "--alpha"),
    ("compare", "--alpha"), ("smooth", "--alpha"), ("export", "--alpha"),
    # gauge cadence screening, the fix that keeps daily gauges out of i15 stats
    ("gauges", "--max-report-min"), ("vgauge", "--max-report-min"),
    ("run", "--max-report-min"), ("compare", "--max-report-min"),
    ("smooth", "--max-report-min"),
    # export products
    ("export", "--layers"), ("export", "--streams"), ("export", "--pdf-crs"),
    ("export", "--streams-format"), ("export", "--no-figures"),
    # climate / smoothing knobs
    ("climate", "--obs-smooth"), ("zoom", "--obs-smooth"),
    ("smooth", "--gauge-analysis"), ("smooth", "--write"),
    # DEM warp resampling — on every command that fetches a DEM
    ("dem", "--resampling"), ("run", "--resampling"), ("zoom", "--resampling"),
    # product layout
    ("dem", "--flat"), ("run", "--flat"), ("gauges", "--flat"),
])
def test_documented_flag_is_present(cmd, flag, capsys):
    with pytest.raises(SystemExit):
        cli.main([cmd, "--help"])
    assert flag in capsys.readouterr().out, f"{cmd} lost {flag}"


@pytest.mark.parametrize("cmd", ["climate", "smooth", "export", "recurrence",
                                 "zoom"])
def test_reuse_commands_require_a_source_event(cmd, capsys):
    """The re-render commands read an already-processed event, so --from-dir and
    --from-key are mandatory; forgetting them must fail fast, not half-run."""
    with pytest.raises(SystemExit):
        cli.main([cmd, "--help"])
    out = capsys.readouterr().out
    assert "--from-dir" in out and "--from-key" in out


def test_aoi_is_required_when_no_bbox_or_aoi_given(monkeypatch):
    """_aoi_from_args exits rather than silently defaulting to somewhere."""
    class Args:
        bbox = None
        aoi = None
    with pytest.raises(SystemExit):
        cli._aoi_from_args(Args())


def test_aoi_from_args_prefers_bbox_then_aoi():
    class WithBbox:
        bbox = [-119.7, 39.3, -119.4, 39.9]
        aoi = "ignored.kmz"
    assert cli._aoi_from_args(WithBbox()) == (-119.7, 39.3, -119.4, 39.9)

    class WithAoi:
        bbox = None
        aoi = "x.kmz"
    assert cli._aoi_from_args(WithAoi()) == "x.kmz"


def test_render_px_scales_with_dpi_but_has_a_floor():
    """The hillshade render cap must never under-sample a high-dpi figure, and must
    stay high enough at low dpi to cover the widest panel."""
    assert cli._render_px(200) >= 2500
    assert cli._render_px(600) > cli._render_px(200)
    assert cli._render_px(None) >= 2500


def test_find_event_aoi_locates_a_saved_aoi(tmp_path):
    """climate/export auto-match their extent to the event AOI written by dem/i15."""
    import geopandas as gpd
    from stormscape.aoi import bbox_polygon
    gpd.GeoDataFrame(geometry=[bbox_polygon((-119.7, 39.3, -119.4, 39.9))],
                     crs=4326).to_file(tmp_path / "ev_aoi.geojson", driver="GeoJSON")
    assert cli._find_event_aoi(str(tmp_path), "ev") is not None


def test_find_event_aoi_returns_none_when_absent(tmp_path):
    assert cli._find_event_aoi(str(tmp_path), "nothing_here") is None


def test_labels_normalises_the_off_switches():
    class A:
        basemap_labels = "none"
    assert cli._labels(A()) is None

    class B:
        basemap_labels = "CartoDB.PositronOnlyLabels"
    assert cli._labels(B()) == "CartoDB.PositronOnlyLabels"


# --------------------------------------------------------------------------- #
# --resampling: it has to REACH the DEM fetch, not merely parse
# --------------------------------------------------------------------------- #
def _dem_kwargs(monkeypatch, extra):
    """Run `dem` with the fetch stubbed; return the kwargs it was called with."""
    import stormscape.dem as d
    seen = {}

    def fake(aoi, **kw):
        seen.update(kw)
        raise SystemExit(0)

    monkeypatch.setattr(d, "fetch_dem_and_hillshade", fake)
    with pytest.raises(SystemExit):
        cli.main(["dem", "--bbox", "-119.6", "39.7", "-119.5", "39.8",
                  "--out-dir", "/tmp/_ss_test"] + extra)
    return seen


def test_resampling_defaults_to_none_so_the_library_owns_the_default(monkeypatch):
    """The parser must not pin the default — dem.DEFAULT_RESAMPLING does."""
    assert _dem_kwargs(monkeypatch, [])["resampling"] is None


@pytest.mark.parametrize("name", ["nearest", "bilinear", "cubic", "lanczos"])
def test_resampling_is_forwarded_to_the_dem_fetch(monkeypatch, name):
    assert _dem_kwargs(monkeypatch, ["--resampling", name])["resampling"] == name


def test_unknown_resampling_is_rejected_by_the_parser(capsys):
    with pytest.raises(SystemExit) as e:
        cli.main(["dem", "--bbox", "-119.6", "39.7", "-119.5", "39.8",
                  "--resampling", "banana"])
    assert e.value.code != 0
    assert "invalid choice" in capsys.readouterr().err


def test_resampling_offers_nearest_so_the_artefact_stays_reproducible(capsys):
    """Reproducing the corduroy hatch is the control that pins it on *nearest*
    rather than on the number of warps — keep it reachable from the CLI."""
    with pytest.raises(SystemExit):
        cli.main(["dem", "--help"])
    assert "nearest" in capsys.readouterr().out


# --------------------------------------------------------------------------- #
# --start/--end scope the radar stack, not only the gauges
# --------------------------------------------------------------------------- #
class _Args:
    def __init__(self, **kw):
        self.date = self.start = self.end = None
        self.__dict__.update(kw)


def test_date_alone_means_the_storm_day_scan():
    assert cli._stack_window(_Args(date="20260814")) is None


def test_start_end_alone_is_accepted_and_used():
    w = cli._stack_window(_Args(start="202608142000", end="202608150400"))
    assert (w[0].hour, w[1].hour) == (20, 4)


def test_start_end_override_date_for_the_stack():
    """One pair of flags means "the analysis window" everywhere. They used to
    scope the gauges while the radar quietly stacked the whole day."""
    w = cli._stack_window(_Args(date="20260101", start="202608142000",
                                end="202608150400"))
    assert w[0].strftime("%Y%m%d") == "20260814"


def test_neither_date_nor_window_is_rejected():
    with pytest.raises(SystemExit):
        cli._stack_window(_Args())


def test_event_label_falls_back_to_the_window_start():
    assert cli._event_label(_Args(start="202608142000",
                                  end="202608150400")) == "20260814"


@pytest.mark.parametrize("cmd", ["i15", "run"])
def test_date_is_no_longer_mandatory(cmd, capsys):
    """--date must be optional so a window alone works; the parser should not
    reject it before _stack_window can validate the combination."""
    with pytest.raises(SystemExit):
        cli.main([cmd, "--help"])
    out = capsys.readouterr().out
    assert "--date" in out and "optional if --start/--end" in out


# --------------------------------------------------------------------------- #
# nexrad --pad-deg reaches the stack (regression: accepted, silently ignored)
# --------------------------------------------------------------------------- #
class _Probe(Exception):
    """Stops _cmd_nexrad right after the call under test records its kwargs."""


def _run_nexrad_intensity(monkeypatch, tmp_path, extra):
    from stormscape import nexrad
    record = {}

    def fake_stack(aoi, start, end, **kw):
        record.update(kw)
        raise _Probe

    monkeypatch.setattr(nexrad, "intensity_stack", fake_stack)
    with pytest.raises(_Probe):
        cli.main(["nexrad", "--intensity",
                  "--bbox", "-120", "39", "-119", "40",
                  "--start", "202606192000", "--end", "202606200200",
                  "--radar", "KRGX", "--out-dir", str(tmp_path)] + extra)
    return record


def test_nexrad_pad_deg_is_forwarded_to_the_stack(monkeypatch, tmp_path):
    kw = _run_nexrad_intensity(monkeypatch, tmp_path, ["--pad-deg", "0.06"])
    assert kw["pad_deg"] == pytest.approx(0.06)


def test_nexrad_default_pad_deg_reaches_the_stack(monkeypatch, tmp_path):
    kw = _run_nexrad_intensity(monkeypatch, tmp_path, [])
    assert kw["pad_deg"] == pytest.approx(0.05)
