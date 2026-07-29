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
