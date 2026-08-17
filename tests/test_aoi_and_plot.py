"""AOI parsing and the small pure helpers behind the figures."""
from __future__ import annotations

import numpy as np
import pytest

from stormscape import aoi as aoi_mod
from stormscape import plot


# --------------------------------------------------------------------------- #
# AOI parsing
# --------------------------------------------------------------------------- #
def test_pad_bounds_expands_symmetrically():
    assert aoi_mod.pad_bounds((-119.7, 39.3, -119.4, 39.9), 0.1) == pytest.approx(
        (-119.8, 39.2, -119.3, 40.0))


def test_pad_bounds_zero_is_a_no_op():
    b = (-119.7, 39.3, -119.4, 39.9)
    assert aoi_mod.pad_bounds(b, 0.0) == pytest.approx(b)


def test_bbox_polygon_matches_its_bounds():
    b = (-119.7, 39.3, -119.4, 39.9)
    assert aoi_mod.bbox_polygon(b).bounds == pytest.approx(b)


def test_load_aoi_from_a_bbox_tuple_has_no_polygon():
    """A --bbox AOI carries no clip geometry, which is why clip_to_aoi is a no-op
    for bbox runs."""
    bounds, geom = aoi_mod.load_aoi((-119.7, 39.3, -119.4, 39.9))
    assert bounds == pytest.approx((-119.7, 39.3, -119.4, 39.9))
    assert geom is None


def test_load_aoi_applies_the_pad():
    bounds, _ = aoi_mod.load_aoi((-119.7, 39.3, -119.4, 39.9), pad_deg=0.05)
    assert bounds == pytest.approx((-119.75, 39.25, -119.35, 39.95))


def test_load_aoi_from_a_shapely_geometry_keeps_the_polygon():
    poly = aoi_mod.bbox_polygon((-119.7, 39.3, -119.4, 39.9))
    bounds, geom = aoi_mod.load_aoi(poly)
    assert geom is not None
    assert bounds == pytest.approx(poly.bounds)


def test_load_aoi_from_a_vector_file_roundtrips(tmp_path):
    import geopandas as gpd
    poly = aoi_mod.bbox_polygon((-119.7, 39.3, -119.4, 39.9))
    p = tmp_path / "aoi.geojson"
    gpd.GeoDataFrame(geometry=[poly], crs=4326).to_file(p, driver="GeoJSON")
    bounds, geom = aoi_mod.load_aoi(str(p))
    assert geom is not None
    assert bounds == pytest.approx(poly.bounds, abs=1e-9)


def test_load_aoi_reprojects_a_non_4326_vector(tmp_path):
    """AOIs arrive in whatever CRS the user drew them in; bounds must come back in
    lon/lat degrees regardless."""
    import geopandas as gpd
    poly = aoi_mod.bbox_polygon((-119.7, 39.3, -119.4, 39.9))
    p = tmp_path / "aoi_utm.geojson"
    gpd.GeoDataFrame(geometry=[poly], crs=4326).to_crs(32611).to_file(
        p, driver="GeoJSON")
    bounds, _ = aoi_mod.load_aoi(str(p))
    assert bounds == pytest.approx((-119.7, 39.3, -119.4, 39.9), abs=1e-4)


# --------------------------------------------------------------------------- #
# gauge-tempo label
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("value", [None, float("nan"), "not-a-number"])
def test_tempo_label_is_omitted_when_unmeasurable(value):
    assert plot._tempo_label(value) is None


def test_tempo_label_names_the_hourly_and_daily_cases():
    assert "hourly" in plot._tempo_label(60)
    assert "daily" in plot._tempo_label(1440)


@pytest.mark.parametrize("m", [5, 10, 15])
def test_tempo_label_reports_subhourly_intervals(m):
    assert f"{m}-min" in plot._tempo_label(m)


# --------------------------------------------------------------------------- #
# hillshade NaN fill
# --------------------------------------------------------------------------- #
def test_fill_hillshade_nan_replaces_corners_with_the_mean_tone():
    """Reprojecting 5070 -> UTM rotates the raster, leaving NaN corners. Filling
    them with the mean terrain tone makes them blend in; leaving them NaN gives
    either stark white wedges or matplotlib's fragile layout-dependent gray box."""
    hsv = np.full((5, 5), 120.0)
    hsv[0, 0] = np.nan
    out = plot._fill_hillshade_nan(hsv)
    assert np.isfinite(out).all()
    assert out[0, 0] == pytest.approx(120.0)


def test_fill_hillshade_nan_leaves_a_complete_array_untouched():
    hsv = np.arange(9, dtype=float).reshape(3, 3)
    assert np.allclose(plot._fill_hillshade_nan(hsv), hsv)


def test_fill_hillshade_nan_passes_an_all_nan_array_through():
    """Nothing to average against — must not raise."""
    out = plot._fill_hillshade_nan(np.full((3, 3), np.nan))
    assert np.isnan(out).all()


def test_default_field_alpha_is_the_documented_project_value():
    """One project-wide drape opacity; every map resolves None to it."""
    assert plot.DEFAULT_FIELD_ALPHA == pytest.approx(0.32)


# --------------------------------------------------------------------------- #
# label placement
# --------------------------------------------------------------------------- #
@pytest.fixture()
def label_ax():
    """A 6x6 inch axis spanning 0-100 in both directions, 100 dpi -> 600 px."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 6), dpi=100)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    fig.canvas.draw()
    yield ax
    plt.close(fig)


def test_shorten_cuts_at_a_word_boundary():
    assert plot.shorten("North Fork Little Humboldt River", 26) == \
        "North Fork Little…"
    assert plot.shorten("Sheep Creek", 26) == "Sheep Creek"


def test_shorten_falls_back_to_a_hard_cut_for_one_long_word():
    """No boundary to cut at — better a truncated word than an overflowing one."""
    assert plot.shorten("Supercalifragilistic", 10) == "Supercali…"


def test_interior_point_prefers_a_vertex_away_from_the_frame():
    from shapely.geometry import LineString

    # runs from the very edge of the window into its middle
    line = LineString([(0.5, 50), (20, 50), (50, 50)])
    x, y = plot.interior_point(line, (0, 0, 100, 100))
    assert (x, y) == (50, 50)


def test_interior_point_ignores_vertices_outside_the_window():
    from shapely.geometry import LineString

    line = LineString([(-40, 50), (60, 50)])
    assert plot.interior_point(line, (0, 0, 100, 100)) == (60, 50)


def test_interior_point_of_a_line_wholly_outside_is_none():
    from shapely.geometry import LineString

    assert plot.interior_point(LineString([(-9, -9), (-5, -5)]),
                               (0, 0, 100, 100)) is None


def test_labels_do_not_land_on_each_other(label_ax):
    lab = plot.Labeller(label_ax)
    assert lab.label(50, 50, "Alpha", fontsize=8)
    assert lab.label(50, 50, "Beta", fontsize=8)     # same point, must move

    (a, b) = lab.taken
    overlap = not (a[2] < b[0] or a[0] > b[2] or a[3] < b[1] or a[1] > b[3])
    assert not overlap


def test_a_label_never_leaves_the_panel(label_ax):
    """The frame test replaces clipping: a label that would run off is moved,
    not cut in half at the edge."""
    lab = plot.Labeller(label_ax)
    assert lab.label(99.5, 99.5, "Corner", fontsize=8)

    ab = label_ax.get_window_extent()
    (x0, y0, x1, y1) = lab.taken[0]
    assert x0 >= ab.x0 and x1 <= ab.x1
    assert y0 >= ab.y0 and y1 <= ab.y1


def test_label_reports_failure_rather_than_stacking(label_ax):
    """A label on top of another is worse than a missing one, so the caller is
    told instead."""
    lab = plot.Labeller(label_ax)
    # wall off every candidate ring around the point
    for dx in range(-90, 91, 6):
        for dy in range(-90, 91, 6):
            lab.taken.append((300 + dx, 300 + dy, 306 + dx, 306 + dy))
    assert lab.label(50, 50, "Nowhere", fontsize=8) is False


def test_a_displaced_label_gets_a_leader_line(label_ax):
    """Far enough to be ambiguous means far enough to need connecting."""
    lab = plot.Labeller(label_ax)
    lab.label(50, 50, "First", fontsize=8)
    before = len(label_ax.lines)
    lab.label(50, 50, "Second", fontsize=8, force_leader=True)

    assert len(label_ax.lines) == before + 1


def test_an_adjacent_label_draws_no_leader(label_ax):
    """A 2 px stub beside its own dot is noise, not information."""
    lab = plot.Labeller(label_ax)
    lab.label(50, 50, "Close", fontsize=8)
    assert not label_ax.lines


def test_forced_leaders_clear_the_label_box(label_ax):
    """The reason annotate's arrows looked absent: an offset shorter than the
    text's own half-width leaves no visible line."""
    lab = plot.Labeller(label_ax)
    lab.label(50, 50, "12", fontsize=8.5, force_leader=True)
    (x0, y0, x1, y1) = lab.taken[0]
    cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    px, py = label_ax.transData.transform((50, 50))

    assert np.hypot(cx - px, cy - py) >= plot.LEADER_FROM
    assert np.hypot(cx - px, cy - py) > (x1 - x0) / 2


def test_crowded_sees_a_neighbouring_marker(label_ax):
    lab = plot.Labeller(label_ax)
    lab.block_many([50, 51], [50, 50])       # two dots ~6 px apart
    assert lab.crowded(50, 50)
    lab2 = plot.Labeller(label_ax)
    lab2.block_many([10, 90], [10, 90])      # far apart
    assert not lab2.crowded(10, 10)


def test_labels_avoid_registered_markers(label_ax):
    lab = plot.Labeller(label_ax)
    lab.block(50, 50, radius_px=20)
    assert lab.label(50, 50, "X", fontsize=8)
    (x0, y0, x1, y1) = lab.taken[-1]
    px, py = label_ax.transData.transform((50, 50))
    blocked = (px - 20, py - 20, px + 20, py + 20)
    assert (x1 < blocked[0] or x0 > blocked[2]
            or y1 < blocked[1] or y0 > blocked[3])
