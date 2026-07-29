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
