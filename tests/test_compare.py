"""Radar-vs-gauge sampling and skill statistics, including the quality and
cadence screens that decide which gauges count."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from stormscape import compare


# --------------------------------------------------------------------------- #
# point sampling
# --------------------------------------------------------------------------- #
def test_sample_raster_at_points_reads_the_containing_cell(field_tif, points_gdf):
    """Cell (row, col) holds value row*10+col; sample the centre of a known cell."""
    arr = np.array([[float(r * 10 + c) for c in range(4)] for r in range(3)])
    tif = field_tif(arr, west=-119.75, north=39.60, res=0.01)
    # centre of row 1, col 2 -> lon = -119.75 + 2.5*0.01, lat = 39.60 - 1.5*0.01
    pts = points_gdf([(-119.75 + 0.025, 39.60 - 0.015)])
    assert compare.sample_raster_at_points(pts, tif)[0] == pytest.approx(12.0)


def test_sample_raster_at_points_off_grid_is_nan(field_tif, points_gdf):
    tif = field_tif(np.ones((3, 3)))
    got = compare.sample_raster_at_points(points_gdf([(-100.0, 20.0)]), tif)
    assert np.isnan(got[0])


def test_sample_raster_at_points_reprojects_non_4326_points(field_tif, points_gdf):
    """Gauges may arrive in any CRS; the sampler must reproject, not mis-index."""
    arr = np.full((5, 5), 42.0)
    tif = field_tif(arr, west=-119.75, north=39.60, res=0.01)
    pts = points_gdf([(-119.73, 39.58)]).to_crs(3857)
    assert compare.sample_raster_at_points(pts, tif)[0] == pytest.approx(42.0)


def test_sampler_parity_with_the_smoothing_module(field_tif, points_gdf):
    """smoothing.gauge_skill_sweep has its own in-memory sampler; it must index
    cells identically (rowcol with op=floor). Using round instead is off by half a
    cell and silently corrupts every sweep result."""
    import rasterio
    from stormscape import smoothing
    arr = np.array([[float(r * 10 + c) for c in range(6)] for r in range(6)])
    tif = field_tif(arr, west=-119.75, north=39.60, res=0.01)
    lonlats = [(-119.75 + 0.01 * (c + 0.5), 39.60 - 0.01 * (r + 0.5))
               for r in (0, 2, 5) for c in (0, 3, 5)]
    pts = points_gdf(lonlats)
    via_compare = compare.sample_raster_at_points(pts, tif)
    with rasterio.open(tif) as ds:
        band = ds.read(1, masked=True).filled(np.nan)
        tr, crs = ds.transform, ds.crs
    via_smoothing = smoothing._sample_array_at_points(band, tr, crs, pts)
    assert np.allclose(via_compare, via_smoothing, equal_nan=True)


# --------------------------------------------------------------------------- #
# skill statistics
# --------------------------------------------------------------------------- #
def _table(gauge, radar, **extra):
    """A radar_vs_gauge-shaped frame for one metric (i15)."""
    g = np.asarray(gauge, float)
    r = np.asarray(radar, float)
    return pd.DataFrame({"i15_mmph": g, "radar_i15max": r, "resid_i15max": r - g,
                         **extra})


def _row(stats, metric="i15max"):
    """One metric's stats as a plain dict.

    Deliberately not a Series: ``corr`` collides with ``Series.corr``, so
    attribute access would silently compare a bound method instead of the value.
    """
    return stats.set_index("metric").loc[metric].to_dict()


def test_perfect_agreement_gives_zero_bias_unit_ratio_unit_corr():
    s = _row(compare.comparison_stats(_table([10, 20, 30], [10, 20, 30]),
                                      pairs={"i15_mmph": "i15max"}))
    assert s["n"] == 3
    assert s["bias"] == pytest.approx(0.0)
    assert s["rmse"] == pytest.approx(0.0)
    assert s["ratio"] == pytest.approx(1.0)
    assert s["corr"] == pytest.approx(1.0)


def test_uniform_over_read_shows_up_as_the_ratio():
    """The headline Hidden Valley finding is a bias *ratio* (~2-4x), so this number
    has to be right."""
    s = _row(compare.comparison_stats(_table([10, 20, 30], [25, 50, 75]),
                                      pairs={"i15_mmph": "i15max"}))
    assert s["ratio"] == pytest.approx(2.5)
    assert s["bias"] == pytest.approx(30.0)
    assert s["corr"] == pytest.approx(1.0)   # perfectly correlated, badly biased


def test_rmse_and_mae_are_computed_correctly():
    s = _row(compare.comparison_stats(_table([10, 10], [12, 6]),
                                      pairs={"i15_mmph": "i15max"}))
    assert s["mae"] == pytest.approx(3.0)                     # |2|, |-4|
    assert s["rmse"] == pytest.approx(np.sqrt((4 + 16) / 2))


def test_nan_pairs_are_dropped_from_n():
    s = _row(compare.comparison_stats(_table([10, np.nan, 30], [10, 20, np.nan]),
                                      pairs={"i15_mmph": "i15max"}))
    assert s["n"] == 1


def test_empty_sample_yields_nan_stats_not_an_error():
    s = _row(compare.comparison_stats(_table([np.nan], [np.nan]),
                                      pairs={"i15_mmph": "i15max"}))
    assert s["n"] == 0 and np.isnan(s["ratio"])


def test_rqi_screen_drops_low_quality_gauges():
    """Beyond radar range or behind terrain the beam overshoots low rain, so RQI
    filtering is required for quantitative work."""
    t = _table([10, 20, 30], [10, 20, 300], rqi=[1.0, 0.95, 0.2])
    assert _row(compare.comparison_stats(t, pairs={"i15_mmph": "i15max"}))["n"] == 3
    screened = _row(compare.comparison_stats(t, pairs={"i15_mmph": "i15max"},
                                             rqi_min=0.8))
    assert screened["n"] == 2
    assert screened["ratio"] == pytest.approx(1.0)   # the bad cell is gone


def test_cadence_screen_applies_to_subhourly_metrics():
    """Coarse reporters have their bursts smeared low by the 1-minute
    interpolation, so they must be excluded from i15/i30/i60."""
    t = _table([10, 20, 5], [10, 20, 60], report_min=[5.0, 15.0, 1440.0])
    screened = _row(compare.comparison_stats(t, pairs={"i15_mmph": "i15max"},
                                             max_report_min=60))
    assert screened["n"] == 2
    assert screened["ratio"] == pytest.approx(1.0)


def test_cadence_screen_does_not_apply_to_the_storm_total():
    """A storm total is cadence-insensitive — a daily gauge still measures it — so
    screening it away would throw out good data."""
    g, r = [10.0, 20.0, 30.0], [10.0, 20.0, 30.0]
    t = pd.DataFrame({"total_mm": g, "radar_total": r,
                      "resid_total": np.subtract(r, g),
                      "report_min": [5.0, 15.0, 1440.0]})
    s = _row(compare.comparison_stats(t, pairs={"total_mm": "total"},
                                      max_report_min=60), metric="total")
    assert s["n"] == 3


def test_missing_metric_columns_are_skipped():
    """Asking for i30/i60 when only i15 was sampled must not fabricate rows."""
    stats = compare.comparison_stats(_table([10], [10]),
                                     pairs={"i15_mmph": "i15max",
                                            "i30_mmph": "i30max"})
    assert set(stats["metric"]) == {"i15max"}


# --------------------------------------------------------------------------- #
# end-to-end radar_vs_gauge against a real raster
# --------------------------------------------------------------------------- #
def test_radar_vs_gauge_builds_radar_and_resid_columns(field_tif, points_gdf):
    tif = field_tif(np.full((6, 6), 50.0), west=-119.75, north=39.60, res=0.01)
    pts = points_gdf([(-119.74, 39.59), (-119.72, 39.57)], i15_mmph=[25.0, 20.0])
    out = compare.radar_vs_gauge(pts, {"i15max": tif},
                                 pairs={"i15_mmph": "i15max"})
    assert out["radar_i15max"].tolist() == pytest.approx([50.0, 50.0])
    assert out["resid_i15max"].tolist() == pytest.approx([25.0, 30.0])
