"""Radar-gauge bias correction and merging."""
from __future__ import annotations

import numpy as np
import pytest
from rasterio.transform import from_origin

from stormscape import merge

TR = from_origin(-119.75, 39.60, 0.01, 0.01)


def _gauge_points(n=5):
    """n gauge locations inside a 20x20 cell grid rooted at TR."""
    lons = [-119.75 + 0.01 * (c + 0.5) for c in range(2, 2 + n)]
    lats = [39.60 - 0.01 * (r + 0.5) for r in range(2, 2 + n)]
    return np.array(lons), np.array(lats)


def test_sample_field_reads_the_containing_cell():
    field = np.array([[float(r * 10 + c) for c in range(5)] for r in range(5)])
    lon = -119.75 + 0.01 * 2.5
    lat = 39.60 - 0.01 * 1.5
    assert merge.sample_field(field, TR, [lon], [lat])[0] == pytest.approx(12.0)


def test_sample_field_off_grid_is_nan():
    field = np.ones((4, 4))
    assert np.isnan(merge.sample_field(field, TR, [-100.0], [20.0])[0])


def test_mean_field_bias_factor_is_gauge_over_radar():
    """Radar reads 4x the gauges, so the correction factor must be 0.25 and the
    corrected field must match the gauges."""
    field = np.full((20, 20), 40.0)
    lons, lats = _gauge_points()
    vals = np.full(len(lons), 10.0)
    corrected, factor = merge.mean_field_bias(field, TR, lons, lats, vals)
    assert factor == pytest.approx(0.25)
    assert np.allclose(corrected, 10.0)


def test_mean_field_bias_is_identity_when_unbiased():
    field = np.full((20, 20), 10.0)
    lons, lats = _gauge_points()
    corrected, factor = merge.mean_field_bias(field, TR, lons, lats,
                                             np.full(len(lons), 10.0))
    assert factor == pytest.approx(1.0)
    assert np.allclose(corrected, field)


def test_mean_field_bias_never_returns_negatives():
    field = np.full((10, 10), 5.0)
    lons, lats = _gauge_points(3)
    corrected, _ = merge.mean_field_bias(field, TR, lons, lats, np.zeros(3))
    assert (corrected >= 0).all()


def test_mean_field_bias_falls_back_to_unity_without_usable_gauges():
    field = np.full((10, 10), 5.0)
    lons, lats = _gauge_points(3)
    _, factor = merge.mean_field_bias(field, TR, lons, lats,
                                      np.full(3, np.nan))
    assert factor == pytest.approx(1.0)


def test_local_bias_preserves_shape_and_finiteness():
    rng = np.random.RandomState(0)
    field = rng.rand(20, 20) * 50 + 1
    lons, lats = _gauge_points(6)
    out = merge.local_bias(field, TR, lons, lats,
                           merge.sample_field(field, TR, lons, lats) * 0.5)
    assert out.shape == field.shape
    assert np.isfinite(out).all()
    assert (out >= 0).all()


def test_local_bias_pulls_the_field_toward_the_gauges():
    """With every gauge reading half the radar, the corrected field should drop."""
    field = np.full((20, 20), 40.0)
    lons, lats = _gauge_points(6)
    out = merge.local_bias(field, TR, lons, lats, np.full(len(lons), 20.0))
    assert np.nanmean(out) < np.nanmean(field)
    assert np.nanmean(out) == pytest.approx(20.0, rel=0.15)


def test_local_bias_falls_back_to_mean_field_bias_with_too_few_gauges():
    """Fewer than 3 wet gauges cannot support a spatial correction."""
    field = np.full((10, 10), 40.0)
    lons, lats = _gauge_points(2)
    out = merge.local_bias(field, TR, lons, lats, np.full(2, 10.0))
    assert np.allclose(out, 10.0)


def test_local_bias_keeps_nan_cells_nan():
    field = np.full((12, 12), 20.0)
    field[0, 0] = np.nan
    lons, lats = _gauge_points(5)
    out = merge.local_bias(field, TR, lons, lats, np.full(5, 10.0))
    assert np.isnan(out[0, 0])


def test_conditional_merge_returns_a_finite_field_of_the_same_shape():
    rng = np.random.RandomState(1)
    field = rng.rand(18, 18) * 30 + 1
    lons, lats = _gauge_points(6)
    vals = merge.sample_field(field, TR, lons, lats) * 1.5
    out = merge.conditional_merge(field, TR, lons, lats, vals)
    assert out.shape == field.shape
    assert np.isfinite(out).sum() > 0
    assert np.nanmin(out) >= 0


def test_loo_cross_validate_predicts_every_held_out_gauge():
    """Leave-one-out is how we check a merge generalizes rather than memorizes: one
    held-out prediction per gauge, per method (raw / mean-field / local / conditional)."""
    rng = np.random.RandomState(2)
    field = rng.rand(20, 20) * 30 + 1
    lons, lats = _gauge_points(6)
    vals = merge.sample_field(field, TR, lons, lats) * 0.8
    table, stats = merge.loo_cross_validate(field, TR, lons, lats, vals)
    assert len(table) == len(lons)
    for col in ("gauge", "raw", "mfb", "lbc", "cm"):
        assert col in table.columns
    assert table["gauge"].to_numpy() == pytest.approx(vals)
    assert stats is not None


def test_loo_cross_validate_corrections_beat_the_raw_field_on_bias():
    """With the radar uniformly over-reading, every correction should land closer to
    the gauges than the raw field does — the whole point of merging."""
    field = np.full((20, 20), 40.0)
    lons, lats = _gauge_points(6)
    vals = np.full(len(lons), 10.0)
    table, _ = merge.loo_cross_validate(field, TR, lons, lats, vals)
    raw_err = float(np.abs(table["raw"] - table["gauge"]).mean())
    for col in ("mfb", "lbc", "cm"):
        assert float(np.abs(table[col] - table["gauge"]).mean()) < raw_err
