"""MRMS stacking: the i15 estimator, date parsing, wet-hour grouping, AOI window."""
from __future__ import annotations

import datetime as dt

import numpy as np
import pytest

from stormscape import mrms


# --------------------------------------------------------------------------- #
# the i15 estimator — the scientific heart of the package
# --------------------------------------------------------------------------- #
def _steps(rate_mmph, n=10, shape=(2, 2)):
    """n two-minute accumulations (mm) for a constant rate (mm/h)."""
    return [np.full(shape, rate_mmph * 2.0 / 60.0) for _ in range(n)]


@pytest.mark.parametrize("rate", [0.0, 12.0, 60.0, 102.3])
def test_compute_i15_recovers_a_constant_rate(rate):
    """i15 = mean(i16, i14); under a constant rate both windows equal that rate."""
    out = mrms.compute_i15(_steps(rate))
    assert np.allclose(out, rate, atol=1e-9)


def test_compute_i15_uses_only_the_trailing_eight_steps():
    """The estimator is a trailing 16-minute window: earlier steps must not leak."""
    stack = _steps(0.0, n=20) + _steps(60.0, n=8)      # dry history, then wet
    assert np.allclose(mrms.compute_i15(stack), 60.0, atol=1e-9)
    # and a wet history followed by dry must read dry
    stack = _steps(60.0, n=20) + _steps(0.0, n=8)
    assert np.allclose(mrms.compute_i15(stack), 0.0, atol=1e-9)


def test_compute_i15_matches_the_hand_computed_formula():
    """Known answer for an uneven stack, from i16 = Σ(8)·60/16 and i14 = Σ(7)·60/14."""
    vals = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]     # mm per 2-min step
    stack = [np.full((1, 1), v) for v in vals]
    i16 = sum(vals) * 60.0 / 16.0
    i14 = sum(vals[1:]) * 60.0 / 14.0
    assert mrms.compute_i15(stack)[0, 0] == pytest.approx((i16 + i14) / 2)


def test_compute_i15_treats_missing_steps_as_zero_not_nan():
    """NaN cells must not poison the sum — a gap is 'no rain measured', not NaN out."""
    stack = _steps(60.0, n=8)
    stack[3] = np.full((2, 2), np.nan)
    out = mrms.compute_i15(stack)
    assert np.isfinite(out).all()
    assert (out < 60.0).all()          # the missing step lowers, not nullifies


def test_compute_i15_is_spatially_independent():
    """Each cell is reduced on its own; no cross-cell mixing."""
    stack = [np.array([[2.0, 0.0]]) for _ in range(8)]
    out = mrms.compute_i15(stack)
    assert out[0, 0] == pytest.approx(60.0)
    assert out[0, 1] == pytest.approx(0.0)


# --------------------------------------------------------------------------- #
# date parsing
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("spec", ["20260619", "2026-06-19", dt.date(2026, 6, 19),
                                  dt.datetime(2026, 6, 19, 13, 30)])
def test_parse_date_accepts_the_documented_forms(spec):
    assert mrms.parse_date(spec) == dt.date(2026, 6, 19)


def test_parse_date_rejects_nonsense():
    with pytest.raises(Exception):
        mrms.parse_date("not-a-date")


# --------------------------------------------------------------------------- #
# wet-hour grouping
# --------------------------------------------------------------------------- #
def test_contiguous_runs_groups_adjacent_hours():
    h = [dt.datetime(2026, 6, 19, x) for x in (20, 21, 22)]
    assert mrms.contiguous_runs(h) == [h]


def test_contiguous_runs_splits_on_a_gap():
    h = [dt.datetime(2026, 6, 19, x) for x in (20, 21, 23, 24 - 1)]
    runs = mrms.contiguous_runs(sorted(h))
    assert len(runs) == 2
    assert [len(r) for r in runs] == [2, 2]


def test_contiguous_runs_single_hour():
    h = [dt.datetime(2026, 6, 19, 22)]
    assert mrms.contiguous_runs(h) == [h]


# --------------------------------------------------------------------------- #
# AOI window on the MRMS CONUS grid
# --------------------------------------------------------------------------- #
def test_aoi_window_covers_the_requested_bounds():
    """The window's transform must place the AOI inside the returned raster."""
    bounds = (-119.75, 39.33, -119.45, 39.88)
    win, tr = mrms.aoi_window(bounds)
    assert win.width > 0 and win.height > 0
    west, north = tr * (0, 0)
    east, south = tr * (win.width, win.height)
    assert west <= bounds[0] and east >= bounds[2]
    assert south <= bounds[1] and north >= bounds[3]


def test_aoi_window_resolution_is_the_mrms_grid():
    _, tr = mrms.aoi_window((-119.75, 39.33, -119.45, 39.88))
    assert tr.a == pytest.approx(mrms.G_RES)
    assert abs(tr.e) == pytest.approx(mrms.G_RES)


def test_aoi_window_clamps_at_the_grid_edge():
    """An AOI straddling the grid's western edge must clamp to it rather than
    producing a negative column offset."""
    win, _ = mrms.aoi_window((mrms.G_W - 1.0, 39.0, mrms.G_W + 1.0, 40.0))
    assert win.col_off == 0
    assert win.width > 0


def test_aoi_window_clamps_at_the_northern_edge():
    win, _ = mrms.aoi_window((-119.7, mrms.G_N - 1.0, -119.4, mrms.G_N + 1.0))
    assert win.row_off == 0
    assert win.height > 0
