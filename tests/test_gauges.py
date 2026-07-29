"""Rain-gauge reduction: cadence measurement, intensity metrics, storm window,
and the filename-collision fix. All offline — no Synoptic token needed."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from stormscape import gauges


# --------------------------------------------------------------------------- #
# native reporting cadence — measured, never taken from metadata
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("step", [5, 10, 15, 60])
def test_report_min_measures_a_uniform_cadence(step, gauge_obs):
    obs = gauge_obs([1.0] * 10, step_min=step)
    assert gauges._report_min(obs) == pytest.approx(step)


def test_report_min_ignores_rows_without_precip(gauge_obs):
    """The multi-variable ASOS trap: a station may report temperature every 5
    minutes but precipitation only at hourly METAR times. Counting *all* rows
    would label it 5-min and wrongly pass a `--max-report-min` screen, so the
    cadence must come from precip-bearing rows only."""
    inc = [None] * 25
    for i in (0, 12, 24):                     # precip only every 60 min
        inc[i] = 2.0
    obs = gauge_obs(inc, step_min=5)
    assert gauges._report_min(obs) == pytest.approx(60.0)


def test_report_min_is_nan_without_a_precip_column():
    obs = pd.DataFrame({"date_time": pd.date_range("2026-06-19", periods=5, freq="5min"),
                        "air_temp_set_1": [10.0] * 5})
    assert np.isnan(gauges._report_min(obs))


def test_report_min_is_nan_with_a_single_observation(gauge_obs):
    """One point defines no interval."""
    assert np.isnan(gauges._report_min(gauge_obs([1.0])))


# --------------------------------------------------------------------------- #
# per-gauge metrics
# --------------------------------------------------------------------------- #
def test_gauge_intensities_totals_the_increments(gauge_obs):
    obs = gauge_obs([1.0, 2.0, 3.0, 4.0], step_min=5)
    assert gauges.gauge_intensities(obs)["total_mm"] == pytest.approx(10.0)


def test_gauge_intensities_recovers_a_constant_rate(gauge_obs):
    """5 mm per 5 min = 60 mm/h; after 1-minute interpolation every peak-intensity
    window should read ~60, and the gauge-side i15 uses the same (i16+i14)/2
    estimator as the radar so the two are directly comparable."""
    obs = gauge_obs([5.0] * 24, step_min=5)          # 2 hours of steady rain
    m = gauges.gauge_intensities(obs)
    assert m["total_mm"] == pytest.approx(120.0)
    for d in (15, 30, 60):
        assert m[f"i{d}_mmph"] == pytest.approx(60.0, rel=0.05)


def test_gauge_intensities_orders_the_peak_durations(gauge_obs):
    """A short burst must satisfy i15 >= i30 >= i60 — the same cell-wise ordering
    the radar fields obey."""
    inc = [0.0] * 6 + [6.0] * 4 + [0.0] * 20         # a 20-minute burst
    m = gauges.gauge_intensities(gauge_obs(inc, step_min=5))
    assert m["i15_mmph"] >= m["i30_mmph"] >= m["i60_mmph"] > 0


def test_gauge_intensities_all_nan_without_a_precip_variable():
    """Stations reporting no precip must yield NaN metrics, not crash."""
    obs = pd.DataFrame({"date_time": pd.date_range("2026-06-19", periods=5, freq="5min"),
                        "wind_speed_set_1": [3.0] * 5})
    m = gauges.gauge_intensities(obs)
    assert np.isnan(m["total_mm"])
    assert all(np.isnan(m[f"i{d}_mmph"]) for d in (15, 30, 60))


def test_gauge_intensities_dry_gauge_is_zero_not_nan(gauge_obs):
    """A gauge that reported and measured nothing is dry (0), which is data —
    distinct from a gauge with no precip sensor (NaN)."""
    m = gauges.gauge_intensities(gauge_obs([0.0] * 20, step_min=5))
    assert m["total_mm"] == pytest.approx(0.0)


def test_negative_increments_do_not_reduce_the_total(gauge_obs):
    """Gauge resets/spurious negatives are clipped, not subtracted."""
    m = gauges.gauge_intensities(gauge_obs([2.0, -5.0, 2.0], step_min=5))
    assert m["total_mm"] == pytest.approx(4.0)


# --------------------------------------------------------------------------- #
# 1-minute interpolation
# --------------------------------------------------------------------------- #
def test_precipitation_per_minute_conserves_the_depth_it_can_span(gauge_obs):
    """The interpolated series distributes each increment back over the interval
    that produced it, so it conserves every increment *except the first* — the
    first reading has no preceding interval to spread across. The reported storm
    total still counts it, which is why `total_mm` exceeds the series sum."""
    obs = gauge_obs([5.0] * 12, step_min=5)
    _, _, per_min = gauges.precipitation_per_minute(obs)
    assert float(np.nansum(per_min)) == pytest.approx(55.0, rel=1e-6)
    assert gauges.gauge_intensities(obs)["total_mm"] == pytest.approx(60.0)


def test_precipitation_per_minute_spreads_a_burst_evenly(gauge_obs):
    """A single 5 mm tip over a 5-minute interval becomes 1 mm/min, not a spike in
    one minute — this is what smears a coarse reporter's peak intensity."""
    obs = gauge_obs([0.0, 5.0], step_min=5)
    _, _, per_min = gauges.precipitation_per_minute(obs)
    nz = np.asarray(per_min)[np.asarray(per_min) > 0]
    assert len(nz) == 5
    assert np.allclose(nz, 1.0)


def test_precipitation_per_minute_expands_to_minute_resolution(gauge_obs):
    obs = gauge_obs([5.0] * 12, step_min=5)
    _, _, per_min = gauges.precipitation_per_minute(obs)
    assert len(per_min) >= 55            # ~60 minutes of 1-minute bins


# --------------------------------------------------------------------------- #
# storm window — mass-weighted, so stray tips cannot stretch it
# --------------------------------------------------------------------------- #
def test_storm_window_brackets_the_rain(minute_series):
    rates = np.zeros(600)
    rates[200:260] = 60.0
    win = gauges.storm_window({"g": minute_series(rates)}, pad_min=30)
    assert win is not None
    start, end = win
    idx = minute_series(rates).index
    assert start <= idx[200] and end >= idx[259]


def test_storm_window_excludes_a_distant_stray_tip(minute_series):
    """Mass coverage, not any-nonzero: one late 0.001 mm tip must not stretch the
    window by hours (that bug turned a 5-hour storm into a 27-hour window)."""
    rates = np.zeros(600)
    rates[200:260] = 60.0                 # 60 mm of real storm
    rates[590] = 0.06                     # a single stray 0.001 mm tip
    _, end = gauges.storm_window({"g": minute_series(rates)}, pad_min=30)
    idx = minute_series(rates).index
    assert end < idx[400], "the stray tip stretched the window"


def test_storm_window_is_none_when_dry(minute_series):
    assert gauges.storm_window({"g": minute_series(np.zeros(120))}) is None


def test_storm_window_combines_multiple_gauges(minute_series):
    """The window is taken over the summed mass of all gauges, so it must reach the
    second gauge's burst even though the first gauge is dry by then. (The exact
    edges sit a minute inside each burst because the central 99% of mass trims
    0.5% from each tail.)"""
    a, b = np.zeros(400), np.zeros(400)
    a[100:120] = 30.0
    b[300:320] = 30.0
    start, end = gauges.storm_window({"a": minute_series(a), "b": minute_series(b)},
                                     pad_min=0)
    idx = minute_series(a).index
    assert start <= idx[101]
    assert end >= idx[300], "the window never reached the second gauge"


def test_storm_window_pad_widens_the_window(minute_series):
    rates = np.zeros(300)
    rates[100:140] = 40.0
    narrow = gauges.storm_window({"g": minute_series(rates)}, pad_min=0)
    wide = gauges.storm_window({"g": minute_series(rates)}, pad_min=45)
    assert wide[0] < narrow[0] and wide[1] > narrow[1]


# --------------------------------------------------------------------------- #
# case-insensitive filename collisions (macOS APFS/HFS+)
# --------------------------------------------------------------------------- #
def test_unique_safe_names_separates_case_only_duplicates():
    """'Virginia City' and 'VIRGINIA CITY' are two real, distinct stations whose
    sanitized stems collide on a case-insensitive filesystem — one gauge's CSV and
    figure silently overwrote the other's."""
    m = gauges.unique_safe_names(["Virginia City", "VIRGINIA CITY"])
    assert len(set(v.lower() for v in m.values())) == 2


def test_unique_safe_names_is_order_independent():
    """Reader and writer must derive the same mapping regardless of input order,
    or a store written in one pass cannot be read back in another."""
    names = ["Virginia City", "VIRGINIA CITY", "Six Mile Canyon"]
    assert gauges.unique_safe_names(names) == gauges.unique_safe_names(names[::-1])


def test_unique_safe_names_leaves_non_colliding_names_unsuffixed():
    """Old stores must keep resolving, so only colliding stems may change."""
    m = gauges.unique_safe_names(["Six Mile Canyon", "Hidden Valley"])
    assert not any(v.endswith("_2") for v in m.values())


def test_unique_safe_names_covers_every_input():
    names = ["A b", "a B", "a b", "Unique One"]
    m = gauges.unique_safe_names(names)
    assert set(m) == set(names)
    assert len(set(v.lower() for v in m.values())) == 4
