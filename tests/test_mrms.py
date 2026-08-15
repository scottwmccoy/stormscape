"""MRMS stacking: the i15 estimator, date parsing, wet-hour grouping, AOI window."""
from __future__ import annotations

import datetime as dt
import warnings

import numpy as np
import pandas as pd
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


# --------------------------------------------------------------------------- #
# wet_window — the cheap radar probe that bounds an expensive NEXRAD fetch
# --------------------------------------------------------------------------- #
def _wet_window_over(monkeypatch, wet_hours, start, end, **kw):
    """Run wet_window with fetch_many faked: `wet_hours` return a rainy grid."""
    monkeypatch.setattr(mrms, "load_aoi", lambda *a, **k: ((-119.8, 39.4, -119.5, 39.7), None))
    def fake_fetch_many(product, times, win, workers=None):
        return {t: (np.full((2, 2), 20.0) if t in wet_hours else np.zeros((2, 2)))
                for t in times}
    monkeypatch.setattr(mrms, "fetch_many", fake_fetch_many)
    return mrms.wet_window("aoi", start, end, **kw)


def test_wet_window_brackets_only_the_wet_hours(monkeypatch):
    """A 30 h storm-DAY window holding a 3 h storm must come back ~hours, not ~days."""
    start, end = dt.datetime(2026, 8, 12, 4), dt.datetime(2026, 8, 13, 10)
    wet = [dt.datetime(2026, 8, 12, h) for h in (21, 22, 23)]
    got = _wet_window_over(monkeypatch, wet, start, end, pad_min=30)
    assert got is not None
    lo, hi = got
    # hourly QPE at HH covers the hour ENDING at HH, so the window opens before 21Z
    assert lo == dt.datetime(2026, 8, 12, 19, 30)
    assert hi == dt.datetime(2026, 8, 12, 23, 30)
    assert (hi - lo).total_seconds() / 3600 < 5           # not the full 30 h


def test_wet_window_is_none_when_every_hour_is_dry(monkeypatch):
    start, end = dt.datetime(2026, 8, 12, 4), dt.datetime(2026, 8, 13, 10)
    assert _wet_window_over(monkeypatch, [], start, end) is None


def test_wet_window_spans_two_separate_cells(monkeypatch):
    """Two cells split by a dry gap yield ONE window covering both, not the first."""
    start, end = dt.datetime(2026, 8, 13, 4), dt.datetime(2026, 8, 14, 10)
    wet = [dt.datetime(2026, 8, 13, 20), dt.datetime(2026, 8, 14, 2)]
    lo, hi = _wet_window_over(monkeypatch, wet, start, end, pad_min=0)
    assert lo == dt.datetime(2026, 8, 13, 19)
    assert hi == dt.datetime(2026, 8, 14, 2)


def test_wet_window_never_escapes_the_requested_span(monkeypatch):
    """Padding must clamp to [start, end] so the caller's bounds are respected."""
    start, end = dt.datetime(2026, 8, 12, 20), dt.datetime(2026, 8, 12, 23)
    wet = [dt.datetime(2026, 8, 12, h) for h in (20, 21, 22, 23)]
    lo, hi = _wet_window_over(monkeypatch, wet, start, end, pad_min=600)
    assert lo == start and hi == end


def test_wet_window_honours_the_qpe_threshold(monkeypatch):
    """An hour at/below qpe_thresh is dry — the threshold is what keeps drizzle
    from stretching the window back to a full day."""
    start, end = dt.datetime(2026, 8, 12, 4), dt.datetime(2026, 8, 12, 8)
    wet = [dt.datetime(2026, 8, 12, 6)]
    assert _wet_window_over(monkeypatch, wet, start, end, qpe_thresh=50.0) is None
    assert _wet_window_over(monkeypatch, wet, start, end, qpe_thresh=5.0) is not None


# --------------------------------------------------------------------------- #
# window_hours — an explicit analysis window, or the storm-day span
# --------------------------------------------------------------------------- #
def test_storm_day_span_covers_the_local_day():
    hours = mrms.window_hours(dt.date(2026, 8, 14))
    assert len(hours) == 30                       # [04Z, next-day 10Z]
    assert hours[0] == dt.datetime(2026, 8, 14, 4)
    assert hours[-1] == dt.datetime(2026, 8, 15, 9)


def test_explicit_window_is_used_verbatim():
    w = (dt.datetime(2026, 8, 14, 20), dt.datetime(2026, 8, 15, 4))
    hours = mrms.window_hours(window=w)
    assert hours[0] == w[0] and hours[-1] == w[1]
    assert len(hours) == 9


def test_explicit_window_overrides_the_date():
    """Back-to-back evening storms share a storm-day span; the window is how a
    run says "only tonight's cells, not last night's tail"."""
    w = (dt.datetime(2026, 8, 14, 20), dt.datetime(2026, 8, 15, 4))
    assert mrms.window_hours(dt.date(2026, 1, 1), window=w)[0] == w[0]


def test_window_start_is_floored_to_the_hour():
    hours = mrms.window_hours(window=(dt.datetime(2026, 8, 14, 20, 37),
                                      dt.datetime(2026, 8, 14, 23)))
    assert hours[0] == dt.datetime(2026, 8, 14, 20)


def test_a_date_or_a_window_is_required():
    with pytest.raises(ValueError):
        mrms.window_hours()


# --------------------------------------------------------------------------- #
# the stacked span — hourly QPE is stamped at the END of the hour it describes
# --------------------------------------------------------------------------- #
def _stacked_span(monkeypatch, wet_stamps):
    """The (t0, t1) i15_storm_day actually fetches PrecipRate over."""
    seen = {}
    monkeypatch.setattr(mrms, "load_aoi",
                        lambda *a, **k: ((-119.8, 39.4, -119.5, 39.7), None))
    monkeypatch.setattr(mrms, "find_wet_hours",
                        lambda *a, **k: (pd.DataFrame({"t": wet_stamps}),
                                         pd.DataFrame({"t": wet_stamps,
                                                       "qmax": [9.0] * len(wet_stamps)})))

    def fake_fetch_many(product, times, win, workers=None):
        times = list(times)
        seen.setdefault(product, []).append((times[0], times[-1]))
        return {}                       # empty -> stack resets, loop still runs

    monkeypatch.setattr(mrms, "fetch_many", fake_fetch_many)
    # RQI/SHSR go through `fetch` (singular), not fetch_many — without this the
    # test makes real S3 requests and CI (-m "not network") would break.
    monkeypatch.setattr(mrms, "fetch",
                        lambda product, t, win: np.ones((2, 2), np.float32))
    mrms.i15_storm_day("aoi", dt.date(2026, 8, 13), verbose=False)
    return seen["PrecipRate"][0]


def test_stack_covers_the_hour_each_wet_stamp_describes(monkeypatch):
    """QPE(HH) is the rain in [HH-1, HH] — verified against the 2-min
    accumulation (ratio 1.000 vs 0.13-0.33 for [HH, HH+1]). So wet stamps
    [20Z..23Z] must be stacked over [19Z, 23Z], plus the 14-min i15 lead."""
    wet = [dt.datetime(2026, 8, 13, h) for h in (20, 21, 22, 23)]
    t0, t1 = _stacked_span(monkeypatch, wet)
    assert t0 == dt.datetime(2026, 8, 13, 18, 46)     # 20Z - 1h - 14min
    assert t1 == dt.datetime(2026, 8, 13, 23)         # last wet stamp


def test_stack_does_not_run_past_the_last_wet_hour(monkeypatch):
    """The old span ran to hn+1h, spending a fetch on a usually-dry hour while
    missing the front of the storm."""
    wet = [dt.datetime(2026, 8, 13, 20)]
    t0, t1 = _stacked_span(monkeypatch, wet)
    assert t1 == dt.datetime(2026, 8, 13, 20)
    assert (t1 - t0) == dt.timedelta(hours=1, minutes=14)


def test_a_single_wet_hour_still_gets_the_i15_lead_in(monkeypatch):
    """8 two-minute steps are needed before i15 is defined; the lead must be
    at least that (14 min = 7 steps before the hour's first minute)."""
    wet = [dt.datetime(2026, 8, 13, 20)]
    t0, _ = _stacked_span(monkeypatch, wet)
    assert t0 <= dt.datetime(2026, 8, 13, 19) - dt.timedelta(minutes=14)


# --------------------------------------------------------------------------- #
# the max_wet_hours cap — it ranks by INTENSITY, so it drops the storm's tails
# --------------------------------------------------------------------------- #
def _find_wet_over(monkeypatch, qmax_by_hour, window, **kw):
    """Run find_wet_hours with fetch_many faked from {hour: areal-max mm}."""
    def fake_fetch_many(product, times, win, workers=None):
        return {t: np.full((2, 2), qmax_by_hour.get(t, 0.0)) for t in times}
    monkeypatch.setattr(mrms, "fetch_many", fake_fetch_many)
    return mrms.find_wet_hours(None, object(), window=window, **kw)


# the real 13 Aug 2026 Bug/Stallion storm: 9 wet hours, one over the default cap
_AUG13 = {dt.datetime(2026, 8, 13, 20): 21.7,
          dt.datetime(2026, 8, 13, 21): 38.7,
          dt.datetime(2026, 8, 13, 22): 38.0,
          dt.datetime(2026, 8, 13, 23): 13.7,
          dt.datetime(2026, 8, 14, 0): 10.0,
          dt.datetime(2026, 8, 14, 1): 25.2,
          dt.datetime(2026, 8, 14, 2): 17.2,
          dt.datetime(2026, 8, 14, 3): 8.4,
          dt.datetime(2026, 8, 14, 4): 6.0}
_AUG13_WINDOW = (dt.datetime(2026, 8, 13, 18), dt.datetime(2026, 8, 14, 5))


def test_find_wet_hours_warns_when_the_cap_truncates(monkeypatch):
    """9 wet hours against a cap of 8 must say so — this went unnoticed once."""
    with pytest.warns(UserWarning, match="9 wet hours found"):
        wet, _ = _find_wet_over(monkeypatch, _AUG13, _AUG13_WINDOW,
                                max_wet_hours=8)
    assert len(wet) == 8


def test_the_truncation_warning_names_the_dropped_hour(monkeypatch):
    """Naming the hour is the point: it tells you which end of the storm was lost."""
    with pytest.warns(UserWarning, match=r"08-14 04Z") as rec:
        _find_wet_over(monkeypatch, _AUG13, _AUG13_WINDOW, max_wet_hours=8)
    msg = str(rec[0].message)
    assert "6.0 mm" in msg                    # strongest dropped
    assert "max_wet_hours" in msg             # and how to fix it


def test_the_cap_drops_the_weakest_hour_not_the_last(monkeypatch):
    """It ranks by intensity: 04Z (6.0 mm) goes, even though 03Z (8.4) is later."""
    with pytest.warns(UserWarning):
        wet, _ = _find_wet_over(monkeypatch, _AUG13, _AUG13_WINDOW,
                                max_wet_hours=8)
    kept = list(wet.t)
    assert dt.datetime(2026, 8, 14, 4) not in kept
    assert dt.datetime(2026, 8, 14, 3) in kept


def test_truncation_shortens_the_stacked_span(monkeypatch):
    """The real damage: the run ends an hour early, so `total` loses that rain."""
    with pytest.warns(UserWarning):
        wet8, _ = _find_wet_over(monkeypatch, _AUG13, _AUG13_WINDOW,
                                 max_wet_hours=8)
    wet12, _ = _find_wet_over(monkeypatch, _AUG13, _AUG13_WINDOW,
                              max_wet_hours=12)
    assert mrms.contiguous_runs(list(wet8.t))[0][-1] == dt.datetime(2026, 8, 14, 3)
    assert mrms.contiguous_runs(list(wet12.t))[0][-1] == dt.datetime(2026, 8, 14, 4)


def test_find_wet_hours_is_quiet_when_the_cap_is_not_reached(monkeypatch):
    """A 4-wet-hour storm (14 Aug) under the default cap must not warn."""
    q = {dt.datetime(2026, 8, 14, 23): 11.9,
         dt.datetime(2026, 8, 15, 0): 32.9,
         dt.datetime(2026, 8, 15, 1): 34.8,
         dt.datetime(2026, 8, 15, 2): 3.4}
    with warnings.catch_warnings():
        warnings.simplefilter("error")        # any warning fails the test
        wet, _ = _find_wet_over(monkeypatch, q,
                                (dt.datetime(2026, 8, 14, 20),
                                 dt.datetime(2026, 8, 15, 3)))
    assert len(wet) == 4


def test_the_all_dry_fallback_never_warns(monkeypatch):
    """A dry span falls back to the single best hour — one row, never truncation."""
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        wet, _ = _find_wet_over(monkeypatch, {},
                                (dt.datetime(2026, 8, 11, 0),
                                 dt.datetime(2026, 8, 11, 6)),
                                max_wet_hours=1)
    assert len(wet) == 1
