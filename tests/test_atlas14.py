"""NOAA Atlas 14: recurrence-interval inversion, region lookup, grid URLs."""
from __future__ import annotations

import numpy as np
import pytest

from stormscape import atlas14

# a representative Atlas 14 partial-duration 15-min intensity curve (mm/h)
ARIS = [1, 2, 5, 10, 25, 50, 100, 200, 500, 1000]
CURVE = [18.0, 24.0, 32.0, 40.0, 52.0, 62.0, 73.0, 85.0, 103.0, 118.0]


@pytest.mark.parametrize("i, expected", list(zip(CURVE, ARIS)))
def test_tabulated_intensities_invert_to_their_own_ari(i, expected):
    """Every quantile on the curve must map back to its own recurrence interval."""
    assert atlas14.recurrence_interval(i, ARIS, CURVE) == pytest.approx(expected, rel=1e-6)


def test_interpolation_between_tabulated_points_is_bracketed():
    ri = atlas14.recurrence_interval(36.0, ARIS, CURVE)      # between 5 and 10 yr
    assert 5.0 < ri < 10.0


def test_sub_annual_intensities_return_a_number_not_nan():
    """Partial-duration series admits >1 event/year, so below the 1-yr quantile the
    RI is a real sub-annual number. It previously returned NaN and the CLI printed
    a bare '<1', losing information."""
    for obs in (15.0, 10.0, 5.0):
        ri = atlas14.recurrence_interval(obs, ARIS, CURVE)
        assert np.isfinite(ri), f"{obs} mm/h gave a non-finite RI"
        assert 0.0 < ri < 1.0


def test_sub_annual_branch_is_continuous_across_one_year():
    """No discontinuity at the 1-yr knot where the extrapolation takes over."""
    just_below = atlas14.recurrence_interval(17.99, ARIS, CURVE)
    at_one = atlas14.recurrence_interval(18.0, ARIS, CURVE)
    assert just_below < at_one
    assert at_one - just_below < 0.02


def test_recurrence_is_monotone_increasing_in_intensity():
    ris = [atlas14.recurrence_interval(x, ARIS, CURVE)
           for x in (5, 15, 18, 24, 40, 73, 118)]
    assert all(b > a for a, b in zip(ris, ris[1:]))


def test_above_the_top_ari_is_infinite_not_extrapolated():
    """Beyond the published curve we decline to invent a number."""
    assert atlas14.recurrence_interval(500.0, ARIS, CURVE) == float("inf")


@pytest.mark.parametrize("bad", [0.0, -3.0, float("nan")])
def test_non_positive_or_missing_intensity_is_nan(bad):
    assert np.isnan(atlas14.recurrence_interval(bad, ARIS, CURVE))


# --------------------------------------------------------------------------- #
# region lookup (bundled table, no network)
# --------------------------------------------------------------------------- #
def test_nevada_maps_to_the_southwest_volume():
    """Volume 1 (Southwest) covers NV/CA/AZ/UT/NM — the Hidden Valley AOI."""
    assert atlas14.region_for_bounds((-119.75, 39.33, -119.45, 39.88)) == "sw"


def test_region_lookup_returns_a_code_for_other_conus_aois():
    """Texas and the Ohio River basin must resolve to some region, not crash."""
    for bounds in ((-97.9, 30.1, -97.5, 30.5), (-84.6, 39.0, -84.3, 39.3)):
        code = atlas14.region_for_bounds(bounds)
        assert isinstance(code, str) and code


# --------------------------------------------------------------------------- #
# grid URL construction
# --------------------------------------------------------------------------- #
def test_grid_url_encodes_region_ari_duration_and_stat():
    url = atlas14.grid_url("sw", 1, 15)
    assert url.startswith("http")
    assert "sw" in url and url.endswith(".zip")


def test_grid_url_distinguishes_ari_and_duration():
    assert atlas14.grid_url("sw", 1, 15) != atlas14.grid_url("sw", 100, 15)
    assert atlas14.grid_url("sw", 1, 15) != atlas14.grid_url("sw", 1, 60)


def test_grid_url_distinguishes_the_ci_bounds():
    """mean / 90% lower / 90% upper are three different products."""
    urls = {atlas14.grid_url("sw", 1, 15, stat=s) for s in ("a", "al", "au")}
    assert len(urls) == 3
