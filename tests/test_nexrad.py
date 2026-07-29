"""Single-radar NEXRAD: the Z-R relation, the hail cap, and the site table.

Only the parts that need no Py-ART/nexradaws — those are lazily imported, and the
volume-reading paths need network access anyway.
"""
from __future__ import annotations

import numpy as np
import pytest

from stormscape import nexrad


# --------------------------------------------------------------------------- #
# Z-R: R = (z/a)^(1/b) with z = 10^(dBZ/10)
# --------------------------------------------------------------------------- #
def test_z_to_rate_matches_the_closed_form():
    dbz = np.array([20.0, 35.0, 50.0])
    z = 10.0 ** (dbz / 10.0)
    assert np.allclose(nexrad.z_to_rate(dbz), (z / 300.0) ** (1.0 / 1.4))


def test_z_to_rate_is_monotone_in_reflectivity():
    r = nexrad.z_to_rate(np.array([10.0, 25.0, 40.0, 55.0, 70.0]))
    assert all(b > a for a, b in zip(r, r[1:]))


def test_z_to_rate_honours_alternative_coefficients():
    """The Beard-Chuang / Marshall-Palmer style swaps must actually take effect."""
    dbz = np.array([45.0])
    assert nexrad.z_to_rate(dbz, a=200.0, b=1.6)[0] != pytest.approx(
        nexrad.z_to_rate(dbz)[0])


def test_dbz_cap_limits_the_rate_at_hail_reflectivities():
    """The operational hail cap: without it, 70 dBZ converts to a physically absurd
    rain rate. The gap between capped and uncapped at high dBZ *is* the hail
    over-estimation this project quantified."""
    capped = nexrad.z_to_rate(np.array([70.0]), dbz_cap=53.0)[0]
    at_cap = nexrad.z_to_rate(np.array([53.0]))[0]
    uncapped = nexrad.z_to_rate(np.array([70.0]))[0]
    assert capped == pytest.approx(at_cap)
    assert uncapped > 3 * capped


def test_dbz_cap_leaves_light_rain_untouched():
    dbz = np.array([25.0, 35.0])
    assert np.allclose(nexrad.z_to_rate(dbz, dbz_cap=53.0),
                       nexrad.z_to_rate(dbz))


def test_z_to_rate_propagates_nan():
    out = nexrad.z_to_rate(np.array([np.nan, 40.0]))
    assert np.isnan(out[0]) and np.isfinite(out[1])


def test_z_to_rate_accepts_a_scalar():
    assert float(np.asarray(nexrad.z_to_rate(40.0))) > 0


# --------------------------------------------------------------------------- #
# bundled WSR-88D site table (no network)
# --------------------------------------------------------------------------- #
def test_nearest_radar_to_the_reno_aoi_is_krgx():
    """KRGX is the mountaintop radar serving the Hidden Valley AOI.
    Return contract is ``(id, dist_km, lat, lon)``."""
    rid, dist_km, lat, lon = nexrad.nearest_radar((-119.75, 39.33, -119.45, 39.88))
    assert rid == "KRGX"
    assert 0 < dist_km < 60                      # ~20 km from the AOI centroid
    assert 39.0 < lat < 40.5 and -120.5 < lon < -119.0


@pytest.mark.parametrize("aoi, expected", [
    ((-104.8, 39.6, -104.5, 39.9), "KFTG"),    # Denver
    ((-97.6, 35.3, -97.3, 35.6), "KTLX"),      # Oklahoma City
])
def test_nearest_radar_for_other_known_cities(aoi, expected):
    assert nexrad.nearest_radar(aoi)[0] == expected


def test_nearest_radar_returns_plausible_coordinates_anywhere_in_conus():
    rid, dist_km, lat, lon = nexrad.nearest_radar((-90.3, 38.5, -90.1, 38.7))
    assert isinstance(rid, str) and len(rid) == 4
    assert 24.0 < lat < 50.0 and -125.0 < lon < -66.0
    assert dist_km >= 0


def test_nearest_radar_distance_grows_with_offset():
    """Sanity on the distance term: an AOI further from KRGX must report further."""
    near = nexrad.nearest_radar((-119.50, 39.70, -119.45, 39.80))[1]
    far = nexrad.nearest_radar((-119.75, 39.33, -119.70, 39.38))[1]
    assert far > near


def test_radar_location_agrees_with_nearest_radar():
    """``radar_location`` returns ``(lat, lon, elev_m)`` for the same site."""
    rid, _, lat, lon = nexrad.nearest_radar((-119.75, 39.33, -119.45, 39.88))
    rlat, rlon, relev = nexrad.radar_location(rid)
    assert (rlat, rlon) == pytest.approx((lat, lon))
    assert relev > 2000                          # KRGX sits at ~2559 m


def test_radar_location_rejects_an_unknown_id():
    with pytest.raises(Exception):
        nexrad.radar_location("ZZZZ")
