"""Hydro blend: the per-gate relation selector and its envelope."""
from __future__ import annotations

import numpy as np
import pytest

from stormscape.nexrad import (HYDRO_CENSORED, HYDRO_DRY, HYDRO_KDP_HAIL,
                               HYDRO_ZR, HYDRO_ZZDR, _KDP_GUARD, hydro_select,
                               z_to_rate, zzdr_to_rate)


def _sel(dbz, zdr, rhohv=None, kdp_rate=None):
    """Run hydro_select on one cell with self-consistent member rates."""
    dbz = np.array([dbz], dtype=float)
    zdr = np.array([zdr], dtype=float) if zdr is not None else None
    r_z = z_to_rate(dbz, dbz_cap=53.0)
    r_zzdr = (zzdr_to_rate(dbz, zdr, dbz_cap=53.0) if zdr is not None
              else np.full(1, np.nan))
    r_kdp = np.array([kdp_rate], dtype=float) if kdp_rate is not None \
        else np.full(1, np.nan)
    cc = np.array([rhohv], dtype=float) if rhohv is not None else None
    rate, lo, hi, cls = hydro_select(dbz, zdr, cc, r_z, r_kdp, r_zzdr)
    return float(rate[0]), float(lo[0]), float(hi[0]), int(cls[0]), float(r_z[0])


def test_light_rain_uses_plain_zr():
    rate, lo, hi, cls, r_z = _sel(28.0, 1.0)
    assert cls == HYDRO_ZR and rate == pytest.approx(r_z)
    assert lo == hi == pytest.approx(r_z)          # single-member envelope


def test_big_drop_rain_uses_zzdr():
    """The 13 Aug core signature: 50 dBZ, ZDR 2.3 -> R(Z,ZDR), below Z-R."""
    rate, lo, hi, cls, r_z = _sel(50.0, 2.3)
    assert cls == HYDRO_ZZDR
    assert rate < r_z
    assert lo <= rate <= hi and hi == pytest.approx(r_z)


def test_hail_routes_to_kdp_and_excludes_zzdr_from_envelope():
    """The 12 Aug core signature: 55 dBZ, ZDR 0.2 -> R(Kdp); R(Z,ZDR) would
    blow up there and must not set the envelope ceiling."""
    rate, lo, hi, cls, r_z = _sel(55.0, 0.2, kdp_rate=40.0)
    assert cls == HYDRO_KDP_HAIL and rate == pytest.approx(40.0)
    r_zzdr_would_be = float(zzdr_to_rate(np.array([55.0]), np.array([0.2]),
                                         dbz_cap=53.0)[0])
    assert r_zzdr_would_be > r_z                   # the blow-up is real...
    assert hi <= r_z + 1e-9                        # ...and excluded from hi


def test_kdp_guard_caps_delta_blowups_in_hail_branch():
    rate, lo, hi, cls, r_z = _sel(55.0, 0.2, kdp_rate=10.0 * 104.0)
    assert cls == HYDRO_KDP_HAIL
    assert rate == pytest.approx(_KDP_GUARD * r_z)


def test_hail_with_invalid_kdp_falls_to_zr():
    rate, lo, hi, cls, r_z = _sel(55.0, 0.2, kdp_rate=None)
    assert cls == HYDRO_KDP_HAIL and rate == pytest.approx(r_z)


def test_marginal_z_never_selects_kdp():
    """The 13 Aug West Valley failure: raging Kdp at ~41 dBZ with real ZDR --
    the rain branch must ignore it entirely."""
    rate, lo, hi, cls, _ = _sel(41.0, 1.5, kdp_rate=73.0)
    assert cls == HYDRO_ZZDR and rate < 30.0


def test_low_rhohv_censors_to_zero():
    rate, lo, hi, cls, _ = _sel(45.0, 1.0, rhohv=0.70)
    assert cls == HYDRO_CENSORED and rate == lo == hi == 0.0


def test_nan_reflectivity_stays_nan():
    rate, lo, hi, cls, _ = _sel(np.nan, 1.0)
    assert cls == HYDRO_DRY and np.isnan(rate) and np.isnan(lo)


def test_single_pol_degrades_to_zr():
    rate, lo, hi, cls, r_z = _sel(50.0, None, kdp_rate=None)
    assert cls == HYDRO_ZR and rate == pytest.approx(r_z)


def test_envelope_brackets_the_blend_everywhere():
    rng = np.random.default_rng(7)
    dbz = rng.uniform(10, 60, 500)
    zdr = rng.uniform(-0.5, 4.5, 500)
    cc = rng.uniform(0.8, 1.0, 500)
    r_z = z_to_rate(dbz, dbz_cap=53.0)
    r_zzdr = zzdr_to_rate(dbz, zdr, dbz_cap=53.0)
    r_kdp = rng.uniform(0, 150, 500)
    rate, lo, hi, cls = hydro_select(dbz, zdr, cc, r_z, r_kdp, r_zzdr)
    ok = np.isfinite(rate)
    assert (lo[ok] <= rate[ok] + 1e-9).all()
    assert (rate[ok] <= hi[ok] + 1e-9).all()
