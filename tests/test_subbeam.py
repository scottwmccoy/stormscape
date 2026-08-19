"""Sub-beam evaporation: model limits, sounding parse, raster wrapper."""
from __future__ import annotations

import datetime as dt

import numpy as np
import pytest

from stormscape import subbeam
from stormscape.subbeam import evap_factor, mean_rh, parse_wyoming


def test_saturated_air_loses_nothing():
    assert evap_factor(20.0, 1.0, rh=1.0) == pytest.approx(1.0)


def test_zero_depth_loses_nothing():
    assert evap_factor(20.0, 0.0, rh=0.2) == pytest.approx(1.0)


def test_drier_air_loses_more():
    f_dry = evap_factor(10.0, 1.0, rh=0.2)
    f_moist = evap_factor(10.0, 1.0, rh=0.7)
    assert f_dry < f_moist < 1.0


def test_light_rain_loses_proportionally_more_than_heavy():
    assert evap_factor(2.0, 1.0, rh=0.2) < evap_factor(50.0, 1.0, rh=0.2)


def test_calibration_anchor():
    """~10 mm/h through 1 km of RH=20% air loses ~35% (Rosenfeld & Mintz 88)."""
    f = evap_factor(10.0, 1.0, rh=0.2)
    assert 0.55 < f < 0.75


def test_complete_evaporation_floors_at_zero():
    assert evap_factor(0.3, 3.0, rh=0.05) == 0.0


_FIXTURE = """
-----------------------------------------------------------------------------
   PRES   HGHT   TEMP   DWPT   RELH   MIXR   DRCT   SKNT   THTA   THTE   THTV
    hPa     m      C      C      %    g/kg    deg   knot     K      K      K
-----------------------------------------------------------------------------
 1000.0    111
  850.0   1509   28.0    2.0   19.0   5.50    180     10  320.0  339.0  321.0
  700.0   3192   12.0   -3.0   35.0   4.20    200     15  322.0  337.0  323.0
  500.0   5880   -8.0  -18.0   45.0   1.80    220     25  325.0  332.0  325.0
"""


def test_parse_wyoming_fixture():
    df = parse_wyoming(_FIXTURE)
    assert len(df) == 3                       # the 2-column line is skipped
    assert df.z_m.tolist() == [1509.0, 3192.0, 5880.0]
    assert df.rh.iloc[0] == pytest.approx(0.19)


def test_parse_wyoming_rejects_garbage():
    with pytest.raises(ValueError, match="TEXT:LIST"):
        parse_wyoming("<html>not a sounding</html>")


def test_mean_rh_layer():
    df = parse_wyoming(_FIXTURE)
    assert mean_rh(df, 1500, 3200) == pytest.approx((0.19 + 0.35) / 2)


def _write(path, arr, res_deg, crs="EPSG:4326"):
    import rasterio
    tr = rasterio.transform.from_origin(-119.6, 39.9, res_deg, res_deg)
    with rasterio.open(path, "w", driver="GTiff", height=arr.shape[0],
                       width=arr.shape[1], count=1, dtype="float32",
                       crs=crs, transform=tr, nodata=np.nan) as ds:
        ds.write(arr.astype("float32"), 1)
    return tr


def test_subbeam_correct_geometry_and_output(tmp_path):
    """Loss grows with range (deeper sub-beam column) at fixed RH."""
    import rasterio
    rate = np.full((20, 20), 20.0)
    dem = np.full((20, 20), 1500.0)
    _write(tmp_path / "rate.tif", rate, 0.01)
    _write(tmp_path / "dem.tif", dem, 0.01)
    res = subbeam.subbeam_correct(str(tmp_path / "rate.tif"),
                                  str(tmp_path / "dem.tif"),
                                  radar_lonlat=(-119.6, 39.9),
                                  radar_elev_m=2500.0,
                                  out_dir=str(tmp_path), key="t",
                                  rh=0.25, tilt_deg=0.0)
    with rasterio.open(res["subbeam_tif"]) as ds:
        out = ds.read(1)
        assert ds.tags()["MODEL"].startswith("dR/dz")
    # radar at the NW origin cell: nearby cells lose less than the far corner
    assert out[0, 0] > out[19, 19]
    assert (out <= rate + 1e-6).all()
    assert res["loss_pct_wet_cells"]["median"] > 0
    assert "rasters" in res["subbeam_tif"]


@pytest.mark.network
def test_fetch_sounding_live():
    df = subbeam.fetch_sounding("REV", dt.datetime(2026, 8, 14, 22))
    assert len(df) > 20 and df.rh.between(0, 1).all()


# --------------------------------------------------------------------------- #
# R(Z,ZDR) -- the third leg of the estimator triangle (lives in nexrad.py)
# --------------------------------------------------------------------------- #
def test_zzdr_deflates_big_drop_cores():
    """At 50 dBZ, ZDR 2.3 dB (the observed Stallion core) must read well below
    the fixed convective Z-R -- that differential is the whole point."""
    from stormscape.nexrad import z_to_rate, zzdr_to_rate
    r_zr = z_to_rate(50.0, dbz_cap=53.0)
    r_zzdr = zzdr_to_rate(50.0, 2.3, dbz_cap=53.0)
    assert r_zzdr < 0.75 * r_zr


def test_zzdr_boosts_small_drop_rain():
    """Low ZDR at a given Z = many small drops = MORE rain than Z-R assumes."""
    from stormscape.nexrad import z_to_rate, zzdr_to_rate
    assert zzdr_to_rate(40.0, 0.2) > z_to_rate(40.0)


def test_zzdr_monotone_decreasing_in_zdr():
    from stormscape.nexrad import zzdr_to_rate
    r = [float(zzdr_to_rate(45.0, z)) for z in (0.0, 1.0, 2.0, 3.0, 4.0)]
    assert all(a > b for a, b in zip(r, r[1:]))


def test_zzdr_nan_zdr_yields_nan_for_fallback():
    from stormscape.nexrad import zzdr_to_rate
    out = zzdr_to_rate(np.array([45.0, 45.0]), np.array([1.0, np.nan]))
    assert np.isfinite(out[0]) and np.isnan(out[1])


def test_zzdr_clips_nonphysical_zdr():
    """ZDR below 0 / above 4 dB is noise or calibration: clipped, not obeyed."""
    from stormscape.nexrad import zzdr_to_rate
    assert zzdr_to_rate(45.0, -2.0) == pytest.approx(zzdr_to_rate(45.0, 0.0))
    assert zzdr_to_rate(45.0, 9.0) == pytest.approx(zzdr_to_rate(45.0, 4.0))


def test_zzdr_hail_cap_applies():
    from stormscape.nexrad import zzdr_to_rate
    assert zzdr_to_rate(63.0, 2.0, dbz_cap=53.0) == \
        pytest.approx(zzdr_to_rate(53.0, 2.0))


def test_rh_from_t_td_limits():
    from stormscape.subbeam import rh_from_t_td
    assert rh_from_t_td(20.0, 20.0) == pytest.approx(1.0)      # saturated
    assert 0.15 < rh_from_t_td(28.0, 0.0) < 0.20               # hot + dry
