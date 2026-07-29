"""NaN-aware field smoothing: the invariants that keep a smoothed field honest."""
from __future__ import annotations

import numpy as np
import pytest

from stormscape import smoothing

METHODS = list(smoothing.METHODS)


@pytest.fixture
def peaky():
    rng = np.random.RandomState(0)
    return rng.rand(41, 41).astype(float) * 100.0


@pytest.mark.parametrize("method", METHODS)
def test_radius_zero_is_the_identity(method, peaky):
    """radius_km=0 must return the field untouched, for every method."""
    assert np.allclose(smoothing.smooth_array(peaky, 1.0, method, 0.0), peaky)


@pytest.mark.parametrize("method", METHODS)
def test_smoothing_never_raises_the_peak(method, peaky):
    """Low-pass filters cannot amplify: this underpins the 'smoothing lowers the
    radar's positive i15 bias mechanically' caveat in the docs."""
    out = smoothing.smooth_array(peaky, 1.0, method, 2.0)
    assert np.nanmax(out) <= np.nanmax(peaky) + 1e-9


@pytest.mark.parametrize("method", METHODS)
def test_peak_reduction_is_monotone_in_radius(method, peaky):
    peaks = [np.nanmax(smoothing.smooth_array(peaky, 1.0, method, r))
             for r in (0.0, 1.0, 2.0, 4.0)]
    assert all(b <= a + 1e-9 for a, b in zip(peaks, peaks[1:]))


@pytest.mark.parametrize("method", METHODS)
def test_constant_field_with_a_nan_hole_does_not_zero_bleed(method):
    """The normalized-convolution guarantee: a masked cell must not drag its
    neighbours toward zero. This is the bug that silently deflates fields."""
    a = np.full((21, 21), 7.0)
    a[10, 10] = np.nan
    out = smoothing.smooth_array(a, 1.0, method, 2.0)
    assert np.allclose(out[np.isfinite(out)], 7.0, atol=1e-6)


@pytest.mark.parametrize("method", METHODS)
def test_all_nan_field_stays_all_nan(method):
    out = smoothing.smooth_array(np.full((9, 9), np.nan), 1.0, method, 2.0)
    assert np.isnan(out).all()


@pytest.mark.parametrize("method", METHODS)
def test_mean_is_approximately_preserved(method):
    """A low-pass filter redistributes mass; it should not create or destroy much."""
    a = np.full((31, 31), 5.0)
    out = smoothing.smooth_array(a, 1.0, method, 2.0)
    assert np.nanmean(out) == pytest.approx(5.0, rel=1e-6)


@pytest.mark.parametrize("method", ["gaussian", "uniform", "idw"])
def test_delta_spike_smooths_symmetrically(method):
    """Catches an even-window origin shift: the response to a centred spike must
    be symmetric about that centre."""
    a = np.zeros((21, 21))
    a[10, 10] = 100.0
    out = smoothing.smooth_array(a, 1.0, method, 2.0)
    assert np.allclose(out, out[::-1, :], atol=1e-9)      # up/down
    assert np.allclose(out, out[:, ::-1], atol=1e-9)      # left/right


def test_median_removes_an_isolated_spike():
    """The rank filter's purpose: despeckle without smearing."""
    a = np.full((15, 15), 4.0)
    a[7, 7] = 500.0
    out = smoothing.smooth_array(a, 1.0, "median", 2.0)
    assert out[7, 7] == pytest.approx(4.0)


def test_unknown_method_is_rejected():
    with pytest.raises((KeyError, ValueError)):
        smoothing.smooth_array(np.zeros((5, 5)), 1.0, "bicubic-wizardry", 1.0)


# --------------------------------------------------------------------------- #
# km -> pixel conversion
# --------------------------------------------------------------------------- #
def test_cell_size_km_on_the_mrms_grid():
    """0.01° at ~39.5°N is ~0.99 km — the documented MRMS cell size."""
    from rasterio.transform import from_origin
    tr = from_origin(-119.75, 39.60, 0.01, 0.01)
    assert smoothing.cell_size_km(tr, 39.5) == pytest.approx(0.99, abs=0.06)


def test_cell_size_km_shrinks_with_latitude():
    """Longitude degrees converge poleward, so the same grid has smaller cells."""
    from rasterio.transform import from_origin
    tr = from_origin(-119.75, 39.60, 0.01, 0.01)
    assert smoothing.cell_size_km(tr, 70.0) < smoothing.cell_size_km(tr, 10.0)


# --------------------------------------------------------------------------- #
# raster / DataArray entry points
# --------------------------------------------------------------------------- #
def test_smooth_tif_preserves_grid_and_crs(field_tif, tmp_path):
    import rasterio
    src = field_tif(np.random.RandomState(1).rand(20, 15) * 50)
    out = smoothing.smooth_tif(src, str(tmp_path / "sm.tif"), "gaussian", 1.0)
    with rasterio.open(src) as a, rasterio.open(out) as b:
        assert (a.width, a.height) == (b.width, b.height)
        assert a.crs == b.crs
        assert a.transform == b.transform


def test_smooth_dataarray_keeps_the_input_dim_order(field_tif):
    """Must stay (band, y, x): a squeezed 2-D result breaks rio.to_raster later
    (the InvalidDimensionOrder failure that only surfaced at write time)."""
    import rioxarray
    src = field_tif(np.random.RandomState(2).rand(12, 9) * 30)
    with rioxarray.open_rasterio(src, masked=True) as da:
        da = da.load()
    out = smoothing.smooth_dataarray(da, "gaussian", 1.0)
    assert out.dims == da.dims
    assert out.shape == da.shape


def test_smooth_dataarray_accepts_a_path(field_tif):
    src = field_tif(np.full((10, 10), 3.0))
    out = smoothing.smooth_dataarray(src, "gaussian", 1.0)
    assert float(np.nanmean(out.values)) == pytest.approx(3.0, rel=1e-6)
