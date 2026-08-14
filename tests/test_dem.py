"""3DEP fetch routing and the single-warp contract.

These pin the fix for the nearest-neighbour resampling artefact: at 10/30/60 m
``py3dep.get_dem`` nearest-warps the seamless EPSG:4269 VRT to 5070 (a rotation
plus a non-integer scale, so nearest aliases along diagonals), and a second warp
to the requested resolution then drops a row and column every ~15 cells. The
contract is: read the VRT on its native grid, warp exactly once, interpolate.

No network -- ``py3dep`` is replaced by a recording double.
"""
from __future__ import annotations

import numpy as np
import pytest
import rioxarray  # noqa: F401  (registers the .rio accessor)
import xarray as xr
from rasterio.enums import Resampling
from rioxarray.raster_array import RasterArray

from stormscape import dem as dem_mod


# --------------------------------------------------------------------------- #
# a py3dep double that records which endpoint was asked for
# --------------------------------------------------------------------------- #
def _grid(west, north, res, n=24, crs="EPSG:4269"):
    """A tiny DEM on a regular grid, georeferenced like the real thing."""
    x = west + res * (np.arange(n) + 0.5)
    y = north - res * (np.arange(n) + 0.5)
    z = xr.DataArray(np.linspace(1000, 1100, n * n).reshape(n, n).astype("f4"),
                     dims=("y", "x"), coords={"y": y, "x": x})
    z.rio.write_crs(crs, inplace=True)
    return z


class FakePy3dep:
    def __init__(self, static_grid=None, dynamic_grid=None):
        self.calls = []
        self._static = static_grid
        self._dynamic = dynamic_grid

    def static_3dep_dem(self, geometry, crs, resolution=10):
        self.calls.append(("static_3dep_dem", resolution))
        # the real one returns the VRT's native 1/3" EPSG:4269 grid
        return (self._static if self._static is not None
                else _grid(-111.4, 33.85, 1 / 10800.0))

    def get_dem(self, geometry, resolution, crs=4326):
        self.calls.append(("get_dem", resolution))
        # the dynamic service resamples server-side and returns 5070 metres
        return (self._dynamic if self._dynamic is not None
                else _grid(-1408797.0, 1310881.0, float(resolution),
                           crs="EPSG:5070"))


@pytest.fixture
def fake3dep(monkeypatch):
    def install(**kw):
        fake = FakePy3dep(**kw)
        monkeypatch.setattr(dem_mod, "_py3dep", lambda: fake)
        return fake
    return install


AOI = (-111.40, 33.80, -111.35, 33.85)


# --------------------------------------------------------------------------- #
# routing: which 3DEP endpoint gets used
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("res", dem_mod.STATIC_RESOLUTIONS)
def test_static_resolutions_bypass_py3dep_get_dem(fake3dep, res):
    """10/30/60 m must read the VRT directly -- py3dep.get_dem is the call that
    hides a nearest warp, so it must not appear."""
    fake = fake3dep()
    dem_mod.get_dem(AOI, resolution=res)
    assert [c[0] for c in fake.calls] == ["static_3dep_dem"]


@pytest.mark.parametrize("res", [1, 3, 5])
def test_other_resolutions_use_the_dynamic_service(fake3dep, res):
    """Outside 10/30/60 py3dep routes to the dynamic image service, which
    resamples server-side and never had the artefact."""
    fake = fake3dep()
    dem_mod.get_dem(AOI, resolution=res)
    assert [c[0] for c in fake.calls] == ["get_dem"]


# --------------------------------------------------------------------------- #
# the single-warp contract
# --------------------------------------------------------------------------- #
def test_static_path_lands_on_the_requested_resolution(fake3dep):
    fake3dep()
    out = dem_mod.get_dem(AOI, resolution=10)
    assert out.rio.crs.to_epsg() == 5070
    assert abs(out.rio.resolution()[0]) == pytest.approx(10.0)


def test_warp_is_never_nearest_by_default(fake3dep, monkeypatch):
    """Nearest is only correct for categorical rasters; on elevation it is what
    produced the diagonal hatch."""
    fake3dep()
    seen = {}
    orig = RasterArray.reproject

    def spy(self, *a, **kw):
        seen["resampling"] = kw.get("resampling")
        return orig(self, *a, **kw)

    monkeypatch.setattr(RasterArray, "reproject", spy)
    dem_mod.get_dem(AOI, resolution=10)
    assert seen["resampling"] is Resampling.bilinear


def test_no_warp_when_already_on_the_target_grid(fake3dep, monkeypatch):
    """The dynamic service already returns 5070 at the requested spacing.
    Reprojecting a raster onto the grid it is already on is the exact mistake
    this module exists to avoid, so it must be skipped."""
    fake3dep()
    calls = []
    orig = RasterArray.reproject

    def spy(self, *a, **kw):
        calls.append(kw.get("resampling"))
        return orig(self, *a, **kw)

    monkeypatch.setattr(RasterArray, "reproject", spy)
    dem_mod.get_dem(AOI, resolution=5)          # dynamic double returns 5 m/5070
    assert calls == []


def test_warp_happens_when_the_crs_differs(fake3dep):
    fake3dep()
    out = dem_mod.get_dem(AOI, resolution=5, dst_crs="EPSG:26912")
    assert out.rio.crs.to_epsg() == 26912


def test_dst_crs_none_leaves_the_native_grid_alone(fake3dep):
    fake3dep()
    out = dem_mod.get_dem(AOI, resolution=5, dst_crs=None)
    assert out.rio.crs.to_epsg() == 5070         # whatever it arrived on


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def test_resampling_accepts_a_name_or_a_member():
    assert dem_mod._resolve_resampling("cubic") is Resampling.cubic
    assert dem_mod._resolve_resampling(Resampling.average) is Resampling.average


def test_unknown_resampling_is_rejected_by_name():
    with pytest.raises(ValueError, match="unknown resampling"):
        dem_mod._resolve_resampling("bilinearish")


def test_needs_warp_is_false_only_on_a_matching_grid():
    g = _grid(-1408797.0, 1310881.0, 10.0, crs="EPSG:5070")
    assert not dem_mod._needs_warp(g, "EPSG:5070", 10)
    assert dem_mod._needs_warp(g, "EPSG:5070", 30)      # wrong spacing
    assert dem_mod._needs_warp(g, "EPSG:26912", 10)     # wrong CRS
    assert not dem_mod._needs_warp(g, None, 10)         # asked to leave it


def test_needs_warp_tolerates_sub_percent_spacing_drift():
    """Float round-trips through a transform, so an exact == would warp for no
    reason; 1% is far tighter than the 6.6% mismatch that caused the artefact."""
    assert not dem_mod._needs_warp(
        _grid(-1408797.0, 1310881.0, 10.004, crs="EPSG:5070"), "EPSG:5070", 10)
    assert dem_mod._needs_warp(
        _grid(-1408797.0, 1310881.0, 9.3817, crs="EPSG:5070"), "EPSG:5070", 10)


# --------------------------------------------------------------------------- #
# hillshade, unchanged but worth pinning alongside
# --------------------------------------------------------------------------- #
def test_hillshade_is_0_255_on_the_dem_grid():
    g = _grid(-1408797.0, 1310881.0, 10.0, crs="EPSG:5070")
    hs = dem_mod.hillshade(g)
    assert hs.shape == g.shape
    assert hs.rio.resolution() == g.rio.resolution()
    assert float(hs.min()) >= 0.0 and float(hs.max()) <= 255.0
