"""Shaded relief: the four rules in stormscape.relief's docstring, pinned.

Every failure mode here was observed in a real statewide figure on 2026-08-14:
mushy terrain from shading at display resolution, a grey block where tile
coverage ended, and — nearly shipped — banding seams from LightSource's
per-array contrast stretch.

Offline: synthetic geographic tiles written to tmp_path stand in for 3DEP.
"""

import numpy as np
import pytest
import rasterio
from matplotlib.colors import LightSource
from rasterio.transform import from_origin

from stormscape import relief

RES = relief.SHADE_RES_M

#: Synthetic tile pixel size. ~1 arcsecond (~26 m at 40 N) -- coarser than real
#: 3DEP 1/3-arcsecond tiles, but still *finer* than the shading grid, which is
#: the case that matters: shading must downsample its source, never upsample it.
TILE_DEG = 0.0003


def _write(path, arr, *, west, north, res_deg, nodata=None):
    """A synthetic 3DEP-like tile: geographic CRS, north-up."""
    h, w = arr.shape
    with rasterio.open(
        path, "w", driver="GTiff", height=h, width=w, count=1,
        dtype="float32", crs="EPSG:4326", nodata=nodata,
        transform=from_origin(west, north, res_deg, res_deg),
    ) as ds:
        ds.write(arr.astype("float32"), 1)
    return path


def _rough(shape, seed=0, amp=8.0):
    """A rough but smooth-ish surface with texture at every scale."""
    rng = np.random.default_rng(seed)
    z = rng.standard_normal(shape).cumsum(0).cumsum(1) * 0.02
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    return (z + amp * np.sin(xx / 9.0) * np.cos(yy / 11.0)).astype("float32")


# --- rule 4: a fixed 0-1 map, so splitting the work cannot change the result -

def test_flat_ground_is_sin_altitude():
    hs = relief.shade(np.zeros((6, 6), "float32"), RES)
    assert np.allclose(hs, np.sin(np.radians(relief.ALTITUDE_DEG)))


def test_planar_slope_matches_hand_computation():
    """Known answer: on a plane, shading is closed-form."""
    a, b, res = 0.10, -0.04, 30.0            # dz/dx, dz/dy
    h, w = 12, 15
    rr, cc = np.mgrid[0:h, 0:w]
    z = (a * (cc * res) + b * ((h - 1 - rr) * res)).astype("float32")

    d = relief.light_direction()
    v = relief.VERT_EXAG
    expect = (-v * a * d[0] - v * b * d[1] + d[2]) / np.sqrt((v * a) ** 2
                                                            + (v * b) ** 2 + 1)
    hs = relief.shade(z, res)
    assert np.allclose(hs[1:-1, 1:-1], expect, atol=1e-6)


@pytest.mark.parametrize("band_rows", [1, 7, 40, 4000])
def test_shade_is_band_split_invariant(band_rows):
    """The banding regression. Shading in bands with a one-row halo must equal
    shading whole -- exactly, not approximately."""
    z = _rough((97, 31), seed=3)
    whole = relief.shade(z, RES)
    out = np.empty_like(whole)
    for y0 in range(0, z.shape[0], band_rows):
        y1 = min(z.shape[0], y0 + band_rows)
        top = max(0, y0 - 1)
        band = relief.shade(z[top:min(z.shape[0], y1 + 1)], RES)
        out[y0:y1] = band[y0 - top:y0 - top + (y1 - y0)]
    assert np.array_equal(whole, out)


def test_lightsource_is_not_band_split_invariant():
    """Why rule 4 exists: the same split through LightSource does *not* agree,
    because its contrast stretch is rescaled per array."""
    z = _rough((97, 31), seed=3)
    ls = LightSource(azdeg=relief.AZIMUTH_DEG, altdeg=relief.ALTITUDE_DEG)
    kw = dict(vert_exag=relief.VERT_EXAG, dx=RES, dy=RES)
    whole = ls.hillshade(z, **kw)
    top_half = ls.hillshade(z[:49], **kw)
    assert not np.allclose(whole[:48], top_half[:48], atol=1e-3)


# --- rule 2: mosaic elevations before shading, so tile edges leave no trace ---

def test_tile_seam_is_invisible(tmp_path):
    """Two abutting tiles must give the same relief as one tile of the whole
    surface. Shading tiles separately would print a line at the join."""
    z = _rough((120, 240), seed=11, amp=25.0)
    d = TILE_DEG
    one = _write(tmp_path / "one.tif", z, west=-119.0, north=40.0, res_deg=d)
    left = _write(tmp_path / "l.tif", z[:, :120], west=-119.0, north=40.0, res_deg=d)
    right = _write(tmp_path / "r.tif", z[:, 120:], west=-119.0 + 120 * d,
                   north=40.0, res_deg=d)

    # display grid, deliberately coarser than the shading grid
    tr = from_origin(-119.0, 40.0, 5 * d, 5 * d)
    shape = (24, 48)
    kw = dict(crs="EPSG:4326", shade_res_m=RES, band_rows=7)
    whole = relief.shaded_relief([one], tr, shape, **kw)
    split = relief.shaded_relief([left, right], tr, shape, **kw)

    assert np.abs(whole - split).max() < 0.02
    # and no localised spike at the join column
    col = np.abs(whole - split).mean(axis=0)
    assert col.max() < 5 * np.median(col) + 0.01


# --- rule 3: shade fine, then average down, or the texture is gone ----------

def test_fine_shading_keeps_more_texture(tmp_path):
    z = _rough((240, 240), seed=5, amp=25.0)
    d = TILE_DEG
    t = _write(tmp_path / "t.tif", z, west=-119.0, north=40.0, res_deg=d)
    tr = from_origin(-119.0, 40.0, 8 * d, 8 * d)
    shape = (30, 30)

    fine = relief.shaded_relief([t], tr, shape, shade_res_m=40.0, band_rows=64)
    coarse = relief.shaded_relief([t], tr, shape, shade_res_m=400.0, band_rows=64)

    def roughness(hs):
        lap = (hs[1:-1, 2:] + hs[1:-1, :-2] + hs[2:, 1:-1] + hs[:-2, 1:-1]
               - 4 * hs[1:-1, 1:-1])
        return float(np.abs(lap).mean())

    assert roughness(fine) > 1.5 * roughness(coarse)


# --- coverage gaps read as flat ground, never as a block or a cliff ---------

def test_uncovered_ground_is_flat_value(tmp_path):
    """A tile covering only the west half; the east half must come back at
    exactly sin(altitude), with no cliff at the boundary."""
    z = _rough((200, 200), seed=7, amp=25.0)
    d = TILE_DEG
    t = _write(tmp_path / "half.tif", z, west=-119.0, north=40.0, res_deg=d)
    tr = from_origin(-119.0, 40.0, 5 * d, 5 * d)
    hs = relief.shaded_relief([t], tr, (40, 80), shade_res_m=RES, band_rows=9)

    flat = np.float32(np.sin(np.radians(relief.ALTITUDE_DEG)))
    assert np.allclose(hs[:, 55:], flat, atol=1e-6)   # well past the edge
    assert hs.min() >= 0.0 and hs.max() <= 1.0
    # the covered half must still carry real relief
    assert hs[:, :30].std() > 0.02


def test_no_tiles_returns_flat(tmp_path):
    tr = from_origin(-119.0, 40.0, 0.01, 0.01)
    hs = relief.shaded_relief([], tr, (8, 8))
    assert np.allclose(hs, np.sin(np.radians(relief.ALTITUDE_DEG)))


def test_tile_outside_grid_is_skipped(tmp_path):
    z = _rough((40, 40), seed=1)
    far = _write(tmp_path / "far.tif", z, west=-100.0, north=35.0, res_deg=0.002)
    tr = from_origin(-119.0, 40.0, 0.008, 0.008)
    hs = relief.shaded_relief([far], tr, (10, 10))
    assert np.allclose(hs, np.sin(np.radians(relief.ALTITUDE_DEG)))


def test_source_nodata_does_not_survive_as_elevation(tmp_path):
    """Explicit nodata must become NaN before the warp -- otherwise -9999 warps
    in as a canyon and the hillshade grows a black gash."""
    z = _rough((200, 200), seed=2, amp=25.0) + 1500.0
    z[60:140, 60:140] = -9999.0
    d = TILE_DEG
    t = _write(tmp_path / "nd.tif", z, west=-119.0, north=40.0,
               res_deg=d, nodata=-9999.0)
    tr = from_origin(-119.0, 40.0, 5 * d, 5 * d)
    hs = relief.shaded_relief([t], tr, (40, 40), shade_res_m=RES, band_rows=8)

    flat = np.float32(np.sin(np.radians(relief.ALTITUDE_DEG)))
    assert np.allclose(hs[18:22, 18:22], flat, atol=1e-6)
    assert hs[:10, :10].std() > 0.01        # real terrain outside the hole
