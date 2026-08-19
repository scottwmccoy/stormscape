"""Shaded relief for map backdrops, built from a local 3DEP tile store.

Lives here rather than in a downstream project because a terrain backdrop is
not specific to any one science question: firescape draws hazard over it and
tracescape draws change over it, and neither should own the other's hillshade.
:func:`stormscape.dem.hillshade` shades a DEM you already hold; this shades an
arbitrary *display* grid straight from the tile store, which is what a map
figure actually needs.

Hillshading differentiates elevation, so it is the loudest possible display of
a resampling mistake. Four rules, each one paid for:

1. **Warp each tile exactly once, from its native grid, with bilinear.** No
   decimating pre-pass. Reading tiles at 1/4 with ``Resampling.average`` and
   then warping cost 2.2 m elevation RMS against the single-warp path — the
   same error class this package removed on 2026-07-31 (see CLAUDE.md: a
   nearest warp inside ``py3dep.get_dem`` was biasing slope ~1.8 deg and
   halving plan curvature). firescape's ``statewide.dem_for_unit`` already
   follows this for the science grids; so does this module for the display grid.

2. **Mosaic elevations into one continuous surface before shading.** Shading
   tiles individually and mosaicking the shaded results prints a hard line at
   every tile boundary, because the gradient kernel at a tile edge has no data
   beyond it. Source windows are read with a 2-cell margin for the same
   reason, mirroring :mod:`stormscape.plot`'s "+20 cells so the kernel never
   reaches past the data".

3. **Shade finer than the display grid, then average the shaded array down.**
   Shading at display resolution throws away every landform smaller than a
   display pixel: measured Laplacian roughness 0.045 vs 0.088 at 150 m. This
   mirrors ``plot._prepare_hillshade``, which downsamples the *hillshade*
   with ``Resampling.average``.

4. **Use a fixed 0-1 shading map** (:func:`shade`), never
   ``matplotlib.colors.LightSource.hillshade``, whose contrast stretch is
   rescaled by the min/max of whatever array it is handed. Statewide DEMs must
   be shaded in bands to bound memory, and a per-array stretch would give each
   band its own contrast — printing exactly the seams rule 2 removes. Here
   flat ground is always ``sin(altitude)`` however the work is split.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

WORK_CRS = "EPSG:5070"

#: Sun position and vertical exaggeration (matches :func:`dem.hillshade`).
AZIMUTH_DEG = 315.0
ALTITUDE_DEG = 45.0
VERT_EXAG = 1.5

#: Metric grid the shading is computed on, before averaging to the display grid.
SHADE_RES_M = 50.0

#: Rows per shading band. Bounds peak memory; results are band-count invariant.
BAND_ROWS = 4000


def light_direction(azimuth: float = AZIMUTH_DEG,
                    altitude: float = ALTITUDE_DEG) -> np.ndarray:
    """Unit vector toward the light source (matplotlib's LightSource convention:
    azimuth in degrees clockwise from north, altitude above the horizon)."""
    az = np.radians(90.0 - azimuth)
    alt = np.radians(altitude)
    return np.array([np.cos(az) * np.cos(alt),
                     np.sin(az) * np.cos(alt),
                     np.sin(alt)])


def shade(elevation, res_m: float, *, direction=None,
          vert_exag: float = VERT_EXAG) -> np.ndarray:
    """Lambertian shading in [0, 1] on a metric grid, with a **fixed** mapping.

    Flat ground always maps to ``sin(altitude)`` and self-shadowed slopes clip
    at 0, independent of what else is in the array — so shading an array in
    pieces gives bit-identical results to shading it whole (given a one-row
    halo on each interior band, which :func:`shaded_relief` supplies).

    ``LightSource.hillshade`` instead rescales by the array's own min/max; it
    is deliberately not used here. See the module docstring, rule 4.
    """
    if direction is None:
        direction = light_direction()
    z = np.asarray(elevation, dtype="float32")
    # Row spacing is negative: row index increases southward.
    e_dy, e_dx = np.gradient(vert_exag * z, -res_m, res_m)
    inten = -e_dx * direction[0] - e_dy * direction[1] + direction[2]
    inten /= np.sqrt(e_dx * e_dx + e_dy * e_dy + 1.0)
    return np.clip(inten, 0.0, 1.0, out=inten)


def _elevation_mosaic(tiles, bounds5070, res_m: float):
    """Continuous float32 elevation mosaic on an EPSG:5070 grid, NaN where no
    tile covers. Each tile is warped exactly once into only the destination
    block it touches, so peak memory stays near one output array even when the
    grid is statewide.

    Kept separate from firescape's ``statewide._mosaic_to_grid``, which serves
    the science grids: that one snaps/pads a unit-sized grid and also handles
    integer class rasters, while this one honours a caller's grid and is
    destination-windowed for arrays two orders of magnitude larger.
    """
    import rasterio
    from rasterio.transform import from_origin
    from rasterio.warp import Resampling, reproject, transform_bounds
    from rasterio.windows import Window
    from rasterio.windows import from_bounds as window_from_bounds

    w5, s5, e5, n5 = bounds5070
    width = int(np.ceil((e5 - w5) / res_m))
    height = int(np.ceil((n5 - s5) / res_m))
    if width <= 0 or height <= 0:
        raise ValueError(f"empty relief grid: {width}x{height}")
    dst = np.full((height, width), np.nan, dtype="float32")

    for p in tiles:
        with rasterio.open(p) as ds:
            try:
                tb = transform_bounds(ds.crs, WORK_CRS, *ds.bounds)
            except Exception:
                continue
            if not (tb[0] < e5 and tb[2] > w5 and tb[1] < n5 and tb[3] > s5):
                continue
            c0 = max(0, int(np.floor((tb[0] - w5) / res_m)) - 2)
            c1 = min(width, int(np.ceil((tb[2] - w5) / res_m)) + 2)
            r0 = max(0, int(np.floor((n5 - tb[3]) / res_m)) - 2)
            r1 = min(height, int(np.ceil((n5 - tb[1]) / res_m)) + 2)
            if c1 <= c0 or r1 <= r0:
                continue
            sub_w, sub_n = w5 + c0 * res_m, n5 - r0 * res_m
            sub_e, sub_s = w5 + c1 * res_m, n5 - r1 * res_m

            # Read the source window covering this block, plus a 2-cell margin
            # so the bilinear kernel never reaches past what was read.
            sw, ss, se, sn = transform_bounds(WORK_CRS, ds.crs,
                                              sub_w, sub_s, sub_e, sub_n)
            win = window_from_bounds(sw, ss, se, sn, ds.transform)
            win = win.round_offsets().round_lengths()
            win = Window(win.col_off - 2, win.row_off - 2,
                         win.width + 4, win.height + 4)
            try:
                win = win.intersection(Window(0, 0, ds.width, ds.height))
            except rasterio.errors.WindowError:
                continue
            if win.width <= 0 or win.height <= 0:
                continue

            arr = ds.read(1, window=win).astype("float32")
            if ds.nodata is not None:
                arr[arr == ds.nodata] = np.nan
            piece = np.full((r1 - r0, c1 - c0), np.nan, dtype="float32")
            reproject(arr, piece,
                      src_transform=ds.window_transform(win), src_crs=ds.crs,
                      src_nodata=np.nan,
                      dst_transform=from_origin(sub_w, sub_n, res_m, res_m),
                      dst_crs=WORK_CRS, dst_nodata=np.nan,
                      resampling=Resampling.bilinear)
            view = dst[r0:r1, c0:c1]           # basic slice -> writes through
            put = np.isnan(view) & ~np.isnan(piece)
            view[put] = piece[put]
            del arr, piece
    return dst


def shaded_relief(tiles, transform, shape, *, crs: str = "EPSG:4326",
                  shade_res_m: float = SHADE_RES_M,
                  azimuth: float = AZIMUTH_DEG,
                  altitude: float = ALTITUDE_DEG,
                  vert_exag: float = VERT_EXAG,
                  band_rows: int = BAND_ROWS) -> np.ndarray:
    """Hillshade on an arbitrary display grid, from 3DEP tiles.

    Parameters
    ----------
    tiles
        Paths of the DEM tiles to consider; those missing the grid are skipped.
    transform, shape
        Destination grid, as an affine transform and ``(rows, cols)``. Typically
        a geographic grid so map axes can carry degree labels — shading still
        happens on the metric grid, so illumination stays geometrically true.
    crs
        CRS of the destination grid. Default EPSG:4326.
    shade_res_m
        Metric resolution to shade at. Should be **finer** than the display
        grid (rule 3); the shaded array is then averaged down.
    band_rows
        Rows per shading band. Results do not depend on this — it only bounds
        peak memory. Each band carries a one-row halo so interior rows see both
        neighbours.

    Returns
    -------
    float32 array of ``shape``, in [0, 1]. Ground not covered by any tile is
    set to the flat-ground value, ``sin(altitude)``, so it reads as level
    terrain rather than as a grey block or a cliff.
    """
    from rasterio.transform import from_origin
    from rasterio.warp import Resampling, reproject, transform_bounds
    from scipy.ndimage import binary_dilation

    tiles = [Path(t) for t in tiles]
    direction = light_direction(azimuth, altitude)
    flat = np.float32(direction[2])
    rows, cols = int(shape[0]), int(shape[1])

    west, north = transform.c, transform.f
    east = west + cols * transform.a
    south = north + rows * transform.e          # transform.e is negative
    w5, s5, e5, n5 = transform_bounds(crs, WORK_CRS, west, south, east, north)

    dem = _elevation_mosaic(tiles, (w5, s5, e5, n5), shade_res_m)
    height, width = dem.shape
    tr5 = from_origin(w5, n5, shade_res_m, shade_res_m)

    void = np.isnan(dem)
    if void.all():
        return np.full((rows, cols), flat, dtype="float32")
    # Any finite fill works: the plateau it creates, and the cliff at its rim,
    # are both overwritten below. The median just keeps the array debuggable.
    np.nan_to_num(dem, copy=False, nan=float(np.nanmedian(dem)))

    hs5 = np.empty((height, width), dtype="float32")
    step = max(1, int(band_rows))
    for y0 in range(0, height, step):
        y1 = min(height, y0 + step)
        top = max(0, y0 - 1)                    # one-row halo each side
        band = shade(dem[top:min(height, y1 + 1)], shade_res_m,
                     direction=direction, vert_exag=vert_exag)
        hs5[y0:y1] = band[y0 - top:y0 - top + (y1 - y0)]
    del dem

    if void.any():
        # Dilate past the kernel's reach so the rim cliff goes too.
        hs5[binary_dilation(void, np.ones((5, 5), bool))] = flat
    del void

    out = np.full((rows, cols), flat, dtype="float32")
    reproject(hs5, out, src_transform=tr5, src_crs=WORK_CRS,
              dst_transform=transform, dst_crs=crs,
              resampling=Resampling.average)
    return out
