"""USGS 3DEP / The National Map DEMs and hillshades.

Thin, reusable wrappers around :mod:`py3dep` (the HyRiver 3DEP client) for the
two things this toolkit needs:

  * download a bare-earth DEM for an AOI at a chosen resolution (1 m lidar,
    10 m, 30 m, ...), reprojected to a working CRS; and
  * render a matplotlib hillshade on the same grid to use as a map backdrop.

3DEP is the public elevation layer behind The National Map, so this fetches
the same data the interactive downloader serves, but scripted and AOI-clipped.
"""

from __future__ import annotations

import time
import warnings

import numpy as np
import rioxarray  # noqa: F401  (registers the .rio accessor on xarray objects)
import xarray as xr
from matplotlib.colors import LightSource

from .aoi import bbox_polygon, load_aoi

warnings.filterwarnings("ignore")

# Resolutions 3DEP publishes, coarse->fine. py3dep wants an int (metres).
STD_RESOLUTIONS = (60, 30, 10, 5, 3, 1)


def _py3dep():
    """Import py3dep on first use.

    py3dep pulls in the whole HyRiver stack, but it is only needed for the
    *network* 3DEP fetches in this module -- so importing it lazily keeps
    ``import stormscape`` (and everything that does not touch 3DEP: MRMS,
    gauges, smoothing, export, plotting) working without it, mirroring how
    :mod:`stormscape.nexrad` defers ``arm_pyart``/``nexradaws``.
    """
    try:
        import py3dep
    except ImportError as e:                              # pragma: no cover
        raise ImportError(
            "py3dep is required for 3DEP DEM access. Install it with\n"
            "  conda install -c conda-forge py3dep") from e
    return py3dep


def dem_sources(aoi, res=None):
    """Query 3DEP source footprints over an AOI (which lidar projects cover it).

    Returns the raw GeoDataFrame from ``py3dep.query_3dep_sources`` (columns
    include ``dem_res`` like ``'1m'``/``'10m'``). ``res`` optionally filters to
    one or more resolution strings. Useful before committing to a 1 m fetch.
    """
    bounds, _ = load_aoi(aoi)
    return _py3dep().query_3dep_sources(bounds, crs=4326, res=res)


def coverage_fraction(aoi, res="1m"):
    """Fraction of the AOI covered by 3DEP sources at resolution ``res``.

    Mirrors the 1 m-availability check used to decide whether a basin can be
    fetched at lidar resolution without seam artefacts from 10 m fill.
    """
    bounds, geom = load_aoi(aoi)
    if geom is None:
        geom = bbox_polygon(bounds)
    try:
        src = _py3dep().query_3dep_sources(bounds, crs=4326, res=res)
    except Exception:                                  # noqa: BLE001
        return 0.0
    if src is None or not len(src):
        return 0.0
    src = src.to_crs(5070)
    g = (__import__("geopandas").GeoSeries([geom], crs=4326)
         .to_crs(5070).iloc[0])
    inter = src.geometry.union_all().intersection(g)
    return float(inter.area / g.area) if g.area else 0.0


def get_dem(aoi, resolution=10, dst_crs="EPSG:5070", pad_deg=0.02,
            clip=False, drop_fill=True, retries=2, retry_wait=5):
    """Download a 3DEP DEM for an AOI and reproject it.

    Parameters
    ----------
    aoi
        Anything :func:`stormscape.aoi.load_aoi` accepts (bbox tuple, vector
        path, or shapely geometry in EPSG:4326).
    resolution
        Target resolution in metres (1, 3, 5, 10, 30, 60). 1 m requires lidar
        coverage -- check :func:`coverage_fraction` first.
    dst_crs
        CRS to reproject into. Default EPSG:5070 (CONUS Albers, equal-area,
        metres) so pixels are square metres for slope/hillshade.
    pad_deg
        Degrees of padding around the AOI bounds before fetching.
    clip
        If a clip polygon is available (vector AOI), mask outside it to NaN.
    drop_fill
        Replace 3DEP fill/no-data (< -1e4) with NaN.
    retries, retry_wait
        The 3DEP dynamic WMS can be slow/flaky (especially at 1 m, where a large
        tile can exceed py3dep's 120 s request timeout). Retry the fetch up to
        ``retries`` times with a linear backoff of ``retry_wait`` seconds.

    Returns
    -------
    xarray.DataArray
        Elevation (m) on the ``dst_crs`` grid with ``.rio`` georeferencing.
    """
    bounds, geom = load_aoi(aoi, pad_deg=pad_deg)
    for attempt in range(retries + 1):
        try:
            dem = _py3dep().get_dem(bbox_polygon(bounds),
                                    resolution=resolution, crs=4326)
            break
        except Exception:                      # noqa: BLE001  (WMS timeout, etc.)
            if attempt >= retries:
                raise
            time.sleep(retry_wait * (attempt + 1))
    if dst_crs is not None:
        dem = dem.rio.reproject(dst_crs, resolution=resolution)
    if drop_fill:
        dem = dem.where(dem > -1e4)
    if clip and geom is not None:
        import geopandas as gpd
        gseries = gpd.GeoSeries([geom], crs=4326).to_crs(dem.rio.crs)
        dem = dem.rio.clip(gseries.geometry, drop=True)
    dem.name = "elevation"
    return dem


def hillshade(dem: xr.DataArray, azimuth=315, altitude=45, vert_exag=1.5):
    """Matplotlib hillshade (0-255 float) on the DEM's own grid.

    Pixel spacing is read from the DEM georeferencing, so this is correct for
    any projected CRS/resolution. NaNs (outside-AOI / fill) are filled to the
    DEM minimum before shading and remain visually flat.
    """
    z = dem.values.astype("float64")
    if z.ndim == 3:                                    # (band, y, x) -> (y, x)
        z = z[0]
    dx, dy = dem.rio.resolution()
    zmin = np.nanmin(z)
    hs = LightSource(azdeg=azimuth, altdeg=altitude).hillshade(
        np.nan_to_num(z, nan=zmin), vert_exag=vert_exag,
        dx=abs(dx), dy=abs(dy))
    hs = (hs * 255.0).astype("float32")
    out = dem.squeeze().copy(data=hs)
    out.name = "hillshade"
    return out


def fetch_dem_and_hillshade(aoi, resolution=10, dst_crs="EPSG:5070",
                            dem_path=None, hillshade_path=None, pad_deg=0.02,
                            azimuth=315, altitude=45, vert_exag=1.5,
                            clip=False, retries=2, retry_wait=5):
    """Convenience: download a DEM, build its hillshade, optionally write both.

    Returns ``(dem, hillshade)`` DataArrays. GeoTIFFs are LZW-compressed.
    ``retries``/``retry_wait`` are forwarded to :func:`get_dem` (3DEP can be
    slow at 1 m).
    """
    dem = get_dem(aoi, resolution=resolution, dst_crs=dst_crs, pad_deg=pad_deg,
                  clip=clip, retries=retries, retry_wait=retry_wait)
    hs = hillshade(dem, azimuth=azimuth, altitude=altitude,
                   vert_exag=vert_exag)
    if dem_path:
        dem.rio.to_raster(dem_path, compress="LZW")
    if hillshade_path:
        hs.rio.to_raster(hillshade_path, compress="LZW")
    return dem, hs
