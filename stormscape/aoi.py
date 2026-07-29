"""Area-of-interest helpers: turn anything into a lon/lat bounding box.

The DEM and MRMS engines only need an AOI expressed as a geographic
(EPSG:4326) bounding box ``(west, south, east, north)`` plus, optionally, a
clip polygon. This module accepts the three things a user is likely to have
on hand -- a bbox tuple, a shapely geometry, or a vector file (GeoJSON /
Shapefile / GeoPackage) -- and normalises them, so nothing downstream is
tied to one project's data layout.
"""

from __future__ import annotations

import os
from typing import Iterable

import geopandas as gpd
from shapely.geometry import box
from shapely.geometry.base import BaseGeometry


def pad_bounds(bounds: tuple, pad_deg: float) -> tuple:
    """Expand a (W, S, E, N) bbox by ``pad_deg`` degrees on every side."""
    w, s, e, n = bounds
    return (w - pad_deg, s - pad_deg, e + pad_deg, n + pad_deg)


def load_aoi(spec, layer=None, pad_deg: float = 0.0):
    """Normalise an AOI specification to ``(bounds, geometry)`` in EPSG:4326.

    Parameters
    ----------
    spec
        One of:
          * a 4-tuple/list ``(west, south, east, north)`` in lon/lat degrees;
          * a path to a vector file readable by GeoPandas (``.geojson``,
            ``.shp``, ``.gpkg`` ...); all features are dissolved into one
            geometry and reprojected to 4326;
          * a shapely geometry (assumed already in EPSG:4326).
    layer
        Optional layer name for multi-layer sources (e.g. GeoPackage).
    pad_deg
        Degrees to pad the returned bounds on each side (the clip geometry is
        not padded). Use a small pad (~0.02-0.05) so edge derivatives such as
        slope/hillshade and the rolling i15 window are valid to the border.

    Returns
    -------
    (bounds, geometry)
        ``bounds`` is the padded ``(W, S, E, N)`` tuple; ``geometry`` is the
        unpadded clip polygon in EPSG:4326, or ``None`` for a bare bbox.
    """
    geom = None
    if isinstance(spec, (tuple, list)) and len(spec) == 4 and all(
            isinstance(v, (int, float)) for v in spec):
        bounds = tuple(float(v) for v in spec)
    elif isinstance(spec, BaseGeometry):
        geom = spec
        bounds = geom.bounds
    elif isinstance(spec, str) and os.path.exists(spec):
        gdf = gpd.read_file(spec, layer=layer)
        if gdf.crs is None:
            raise ValueError(f"{spec} has no CRS; cannot place it on Earth.")
        gdf = gdf.to_crs(4326)
        geom = gdf.geometry.union_all()
        bounds = tuple(gdf.total_bounds)
    else:
        raise TypeError(
            "aoi must be a (W,S,E,N) tuple, a shapely geometry, or a path to "
            f"a vector file; got {spec!r}")
    return pad_bounds(bounds, pad_deg), geom


def bbox_polygon(bounds: tuple) -> BaseGeometry:
    """Shapely box for a (W, S, E, N) bbox (EPSG:4326)."""
    return box(*bounds)


def read_overlay(spec, dst_crs="EPSG:5070", layer=None):
    """Read a vector overlay (perimeter, basins, points) for plotting.

    Accepts a path, a GeoDataFrame, or ``None`` (passthrough). Reprojects to
    ``dst_crs`` so it lines up with a projected hillshade.
    """
    if spec is None:
        return None
    if isinstance(spec, gpd.GeoDataFrame):
        gdf = spec
    else:
        gdf = gpd.read_file(spec, layer=layer)
    if gdf.crs is not None and dst_crs is not None:
        gdf = gdf.to_crs(dst_crs)
    return gdf
