"""stormscape.export -- georeferenced exports for GIS / CalTopo.

Two products, both derived from an already-processed event directory:

1. **EPSG:3857 GeoTIFFs** of the rainfall fields (Web-Mercator, CalTopo's native
   projection). For each field both a **raw single-band float** GeoTIFF (the data
   values, for analysis / restyling) and a **colorized RGBA** GeoTIFF (the
   project colormap with transparent dry / no-data cells, so it draws correctly
   as a CalTopo overlay without the viewer having to style a float raster).

2. **Georeferenced PDFs** of the primary map figures (the i15 map and the
   anomaly map). The full styled matplotlib figure is rasterized and wrapped as a
   **GeoPDF**: the map *axes* are registered to their projected coordinates and a
   **neatline** bounds the georeferenced area to the map frame, so the colorbar /
   title / margins are not mis-georeferenced. Readable as a georeferenced layer
   by CalTopo / Avenza / QGIS / ArcGIS.

The GeoPDF path needs GDAL's PDF driver, a separate conda-forge plugin::

    conda install -c conda-forge libgdal-pdf

It is imported lazily, so ``import stormscape`` and the GeoTIFF export work
without it; only :func:`figure_to_geopdf` requires it. The GeoTIFF export needs
only rioxarray (already a core dependency).
"""
from __future__ import annotations

import os
from typing import Optional, Sequence, Union

import numpy as np

from .layout import find, out_path

WEB_MERCATOR = "EPSG:3857"

# field suffixes (after "<key>_") exported to EPSG:3857 by default -- the
# anomaly and the peak-i15 field, the two layers a CalTopo overlay wants.
DEFAULT_EXPORT_FIELDS = ("anom_i15", "i15max")

# field-name substrings that are categorical / quality flags: reproject by
# nearest neighbour and never colour-ramp them like a continuous field.
_CATEGORICAL = ("tpki15", "rqi", "shsr", "cbb", "flag", "severity")


def _is_categorical(field: str) -> bool:
    f = str(field).lower()
    return any(t in f for t in _CATEGORICAL)


def _crs_tag(dst_crs) -> str:
    """Short filename tag for a CRS (``EPSG:3857`` -> ``3857``)."""
    s = str(dst_crs).upper().replace("EPSG:", "").strip()
    return s if s.isdigit() else "proj"


def _open_rio(path):
    """Open a raster as a masked (NaN-nodata) DataArray, loaded + closed promptly
    (the lazy handle otherwise lingers to GC and trips a benign excepthook at
    interpreter exit)."""
    import rioxarray
    with rioxarray.open_rasterio(path, masked=True) as da:
        return da.load()


# --------------------------------------------------------------------------- #
# 1) EPSG:3857 GeoTIFF export (raw float + colorized RGBA)
# --------------------------------------------------------------------------- #
def reproject_geotiff(in_tif: str, out_tif: str, dst_crs: str = WEB_MERCATOR,
                      resampling=None, compress: str = "deflate") -> str:
    """Reproject one GeoTIFF to ``dst_crs`` (default EPSG:3857), preserving the
    NaN nodata. ``resampling`` defaults to **bilinear** for continuous fields and
    **nearest** for categorical ones (inferred from the filename); pass a
    :class:`rasterio.enums.Resampling` or its name to override."""
    import rioxarray  # noqa: F401  (registers the .rio accessor)
    from rasterio.enums import Resampling
    if resampling is None:
        resampling = (Resampling.nearest if _is_categorical(os.path.basename(in_tif))
                      else Resampling.bilinear)
    elif isinstance(resampling, str):
        resampling = getattr(Resampling, resampling)
    da = _open_rio(in_tif)
    out = da.rio.reproject(dst_crs, resampling=resampling)
    os.makedirs(os.path.dirname(os.path.abspath(out_tif)), exist_ok=True)
    out.rio.to_raster(out_tif, compress=compress)
    return out_tif


def _field_style(field: str, arr, *, cmap: str = "YlGnBu",
                 anomaly_cmap: str = "cmc.vik", wet_min: float = 5.0,
                 vmax: Optional[float] = None):
    """``(cmap, vmin, vmax, norm, mask_below)`` for colorizing ``field``, matching
    the figures: anomalies use a diverging map centred at 1.0 (``TwoSlopeNorm``);
    intensity / depth fields use the sequential ``cmap`` with the dry tail masked.
    ``arr`` supplies the data for the 99th-percentile autoscale."""
    from matplotlib.colors import TwoSlopeNorm
    finite = arr[np.isfinite(arr)]
    f = str(field).lower()
    if f.startswith("anom"):
        vmx = vmax if vmax is not None else (
            max(2.0, float(np.percentile(finite, 99))) if finite.size else 2.0)
        return anomaly_cmap, 0.0, vmx, TwoSlopeNorm(vcenter=1.0, vmin=0.0,
                                                    vmax=vmx), None
    if f.startswith("dnbr") or "severity" in f:
        # burn severity runs 0-1 (dNBR) or 0-N (classes), not 0-100+ like an
        # intensity, so it needs its own scale; below the unburned break stays
        # transparent so an unburned hillside shows the basemap, not pale paint
        if "severity" in f:
            vmx = vmax if vmax is not None else (
                float(np.nanmax(finite)) if finite.size else 4.0)
            return "YlOrRd", 0.0, max(vmx, 1.0), None, 0.5
        vmx = vmax if vmax is not None else 1.0
        return "YlOrRd", 0.0, vmx, None, 0.10
    vmx = vmax if vmax is not None else (
        max(10.0, float(np.percentile(finite, 99))) if finite.size else 10.0)
    intensity = any(t in f for t in ("i15", "i30", "i60", "i2max", "peakrate"))
    return cmap, 0.0, vmx, None, (wet_min if intensity else None)


def _colormap_rgba(arr, cmap, vmin, vmax, norm=None, mask_below=None):
    """Map a 2-D float array to an ``(H, W, 4)`` uint8 RGBA image via a matplotlib
    colormap. NaN cells -- and, if given, cells below ``mask_below`` -- are fully
    transparent (alpha 0) so dry / no-data areas let the CalTopo basemap show."""
    import matplotlib
    from matplotlib.colors import Normalize
    if isinstance(cmap, str) and cmap.startswith("cmc."):
        try:
            import cmcrameri.cm  # noqa: F401  registers the cmc.* colormaps
        except ImportError:
            cmap = "RdBu_r"
    cmap_obj = matplotlib.colormaps[cmap] if isinstance(cmap, str) else cmap
    if norm is None:
        norm = Normalize(vmin=vmin, vmax=vmax)
    a = np.asarray(arr, dtype=float)
    n = np.ma.filled(norm(np.ma.masked_invalid(a)).astype(float), np.nan)
    # Normalizers may return +/-inf for values outside [vmin, vmax] (TwoSlopeNorm
    # does), and NaN for the masked cells. Resolve both explicitly and clip to the
    # unit interval so out-of-range values take the end colours -- what imshow
    # shows -- instead of overflowing inside the colormap lookup.
    n = np.clip(np.nan_to_num(n, nan=0.0, posinf=1.0, neginf=0.0), 0.0, 1.0)
    rgba = (cmap_obj(n) * 255).astype(np.uint8)
    alpha = np.where(np.isfinite(a), 255, 0)
    if mask_below is not None:
        alpha = np.where(a < float(mask_below), 0, alpha)
    rgba[..., 3] = alpha.astype(np.uint8)
    return rgba


def write_rgba(out_tif: str, rgba, transform, crs) -> str:
    """Write an ``(H, W, 4)`` uint8 RGBA array as a 4-band GeoTIFF with the alpha
    band tagged, so viewers honour the transparency.

    Tagging the fourth band ``ColorInterp.alpha`` is the whole point: without it
    CalTopo paints the layer's full rectangle over the basemap instead of
    letting the untouched cells through.
    """
    import rasterio
    from rasterio.enums import ColorInterp
    h, w = rgba.shape[:2]
    prof = dict(driver="GTiff", width=w, height=h, count=4, dtype="uint8",
                crs=crs, transform=transform, compress="deflate")
    os.makedirs(os.path.dirname(os.path.abspath(out_tif)), exist_ok=True)
    with rasterio.open(out_tif, "w", **prof) as dst:
        for b in range(4):
            dst.write(rgba[:, :, b], b + 1)
        dst.colorinterp = [ColorInterp.red, ColorInterp.green,
                           ColorInterp.blue, ColorInterp.alpha]
    return out_tif


def export_geotiffs(in_dir: str, key: str, out_dir: str,
                    fields: Sequence[str] = DEFAULT_EXPORT_FIELDS,
                    dst_crs: str = WEB_MERCATOR, out_key: Optional[str] = None,
                    raw: bool = True, colorize: bool = True,
                    cmap: str = "YlGnBu", anomaly_cmap: str = "cmc.vik",
                    wet_min: float = 5.0, vmax: Optional[float] = None) -> list:
    """Reproject an event's fields to ``dst_crs`` (default EPSG:3857) for GIS /
    CalTopo.

    For each ``<in_dir>/<key>_<field>.tif`` writes, under ``out_dir``:

      * ``<out_key>_<field>_<tag>.tif`` -- the **raw** single-band float field
        (when ``raw``), and
      * ``<out_key>_<field>_<tag>_rgb.tif`` -- a **colorized RGBA** image with
        transparent dry / no-data cells (when ``colorize`` and the field is
        continuous), styled like the figures.

    ``<tag>`` is the EPSG code (``3857``). Missing fields are skipped with a note.
    Returns the list of written paths.
    """
    out_key = out_key or key
    tag = _crs_tag(dst_crs)
    os.makedirs(out_dir, exist_ok=True)
    written = []
    for field in fields:
        src = find(in_dir, f"{key}_{field}.tif")
        if not os.path.exists(src):
            print(f"note: {src} not found; skipping {field}")
            continue
        raw_path = out_path(out_dir, f"{out_key}_{field}_{tag}.tif")
        reproject_geotiff(src, raw_path, dst_crs=dst_crs)
        if colorize and not _is_categorical(field):
            da = _open_rio(raw_path)
            arr = np.asarray(da.squeeze().values, dtype=float)
            sty, vmin, vmx, norm, mask_below = _field_style(
                field, arr, cmap=cmap, anomaly_cmap=anomaly_cmap,
                wet_min=wet_min, vmax=vmax)
            rgba = _colormap_rgba(arr, sty, vmin, vmx, norm=norm,
                                  mask_below=mask_below)
            rgb_path = out_path(out_dir, f"{out_key}_{field}_{tag}_rgb.tif")
            write_rgba(rgb_path, rgba, da.rio.transform(), da.rio.crs)
            written.append(rgb_path)
        if raw:
            written.append(raw_path)
        else:
            os.remove(raw_path)              # only the colorized image was wanted
    return written


# --------------------------------------------------------------------------- #
# 2) full-resolution stream network as a vector layer
# --------------------------------------------------------------------------- #
# OGR write driver per file extension (all present in the conda-forge GDAL).
_VECTOR_DRIVERS = {".geojson": "GeoJSON", ".json": "GeoJSON", ".gpkg": "GPKG",
                   ".shp": "ESRI Shapefile", ".kml": "KML", ".gpx": "GPX"}


def _vector_driver(path: str) -> str:
    return _VECTOR_DRIVERS.get(os.path.splitext(path)[1].lower(), "GeoJSON")


def export_streams(aoi, out_path: str, *, named_only: bool = False,
                   watercourse_only: bool = True, clip_to_aoi: bool = True,
                   driver: Optional[str] = None):
    """Export the **full-resolution NHD stream network** for an AOI as a vector
    layer (default GeoJSON, EPSG:4326 -- the CalTopo-native vector CRS).

    Fetches USGS **NHDPlus HR** flowlines for the AOI (the same dense network the
    figures' reference overlay draws), keeping the *whole* network by default --
    named creeks/rivers **and** unnamed headwater rills. ``clip_to_aoi`` truncates
    the flowlines to the AOI polygon when the AOI is a polygon (a ``--bbox`` AOI
    has no polygon, so the bounding-box network is kept either way); set it False
    to keep whole flowlines over the AOI bounding box. ``named_only`` keeps only
    GNIS-named flowlines; ``watercourse_only`` drops canals / ditches / pipelines.
    ``driver`` is inferred from the file extension (``.geojson`` / ``.gpkg`` /
    ``.shp`` / ``.kml``) unless given. Returns the path, or None if the network
    fetch returned nothing.
    """
    import geopandas as gpd

    from . import refdata
    from .aoi import load_aoi
    g = refdata.streams(aoi, named_only=named_only,
                        watercourse_only=watercourse_only)
    if not len(g):
        print("note: no NHD flowlines returned for the AOI; nothing written")
        return None
    if clip_to_aoi:
        _, geom = load_aoi(aoi)
        if geom is not None:
            g = gpd.clip(g, gpd.GeoDataFrame(geometry=[geom], crs=4326))
            g = g[~g.geometry.is_empty & g.geometry.notna()].reset_index(drop=True)
    keep = [c for c in ("name", "fcode") if c in g.columns]   # tidy attributes
    g = g[keep + ["geometry"]]
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    g.to_file(out_path, driver=driver or _vector_driver(out_path))
    return out_path


# --------------------------------------------------------------------------- #
# 3) georeferenced PDF of a map figure (GeoPDF)
# --------------------------------------------------------------------------- #
def geopdf_supported() -> bool:
    """True if GDAL's PDF write driver is available (needs ``libgdal-pdf``)."""
    try:
        from osgeo import gdal
        return gdal.GetDriverByName("PDF") is not None
    except Exception:                          # noqa: BLE001 (osgeo may be absent)
        return False


def figure_to_geopdf(fig, ax, out_pdf: str, crs, *, dpi: int = 200,
                     neatline: bool = True, author: str = "stormscape",
                     title: Optional[str] = None, subject: Optional[str] = None,
                     geo_encoding: str = "ISO32000") -> str:
    """Write a matplotlib map figure as a **georeferenced PDF**.

    The whole styled figure is rasterized at ``dpi`` and embedded; the **map
    axes** (``ax``) are registered to their projected coordinates. ``crs`` must be
    the CRS the axes are plotted in (i.e. the figure's ``work_crs``) -- the
    data<->page mapping over the axes is then exactly affine. A **neatline**
    (default) limits the georeferenced area to the map frame, so the colorbar /
    title margins are excluded rather than mis-georeferenced.

    Needs GDAL's PDF driver (``conda install -c conda-forge
    libgdal-pdf``). ``geo_encoding`` is ``ISO32000`` (the modern Adobe geospatial
    standard, default) or ``OGC_BP`` (the legacy TerraGo/Avenza flavour).
    """
    from io import BytesIO

    from osgeo import gdal, osr
    gdal.UseExceptions()
    if gdal.GetDriverByName("PDF") is None:
        raise RuntimeError(
            "GDAL has no PDF driver; install it with\n"
            "  conda install -c conda-forge libgdal-pdf")
    try:
        from PIL import Image                  # ships with matplotlib
    except ImportError as e:                   # pragma: no cover
        raise RuntimeError("Pillow is required to rasterize the figure") from e

    # finalize layout, then read the map-axes rectangle (figure fractions, y up)
    fig.canvas.draw()
    pos = ax.get_position()

    # rasterize the full page at dpi (NO bbox_inches='tight', so the axes
    # fraction maps deterministically onto the page pixels)
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=dpi)
    buf.seek(0)
    img = np.asarray(Image.open(buf).convert("RGB"))      # (H, W, 3) uint8
    H, W = img.shape[:2]

    # map-axes rectangle in image pixels (image rows top-down; fig y is bottom-up)
    left, right = pos.x0 * W, pos.x1 * W
    top, bottom = (1.0 - pos.y1) * H, (1.0 - pos.y0) * H

    # axes data limits in `crs`: (left,right) spine x and (bottom,top) spine y --
    # this spine convention is exact even if an axis is inverted.
    xl0, xl1 = ax.get_xlim()
    yl0, yl1 = ax.get_ylim()
    dx = (xl1 - xl0) / (right - left)
    dy = (yl1 - yl0) / (top - bottom)          # top<bottom in rows -> dy<0 north-up
    gt = (xl0 - left * dx, dx, 0.0, yl1 - top * dy, 0.0, dy)

    mem = gdal.GetDriverByName("MEM").Create("", W, H, 3, gdal.GDT_Byte)
    for b in range(3):
        mem.GetRasterBand(b + 1).WriteArray(np.ascontiguousarray(img[:, :, b]))
    mem.SetGeoTransform(gt)
    srs = osr.SpatialReference()
    srs.SetFromUserInput(str(crs))
    srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    mem.SetProjection(srs.ExportToWkt())

    opts = [f"DPI={int(round(dpi))}", f"GEO_ENCODING={geo_encoding}",
            "MARGIN=0", f"AUTHOR={author}", "CREATOR=stormscape"]
    if title:
        opts.append(f"TITLE={title}")
    if subject:
        opts.append(f"SUBJECT={subject}")
    if neatline:
        gx0, gx1 = sorted((xl0, xl1))
        gy0, gy1 = sorted((yl0, yl1))
        opts.append("NEATLINE=POLYGON(({0} {2},{0} {3},{1} {3},{1} {2},{0} {2}))"
                    .format(gx0, gx1, gy0, gy1))
    os.makedirs(os.path.dirname(os.path.abspath(out_pdf)), exist_ok=True)
    out = gdal.GetDriverByName("PDF").CreateCopy(out_pdf, mem, options=opts)
    out.FlushCache()
    out = None
    mem = None
    return out_pdf
