"""Drape the i15 field over a hillshade and add optional vector overlays.

Reproduces the storm-day map style from the debris-flow project, but with
every project-specific layer (fire perimeter, basins, gauges) reduced to a
generic, optional overlay so the same figure works for any AOI.
"""

from __future__ import annotations

import os
import warnings

import numpy as np
import pandas as pd
import rioxarray  # noqa: F401  (registers .rio accessor)
import xarray as xr
from matplotlib.lines import Line2D

from .aoi import read_overlay
from .layout import find

warnings.filterwarnings("ignore")

# Project-wide default opacity for the rainfall/field layer draped over the
# hillshade. Every map function resolves ``alpha=None`` to this value, so the
# whole toolkit (and CLI) shares one default; pass an explicit ``alpha`` (or the
# ``--alpha`` CLI flag) to override on any single figure.
DEFAULT_FIELD_ALPHA = 0.32


def _load(da_or_path):
    """Accept an xarray DataArray or a raster path -> squeezed DataArray."""
    if isinstance(da_or_path, xr.DataArray):
        return da_or_path.squeeze()
    return rioxarray.open_rasterio(da_or_path, masked=True).squeeze()


def _to_crs(da, crs):
    if crs is None or da.rio.crs is None:
        return da
    if da.rio.crs.to_string() == str(crs) or (
            da.rio.crs.to_epsg() is not None
            and str(crs).upper().endswith(str(da.rio.crs.to_epsg()))):
        return da
    return da.rio.reproject(crs)


def _extent(da):
    b = da.rio.bounds()
    return [b[0], b[2], b[1], b[3]]      # left, right, bottom, top (imshow)


def _resolve_provider(name):
    import contextily as cx
    src = cx.providers
    for part in str(name).split("."):
        src = src[part]
    return src


def _add_basemap(ax, crs, provider, zoom, labels=None):
    """Underlay an open-source contextily basemap (roads/rivers/place names).

    The base tiles are drawn at zorder 0 (beneath the semi-transparent
    hillshade and i15 drape). If ``labels`` is given, a transparent
    labels-only layer is drawn on top (zorder 7) so place / road / river
    names stay legible even where the i15 field covers the map. contextily
    is an optional dependency; a helpful error is raised if it is missing.
    """
    try:
        import contextily as cx
    except ImportError as e:                       # pragma: no cover
        raise ImportError(
            "basemap=True needs contextily: "
            "conda install -c conda-forge contextily") from e
    # extra kwargs (zorder) are collected by add_basemap and forwarded to
    # imshow, so the tiles land on the right layer.
    zk = {"zoom": zoom} if zoom is not None else {}
    cx.add_basemap(ax, source=_resolve_provider(provider), crs=str(crs),
                   attribution_size=6, zorder=0, **zk)
    if labels:
        cx.add_basemap(ax, source=_resolve_provider(labels), crs=str(crs),
                       attribution=False, zorder=7, **zk)


def _bounds_4326(da):
    """Lon/lat (W,S,E,N) bbox of a georeferenced DataArray."""
    import geopandas as gpd
    from shapely.geometry import box
    b = da.rio.bounds()
    return tuple(gpd.GeoSeries([box(*b)], crs=da.rio.crs)
                 .to_crs(4326).total_bounds)


def _label_line(ax, geom, text, color, fontsize, italic=False):
    """Place a halo-stroked label along the middle of a line, rotated to it."""
    import math

    import matplotlib.patheffects as pe
    if geom is None or geom.is_empty:
        return
    if geom.geom_type == "MultiLineString":
        geom = max(geom.geoms, key=lambda g: g.length)
    if geom.length == 0:
        return
    mid = geom.interpolate(0.5, normalized=True)
    a = geom.interpolate(0.42, normalized=True)
    b = geom.interpolate(0.58, normalized=True)
    ang = math.degrees(math.atan2(b.y - a.y, b.x - a.x)) % 180
    if ang > 90:
        ang -= 180
    t = ax.text(mid.x, mid.y, text, fontsize=fontsize, color=color,
                style="italic" if italic else "normal", ha="center",
                va="center", rotation=ang, rotation_mode="anchor", zorder=9,
                clip_on=True)
    t.set_path_effects([pe.withStroke(linewidth=2.0, foreground="white")])


def _label_point(ax, x, y, text, fontsize):
    import matplotlib.patheffects as pe
    t = ax.text(x, y, "  " + text, fontsize=fontsize, color="black",
                va="center", ha="left", zorder=9, clip_on=True)
    t.set_path_effects([pe.withStroke(linewidth=1.6, foreground="white")])


def _longest_per_name(gdf):
    """One row per unique name (the longest geometry), for label de-cluttering."""
    g = gdf[gdf.name.notna() & (gdf.name.astype(str).str.len() > 0)].copy()
    if not len(g):
        return g
    g["_len"] = g.geometry.length
    return g.sort_values("_len", ascending=False).drop_duplicates("name")


def add_reference(ax, work_crs, streams=None, roads=None, places=None,
                  label=True, max_stream_labels=22, max_road_labels=16,
                  max_place_labels=24):
    """Draw labelled stream / road / place overlays (reprojected to work_crs).

    Each argument is a GeoDataFrame or vector path (see
    :mod:`stormscape.refdata` for AOI fetchers). Returns legend handles.
    """
    from matplotlib.lines import Line2D

    from .aoi import read_overlay
    handles = []

    s = read_overlay(streams, work_crs)
    if s is not None and len(s):
        if "name" in s.columns:
            named = s[s.name.notna() & (s.name.astype(str).str.len() > 0)]
            minor = s[~s.index.isin(named.index)]
        else:
            named, minor = s, s.iloc[0:0]
        if len(minor):                  # unnamed headwaters: finer + lighter
            minor.plot(ax=ax, color="#9ecae1", lw=0.3, alpha=0.6, zorder=2.4)
        if len(named):
            named.plot(ax=ax, color="#2c7fb8", lw=0.6, alpha=0.85, zorder=2.5)
        handles.append(Line2D([], [], color="#2c7fb8", lw=1.2,
                              label="stream (NHD)"))
        if label:
            for _, row in _longest_per_name(named).head(max_stream_labels).iterrows():
                _label_line(ax, row.geometry, row["name"], "#1b4f72", 5.5,
                            italic=True)

    r = read_overlay(roads, work_crs)
    if r is not None and len(r):
        if "kind" in r.columns:
            prim, sec = r[r.kind == "primary"], r[r.kind == "secondary"]
            loc = r[r.kind == "local"]
        else:
            prim, sec, loc = r, r.iloc[0:0], r.iloc[0:0]
        if len(loc):
            loc.plot(ax=ax, color="0.55", lw=0.5, alpha=0.7, zorder=2.6)
        if len(sec):
            sec.plot(ax=ax, color="0.2", lw=1.0, zorder=2.7)
        if len(prim):
            prim.plot(ax=ax, color="#9b2d6f", lw=1.8, zorder=2.8)
        handles.append(Line2D([], [], color="0.2", lw=1.5,
                              label="road (TIGER)"))
        if label:
            for _, row in _longest_per_name(r).head(max_road_labels).iterrows():
                _label_line(ax, row.geometry, row["name"], "black", 6.5)

    p = read_overlay(places, work_crs)
    if p is not None and len(p):
        p = p.head(max_place_labels).reset_index(drop=True)
        # GNIS mixes points (populated places) and polygons (civil places);
        # representative_point() gives a label anchor for any geometry type.
        reps = p.geometry.representative_point()
        ax.scatter(reps.x, reps.y, s=12, color="black",
                   edgecolor="white", linewidth=0.4, zorder=8.5)
        handles.append(Line2D([], [], marker="o", color="black", mec="white",
                              lw=0, markersize=5, label="place (GNIS)"))
        if label:
            for i, row in p.iterrows():
                if row.get("name"):
                    _label_point(ax, reps.iloc[i].x, reps.iloc[i].y,
                                 row["name"], 5.5)
    return handles


def add_gauges(ax, work_crs, gauges, value=None, cmap="turbo", vmin=None,
               vmax=None, size=120, label="gauge", marker="o"):
    """Overlay rain gauges, optionally coloured by a value/residual column.

    ``gauges`` is a GeoDataFrame or vector path (e.g. from
    :func:`stormscape.gauges.gauge_fields` or
    :func:`stormscape.compare.radar_vs_gauge`); it is reprojected to
    ``work_crs``. With ``value=None`` all gauges are a single-colour marker.
    Naming a column colours the finite-valued gauges by it (the returned
    scatter lets the caller attach a colour bar); gauges with no value are
    drawn as small open grey markers so missing data stays visible.

    Returns ``(legend_handles, scatter_or_None)``.
    """
    from matplotlib.lines import Line2D

    from .aoi import read_overlay
    g = read_overlay(gauges, work_crs)
    handles = []
    if g is None or not len(g):
        return handles, None
    xs, ys = g.geometry.x.to_numpy(), g.geometry.y.to_numpy()
    if value is None or value not in g.columns:
        ax.scatter(xs, ys, marker=marker, s=size, facecolor="yellow",
                   edgecolor="black", linewidth=1.2, zorder=6.5)
        handles.append(Line2D([], [], marker=marker, color="yellow",
                              mec="black", lw=0, markersize=10, label=label))
        return handles, None
    vals = pd.to_numeric(g[value], errors="coerce").to_numpy(dtype="float64")
    finite = np.isfinite(vals)
    if (~finite).any():                          # gauges with no value: open grey
        ax.scatter(xs[~finite], ys[~finite], marker=marker, s=size * 0.55,
                   facecolor="none", edgecolor="0.5", linewidth=1.0, zorder=6.4)
    sc = None
    if finite.any():
        sc = ax.scatter(xs[finite], ys[finite], c=vals[finite], marker=marker,
                        s=size, cmap=cmap, vmin=vmin, vmax=vmax,
                        edgecolor="white", linewidth=1.6, zorder=6.6)
        handles.append(Line2D([], [], marker=marker, color="0.3", mec="white",
                              lw=0, markersize=10, label=label))
    return handles, sc


def _add_north_arrow(ax, x=0.07, y=0.90):
    """A simple up-pointing north arrow with an 'N' label (axes fraction)."""
    import matplotlib.patheffects as pe
    ann = ax.annotate("", xy=(x, y + 0.06), xytext=(x, y - 0.06),
                      xycoords="axes fraction",
                      arrowprops=dict(arrowstyle="-|>", lw=2.2, color="black"))
    if ann.arrow_patch is not None:                  # white halo for contrast
        ann.arrow_patch.set_path_effects(
            [pe.withStroke(linewidth=4, foreground="white")])
    t = ax.text(x, y + 0.075, "N", transform=ax.transAxes, ha="center",
                va="bottom", fontsize=12, fontweight="bold", color="black")
    t.set_path_effects([pe.withStroke(linewidth=2.5, foreground="white")])


def _utm_crs_for(da):
    """UTM EPSG (e.g. 'EPSG:32611') for the centre of a georeferenced array."""
    import geopandas as gpd
    from shapely.geometry import box
    b = da.rio.bounds()
    cen = gpd.GeoSeries([box(b[0], b[1], b[2], b[3])], crs=da.rio.crs) \
        .to_crs(4326).iloc[0].centroid
    zone = int((cen.x + 180) // 6) + 1
    return f"EPSG:{(32600 if cen.y >= 0 else 32700) + zone}"


def _prepare_hillshade(hillshade, work_crs="UTM", max_px=2600):
    """Load, reproject, and downsample a hillshade **once** for reuse across the
    figures of a command, so the expensive prep happens a single time instead of
    per figure.

    The 1 m zoom hillshades are ~170 M cells -- ~25x more than a 200-300 dpi
    figure can resolve -- so rendering them raw is slow and memory-heavy (a
    ~170 M-cell float array is ~1.4 GB, and the rainfall field gets upsampled to
    that same grid to overlay it). This reprojects to the working CRS and
    coarsens to at most ``max_px`` cells on the long side -- an invisible loss at
    figure resolution. Returns ``(hillshade_da, resolved_work_crs)``: pass *both*
    to the figure functions (``drape_i15(hillshade=hs_da, work_crs=wc, ...)``)
    so each function's ``_to_crs`` is a no-op and the field upsamples to the small
    grid. ``max_px=None`` keeps full resolution. The on-disk GeoTIFF is never
    modified -- this only affects what gets rendered.
    """
    from rasterio.enums import Resampling
    from rasterio.warp import calculate_default_transform
    hs0 = _load(hillshade)
    wc = (_utm_crs_for(hs0) if str(work_crs).strip().upper() in ("UTM", "AUTO")
          else work_crs)
    if not max_px:
        return _to_crs(hs0, wc), wc
    left, bottom, right, top = hs0.rio.bounds()
    _, w, h = calculate_default_transform(
        hs0.rio.crs, wc, hs0.rio.width, hs0.rio.height, left, bottom, right, top)
    if max(w, h) <= max_px:                       # already at/under target -> just reproject
        return _to_crs(hs0, wc), wc
    scale = float(max_px) / max(w, h)
    shape = (max(1, round(h * scale)), max(1, round(w * scale)))  # (height, width)
    # warp straight to the downsampled grid (GDAL reads+averages in one pass --
    # never materializes the full-res reprojection); ``average`` anti-aliases.
    return hs0.rio.reproject(wc, shape=shape,
                             resampling=Resampling.average), wc


def _add_scale_ticks(ax, work_crs):
    """Latitude/longitude tick marks, placed through the figure CRS.

    Ticks sit at round lon/lat values, positioned at the projected coordinate of
    each meridian/parallel at the view centre. In a near-conformal local CRS
    (UTM) meridian convergence is small, so the labels are accurate to far
    better than their 0.01-degree precision; in a strongly-convergent CRS
    (e.g. CONUS Albers far from its central meridian) prefer plotting in UTM.
    """
    import math

    from pyproj import Transformer
    fwd = Transformer.from_crs("EPSG:4326", work_crs, always_xy=True)
    inv = Transformer.from_crs(work_crs, "EPSG:4326", always_xy=True)
    x0, x1 = sorted(ax.get_xlim())
    y0, y1 = sorted(ax.get_ylim())
    lls = [inv.transform(xx, yy) for xx in (x0, x1) for yy in (y0, y1)]
    lons, lats = [p[0] for p in lls], [p[1] for p in lls]
    lon0, lon1, lat0, lat1 = min(lons), max(lons), min(lats), max(lats)
    lonm, latm = (lon0 + lon1) / 2, (lat0 + lat1) / 2

    def _step(span):
        return next((s for s in (0.02, 0.05, 0.1, 0.2, 0.25, 0.5, 1, 2, 5)
                     if span / 5 <= s), 10)
    ls, ts = _step(lon1 - lon0), _step(lat1 - lat0)
    lon_t = np.arange(math.ceil(lon0 / ls) * ls, lon1 + 1e-9, ls)
    lat_t = np.arange(math.ceil(lat0 / ts) * ts, lat1 + 1e-9, ts)
    ax.set_xticks([fwd.transform(lo, latm)[0] for lo in lon_t])
    ax.set_yticks([fwd.transform(lonm, la)[1] for la in lat_t])
    ax.set_xticklabels([f"{abs(lo):.2f}°{'W' if lo < 0 else 'E'}"
                        for lo in lon_t])
    ax.set_yticklabels([f"{abs(la):.2f}°{'N' if la >= 0 else 'S'}"
                        for la in lat_t])
    ax.tick_params(direction="out", length=4, labelsize=7)


def drape_i15(hillshade, i15, out_path=None, work_crs="EPSG:5070",
              wet_min=5.0, vmax=None, cmap="YlGnBu", alpha=None, norm=None,
              cbar_ticks=None, cbar_ticklabels=None,
              perimeters=None, basins=None, highlight=None, points=None,
              gauges=None, gauge_value=None, gauge_cmap=None, gauge_vmax=None,
              gauge_label=None,
              title=None, cbar_label=None, figsize=(9, 8.5), dpi=200,
              basemap=False, basemap_provider="USGS.USTopo",
              basemap_labels=None, basemap_zoom=None, hillshade_alpha=None,
              hillshade_vmin=-20, hillshade_vmax=255,
              streams=None, roads=None, places=None, reference=False,
              local_roads=False, label_reference=True,
              clip=None, clip_margin=0.04,
              field_smooth=None, field_smooth_radius_km=0.0, smooth_power=2.0,
              north_arrow=False, scale_ticks=False, legend="all"):
    """Plot an i15 field draped over a hillshade with optional overlays.

    Parameters
    ----------
    hillshade, i15
        DataArrays or raster paths. ``i15`` is reprojected to match the
        hillshade grid; both are shown in ``work_crs``.
    wet_min
        i15 below this (mm/h) is left transparent so the basemap/hillshade
        shows through.
    alpha
        i15 drape opacity. Default ``None`` -> ``DEFAULT_FIELD_ALPHA`` (0.32),
        the project-wide default shared by every map so terrain/basemap reads
        through the rain; pass a value to override.
    vmax
        Colour-scale max; default is the 99th percentile of i15 (>= 10).
    norm, cbar_ticks, cbar_ticklabels
        Override the linear ``0..vmax`` scale with a matplotlib norm, and label
        the colour bar with something other than numbers. Together these draw a
        *classed* map -- a ``BoundaryNorm`` on the class breaks plus the class
        names as tick labels, which is how burn severity is published (see
        :func:`stormscape.burn.severity_colors`).
    perimeters, basins, highlight, points
        Optional vector overlays (path or GeoDataFrame, any CRS):
          * ``perimeters`` -- outlined white-over-black (e.g. fire/AOI border);
          * ``basins`` -- thin grey outlines (e.g. all candidate basins);
          * ``highlight`` -- bold cyan-over-black outlines (features of
            interest, e.g. basins that responded);
          * ``points`` -- yellow triangles (e.g. rain gauges, sites).
    gauges, gauge_value, gauge_cmap, gauge_vmax, gauge_label
        Rain-gauge overlay (GeoDataFrame or path). With ``gauge_value=None``
        the gauges are plain yellow markers; name a column to colour them. A
        value column (e.g. ``"i15_mmph"``) shares the i15 colour scale so the
        gauges and the radar field read alike, while a ``"resid_*"`` column
        (from :func:`stormscape.compare.radar_vs_gauge`) uses a diverging
        radar-minus-gauge scale with its own colour bar. ``gauge_cmap`` /
        ``gauge_vmax`` / ``gauge_label`` override the defaults.
    north_arrow, scale_ticks, legend
        Cartographic touches: ``north_arrow=True`` draws an up-arrow + "N";
        ``scale_ticks=True`` replaces the blank axes with distance tick marks
        (km from the lower-left corner); ``legend`` is ``"all"`` (default, every
        overlay), ``"gauges"`` (only the gauge marker), or ``False``/``None``.
    reference
        If ``True``, auto-fetch and overlay labelled vector reference layers
        for the map extent -- USGS **NHD** streams, **TIGER** roads, and
        **GNIS** place names (see :mod:`stormscape.refdata`) -- drawn directly
        on the hillshade. ``local_roads`` adds residential streets;
        ``label_reference=False`` draws the lines/points without text.
    streams, roads, places
        Explicit reference overlays (GeoDataFrame or path) drawn the same way;
        given these, ``reference`` need not auto-fetch. ``roads`` is styled by
        a ``kind`` column (primary/secondary/local) if present.
    clip, clip_margin
        Tighten the view to a geometry's extent. ``clip`` is a path /
        GeoDataFrame (e.g. the AOI or perimeter); the axes are limited to its
        bounds expanded by ``clip_margin`` (fraction, default 0.04). With
        ``reference``, the vectors are fetched for this tighter extent too.
    field_smooth, field_smooth_radius_km
        Optionally smooth the i15 field for *display* before draping (NaN-aware,
        on its native grid). ``field_smooth`` is a
        :data:`stormscape.smoothing.METHODS` key (e.g. ``"gaussian"``) and
        ``field_smooth_radius_km`` the nominal scale (~Gaussian sigma); the
        default (``None`` / ``0``) draws the raw field. Only the rendered figure
        is affected -- the source raster is untouched.
    hillshade_vmin, hillshade_vmax
        Grey-scale stretch of the hillshade. Defaults (-20, 255) are slightly
        darker than a full lift; lower ``hillshade_vmin`` to darken further.
    dpi
        Output resolution (default 200).
    basemap
        If truthy, underlay an open-source basemap (roads, rivers, place
        names) downloaded as map tiles via contextily. Pass ``True`` to use
        ``basemap_provider``, or a provider key string directly.
    basemap_provider
        contextily provider key for the base tiles. Default ``"USGS.USTopo"``
        -- the USGS National Map topographic basemap (public-domain USGS
        tiles with named creeks/rivers, roads, contours, and place names; the
        same National Map service pfdf draws from). Other good options:
        ``"USGS.USImageryTopo"`` (imagery + labels), ``"OpenStreetMap.Mapnik"``
        (dense roads/waterways), or a label-free base like
        ``"CartoDB.VoyagerNoLabels"`` paired with ``basemap_labels``.
    basemap_labels
        Optional provider key for a transparent labels-only layer drawn *on
        top* of the i15 field (e.g. ``"CartoDB.PositronOnlyLabels"``), useful
        when the base provider has no baked-in labels. Default ``None`` (the
        USGS topo base is already labelled, so the rain is kept lighter to let
        those labels read through instead).
    basemap_zoom
        Tile zoom level (int). ``None`` lets contextily choose for the AOI.
    hillshade_alpha
        Hillshade opacity. Default ``None`` -> opaque (1.0) with no basemap,
        and 0.0 with a basemap (the USGS topo base already conveys relief;
        raise it to blend the hillshade back in over a flat basemap).
    title, cbar_label
        Figure title and colour-bar label (sensible defaults).

    Returns
    -------
    (fig, ax)
        The figure is also written to ``out_path`` if given.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if isinstance(cmap, str) and cmap.startswith("cmc."):
        try:
            import cmcrameri.cm  # noqa: F401  registers the cmc.* colormaps
        except ImportError:
            pass

    def _load_i15():                       # load + (optionally) smooth the field
        da = _load(i15)
        if field_smooth and field_smooth_radius_km and field_smooth_radius_km > 0:
            from .smoothing import smooth_dataarray
            da = smooth_dataarray(da, field_smooth, field_smooth_radius_km,
                                  power=smooth_power)
        return da

    if hillshade is not None:
        hs0 = _load(hillshade)
        if str(work_crs).strip().upper() in ("UTM", "AUTO"):
            work_crs = _utm_crs_for(hs0)     # auto UTM zone -> near north-up
        hs = _to_crs(hs0, work_crs)
        i15da = _to_crs(_load_i15(), work_crs).rio.reproject_match(hs)
    else:                          # no terrain: show the i15 field on its own grid
        i150 = _load_i15()
        if str(work_crs).strip().upper() in ("UTM", "AUTO"):
            work_crs = _utm_crs_for(i150)
        hs = None
        i15da = _to_crs(i150, work_crs)
    i15v = np.nan_to_num(i15da.values)
    if i15v.ndim == 3:
        i15v = i15v[0]
    i15m = np.ma.masked_less(i15v, wet_min)
    if vmax is None and norm is None:
        finite = i15da.values[np.isfinite(i15da.values)]
        vmax = max(float(np.percentile(finite, 99)) if finite.size else 10.0,
                   10.0)
    ext = _extent(hs if hs is not None else i15da)

    # optional tighter view: clip the axes to a geometry's extent (+ margin)
    clip_gdf = read_overlay(clip, work_crs) if clip is not None else None
    clip_bounds = None
    if clip_gdf is not None and len(clip_gdf):
        x0, y0, x1, y1 = clip_gdf.total_bounds
        mx, my = (x1 - x0) * clip_margin, (y1 - y0) * clip_margin
        clip_bounds = (x0 - mx, y0 - my, x1 + mx, y1 + my)

    def _apply_clip():
        if clip_bounds:
            ax.set_xlim(clip_bounds[0], clip_bounds[2])
            ax.set_ylim(clip_bounds[1], clip_bounds[3])

    fig, ax = plt.subplots(figsize=figsize)
    i15_alpha = alpha if alpha is not None else DEFAULT_FIELD_ALPHA
    hs_alpha = (hillshade_alpha if hillshade_alpha is not None
                else (0.0 if basemap else 1.0))
    if hs is not None and hs_alpha > 0:              # skip if no / clear hillshade
        hsv = _fill_hillshade_nan(hs.values[0] if hs.values.ndim == 3
                                  else hs.values)
        ax.imshow(hsv, cmap="gray", extent=ext, origin="upper",
                  vmin=hillshade_vmin, vmax=hillshade_vmax, zorder=0.6,
                  alpha=hs_alpha)
    # a norm (e.g. the BoundaryNorm behind a classed burn-severity map) replaces
    # the linear vmin/vmax scale rather than fighting with it
    scale_kw = dict(norm=norm) if norm is not None else dict(vmin=0, vmax=vmax)
    im = ax.imshow(i15m, cmap=cmap, extent=ext, origin="upper",
                   alpha=i15_alpha, zorder=1,
                   interpolation="nearest", **scale_kw)
    _apply_clip()                      # set view before any basemap tile fetch
    if basemap:
        _add_basemap(ax, work_crs,
                     basemap if isinstance(basemap, str) else basemap_provider,
                     basemap_zoom, labels=basemap_labels)

    handles = []
    gauge_handles = []
    resid_scatter, resid_label = None, None
    bas = read_overlay(basins, work_crs)
    if bas is not None and len(bas):
        bas.boundary.plot(ax=ax, color="0.25", lw=0.7, alpha=0.7, zorder=2)
        handles.append(Line2D([], [], color="0.25", lw=1.0,
                              label="basin"))
    per = read_overlay(perimeters, work_crs)
    if per is not None and len(per):
        per.boundary.plot(ax=ax, color="white", lw=2.4, zorder=3)
        per.boundary.plot(ax=ax, color="black", lw=1.1, zorder=3)
        handles.append(Line2D([], [], color="black", lw=1.5,
                              label="perimeter"))
    hl = read_overlay(highlight, work_crs)
    if hl is not None and len(hl):
        hl.boundary.plot(ax=ax, color="black", lw=3.6, zorder=5)
        hl.boundary.plot(ax=ax, color="#00e5ff", lw=2.0, zorder=5)
        handles.append(Line2D([], [], color="#00e5ff", lw=2.2,
                              label="highlighted"))
    pts = read_overlay(points, work_crs)
    if pts is not None and len(pts):
        ax.scatter(pts.geometry.x, pts.geometry.y, marker="^", s=130,
                   facecolor="yellow", edgecolor="black", lw=1.2, zorder=6)
        handles.append(Line2D([], [], marker="^", color="yellow", mec="black",
                              lw=0, markersize=11, label="point"))

    # rain-gauge overlay: plain markers, value-coloured (shares the i15 scale),
    # or radar-minus-gauge residuals on a diverging scale with its own colour bar
    gg = read_overlay(gauges, work_crs) if gauges is not None else None
    if gg is not None and len(gg):
        if gauge_value is None:
            gh, _ = add_gauges(ax, work_crs, gg, label=gauge_label or "gauge")
        elif str(gauge_value).startswith("resid"):
            col = pd.to_numeric(gg.get(gauge_value), errors="coerce")
            m = (gauge_vmax if gauge_vmax is not None
                 else (float(np.nanmax(np.abs(col))) if col.notna().any()
                       else 1.0)) or 1.0
            gh, gsc = add_gauges(ax, work_crs, gg, value=gauge_value,
                                 cmap=gauge_cmap or "RdBu_r", vmin=-m, vmax=m,
                                 label=gauge_label or "gauge (radar - gauge)")
            resid_scatter = gsc          # colour bar built tight to the map below
            resid_label = gauge_label or "radar - gauge  (mm h$^{-1}$)"
        else:                              # value mode: share the field's scale
            gh, _ = add_gauges(ax, work_crs, gg, value=gauge_value,
                               cmap=gauge_cmap or cmap, vmin=0,
                               vmax=gauge_vmax or vmax,
                               label=gauge_label or "gauge")
        handles += gh
        gauge_handles = gh

    # vector reference layers (NHD streams / TIGER roads / GNIS places)
    if reference and streams is None and roads is None and places is None:
        from . import refdata
        # fetch over the (tighter) clip extent if clipping, else the full map
        if clip_gdf is not None and len(clip_gdf):
            bnds = tuple(clip_gdf.to_crs(4326).total_bounds)
        else:
            bnds = _bounds_4326(hs)
        # named watercourses only -> the streams a topo map labels (NHD HR
        # carries every headwater rill, which would bury the i15 field)
        streams = refdata.streams(bnds, named_only=False)   # incl. minor/unnamed
        roads = refdata.roads(bnds, local=local_roads)
        places = refdata.places(bnds)
    if streams is not None or roads is not None or places is not None:
        handles += add_reference(ax, work_crs, streams=streams, roads=roads,
                                 places=places, label=label_reference)

    # colour bars pinned tight to the map (append_axes hugs the axes regardless
    # of the figure aspect, so a tall narrow UTM map doesn't strand them)
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    div = make_axes_locatable(ax)
    cax = div.append_axes("right", size="4%", pad=0.12)
    cb = fig.colorbar(im, cax=cax,
                      **(dict(ticks=cbar_ticks) if cbar_ticks is not None
                         else {}))
    if cbar_ticklabels is not None:      # e.g. severity class names, not numbers
        cb.ax.set_yticklabels(cbar_ticklabels, fontsize=9)
    cb.set_label(cbar_label or
                 "peak 15-min rainfall intensity  i$_{15}$  (mm h$^{-1}$)",
                 fontsize=10)
    if resid_scatter is not None:        # diverging radar-minus-gauge bar, left
        caxl = div.append_axes("left", size="4%", pad=0.75)
        cbl = fig.colorbar(resid_scatter, cax=caxl)
        caxl.yaxis.set_ticks_position("left")
        caxl.yaxis.set_label_position("left")
        cbl.set_label(resid_label, fontsize=9)
    if title:
        ax.set_title(title, fontsize=12)
    # legend: "all" (every overlay), "gauges" (only the gauge marker), or off
    if legend and str(legend).lower() not in ("none", "false"):
        lh = (gauge_handles if str(legend).lower() in ("gauge", "gauges")
              else handles)
        if lh:
            ax.legend(handles=lh, loc="upper right", fontsize=9, framealpha=0.9)
    # frame on the clip extent, else the radar field, so stray off-field
    # gauges / overlays cannot stretch the axes
    if clip_bounds:
        ax.set_xlim(clip_bounds[0], clip_bounds[2])
        ax.set_ylim(clip_bounds[1], clip_bounds[3])
    else:
        ax.set_xlim(ext[0], ext[1])
        ax.set_ylim(ext[2], ext[3])
    if scale_ticks:                    # lat/long coordinate tick marks
        _add_scale_ticks(ax, work_crs)
    else:
        ax.set_xticks([])
        ax.set_yticks([])
    if north_arrow:
        _add_north_arrow(ax)
    if out_path:
        fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    return fig, ax


# field name -> (title, cmap, vmin, vmax, unit); vmin/vmax None -> data range.
# Colormaps are colourblind-safe (per the scientific-plots convention); never
# jet/turbo as the original MRMS_stack used.
_PANEL_SPECS = {
    "tpki15": ("Time of peak i15", "viridis", None, None, "UTC hour"),
    "total": ("QPE storm total", "YlGnBu", 0.0, None, "mm"),
    "total_mm": ("Storm total", "YlGnBu", 0.0, None, "mm"),
    "rqi": ("Radar quality index", "viridis", 0.0, 1.0, "RQI [-]"),
    "shsr": ("Beam height (SHSR)", "viridis_r", 0.0, None, "km AGL"),
    "cbb": ("Beam blockage (CBB)", "inferno", 0.0, 1.0, "fraction"),
    "i15max": ("Peak i15", "YlGnBu", 0.0, None, "mm h$^{-1}$"),
    "i30max": ("Peak i30", "YlGnBu", 0.0, None, "mm h$^{-1}$"),
    "i60max": ("Peak i60", "YlGnBu", 0.0, None, "mm h$^{-1}$"),
    "i2max": ("Peak i2", "YlGnBu", 0.0, None, "mm h$^{-1}$"),
    "peakrate_mmph": ("Peak rate", "YlGnBu", 0.0, None, "mm h$^{-1}$"),
    # burn severity (stormscape.burn). Not in _MASK_DRY: BRISK dNBR is already
    # NaN outside the burn, and that mask's 0.5 cut would erase everything below
    # high severity -- dNBR runs about -0.3 to 1.
    "dnbr": ("Burn severity (dNBR)", "YlOrRd", 0.0, 1.0, "dNBR [-]"),
    "severity": ("Burn severity class", "YlOrRd", 0.0, None, "class"),
}


# fields where 0 = dry: mask them so the hillshade shows through (like the i15
# drape). Others (tpki15/rqi/shsr/cbb) are shown wherever they have data.
_MASK_DRY = {"total", "total_mm", "i15max", "i30max", "i60max", "i2max",
             "peakrate_mmph", "mstotal"}


def _fill_hillshade_nan(hsv):
    """Replace the rotated 5070->UTM NaN corners of a reprojected hillshade with the
    mean terrain tone, so the terrain visually fills the view frame instead of the
    corners rendering as stark white wedges or matplotlib's fragile, layout-
    dependent resize gray-fill (the "rotated gray box" relic behind the map)."""
    hsv = np.asarray(hsv, dtype=float)
    finite = np.isfinite(hsv)
    if finite.any() and not finite.all():
        hsv = np.where(finite, hsv, float(np.nanmean(hsv)))
    return hsv


def _draw_field(ax, arr, ext, wc, *, hs=None, cmap="viridis", vmin=None,
                vmax=None, norm=None, mask_dry=False, alpha=None,
                hillshade_alpha=1.0, streams=None, roads=None, places=None,
                reference=False, label_reference=True, perim=None, gauges=None,
                clip_bounds=None, scale_ticks=True, north_arrow=True):
    """Draw one field array over the (optional) hillshade with the toolkit's
    shared map context -- reference vectors, perimeter, gauges, clip limits,
    lat/long ticks, north arrow. ``arr`` must already be on the hillshade grid
    (the caller reprojects). Returns the field's imshow handle so the caller can
    attach its own colorbar/title. The per-panel core of :func:`drape_i15` /
    :func:`diagnostic_panels`, shared by the climatology figures."""
    if hs is not None and hillshade_alpha > 0:
        hsv = hs.values[0] if hs.values.ndim == 3 else hs.values
        # Fill the rotated 5070->UTM NaN corners with the mean terrain tone so the
        # terrain fills the frame (like the production i15 maps / drape_i15) instead
        # of rendering as white wedges or a fragile, layout-dependent resize
        # gray-fill -- the "rotated gray box" relic. An explicit fill is reliable.
        hsv = _fill_hillshade_nan(hsv)
        ax.imshow(hsv, cmap="gray", extent=ext, origin="upper",
                  vmin=-20, vmax=255, zorder=0.6, alpha=hillshade_alpha)
    show = (np.ma.masked_less(np.nan_to_num(arr, nan=-1.0), 0.5)
            if mask_dry else arr)
    kw = dict(norm=norm) if norm is not None else dict(vmin=vmin, vmax=vmax)
    a = alpha if alpha is not None else DEFAULT_FIELD_ALPHA
    im = ax.imshow(show, extent=ext, origin="upper", cmap=cmap, zorder=1,
                   interpolation="nearest",
                   alpha=a if hs is not None else 1.0, **kw)
    if reference:
        add_reference(ax, wc, streams=streams, roads=roads, places=places,
                      label=label_reference, max_stream_labels=8,
                      max_road_labels=8, max_place_labels=10)
    if perim is not None and len(perim):
        perim.boundary.plot(ax=ax, color="white", lw=2.0, zorder=3)
        perim.boundary.plot(ax=ax, color="black", lw=0.9, zorder=3)
    if gauges is not None and len(gauges):
        ax.scatter(gauges.geometry.x, gauges.geometry.y, s=16, c="black",
                   edgecolor="white", linewidth=0.5, zorder=6)
    if clip_bounds:
        ax.set_xlim(clip_bounds[0], clip_bounds[1])
        ax.set_ylim(clip_bounds[2], clip_bounds[3])
    if scale_ticks:
        _add_scale_ticks(ax, wc)
    else:
        ax.set_xticks([])
        ax.set_yticks([])
    if north_arrow:
        _add_north_arrow(ax)
    return im


def diagnostic_panels(radar_dir, key, which=("tpki15", "total", "rqi", "shsr"),
                      out_path=None, work_crs="UTM", hillshade=None,
                      perimeters=None, gauges=None, reference=False,
                      local_roads=False, label_reference=True, alpha=None,
                      hillshade_alpha=1.0, north_arrow=True, scale_ticks=True,
                      clip=None, clip_margin=0.04, title=None, ncols=2,
                      field_smooth=None, field_smooth_radius_km=0.0,
                      smooth_power=2.0, figsize=None, dpi=200):
    """Multi-panel diagnostic map of stacked radar fields.

    Reads ``<radar_dir>/<key>_<field>.tif`` for each name in ``which`` (default:
    time-of-peak i15, QPE total, RQI, beam-height SHSR) and tiles them, sharing
    the **same context as the toolkit's main maps** -- optional hillshade
    underlay, labelled reference vectors (NHD streams / TIGER roads / GNIS
    places), AOI perimeter, gauge markers, north arrow and lat/long scale ticks.
    Colourblind-safe colormaps; time-of-peak is unwrapped across midnight. A port
    (and upgrade) of the multi-map overview in D. Cavagna's ``MRMS_stack``.

    If ``field_smooth`` (a :data:`stormscape.smoothing.METHODS` key) is set, the
    **rainfall** panels (intensities + depths, i.e. those in :data:`_MASK_DRY`)
    are smoothed for display at ``field_smooth_radius_km`` -- the categorical /
    quality panels (tpki15 / rqi / shsr / cbb) are always left raw. Display only;
    the source rasters are untouched.
    """
    import math

    import matplotlib.pyplot as plt
    # smooth only the rainfall fields (those in _MASK_DRY: intensities + depths);
    # leave categorical/quality fields (tpki15/rqi/shsr/cbb) raw
    smoothing_on = (bool(field_smooth) and field_smooth_radius_km
                    and field_smooth_radius_km > 0)
    smooth_label = ""
    if smoothing_on:
        from .smoothing import METHODS as _SM
        from .smoothing import smooth_dataarray
        smooth_label = f"{_SM.get(field_smooth, field_smooth)} {field_smooth_radius_km:g} km"
    das = {}
    for name in which:
        p = find(radar_dir, f"{key}_{name}.tif")
        if os.path.exists(p):
            da = _load(p)
            if smoothing_on and name in _MASK_DRY:
                da = smooth_dataarray(da, field_smooth, field_smooth_radius_km,
                                      power=smooth_power)
            das[name] = da
    if not das:
        raise FileNotFoundError(
            f"no {tuple(which)} field tifs for key '{key}' in {radar_dir}")

    base0 = _load(hillshade) if hillshade is not None else next(iter(das.values()))
    wc = (_utm_crs_for(base0)
          if str(work_crs).strip().upper() in ("UTM", "AUTO") else work_crs)
    hs = _to_crs(base0, wc) if hillshade is not None else None
    perim = read_overlay(perimeters, wc) if perimeters is not None else None
    g = read_overlay(gauges, wc) if gauges is not None else None

    clip_gdf = read_overlay(clip, wc) if clip is not None else None
    clip_bounds = None
    if clip_gdf is not None and len(clip_gdf):
        x0, y0, x1, y1 = clip_gdf.total_bounds
        mx, my = (x1 - x0) * clip_margin, (y1 - y0) * clip_margin
        clip_bounds = (x0 - mx, x1 + mx, y0 - my, y1 + my)

    streams = roads = places = None
    if reference:
        from . import refdata
        bnds = (tuple(clip_gdf.to_crs(4326).total_bounds)
                if clip_gdf is not None and len(clip_gdf)
                else _bounds_4326(hs if hs is not None else base0))
        streams = refdata.streams(bnds, named_only=False)
        roads = refdata.roads(bnds, local=local_roads)
        places = refdata.places(bnds)

    n = len(das)
    ncols = min(ncols, n)
    nrows = math.ceil(n / ncols)
    # size the figure from the map aspect so the equal-aspect panels don't float
    # inside over-wide cells -- that float is the dead white space between
    # columns for tall, narrow AOIs (e.g. Hidden Valley ~0.4 wide:tall).
    base_wc = hs if hs is not None else _to_crs(base0, wc)
    if clip_bounds:
        map_w, map_h = (clip_bounds[1] - clip_bounds[0],
                        clip_bounds[3] - clip_bounds[2])
    else:
        le, ri, bo, to = _extent(base_wc)
        map_w, map_h = ri - le, to - bo
    aspect = (map_w / map_h) if map_h else 1.0
    panel_h = 5.6                                    # map height per row (inches)
    figsize = figsize or ((panel_h * aspect + 1.5) * ncols,
                          panel_h * nrows + (0.5 if title else 0.0))
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
    axes = axes.flatten()
    for ax, (name, da0) in zip(axes, das.items()):
        da = _to_crs(da0, wc)
        if hs is not None:
            da = da.rio.reproject_match(hs)
        arr = np.asarray(da.squeeze().values, dtype=float)
        ttl, cmap, vmin, vmax, unit = _PANEL_SPECS.get(
            name, (name, "viridis", None, None, ""))
        if smoothing_on and name in _MASK_DRY:                 # note the smoothing
            ttl = f"{ttl} ({smooth_label})"
        if name == "tpki15" and np.isfinite(arr).any():       # unwrap midnight
            arr = np.where(arr < np.nanmedian(arr) - 12, arr + 24, arr)
        finite = arr[np.isfinite(arr)]
        vmx = vmax if vmax is not None else (
            float(np.nanpercentile(finite, 99)) if finite.size else 1.0)
        vmn = vmin if vmin is not None else (
            float(np.nanmin(finite)) if finite.size else 0.0)
        ext = _extent(hs if hs is not None else da)
        im = _draw_field(ax, arr, ext, wc, hs=hs, cmap=cmap, vmin=vmn, vmax=vmx,
                         mask_dry=(name in _MASK_DRY), alpha=alpha,
                         hillshade_alpha=hillshade_alpha, streams=streams,
                         roads=roads, places=places, reference=reference,
                         label_reference=label_reference, perim=perim, gauges=g,
                         clip_bounds=clip_bounds, scale_ticks=scale_ticks,
                         north_arrow=north_arrow)
        cb = fig.colorbar(im, ax=ax, shrink=0.82, pad=0.02)
        cb.ax.tick_params(labelsize=8)
        cb.set_label(unit, fontsize=8)
        ax.set_title(ttl, fontsize=10)
    for ax in axes[n:]:
        ax.axis("off")
    if title:
        fig.suptitle(title, fontsize=12)
    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
    return out_path or (fig, axes)


def _setup_context(base0, hs_da, work_crs, perimeters, gauges, clip,
                   clip_margin, reference, local_roads):
    """Resolve the shared map context for the climatology figures -- work CRS,
    reprojected hillshade, vector overlays, clip box, and (optionally) fetched
    reference vectors. Mirrors the setup block in :func:`diagnostic_panels`."""
    wc = (_utm_crs_for(base0)
          if str(work_crs).strip().upper() in ("UTM", "AUTO") else work_crs)
    hs = _to_crs(hs_da, wc) if hs_da is not None else None
    perim = read_overlay(perimeters, wc) if perimeters is not None else None
    g = read_overlay(gauges, wc) if gauges is not None else None
    clip_gdf = read_overlay(clip, wc) if clip is not None else None
    clip_bounds = None
    if clip_gdf is not None and len(clip_gdf):
        x0, y0, x1, y1 = clip_gdf.total_bounds
        mx, my = (x1 - x0) * clip_margin, (y1 - y0) * clip_margin
        clip_bounds = (x0 - mx, x1 + mx, y0 - my, y1 + my)
    # NOTE: the hillshade is intentionally NOT clipped here -- the unclipped
    # terrain fills the view frame (like drape_i15), so the maps don't float in
    # white space. The rotated-gray-box relic is prevented in _draw_field (the
    # hillshade is drawn as explicit RGBA with alpha=0 where masked, so the
    # masked corners stay transparent instead of being gray-filled on resize).
    streams = roads = places = None
    if reference:
        from . import refdata
        bnds = (tuple(clip_gdf.to_crs(4326).total_bounds)
                if clip_gdf is not None and len(clip_gdf)
                else _bounds_4326(hs if hs is not None else base0))
        streams = refdata.streams(bnds, named_only=False)
        roads = refdata.roads(bnds, local=local_roads)
        places = refdata.places(bnds)
    return wc, hs, perim, g, clip_bounds, streams, roads, places


def climatology_comparison(clim_dir, obs_dir, key, durations=(15, 30, 60),
                           obs_key=None, out_path=None, work_crs="UTM",
                           hillshade=None, perimeters=None, gauges=None,
                           reference=False, local_roads=False,
                           label_reference=True, alpha=None, hillshade_alpha=1.0,
                           north_arrow=True, scale_ticks=True, clip=None,
                           clip_margin=0.04, cmap="YlGnBu",
                           shared_row_scale=False, ari=1, title=None,
                           obs_smooth=None, obs_smooth_radius_km=0.0,
                           smooth_power=2.0, figsize=None, dpi=200):
    """3xN figure: observed peak intensity vs NOAA Atlas 14 climatology.

    Rows are durations (15/30/60 min); **column 1** is the Atlas 14 ``ari``-year
    climatology (``<clim_dir>/<key>_clim_i{d}.tif``), **column 2** the observed
    field (``<obs_dir>/<key>_i{d}max.tif``, MRMS or NEXRAD). Each panel carries
    the toolkit's full map context (hillshade, labelled reference vectors,
    perimeter, gauges, ticks). By default every panel has its **own** sequential
    colour scale -- the climatology is smooth and small next to the peaky
    observed field, so a shared scale would wash it out; the magnitude
    comparison is the job of :func:`anomaly_map`. ``shared_row_scale=True``
    shares vmin/vmax within a row instead.

    If ``obs_smooth`` (a :data:`stormscape.smoothing.METHODS` key) is set, the
    observed field is spatially smoothed (NaN-aware) at ``obs_smooth_radius_km``
    on its native grid before display -- so the peaky ~1 km radar field reads
    against the smooth ~800 m climatology rather than dominating with speckle.
    """
    import matplotlib.pyplot as plt
    okey = obs_key or key
    smoothing_on = bool(obs_smooth) and obs_smooth_radius_km and obs_smooth_radius_km > 0
    rows = []
    for d in durations:
        cp = find(clim_dir, f"{key}_clim_i{d}.tif")
        op = find(obs_dir, f"{okey}_i{d}max.tif")
        if os.path.exists(cp) and os.path.exists(op):
            obs_da = _load(op)
            if smoothing_on:
                from .smoothing import smooth_dataarray
                obs_da = smooth_dataarray(obs_da, obs_smooth,
                                          obs_smooth_radius_km, power=smooth_power)
            rows.append((d, _load(cp), obs_da))
    if not rows:
        raise FileNotFoundError(
            f"no clim/observed i<d>.tif pairs for key '{key}' "
            f"(clim_dir={clim_dir}, obs_dir={obs_dir})")

    hs_da = _load(hillshade) if hillshade is not None else None
    base0 = hs_da if hs_da is not None else rows[0][1]
    wc, hs, perim, g, clip_bounds, streams, roads, places = _setup_context(
        base0, hs_da, work_crs, perimeters, gauges, clip, clip_margin,
        reference, local_roads)

    base_wc = hs if hs is not None else _to_crs(base0, wc)
    if clip_bounds:
        map_w, map_h = (clip_bounds[1] - clip_bounds[0],
                        clip_bounds[3] - clip_bounds[2])
    else:
        le, ri, bo, to = _extent(base_wc)
        map_w, map_h = ri - le, to - bo
    aspect = (map_w / map_h) if map_h else 1.0
    panel_h = 4.6
    nrows = len(rows)
    figsize = figsize or ((panel_h * aspect + 1.6) * 2, panel_h * nrows + 0.6)
    fig, axes = plt.subplots(nrows, 2, figsize=figsize, squeeze=False)

    obs_title = "Observed"
    if smoothing_on:
        from .smoothing import METHODS as _SM
        obs_title = (f"Observed ({_SM.get(obs_smooth, obs_smooth)}, "
                     f"{obs_smooth_radius_km:g} km)")
    col_titles = (f"NOAA Atlas 14 — {ari}-yr", obs_title)
    for r, (d, clim_da, obs_da) in enumerate(rows):
        arrs = []
        for da0 in (clim_da, obs_da):
            da = _to_crs(da0, wc)
            if hs is not None:
                da = da.rio.reproject_match(hs)
            arrs.append(np.asarray(da.squeeze().values, dtype=float))
        ext = _extent(hs if hs is not None else _to_crs(clim_da, wc))

        def _lohi(vals):
            # 2nd-99th percentile -> each panel uses its own dynamic range, so
            # the smooth/narrow climatology field shows its gradient instead of
            # sitting flat in the middle of a 0-anchored scale.
            v = vals[np.isfinite(vals)]
            if not v.size:
                return (0.0, 1.0)
            return (float(np.nanpercentile(v, 2)), float(np.nanpercentile(v, 99)))
        if shared_row_scale:
            both = np.concatenate([a[np.isfinite(a)] for a in arrs]
                                  + [np.array([0.0])])
            scales = [_lohi(both)] * 2
        else:
            scales = [_lohi(a) for a in arrs]
        for c in range(2):
            ax = axes[r][c]
            vmn, vmx = scales[c]
            im = _draw_field(
                ax, arrs[c], ext, wc, hs=hs, cmap=cmap, vmin=vmn, vmax=vmx,
                mask_dry=True, alpha=alpha, hillshade_alpha=hillshade_alpha,
                streams=streams, roads=roads, places=places,
                reference=reference, label_reference=label_reference,
                perim=perim, gauges=g, clip_bounds=clip_bounds,
                scale_ticks=scale_ticks,
                north_arrow=(north_arrow and r == 0 and c == 0))
            cb = fig.colorbar(im, ax=ax, shrink=0.82, pad=0.02)
            cb.ax.tick_params(labelsize=8)
            cb.set_label("mm h$^{-1}$", fontsize=8)
            if r == 0:
                ax.set_title(col_titles[c], fontsize=11, fontweight="bold")
        axes[r][0].set_ylabel(f"I$_{{{d}}}$", fontsize=13, fontweight="bold")
    if title:
        fig.suptitle(title, fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.98] if title else None)
    if out_path:
        fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
    return out_path or (fig, axes)


def anomaly_map(anomaly, hillshade=None, out_path=None, work_crs="UTM",
                duration=15, ari=1, cmap="cmc.vik", vmax=None,
                contour_levels=None, perimeters=None, gauges=None,
                reference=False, local_roads=False, label_reference=True,
                alpha=None, hillshade_alpha=1.0, north_arrow=True,
                scale_ticks=True, clip=None, clip_margin=0.04, title=None,
                cbar_label=None, figsize=(9, 8.5), dpi=200):
    """Single-duration rainfall anomaly map: observed / climatology.

    ``anomaly`` is a ratio DataArray or tif (from
    :func:`stormscape.atlas14.anomaly`). Drawn over the hillshade with a
    **diverging** colormap centred at **1.0** (``TwoSlopeNorm``; per the
    scientific-plots default, ``cmc.vik`` if cmcrameri is present, else
    ``RdBu_r``) so cells below the ``ari``-year storm read cool and cells above
    read warm. **Integer contours** (1x, 2x, 3x, ...) are drawn to emphasize the
    recurrence-multiple breaks. Full map context + north arrow + lat/long ticks.
    """
    import math

    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    if isinstance(cmap, str) and cmap.startswith("cmc."):
        try:
            import cmcrameri.cm  # noqa: F401  registers the cmc.* colormaps
        except ImportError:
            cmap = "RdBu_r"

    da0 = _load(anomaly)
    hs_da = _load(hillshade) if hillshade is not None else None
    base0 = hs_da if hs_da is not None else da0
    wc, hs, perim, g, clip_bounds, streams, roads, places = _setup_context(
        base0, hs_da, work_crs, perimeters, gauges, clip, clip_margin,
        reference, local_roads)
    da = _to_crs(da0, wc)
    if hs is not None:
        da = da.rio.reproject_match(hs)
    arr = np.asarray(da.squeeze().values, dtype=float)
    ext = _extent(hs if hs is not None else da)
    finite = arr[np.isfinite(arr)]
    vmx = vmax if vmax is not None else (
        max(2.0, float(np.nanpercentile(finite, 99))) if finite.size else 2.0)
    norm = TwoSlopeNorm(vcenter=1.0, vmin=0.0, vmax=vmx)

    fig, ax = plt.subplots(figsize=figsize)
    im = _draw_field(ax, arr, ext, wc, hs=hs, cmap=cmap, norm=norm,
                     mask_dry=False, alpha=alpha,
                     hillshade_alpha=hillshade_alpha, streams=streams,
                     roads=roads, places=places, reference=reference,
                     label_reference=label_reference, perim=perim, gauges=g,
                     clip_bounds=clip_bounds, scale_ticks=scale_ticks,
                     north_arrow=north_arrow)
    levels = contour_levels or list(range(1, int(math.ceil(vmx)) + 1))
    if levels and finite.size:
        cs = ax.contour(arr, levels=levels, extent=ext, origin="upper",
                        colors="black", linewidths=0.8, zorder=4, alpha=0.85)
        ax.clabel(cs, fmt={lev: f"{lev:g}×" for lev in levels},
                  fontsize=8, inline=True)
    cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02,
                      ticks=list(range(0, int(math.ceil(vmx)) + 1)))
    cb.set_label(cbar_label or
                 f"observed I$_{{{duration}}}$ ÷ {ari}-yr climatology  "
                 f"(×)", fontsize=10)
    ax.set_title(title or
                 f"I$_{{{duration}}}$ anomaly  (observed / {ari}-yr climatology)",
                 fontsize=12)
    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
    return out_path or (fig, ax)


def smoothing_comparison(in_dir, key, field="i15max",
                         methods=("gaussian", "uniform", "median", "idw"),
                         radii_km=(0, 1, 2, 4), out_path=None, work_crs="UTM",
                         hillshade=None, perimeters=None, gauges=None,
                         reference=False, local_roads=False,
                         label_reference=True, alpha=None, hillshade_alpha=1.0,
                         north_arrow=True, scale_ticks=True, clip=None,
                         clip_margin=0.04, cmap="YlGnBu", shared_scale=True,
                         power=2.0, title=None, figsize=None, dpi=200):
    """Grid of smoothing methods × radii for one radar field.

    Reads ``<in_dir>/<key>_<field>.tif`` (default the peak-i15 field) and tiles a
    **methods × radii** grid: rows are smoothing methods
    (:data:`stormscape.smoothing.METHODS`), columns are ``radii_km`` with the
    first column (``radius 0``) the raw field. Each panel is smoothed on the
    native grid (before any reprojection warp) then draped over the hillshade with
    the toolkit's full map context. A **single shared colour scale** (``vmin=0``,
    ``vmax`` = 99th percentile of the *raw* field) spans all panels so the
    peak-flattening from smoothing is visible -- unlike
    :func:`climatology_comparison`'s per-panel default. ``shared_scale=False``
    scales (and colour-bars) each panel on its own. A single method gives a 1×N
    scale strip.
    """
    import matplotlib.pyplot as plt

    from . import smoothing as _sm
    if hillshade_alpha is None:                        # opaque base (no basemap here)
        hillshade_alpha = 1.0
    fpath = find(in_dir, f"{key}_{field}.tif")
    if not os.path.exists(fpath):
        raise FileNotFoundError(f"field tif not found: {fpath}")
    da_native = _load(fpath)
    native = np.asarray(da_native.squeeze().values, dtype=float)
    tr = da_native.rio.transform()
    b = da_native.rio.bounds()
    if da_native.rio.crs is not None and da_native.rio.crs.is_geographic:
        cell_km = _sm.cell_size_km(tr, 0.5 * (b[1] + b[3]))
    else:
        cell_km = 0.5 * (abs(tr.a) + abs(tr.e)) / 1000.0

    methods, radii = list(methods), list(radii_km)
    nrows, ncols = len(methods), len(radii)
    hs_da = _load(hillshade) if hillshade is not None else None
    base0 = hs_da if hs_da is not None else da_native
    wc, hs, perim, g, clip_bounds, streams, roads, places = _setup_context(
        base0, hs_da, work_crs, perimeters, gauges, clip, clip_margin,
        reference, local_roads)

    raw_fin = native[np.isfinite(native)]
    vmax_shared = float(np.nanpercentile(raw_fin, 99)) if raw_fin.size else 1.0
    spec = _PANEL_SPECS.get(field, (field, cmap, 0.0, None, "mm h$^{-1}$"))
    unit = spec[4]

    base_wc = hs if hs is not None else _to_crs(base0, wc)
    if clip_bounds:
        map_w, map_h = (clip_bounds[1] - clip_bounds[0],
                        clip_bounds[3] - clip_bounds[2])
    else:
        le, ri, bo, to = _extent(base_wc)
        map_w, map_h = ri - le, to - bo
    aspect = (map_w / map_h) if map_h else 1.0
    panel_h = 3.6
    figsize = figsize or (panel_h * aspect * ncols + 2.2, panel_h * nrows + 1.0)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False,
                             constrained_layout=True)

    mask_dry = field in _MASK_DRY
    last_im = None
    for ri_, method in enumerate(methods):
        for ci, r in enumerate(radii):
            ax = axes[ri_][ci]
            sm_arr = _sm.smooth_array(native, cell_km, method, r, power=power)
            sm_da = da_native.copy(data=sm_arr).rio.write_crs(da_native.rio.crs)
            da = _to_crs(sm_da, wc)
            if hs is not None:
                da = da.rio.reproject_match(hs)
            arr = np.asarray(da.squeeze().values, dtype=float)
            ext = _extent(hs if hs is not None else da)
            if shared_scale:
                vmn, vmx = 0.0, vmax_shared
            else:
                fin = arr[np.isfinite(arr)]
                vmn, vmx = 0.0, (float(np.nanpercentile(fin, 99))
                                 if fin.size else 1.0)
            last_im = _draw_field(
                ax, arr, ext, wc, hs=hs, cmap=cmap, vmin=vmn, vmax=vmx,
                mask_dry=mask_dry, alpha=alpha, hillshade_alpha=hillshade_alpha,
                streams=streams, roads=roads, places=places, reference=reference,
                label_reference=label_reference, perim=perim, gauges=g,
                clip_bounds=clip_bounds,
                scale_ticks=(scale_ticks and ri_ == nrows - 1 and ci == 0),
                north_arrow=(north_arrow and ri_ == 0 and ci == 0))
            if not shared_scale:
                cb = fig.colorbar(last_im, ax=ax, shrink=0.82, pad=0.02)
                cb.ax.tick_params(labelsize=7)
            if ri_ == 0:
                ax.set_title("raw" if not r else f"{r:g} km", fontsize=11,
                             fontweight="bold")
        axes[ri_][0].set_ylabel(_sm.METHODS.get(method, method), fontsize=11,
                                fontweight="bold")
    if shared_scale and last_im is not None:
        cb = fig.colorbar(last_im, ax=axes.ravel().tolist(), shrink=0.6, pad=0.02)
        cb.set_label(unit, fontsize=9)
    fig.suptitle(title or f"{key}  —  {spec[0]} smoothing comparison",
                 fontsize=13)
    if out_path:
        fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
    return out_path or (fig, axes)


def smoothing_skill_plot(sweep_df, durations=(15, 30, 60),
                         metrics=("corr", "rmse", "ratio"), out_path=None,
                         figsize=None, dpi=200):
    """Radar–gauge skill vs smoothing radius, per method and duration.

    Takes the tidy DataFrame from
    :func:`stormscape.smoothing.gauge_skill_sweep`. Rows are ``metrics`` (default
    correlation / RMSE / bias-ratio), columns are ``durations``. Each panel plots
    the metric vs smoothing radius with one line per method; the raw (radius 0)
    value is a hollow square and the optimum is starred (correlation max / RMSE
    min). Because i15 is a *peak* metric, smoothing drives the bias **ratio**
    through 1.0 as a mechanical side-effect -- the ratio panel marks 1.0 but is
    **not** the skill criterion; trust correlation (up) and RMSE (down).
    """
    import matplotlib.pyplot as plt

    from . import smoothing as _sm
    palette = {"gaussian": "#0072B2", "uniform": "#E69F00",
               "median": "#009E73", "idw": "#D55E00"}
    mlabel = {"corr": "correlation r", "rmse": "RMSE (mm h$^{-1}$)",
              "mae": "MAE (mm h$^{-1}$)", "bias": "bias (mm h$^{-1}$)",
              "ratio": "bias ratio  Σradar / Σgauge"}
    present = set(sweep_df["duration"].unique())
    durations = [d for d in durations if d in present]
    metrics = list(metrics)
    methods = list(dict.fromkeys(sweep_df["method"]))
    nrows, ncols = len(metrics), max(1, len(durations))
    figsize = figsize or (3.5 * ncols + 0.8, 2.7 * nrows + 1.0)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False,
                             sharex=True)
    for mi, metric in enumerate(metrics):
        for di, d in enumerate(durations):
            ax = axes[mi][di]
            sub_d = sweep_df[sweep_df["duration"] == d]
            for method in methods:
                s = sub_d[sub_d["method"] == method].sort_values("radius_km")
                if not len(s) or metric not in s.columns:
                    continue
                color = palette.get(method)
                ax.plot(s["radius_km"], s[metric], "-o", ms=3, lw=1.3,
                        color=color, label=_sm.METHODS.get(method, method))
                raw = s[s["radius_km"] == 0]
                if len(raw):
                    ax.plot(raw["radius_km"], raw[metric], "s", ms=6,
                            mfc="white", mec=color, zorder=5)
                if metric in ("corr", "rmse", "mae"):
                    pos = s[s["radius_km"] > 0].dropna(subset=[metric])
                    if len(pos):
                        idx = (pos[metric].idxmin() if metric in ("rmse", "mae")
                               else pos[metric].idxmax())
                        ax.plot(pos.loc[idx, "radius_km"], pos.loc[idx, metric],
                                "*", ms=14, color=color, mec="black", mew=0.4,
                                zorder=6)
            if metric == "ratio":
                ax.axhline(1.0, color="0.4", ls="--", lw=1.0)
            if mi == 0:
                ax.set_title(f"I$_{{{d}}}$", fontsize=12, fontweight="bold")
            if di == 0:
                ax.set_ylabel(mlabel.get(metric, metric), fontsize=10)
            if mi == nrows - 1:
                ax.set_xlabel("smoothing radius (km)", fontsize=10)
            ax.grid(True, alpha=0.3)
    handles, lbls = axes[0][0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, lbls, loc="lower center", ncol=len(handles),
                   fontsize=10, frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Radar–gauge skill vs smoothing radius "
                 "(★ optimum; ▫ raw)", fontsize=13)
    fig.tight_layout(rect=[0, 0.02, 1, 0.96])
    if out_path:
        fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
    return out_path or (fig, axes)


def plot_virtual_gauge(df, name="virtual gauge", out_path=None,
                       durations=(5, 15, 30, 60), figsize=(7, 5.5), dpi=200):
    """Two-panel virtual-gauge figure (radar rainfall time series at a point).

    Top: trailing-window intensities ``i{d}`` from PrecipRate plus the hourly
    gauge-corrected ``i60`` (QPE). Bottom: cumulative total, PrecipRate vs QPE.
    Takes one DataFrame from :func:`stormscape.mrms.virtual_gauge_timeseries`.
    A port of ``VirtualGage.plot_data`` in D. Cavagna's ``MRMS_stack``.
    """
    import matplotlib.pyplot as plt
    styles = {2: (":", "0.4"), 5: (":", "0.4"), 15: ("-", "red"),
              30: ("-", "blue"), 60: ("-", "teal")}
    cyc = ["red", "blue", "teal", "purple", "green"]
    fig, ax = plt.subplots(2, 1, figsize=figsize, sharex=True)
    for i, d in enumerate(durations):
        col = f"i{d}_mmph"
        if col not in df:
            continue
        ls, color = styles.get(d, ("-", cyc[i % len(cyc)]))
        ax[0].plot(df.index, df[col], ls, color=color,
                   lw=0.8 if ls == ":" else 1.0, label=f"I$_{{{d}}}$")
    if "i60_qpe_mmph" in df:
        m = df["i60_qpe_mmph"].notna()
        ax[0].plot(df.index[m], df["i60_qpe_mmph"][m], ".--", color="teal",
                   lw=1, label="I$_{60}$ (QPE)")
    ax[0].set_ylabel("I [mm h$^{-1}$]")
    ax[0].legend(fontsize=8)
    ax[0].grid(lw=0.5, alpha=0.5)

    ax[1].plot(df.index, df["total_mm"], "-k", lw=1.2, label="PrecipRate")
    if "total_qpe_mm" in df:
        m = df["total_qpe_mm"].notna()
        ax[1].plot(df.index[m], df["total_qpe_mm"][m], ".--k", lw=1, label="QPE")
    ax[1].set_ylabel("Total [mm]")
    ax[1].legend(title="Method", fontsize=8)
    ax[1].grid(lw=0.5, alpha=0.5)
    for lbl in ax[1].get_xticklabels():
        lbl.set_rotation(45)
        lbl.set_ha("right")
    fig.suptitle(name)
    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
    return out_path or (fig, ax)


def _overlay_vg_series(ax, sources, real_series, name, value, colors,
                       lw=1.0, real_lw=1.1):
    """Plot one gauge's ``value`` column on ``ax``: each radar source solid +
    coloured (atlas palette), the matching real gauge dashed black. Shared by
    :func:`virtual_gauge_atlas` and :func:`virtual_gauge_detail`. Returns the
    max y value drawn (for setting a non-degenerate ylim)."""
    ymax = 0.0
    for i, (_lbl, s) in enumerate(sources.items()):
        df = s.get(name)
        if df is not None and value in df and df[value].notna().any():
            ax.plot(df.index, df[value], "-", color=colors[i % len(colors)],
                    lw=lw)
            ymax = max(ymax, float(np.nanmax(df[value])))
    rs = (real_series or {}).get(name)
    if rs is not None and value in rs and rs[value].notna().any():
        ax.plot(rs.index, rs[value], "--", color="black", lw=real_lw)
        ymax = max(ymax, float(np.nanmax(rs[value])))
    return ymax


def virtual_gauge_atlas(sources, real_series=None, value="i15_mmph",
                        out_path=None, ncols=None, title=None, figsize=None,
                        dpi=200):
    """Atlas of virtual-gauge time series -- one small panel per gauge.

    ``sources`` is ``{label: {gauge_name: DataFrame}}`` for one or more radar
    sources (e.g. ``{"MRMS": ..., "NEXRAD": ...}``), each drawn as a solid
    coloured line; a plain ``{gauge_name: DataFrame}`` dict is also accepted (one
    source). Where a real gauge of the same name exists in ``real_series`` it is
    overlaid dashed black. Gauge panels are the union of names across sources +
    real. Inputs are the dicts from
    :func:`stormscape.mrms.virtual_gauge_timeseries` /
    :func:`stormscape.nexrad.virtual_gauge_timeseries` /
    :func:`stormscape.gauges.gauge_timeseries`.
    """
    import math

    import matplotlib.pyplot as plt
    if sources and not isinstance(next(iter(sources.values())), dict):
        sources = {"radar VG": sources}              # back-compat: single source
    names = []
    for s in (*sources.values(), real_series or {}):
        for nm in s:
            if nm not in names:
                names.append(nm)
    if not names:
        raise ValueError("no virtual-gauge series to plot")
    colors = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e"]
    n = len(names)
    ncols = ncols or min(5, max(1, int(round(n ** 0.5))))
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols, sharex=True, squeeze=False,
                             figsize=figsize or (3.1 * ncols, 2.2 * nrows))
    axes = axes.flatten()
    label = {"i5_mmph": "I$_5$", "i15_mmph": "I$_{15}$", "i30_mmph": "I$_{30}$",
             "i60_mmph": "I$_{60}$", "total_mm": "total"}.get(value, value)
    unit = "mm" if value == "total_mm" else "mm h$^{-1}$"
    for ax, name in zip(axes, names):
        ymax = _overlay_vg_series(ax, sources, real_series, name, value, colors)
        ax.set_ylim(0, max(ymax * 1.15, 1.0))        # avoid ~0 dry-panel autoscale
        ax.set_title(str(name)[:24], fontsize=7)
        ax.tick_params(labelsize=6)
        ax.grid(lw=0.3, alpha=0.4)
        for lab in ax.get_xticklabels():
            lab.set_rotation(45)
            lab.set_ha("right")
    for ax in axes[n:]:
        ax.axis("off")
    handles = [Line2D([], [], color=colors[i % len(colors)], lw=1.4, label=lbl)
               for i, lbl in enumerate(sources)]
    if real_series:
        handles.append(Line2D([], [], color="black", ls="--", lw=1.4,
                              label="gauge"))
    fig.legend(handles=handles, loc="upper right", fontsize=8)
    fig.suptitle(title or f"Virtual gauges — {label} [{unit}]", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])           # leave room for suptitle
    if out_path:
        fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
    return out_path or (fig, axes)


def _tempo_label(report_min):
    """Short label for a gauge's native reporting cadence in minutes (the
    ``report_min`` field from :func:`stormscape.gauges._report_min`); ``None`` when
    unknown, so callers can omit it. E.g. ``5 -> "~5-min reporting"``,
    ``60 -> "~hourly (60-min)"``, ``1440 -> "~daily (1440-min)"``."""
    try:
        m = float(report_min)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(m):
        return None
    if m >= 1440:
        return f"~daily ({m:g}-min)"
    if m == 60:
        return "~hourly (60-min)"
    return f"~{m:g}-min reporting"


def virtual_gauge_detail(sources, name, real_series=None, out_path=None,
                         figsize=(8.5, 10.0), dpi=200, title=None,
                         report_min=None):
    """Detailed single-gauge comparison — a 4-row vertical stack of cumulative
    rainfall, then I60 / I30 / I15 intensity, each overlaying every radar
    source (solid, same palette as :func:`virtual_gauge_atlas` — MRMS blue,
    NEXRAD red) and the matching real gauge (dashed black).

    ``sources`` is ``{label: {gauge_name: DataFrame}}`` (a bare
    ``{gauge_name: DataFrame}`` dict is also accepted); ``name`` selects which
    gauge to draw. The big-figure counterpart of one atlas panel.

    The ground gauge's native reporting **tempo** is annotated in the title:
    ``report_min`` (minutes) is taken from the explicit ``report_min`` argument,
    else the real series' ``df.attrs['report_min']`` (set by
    :func:`stormscape.gauges.fetch_gauge_event` / ``gauge_timeseries`` /
    ``load_event_series``). Coarse (hourly/daily) reporters smear bursts under the
    1-min interpolation -> their i15/i30 read low, so the tempo flags how far to
    trust the peaks. Omitted when unknown.
    """
    import matplotlib.pyplot as plt
    if sources and not isinstance(next(iter(sources.values())), dict):
        sources = {"radar VG": sources}              # back-compat: single source
    colors = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e"]
    panels = [("total_mm", "Cumulative rainfall", "mm"),
              ("i60_mmph", "I$_{60}$", "mm h$^{-1}$"),
              ("i30_mmph", "I$_{30}$", "mm h$^{-1}$"),
              ("i15_mmph", "I$_{15}$", "mm h$^{-1}$")]
    fig, axes = plt.subplots(4, 1, sharex=True, figsize=figsize)
    for ax, (value, lab, unit) in zip(axes, panels):
        ymax = _overlay_vg_series(ax, sources, real_series, name, value,
                                  colors, lw=1.8, real_lw=2.0)
        ax.set_ylim(0, max(ymax * 1.1, 1.0))
        ax.set_ylabel(f"{lab} [{unit}]", fontsize=10)
        ax.grid(lw=0.3, alpha=0.4)
        ax.tick_params(labelsize=8)
    axes[-1].set_xlabel("time (UTC)", fontsize=10)
    for t in axes[-1].get_xticklabels():
        t.set_rotation(30)
        t.set_ha("right")
    handles = [Line2D([], [], color=colors[i % len(colors)], lw=2.0, label=lbl)
               for i, lbl in enumerate(sources)]
    if (real_series or {}).get(name) is not None:
        handles.append(Line2D([], [], color="black", ls="--", lw=2.0,
                              label="gauge"))
    axes[0].legend(handles=handles, loc="upper left", fontsize=9)
    if report_min is None and real_series and name in real_series:
        report_min = (real_series[name].attrs or {}).get("report_min")
    tempo = _tempo_label(report_min)
    base = f"{name} — virtual gauge vs gauge"
    fig.suptitle(title or (f"{base}   (ground gauge: {tempo})" if tempo else base),
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    if out_path:
        fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
    return out_path or (fig, axes)


_PICKER_HTML = r"""<!doctype html>
<html><head><meta charset="utf-8"><title>__TITLE__</title>
<style>
 body{font-family:system-ui,Arial,sans-serif;margin:0;background:#1e1e1e;color:#eee}
 header{padding:10px 14px;background:#2a2a2a;border-bottom:1px solid #444}
 header b{color:#fff} .hint{color:#aaa;font-size:13px}
 #wrap{position:relative;display:inline-block;margin:12px;line-height:0}
 #map{max-width:70vw;height:auto;display:block;user-select:none;-webkit-user-drag:none}
 #ov{position:absolute;left:0;top:0;cursor:crosshair}
 #panel{position:fixed;right:14px;top:64px;width:380px;background:#2a2a2a;
        border:1px solid #444;border-radius:8px;padding:14px;font-size:13px}
 .row{margin:10px 0} .lab{color:#9cf;font-weight:600;margin-bottom:3px}
 code{background:#111;color:#7fd;padding:6px 8px;border-radius:4px;display:block;
      white-space:pre-wrap;word-break:break-all;line-height:1.45;font-size:12px}
 button{background:#3a6;color:#fff;border:0;border-radius:5px;padding:6px 10px;
        cursor:pointer;font-size:12px;margin-top:5px} button:hover{background:#4b7}
 .muted{color:#999;font-size:12px}
</style></head>
<body>
<header><b>stormscape — pick a zoom bbox</b>
 <span class="hint">&nbsp;drag a rectangle on the map, then copy the command (or just the bbox).</span></header>
<div id="wrap"><img id="map" src="data:image/png;base64,__IMG__">
 <canvas id="ov"></canvas></div>
<div id="panel">
 <div class="row muted">cursor (lon, lat): <span id="cur">—</span></div>
 <div class="row"><div class="lab">--bbox  W S E N</div>
  <code id="bbox">drag to select…</code>
  <button onclick="cp('bbox')">copy bbox</button></div>
 <div class="row"><div class="lab">zoom command</div>
  <code id="cmd">—</code>
  <button onclick="cp('cmd')">copy command</button></div>
 <div class="row muted">tip: lat/long ticks line the map edges; the box maps to exact lon/lat.</div>
</div>
<script>
 const PRE=__PREFIX_JSON__, SUF=__SUFFIX_JSON__;
 const AX={l:__AXL__,r:__AXR__,t:__AXT__,b:__AXB__};            // map axes, image fraction
 const C={tl:[__LONTL__,__LATTL__],tr:[__LONTR__,__LATTR__],     // corner lon/lat
          bl:[__LONBL__,__LATBL__],br:[__LONBR__,__LATBR__]};
 const img=document.getElementById('map'), ov=document.getElementById('ov');
 const ctx=ov.getContext('2d'); let drag=null, box=null;
 function fit(){ ov.width=img.clientWidth; ov.height=img.clientHeight;
   ov.style.width=img.clientWidth+'px'; ov.style.height=img.clientHeight+'px'; redraw(); }
 img.onload=fit; window.addEventListener('resize',fit); if(img.complete) fit();
 function pos(e){ const r=ov.getBoundingClientRect();
   return {x:Math.min(Math.max(e.clientX-r.left,0),r.width),
           y:Math.min(Math.max(e.clientY-r.top,0),r.height)}; }
 function ll(x,y){ const w=ov.clientWidth,h=ov.clientHeight;       // -> {lon,lat}
   let fx=((x/w)-AX.l)/(AX.r-AX.l), fy=((y/h)-AX.t)/(AX.b-AX.t);
   fx=Math.min(Math.max(fx,0),1); fy=Math.min(Math.max(fy,0),1);   // bilinear over corners
   const B=p=>(1-fx)*(1-fy)*C.tl[p]+fx*(1-fy)*C.tr[p]+(1-fx)*fy*C.bl[p]+fx*fy*C.br[p];
   return {lon:B(0), lat:B(1)}; }
 ov.addEventListener('mousedown',e=>{drag=pos(e); box=null;});
 ov.addEventListener('mousemove',e=>{ const p=pos(e), g=ll(p.x,p.y);
   document.getElementById('cur').textContent=g.lon.toFixed(4)+', '+g.lat.toFixed(4);
   if(drag){ box={a:drag,b:p}; redraw(); } });
 window.addEventListener('mouseup',e=>{ if(drag){ box={a:drag,b:pos(e)}; drag=null; redraw(); update(); } });
 function redraw(){ ctx.clearRect(0,0,ov.width,ov.height); if(!box) return;
   const x=Math.min(box.a.x,box.b.x), y=Math.min(box.a.y,box.b.y),
         w=Math.abs(box.a.x-box.b.x), h=Math.abs(box.a.y-box.b.y);
   ctx.fillStyle='rgba(0,200,255,0.15)'; ctx.fillRect(x,y,w,h);
   ctx.strokeStyle='#0cf'; ctx.lineWidth=2; ctx.strokeRect(x,y,w,h); }
 function update(){ if(!box) return;                              // 4 corners -> min/max
   const pts=[[box.a.x,box.a.y],[box.b.x,box.a.y],[box.a.x,box.b.y],[box.b.x,box.b.y]]
             .map(c=>ll(c[0],c[1]));
   const lo=pts.map(p=>p.lon), la=pts.map(p=>p.lat);
   const w=Math.min(...lo).toFixed(4), e=Math.max(...lo).toFixed(4),
         s=Math.min(...la).toFixed(4), n=Math.max(...la).toFixed(4);
   const bb=w+' '+s+' '+e+' '+n;
   document.getElementById('bbox').textContent=bb;
   document.getElementById('cmd').textContent=(PRE+' --bbox '+bb+(SUF?(' '+SUF):'')).trim(); }
 function cp(id){ const t=document.getElementById(id).textContent;
   navigator.clipboard&&navigator.clipboard.writeText(t); }
</script>
</body></html>
"""


def bbox_picker(i15, hillshade=None, out_path=None, cmap="YlGnBu", wet_min=5.0,
                reference=False, local_roads=False, label_reference=True,
                perimeters=None, gauges=None, cmd_prefix="", cmd_suffix="",
                title="stormscape bbox picker", dpi=160, max_px=1700):
    """Write a self-contained HTML bbox picker for choosing a zoom sub-AOI.

    Renders the event's i15-over-hillshade through :func:`drape_i15` (the **full
    main-map context** -- labelled NHD streams / TIGER roads / GNIS places, AOI
    perimeter, gauges, north arrow, lat/long ticks, colorbar -- in the same UTM
    projection as the production maps), embeds it as a base64 PNG, and overlays a
    vanilla-JS canvas. The map *axes* is a known rectangle of the figure with
    known corner coordinates, so a dragged rectangle maps back to an accurate
    ``--bbox W S E N`` (bilinear over the corners; exact for an AOI this size)
    plus a ready-to-run ``zoom`` command. No server, no GUI toolkit, no internet
    -> opens in any browser on Windows/macOS/Linux; rendering is headless
    (Agg/PIL), so generating it never needs a display. Returns the HTML string
    (written to ``out_path`` if given).
    """
    import base64
    import io
    import json

    import matplotlib.pyplot as plt
    from PIL import Image
    from pyproj import Transformer

    import geopandas as gpd

    from .aoi import bbox_polygon
    ref = _load(hillshade if hillshade is not None else i15)
    wc = _utm_crs_for(ref)                              # near north-up UTM zone
    ib = tuple(float(v) for v in _load(i15).rio.bounds())   # i15 extent (4326)
    clip_gdf = gpd.GeoDataFrame(geometry=[bbox_polygon(ib)], crs="EPSG:4326")
    cb = clip_gdf.to_crs(wc).total_bounds              # AOI extent in UTM metres
    aspect = (cb[2] - cb[0]) / (cb[3] - cb[1]) if (cb[3] - cb[1]) else 1.0
    ph = 7.6
    # frame tightly to the AOI (like the production maps) so the Albers->UTM
    # hillshade rotation doesn't show as gray wedges around the edges.
    fig, ax = drape_i15(
        hillshade, i15, out_path=None, work_crs=wc, cmap=cmap, wet_min=wet_min,
        alpha=None, reference=reference, local_roads=local_roads,
        label_reference=label_reference, perimeters=perimeters, gauges=gauges,
        gauge_value=("i15_mmph" if gauges is not None else None),
        legend=("gauges" if gauges is not None else "all"), title=title,
        north_arrow=True, scale_ticks=True, clip=clip_gdf, clip_margin=0.02,
        figsize=(ph * aspect + 1.9, ph), dpi=dpi)
    fig.canvas.draw()                                  # finalise layout

    pos = ax.get_position()                            # map axes, figure fraction
    xa, xb = ax.get_xlim()
    ya, yb = ax.get_ylim()                             # UTM metres; yb = top (north)
    t = Transformer.from_crs(wc, "EPSG:4326", always_xy=True)
    (lontl, lattl), (lontr, lattr) = t.transform(xa, yb), t.transform(xb, yb)
    (lonbl, latbl), (lonbr, latbr) = t.transform(xa, ya), t.transform(xb, ya)
    axl, axr, axt, axb = pos.x0, pos.x1, 1 - pos.y1, 1 - pos.y0  # image fraction

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi)            # NOT bbox_inches=tight (keep pos)
    plt.close(fig)
    im = Image.open(io.BytesIO(buf.getvalue()))
    if max(im.size) > max_px:                          # keep the HTML lean
        im = im.convert("RGB").resize(
            (round(im.width * max_px / max(im.size)),
             round(im.height * max_px / max(im.size))))
    obuf = io.BytesIO()
    im.save(obuf, format="PNG")
    b64 = base64.b64encode(obuf.getvalue()).decode()

    repl = {"__IMG__": b64, "__TITLE__": title,
            "__PREFIX_JSON__": json.dumps(cmd_prefix),
            "__SUFFIX_JSON__": json.dumps(cmd_suffix),
            "__AXL__": f"{axl:.5f}", "__AXR__": f"{axr:.5f}",
            "__AXT__": f"{axt:.5f}", "__AXB__": f"{axb:.5f}",
            "__LONTL__": f"{lontl:.6f}", "__LATTL__": f"{lattl:.6f}",
            "__LONTR__": f"{lontr:.6f}", "__LATTR__": f"{lattr:.6f}",
            "__LONBL__": f"{lonbl:.6f}", "__LATBL__": f"{latbl:.6f}",
            "__LONBR__": f"{lonbr:.6f}", "__LATBR__": f"{latbr:.6f}"}
    html = _PICKER_HTML
    for k, v in repl.items():
        html = html.replace(k, v)
    if out_path:
        with open(out_path, "w") as f:
            f.write(html)
    return html
