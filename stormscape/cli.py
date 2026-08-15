"""Command-line interface: ``python -m stormscape <command>``.

Commands
--------
dem      download a DEM + hillshade for an AOI
i15      stack MRMS into peak-i15/i30/i60 fields for an AOI + storm-day
map      drape an existing i15 GeoTIFF over an existing hillshade GeoTIFF
run      DEM -> hillshade -> i15 -> figure in one shot (+ --gauges/--compare)
gauges   fetch Synoptic rain-gauge total + peak 15/30/60-min intensities
compare  sample radar rasters at gauges -> residuals + skill stats (+ map)
nexrad   single-radar NEXRAD Level II reflectivity tilt -> GeoTIFF + figure
panels   multi-panel diagnostic map (time-of-peak i15, QPE total, RQI, SHSR)
vgauge   virtual rain gauges: radar rainfall time series at point(s)
zoom     re-render an existing event's figures clipped to a sub-AOI (no re-run)
pick     interactive browser bbox picker: drag a zoom box on an event's map
climate  NOAA Atlas 14 climatology vs observed i15/i30/i60 + anomaly maps
export   georeferenced exports for GIS/CalTopo: EPSG:3857 GeoTIFFs + GeoPDFs

The AOI is given as either ``--bbox W S E N`` (lon/lat degrees) or
``--aoi path`` (any vector file GeoPandas can read).
"""

from __future__ import annotations

import argparse
import json
import os
import sys

from .export import DEFAULT_EXPORT_FIELDS
from .layout import find, find_subdir, out_path, subdir


def _aoi_from_args(args):
    if args.bbox:
        return tuple(args.bbox)
    if args.aoi:
        return args.aoi
    sys.exit("error: provide --bbox W S E N or --aoi <vector file>")


def _find_event_aoi(from_dir, from_key):
    """Locate a processed event's AOI vector, for auto-matching a reused figure's
    extent to the i15 maps (which clip to the AOI, not the wider padded radar
    footprint). Checks ``<from_key>_aoi.geojson`` (saved by dem/i15/run/nexrad) then
    a user-placed ``<from_key>_AOI.{kmz,geojson,gpkg,shp}``. Returns a path or None."""
    if not (from_dir and from_key):
        return None
    names = [f"{from_key}_aoi.geojson"] + [
        f"{from_key}_AOI.{e}" for e in ("kmz", "geojson", "gpkg", "shp")]
    for n in names:
        p = find(from_dir, n)
        if os.path.exists(p):
            return p
    return None


def _save_event_aoi(args, out_dir, key):
    """Persist the resolved AOI as ``<key>_aoi.geojson`` next to the event rasters,
    so later commands (climate, ...) can auto-match the figure extent to the i15
    maps. Best-effort -- warns but never fails the run."""
    if not (getattr(args, "bbox", None) or getattr(args, "aoi", None)):
        return None
    try:
        import geopandas as gpd
        from .aoi import bbox_polygon, load_aoi
        bounds, geom = load_aoi(_aoi_from_args(args))
        g = gpd.GeoDataFrame(geometry=[geom or bbox_polygon(bounds)], crs=4326)
        path = out_path(out_dir, f"{key}_aoi.geojson")
        g.to_file(path, driver="GeoJSON")
        return path
    except Exception as e:                           # noqa: BLE001 (best-effort)
        print(f"note: could not save event AOI ({e})")
        return None


def _labels(args):
    """Normalise --basemap-labels ('none'/'off'/empty -> None)."""
    v = getattr(args, "basemap_labels", None)
    if v is None or str(v).strip().lower() in ("", "none", "off"):
        return None
    return v


def _cmd_dem(args):
    from .dem import fetch_dem_and_hillshade
    os.makedirs(args.out_dir, exist_ok=True)
    key = args.key or "aoi"
    dem_path = out_path(args.out_dir, f"{key}_dem.tif")
    hs_path = out_path(args.out_dir, f"{key}_hillshade.tif")
    fetch_dem_and_hillshade(_aoi_from_args(args), resolution=args.resolution,
                            dst_crs=args.dst_crs, dem_path=dem_path,
                            hillshade_path=hs_path, pad_deg=args.pad_deg,
                            clip=args.clip_dem,
                            resampling=getattr(args, "resampling", None))
    print(f"wrote {dem_path}\nwrote {hs_path}")
    return dem_path, hs_path


def _cmd_i15(args):
    from .mrms import i15_storm_day, save_fields
    os.makedirs(args.out_dir, exist_ok=True)
    key = args.key or f"i15_{_event_label(args)}"
    _save_event_aoi(args, args.out_dir, key)        # for climate auto-extent-match
    # --start/--end scope the RADAR STACK too, not just the gauges. A --date
    # window spans ~30 h so the local day is covered, which silently merges
    # back-to-back evening storms: on 2026-08-14 the 04Z hour was the tail of
    # the 13 Aug storm and landed in "today's" peak maps.
    window = _stack_window(args)
    res = i15_storm_day(_aoi_from_args(args), args.date, pad_deg=args.pad_deg,
                        qpe_thresh=args.qpe_thresh,
                        max_wet_hours=args.max_wet_hours, workers=args.workers,
                        window=window)
    paths = save_fields(res, args.out_dir, key)
    print(json.dumps(res["meta"], indent=2))
    for p in paths:
        print(f"wrote {p}")
    if getattr(args, "multisensor", False):
        from .mrms import multisensor_total
        ms = multisensor_total(_aoi_from_args(args), args.date,
                               pad_deg=args.pad_deg, workers=args.workers,
                               window=window)
        for p in save_fields(ms, args.out_dir, key):
            print(f"wrote {p}")
    return out_path(args.out_dir, f"{key}_i15max.tif")


def _cmd_map(args):
    from .plot import drape_i15
    out = args.out or out_path(args.out_dir or ".",
                               (args.key or "i15_map") + ".png")
    drape_i15(args.hillshade, args.i15, out_path=out, work_crs="UTM",
              cmap=args.cmap, wet_min=args.wet_min, perimeters=args.perimeters,
              basins=args.basins, highlight=args.highlight, points=args.points,
              title=args.title, basemap=args.basemap,
              basemap_provider=args.basemap_provider,
              basemap_labels=_labels(args), basemap_zoom=args.basemap_zoom,
              hillshade_alpha=args.hillshade_alpha,
              alpha=args.alpha,
              reference=args.reference, local_roads=args.local_roads,
              label_reference=not args.no_reference_labels,
              north_arrow=True, scale_ticks=True,
              clip=(args.perimeters if args.clip else None),
              clip_margin=args.clip_margin, dpi=args.dpi)
    print(f"wrote {out}")
    return out


def _cmd_run(args):
    _, hs_path = _cmd_dem(args)
    i15_path = _cmd_i15(args)
    from .plot import drape_i15
    key = args.key or "run"
    _save_event_aoi(args, args.out_dir, key)        # for climate auto-extent-match
    out = out_path(args.out_dir, f"{key}.png")
    title = args.title or f"{key}  -  i15 over terrain  ({_event_label(args)})"
    clip = (args.perimeters or args.aoi) if args.clip else None

    # optional rain gauges (--gauges) and radar-vs-gauge comparison (--compare,
    # which implies the fetch). Gauges are overlaid on the main figure coloured
    # on the i15 scale; --compare also writes a CSV + a residual map.
    gauges_gdf = None
    if args.gauges or args.compare:
        g, start, end = _fetch_gauges(args)
        if len(g):
            gpath = out_path(args.out_dir, f"{key}_gauges.geojson")
            g.to_file(gpath, driver="GeoJSON")
            print(f"wrote {gpath}  ({len(g)} gauges)")
            gauges_gdf = g
        else:
            print(f"no gauges found in AOI for {start} to {end}")

    drape_i15(hs_path, i15_path, out_path=out, work_crs="UTM",
              cmap=args.cmap, wet_min=args.wet_min, perimeters=args.perimeters,
              basins=args.basins, highlight=args.highlight, points=args.points,
              gauges=gauges_gdf,
              gauge_value=("i15_mmph" if gauges_gdf is not None else None),
              title=title, basemap=args.basemap,
              basemap_provider=args.basemap_provider,
              basemap_labels=_labels(args), basemap_zoom=args.basemap_zoom,
              hillshade_alpha=args.hillshade_alpha,
              alpha=args.alpha,
              reference=args.reference, local_roads=args.local_roads,
              label_reference=not args.no_reference_labels,
              north_arrow=True, scale_ticks=True,
              legend=("gauges" if gauges_gdf is not None else "all"),
              clip=clip, clip_margin=args.clip_margin, dpi=args.dpi)
    print(f"wrote {out}")

    if args.compare and gauges_gdf is not None:
        from .compare import compare_storm
        table, stats = compare_storm(gauges_gdf, args.out_dir, key,
                                     rqi_min=args.rqi_min,
                                     max_report_min=args.max_report_min,
                                     multisensor=args.multisensor)
        csv = out_path(args.out_dir, f"{key}_compare.csv")
        table.drop(columns="geometry").to_csv(csv, index=False)
        print(stats.to_string(index=False))
        print(f"wrote {csv}")
        cmp_png = out_path(args.out_dir, f"{key}_compare.png")
        drape_i15(hs_path, i15_path, out_path=cmp_png, work_crs="UTM",
                  cmap=args.cmap, wet_min=args.wet_min,
                  perimeters=args.perimeters,
                  gauges=table, gauge_value="resid_i15max",
                  reference=args.reference, local_roads=args.local_roads,
                  label_reference=not args.no_reference_labels,
                  alpha=args.alpha, north_arrow=True, scale_ticks=True,
                  legend="gauges",
                  clip=clip, clip_margin=args.clip_margin, dpi=args.dpi,
                  title=f"{key}  -  radar - gauge  i15  ({_event_label(args)})")
        print(f"wrote {cmp_png}")


def _gauge_window(args):
    """(start, end) UTC datetimes from --date (local-day scan) or --start/--end."""
    import datetime as dt

    from .mrms import SCAN_PAD_H, parse_date
    if getattr(args, "start", None) and getattr(args, "end", None):
        f = "%Y%m%d%H%M"
        return (dt.datetime.strptime(args.start, f),
                dt.datetime.strptime(args.end, f))
    if not getattr(args, "date", None):
        sys.exit("error: provide --date YYYYMMDD or both --start and --end")
    d = parse_date(args.date)
    start = dt.datetime(d.year, d.month, d.day, SCAN_PAD_H[0])
    end = (dt.datetime(d.year, d.month, d.day)
           + dt.timedelta(days=1, hours=SCAN_PAD_H[1]))
    return start, end


def _event_label(args):
    """``YYYYMMDD`` for keys/titles, from --date or from the window start."""
    if getattr(args, "date", None):
        return str(args.date).replace("-", "")
    win = _stack_window(args)
    return f"{win[0]:%Y%m%d}" if win else "event"


def _stack_window(args):
    """``(start, end)`` for the MRMS stack when the user pinned one, else None.

    Deliberately the SAME ``--start``/``--end`` the gauge side uses: having one
    pair of flags mean "the analysis window" everywhere is less surprising than
    having them scope the gauges while the radar quietly stacks a whole day.
    """
    import datetime as dt

    if not (getattr(args, "start", None) and getattr(args, "end", None)):
        if not getattr(args, "date", None):
            sys.exit("error: provide --date YYYYMMDD, or both --start and --end")
        return None
    f = "%Y%m%d%H%M"
    return (dt.datetime.strptime(args.start, f), dt.datetime.strptime(args.end, f))


def _explicit_window(args):
    """True when the user pinned the window with --start/--end (never narrow it)."""
    return bool(getattr(args, "start", None) and getattr(args, "end", None))


def _narrow_to_storm(args, start, end, series=None, reason="NEXRAD"):
    """Shrink a storm-DAY window to the storm's wet period *before* an expensive
    fetch, using only cheap products.

    A ``--date`` window spans ~30 h (:data:`stormscape.mrms.SCAN_PAD_H`) so the
    local day is fully covered, but a NEXRAD Level II fetch over that window
    pulls ~300 volumes / ~2 GB when the storm itself lasted a few hours. Order of
    preference:

    1. an explicit ``--start``/``--end`` -- the user's window always wins;
    2. the gauge series' mass-weighted :func:`stormscape.gauges.storm_window`;
    3. the radar-side :func:`stormscape.mrms.wet_window` (hourly QPE, a few KB
       per hour) -- which also covers the no-gauge and dry-gauge cases, where
       (2) returns ``None`` yet the radar clearly saw rain.

    Falls back to the original window if everything reads dry, so a genuinely
    unknown storm is never silently truncated to nothing.
    """
    if _explicit_window(args):
        return start, end
    span_h = (end - start).total_seconds() / 3600.0
    wins = []
    if series:
        from .gauges import storm_window
        win = storm_window(series, pad_min=getattr(args, "pad_min", 30))
        if win:
            wins.append(("gauge rain mass", win))
    try:
        from .mrms import wet_window
        win = wet_window(_aoi_from_args(args), start, end,
                         pad_deg=getattr(args, "pad_deg", 0.05),
                         pad_min=getattr(args, "pad_min", 30))
        if win:
            wins.append(("MRMS hourly QPE", win))
    except Exception as e:                                     # noqa: BLE001
        print(f"  note: MRMS wet-window probe failed ({repr(e)[:60]})")
    if not wins:
        print(f"  note: no wet period found in the {span_h:.0f} h window; keeping it whole")
        return start, end
    # UNION, not first-match. Gauges are point samples: on 2026-08-13 a single
    # wet gauge put the rain mass in 19:01-20:29 while MRMS saw the AOI wet
    # 19-04Z, and trusting the gauge window alone would have cut the storm's
    # peak and its second cell out of the comparison entirely. The radar bounds
    # when it rained *anywhere in the AOI*; the gauges can only extend that.
    lo = min(w[0] for _, w in wins)
    hi = max(w[1] for _, w in wins)
    src = " + ".join(n for n, _ in wins)
    print(f"  {reason} window from {src}: {lo:%m-%d %H:%M}-{hi:%m-%d %H:%M} "
          f"({(hi - lo).total_seconds()/3600:.1f} h of {span_h:.0f} h)")
    return lo, hi


def _fetch_gauges(args):
    from .gauges import SynopticError, gauge_fields
    start, end = _gauge_window(args)
    try:
        g = gauge_fields(_aoi_from_args(args), start, end, token=args.token,
                         durations=tuple(args.durations), pad_deg=args.pad_deg)
    except SynopticError as e:
        sys.exit(f"Synoptic error: {e}")
    return g, start, end


def _cmd_gauges(args):
    """Unified gauge pipeline. In ONE Synoptic draw: (1) build the canonical gauge
    store -- ``<key>_gauges.geojson`` (coords + peak metrics + I15 time-of-peak) +
    per-gauge series CSVs in ``RainGaugeData/`` (reused by compare / recurrence /
    vgauge) -- then, by default, (2) drop virtual gauges at every wet near-AOI
    station from MRMS (+NEXRAD, 3-way) and render the rainfall comparison atlas and
    (3) the per-gauge detail figures (``VirtualGaugeFigures/``). The store keeps the
    full storm-DAY record; the products use the tight rain window. The atlas and the
    detail figures share one filtered (wet) station set, so their gauges are uniform.
    ``--store-only`` reverts to just the store."""
    import pandas as pd
    from .gauges import SynopticError, fetch_gauge_event
    os.makedirs(args.out_dir, exist_ok=True)
    key = args.key or "gauges"
    start, end = _gauge_window(args)
    try:
        g, series = fetch_gauge_event(
            _aoi_from_args(args), start, end, args.out_dir, key,
            token=args.token, durations=tuple(args.durations),
            pad_deg=args.pad_deg, write_series=not args.no_series)
    except SynopticError as e:
        sys.exit(f"Synoptic error: {e}")
    out = out_path(args.out_dir, f"{key}_gauges.geojson")
    if not len(g):
        print(f"no gauges found in AOI for {start} to {end}")
        return out
    wet = int((pd.to_numeric(g.get("i15_mmph"), errors="coerce").fillna(0)
               > 0).sum()) if "i15_mmph" in g.columns else 0
    print(f"wrote {out}  ({len(g)} gauges, {wet} wet; "
          f"{len(series)} series -> {args.out_dir}/RainGaugeData/; "
          f"{start:%Y-%m-%d %HZ} to {end:%m-%d %HZ})")

    # default: also build the virtual-gauge comparison products (atlas + detail)
    # from the SAME draw -- the full pipeline in one command. Clip the full-day
    # store to the storm's rain window, then filter to wet (default) near-AOI
    # stations; that one station set feeds both the atlas and the detail figures.
    if getattr(args, "store_only", False) or not series:
        return out
    # Bound the products to the storm. Uses gauge rain mass UNIONed with the MRMS
    # wet window -- a sparse or unlucky gauge set can otherwise put the window on
    # a fraction of a storm the radar saw across the whole AOI.
    vstart = min(d.index.min() for d in series.values())
    vend = max(d.index.max() for d in series.values())
    vstart, vend = _narrow_to_storm(args, vstart, vend, series=series)
    real = {}
    for n, df in series.items():                # full storm-day -> rain window
        d = df.loc[vstart:vend]
        d.attrs = dict(df.attrs)
        real[n] = d
    real = _filter_stations(real, args)         # --wet-only (default) / --max-dist-km
    if not real:
        print("note: no stations left after filtering; skipping atlas/detail")
        return out
    points = [(n, df.attrs["lon"], df.attrs["lat"]) for n, df in real.items()]
    print(f"virtual gauges over rain window {vstart:%m-%d %H:%M}-{vend:%m-%d %H:%M} "
          f"({len(real)} stations; sources: MRMS"
          + ("+NEXRAD" if getattr(args, 'nexrad', False) else "") + "+real)")
    _virtual_gauge_products(
        args, points=points, real=real, reusing=True, user_points=[],
        start=vstart, end=vend, key=key,
        make_atlas=not args.no_atlas, make_detail=not args.no_detail)
    return out


def _cmd_compare(args):
    from .compare import compare_storm
    radar_dir = args.radar_dir or args.out_dir or "."
    key = args.key or "i15"
    if args.gauges:
        import geopandas as gpd
        gauges = gpd.read_file(args.gauges)
    else:
        gauges, _, _ = _fetch_gauges(args)
    if not len(gauges):
        sys.exit("no gauges to compare")
    table, stats = compare_storm(gauges, radar_dir, key, rqi_min=args.rqi_min,
                                 max_report_min=args.max_report_min,
                                 multisensor=args.multisensor)
    out_csv = args.out or out_path(args.out_dir or ".", f"{key}_compare.csv")
    os.makedirs(os.path.dirname(os.path.abspath(out_csv)), exist_ok=True)
    table.drop(columns="geometry").to_csv(out_csv, index=False)
    print(stats.to_string(index=False))
    print(f"wrote {out_csv}")
    if args.map:
        from .plot import drape_i15
        i15_tif = find(radar_dir, f"{key}_i15max.tif")
        drape_i15(args.hillshade, i15_tif, out_path=args.map,
                  work_crs="UTM", cmap=args.cmap, gauges=table,
                  gauge_value=f"resid_{args.map_metric}",
                  perimeters=args.perimeters,
                  clip=(args.perimeters if args.clip else None),
                  reference=args.reference, alpha=args.alpha,
                  north_arrow=True, scale_ticks=True, legend="gauges",
                  title=args.title or f"radar - gauge  ({args.map_metric})")
        print(f"wrote {args.map}")
    return out_csv


def _cmd_panels(args):
    from .plot import diagnostic_panels
    out = args.out or out_path(args.radar_dir or ".", f"{args.key}_panels.png")
    diagnostic_panels(args.radar_dir or ".", args.key, which=args.fields,
                      out_path=out, hillshade=args.hillshade,
                      perimeters=args.perimeters, gauges=args.gauges,
                      reference=args.reference, local_roads=args.local_roads,
                      label_reference=not args.no_reference_labels,
                      alpha=args.alpha,
                      clip=(args.perimeters if args.clip else None),
                      clip_margin=args.clip_margin, title=args.title, dpi=args.dpi)
    print(f"wrote {out}")
    return out


def _crop_rasters(src_dir, src_key, out_dir, key, bounds, pad_deg, skip=()):
    """Write cropped copies of an event's GeoTIFFs to a new folder (a
    self-contained zoom). CRS-aware: the lon/lat zoom box is reprojected to each
    raster's own CRS before clipping (MRMS fields are 4326, terrain is 5070)."""
    import geopandas as gpd
    import rioxarray  # noqa: F401  (registers the .rio accessor)
    from rioxarray import open_rasterio
    from shapely.geometry import box

    from .aoi import pad_bounds
    w, s, e, n = pad_bounds(bounds, pad_deg)
    fields = ("dem", "hillshade", "i15max", "i30max", "i60max", "i2max",
              "total", "tpki15", "rqi", "shsr", "peakrate_mmph", "mstotal")
    wrote = 0
    for f in fields:
        if f in skip:
            continue
        p = os.path.join(src_dir, f"{src_key}_{f}.tif")
        if not os.path.exists(p):
            continue
        with open_rasterio(p) as da:               # close the handle promptly
            bb = gpd.GeoSeries([box(w, s, e, n)], crs=4326).to_crs(da.rio.crs)
            try:
                cropped = da.rio.clip_box(*bb.total_bounds).load()
            except Exception as ex:                # zoom box outside raster, etc.
                print(f"note: could not crop {f}: {ex}")
                continue
        cropped.rio.to_raster(out_path(out_dir, f"{key}_{f}.tif"))
        wrote += 1
    print(f"cropped {wrote} rasters -> {out_dir}")


def _cmd_zoom(args):
    """Re-render an already-processed event's figures clipped to a sub-AOI,
    reusing its rasters (no MRMS re-download -- MRMS has no finer resolution).
    Only --refine-dem re-fetches, and only the terrain (the one product that
    benefits from zooming)."""
    import geopandas as gpd

    from .aoi import bbox_polygon, load_aoi
    from .plot import _prepare_hillshade, diagnostic_panels, drape_i15
    os.makedirs(args.out_dir, exist_ok=True)
    key = args.key or f"{args.from_key}_zoom"
    src = args.from_dir

    # zoom sub-AOI -> a clip GeoDataFrame (EPSG:4326) the plot fns view-clip to
    spec = _aoi_from_args(args)
    bounds, geom = load_aoi(spec)
    clip_gdf = gpd.GeoDataFrame(geometry=[geom or bbox_polygon(bounds)], crs=4326)

    # one drape opacity for every zoom figure (map + panels); None falls back to
    # the project-wide DEFAULT_FIELD_ALPHA (0.32) in plot.py, --alpha overrides.
    zoom_alpha = args.alpha
    # smooth the radar fields for display by default (Gaussian 1 km) -- the same
    # --obs-smooth knob the climate maps use, applied to the i15 map + rainfall
    # panels too so every rainfall map in the zoom folder is smoothed alike
    osm = None if args.obs_smooth == "none" else args.obs_smooth
    osr = args.obs_smooth_radius
    sm_note = f", {osm} {osr:g} km" if (osm and osr > 0) else ""

    # terrain backdrop: reuse the source hillshade, or re-fetch it finer for the
    # zoom extent (3DEP has 3 m / 1 m tiers; MRMS does not -> reuse those). The
    # 3DEP 1 m WMS can be slow/flaky, so retry; if it still fails, fall back to
    # the existing hillshade so the zoom still renders (just at coarser terrain).
    refined = False
    src_hs = find(src, f"{args.from_key}_hillshade.tif")
    if args.refine_dem:
        from .dem import fetch_dem_and_hillshade
        dem_path = out_path(args.out_dir, f"{key}_dem.tif")
        hs_path = out_path(args.out_dir, f"{key}_hillshade.tif")
        # cap the DEM pad small: the view is clipped to the box, so the hillshade
        # only needs a sliver of edge padding -- and at 1 m the default ~0.05 deg
        # (~5.5 km) pad would ~5x the area and blow past the WMS 120 s timeout.
        refine_pad = min(args.pad_deg, 0.003)
        print(f"fetching {args.resolution} m DEM for the zoom extent "
              f"(3DEP can be slow at 1 m; will retry, then fall back)...")
        try:
            fetch_dem_and_hillshade(spec, resolution=args.resolution,
                                    dst_crs=args.dst_crs, dem_path=dem_path,
                                    hillshade_path=hs_path, pad_deg=refine_pad,
                                    clip=args.clip_dem, retries=2,
                                    resampling=getattr(args, "resampling", None))
            print(f"wrote {dem_path}\nwrote {hs_path}")
            refined = True
        except Exception as ex:                # noqa: BLE001 (3DEP timeout, etc.)
            print(f"warning: --refine-dem {args.resolution} m fetch failed "
                  f"({type(ex).__name__}); falling back to the existing "
                  f"{args.from_key} hillshade so the zoom still renders.")
            hs_path = src_hs if os.path.exists(src_hs) else None
    else:
        hs_path = src_hs if os.path.exists(src_hs) else None

    if args.crop_rasters:
        # if the finer DEM/hillshade were written, don't overwrite them with
        # cropped source copies; otherwise crop the source terrain too.
        _crop_rasters(src, args.from_key, args.out_dir, key, bounds,
                      args.pad_deg,
                      skip=("dem", "hillshade") if refined else ())

    gpath = find(src, f"{args.from_key}_gauges.geojson")
    have_gauges = args.gauges and os.path.exists(gpath)
    if args.gauges and not have_gauges:
        print(f"note: {gpath} not found; skipping gauge overlay")

    # prep the terrain ONCE: reproject to UTM + downsample to render resolution,
    # then reuse the same array across the map, panels, and climate figures. The
    # 1 m zoom hillshade is ~170 M cells (~25x a 200-300 dpi figure can show), so
    # rendering it raw per figure is the dominant cost; this is invisible in the
    # output and leaves the on-disk 1 m DEM/hillshade untouched.
    hs_da, hs_wc = ((None, "UTM") if hs_path is None
                    else _prepare_hillshade(hs_path, "UTM", _render_px(args.dpi)))

    if not args.no_map:
        i15_tif = find(src, f"{args.from_key}_i15max.tif")
        if not os.path.exists(i15_tif):
            sys.exit(f"error: {i15_tif} not found (need the source i15max field)")
        gauges_gdf = gpd.read_file(gpath) if have_gauges else None
        out = out_path(args.out_dir, f"{key}.png")
        drape_i15(hs_da, i15_tif, out_path=out, work_crs=hs_wc,
                  cmap=args.cmap, wet_min=args.wet_min,
                  perimeters=args.perimeters, basins=args.basins,
                  highlight=args.highlight, points=args.points,
                  gauges=gauges_gdf,
                  gauge_value=("i15_mmph" if gauges_gdf is not None else None),
                  title=args.title or f"{key}  -  peak i15 (zoom{sm_note})",
                  basemap=args.basemap, basemap_provider=args.basemap_provider,
                  basemap_labels=_labels(args), basemap_zoom=args.basemap_zoom,
                  hillshade_alpha=args.hillshade_alpha,
                  alpha=zoom_alpha,
                  field_smooth=osm, field_smooth_radius_km=osr,
                  reference=args.reference, local_roads=args.local_roads,
                  label_reference=not args.no_reference_labels,
                  north_arrow=True, scale_ticks=True,
                  legend=("gauges" if gauges_gdf is not None else "all"),
                  clip=clip_gdf, clip_margin=args.clip_margin, dpi=args.dpi)
        print(f"wrote {out}")

    if not args.no_panels:
        out = out_path(args.out_dir, f"{key}_panels.png")
        diagnostic_panels(src, args.from_key, which=args.fields, out_path=out,
                          hillshade=hs_da, work_crs=hs_wc,
                          perimeters=args.perimeters,
                          gauges=(gpath if have_gauges else None),
                          reference=args.reference, local_roads=args.local_roads,
                          label_reference=not args.no_reference_labels,
                          alpha=zoom_alpha,
                          field_smooth=osm, field_smooth_radius_km=osr,
                          clip=clip_gdf, clip_margin=args.clip_margin,
                          title=args.title, dpi=args.dpi)
        print(f"wrote {out}")

    # NOAA Atlas 14 climatology comparison + anomaly for the zoom sub-AOI
    # (clim fetched fresh for the tighter extent; observed reused from the source
    # event -- MRMS has no finer resolution). Guarded: an Atlas 14 fetch failure
    # must not lose the map/panels already written above.
    if not args.no_climate:
        try:
            _run_climate(args, src=src, from_key=args.from_key,
                         out_dir=args.out_dir, key=key, bounds=bounds,
                         clip_gdf=clip_gdf, hillshade=hs_da, work_crs=hs_wc,
                         gauges_arg=(gpath if have_gauges else None))
        except Exception as ex:                # noqa: BLE001 (Atlas 14 fetch etc.)
            print(f"warning: --climate step failed ({type(ex).__name__}: {ex}); "
                  f"skipping climate maps (the zoom map/panels were still written)")


def _cmd_pick(args):
    """Write a self-contained HTML bbox picker for an event and open it in the
    browser. Drag a rectangle -> copy the ready-to-run `zoom` command. Browser
    based so it works on Windows / macOS / Linux with no GUI-toolkit dependency."""
    from .plot import bbox_picker
    src = args.from_dir
    i15 = args.i15 or find(src, f"{args.from_key}_i15max.tif")
    if not os.path.exists(i15):
        sys.exit(f"error: {i15} not found (need the source i15max field)")
    hs = args.hillshade or find(src, f"{args.from_key}_hillshade.tif")
    hs = hs if os.path.exists(hs) else None
    gauges = None
    if args.gauges:
        import geopandas as gpd
        gp = find(src, f"{args.from_key}_gauges.geojson")
        if os.path.exists(gp):
            gauges = gpd.read_file(gp)
        else:
            print(f"note: {gp} not found; skipping gauge overlay")
    out = args.out or find(src, f"{args.from_key}_pick.html")
    zoom_dir = os.path.join(src, f"{args.from_key}_zoom")
    prefix = (f'python -m stormscape zoom --from-dir "{src}" '
              f'--from-key {args.from_key}')
    suffix = (f'--out-dir "{zoom_dir}" --key {args.from_key}_zoom '
              f'--reference --crop-rasters')
    bbox_picker(i15, hillshade=hs, out_path=out, cmap=args.cmap,
                wet_min=args.wet_min, reference=not args.no_reference,
                local_roads=args.local_roads,
                label_reference=not args.no_reference_labels,
                perimeters=args.perimeters, gauges=gauges,
                cmd_prefix=prefix, cmd_suffix=suffix,
                title=f"{args.from_key} — pick a zoom bbox")
    print(f"wrote {out}")
    if not args.no_open:
        import webbrowser
        webbrowser.open(f"file://{os.path.abspath(out)}")
        print("opened in your browser — drag a rectangle, then copy the command")
    return out


def _safe_name(name):
    return "".join(ch if ch.isalnum() else "_" for ch in str(name))


def _vgauge_user_points(args):
    """Explicit user virtual-gauge points -> list of (name, lon, lat) (or [])."""
    if args.points_file:
        import geopandas as gpd
        from .mrms import _as_points
        return _as_points(gpd.read_file(args.points_file))
    if args.point:
        pts = []
        for i, s in enumerate(args.point):
            parts = s.split(",")
            name = parts[2] if len(parts) > 2 else f"VG{i + 1}"
            pts.append((name, float(parts[0]), float(parts[1])))
        return pts
    return []


def _dist_to_aoi_km(lon, lat, bounds):
    """Great-circle-ish km from a point to an AOI bbox ``(W,S,E,N)`` (0 if inside)."""
    import math
    w, s, e, n = bounds
    dx = max(w - lon, 0.0, lon - e) * 111.32 * math.cos(math.radians(lat))
    dy = max(s - lat, 0.0, lat - n) * 111.32
    return (dx * dx + dy * dy) ** 0.5


def _vgauge_aoi_bounds(args):
    """AOI bbox for the vgauge distance filter: explicit --bbox/--aoi, else the
    reused event's i15max footprint (``<from-dir>/<from-key>_i15max.tif``)."""
    if args.bbox or args.aoi:
        from .aoi import load_aoi
        return load_aoi(_aoi_from_args(args))[0]
    from_dir, from_key = getattr(args, "from_dir", None), getattr(args, "from_key", None)
    if from_dir and from_key:
        tif = find(from_dir, f"{from_key}_i15max.tif")
        if os.path.exists(tif):
            import rioxarray
            with rioxarray.open_rasterio(tif) as d:
                b = d.rio.bounds()
            return (b[0], b[1], b[2], b[3])
    return None


def _filter_stations(real, args):
    """Filter the real-gauge dict for the vgauge atlas/detail: drop stations farther
    than ``--max-dist-km`` from the AOI, (with ``--wet-only``) those that recorded no
    rain in the atlas metric, and (with ``--max-report-min``) coarse reporters whose
    native cadence exceeds that interval -- e.g. ``--max-report-min 60`` keeps only
    gauges reporting hourly or finer, dropping daily/coarse gauges whose i15/i30 are
    smeared low by the 1-min interpolation. A station whose cadence can't be measured
    (no/sparse precip obs -> NaN ``report_min``) is left to the ``--wet-only`` screen,
    not dropped on cadence. Explicit ``--point``s are unaffected (filtered elsewhere).
    Returns the kept ``{name: df}``."""
    import pandas as pd
    max_d = getattr(args, "max_dist_km", None)
    wet = getattr(args, "wet_only", False)
    max_report = getattr(args, "max_report_min", None)
    if max_d is None and not wet and max_report is None:
        return real
    bounds = _vgauge_aoi_bounds(args) if max_d is not None else None
    if max_d is not None and bounds is None:
        print("note: --max-dist-km set but no AOI found (pass --bbox/--aoi, or a "
              "--from-dir with <from-key>_i15max.tif); skipping distance filter")
    metric = getattr(args, "atlas_metric", 15)
    wet_min = getattr(args, "wet_min", 0.5)      # floor: traces below this = "dry"

    def _wet(df):
        for col in (f"i{metric}_mmph", "rate_mmph", "total_mm"):
            if col in df.columns:
                peak = float(pd.to_numeric(df[col], errors="coerce").fillna(0).max())
                return peak > wet_min
        return True                              # no rain column -> don't drop

    kept = {}
    for n, df in real.items():
        lon, lat = df.attrs.get("lon"), df.attrs.get("lat")
        if bounds is not None and lon is not None and \
                _dist_to_aoi_km(lon, lat, bounds) > max_d:
            continue
        if wet and not _wet(df):
            continue
        rm = df.attrs.get("report_min")          # native gauge tempo (min)
        if max_report is not None and pd.notna(rm) and float(rm) > max_report:
            continue
        kept[n] = df
    if len(kept) < len(real):
        crit = " + ".join(filter(None, [
            f"<={max_d:g} km from AOI" if bounds is not None else "",
            "wet-only" if wet else "",
            f"report<={max_report:g} min" if max_report is not None else ""]))
        print(f"filtered stations: {len(real)} -> {len(kept)}  ({crit})")
    return kept


def _virtual_gauge_products(args, *, points, real, reusing, user_points,
                            start, end, key, make_atlas, make_detail):
    """Shared virtual-gauge product builder used by BOTH ``gauges`` (default) and
    ``vgauge``: sample the radar source(s) (MRMS and/or NEXRAD) at the gauge
    points, write the per-gauge CSVs, and render the rainfall comparison atlas +
    the per-gauge detail figures. The atlas and the detail figures consume the
    SAME ``sources``/``real`` (the already-filtered station set), so their
    wet-gauge content is identical by construction."""
    from .plot import (plot_virtual_gauge, virtual_gauge_atlas,
                       virtual_gauge_detail)
    rgd = subdir(args.out_dir, "RainGaugeData")
    os.makedirs(rgd, exist_ok=True)
    # virtual-gauge series from one or more radar sources (MRMS and/or NEXRAD;
    # --nexrad adds NEXRAD alongside MRMS in the post-2020 both-available epoch)
    sources = {}
    if getattr(args, "source", "mrms") == "mrms":
        from .mrms import virtual_gauge_timeseries as _mrms_vg
        sources["MRMS"] = _mrms_vg(points, start, end,
                                   durations=tuple(args.durations),
                                   multisensor=not args.no_multisensor,
                                   pad_deg=args.pad_deg)
    if getattr(args, "source", "mrms") == "nexrad" or getattr(args, "nexrad", False):
        from .nexrad import virtual_gauge_timeseries as _nex_vg
        sources["NEXRAD"] = _nex_vg(points, start, end, radar=args.radar,
                                    method=args.method,
                                    durations=tuple(args.durations),
                                    cache_dir=args.cache_dir
                                    or os.path.join(args.out_dir, "nexrad_cache"))

    # case-insensitive-unique filename stems so two stations whose names differ
    # only in case (e.g. 'Virginia City' / 'VIRGINIA CITY') don't collide into one
    # file on macOS -- which would silently drop a gauge from the per-gauge CSVs /
    # detail figures while it still appears in the atlas (non-uniform output)
    from .gauges import unique_safe_names
    allnames = (set(real) | {nm for vg in sources.values() for nm in vg}
                | {n for n, _, _ in user_points})
    stems = unique_safe_names(allnames)

    # per-gauge CSVs in RainGaugeData/ (one file per gauge, Cavagna-style)
    for label, vg in sources.items():
        for name, df in vg.items():
            df.to_csv(os.path.join(
                rgd, f"{key}_vgauge_{label.lower()}_{stems[name]}.csv"))
    if not reusing:                       # reuse leaves the full-day store intact
        for name, df in real.items():
            df.to_csv(os.path.join(rgd, f"{key}_gauge_{stems[name]}.csv"))
    print("wrote CSVs [" + ", ".join(f"{l} {len(v)}" for l, v in sources.items())
          + (f", real {len(real)}" if not reusing else f", reused real {len(real)}")
          + f"] -> {rgd}")

    # individual 2-panel figures for explicit user points (MRMS primary)
    primary = sources.get("MRMS") or next(iter(sources.values()), {})
    for name, _, _ in user_points:
        if name in primary:
            png = out_path(args.out_dir, f"{key}_{stems[name]}.png")
            plot_virtual_gauge(primary[name], name=name, out_path=png,
                               durations=tuple(args.durations), dpi=args.dpi)
            print(f"wrote {png}")

    # atlas: every radar source + real overlay where matched
    if make_atlas:
        atlas = out_path(args.out_dir, f"{key}_vg_atlas.png")
        virtual_gauge_atlas(sources, real_series=real,
                            value=f"i{args.atlas_metric}_mmph",
                            out_path=atlas, dpi=args.dpi)
        print(f"wrote {atlas}")

    # detail: one big 4-row figure per gauge (cumulative + I60/I30/I15), same
    # line styles as the atlas, into VirtualGaugeFigures/ -- SAME gauge set as the
    # atlas (uniform wet gauges across both products)
    if make_detail:
        fdir = subdir(args.out_dir, "VirtualGaugeFigures")
        os.makedirs(fdir, exist_ok=True)
        names = []
        for s in (*sources.values(), real):
            for nm in s:
                if nm not in names:
                    names.append(nm)
        for nm in names:
            fp = os.path.join(fdir, f"{key}_vgdetail_{stems[nm]}.png")
            rm = (real[nm].attrs.get("report_min")     # native gauge tempo -> title
                  if real.get(nm) is not None else None)
            virtual_gauge_detail(sources, nm, real_series=real, report_min=rm,
                                 out_path=fp, dpi=args.dpi)
        print(f"wrote {len(names)} detail figures -> {fdir}")
    return sources


def _cmd_vgauge(args):
    os.makedirs(args.out_dir, exist_ok=True)
    # analysis window (lazy): explicit --date / --start+--end now; else, for a
    # reused store, derived from the store's storm rain window below
    start = end = None
    if getattr(args, "date", None) or (args.start and args.end):
        start, end = _gauge_window(args)
    key = args.key or "vgauge"

    user_points = _vgauge_user_points(args)
    points = list(user_points)

    # --gauges: drop a VG at every real Synoptic station + pull the real series.
    # Reuse the canonical store (the full storm-DAY record written by `gauges`) if
    # --from-dir is given, trimming to the storm's rain window; else fetch live.
    real, reusing = {}, False
    if args.gauges:
        store_gj = (find(args.from_dir, f"{args.from_key}_gauges.geojson")
                    if args.from_dir and args.from_key else None)
        if store_gj and os.path.exists(store_gj) and not args.refetch:
            import geopandas as gpd
            from .gauges import load_event_series
            gj = gpd.read_file(store_gj)
            real = load_event_series(find_subdir(args.from_dir, "RainGaugeData"),
                                     args.from_key, gj)
            # Narrow the storm DAY to the storm itself. `start` is already set
            # when --date was given, so gating this on `start is None` (as it
            # once was) silently skipped the trim for every --date run and handed
            # the full ~30 h to the NEXRAD fetch. Only an explicit --start/--end
            # should survive untouched.
            if start is None or not _explicit_window(args):
                if start is None and real:
                    start = min(d.index.min() for d in real.values())
                    end = max(d.index.max() for d in real.values())
                start, end = _narrow_to_storm(args, start, end, series=real)
            clipped = {}
            for n, df in real.items():              # trim to the rain window
                d = df.loc[start:end]
                d.attrs = dict(df.attrs)
                clipped[n] = d
            real, reusing = clipped, True
            print(f"reusing gauge store: {len(real)} series from {args.from_dir}/"
                  f"RainGaugeData, rain window {start:%m-%d %H:%M}-{end:%m-%d %H:%M}")
        else:
            if start is None:                       # a live fetch needs a window
                start, end = _gauge_window(args)
            from .gauges import gauge_timeseries
            real, _ = gauge_timeseries(_aoi_from_args(args), start, end,
                                       token=args.token,
                                       durations=tuple(args.durations),
                                       pad_deg=args.pad_deg)
        real = _filter_stations(real, args)         # --max-dist-km / --wet-only
        points += [(n, df.attrs["lon"], df.attrs["lat"]) for n, df in real.items()]
    if not points:
        sys.exit("error: provide --point / --points-file and/or --gauges")
    if start is None:                               # explicit points w/o a window
        start, end = _gauge_window(args)
        # explicit --point runs have no gauge series to derive a window from, so
        # probe the radar before a NEXRAD fetch rather than pulling the whole day
        if getattr(args, "nexrad", False) or args.source == "nexrad":
            start, end = _narrow_to_storm(args, start, end)

    # atlas with --atlas or the all-stations --gauges mode; detail with --detail
    _virtual_gauge_products(
        args, points=points, real=real, reusing=reusing, user_points=user_points,
        start=start, end=end, key=key,
        make_atlas=(args.atlas or args.gauges), make_detail=args.detail)


def _dbz_label(field):
    return {"reflectivity": "reflectivity (dBZ)",
            "velocity": "radial velocity (m s$^{-1}$)",
            "cross_correlation_ratio": "correlation coefficient"}.get(field, field)


def _nexrad_when(args):
    """UTC datetime of the scan to grid: --time on --date, else window midpoint."""
    import datetime as dt

    from .mrms import parse_date
    if getattr(args, "time", None):
        if not args.date:
            sys.exit("error: --time needs --date")
        d = parse_date(args.date)
        return dt.datetime(d.year, d.month, d.day,
                           int(args.time[:2]), int(args.time[2:]))
    start, end = _gauge_window(args)
    return start + (end - start) / 2


def _cmd_nexrad(args):
    from .mrms import save_fields
    from .nexrad import (intensity_stack, nearest_radar, reflectivity_composite,
                         reflectivity_field)
    from .plot import drape_i15
    os.makedirs(args.out_dir, exist_ok=True)
    aoi = _aoi_from_args(args)
    key = args.key or "nexrad"
    _save_event_aoi(args, args.out_dir, key)        # for climate auto-extent-match
    cache = args.cache_dir or os.path.join(args.out_dir, "nexrad_cache")
    radar = (args.radar.upper() if args.radar else nearest_radar(aoi)[0])

    if args.intensity:                       # i15/i30/i60 stack (Z-R rain rate)
        start, end = _gauge_window(args)
        res = intensity_stack(aoi, start, end, radar=radar, a=args.zr_a,
                              b=args.zr_b,
                              dbz_cap=(None if args.no_hail_cap else args.dbz_cap),
                              sweep=args.sweep, res_m=args.res_m, cache_dir=cache,
                              method=args.method, z_blend=args.z_blend,
                              rate_cap=args.rate_cap,
                              blockage_dem=args.blockage_dem, cbb_max=args.cbb_max)
        paths = save_fields(res, args.out_dir, key)
        field_tif = out_path(args.out_dir, f"{key}_i15max.tif")
        cbar = "peak 15-min intensity  i15  (mm h$^{-1}$)"
        deftitle = f"{radar}  L2 i15 [{args.method}]  ({start:%Y-%m-%d})"
    elif args.composite:                     # storm-peak reflectivity field
        start, end = _gauge_window(args)
        res = reflectivity_composite(aoi, start, end, radar=radar,
                                     field=args.field, sweep=(args.sweep or 0),
                                     res_m=args.res_m, cache_dir=cache,
                                     max_scans=args.max_scans)
        paths = save_fields(res, args.out_dir, key)
        field_tif = paths[0]
        cbar = _dbz_label(args.field)
        deftitle = f"{radar}  {args.field}  (storm-peak)"
    else:                                    # single scan nearest a time
        res = reflectivity_field(aoi, _nexrad_when(args), radar=radar,
                                 field=args.field, sweep=(args.sweep or 0),
                                 res_m=args.res_m, cache_dir=cache)
        paths = save_fields(res, args.out_dir, key)
        field_tif = paths[0]
        cbar = _dbz_label(args.field)
        deftitle = f"{radar}  {args.field}  ({res['meta'].get('scan_time', '')})"
    print(json.dumps(res["meta"], indent=2))
    for p in paths:
        print(f"wrote {p}")

    # optional rain-gauge overlay, coloured by the radar value sampled beneath
    # each gauge (same scale as the field) -- mirrors the run/compare maps.
    gauges_gdf, gv = None, None
    if args.gauges_file:
        import geopandas as gpd
        gauges_gdf = gpd.read_file(args.gauges_file)
    elif args.gauges:
        gauges_gdf, _, _ = _fetch_gauges(args)
    if gauges_gdf is not None and len(gauges_gdf):
        from .compare import sample_raster_at_points
        gauges_gdf = gauges_gdf.copy()
        gauges_gdf["radar_val"] = sample_raster_at_points(gauges_gdf, field_tif)
        gv = "radar_val"

    out_png = out_path(args.out_dir, f"{key}_nexrad.png")
    drape_i15(args.hillshade, field_tif, out_path=out_png, work_crs="UTM",
              cmap=args.cmap, wet_min=args.wet_min, cbar_label=cbar,
              perimeters=args.perimeters, reference=args.reference,
              local_roads=args.local_roads,
              label_reference=not args.no_reference_labels,
              gauges=gauges_gdf, gauge_value=gv,
              alpha=args.alpha,
              north_arrow=True, scale_ticks=True,
              legend=("gauges" if gauges_gdf is not None else "all"),
              clip=((args.perimeters or args.aoi) if args.clip else None),
              clip_margin=args.clip_margin, dpi=args.dpi,
              title=args.title or deftitle)
    print(f"wrote {out_png}")
    return field_tif


def _render_px(dpi):
    """Max hillshade pixels (long side) to render at ``dpi``. A couple thousand
    covers the widest panel (~9-12 in) with margin, so downsampling the terrain
    to this is invisible in the output while avoiding the slow/memory-heavy
    full-resolution render of fine (e.g. 1 m) hillshades."""
    return max(2500, int(round((dpi or 200) * 12)))


def _run_climate(args, *, src, from_key, out_dir, key, bounds, clip_gdf,
                 hillshade, work_crs="UTM", gauges_arg, title=None):
    """Fetch the NOAA Atlas 14 climatology for ``bounds`` and write the
    observed-vs-climatology comparison figure + per-duration anomaly maps for an
    already-resolved AOI / clip / terrain.

    Shared by the ``climate`` command and ``zoom --climate``. The observed fields
    are read from ``<src>/<from_key>_i*max.tif``; the climatology and the figures
    are written under ``out_dir`` with ``key``. The observed field is smoothed
    (default Gaussian 1 km) before both the comparison and the anomaly."""
    from . import atlas14
    from .mrms import save_fields
    from .plot import anomaly_map, climatology_comparison
    durations = tuple(args.durations)
    # observed-field smoothing (default Gaussian 1 km) so the peaky ~1 km radar
    # field reads against the smooth ~800 m climatology
    osm = None if args.obs_smooth == "none" else args.obs_smooth
    osr = args.obs_smooth_radius
    if osm and osr > 0:
        print(f"smoothing observed fields for climate plots: {osm}, {osr:g} km")

    # 1) climatology fields
    cache = os.path.join(out_dir, "atlas14_cache")
    clim = atlas14.climatology_field(bounds, durations=durations, ari=args.ari,
                                     region=args.region, stat=args.stat,
                                     cache_dir=cache, pad_deg=args.pad_deg)
    save_fields(clim, out_dir, f"{key}_clim")
    print(f"wrote {out_dir}/{key}_clim_i*.tif  (NOAA Atlas 14 region "
          f"'{clim['meta']['region']}', {args.ari}-yr {args.stat})")

    # 2) observed-vs-climatology comparison figure
    if not args.no_comparison:
        out = out_path(out_dir, f"{key}_climate_compare.png")
        climatology_comparison(
            out_dir, src, key, durations=durations, obs_key=from_key,
            hillshade=hillshade, work_crs=work_crs,
            perimeters=args.perimeters, gauges=gauges_arg,
            reference=args.reference, local_roads=args.local_roads,
            label_reference=not args.no_reference_labels, alpha=args.alpha,
            cmap=args.cmap, shared_row_scale=args.shared_row_scale, ari=args.ari,
            obs_smooth=osm, obs_smooth_radius_km=osr,
            clip=clip_gdf, clip_margin=args.clip_margin,
            title=title or f"{key}  -  observed vs NOAA Atlas 14 {args.ari}-yr",
            dpi=args.dpi, out_path=out)
        print(f"wrote {out}")

    # 3) per-duration anomaly maps (observed / climatology)
    if not args.no_anomaly:
        for d in durations:
            clim_tif = out_path(out_dir, f"{key}_clim_i{d}.tif")
            obs_tif = find(src, f"{from_key}_i{d}max.tif")
            if not (os.path.exists(clim_tif) and os.path.exists(obs_tif)):
                continue
            obs_in = obs_tif
            if osm and osr > 0:                      # smooth obs before dividing
                from .smoothing import smooth_dataarray
                obs_in = smooth_dataarray(obs_tif, osm, osr)
            ratio = atlas14.anomaly(obs_in, clim_tif, eps=args.eps)
            atif = out_path(out_dir, f"{key}_anom_i{d}.tif")
            ratio.rio.to_raster(atif)
            out = out_path(out_dir, f"{key}_anom_i{d}.png")
            anomaly_map(ratio, hillshade=hillshade, work_crs=work_crs,
                        duration=d, ari=args.ari,
                        cmap=args.anomaly_cmap, perimeters=args.perimeters,
                        gauges=gauges_arg, reference=args.reference,
                        local_roads=args.local_roads,
                        label_reference=not args.no_reference_labels,
                        alpha=args.alpha, clip=clip_gdf,
                        clip_margin=args.clip_margin, dpi=args.dpi, out_path=out)
            print(f"wrote {atif}\nwrote {out}")


def _cmd_climate(args):
    """NOAA Atlas 14 rainfall climatology vs an event's observed i15/i30/i60.

    Re-uses an already-processed event's observed fields (``<from-key>_i*max.tif``
    from MRMS or NEXRAD) -- no radar re-run. Fetches the matching Atlas 14
    gridded climatology, writes ``<key>_clim_i*.tif``, a 3xN observed-vs-clim
    comparison figure, and per-duration anomaly (observed / climatology) maps."""
    import geopandas as gpd

    from .aoi import bbox_polygon, load_aoi
    from .plot import _prepare_hillshade
    os.makedirs(args.out_dir, exist_ok=True)
    src, from_key = args.from_dir, args.from_key
    key = args.key or from_key

    obs_i15 = find(src, f"{from_key}_i15max.tif")
    if not os.path.exists(obs_i15):
        sys.exit(f"error: {obs_i15} not found (need the observed i15max field)")

    # AOI extent: explicit --aoi/--bbox; else auto-match the event's saved/placed
    # AOI so the figures frame exactly like the i15 maps (which clip to the AOI, not
    # the wider i15max footprint = AOI + MRMS pad); else fall back to that footprint.
    if args.bbox or args.aoi:
        bounds, geom = load_aoi(_aoi_from_args(args))
        if not args.perimeters and args.aoi:        # draw the AOI outline by default
            args.perimeters = args.aoi
    else:
        auto_aoi = _find_event_aoi(src, from_key)
        if auto_aoi:
            bounds, geom = load_aoi(auto_aoi)
            if not args.perimeters:
                args.perimeters = auto_aoi
            print(f"auto-matched figure extent to event AOI: "
                  f"{os.path.basename(auto_aoi)}")
        else:
            import rioxarray
            with rioxarray.open_rasterio(obs_i15) as d:
                b = d.rio.bounds()
            bounds, geom = (b[0], b[1], b[2], b[3]), None
            print("note: no event AOI found; framing to the i15max footprint "
                  "(pass --aoi/--bbox to match the i15 maps)")
    clip_gdf = gpd.GeoDataFrame(geometry=[geom or bbox_polygon(bounds)], crs=4326)

    hs_path = args.hillshade or find(src, f"{from_key}_hillshade.tif")
    if not os.path.exists(hs_path):
        print(f"note: hillshade {hs_path} not found; rendering without terrain")
        hs_path = None
    # prep the hillshade once (reproject + downsample to render resolution) so the
    # 3xN comparison + N anomaly figures don't each re-process a big raster
    hs_da, hs_wc = ((None, "UTM") if hs_path is None
                    else _prepare_hillshade(hs_path, "UTM", _render_px(args.dpi)))
    gpath = find(src, f"{from_key}_gauges.geojson")
    have_gauges = args.gauges and os.path.exists(gpath)
    if args.gauges and not have_gauges:
        print(f"note: {gpath} not found; skipping gauge overlay")
    gauges_arg = gpath if have_gauges else None

    _run_climate(args, src=src, from_key=from_key, out_dir=args.out_dir, key=key,
                 bounds=bounds, clip_gdf=clip_gdf, hillshade=hs_da, work_crs=hs_wc,
                 gauges_arg=gauges_arg, title=args.title)


def _burn_class_table(result, bounds, path):
    """Per-class pixel count and **true** ground area of a severity mosaic.

    The grid is EPSG:3857, whose metres shrink with latitude, so a nominal 60 m
    cell is ~46.5 m of ground at 39 deg N. Areas computed on the raw transform
    would be ~66% too large there; scale by cos(lat) first.
    """
    import math

    import numpy as np
    import pandas as pd

    from .burn import SBS_LABELS, SEVERITY_SCHEMES

    cls = result["fields"].get("severity")
    if cls is None:
        return None
    lat = (bounds[1] + bounds[3]) / 2.0
    res = abs(result["transform"].a) * math.cos(math.radians(lat))
    cell_km2 = (res / 1000.0) ** 2
    scheme = result["meta"].get("scheme")
    if scheme == "baer":
        labels = SBS_LABELS
    else:
        labels = dict(enumerate(SEVERITY_SCHEMES[scheme]["labels"]))
    finite = cls[np.isfinite(cls)]
    rows = []
    for c in sorted(set(np.unique(finite).tolist())):
        n = int((finite == c).sum())
        rows.append(dict(cls=int(c), label=labels.get(int(c), f"class {int(c)}"),
                         pixels=n, area_km2=round(n * cell_km2, 4),
                         fraction=round(n / max(finite.size, 1), 4)))
    df = pd.DataFrame(rows)
    df.to_csv(path, index=False)
    return df


def _burn_display_defaults(product):
    """``(wet_min, vmax)`` for the burn drape, which differ by product.

    dNBR is a continuous 0-1 index cut at the 0.10 unburned break; BAER soil
    burn severity is classes 1-4, so it needs a class-height scale and its cut
    just above class 1 ("unburned/very low"). Sharing one default would either
    clip the class map to its lowest class or paint unburned ground."""
    if product == "sbs":
        return 1.5, 4.0
    return 0.10, 1.0


def _cmd_burn(args):
    """Near-real-time burn severity (CIMSS BRISK dNBR) over an AOI.

    Screens the archive for fires intersecting the AOI, caches their GeoTIFFs,
    mosaics them onto one grid, and drapes the result over a hillshade."""
    import geopandas as gpd

    from .aoi import bbox_polygon, load_aoi
    from .burn import (CACHE_DIR, MATURITY_DAYS, SBS_LABELS, burn_severity,
                       find_scenes, register_baer_cmap, severity_colors)
    from .mrms import save_fields
    from .plot import drape_i15

    aoi = _aoi_from_args(args)
    bounds, geom = load_aoi(aoi)
    key = args.key or "burn"
    cache = args.cache_dir or os.path.join(args.out_dir, CACHE_DIR)
    years = args.years or None

    if args.list:                       # cheap: headers only, nothing downloaded
        sc = find_scenes(aoi, product=args.product, date=args.date,
                         since=args.since, latest=not args.all_dates,
                         years=years, min_age_days=args.min_age,
                         cache_dir=cache)
        if not len(sc):
            print("no burn-severity scenes intersect this AOI")
            return None
        print(f"\n{'fire':<28} {'st':<3} {'date':<12} {'age':>5}  overlap")
        for r in sc.itertuples(index=False):
            flag = " *" if r.age_days < MATURITY_DAYS else ""
            print(f"{str(r.fire):<28} {str(r.state or ''):<3} "
                  f"{str(r.date):<12} {r.age_days:>4}d  {r.overlap:.5f}{flag}")
        if (sc.age_days < MATURITY_DAYS).any():
            print(f"  * composite younger than {MATURITY_DAYS} d -- pattern is "
                  f"reliable, magnitude may under-read (see --min-age)")
        return sc

    os.makedirs(args.out_dir, exist_ok=True)
    res = burn_severity(aoi, date=args.date, product=args.product,
                        scheme=args.scheme, since=args.since,
                        fires=args.fire, years=years, cache_dir=cache,
                        min_age_days=args.min_age, pad_deg=args.pad_deg,
                        workers=args.workers)
    if res is None:
        print("nothing to map: no fire in this AOI "
              "(widen --bbox/--aoi, or check --date / --since)")
        return None

    m = res["meta"]
    dates = sorted(set(m["scene_dates"]))
    print(f"\n{len(m['fires'])} fire(s): {', '.join(map(str, m['fires']))}")
    print(f"scene dates: {dates[0]}" + (f" .. {dates[-1]}" if len(dates) > 1 else ""))
    if len(dates) > 1:
        # a mixed-date mosaic is legitimate (fires burn at different times) but
        # it is not one snapshot, so say so rather than let it pass unnoticed
        print("note: scenes span several dates; each fire is shown at its own "
              "latest observation on or before --date")

    paths = save_fields(res, args.out_dir, key, layout=args.layout)
    for p in paths:
        print(f"wrote {p}")

    scenes = find_scenes(aoi, product=args.product, date=args.date,
                         since=args.since, latest=not args.all_dates,
                         years=years, min_age_days=args.min_age,
                         cache_dir=cache, verbose=False)
    if len(scenes):
        gp = out_path(args.out_dir, f"{key}_burn_scenes.geojson",
                      layout=args.layout)
        scenes.assign(date=scenes.date.astype(str)).to_file(gp, driver="GeoJSON")
        print(f"wrote {gp}")

    tpath = out_path(args.out_dir, f"{key}_burn_classes.csv", layout=args.layout)
    tbl = _burn_class_table(res, bounds, tpath)
    if tbl is not None:
        print(f"wrote {tpath}")
        print(tbl.to_string(index=False))

    if args.no_map:
        return res

    hs_path = args.hillshade or find(args.out_dir, f"{key}_hillshade.tif")
    if args.dem and not os.path.exists(hs_path):
        from .dem import fetch_dem_and_hillshade
        dem_path = out_path(args.out_dir, f"{key}_dem.tif", layout=args.layout)
        hs_path = out_path(args.out_dir, f"{key}_hillshade.tif",
                           layout=args.layout)
        print(f"fetching {args.resolution} m DEM for the AOI")
        fetch_dem_and_hillshade(aoi, resolution=args.resolution,
                                dst_crs=args.dst_crs, dem_path=dem_path,
                                hillshade_path=hs_path, pad_deg=args.pad_deg,
                                resampling=getattr(args, "resampling", None))
    if not os.path.exists(hs_path):
        print(f"note: no hillshade at {hs_path}; rendering without terrain "
              f"(pass --hillshade, or --dem to fetch one)")
        hs_path = None

    d_wet, d_vmax = _burn_display_defaults(args.product)
    wet_min = d_wet if args.wet_min is None else args.wet_min
    vmax = d_vmax if args.vmax is None else args.vmax
    field = "severity" if args.product == "sbs" else "dnbr"

    # Classed by default, in the BAER class colours -- how burn severity is
    # actually published. --continuous falls back to the smooth ramp (same
    # palette), and an explicit --cmap opts out of the class colours entirely.
    norm = cbar_ticks = cbar_ticklabels = None
    cmap = args.cmap
    if cmap in (None, "baer", "brisk"):
        if args.continuous or args.product == "sbs":
            cmap = register_baer_cmap().name
        else:
            cmap, norm, cbar_ticks, cbar_ticklabels = severity_colors(args.scheme)
            vmax = None                   # the norm owns the scale now
    if args.product == "sbs":             # already classed 1-4; label them
        cbar_ticks = sorted(SBS_LABELS)
        cbar_ticklabels = [SBS_LABELS[k] for k in cbar_ticks]
    fpath = find(args.out_dir, f"{key}_{field}.tif")
    png = out_path(args.out_dir, f"{key}_burn.png", layout=args.layout)
    label = ("soil burn severity (BAER)" if args.product == "sbs"
             else "burn severity  (dNBR class)" if norm is not None
             else "dNBR [-]")
    title = args.title or (
        f"Burn severity ({'BAER SBS' if args.product == 'sbs' else 'BRISK dNBR'})"
        f" - {', '.join(map(str, m['fires']))[:60]}  {dates[-1]}")
    clip_gdf = (gpd.GeoDataFrame(geometry=[geom or bbox_polygon(bounds)], crs=4326)
                if args.clip else None)
    drape_i15(hs_path, fpath, out_path=png, work_crs="UTM",
              wet_min=wet_min, vmax=vmax, cmap=cmap, norm=norm,
              cbar_ticks=cbar_ticks, cbar_ticklabels=cbar_ticklabels,
              alpha=args.alpha, title=title, cbar_label=label,
              perimeters=args.perimeters, basins=args.basins,
              highlight=args.highlight, points=args.points,
              reference=args.reference, local_roads=args.local_roads,
              label_reference=not args.no_reference_labels,
              basemap=args.basemap, basemap_provider=args.basemap_provider,
              basemap_labels=args.basemap_labels, basemap_zoom=args.basemap_zoom,
              hillshade_alpha=args.hillshade_alpha,
              clip=clip_gdf, clip_margin=args.clip_margin,
              north_arrow=True, scale_ticks=True, dpi=args.dpi)
    print(f"wrote {png}")
    _save_event_aoi(args, args.out_dir, key)
    return res


def _cmd_export(args):
    """Georeferenced exports for GIS / CalTopo from an already-processed event:
    EPSG:3857 GeoTIFFs of the rainfall fields (raw float + colorized RGBA) and
    GeoPDFs of the primary map figures (i15 map + anomaly map). Re-uses the
    event's rasters -- no radar re-run."""
    import geopandas as gpd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from . import export as _export
    from .aoi import bbox_polygon, load_aoi
    from .plot import _prepare_hillshade, anomaly_map, drape_i15
    os.makedirs(args.out_dir, exist_ok=True)
    src, from_key = args.from_dir, args.from_key
    key = args.key or from_key

    # 1) EPSG:3857 GeoTIFF layers (raw float + colorized RGBA) ----------------
    if not args.no_geotiff:
        for p in _export.export_geotiffs(
                src, from_key, args.out_dir, fields=args.layers, dst_crs=args.crs,
                out_key=key, raw=not args.no_raw, colorize=not args.no_rgb,
                cmap=args.cmap, anomaly_cmap=args.anomaly_cmap,
                wet_min=args.wet_min):
            print(f"wrote {p}")

    # 2) full-resolution NHD stream network as a vector layer ----------------
    if args.streams:
        aoi_spec = (_aoi_from_args(args) if (args.bbox or args.aoi)
                    else _find_event_aoi(src, from_key))
        if not aoi_spec:
            print("note: no AOI for --streams (pass --aoi/--bbox, or place "
                  f"{from_key}_aoi.geojson in --from-dir); skipping streams")
        else:
            out_s = out_path(args.out_dir,
                             f"{key}_streams.{args.streams_format}")
            if _export.export_streams(aoi_spec, out_s,
                                      named_only=args.streams_named_only,
                                      clip_to_aoi=not args.streams_bbox):
                print(f"wrote {out_s}")

    # 3) georeferenced PDFs of the primary figures ---------------------------
    if args.no_figures:
        return
    if not _export.geopdf_supported():
        print("note: GDAL has no PDF driver; skipping the GeoPDF figures "
              "(install: conda install -c conda-forge libgdal-pdf)")
        return

    # frame exactly like the i15 maps: explicit --aoi/--bbox, else auto-match the
    # event AOI (the maps clip to the AOI, not the wider padded radar footprint)
    if args.bbox or args.aoi:
        bounds, geom = load_aoi(_aoi_from_args(args))
        if not args.perimeters and args.aoi:
            args.perimeters = args.aoi
    else:
        auto = _find_event_aoi(src, from_key)
        if auto:
            bounds, geom = load_aoi(auto)
            if not args.perimeters:
                args.perimeters = auto
            print(f"auto-matched figure extent to event AOI: "
                  f"{os.path.basename(auto)}")
        else:
            bounds, geom = None, None
            print("note: no event AOI found; framing to the field footprint "
                  "(pass --aoi/--bbox to match the i15 maps)")
    clip_gdf = (gpd.GeoDataFrame(geometry=[geom or bbox_polygon(bounds)], crs=4326)
                if bounds is not None else None)

    # terrain backdrop + the figure CRS. Default --pdf-crs UTM renders the GeoPDF
    # identically to the PNG deliverables; the georeferencing is exact in whatever
    # projected CRS we render in (pass EPSG:3857 to match the GeoTIFF layers).
    hs_path = args.hillshade or find(src, f"{from_key}_hillshade.tif")
    if not os.path.exists(hs_path):
        print(f"note: hillshade {hs_path} not found; rendering without terrain")
        hs_path = None
    hs_da, hs_wc = ((None, args.pdf_crs) if hs_path is None
                    else _prepare_hillshade(hs_path, args.pdf_crs,
                                            _render_px(args.dpi)))

    gpath = find(src, f"{from_key}_gauges.geojson")
    have_gauges = args.gauges and os.path.exists(gpath)
    if args.gauges and not have_gauges:
        print(f"note: {gpath} not found; skipping gauge overlay")
    gauges_gdf = gpd.read_file(gpath) if have_gauges else None

    # i15-map GeoPDF
    if not args.no_i15:
        i15_tif = find(src, f"{from_key}_i15max.tif")
        if not os.path.exists(i15_tif):
            print(f"note: {i15_tif} not found; skipping the i15 GeoPDF")
        else:
            fig, ax = drape_i15(
                hs_da, i15_tif, out_path=None, work_crs=hs_wc, cmap=args.cmap,
                wet_min=args.wet_min, perimeters=args.perimeters,
                basins=args.basins, highlight=args.highlight, points=args.points,
                gauges=gauges_gdf,
                gauge_value=("i15_mmph" if gauges_gdf is not None else None),
                title=args.title or f"{key}  -  peak i15 (georeferenced)",
                hillshade_alpha=args.hillshade_alpha, alpha=args.alpha,
                reference=args.reference, local_roads=args.local_roads,
                label_reference=not args.no_reference_labels,
                north_arrow=True, scale_ticks=True,
                legend=("gauges" if gauges_gdf is not None else "all"),
                clip=clip_gdf, clip_margin=args.clip_margin, dpi=args.dpi)
            out_pdf = out_path(args.out_dir, f"{key}.pdf")
            _export.figure_to_geopdf(fig, ax, out_pdf, crs=hs_wc, dpi=args.dpi,
                                     title=f"{key} peak i15",
                                     subject="peak 15-min rainfall intensity (mm/h)")
            plt.close(fig)
            print(f"wrote {out_pdf}")

    # anomaly-map GeoPDF (needs <from-key>_anom_i<d>.tif from `climate`)
    if not args.no_anom:
        d = args.anom_duration
        anom_tif = find(src, f"{from_key}_anom_i{d}.tif")
        if not os.path.exists(anom_tif):
            print(f"note: {anom_tif} not found; skipping the anomaly GeoPDF "
                  f"(run `climate` first to write it)")
        else:
            fig, ax = anomaly_map(
                anom_tif, hillshade=hs_da, out_path=None, work_crs=hs_wc,
                duration=d, ari=args.ari, cmap=args.anomaly_cmap,
                perimeters=args.perimeters, gauges=gauges_gdf,
                reference=args.reference, local_roads=args.local_roads,
                label_reference=not args.no_reference_labels, alpha=args.alpha,
                clip=clip_gdf, clip_margin=args.clip_margin, dpi=args.dpi)
            out_pdf = out_path(args.out_dir, f"{key}_anom_i{d}.pdf")
            _export.figure_to_geopdf(fig, ax, out_pdf, crs=hs_wc, dpi=args.dpi,
                                     title=f"{key} I{d} anomaly",
                                     subject=f"observed I{d} / {args.ari}-yr "
                                             f"NOAA Atlas 14")
            plt.close(fig)
            print(f"wrote {out_pdf}")


def _analysis_gauges(args, src, from_key):
    """Resolve the gauge GeoDataFrame for ``smooth --gauge-analysis``:
    ``--gauges-file`` overrides; else the source ``<from-key>_gauges.geojson``;
    else a live Synoptic fetch (needs --aoi/--bbox + --date/--start/--end)."""
    import geopandas as gpd
    if getattr(args, "gauges_file", None):
        return gpd.read_file(args.gauges_file)
    gpath = find(src, f"{from_key}_gauges.geojson")
    if os.path.exists(gpath):
        return gpd.read_file(gpath)
    g, _, _ = _fetch_gauges(args)
    return g


def _gauge_peak_times(rgd_dir, key, names, duration=15):
    """``{gauge_name: timestamp}`` of each gauge's I{duration} peak, read from the
    saved RainGaugeData per-gauge series CSVs (offline -- no Synoptic token).
    Missing/unreadable files are skipped."""
    import pandas as pd
    out = {}
    if not rgd_dir or not os.path.isdir(rgd_dir):
        return out
    col = f"i{duration}_mmph"
    for n in names:
        p = os.path.join(rgd_dir, f"{key}_gauge_{_safe_name(n)}.csv")
        if not os.path.exists(p):
            continue
        try:
            s = pd.read_csv(p, index_col=0, parse_dates=True)[col]
            if s.notna().any():
                out[n] = s.idxmax()
        except Exception:                                # noqa: BLE001
            pass
    return out


def _ri_str(ri):
    """Format a recurrence interval (years) for display: numeric incl. sub-1-yr
    (``>=1`` to 1 dp, below to 2 sig figs), ``>1000`` above the top tabulated ARI,
    blank if missing."""
    import math
    if ri is None or (isinstance(ri, float) and math.isnan(ri)):
        return ""
    if math.isinf(ri):
        return ">1000"
    return f"{ri:.1f}" if ri >= 1 else f"{ri:.2g}"


def _write_recurrence_md(df, durations, path, with_time):
    import pandas as pd
    def f1(x): return "" if (x is None or pd.isna(x)) else f"{x:.1f}"
    with_dist = "dist_to_aoi_km" in df.columns
    with_tempo = "report_min" in df.columns
    cols = ["Gauge"]
    if with_tempo:
        cols.append("tempo min")
    if with_time:
        cols.append("I15 peak (UTC)")
    if with_dist:
        cols.append("dist km")
    cols += ([f"I{d}" for d in durations] + [f"anom{d}" for d in durations]
             + [f"RI{d}" for d in durations])
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        cells = [str(r.get("gauge"))[:28]]
        if with_tempo:
            rm = r.get("report_min")
            cells.append("" if (rm is None or pd.isna(rm)) else f"{float(rm):g}")
        if with_time:
            t = r.get("i15_peak_time")
            cells.append("" if pd.isna(t) else pd.Timestamp(t).strftime("%m-%d %H:%M"))
        if with_dist:
            d0 = r.get("dist_to_aoi_km")
            cells.append("in" if (pd.notna(d0) and d0 <= 0) else f1(d0))
        cells += [f1(r.get(f"i{d}")) for d in durations]
        cells += [f1(r.get(f"anom_i{d}")) for d in durations]
        cells += [_ri_str(r.get(f"RI_i{d}")) for d in durations]
        lines.append("| " + " | ".join(cells) + " |")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def _cmd_recurrence(args):
    """Wet-gauge recurrence table: peak I15/I30/I60, the time of the I15 peak,
    the anomaly vs the 1-yr NOAA Atlas 14 climatology, and the recurrence interval
    of each observed peak (climatology from the NOAA PFDS point service)."""
    import geopandas as gpd
    import pandas as pd
    from .compare import gauge_recurrence_table
    os.makedirs(args.out_dir, exist_ok=True)
    src, from_key = args.from_dir, args.from_key
    key = args.key or from_key
    durations = tuple(args.durations)
    gpath = args.gauges_file or find(src, f"{from_key}_gauges.geojson")
    if not os.path.exists(gpath):
        sys.exit(f"error: {gpath} not found (need the event gauges geojson)")
    gauges = gpd.read_file(gpath)

    # I15 time-of-peak: the canonical store now carries it in the geojson; only
    # fall back to scanning the series CSVs if the geojson lacks the column.
    peak_times = None
    if not args.no_peak_time and "i15_peak_time" not in gauges.columns:
        rgd = args.raingauge_dir or find_subdir(src, "RainGaugeData")
        names = list(gauges["name"]) if "name" in gauges.columns else []
        peak_times = _gauge_peak_times(rgd, from_key, names, duration=15)
        if not peak_times:
            print(f"note: no I15 time-of-peak in the geojson and no per-gauge "
                  f"series CSVs in {rgd}; time-of-peak omitted")

    # AOI for the in/near-AOI flag: --bbox/--aoi if given, else the event footprint
    aoi_bounds = None
    if args.bbox or args.aoi:
        from .aoi import load_aoi
        aoi_bounds, _ = load_aoi(_aoi_from_args(args))
    else:
        i15tif = find(src, f"{from_key}_i15max.tif")
        if os.path.exists(i15tif):
            import rioxarray
            with rioxarray.open_rasterio(i15tif) as d:
                b = d.rio.bounds()
            aoi_bounds = (b[0], b[1], b[2], b[3])

    print(f"querying NOAA PFDS ({args.series}, {args.stat}) per wet gauge...")
    df = gauge_recurrence_table(gauges, durations=durations, stat=args.stat,
                                series=args.series, peak_times=peak_times,
                                aoi_bounds=aoi_bounds)
    if not len(df):
        sys.exit("no wet gauges (i15 > 0) in the gauges file")
    if args.max_dist_km is not None and "dist_to_aoi_km" in df.columns:
        before = len(df)
        df = df[df["dist_to_aoi_km"] <= args.max_dist_km].reset_index(drop=True)
        if len(df) < before:
            print(f"  dropped {before - len(df)} gauge(s) > {args.max_dist_km} km "
                  f"from the AOI (regional outliers)")

    with_time = "i15_peak_time" in df.columns
    csv = out_path(args.out_dir, f"{key}_gauge_recurrence.csv")
    df.to_csv(csv, index=False)
    md = out_path(args.out_dir, f"{key}_gauge_recurrence.md")
    _write_recurrence_md(df, durations, md, with_time)
    print(f"wrote {csv}\nwrote {md}")
    rare = df[pd.to_numeric(df.get("RI_i15"), errors="coerce").fillna(0) >= 1]
    print(f"{len(df)} wet gauges; {len(rare)} reached >=1-yr I15 recurrence "
          f"(max {df['gauge'].iloc[0]}: I15 {df['i15'].iloc[0]:.0f} mm/h, "
          f"RI {_ri_str(df['RI_i15'].iloc[0])} yr)")


def _cmd_smooth(args):
    """Smooth an event's radar fields: a methods×radii comparison figure, an
    optional radar-gauge skill sweep, and optional writing of smoothed fields.

    Reuses an already-processed event's observed fields
    (``<from-key>_<field>.tif`` from MRMS or NEXRAD) -- no radar re-run."""
    from .plot import smoothing_comparison, smoothing_skill_plot
    from .smoothing import best_radius, gauge_skill_sweep, smooth_event_fields
    os.makedirs(args.out_dir, exist_ok=True)
    src, from_key = args.from_dir, args.from_key
    key = args.key or from_key

    field_tif = find(src, f"{from_key}_{args.field}.tif")
    if not os.path.exists(field_tif):
        sys.exit(f"error: {field_tif} not found (need the observed "
                 f"{args.field} field)")

    hs_path = args.hillshade or find(src, f"{from_key}_hillshade.tif")
    if not os.path.exists(hs_path):
        print(f"note: hillshade {hs_path} not found; rendering without terrain")
        hs_path = None
    # gauge overlay: --gauges-file overrides the auto-located <from-key>_gauges
    # (so a NEXRAD field dir can borrow the event's MRMS-keyed gauges file)
    overlay_path = getattr(args, "gauges_file", None) or \
        find(src, f"{from_key}_gauges.geojson")
    have_gauges = args.gauges and os.path.exists(overlay_path)
    if args.gauges and not have_gauges:
        print(f"note: {overlay_path} not found; skipping gauge overlay")
    gauges_arg = overlay_path if have_gauges else None

    clip_arg = None                                   # avoid Albers->UTM wedges
    if args.clip:
        if args.perimeters:
            clip_arg = args.perimeters
        else:                                         # clip to the field footprint
            import geopandas as gpd
            import rioxarray
            from .aoi import bbox_polygon
            with rioxarray.open_rasterio(field_tif) as d:
                b = d.rio.bounds()
            clip_arg = gpd.GeoDataFrame(
                geometry=[bbox_polygon((b[0], b[1], b[2], b[3]))], crs=4326)

    # 1) methods x radii comparison figure
    if not args.no_comparison:
        out = out_path(args.out_dir, f"{key}_smoothing_compare.png")
        smoothing_comparison(
            src, from_key, field=args.field, methods=tuple(args.methods),
            radii_km=tuple(args.radii), hillshade=hs_path,
            perimeters=args.perimeters, gauges=gauges_arg,
            reference=args.reference, local_roads=args.local_roads,
            label_reference=not args.no_reference_labels, alpha=args.alpha,
            hillshade_alpha=args.hillshade_alpha, cmap=args.cmap,
            shared_scale=not args.no_shared_scale, power=args.power,
            clip=clip_arg, clip_margin=args.clip_margin, title=args.title,
            out_path=out, dpi=args.dpi)
        print(f"wrote {out}")

    # 2) gauge-skill sweep: does smoothing improve radar-gauge agreement?
    if args.gauge_analysis:
        gauges = _analysis_gauges(args, src, from_key)
        if not len(gauges):
            sys.exit("no gauges available for --gauge-analysis")
        sweep = gauge_skill_sweep(
            gauges, src, from_key, methods=tuple(args.methods),
            radii_km=tuple(args.sweep), durations=tuple(args.durations),
            rqi_min=args.rqi_min, max_report_min=args.max_report_min,
            power=args.power)
        csv = out_path(args.out_dir, f"{key}_smoothing_skill.csv")
        sweep.to_csv(csv, index=False)
        png = out_path(args.out_dir, f"{key}_smoothing_skill.png")
        smoothing_skill_plot(sweep, durations=tuple(args.durations),
                             out_path=png, dpi=args.dpi)
        print(f"wrote {csv}\nwrote {png}")
        for d in args.durations:
            br = best_radius(sweep, by="rmse", duration=d)
            bc = best_radius(sweep, by="corr", duration=d)
            print(f"  I{d} optimum radius  min-RMSE: " + ", ".join(
                f"{m} {r:.1f}km" for m, (r, _) in br.items()))
            print(f"  I{d} optimum radius  max-corr: " + ", ".join(
                f"{m} {r:.1f}km" for m, (r, _) in bc.items()))

    # 3) optionally write smoothed fields at a chosen method/radius
    if args.write:
        if args.write_radius is None:
            sys.exit("--write requires --write-radius KM")
        paths = smooth_event_fields(src, from_key, args.out_dir, key,
                                    args.write, args.write_radius,
                                    power=args.power)
        print(f"wrote {len(paths)} smoothed field(s) ({args.write}, "
              f"r={args.write_radius:g}km) to {args.out_dir}/")


def _add_window_opts(p):
    """``--start``/``--end``: the analysis window, for commands without gauge
    options. Where a command has both, these are the *same* flags -- one pair
    meaning "the analysis window" for the radar stack and the gauges alike.
    """
    p.add_argument("--start", help="UTC window start YYYYMMDDHHMM (with --end, "
                                   "overrides --date; use when a storm does not "
                                   "line up with a local day)")
    p.add_argument("--end", help="UTC window end YYYYMMDDHHMM")


def _add_gauge_opts(p):
    p.add_argument("--token", help="Synoptic API token (else $SYNOPTIC_TOKEN)")
    p.add_argument("--durations", type=int, nargs="+", default=[15, 30, 60],
                   metavar="MIN",
                   help="peak-intensity durations, minutes (default 15 30 60)")
    p.add_argument("--start", help="UTC window start YYYYMMDDHHMM (with --end, "
                                   "overrides --date)")
    p.add_argument("--end", help="UTC window end YYYYMMDDHHMM")


def _add_gauges_pipeline_opts(p):
    """Virtual-gauge pipeline options for the unified ``gauges`` command (fetch ->
    store -> rainfall comparison atlas + per-gauge detail figures). Defaults run the
    full 3-way (MRMS+NEXRAD+real) pipeline on the wet, near-AOI stations."""
    p.add_argument("--store-only", action="store_true",
                   help="fetch + write only the gauge store; skip the virtual-gauge "
                        "atlas + detail figures (the pre-2026-06-29 behaviour)")
    p.add_argument("--no-atlas", action="store_true",
                   help="skip the rainfall comparison atlas (<key>_vg_atlas.png)")
    p.add_argument("--no-detail", action="store_true",
                   help="skip the per-gauge detail figures (VirtualGaugeFigures/)")
    p.add_argument("--no-nexrad", dest="nexrad", action="store_false",
                   help="MRMS + real only; default adds single-radar NEXRAD (3-way)")
    p.add_argument("--all-gauges", dest="wet_only", action="store_false",
                   help="include dry stations in the atlas/detail "
                        "(default: wet stations only)")
    p.add_argument("--wet-min", type=float, default=0.5, metavar="MMPH",
                   help="peak intensity (mm/h) below which a station counts as dry "
                        "(default 0.5; drops traces)")
    p.add_argument("--max-dist-km", type=float, default=None,
                   help="also drop stations farther than this many km from the AOI")
    p.add_argument("--max-report-min", type=float, default=None, metavar="MIN",
                   help="keep only gauges whose native reporting interval is <= this "
                        "many minutes (e.g. 60 = hourly or finer); drops daily/coarse "
                        "reporters whose i15/i30 are smeared by the 1-min interpolation")
    p.add_argument("--atlas-metric", type=int, default=15,
                   help="intensity duration (min) shown in the atlas panels (default 15)")
    p.add_argument("--no-multisensor", action="store_true",
                   help="skip the hourly gauge-corrected MultiSensor QPE (MRMS) overlay")
    p.add_argument("--radar", help="NEXRAD radar id (nearest to the AOI if unset)")
    p.add_argument("--method", choices=["za", "kdp"], default="kdp",
                   help="NEXRAD rate recipe for the VG series (default kdp)")
    p.add_argument("--cache-dir", help="NEXRAD volume cache dir "
                                       "(default <out-dir>/nexrad_cache)")
    p.add_argument("--pad-min", type=int, default=30,
                   help="padding (min) around the auto-detected storm rain window")
    p.add_argument("--dpi", type=int, default=200)
    # full 3-way pipeline on wet stations, MRMS primary; 5/15/30/60-min VG intervals
    p.set_defaults(nexrad=True, wet_only=True, source="mrms",
                   durations=[5, 15, 30, 60])


def _add_aoi(p):
    p.add_argument("--bbox", type=float, nargs=4,
                   metavar=("W", "S", "E", "N"),
                   help="AOI bounding box in lon/lat degrees")
    p.add_argument("--aoi", help="AOI vector file (GeoJSON/SHP/GPKG)")
    p.add_argument("--pad-deg", type=float, default=0.05,
                   help="degrees to pad the AOI bounds (default 0.05)")
    p.add_argument("--dst-crs", default="EPSG:5070",
                   help="working CRS for DEM/figure (default EPSG:5070)")
    p.add_argument("--out-dir", default="stormscape_out", help="output directory")
    p.add_argument("--key", help="filename stem for outputs")
    _add_layout(p)


def _add_dem_opts(p):
    """DEM-fetch knobs shared by `dem`, `run` and `zoom --refine-dem`.

    The resampling choices are the ones that make sense on a continuous
    elevation surface; `mode`/`sum`/`q1` and friends are valid rasterio
    resamplers but meaningless for a DEM, so they are not offered. `nearest` IS
    offered despite being the artefact this project went to some trouble to
    remove -- reproducing the corduroy hatch is the control experiment that
    pins the blame on nearest rather than on the number of warps.
    """
    p.add_argument("--resampling", default=None,
                   choices=["nearest", "bilinear", "cubic", "cubic_spline",
                            "lanczos", "average"],
                   help="resampling for the 3DEP warp (default bilinear; cubic "
                        "overshoots at cliffs and hillshading turns that into "
                        "bright rims, nearest aliases into a ~45 deg hatch)")


def _add_layout(p):
    """``--flat`` opts out of the sorted figures/rasters/tables/vectors layout.

    Reading is unaffected -- every ``--from-dir`` / ``--radar-dir`` resolves
    both layouts -- so this only changes where a run *writes*.
    """
    p.add_argument("--flat", dest="layout", action="store_const",
                   const="flat", default=None,
                   help="write products straight into --out-dir instead of "
                        "sorting them into figures/ rasters/ tables/ vectors/ "
                        "(also settable with $STORMSCAPE_LAYOUT=flat)")


def _add_overlays(p):
    p.add_argument("--perimeters", help="vector overlay outlined as a border")
    p.add_argument("--basins", help="vector overlay drawn as thin outlines")
    p.add_argument("--highlight", help="vector overlay drawn bold (cyan)")
    p.add_argument("--points", help="point vector overlay (triangles)")
    p.add_argument("--title", help="figure title")
    p.add_argument("--cmap", default="YlGnBu",
                   help="colormap for the i15 field (default YlGnBu, "
                        "colorblind-safe; e.g. inferno, cmc.lajolla, cmc.oslo)")
    p.add_argument("--wet-min", type=float, default=5.0,
                   help="mask i15 below this mm/h in the drape (default 5)")
    p.add_argument("--basemap", action="store_true",
                   help="underlay an open-source basemap (named creeks/rivers,"
                        " roads, place names) via contextily")
    p.add_argument("--basemap-provider", default="USGS.USTopo",
                   help="contextily provider key (default USGS.USTopo = USGS "
                        "National Map topo; try USGS.USImageryTopo or "
                        "OpenStreetMap.Mapnik)")
    p.add_argument("--basemap-labels", default=None,
                   help="optional labels-only provider drawn on top (e.g. "
                        "CartoDB.PositronOnlyLabels) for a label-free base; "
                        "'none' to force off")
    p.add_argument("--basemap-zoom", type=int, default=None,
                   help="basemap tile zoom level (default: auto)")
    p.add_argument("--hillshade-alpha", type=float, default=None,
                   help="hillshade opacity (default 1.0; 0.0 when --basemap)")
    p.add_argument("--alpha", type=float, default=None,
                   help="field drape opacity over terrain/basemap (default 0.32)")
    p.add_argument("--reference", action="store_true",
                   help="auto-fetch + overlay labelled vector streams (NHD), "
                        "roads (TIGER), and place names (GNIS) on the hillshade")
    p.add_argument("--local-roads", action="store_true",
                   help="include residential/local roads in --reference")
    p.add_argument("--no-reference-labels", action="store_true",
                   help="draw reference lines/points without text labels")
    p.add_argument("--clip", action="store_true",
                   help="clip the figure tightly to the AOI / --perimeters extent")
    p.add_argument("--clip-margin", type=float, default=0.04,
                   help="fractional margin around the clip extent (default 0.04)")
    p.add_argument("--dpi", type=int, default=200,
                   help="output figure resolution (default 200)")


def _add_climate_opts(p):
    """NOAA Atlas 14 climate-comparison knobs, shared by `climate` and `zoom`."""
    p.add_argument("--ari", type=int, default=1,
                   help="recurrence interval in years (default 1)")
    p.add_argument("--durations", type=int, nargs="+", default=[15, 30, 60],
                   metavar="MIN",
                   help="intensity durations, minutes (default 15 30 60)")
    p.add_argument("--region",
                   help="NOAA Atlas 14 region code override (else auto from the "
                        "AOI; e.g. sw, tx, se, mw, ne, orb, inw)")
    p.add_argument("--stat", default="mean", choices=["mean", "lower", "upper"],
                   help="grid statistic: mean (default) or 90%% CI lower/upper")
    p.add_argument("--anomaly-cmap", default="cmc.vik",
                   help="diverging colormap for the anomaly maps (default cmc.vik)")
    p.add_argument("--shared-row-scale", action="store_true",
                   help="share the colour scale within each comparison row "
                        "(default: each panel scaled independently)")
    p.add_argument("--obs-smooth", default="gaussian",
                   choices=["gaussian", "uniform", "median", "idw", "none"],
                   help="smooth the radar field for display (default gaussian; "
                        "'none' = raw peaky field). climate: the observed field "
                        "in the comparison + anomaly; zoom: also the i15 map + "
                        "rainfall panels")
    p.add_argument("--obs-smooth-radius", type=float, default=1.0, metavar="KM",
                   help="display smoothing radius, km (default 1.0; nominal "
                        "scale ~ Gaussian sigma)")
    p.add_argument("--eps", type=float, default=0.1,
                   help="climatology floor (mm/h) below which the anomaly is masked")
    p.add_argument("--no-comparison", action="store_true",
                   help="skip the observed-vs-climatology comparison figure")
    p.add_argument("--no-anomaly", action="store_true",
                   help="skip the per-duration anomaly maps")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="stormscape",
                                 description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    pd_ = sub.add_parser("dem", help="download DEM + hillshade")
    _add_aoi(pd_)
    pd_.add_argument("--resolution", type=int, default=10)
    pd_.add_argument("--clip-dem", action="store_true",
                     help="mask DEM outside the AOI polygon")
    _add_dem_opts(pd_)
    pd_.set_defaults(func=_cmd_dem)

    pi = sub.add_parser("i15", help="stack MRMS into a peak-i15 field")
    _add_aoi(pi)
    pi.add_argument("--date", help="storm-day YYYYMMDD (scans [04Z, next-day "
                                   "10Z]); optional if --start/--end are given")
    pi.add_argument("--qpe-thresh", type=float, default=2.5)
    pi.add_argument("--max-wet-hours", type=int, default=8)
    pi.add_argument("--workers", type=int, default=12)
    _add_window_opts(pi)
    pi.add_argument("--multisensor", action="store_true",
                    help="also fetch gauge-corrected MultiSensor QPE total "
                         "-> <key>_mstotal.tif")
    pi.set_defaults(func=_cmd_i15)

    pm = sub.add_parser("map", help="drape an i15 tif over a hillshade tif")
    pm.add_argument("--hillshade", required=True)
    pm.add_argument("--i15", required=True)
    pm.add_argument("--out", help="output PNG path")
    pm.add_argument("--out-dir", default=".")
    pm.add_argument("--key")
    pm.add_argument("--dst-crs", default="EPSG:5070")
    _add_layout(pm)                     # only bites when --out is not given
    _add_overlays(pm)
    pm.set_defaults(func=_cmd_map)

    pr = sub.add_parser("run", help="DEM -> hillshade -> i15 -> figure")
    _add_aoi(pr)
    pr.add_argument("--date", help="storm-day YYYYMMDD (scans [04Z, next-day "
                                   "10Z]); optional if --start/--end are given")
    pr.add_argument("--resolution", type=int, default=10)
    pr.add_argument("--clip-dem", action="store_true",
                    help="mask DEM outside the AOI polygon")
    _add_dem_opts(pr)
    pr.add_argument("--qpe-thresh", type=float, default=2.5)
    pr.add_argument("--max-wet-hours", type=int, default=8)
    pr.add_argument("--workers", type=int, default=12)
    pr.add_argument("--multisensor", action="store_true",
                    help="also fetch gauge-corrected MultiSensor QPE total")
    _add_overlays(pr)
    _add_gauge_opts(pr)
    pr.add_argument("--gauges", action="store_true",
                    help="fetch Synoptic rain gauges + overlay them on the figure")
    pr.add_argument("--compare", action="store_true",
                    help="also compute radar-vs-gauge stats + a residual map "
                         "(implies --gauges)")
    pr.add_argument("--rqi-min", type=float, default=None,
                    help="ignore gauges whose radar RQI is below this (compare)")
    pr.add_argument("--max-report-min", type=float, default=None,
                    help="for i15/i30/i60: ignore gauges coarser than this "
                         "native reporting interval, min (compare)")
    pr.set_defaults(func=_cmd_run)

    pg = sub.add_parser("gauges",
                        help="Synoptic gauges, full pipeline: canonical store "
                             "(<key>_gauges.geojson + RainGaugeData/) + virtual-gauge "
                             "comparison atlas + per-gauge detail figures")
    _add_aoi(pg)
    pg.add_argument("--date", help="storm-day YYYYMMDD (sets the UTC window)")
    pg.add_argument("--no-series", action="store_true",
                    help="write only the geojson, skip the per-gauge series CSVs")
    _add_gauge_opts(pg)
    _add_gauges_pipeline_opts(pg)
    pg.set_defaults(func=_cmd_gauges)

    pc = sub.add_parser("compare",
                        help="radar rasters vs gauges: residuals + skill stats")
    _add_aoi(pc)
    pc.add_argument("--date", help="storm-day YYYYMMDD (for live gauge fetch)")
    _add_gauge_opts(pc)
    pc.add_argument("--gauges",
                    help="precomputed gauges file (GeoJSON); skips live fetch")
    pc.add_argument("--radar-dir",
                    help="dir with <key>_*.tif radar rasters (default --out-dir)")
    pc.add_argument("--rqi-min", type=float, default=None,
                    help="ignore gauges whose radar RQI is below this")
    pc.add_argument("--max-report-min", type=float, default=None,
                    help="for i15/i30/i60 only: ignore gauges whose native "
                         "reporting interval (min) exceeds this, e.g. 15 "
                         "(coarse reporters smear short-duration peaks)")
    pc.add_argument("--multisensor", action="store_true",
                    help="add a gauge-corrected MultiSensor total row "
                         "(needs <radar-dir>/<key>_mstotal.tif)")
    pc.add_argument("--out", help="output CSV path for the per-gauge table")
    pc.add_argument("--map", help="also render a residual-map PNG to this path")
    pc.add_argument("--map-metric", default="i15max",
                    help="radar field for --map residuals (default i15max)")
    pc.add_argument("--hillshade",
                    help="hillshade tif for the --map background (optional)")
    pc.add_argument("--perimeters", help="AOI/perimeter overlay for --map")
    pc.add_argument("--reference", action="store_true",
                    help="labelled vector reference overlay on --map")
    pc.add_argument("--clip", action="store_true",
                    help="clip --map tightly to --perimeters")
    pc.add_argument("--title", help="title for --map")
    pc.add_argument("--cmap", default="YlGnBu",
                    help="colormap for the i15 field background (default YlGnBu)")
    pc.add_argument("--alpha", type=float, default=None,
                    help="field drape opacity on --map (default 0.32)")
    pc.set_defaults(func=_cmd_compare)

    pn = sub.add_parser("nexrad",
                        help="single-radar NEXRAD Level II reflectivity field")
    _add_aoi(pn)
    pn.add_argument("--date", help="storm-day YYYYMMDD")
    pn.add_argument("--time", help="UTC scan time HHMM (single scan; "
                                   "default: storm-window midpoint)")
    pn.add_argument("--start", help="UTC window start YYYYMMDDHHMM (--composite)")
    pn.add_argument("--end", help="UTC window end YYYYMMDDHHMM (--composite)")
    pn.add_argument("--radar", help="four-letter radar id (default: nearest to AOI)")
    pn.add_argument("--field", default="reflectivity",
                    help="pyart field: reflectivity, velocity, "
                         "cross_correlation_ratio, ... (default reflectivity)")
    pn.add_argument("--sweep", type=int, default=None,
                    help="elevation-tilt index; default = lowest (0 for a "
                         "reflectivity field; all SAILS low cuts for --intensity)")
    pn.add_argument("--composite", action="store_true",
                    help="per-cell max over all scans in the window "
                         "(storm-peak reflectivity) instead of one scan")
    pn.add_argument("--intensity", action="store_true",
                    help="build the i15/i30/i60 rainfall-intensity stack "
                         "(capped convective Z-R) instead of a reflectivity "
                         "field -> <key>_i15max.tif etc., like the MRMS i15")
    pn.add_argument("--zr-a", type=float, default=300.0,
                    help="Z-R coefficient a in Z=a*R^b (default 300, convective)")
    pn.add_argument("--zr-b", type=float, default=1.4,
                    help="Z-R exponent b (default 1.4)")
    pn.add_argument("--dbz-cap", type=float, default=53.0,
                    help="hail cap (dBZ) before Z-R conversion (default 53)")
    pn.add_argument("--no-hail-cap", action="store_true",
                    help="disable the dBZ hail cap (raw convective Z-R)")
    pn.add_argument("--method", choices=["za", "kdp"], default="za",
                    help="--intensity rate: za = capped convective Z-R (v1, all "
                         "eras); kdp = dual-pol R(Kdp) blended with Z-R (v2, 2012+)")
    pn.add_argument("--z-blend", type=float, default=35.0,
                    help="--method kdp: use R(Kdp) at/above this dBZ, capped Z-R "
                         "below (default 35)")
    pn.add_argument("--rate-cap", type=float, default=None,
                    help="--intensity: clip every method's per-scan rate to this "
                         "max (mm/h) before stacking (operational hail-cap analogue)")
    pn.add_argument("--blockage-dem",
                    help="--intensity: DEM tif for beam-blockage masking (mask "
                         "cells with cumulative blockage > --cbb-max; adds a cbb field)")
    pn.add_argument("--cbb-max", type=float, default=0.5,
                    help="cumulative beam-blockage fraction above which to mask "
                         "(default 0.5)")
    pn.add_argument("--res-m", type=float, default=500.0,
                    help="output grid resolution, metres (default 500)")
    pn.add_argument("--max-scans", type=int, default=40,
                    help="cap on scans gridded for --composite (default 40)")
    pn.add_argument("--cache-dir",
                    help="dir for downloaded volumes (default <out-dir>/nexrad_cache)")
    pn.add_argument("--hillshade", help="hillshade tif to drape over (optional)")
    pn.add_argument("--gauges", action="store_true",
                    help="fetch Synoptic gauges + overlay, coloured by the radar "
                         "value sampled at each gauge (needs $SYNOPTIC_TOKEN)")
    pn.add_argument("--gauges-file",
                    help="precomputed gauges GeoJSON to overlay (skips live fetch)")
    pn.add_argument("--token", help="Synoptic API token (else $SYNOPTIC_TOKEN)")
    pn.add_argument("--durations", type=int, nargs="+", default=[15, 30, 60],
                    metavar="MIN", help="gauge peak-intensity durations, minutes")
    _add_overlays(pn)
    pn.set_defaults(func=_cmd_nexrad)

    pp = sub.add_parser("panels",
                        help="multi-panel diagnostic map of stacked radar fields")
    pp.add_argument("--radar-dir", default=".",
                    help="dir with <key>_<field>.tif (default .)")
    pp.add_argument("--key", required=True, help="filename stem of the fields")
    pp.add_argument("--fields", nargs="+",
                    default=["tpki15", "total", "rqi", "shsr"],
                    help="field names to panel (default tpki15 total rqi shsr)")
    pp.add_argument("--perimeters", help="AOI/perimeter outline overlay")
    pp.add_argument("--gauges", help="gauges GeoJSON to overlay as points")
    pp.add_argument("--hillshade", help="hillshade tif to drape the panels over")
    pp.add_argument("--reference", action="store_true",
                    help="labelled NHD streams / TIGER roads / GNIS places overlay")
    pp.add_argument("--local-roads", action="store_true",
                    help="include residential roads in --reference")
    pp.add_argument("--no-reference-labels", action="store_true",
                    help="draw reference lines/points without labels")
    pp.add_argument("--clip", action="store_true",
                    help="clip panels tightly to the --perimeters extent")
    pp.add_argument("--clip-margin", type=float, default=0.04)
    pp.add_argument("--alpha", type=float, default=None,
                    help="field drape opacity (default 0.32)")
    pp.add_argument("--out", help="output PNG (default <radar-dir>/<key>_panels.png)")
    pp.add_argument("--title", help="figure title")
    pp.add_argument("--dpi", type=int, default=200)
    pp.set_defaults(func=_cmd_panels)

    pv = sub.add_parser("vgauge",
                        help="virtual rain gauges: radar rainfall time series at point(s)")
    _add_aoi(pv)
    pv.add_argument("--date", help="storm-day YYYYMMDD (sets the UTC window)")
    pv.add_argument("--start", help="UTC window start YYYYMMDDHHMM (with --end)")
    pv.add_argument("--end", help="UTC window end YYYYMMDDHHMM")
    pv.add_argument("--point", action="append", metavar="LON,LAT[,NAME]",
                    help="virtual-gauge location, repeatable (e.g. -119.7,39.5,Site)")
    pv.add_argument("--points-file",
                    help="vector file of gauge points (GeoJSON/SHP/GPKG)")
    pv.add_argument("--durations", type=int, nargs="+", default=[5, 15, 30, 60],
                    metavar="MIN",
                    help="intensity-window durations, minutes (default 5 15 30 60)")
    pv.add_argument("--no-multisensor", action="store_true",
                    help="skip the hourly gauge-corrected MultiSensor QPE overlay")
    pv.add_argument("--gauges", action="store_true",
                    help="also drop a VG at every real Synoptic station + overlay "
                         "the real gauge series (reuses a saved store if --from-dir "
                         "given, else needs $SYNOPTIC_TOKEN)")
    pv.add_argument("--from-dir",
                    help="reuse a canonical gauge store from here (its "
                         "<from-key>_gauges.geojson + RainGaugeData/) instead of "
                         "re-fetching; the full storm-day record is trimmed to the "
                         "rain window")
    pv.add_argument("--from-key", help="key/stem of the saved store (with --from-dir)")
    pv.add_argument("--refetch", action="store_true",
                    help="ignore any saved store and fetch live from Synoptic")
    pv.add_argument("--pad-min", type=int, default=30,
                    help="minutes of padding around the auto-detected rain window "
                         "when trimming a reused store (default 30)")
    pv.add_argument("--max-dist-km", type=float, default=None,
                    help="with --gauges: drop stations farther than this many km "
                         "from the AOI (--bbox/--aoi, else the <from-key>_i15max "
                         "footprint); keeps the atlas to near-AOI gauges")
    pv.add_argument("--wet-only", action="store_true",
                    help="with --gauges: drop stations whose peak atlas-metric "
                         "intensity is below --wet-min (declutters the atlas)")
    pv.add_argument("--wet-min", type=float, default=0.5, metavar="MMPH",
                    help="--wet-only floor: peak intensity (mm/h) below which a "
                         "station counts as dry (default 0.5; drops traces)")
    pv.add_argument("--max-report-min", type=float, default=None, metavar="MIN",
                    help="with --gauges: keep only gauges reporting at this cadence "
                         "or finer (e.g. 60 = hourly or finer); drops coarse reporters")
    pv.add_argument("--token", help="Synoptic API token (else $SYNOPTIC_TOKEN)")
    pv.add_argument("--atlas", action="store_true",
                    help="render an atlas subplot of all virtual gauges (vs real)")
    pv.add_argument("--atlas-metric", type=int, default=15,
                    help="duration (min) shown in the atlas panels (default 15)")
    pv.add_argument("--detail", action="store_true",
                    help="also write a big 4-row figure per gauge (cumulative + "
                         "I60/I30/I15) to VirtualGaugeFigures/")
    pv.add_argument("--source", choices=["mrms", "nexrad"], default="mrms",
                    help="primary rainfall source: mrms (default) or nexrad")
    pv.add_argument("--nexrad", action="store_true",
                    help="also include the single-radar NEXRAD L2 series "
                         "alongside MRMS on the atlas/CSVs (>=2020 events)")
    pv.add_argument("--radar", help="NEXRAD radar id for --source nexrad (nearest if unset)")
    pv.add_argument("--method", choices=["za", "kdp"], default="kdp",
                    help="--source nexrad rate recipe (default kdp)")
    pv.add_argument("--cache-dir", help="NEXRAD volume cache dir (--source nexrad)")
    pv.add_argument("--dpi", type=int, default=200)
    pv.set_defaults(func=_cmd_vgauge)

    pz = sub.add_parser("zoom",
                        help="re-render an existing event's figures clipped to a "
                             "sub-AOI (reuses rasters; no MRMS re-download)")
    _add_aoi(pz)                        # the ZOOM sub-AOI + new --out-dir / --key
    pz.add_argument("--from-dir", required=True,
                    help="directory of the already-processed event (its --out-dir)")
    pz.add_argument("--from-key", required=True,
                    help="key/stem of the already-processed event")
    pz.add_argument("--fields", nargs="+",
                    default=["tpki15", "total", "rqi", "shsr"],
                    help="panel field stems to re-render (default tpki15 total rqi shsr)")
    pz.add_argument("--gauges", action="store_true",
                    help="overlay the source <from-key>_gauges.geojson if present")
    pz.add_argument("--refine-dem", action="store_true",
                    help="re-fetch a finer DEM + hillshade for the zoom extent "
                         "(the only product that benefits; MRMS has no finer res)")
    pz.add_argument("--resolution", type=int, default=10,
                    help="DEM resolution for --refine-dem, metres (default 10)")
    pz.add_argument("--clip-dem", action="store_true",
                    help="with --refine-dem, mask the DEM outside the zoom polygon")
    _add_dem_opts(pz)
    pz.add_argument("--crop-rasters", action="store_true",
                    help="also write cropped copies of the source GeoTIFFs to the "
                         "new folder (a self-contained zoom)")
    pz.add_argument("--no-map", action="store_true", help="skip the main i15 map")
    pz.add_argument("--no-panels", action="store_true",
                    help="skip the diagnostic-panels figure")
    pz.add_argument("--no-climate", action="store_true",
                    help="skip the NOAA Atlas 14 climatology + anomaly maps "
                         "(produced by default for the zoom sub-AOI)")
    _add_climate_opts(pz)              # climate knobs (--ari/--durations/--obs-smooth…)
    _add_overlays(pz)                  # note: zoom always clips to the sub-AOI
    pz.set_defaults(func=_cmd_zoom)

    pk = sub.add_parser("pick",
                        help="interactive browser bbox picker: drag a zoom box "
                             "on an event's map (cross-platform, no GUI toolkit)")
    pk.add_argument("--from-dir", required=True,
                    help="directory of the already-processed event")
    pk.add_argument("--from-key", required=True,
                    help="key/stem of the already-processed event")
    pk.add_argument("--i15",
                    help="i15 GeoTIFF (default <from-dir>/<from-key>_i15max.tif)")
    pk.add_argument("--hillshade",
                    help="hillshade GeoTIFF (default <from-dir>/<from-key>_hillshade.tif)")
    pk.add_argument("--out", help="output HTML (default <from-dir>/<from-key>_pick.html)")
    pk.add_argument("--cmap", default="YlGnBu",
                    help="colormap for the i15 background (default YlGnBu)")
    pk.add_argument("--wet-min", type=float, default=5.0,
                    help="mask i15 below this mm/h in the background (default 5)")
    pk.add_argument("--no-reference", action="store_true",
                    help="skip the labelled NHD/TIGER/GNIS overlay (on by default "
                         "for orientation; disabling skips a network fetch)")
    pk.add_argument("--local-roads", action="store_true",
                    help="include residential/local roads in the reference overlay")
    pk.add_argument("--no-reference-labels", action="store_true",
                    help="draw reference lines/points without text labels")
    pk.add_argument("--perimeters", help="AOI/perimeter outline to draw on the map")
    pk.add_argument("--gauges", action="store_true",
                    help="overlay the source <from-key>_gauges.geojson if present")
    pk.add_argument("--no-open", action="store_true",
                    help="write the HTML but don't open a browser")
    pk.set_defaults(func=_cmd_pick)

    pcl = sub.add_parser("climate",
                         help="NOAA Atlas 14 rainfall climatology vs observed "
                              "i15/i30/i60: comparison + anomaly maps (no re-run)")
    _add_aoi(pcl)              # optional AOI override + --out-dir/--key/--pad-deg
    pcl.add_argument("--from-dir", required=True,
                     help="directory of the already-processed event (its --out-dir)")
    pcl.add_argument("--from-key", required=True,
                     help="key/stem of the already-processed event")
    pcl.add_argument("--hillshade",
                     help="hillshade GeoTIFF (default <from-dir>/<from-key>_hillshade.tif)")
    pcl.add_argument("--gauges", action="store_true",
                     help="overlay the source <from-key>_gauges.geojson if present")
    _add_climate_opts(pcl)
    _add_overlays(pcl)
    pcl.set_defaults(func=_cmd_climate)

    pe = sub.add_parser("export",
                        help="georeferenced exports for GIS/CalTopo: EPSG:3857 "
                             "GeoTIFFs (raw + colorized RGBA) + GeoPDF figures "
                             "(no re-run)")
    _add_aoi(pe)              # optional AOI override + --out-dir/--key/--pad-deg
    pe.add_argument("--from-dir", required=True,
                    help="directory of the already-processed event (its --out-dir)")
    pe.add_argument("--from-key", required=True,
                    help="key/stem of the already-processed event")
    pe.add_argument("--layers", nargs="+", default=list(DEFAULT_EXPORT_FIELDS),
                    metavar="FIELD",
                    help="field suffixes to export to GeoTIFF (default: anom_i15 "
                         "i15max; e.g. add i30max i60max total peakrate_mmph)")
    pe.add_argument("--crs", default="EPSG:3857",
                    help="GeoTIFF export CRS (default EPSG:3857, CalTopo's native)")
    pe.add_argument("--no-geotiff", action="store_true",
                    help="skip the EPSG:3857 GeoTIFF layers")
    pe.add_argument("--no-rgb", action="store_true",
                    help="skip the colorized RGBA GeoTIFFs (raw float only)")
    pe.add_argument("--no-raw", action="store_true",
                    help="skip the raw float GeoTIFFs (colorized RGBA only)")
    pe.add_argument("--no-figures", action="store_true",
                    help="skip the GeoPDF figures (export only the GeoTIFF layers)")
    pe.add_argument("--no-i15", action="store_true",
                    help="skip the i15-map GeoPDF")
    pe.add_argument("--no-anom", action="store_true",
                    help="skip the anomaly-map GeoPDF")
    pe.add_argument("--pdf-crs", default="UTM",
                    help="projected CRS to render the GeoPDF figures in (default "
                         "UTM = identical look to the PNG deliverables; pass "
                         "EPSG:3857 to match the GeoTIFF layers)")
    pe.add_argument("--anom-duration", type=int, default=15, metavar="MIN",
                    help="anomaly duration for the GeoPDF (default 15; needs "
                         "<from-key>_anom_i<MIN>.tif from `climate`)")
    pe.add_argument("--ari", type=int, default=1,
                    help="climatology ARI (years) labelling the anomaly (default 1)")
    pe.add_argument("--anomaly-cmap", default="cmc.vik",
                    help="diverging colormap for the anomaly (default cmc.vik)")
    pe.add_argument("--hillshade",
                    help="hillshade GeoTIFF (default <from-dir>/<from-key>_hillshade.tif)")
    pe.add_argument("--gauges", action="store_true",
                    help="overlay the source <from-key>_gauges.geojson if present")
    pe.add_argument("--streams", action="store_true",
                    help="also export the full-resolution NHD stream network for "
                         "the AOI as a vector layer (CalTopo-ready)")
    pe.add_argument("--streams-format", default="geojson",
                    choices=["geojson", "gpkg", "shp", "kml"],
                    help="vector format for --streams (default geojson)")
    pe.add_argument("--streams-bbox", action="store_true",
                    help="keep whole flowlines over the AOI bounding box "
                         "(default: clip/truncate to the AOI polygon)")
    pe.add_argument("--streams-named-only", action="store_true",
                    help="keep only GNIS-named creeks/rivers (default: the full "
                         "network incl. unnamed headwaters)")
    _add_overlays(pe)
    pe.set_defaults(func=_cmd_export)

    prc = sub.add_parser("recurrence",
                         help="wet-gauge table: peak I15/I30/I60 + time-of-peak, "
                              "anomaly + recurrence interval vs NOAA Atlas 14")
    _add_aoi(prc)             # provides --out-dir/--key (AOI itself is unused;
    prc.add_argument("--from-dir", required=True,   # PFDS locates by gauge lat/lon)
                     help="directory of the processed event (its --out-dir)")
    prc.add_argument("--from-key", required=True,
                     help="key/stem of the processed event")
    prc.add_argument("--gauges-file",
                     help="gauges geojson (default <from-dir>/<from-key>_gauges.geojson)")
    prc.add_argument("--durations", type=int, nargs="+", default=[15, 30, 60],
                     metavar="MIN", help="intensity durations (default 15 30 60)")
    prc.add_argument("--series", default="pds", choices=["pds", "ams"],
                     help="PF series: pds (partial-duration, default; matches the "
                          "maps) or ams (annual-maximum)")
    prc.add_argument("--stat", default="mean", choices=["mean", "lower", "upper"],
                     help="PF estimate: mean (default) or 90%% CI lower/upper")
    prc.add_argument("--raingauge-dir",
                     help="dir of per-gauge series CSVs for the I15 time-of-peak "
                          "(default <from-dir>/RainGaugeData)")
    prc.add_argument("--no-peak-time", action="store_true",
                     help="skip the I15 time-of-peak column")
    prc.add_argument("--max-dist-km", type=float, default=None,
                     help="drop gauges farther than this many km from the AOI "
                          "(default: keep all, flagged with dist_to_aoi_km/in_aoi)")
    prc.set_defaults(func=_cmd_recurrence)

    psm = sub.add_parser("smooth",
                         help="smooth an event's radar fields: methods×radii "
                              "comparison + optional gauge-skill analysis (no re-run)")
    _add_aoi(psm)             # AOI override (gauge fetch) + --out-dir/--key/--pad-deg
    psm.add_argument("--from-dir", required=True,
                     help="directory of the already-processed event (its --out-dir)")
    psm.add_argument("--from-key", required=True,
                     help="key/stem of the already-processed event")
    psm.add_argument("--field", default="i15max",
                     help="field to smooth/compare (default i15max; e.g. i30max, "
                          "i60max, total)")
    psm.add_argument("--methods", nargs="+",
                     default=["gaussian", "uniform", "median", "idw"],
                     choices=["gaussian", "uniform", "median", "idw"],
                     help="smoothing methods to compare (default: all four)")
    psm.add_argument("--radii", type=float, nargs="+", default=[0, 1, 2, 4],
                     metavar="KM",
                     help="smoothing radii (km) for the comparison columns "
                          "(default 0 1 2 4; 0 = raw)")
    psm.add_argument("--power", type=float, default=2.0,
                     help="IDW distance power (default 2)")
    psm.add_argument("--no-shared-scale", action="store_true",
                     help="scale each comparison panel independently (default: one "
                          "shared colour scale so peak-flattening is visible)")
    psm.add_argument("--hillshade",
                     help="hillshade GeoTIFF (default <from-dir>/<from-key>_hillshade.tif)")
    psm.add_argument("--gauge-analysis", action="store_true",
                     help="sweep smoothing radius vs radar-gauge skill "
                          "(corr/RMSE/bias-ratio) -> skill figure + CSV")
    psm.add_argument("--gauges", action="store_true",
                     help="overlay the source <from-key>_gauges.geojson on the maps")
    psm.add_argument("--gauges-file",
                     help="gauges GeoJSON for --gauge-analysis (overrides the "
                          "source <from-key>_gauges.geojson)")
    psm.add_argument("--sweep", type=float, nargs="+",
                     default=[0, 0.5, 1, 2, 3, 4, 6, 8], metavar="KM",
                     help="radii (km) for the --gauge-analysis sweep "
                          "(default 0 0.5 1 2 3 4 6 8)")
    psm.add_argument("--date", help="storm-day YYYYMMDD (only for a live gauge "
                                    "fetch fallback in --gauge-analysis)")
    psm.add_argument("--rqi-min", type=float, default=None,
                     help="ignore gauges whose radar RQI is below this")
    psm.add_argument("--max-report-min", type=float, default=None,
                     help="for i15/i30/i60: ignore gauges whose native reporting "
                          "interval (min) exceeds this (coarse reporters smear peaks)")
    psm.add_argument("--write", choices=["gaussian", "uniform", "median", "idw"],
                     help="write smoothed field tifs at this method (needs "
                          "--write-radius); they flow into compare/map/climate")
    psm.add_argument("--write-radius", type=float, default=None, metavar="KM",
                     help="radius (km) for --write")
    psm.add_argument("--no-comparison", action="store_true",
                     help="skip the methods×radii comparison figure")
    _add_gauge_opts(psm)      # --token/--durations/--start/--end (fetch fallback)
    _add_overlays(psm)
    psm.set_defaults(func=_cmd_smooth)

    pb = sub.add_parser(
        "burn",
        help="near-real-time burn severity (CIMSS BRISK dNBR) over an AOI",
        description="Fetch, cache and map near-real-time burn severity. BRISK "
                    "maps every large US fire daily from a nine-satellite dNBR "
                    "composite -- an INTERIM product: supersede it with BAER "
                    "soil burn severity (--product sbs) when that lands.")
    _add_aoi(pb)
    pb.add_argument("--date", help="scar as of this date (YYYYMMDD); each fire "
                                   "is shown at its latest scene on or before "
                                   "it (default: latest available)")
    pb.add_argument("--since", help="ignore scenes older than YYYYMMDD (e.g. to "
                                    "exclude last season's fires)")
    pb.add_argument("--product", default="dnbr",
                    choices=["dnbr", "sbs", "baer_dnbr"],
                    help="dnbr = BRISK daily composite (default); sbs = BAER "
                         "soil burn severity, authoritative but only for fires "
                         "with a completed assessment; baer_dnbr = the BAER "
                         "teams' own dNBR (2025 only so far)")
    pb.add_argument("--min-age", type=float, default=None, metavar="DAYS",
                    help="require the composite to be at least DAYS old (days "
                         "since the fire entered the archive) and skip fires "
                         "that are not, naming them. A young composite has the "
                         "right pattern but reads LOW; ~14 d is where agreement "
                         "with BAER's own dNBR settles (examples/brisk_vs_baer.py). "
                         "Off by default: an immature scar still beats none, and "
                         "the run says so.")
    pb.add_argument("--scheme", default="usgs", choices=["usgs", "brisk"],
                    help="dNBR severity breaks: usgs = MTBS/USGS 0.10/0.27/0.44/"
                         "0.66 (default), brisk = the portal's 0.10/0.40/0.70")
    pb.add_argument("--fire", nargs="+", metavar="NAME",
                    help="restrict to these fire names (as listed by --list)")
    pb.add_argument("--years", type=int, nargs="+", metavar="Y",
                    help="archive years to search (default: this year and last)")
    pb.add_argument("--all-dates", action="store_true",
                    help="with --list, show every scene per fire instead of "
                         "only the latest")
    pb.add_argument("--list", action="store_true",
                    help="list the fires intersecting the AOI and stop "
                         "(headers only -- downloads nothing)")
    pb.add_argument("--cache-dir",
                    help="scene cache (default <out-dir>/brisk_cache)")
    pb.add_argument("--workers", type=int, default=12)
    pb.add_argument("--vmax", type=float, default=None,
                    help="colour-scale max for the drape (default 1.0 for dNBR, "
                         "4.0 for the sbs class field)")
    pb.add_argument("--hillshade", help="hillshade tif for the map backdrop")
    pb.add_argument("--dem", action="store_true",
                    help="fetch a DEM + hillshade for the AOI if none is found")
    pb.add_argument("--resolution", type=int, default=10,
                    help="DEM resolution for --dem (default 10 m)")
    _add_dem_opts(pb)
    pb.add_argument("--no-map", action="store_true",
                    help="write the rasters and table only, no figure")
    _add_overlays(pb)
    pb.add_argument("--continuous", action="store_true",
                    help="shade dNBR as a continuous ramp instead of banding it "
                         "into the severity classes BAER publishes (same colours)")
    # burn severity is not rainfall: the YlGnBu / 5 mm-h defaults from
    # _add_overlays would mis-scale it. The default colours are the BAER class
    # palette (verified identical across every 2025 BAER product, and the same
    # four colours BRISK publishes); the cut and scale are resolved per product
    # in _burn_display_defaults, so wet_min starts as None rather than 5 mm/h.
    pb.set_defaults(func=_cmd_burn, cmap="baer", wet_min=None)

    args = ap.parse_args(argv)
    # One switch for the whole process rather than threading `layout=` through
    # every writer: `out_path`/`subdir` read $STORMSCAPE_LAYOUT by default, so
    # setting it here makes --flat reach library calls (mrms.save_fields,
    # gauges.fetch_gauge_event, ...) without changing their call sites.
    if getattr(args, "layout", None):
        os.environ["STORMSCAPE_LAYOUT"] = args.layout
    args.func(args)


if __name__ == "__main__":
    main()
