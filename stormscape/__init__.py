"""stormscape -- terrain + radar-rainfall mapping from public US data.

A small, reusable package that:

  1. downloads USGS 3DEP / The National Map DEMs for an AOI and renders a
     hillshade (:mod:`stormscape.dem`);
  2. stacks 2-minute NOAA MRMS ``PrecipRate`` returns into a peak 15-minute
     rainfall-intensity (i15) field for a storm-day (:mod:`stormscape.mrms`);
  3. drapes the i15 field over the hillshade with optional vector overlays
     (:mod:`stormscape.plot`).

Everything is driven by an AOI (bbox tuple, vector file, or shapely geometry)
plus, for rainfall, a date -- no project-specific inputs required.

Quick start
-----------
>>> from stormscape import fetch_dem_and_hillshade, i15_storm_day, drape_i15
>>> aoi = (-105.55, 40.55, -105.25, 40.80)          # W, S, E, N (lon/lat)
>>> dem, hs = fetch_dem_and_hillshade(aoi, resolution=10,
...                                   dem_path="dem.tif",
...                                   hillshade_path="hs.tif")
>>> res = i15_storm_day(aoi, "2021-07-20")
>>> from stormscape.mrms import save_fields
>>> save_fields(res, ".", key="event")
>>> drape_i15("hs.tif", "event_i15max.tif", out_path="event.png")

or from the command line::

    python -m stormscape run --bbox -105.55 40.55 -105.25 40.80 \
        --date 2021-07-20 --resolution 10 --out-dir ./out --key cameron
"""

from .aoi import load_aoi, pad_bounds
from .atlas14 import (anomaly, climatology_field, fetch_grid, grid_url,
                      pf_point, recurrence_interval, region_for_bounds)
from .burn import (SEVERITY_SCHEMES, burn_severity, catalog as burn_catalog,
                   classify as burn_classify, fetch_scene, find_scenes,
                   register_brisk_cmap)
from .compare import (comparison_stats, compare_storm, gauge_recurrence_table,
                      radar_vs_gauge, sample_raster_at_points)
from .merge import (conditional_merge, local_bias, loo_cross_validate,
                    mean_field_bias, sample_field)
from .dem import (coverage_fraction, dem_sources, fetch_dem_and_hillshade,
                  get_dem, hillshade)
from .export import (DEFAULT_EXPORT_FIELDS, export_geotiffs, export_streams,
                     figure_to_geopdf, geopdf_supported, reproject_geotiff)
from .gauges import (fetch_gauge_event, gauge_fields, gauge_intensities,
                     gauge_timeseries, get_rainfall as gauge_rainfall,
                     get_stations as gauge_stations, load_event_series,
                     precipitation_per_minute, storm_window)
from .mines import (GROUPS as MINE_GROUPS, SOURCES as MINE_SOURCES,
                    density_grid as mine_density_grid,
                    density_raster as mine_density_raster,
                    group_counts as mine_group_counts, group_of as mine_group_of,
                    mine_features, register_source as register_mine_source)
from .streamflow import (SOURCES as STREAM_SOURCES, fetch_stream_event,
                         flow_summary, load_event_series as load_stream_series,
                         stream_series, stream_sites)
from .mrms import (compute_i15, fetch, fetch_many, i15_storm_day,
                   multisensor_total, parse_date, save_fields,
                   virtual_gauge_timeseries)
from .nexrad import (available_scans, beam_blockage, download_scans,
                     intensity_stack, lowest_tilt_grid, nearest_radar,
                     nearest_scan, radar_location, read_sweep,
                     reflectivity_composite, reflectivity_field,
                     sample_radar_at_points, z_to_rate)
from .plot import (add_gauges, add_reference, anomaly_map, bbox_picker,
                   climatology_comparison, diagnostic_panels, drape_i15,
                   plot_virtual_gauge, smoothing_comparison,
                   smoothing_skill_plot, virtual_gauge_atlas,
                   virtual_gauge_detail)
from .refdata import places, reference_layers, roads, streams
from .smoothing import (METHODS as SMOOTH_METHODS, best_radius, cell_size_km,
                        gauge_skill_sweep, smooth_array, smooth_dataarray,
                        smooth_event_fields, smooth_tif)

__version__ = "0.1.0"

__all__ = [
    "load_aoi", "pad_bounds",
    "get_dem", "hillshade", "fetch_dem_and_hillshade", "dem_sources",
    "coverage_fraction",
    "i15_storm_day", "save_fields", "compute_i15", "fetch", "fetch_many",
    "parse_date", "multisensor_total",
    "gauge_stations", "gauge_rainfall", "gauge_intensities", "gauge_fields",
    "precipitation_per_minute", "gauge_timeseries", "fetch_gauge_event",
    "load_event_series", "storm_window",
    "radar_vs_gauge", "comparison_stats", "compare_storm",
    "sample_raster_at_points", "gauge_recurrence_table",
    "mean_field_bias", "local_bias", "conditional_merge", "loo_cross_validate",
    "sample_field",
    "nearest_radar", "radar_location", "available_scans", "download_scans",
    "nearest_scan", "read_sweep", "lowest_tilt_grid", "reflectivity_field",
    "reflectivity_composite", "sample_radar_at_points", "z_to_rate",
    "intensity_stack", "beam_blockage",
    "drape_i15", "add_reference", "add_gauges", "diagnostic_panels",
    "virtual_gauge_timeseries", "plot_virtual_gauge", "virtual_gauge_atlas",
    "virtual_gauge_detail", "bbox_picker",
    "climatology_field", "anomaly", "region_for_bounds", "grid_url",
    "fetch_grid", "pf_point", "recurrence_interval",
    "climatology_comparison", "anomaly_map",
    "smooth_array", "smooth_tif", "smooth_event_fields", "smooth_dataarray",
    "gauge_skill_sweep", "best_radius", "cell_size_km", "SMOOTH_METHODS",
    "smoothing_comparison", "smoothing_skill_plot",
    "reproject_geotiff", "export_geotiffs", "export_streams",
    "figure_to_geopdf", "geopdf_supported", "DEFAULT_EXPORT_FIELDS",
    "streams", "roads", "places", "reference_layers",
    "burn_severity", "find_scenes", "burn_catalog", "fetch_scene",
    "burn_classify", "SEVERITY_SCHEMES", "register_brisk_cmap",
    "mine_features", "mine_density_grid", "mine_density_raster",
    "mine_group_counts", "mine_group_of", "register_mine_source",
    "MINE_GROUPS", "MINE_SOURCES",
    "stream_sites", "stream_series", "flow_summary", "fetch_stream_event",
    "load_stream_series", "STREAM_SOURCES",
]
