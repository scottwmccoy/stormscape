"""Shared fixtures for the stormscape test suite.

Every test here runs **offline** — no MRMS/NEXRAD/Synoptic/3DEP/Atlas 14 requests
and no API token. The suite exercises the science math, the raster plumbing, and
the CLI wiring against synthetic known-answer data, which is what regressions
actually show up in. Tests needing an optional dependency (Py-ART, wradlib,
GDAL's PDF driver) or the network are marked and skipped when unavailable.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

# headless plotting for any test that builds a figure
import matplotlib
matplotlib.use("Agg")


@pytest.fixture
def field_tif(tmp_path):
    """Factory: write a small EPSG:4326 float32 GeoTIFF with NaN nodata.

    Returns ``make(array, name=..., west=..., north=..., res=...) -> path``.
    The default origin/resolution mimics the MRMS CONUS grid near Reno, NV
    (0.01°, so ~1 km cells) so cell-size maths land in a realistic range.
    """
    import rasterio
    from rasterio.transform import from_origin

    def make(arr, name="field.tif", west=-119.75, north=39.60, res=0.01):
        arr = np.asarray(arr, dtype="float32")
        path = tmp_path / name
        with rasterio.open(
            path, "w", driver="GTiff", height=arr.shape[0], width=arr.shape[1],
            count=1, dtype="float32", crs="EPSG:4326",
            transform=from_origin(west, north, res, res), nodata=np.nan,
        ) as dst:
            dst.write(arr, 1)
        return str(path)

    return make


@pytest.fixture
def points_gdf():
    """Factory: ``make([(lon, lat), ...], **cols) -> GeoDataFrame`` in EPSG:4326."""
    import geopandas as gpd
    from shapely.geometry import Point

    def make(lonlats, **cols):
        return gpd.GeoDataFrame(
            dict(**cols), geometry=[Point(x, y) for x, y in lonlats], crs=4326)

    return make


@pytest.fixture
def gauge_obs():
    """Factory for one gauge's Synoptic-shaped observation frame.

    ``make(increments_mm, step_min=5, start=...)`` builds a frame with the
    ``date_time`` and ``precip_intervals_set_1d`` columns
    :mod:`stormscape.gauges` expects. Pass ``None`` inside ``increments_mm`` to
    simulate a multi-variable station whose rows carry no precip value.
    """
    from stormscape.gauges import INT_VAR

    def make(increments_mm, step_min=5, start="2026-06-19T20:00:00"):
        t0 = pd.Timestamp(start)
        times = [t0 + pd.Timedelta(minutes=step_min * i)
                 for i in range(len(increments_mm))]
        return pd.DataFrame({"date_time": times, INT_VAR: increments_mm})

    return make


@pytest.fixture
def minute_series():
    """Factory: a 1-minute rate series DataFrame as ``load_event_series`` returns.

    ``make(rates_mmph, start=...)`` -> DataFrame indexed by minute with a
    ``rate_mmph`` column (what :func:`stormscape.gauges.storm_window` reduces).
    """
    def make(rates_mmph, start="2026-06-19T20:00:00"):
        idx = pd.date_range(start, periods=len(rates_mmph), freq="1min")
        return pd.DataFrame({"rate_mmph": np.asarray(rates_mmph, dtype=float)},
                            index=idx)

    return make


def pytest_configure(config):
    config.addinivalue_line("markers", "network: needs internet (skipped by default)")
    config.addinivalue_line("markers", "optional_deps: needs an optional dependency")
