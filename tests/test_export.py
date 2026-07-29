"""Georeferenced export: EPSG:3857 reprojection, colorized RGBA transparency,
and the GeoPDF gate."""
from __future__ import annotations

import os

import numpy as np
import pytest

from stormscape import export


def has_module(name: str) -> bool:
    import importlib.util
    try:
        return importlib.util.find_spec(name) is not None
    except (ImportError, ValueError):
        return False


def test_reproject_geotiff_lands_in_web_mercator(field_tif, tmp_path):
    """CalTopo consumes EPSG:3857, so this is the whole point of the layer export."""
    import rasterio
    src = field_tif(np.random.RandomState(0).rand(12, 10) * 80)
    out = export.reproject_geotiff(src, str(tmp_path / "w.tif"))
    with rasterio.open(out) as ds:
        assert ds.crs.to_epsg() == 3857


def test_reproject_geotiff_honours_an_explicit_crs(field_tif, tmp_path):
    import rasterio
    src = field_tif(np.ones((8, 8)))
    out = export.reproject_geotiff(src, str(tmp_path / "utm.tif"),
                                   dst_crs="EPSG:32611")
    with rasterio.open(out) as ds:
        assert ds.crs.to_epsg() == 32611


def test_reproject_geotiff_roughly_preserves_the_value_range(field_tif, tmp_path):
    import rasterio
    arr = np.random.RandomState(1).rand(20, 20) * 100
    out = export.reproject_geotiff(field_tif(arr), str(tmp_path / "w.tif"))
    with rasterio.open(out) as ds:
        got = ds.read(1, masked=True)
    assert float(got.max()) <= arr.max() + 1e-3
    assert float(got.min()) >= arr.min() - 1e-3


def test_categorical_fields_are_reprojected_by_nearest_neighbour(field_tif, tmp_path):
    """Interpolating a quality flag or a time-of-peak code would invent values that
    mean nothing, so those fields must use nearest neighbour."""
    import rasterio
    arr = np.random.RandomState(2).randint(0, 4, (14, 14)).astype("float32")
    out = export.reproject_geotiff(field_tif(arr, name="event_rqi.tif"),
                                   str(tmp_path / "rqi_3857.tif"))
    with rasterio.open(out) as ds:
        vals = np.unique(ds.read(1, masked=True).compressed())
    assert set(np.round(vals, 6)).issubset({0.0, 1.0, 2.0, 3.0})


def test_is_categorical_recognises_the_flag_fields():
    for f in ("tpki15", "rqi", "shsr", "cbb"):
        assert export._is_categorical(f)
    for f in ("i15max", "i30max", "total", "anom_i15", "peakrate_mmph"):
        assert not export._is_categorical(f)


# --------------------------------------------------------------------------- #
# colorized RGBA — what actually gets dropped on a CalTopo basemap
# --------------------------------------------------------------------------- #
def test_colormap_rgba_makes_nan_fully_transparent():
    arr = np.array([[10.0, np.nan], [30.0, 40.0]])
    rgba = export._colormap_rgba(arr, "YlGnBu", 0.0, 50.0)
    assert rgba.shape == (2, 2, 4)
    assert rgba[0, 1, 3] == 0
    assert rgba[0, 0, 3] == 255


def test_colormap_rgba_masks_the_dry_tail():
    """Dry cells must read through to the basemap rather than painting it pale."""
    arr = np.array([[1.0, 8.0, 60.0]])
    rgba = export._colormap_rgba(arr, "YlGnBu", 0.0, 100.0, mask_below=5.0)
    assert rgba[0, 0, 3] == 0            # 1.0 < wet_min
    assert rgba[0, 1, 3] == 255
    assert rgba[0, 2, 3] == 255


def test_colormap_rgba_is_monotone_in_value():
    """Higher intensity must not map to a visually lower colour-ramp position."""
    arr = np.array([[0.0, 25.0, 50.0, 100.0]])
    rgba = export._colormap_rgba(arr, "YlGnBu", 0.0, 100.0)
    # YlGnBu darkens with value: luminance should fall
    lum = rgba[0, :, :3].astype(float).mean(axis=1)
    assert all(b <= a + 1e-9 for a, b in zip(lum, lum[1:]))


def test_colormap_rgba_handles_values_above_vmax_without_overflow():
    """TwoSlopeNorm maps anything above vmax to +inf. Feeding that to the colormap
    used to overflow (a RuntimeWarning, and only accidentally the right colour);
    out-of-range values must clip to the end colour, as imshow shows them."""
    import warnings
    from matplotlib.colors import TwoSlopeNorm
    arr = np.array([[0.5, 1.0, 2.0, 99.0]])          # 99 is far above vmax
    norm = TwoSlopeNorm(vcenter=1.0, vmin=0.0, vmax=3.0)
    with warnings.catch_warnings():
        warnings.simplefilter("error")               # any warning fails the test
        rgba = export._colormap_rgba(arr, "RdBu_r", 0.0, 3.0, norm=norm)
    top = export._colormap_rgba(np.array([[3.0]]), "RdBu_r", 0.0, 3.0, norm=norm)
    assert np.array_equal(rgba[0, 3, :3], top[0, 0, :3])
    assert rgba[0, 3, 3] == 255


def test_colormap_rgba_handles_values_below_vmin_without_overflow():
    import warnings
    from matplotlib.colors import TwoSlopeNorm
    arr = np.array([[-50.0, 1.0]])
    norm = TwoSlopeNorm(vcenter=1.0, vmin=0.0, vmax=3.0)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        rgba = export._colormap_rgba(arr, "RdBu_r", 0.0, 3.0, norm=norm)
    assert rgba[0, 0, 3] == 255                      # finite, so opaque


def test_field_style_centres_the_anomaly_on_one():
    """An anomaly is a ratio: the diverging colormap must be centred at 1x or the
    'above/below climatology' reading inverts."""
    arr = np.array([[0.5, 1.0, 3.0]])
    cmap, vmin, vmax, norm, mask_below = export._field_style("anom_i15", arr)
    assert norm is not None
    assert norm.vcenter == pytest.approx(1.0)
    assert mask_below is None            # anomalies are not dry-masked


def test_field_style_dry_masks_intensity_fields():
    arr = np.array([[0.0, 50.0]])
    _, _, _, norm, mask_below = export._field_style("i15max", arr, wet_min=5.0)
    assert norm is None
    assert mask_below == pytest.approx(5.0)


# --------------------------------------------------------------------------- #
# batch export
# --------------------------------------------------------------------------- #
def test_export_geotiffs_writes_raw_and_colorized_per_field(field_tif, tmp_path):
    d = tmp_path / "src"
    d.mkdir()
    for f in ("i15max", "anom_i15"):
        arr = np.random.RandomState(3).rand(10, 10) * (100 if f == "i15max" else 4)
        field_tif(arr, name=f"src/ev_{f}.tif")
    out = tmp_path / "out"
    written = export.export_geotiffs(str(d), "ev", str(out),
                                     fields=("i15max", "anom_i15"))
    names = {os.path.basename(p) for p in written}
    for f in ("i15max", "anom_i15"):
        assert f"ev_{f}_3857.tif" in names
        assert f"ev_{f}_3857_rgb.tif" in names


def test_export_geotiffs_rgba_has_a_tagged_alpha_band(field_tif, tmp_path):
    import rasterio
    from rasterio.enums import ColorInterp
    d = tmp_path / "src"
    d.mkdir()
    field_tif(np.random.RandomState(4).rand(10, 10) * 100, name="src/ev_i15max.tif")
    export.export_geotiffs(str(d), "ev", str(tmp_path / "o"), fields=("i15max",))
    with rasterio.open(tmp_path / "o" / "ev_i15max_3857_rgb.tif") as ds:
        assert ds.count == 4
        assert ds.colorinterp[3] == ColorInterp.alpha


def test_export_geotiffs_skips_missing_fields_without_failing(field_tif, tmp_path):
    d = tmp_path / "src"
    d.mkdir()
    field_tif(np.ones((6, 6)), name="src/ev_i15max.tif")
    written = export.export_geotiffs(str(d), "ev", str(tmp_path / "o"),
                                     fields=("i15max", "does_not_exist"))
    assert any("i15max" in p for p in written)
    assert not any("does_not_exist" in p for p in written)


def test_export_geotiffs_can_skip_the_raw_layer(field_tif, tmp_path):
    d = tmp_path / "src"
    d.mkdir()
    field_tif(np.ones((6, 6)) * 20, name="src/ev_i15max.tif")
    written = export.export_geotiffs(str(d), "ev", str(tmp_path / "o"),
                                     fields=("i15max",), raw=False)
    assert all(p.endswith("_rgb.tif") for p in written)
    assert not (tmp_path / "o" / "ev_i15max_3857.tif").exists()


def test_vector_driver_is_inferred_from_the_extension():
    assert export._vector_driver("x.geojson") == "GeoJSON"
    assert export._vector_driver("x.gpkg") == "GPKG"
    assert export._vector_driver("x.shp") == "ESRI Shapefile"
    assert export._vector_driver("x.kml") == "KML"
    assert export._vector_driver("x.unknown") == "GeoJSON"      # safe default


# --------------------------------------------------------------------------- #
# GeoPDF — gated on GDAL's PDF driver
# --------------------------------------------------------------------------- #
def test_geopdf_supported_returns_a_bool():
    assert isinstance(export.geopdf_supported(), bool)


@pytest.mark.optional_deps
@pytest.mark.skipif(not has_module("osgeo") or not export.geopdf_supported(),
                    reason="needs GDAL's PDF driver (conda install -c conda-forge libgdal-pdf)")
def test_figure_to_geopdf_georeferences_the_map_axes(tmp_path):
    """The neatline must bound the *map frame* in the axes CRS, so the colorbar and
    title margins are excluded rather than mis-georeferenced."""
    import matplotlib.pyplot as plt
    from osgeo import gdal, ogr, osr
    gdal.UseExceptions()

    x0, x1, y0, y1 = -13330509.0, -13317151.0, 4786337.0, 4802207.0
    fig, ax = plt.subplots(figsize=(5, 6))
    im = ax.imshow(np.random.RandomState(0).rand(20, 15),
                   extent=(x0, x1, y0, y1), origin="upper")
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    fig.colorbar(im, ax=ax)                      # chrome the neatline must exclude
    out = export.figure_to_geopdf(fig, ax, str(tmp_path / "m.pdf"),
                                  crs="EPSG:3857", dpi=100)
    plt.close(fig)

    ds = gdal.Open(out)
    assert osr.SpatialReference(ds.GetProjection()).GetAuthorityCode(None) == "3857"
    ring = ogr.CreateGeometryFromWkt(ds.GetMetadataItem("NEATLINE")).GetGeometryRef(0)
    xs = [ring.GetPoint(i)[0] for i in range(ring.GetPointCount())]
    ys = [ring.GetPoint(i)[1] for i in range(ring.GetPointCount())]
    assert min(xs) == pytest.approx(x0, abs=1.0)
    assert max(xs) == pytest.approx(x1, abs=1.0)
    assert min(ys) == pytest.approx(y0, abs=1.0)
    assert max(ys) == pytest.approx(y1, abs=1.0)
