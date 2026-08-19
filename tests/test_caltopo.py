"""CalTopo GeoJSON: folder grouping, class colours, and what survives a round trip.

The behaviour under test is the one that failed in the field -- a network
exported as plain GeoJSON arrives as N independent objects with N checkboxes.
Every folder assertion here is really asserting "one toggle, not 1,400".
"""
from __future__ import annotations

import json

import pytest

from stormscape import caltopo

gpd = pytest.importorskip("geopandas")
from shapely.geometry import LineString, Point, Polygon  # noqa: E402


def lines(n=3, crs="EPSG:4326"):
    return gpd.GeoDataFrame(
        {"seg": list(range(n)), "P": [round(0.05 + 0.2 * i, 2) for i in range(n)]},
        geometry=[LineString([(-119.8 + 0.01 * i, 39.6),
                              (-119.79 + 0.01 * i, 39.61)]) for i in range(n)],
        crs=crs)


def folders_of(fc):
    return [f for f in fc["features"] if f["properties"]["class"] == "Folder"]


def shapes_of(fc):
    return [f for f in fc["features"] if f["properties"]["class"] != "Folder"]


def test_every_shape_lands_in_its_layer_folder():
    """The whole reason this module exists: one checkbox per layer."""
    fc = caltopo.build([caltopo.Layer("network", lines(4)),
                        caltopo.Layer("basins", lines(2))])
    fol = folders_of(fc)
    assert [f["properties"]["title"] for f in fol] == ["network", "basins"]
    ids = {f["id"] for f in fol}
    assert len(ids) == 2
    for s in shapes_of(fc):
        assert s["properties"]["folderId"] in ids
    by_folder = {}
    for s in shapes_of(fc):
        by_folder[s["properties"]["folderId"]] = \
            by_folder.get(s["properties"]["folderId"], 0) + 1
    assert sorted(by_folder.values()) == [2, 4]


def test_folders_precede_the_shapes_that_reference_them():
    fc = caltopo.build([caltopo.Layer("a", lines(2))])
    classes = [f["properties"]["class"] for f in fc["features"]]
    assert classes[0] == "Folder"
    assert set(classes[1:]) == {"Shape"}


def test_an_empty_layer_contributes_no_folder():
    """An empty folder in CalTopo is a checkbox that does nothing."""
    empty = lines(1).iloc[:0]
    fc = caltopo.build([caltopo.Layer("gone", empty),
                        caltopo.Layer("kept", lines(2))])
    assert [f["properties"]["title"] for f in folders_of(fc)] == ["kept"]


def test_points_become_markers_and_lines_become_shapes():
    pts = gpd.GeoDataFrame({"n": [1]}, geometry=[Point(-119.8, 39.6)],
                           crs="EPSG:4326")
    fc = caltopo.build([caltopo.Layer("m", pts, color="#112233"),
                        caltopo.Layer("l", lines(1))])
    kinds = {s["properties"]["class"]: s for s in shapes_of(fc)}
    assert set(kinds) == {"Marker", "Shape"}
    assert kinds["Marker"]["properties"]["marker-color"] == "#112233"
    assert "stroke" not in kinds["Marker"]["properties"]
    assert kinds["Shape"]["properties"]["stroke-width"] == 3.0


def test_label_column_becomes_title_not_name():
    """CalTopo labels from `title`; a `name` column imports unlabelled."""
    g = lines(2)
    g["seg_uid"] = ["16050102-7", "16050102-9"]
    fc = caltopo.build([caltopo.Layer("net", g, label="seg_uid")])
    assert sorted(s["properties"]["title"] for s in shapes_of(fc)) == \
        ["16050102-7", "16050102-9"]


def test_fields_travel_into_the_description():
    """The only route a modelled attribute takes into the field."""
    fc = caltopo.build([caltopo.Layer("net", lines(1), fields=("seg", "P"),
                                      description="observed storm")])
    d = shapes_of(fc)[0]["properties"]["description"]
    assert d.splitlines() == ["observed storm", "seg: 0", "P: 0.05"]


def test_description_drops_missing_and_blank_values():
    g = lines(2)
    g["note"] = ["", None]
    fc = caltopo.build([caltopo.Layer("n", g, fields=("note", "absent", "seg"))])
    for s in shapes_of(fc):
        assert s["properties"]["description"].startswith("seg: ")


def test_classify_assigns_colours_on_the_break_boundaries():
    cols = caltopo.classify([0.0, 0.2, 0.4, 0.6, 0.8, 0.95],
                            [0.2, 0.4, 0.6, 0.8])
    assert cols == [caltopo.CLASS_COLORS[0], caltopo.CLASS_COLORS[1],
                    caltopo.CLASS_COLORS[2], caltopo.CLASS_COLORS[3],
                    caltopo.CLASS_COLORS[4], caltopo.CLASS_COLORS[4]]


def test_classify_returns_none_for_nan_so_the_layer_colour_wins():
    assert caltopo.classify([float("nan"), None, "x"], [0.5],
                            ["#000000", "#FFFFFF"]) == [None, None, None]


def test_a_nan_class_paints_the_fallback_not_a_crash():
    g = lines(2)
    cols = caltopo.classify([float("nan"), 0.9], [0.5], ["#000000", "#FFFFFF"])
    fc = caltopo.build([caltopo.Layer("n", g, color=cols,
                                      fallback_color="#808080")])
    assert sorted(s["properties"]["stroke"] for s in shapes_of(fc)) == \
        ["#808080", "#FFFFFF"]


def test_classify_rejects_a_mismatched_colour_count():
    with pytest.raises(ValueError, match="breaks need"):
        caltopo.classify([0.1], [0.2, 0.4], ["#000000", "#FFFFFF"])


def test_class_labels_read_as_a_legend():
    assert caltopo.class_labels([0.2, 0.5], unit="") == \
        ["< 0.2", "0.2-0.5", "≥ 0.5"]


def test_a_wrong_length_colour_list_is_refused():
    with pytest.raises(ValueError, match="colours for"):
        caltopo.build([caltopo.Layer("n", lines(3), color=["#000000"])])


def test_non_wgs84_input_is_reprojected():
    """GeoJSON is WGS84 by spec and CalTopo assumes it."""
    g = lines(2).to_crs("EPSG:5070")
    fc = caltopo.build([caltopo.Layer("n", g)])
    xs = [c[0] for s in shapes_of(fc) for c in s["geometry"]["coordinates"]]
    assert all(-120.5 < x < -119.0 for x in xs), xs


def test_positions_match_caltopos_four_element_form_by_default():
    fc = caltopo.build([caltopo.Layer("n", lines(1))])
    for c in shapes_of(fc)[0]["geometry"]["coordinates"]:
        assert len(c) == 4 and c[2] == 0 and c[3] == 0


def test_strict_geojson_positions_on_request():
    fc = caltopo.build([caltopo.Layer("n", lines(1))], coord_len=2)
    assert all(len(c) == 2 for c in shapes_of(fc)[0]["geometry"]["coordinates"])


def test_polygon_rings_survive_rounding():
    poly = gpd.GeoDataFrame({"a": [1]}, crs="EPSG:4326", geometry=[Polygon(
        [(-119.8, 39.6), (-119.7, 39.6), (-119.7, 39.7), (-119.8, 39.6)])])
    fc = caltopo.build([caltopo.Layer("b", poly, fill="#3388FF")])
    s = shapes_of(fc)[0]
    assert s["geometry"]["type"] == "Polygon"
    assert len(s["geometry"]["coordinates"][0]) == 4
    assert s["properties"]["fill"] == "#3388FF"


def test_fill_is_omitted_unless_asked_for():
    """An unfilled basin outline lets the terrain show through it."""
    fc = caltopo.build([caltopo.Layer("b", lines(1))])
    assert "fill" not in shapes_of(fc)[0]["properties"]


def test_colours_are_normalized_to_caltopos_hex_form():
    fc = caltopo.build([caltopo.Layer("n", lines(1), color="d7191c")])
    assert shapes_of(fc)[0]["properties"]["stroke"] == "#D7191C"


def test_labels_default_off_because_a_network_is_unreadable_with_them_on():
    fc = caltopo.build([caltopo.Layer("n", lines(2))])
    assert folders_of(fc)[0]["properties"]["labelVisible"] is False


def test_written_file_is_valid_geojson_gdal_can_reopen(tmp_path):
    """The file also has to survive a QGIS sanity check, not just CalTopo.

    GDAL reads the 4-element positions, warns that it is dropping the tail past
    three, and keeps the third -- so a CalTopo-form file reopens as *3-D*
    geometry with z=0. Harmless for CalTopo, visible in QGIS; `coord_len=2`
    is the way out, and the bundle's GPKG copies never go through this path.
    """
    p = caltopo.write([caltopo.Layer("net", lines(3), fields=("P",))],
                      str(tmp_path / "bundle.geojson"))
    back = gpd.read_file(p)
    # The folder rides along as a null-geometry row -- GDAL keeps it, so a QGIS
    # check shows one extra empty row per folder rather than dropping them.
    assert len(back) == 4
    assert back.geometry.isna().sum() == 1
    shapes = back[back["class"] == "Shape"]
    assert len(shapes) == 3
    assert json.loads(open(p).read())["type"] == "FeatureCollection"
    assert shapes.geometry.iloc[0].coords[0] == pytest.approx((-119.8, 39.6, 0.0))
    assert shapes.geometry.iloc[0].has_z


def test_strict_positions_reopen_as_plain_2d(tmp_path):
    p = caltopo.write([caltopo.Layer("net", lines(2))],
                      str(tmp_path / "strict.geojson"), coord_len=2)
    back = gpd.read_file(p)
    shapes = back[back["class"] == "Shape"]
    assert len(shapes) == 2
    assert not shapes.geometry.iloc[0].has_z


def test_summary_counts_what_caltopo_will_show(tmp_path):
    p = caltopo.write([caltopo.Layer("net", lines(4)),
                       caltopo.Layer("basins", lines(2))],
                      str(tmp_path / "b.geojson"))
    s = caltopo.summary(p)
    assert s["folders"] == {"net": 4, "basins": 2}
    assert s["ungrouped"] == 0
    assert s["bytes"] > 0


def test_simplify_drops_vertices_but_keeps_the_line():
    dense = gpd.GeoDataFrame({"a": [1]}, crs="EPSG:4326", geometry=[LineString(
        [(-119.8 + 0.0001 * i, 39.6 + 0.00001 * (i % 2)) for i in range(60)])])
    plain = caltopo.build([caltopo.Layer("n", dense)])
    thin = caltopo.build([caltopo.Layer("n", dense, simplify_m=20.0)])
    n0 = len(shapes_of(plain)[0]["geometry"]["coordinates"])
    n1 = len(shapes_of(thin)[0]["geometry"]["coordinates"])
    assert n1 < n0 and n1 >= 2
