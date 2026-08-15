"""Abandoned-mine features: grouping, the ArcGIS paging contract, and density.

Entirely offline -- the USMIN service is replaced by a recording double that
replays canned GeoJSON, so no test here touches the network.

The invariants worth pinning are the ones that failed silently during
development: a hosted FeatureServer hides its "there is more" flag somewhere a
MapServer does not, an ArcGIS error arrives as HTTP 200 and reads as an empty
AOI, and mine waste lives almost entirely in the polygon layer.
"""
from __future__ import annotations

import json

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import Polygon

from stormscape import mines, refdata


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _feature(lon, lat, **props):
    return {"type": "Feature",
            "properties": props,
            "geometry": {"type": "Point", "coordinates": [lon, lat]}}


def _poly_feature(lon, lat, **props):
    d = 0.001
    ring = [[lon, lat], [lon + d, lat], [lon + d, lat + d], [lon, lat + d],
            [lon, lat]]
    return {"type": "Feature", "properties": props,
            "geometry": {"type": "Polygon", "coordinates": [ring]}}


class _FakeResponse:
    def __init__(self, payload):
        self._payload = payload

    def json(self):
        return self._payload


def _install(monkeypatch, pages, record=None):
    """Replace requests.get with a double replaying ``pages`` in order."""
    seq = list(pages)

    def fake_get(url, params=None, timeout=None):
        if record is not None:
            record.append(dict(url=url, params=dict(params or {})))
        return _FakeResponse(seq.pop(0) if seq else {"type": "FeatureCollection",
                                                     "features": []})

    monkeypatch.setattr(refdata.requests, "get", fake_get)


# --------------------------------------------------------------------------- #
# feature-type grouping
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("ftr_type,group", [
    ("Mine Dump", "waste"),
    ("Tailings - Undifferentiated", "waste"),
    ("Slag Pile", "waste"),
    ("Adit", "openings"),
    ("Mine Shaft", "openings"),
    ("Open Pit Mine", "surface"),
    ("Quarry - Limestone", "surface"),
    ("Gravel Pit", "aggregate"),
    ("Prospect Pit", "prospect"),
    ("Mill Site", "other"),
])
def test_group_of_covers_the_usmin_vocabulary(ftr_type, group):
    assert mines.group_of(ftr_type) == group


def test_group_of_is_case_and_space_insensitive():
    assert mines.group_of("  mine DUMP ") == "waste"


def test_group_of_sends_unknown_types_to_other_not_an_error():
    """A type USMIN adds later must land somewhere visible, not crash."""
    assert mines.group_of("Antimatter Pit") == "other"
    assert mines.group_of(None) == "other"
    assert mines.group_of("") == "other"


def test_prefix_families_absorb_future_subtypes():
    """'Tailings - <anything>' is waste without editing the table -- the whole
    point of the prefix rules."""
    assert mines.group_of("Tailings - Something New") == "waste"
    assert mines.group_of("Quarry - Unobtainium") == "surface"


def test_every_known_type_maps_into_a_declared_group():
    for t, g in mines._EXACT.items():
        assert g in mines.GROUPS, f"{t} -> undeclared group {g}"


def test_default_kinds_are_real_groups():
    assert set(mines.DEFAULT_KINDS) <= set(mines.GROUPS)


def test_cli_default_kinds_track_the_module():
    """The parser hard-codes the default for help text; catch the drift here."""
    from stormscape import cli
    assert tuple(cli._MINE_DEFAULT_KINDS) == tuple(mines.DEFAULT_KINDS)


# --------------------------------------------------------------------------- #
# kind resolution and the server-side where clause
# --------------------------------------------------------------------------- #
def test_resolve_kinds_accepts_groups_types_and_strings():
    g, t = mines._resolve_kinds(["waste", "Mine Shaft"])
    assert g == {"waste"} and t == {"mine shaft"}
    g, t = mines._resolve_kinds("waste,openings")
    assert g == {"waste", "openings"} and not t


def test_resolve_kinds_all_means_everything():
    for spec in (None, "all", ["all"], []):
        g, t = mines._resolve_kinds(spec)
        assert g == set(mines.GROUPS) and not t


def test_where_clause_filters_server_side_and_keeps_prefix_families():
    where = mines._where_for({"waste"}, set(), "ftr_type")
    assert "'mine dump'" in where
    assert "LIKE 'tailings%'" in where       # future subtypes without a round trip


def test_where_clause_is_none_when_it_cannot_be_expressed():
    """'other' is a catch-all for *unrecognised* types, so it cannot be
    enumerated in SQL -- the caller must filter locally instead."""
    assert mines._where_for({"other"}, set(), "ftr_type") is None
    assert mines._where_for(set(mines.GROUPS), set(), "ftr_type") is None


def test_where_clause_escapes_quotes():
    w = mines._where_for(set(), {"o'brien pit"}, "ftr_type")
    assert "o''brien pit" in w


# --------------------------------------------------------------------------- #
# the paging contract -- the bug that silently truncated USMIN
# --------------------------------------------------------------------------- #
def test_more_pages_reads_the_mapserver_top_level_flag():
    assert refdata._more_pages({"exceededTransferLimit": True})


def test_more_pages_reads_the_featureserver_nested_flag():
    """A hosted FeatureServer puts it ONLY under 'properties'. Checking the top
    level alone stops after one page -- USMIN returned exactly 2000 of 3134
    features over the Hidden Valley AOI, with no error to notice."""
    assert refdata._more_pages({"properties": {"exceededTransferLimit": True}})


def test_more_pages_false_when_complete():
    assert not refdata._more_pages({"features": [], "properties": {}})
    assert not refdata._more_pages({})


def test_query_pages_through_a_nested_transfer_limit(monkeypatch):
    page1 = {"type": "FeatureCollection",
             "features": [_feature(-119.5, 39.5, ftr_type="Adit")],
             "properties": {"exceededTransferLimit": True}}
    page2 = {"type": "FeatureCollection",
             "features": [_feature(-119.6, 39.6, ftr_type="Adit")],
             "properties": {}}
    calls = []
    _install(monkeypatch, [page1, page2], record=calls)
    g = refdata.arcgis_query("svc", 17, (-120, 39, -119, 40), page=1)
    assert len(g) == 2, "nested flag ignored -> silent truncation"
    assert [c["params"]["resultOffset"] for c in calls] == [0, 1]


def test_query_stops_when_a_page_is_short(monkeypatch):
    _install(monkeypatch, [{"type": "FeatureCollection", "features": []}])
    assert not len(refdata.arcgis_query("svc", 17, (-120, 39, -119, 40)))


# --------------------------------------------------------------------------- #
# a service error must not read as an empty AOI
# --------------------------------------------------------------------------- #
def test_service_error_warns_rather_than_returning_a_quiet_empty(monkeypatch,
                                                                capsys):
    """ArcGIS returns errors as HTTP 200 with an {"error": ...} body. Left
    unchecked that is indistinguishable from an AOI with no features -- which
    is how the NBMG mirror's 'Pagination is not supported' went unnoticed."""
    _install(monkeypatch, [{"error": {"code": 400,
                                      "message": "Pagination is not supported."}}])
    g = refdata.arcgis_query("svc", 0, (-120, 39, -119, 40), what="usmin_nbmg")
    assert not len(g)
    err = capsys.readouterr().err
    assert "usmin_nbmg" in err and "Pagination is not supported" in err


def test_paginates_false_omits_the_offset_params(monkeypatch):
    """A server that rejects resultOffset outright must not be sent one."""
    calls = []
    _install(monkeypatch, [{"type": "FeatureCollection", "features": []}],
             record=calls)
    refdata.arcgis_query("svc", 0, (-120, 39, -119, 40), paginates=False)
    assert "resultOffset" not in calls[0]["params"]
    assert "resultRecordCount" not in calls[0]["params"]


def test_token_is_passed_through_when_given(monkeypatch):
    calls = []
    _install(monkeypatch, [{"type": "FeatureCollection", "features": []}],
             record=calls)
    refdata.arcgis_query("svc", 0, (-120, 39, -119, 40), token="SEKRIT")
    assert calls[0]["params"]["token"] == "SEKRIT"


def test_no_token_key_when_absent(monkeypatch):
    calls = []
    _install(monkeypatch, [{"type": "FeatureCollection", "features": []}],
             record=calls)
    refdata.arcgis_query("svc", 0, (-120, 39, -119, 40))
    assert "token" not in calls[0]["params"]


# --------------------------------------------------------------------------- #
# mine_features
# --------------------------------------------------------------------------- #
def test_mine_features_normalises_columns_and_groups(monkeypatch):
    page = {"type": "FeatureCollection", "features": [
        _feature(-119.5, 39.5, ftr_type="Adit", ftr_name="Devils Gate",
                 state="NV", county="Storey", topo_name="Virginia City",
                 topo_date="1950"),
        _feature(-119.6, 39.6, ftr_type="Mine Shaft", ftr_name="",
                 state="NV", county="Storey", topo_name="Virginia City",
                 topo_date="1950"),
    ], "properties": {}}
    _install(monkeypatch, [page])
    g = mines.mine_features((-120, 39, -119, 40), kinds="openings",
                            geometry="points")
    assert list(g.columns)[-1] == "geometry"
    assert set(g.group) == {"openings"}
    assert g.source.unique().tolist() == ["usmin"]
    assert g.geom_kind.unique().tolist() == ["point"]


def test_blank_names_become_null(monkeypatch):
    """USMIN writes unnamed features as '', not null -- 2990 of 3134 over the
    Hidden Valley AOI. Left as '', .notna() calls every feature named."""
    page = {"type": "FeatureCollection", "features": [
        _feature(-119.5, 39.5, ftr_type="Adit", ftr_name=""),
        _feature(-119.6, 39.6, ftr_type="Adit", ftr_name="   "),
        _feature(-119.7, 39.7, ftr_type="Adit", ftr_name="Devils Gate"),
    ], "properties": {}}
    _install(monkeypatch, [page])
    g = mines.mine_features((-120, 39, -119, 40), kinds="openings",
                            geometry="points")
    assert g.name.notna().sum() == 1


def test_named_only_keeps_just_the_named(monkeypatch):
    page = {"type": "FeatureCollection", "features": [
        _feature(-119.5, 39.5, ftr_type="Adit", ftr_name=""),
        _feature(-119.7, 39.7, ftr_type="Adit", ftr_name="Devils Gate"),
    ], "properties": {}}
    _install(monkeypatch, [page])
    g = mines.mine_features((-120, 39, -119, 40), kinds="openings",
                            geometry="points", named_only=True)
    assert len(g) == 1 and g.name.iloc[0] == "Devils Gate"


def test_kinds_are_reapplied_locally(monkeypatch):
    """The server-side where clause is an optimisation, not the contract: a
    source that ignores it must still honour `kinds`."""
    page = {"type": "FeatureCollection", "features": [
        _feature(-119.5, 39.5, ftr_type="Adit"),
        _feature(-119.6, 39.6, ftr_type="Prospect Pit"),   # not requested
    ], "properties": {}}
    _install(monkeypatch, [page])
    g = mines.mine_features((-120, 39, -119, 40), kinds="openings",
                            geometry="points")
    assert set(g.ftr_type) == {"Adit"}


def test_geometry_both_reads_the_polygon_layer_too(monkeypatch):
    """Mine waste is 14,815 polygons against 413 points nationally, so a
    points-only default would make a waste query look like an empty AOI."""
    pts = {"type": "FeatureCollection",
           "features": [_feature(-119.5, 39.5, ftr_type="Mine Dump")],
           "properties": {}}
    polys = {"type": "FeatureCollection",
             "features": [_poly_feature(-119.6, 39.6, ftr_type="Mine Dump")],
             "properties": {}}
    _install(monkeypatch, [pts, polys])
    g = mines.mine_features((-120, 39, -119, 40), kinds="waste", geometry="both")
    assert sorted(g.geom_kind) == ["area", "point"]


def test_geometry_default_is_both():
    import inspect
    sig = inspect.signature(mines.mine_features)
    assert sig.parameters["geometry"].default == "both"


def test_unknown_source_is_rejected():
    with pytest.raises(ValueError, match="unknown mine source"):
        mines.mine_features((-120, 39, -119, 40), source="nope")


def test_bad_geometry_is_rejected():
    with pytest.raises(ValueError, match="geometry must be"):
        mines.mine_features((-120, 39, -119, 40), geometry="sideways")


def test_gated_source_without_a_token_raises_a_useful_error(monkeypatch):
    monkeypatch.delenv("STORMSCAPE_NDOM_TOKEN", raising=False)
    with pytest.raises(PermissionError, match="needs a token"):
        mines.mine_features((-120, 39, -119, 40), source="ndom")


def test_gated_source_accepts_an_env_token(monkeypatch):
    monkeypatch.setenv("STORMSCAPE_NDOM_TOKEN", "abc")
    _install(monkeypatch, [{"type": "FeatureCollection", "features": []}])
    assert not len(mines.mine_features((-120, 39, -119, 40), source="ndom"))


def test_empty_result_has_the_full_schema(monkeypatch):
    """Downstream code indexes these columns; an empty frame must still carry
    them or a quiet AOI becomes a KeyError."""
    _install(monkeypatch, [{"type": "FeatureCollection", "features": []}])
    g = mines.mine_features((-120, 39, -119, 40))
    assert not len(g)
    for col in ("name", "ftr_type", "group", "geom_kind", "source", "geometry"):
        assert col in g.columns
    assert g.crs is not None


def test_register_source_round_trips():
    try:
        mines.register_source("scratch", label="x", service="s",
                              layers={"points": 0}, fields=mines._USMIN_FIELDS,
                              public=True)
        assert "scratch" in mines.SOURCES
    finally:
        mines.SOURCES.pop("scratch", None)


# --------------------------------------------------------------------------- #
# density
# --------------------------------------------------------------------------- #
def _grid_gdf(n_by_cell):
    """Points clustered so each 1 km cell holds a known count."""
    rows = []
    for i, n in enumerate(n_by_cell):
        for _ in range(n):
            rows.append(_feature(-119.5 + i * 0.05, 39.5, ftr_type="Prospect Pit"))
    g = gpd.GeoDataFrame.from_features(rows, crs=4326)
    g["group"] = "prospect"
    return g


def test_density_grid_counts_and_per_km2():
    g = _grid_gdf([5, 3])
    d = mines.density_grid(g, cell_km=1.0)
    assert d["count"].sum() == 8
    assert set(d["count"]) == {5, 3}
    np.testing.assert_allclose(d["per_km2"], d["count"])   # 1 km cell


def test_density_grid_scales_per_km2_with_cell_size():
    g = _grid_gdf([4])
    d = mines.density_grid(g, cell_km=2.0)
    assert d["count"].sum() == 4
    np.testing.assert_allclose(d["per_km2"].sum(), 1.0)    # 4 / (2 km)^2


def test_density_grid_sorted_by_count_descending():
    d = mines.density_grid(_grid_gdf([2, 9, 5]), cell_km=1.0)
    assert list(d["count"]) == sorted(d["count"], reverse=True)


def test_density_grid_filters_by_group():
    g = _grid_gdf([4])
    g.loc[g.index[:2], "group"] = "waste"
    assert mines.density_grid(g, groups="waste")["count"].sum() == 2
    assert mines.density_grid(g, groups="prospect")["count"].sum() == 2


def test_density_grid_empty_input_returns_typed_empty():
    d = mines.density_grid(_grid_gdf([1]).iloc[0:0])
    assert not len(d) and "count" in d.columns and d.crs is not None


def test_density_grid_rejects_a_nonpositive_cell():
    with pytest.raises(ValueError, match="cell_km"):
        mines.density_grid(_grid_gdf([1]), cell_km=0)


def test_density_uses_polygon_representative_points():
    """A polygon has no .x/.y, and a crescent dump's centroid can fall outside
    it -- representative_point is guaranteed inside."""
    g = gpd.GeoDataFrame(
        {"group": ["waste"]},
        geometry=[Polygon([(-119.5, 39.5), (-119.4, 39.5), (-119.4, 39.6),
                           (-119.5, 39.6)])], crs=4326)
    assert mines.density_grid(g)["count"].sum() == 1


def test_density_raster_totals_match_the_vector_grid():
    """The two density paths must not disagree -- one drives the figure, the
    other the GeoTIFF."""
    g = _grid_gdf([5, 3])
    r = mines.density_raster(g, (-119.6, 39.45, -119.4, 39.55), cell_km=1.0)
    arr = r["fields"]["mine_density"]
    vec = mines.density_grid(g, cell_km=1.0)
    assert arr.sum() == pytest.approx(vec["count"].sum())
    assert arr.max() == pytest.approx(vec["count"].max())


def test_density_raster_profile_is_equal_area_and_writable():
    r = mines.density_raster(_grid_gdf([2]), (-119.6, 39.45, -119.4, 39.55))
    prof = r["profile"]
    assert prof["crs"].to_epsg() == 5070          # equal-area, not 4326/3857
    assert prof["dtype"] == "float32" and prof["count"] == 1
    assert prof["height"] > 0 and prof["width"] > 0


def test_density_raster_survives_an_empty_frame():
    r = mines.density_raster(_grid_gdf([1]).iloc[0:0],
                             (-119.6, 39.45, -119.4, 39.55))
    assert r["fields"]["mine_density"].max() == 0
    assert r["meta"]["n_features"] == 0


def test_density_raster_saves_through_save_fields(tmp_path):
    """It claims to be an mrms-style result dict; hold it to that."""
    from stormscape.mrms import save_fields
    r = mines.density_raster(_grid_gdf([3]), (-119.6, 39.45, -119.4, 39.55))
    paths = save_fields(r, str(tmp_path), "ev", layout="flat")
    assert len(paths) == 1 and paths[0].endswith("ev_mine_density.tif")


# --------------------------------------------------------------------------- #
# tables
# --------------------------------------------------------------------------- #
def test_group_counts_lists_types_present():
    g = _grid_gdf([2])
    t = mines.group_counts(g)
    assert t.n.sum() == 2 and "Prospect Pit" in t.types.iloc[0]


def test_group_counts_empty_frame_has_columns():
    t = mines.group_counts(_grid_gdf([1]).iloc[0:0])
    assert list(t.columns) == ["group", "description", "n", "types"]


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _run_mines_cli(monkeypatch, argv):
    from stormscape import cli
    seen = {}
    monkeypatch.setattr(cli, "_cmd_mines", lambda args: seen.update(vars(args)))
    cli.main(argv)
    return seen


def test_mines_cli_defaults_suit_counts_not_rainfall(monkeypatch):
    got = _run_mines_cli(monkeypatch,
                         ["mines", "--bbox", "-120", "39", "-119", "40"])
    assert got["cmap"] == "YlOrBr"            # not the rainfall YlGnBu
    assert got["wet_min"] == 0.5              # not 5 mm/h
    assert got["geometry"] == "both"          # or the waste goes missing
    assert got["kinds"] == list(mines.DEFAULT_KINDS)
    assert got["source"] == "usmin"


def test_mines_cli_accepts_kinds_and_density_knobs(monkeypatch):
    got = _run_mines_cli(monkeypatch,
                         ["mines", "--bbox", "-120", "39", "-119", "40",
                          "--kinds", "waste", "Mine Dump", "--cell-km", "2",
                          "--density-group", "waste", "--mines-mode", "density"])
    assert got["kinds"] == ["waste", "Mine Dump"]
    assert got["cell_km"] == 2.0
    assert got["density_group"] == "waste"
    assert got["mines_mode"] == "density"


@pytest.mark.parametrize("cmd", ["map", "run", "nexrad", "zoom", "burn",
                                 "export"])
def test_map_commands_expose_the_mine_overlay(cmd, monkeypatch):
    """--mines has to reach every command that renders through drape_i15, or
    the flag silently does nothing on some of them."""
    from stormscape import cli
    seen = {}
    monkeypatch.setattr(cli, f"_cmd_{cmd}", lambda args: seen.update(vars(args)))
    base = {"map": ["--hillshade", "h.tif", "--i15", "i.tif"],
            "run": ["--bbox", "-120", "39", "-119", "40", "--date", "20260619"],
            "nexrad": ["--bbox", "-120", "39", "-119", "40"],
            "zoom": ["--from-dir", ".", "--from-key", "k",
                     "--bbox", "-120", "39", "-119", "40"],
            "burn": ["--bbox", "-120", "39", "-119", "40"],
            "export": ["--from-dir", ".", "--from-key", "k"]}[cmd]
    cli.main([cmd] + base + ["--mines", "--mines-mode", "density"])
    assert seen["mines"] is True and seen["mines_mode"] == "density"


def test_cmd_mines_only_reads_args_the_parser_defines():
    """Every ``args.<name>`` the command touches must exist in the namespace.

    Regression: ``_cmd_mines`` read ``args.vmax``, which only the `burn` parser
    defines -- so the command fetched, wrote three files and then died with
    AttributeError at the very last step. Patching the command out (as the CLI
    tests above do) cannot catch that, and the real path needs the network, so
    check the two statically instead.
    """
    import ast
    import inspect

    from stormscape import cli

    tree = ast.parse(inspect.getsource(cli._cmd_mines))
    used = {n.attr for n in ast.walk(tree)
            if isinstance(n, ast.Attribute) and isinstance(n.value, ast.Name)
            and n.value.id == "args"}
    seen = {}
    monkey = lambda a: seen.update(vars(a))            # noqa: E731
    real, cli._cmd_mines = cli._cmd_mines, monkey
    try:
        cli.main(["mines", "--bbox", "-120", "39", "-119", "40"])
    finally:
        cli._cmd_mines = real
    missing = sorted(used - set(seen))
    assert not missing, f"_cmd_mines reads undefined args: {missing}"


def test_mine_kwargs_is_empty_when_the_flag_is_off():
    from argparse import Namespace

    from stormscape import cli
    assert cli._mine_kwargs(Namespace(mines=False)) == {}
    assert cli._mine_kwargs(Namespace()) == {}       # command without the flags


def test_mine_kwargs_passes_the_knobs_through():
    from argparse import Namespace

    from stormscape import cli
    kw = cli._mine_kwargs(Namespace(mines=True, mines_kinds=["waste"],
                                    mines_mode="density", mines_cell_km=2.0,
                                    mines_group="waste", mines_labels=True))
    assert kw["mines"] is True and kw["mines_cell_km"] == 2.0
    assert kw["mines_groups"] == "waste" and kw["mines_kinds"] == ["waste"]


# --------------------------------------------------------------------------- #
# plotting styles
# --------------------------------------------------------------------------- #
def test_every_group_has_a_plot_style():
    from stormscape import plot
    for g in mines.GROUPS:
        assert g in plot.MINE_STYLE, f"no style for group {g}"


def test_add_mines_on_an_empty_frame_is_a_no_op():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from stormscape import plot
    fig, ax = plt.subplots()
    empty = _grid_gdf([1]).iloc[0:0]
    assert plot.add_mines(ax, "EPSG:4326", empty) == []
    plt.close(fig)


def test_add_mines_points_returns_one_handle_per_group():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from stormscape import plot
    g = _grid_gdf([3])
    fig, ax = plt.subplots()
    handles = plot.add_mines(ax, "EPSG:4326", g, mode="points")
    assert len(handles) == 1 and "prospect" in handles[0].get_label()
    plt.close(fig)


def test_add_mines_density_draws_a_size_legend():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from stormscape import plot
    fig, ax = plt.subplots()
    handles = plot.add_mines(ax, "EPSG:4326", _grid_gdf([6, 2]), mode="density")
    assert handles, "a graduated symbol is unreadable without a key"
    assert any("per" in h.get_label() for h in handles)
    plt.close(fig)


def test_drape_renders_terrain_only_with_no_field(field_tif):
    """A map whose subject is the vectors must not also drape a raster of the
    same counts -- drawn both ways at once the reader sees one quantity twice.
    So `i15=None` is a supported mode: terrain + overlays, and no colour bar."""
    from stormscape.plot import drape_i15

    hs = field_tif(np.linspace(0, 255, 36).reshape(6, 6), name="hs.tif")
    fig, ax = drape_i15(hs, None, work_crs="EPSG:4326",
                        mines=_grid_gdf([2]), mines_mode="points")
    # the field imshow and its colour bar are both absent; the hillshade remains
    assert len(fig.axes) == 1, "a colour bar was drawn for a non-existent field"
    assert len(ax.images) == 1
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_drape_still_draws_the_colourbar_when_a_field_is_given(field_tif):
    from stormscape.plot import drape_i15

    hs = field_tif(np.linspace(0, 255, 36).reshape(6, 6), name="hs.tif")
    fld = field_tif(np.full((6, 6), 20.0), name="f.tif")
    fig, ax = drape_i15(hs, fld, work_crs="EPSG:4326")
    assert len(fig.axes) == 2                     # map + colour bar
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_drape_needs_a_hillshade_or_a_field():
    from stormscape.plot import drape_i15
    with pytest.raises(ValueError, match="hillshade, a field, or both"):
        drape_i15(None, None)


def test_density_map_is_off_by_default(monkeypatch):
    got = _run_mines_cli(monkeypatch,
                         ["mines", "--bbox", "-120", "39", "-119", "40"])
    assert got["density_map"] is False


def test_add_mines_auto_switches_to_density_when_crowded(monkeypatch):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from stormscape import plot
    monkeypatch.setattr(plot, "MINE_DENSITY_SWITCH", 2)
    fig, ax = plt.subplots()
    handles = plot.add_mines(ax, "EPSG:4326", _grid_gdf([5]), mode="auto")
    assert any("per" in h.get_label() for h in handles)   # density, not points
    plt.close(fig)
