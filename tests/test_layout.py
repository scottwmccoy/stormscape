"""Product layout: sorting on write, and reading BOTH layouts.

The backward-compatibility rule is the important one — every event folder
written before the sorted layout existed is flat, and those must keep working
through `--from-dir` / `--radar-dir` without migration.
"""
from __future__ import annotations

import os

import pytest

from stormscape import layout


# --------------------------------------------------------------------------- #
# which subdirectory a product belongs in
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("name,sub", [
    ("ev_i15max.tif", "rasters"),
    ("ev_dem.tiff", "rasters"),
    ("ev.png", "figures"),
    ("ev_anom_i15.pdf", "figures"),
    ("ev_compare.csv", "tables"),
    ("ev_gauge_recurrence.md", "tables"),
    ("ev_gauges.geojson", "vectors"),
    ("ev_streams.gpkg", "vectors"),
])
def test_products_sort_by_extension(name, sub):
    assert layout.subdir_for(name) == sub


def test_unknown_extensions_stay_at_the_event_root():
    assert layout.subdir_for("ev_pick.html") is None
    assert layout.subdir_for("notes.rst") is None


def test_readme_stays_at_the_event_root():
    """README.md describes the whole folder — filing it under tables/ with the
    CSVs would bury it."""
    assert layout.subdir_for("README.md") is None
    assert layout.subdir_for("ev_gauge_recurrence.md") == "tables"


# --------------------------------------------------------------------------- #
# writing
# --------------------------------------------------------------------------- #
def test_out_path_sorts_and_creates_the_subdirectory(tmp_path):
    p = layout.out_path(str(tmp_path), "ev_i15max.tif")
    assert p == str(tmp_path / "rasters" / "ev_i15max.tif")
    assert os.path.isdir(tmp_path / "rasters")


def test_flat_layout_writes_straight_to_the_event_root(tmp_path):
    p = layout.out_path(str(tmp_path), "ev_i15max.tif", layout="flat")
    assert p == str(tmp_path / "ev_i15max.tif")
    assert not os.path.exists(tmp_path / "rasters")


def test_env_var_switches_the_default(tmp_path, monkeypatch):
    monkeypatch.setenv("STORMSCAPE_LAYOUT", "flat")
    assert layout.out_path(str(tmp_path), "a.tif") == str(tmp_path / "a.tif")
    monkeypatch.setenv("STORMSCAPE_LAYOUT", "sorted")
    assert layout.out_path(str(tmp_path), "a.tif") == str(tmp_path / "rasters" / "a.tif")


def test_explicit_layout_beats_the_env_var(tmp_path, monkeypatch):
    monkeypatch.setenv("STORMSCAPE_LAYOUT", "flat")
    p = layout.out_path(str(tmp_path), "a.tif", layout="sorted")
    assert p == str(tmp_path / "rasters" / "a.tif")


# --------------------------------------------------------------------------- #
# reading — the backward-compatibility contract
# --------------------------------------------------------------------------- #
def test_find_reads_a_legacy_flat_event_folder(tmp_path):
    """Folders written before the sorted layout must keep working untouched."""
    (tmp_path / "ev_i15max.tif").write_text("x")
    assert layout.find(str(tmp_path), "ev_i15max.tif") == str(tmp_path / "ev_i15max.tif")


def test_find_reads_a_sorted_event_folder(tmp_path):
    (tmp_path / "rasters").mkdir()
    (tmp_path / "rasters" / "ev_i15max.tif").write_text("x")
    assert layout.find(str(tmp_path), "ev_i15max.tif") == str(
        tmp_path / "rasters" / "ev_i15max.tif")


def test_find_prefers_the_sorted_copy_when_both_exist(tmp_path):
    (tmp_path / "rasters").mkdir()
    (tmp_path / "ev.tif").write_text("old")
    (tmp_path / "rasters" / "ev.tif").write_text("new")
    assert layout.find(str(tmp_path), "ev.tif") == str(tmp_path / "rasters" / "ev.tif")


def test_find_points_at_the_sorted_path_when_the_file_is_absent(tmp_path):
    """A missing-file error should name where a fresh run would have written it."""
    assert layout.find(str(tmp_path), "ev.tif") == str(tmp_path / "rasters" / "ev.tif")


def test_find_accepts_a_full_path_and_uses_only_the_basename(tmp_path):
    (tmp_path / "ev.tif").write_text("x")
    assert layout.find(str(tmp_path), "/somewhere/else/ev.tif") == str(tmp_path / "ev.tif")


# --------------------------------------------------------------------------- #
# product directories
# --------------------------------------------------------------------------- #
def test_figure_directories_nest_under_figures(tmp_path):
    d = layout.subdir(str(tmp_path), "VirtualGaugeFigures")
    assert d == str(tmp_path / "figures" / "VirtualGaugeFigures")


def test_stores_and_caches_stay_at_the_event_root(tmp_path):
    """RainGaugeData is a store and the caches are inputs — not sorted products,
    and nexrad_cache in particular should stay obvious to delete."""
    for name in ("RainGaugeData", "nexrad_cache", "atlas14_cache"):
        assert layout.subdir(str(tmp_path), name) == str(tmp_path / name)


def test_flat_layout_keeps_figure_directories_at_the_root(tmp_path):
    d = layout.subdir(str(tmp_path), "VirtualGaugeFigures", layout="flat")
    assert d == str(tmp_path / "VirtualGaugeFigures")


def test_find_subdir_resolves_both_layouts(tmp_path):
    (tmp_path / "VirtualGaugeFigures").mkdir()
    assert layout.find_subdir(str(tmp_path), "VirtualGaugeFigures") == str(
        tmp_path / "VirtualGaugeFigures")
    nested = tmp_path / "figures" / "VirtualGaugeFigures"
    nested.mkdir(parents=True)
    assert layout.find_subdir(str(tmp_path), "VirtualGaugeFigures") == str(nested)


# --------------------------------------------------------------------------- #
# round-trip: what a writer wrote, a reader finds
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("lay", ["sorted", "flat"])
def test_write_then_find_round_trips_in_either_layout(tmp_path, lay):
    for name in ("ev_i15max.tif", "ev.png", "ev_compare.csv", "ev_gauges.geojson"):
        p = layout.out_path(str(tmp_path), name, layout=lay)
        open(p, "w").write("x")
        assert layout.find(str(tmp_path), name) == p


# --------------------------------------------------------------------------- #
# CLI default output paths must honour the layout too
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("name, sub", [
    ("ev_panels.png", "figures"),      # panels: default was <radar-dir>/<key>_panels.png
    ("ev.png", "figures"),             # map:    default was <out-dir>/<key>.png
    ("ev_compare.csv", "tables"),      # compare
    ("ev_streams.geojson", "vectors"), # export --streams
])
def test_cli_default_output_names_route_into_subdirs(tmp_path, name, sub):
    """These four defaults were built with os.path.join and bypassed the sort,
    dropping figures next to the GeoTIFFs. Caught on a live run, not by a test —
    the commands that write them are network-bound and are not exercised here."""
    assert layout.out_path(str(tmp_path), name) == str(tmp_path / sub / name)


# --------------------------------------------------------------------------- #
# pointing --from-dir at a layout subdirectory still resolves the whole event
# --------------------------------------------------------------------------- #
def _event(tmp_path):
    """A sorted event folder: a raster, the AOI, and a gauge store."""
    for sub, name in [("rasters", "ev_i15max.tif"), ("vectors", "ev_aoi.geojson"),
                      ("vectors", "ev_gauges.geojson")]:
        (tmp_path / sub).mkdir(exist_ok=True)
        (tmp_path / sub / name).write_text("x")
    (tmp_path / "RainGaugeData").mkdir()
    return tmp_path


def test_from_dir_at_rasters_still_finds_the_event_aoi(tmp_path):
    """`--from-dir <event>/rasters` is a natural thing to type -- it is where
    the GeoTIFFs are -- but the AOI and gauge store live in vectors/. Anchoring
    there used to miss them, and the climate figures silently fell back to the
    i15 footprint instead of the event AOI."""
    ev = _event(tmp_path)
    for name in ("ev_i15max.tif", "ev_aoi.geojson", "ev_gauges.geojson"):
        got = layout.find(str(ev / "rasters"), name)
        assert os.path.exists(got), f"{name} not found from rasters/"


def test_from_dir_at_rasters_still_finds_raingaugedata(tmp_path):
    ev = _event(tmp_path)
    assert layout.find_subdir(str(ev / "rasters"), "RainGaugeData") == \
        str(ev / "RainGaugeData")


def test_a_local_hit_still_beats_the_parent(tmp_path):
    """Climbing is only ever a fallback -- a file in the given directory wins."""
    ev = _event(tmp_path)
    (ev / "rasters" / "ev_aoi.geojson").write_text("local")
    assert layout.find(str(ev / "rasters"), "ev_aoi.geojson") == \
        str(ev / "rasters" / "ev_aoi.geojson")


def test_climbing_only_happens_from_a_layout_subdirectory(tmp_path):
    """An ordinary directory must not start reading its parent."""
    (tmp_path / "vectors").mkdir()
    (tmp_path / "vectors" / "ev_aoi.geojson").write_text("x")
    plain = tmp_path / "scratch"
    plain.mkdir()
    got = layout.find(str(plain), "ev_aoi.geojson")
    assert not os.path.exists(got)
