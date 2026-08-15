"""Near-real-time burn severity: archive parsing, AOI screening, severity
classification, and the mosaic rules that keep neighbouring fires from erasing
each other.

Entirely offline — the BRISK archive is replaced by synthetic listings and
locally written GeoTIFFs.
"""
from __future__ import annotations

import datetime as dt
import json
import os

import numpy as np
import pytest

from stormscape import burn


@pytest.fixture
def scene_tif(tmp_path):
    """Factory: a BRISK-shaped EPSG:3857 float32 scene with a NaN surround.

    Mirrors the real files: 60 m cells on the Web-Mercator grid, no nodata tag
    (BRISK marks the outside with bare NaN), origin snapped to a 60 m multiple.
    """
    import rasterio
    from rasterio.transform import from_origin

    def make(arr, name="scene.tif", left=-12807600.0, top=4759020.0, res=60.0,
             crs="EPSG:3857", dtype="float32", nodata=None):
        arr = np.asarray(arr, dtype=dtype)
        path = tmp_path / name
        with rasterio.open(path, "w", driver="GTiff", height=arr.shape[0],
                           width=arr.shape[1], count=1, dtype=dtype, crs=crs,
                           transform=from_origin(left, top, res, res),
                           nodata=nodata) as dst:
            dst.write(arr, 1)
        return str(path)

    return make


# --------------------------------------------------------------------------- #
# archive filenames
# --------------------------------------------------------------------------- #
def test_parse_name_splits_fire_state_and_date():
    r = burn.parse_name("Hidden-Valley-NV-dNBR_20260814_235959.tif")
    assert r["fire"] == "Hidden-Valley"
    assert r["state"] == "NV"
    assert r["date"] == dt.date(2026, 8, 14)


def test_parse_name_keeps_incident_numbered_fires():
    """Some scenes are named by incident number with no fire name at all."""
    r = burn.parse_name("0231-OR-dNBR_20260528_235959.tif")
    assert r["fire"] == "0231" and r["state"] == "OR"


def test_parse_name_tolerates_a_missing_state():
    """The state is informational, so a name without one must still parse."""
    r = burn.parse_name("Somefire-dNBR_20260814_235959.tif")
    assert r is not None
    assert r["fire"] == "Somefire" and r["state"] is None


def test_parse_name_rejects_unrelated_files():
    for n in ("qgis_BRISK_dNBR_colorscale_v2.txt", "notes.pdf",
              "Ward-NV-dNBR_2026081_235959.tif"):
        assert burn.parse_name(n) is None


def test_parse_name_reads_the_soil_burn_severity_token():
    r = burn.parse_name("Cottonwood-Peak-NV-sbs_20260217_235959.tif", "sbs")
    assert r["fire"] == "Cottonwood-Peak" and r["date"] == dt.date(2026, 2, 17)
    # the dNBR parser must not claim an sbs file
    assert burn.parse_name("Cottonwood-Peak-NV-sbs_20260217_235959.tif") is None


def test_as_date_accepts_the_usual_spellings():
    want = dt.date(2026, 8, 14)
    for v in ("20260814", "2026-08-14", dt.date(2026, 8, 14),
              dt.datetime(2026, 8, 14, 23, 59)):
        assert burn._as_date(v) == want


# --------------------------------------------------------------------------- #
# severity classification
# --------------------------------------------------------------------------- #
def test_classify_puts_values_in_the_documented_usgs_classes():
    """USGS/MTBS breaks are published x1000 (100/270/440/660); the archive ships
    unscaled dNBR, so the breaks must be the /1000 form."""
    arr = np.array([[0.0, 0.09, 0.10, 0.26, 0.27, 0.43, 0.44, 0.65, 0.66, 0.9]])
    got = burn.classify(arr, "usgs")
    assert got.tolist() == [[0, 0, 1, 1, 2, 2, 3, 3, 4, 4]]


def test_classify_is_monotone_non_decreasing():
    a = np.linspace(-0.5, 1.5, 200).reshape(1, -1)
    c = burn.classify(a, "usgs")[0]
    assert all(b >= x for x, b in zip(c, c[1:]))


def test_classify_keeps_unburned_distinct_from_unobserved():
    """Class 0 means 'we looked and it did not burn'; NaN means 'no observation'.
    Collapsing the two would silently inflate the unburned area."""
    arr = np.array([[0.02, np.nan]])
    got = burn.classify(arr)
    assert got[0, 0] == 0.0
    assert np.isnan(got[0, 1])


def test_classify_handles_negative_dnbr():
    """Enhanced regrowth is legitimately negative and is not a burn class."""
    assert burn.classify(np.array([[-0.35]]))[0, 0] == 0.0


def test_classify_rejects_an_unknown_scheme():
    with pytest.raises(ValueError, match="scheme must be one of"):
        burn.classify(np.zeros((2, 2)), "nonsense")


@pytest.mark.parametrize("scheme", sorted(burn.SEVERITY_SCHEMES))
def test_every_scheme_labels_each_class_it_can_produce(scheme):
    """N breaks produce N+1 classes, so a missing label would make the class
    table fall back to 'class 3' for a real severity level."""
    spec = burn.SEVERITY_SCHEMES[scheme]
    assert len(spec["labels"]) == len(spec["breaks"]) + 1
    assert list(spec["breaks"]) == sorted(spec["breaks"])
    top = burn.classify(np.array([[9.0]]), scheme)[0, 0]
    assert int(top) == len(spec["breaks"])


# --------------------------------------------------------------------------- #
# mosaic
# --------------------------------------------------------------------------- #
def _bounds_4326(path):
    import rasterio
    from rasterio.warp import transform_bounds
    with rasterio.open(path) as ds:
        return transform_bounds(ds.crs, "EPSG:4326", *ds.bounds)


def test_mosaic_of_one_scene_reproduces_its_values(scene_tif):
    """Same CRS and an aligned grid: the mosaic must be a paste, not a resample —
    otherwise every single-fire map is quietly interpolated."""
    arr = np.arange(48, dtype="float32").reshape(6, 8) / 48.0
    p = scene_tif(arr)
    out, _, crs = burn.mosaic([p], _bounds_4326(p), dst_crs="EPSG:3857")
    assert str(crs) == "EPSG:3857"
    inner = out[np.isfinite(out)]
    assert np.isclose(inner.max(), arr.max())
    assert np.isclose(inner.min(), arr.min())
    # every source value survives unchanged
    assert set(np.round(inner, 6)).issuperset(set(np.round(arr.ravel(), 6)))


def test_mosaic_does_not_let_one_scenes_nan_surround_erase_a_neighbour(scene_tif):
    """The failure this guards against: scenes are mostly NaN outside their own
    fire, so a plain 'last one wins' merge blanks whichever fire is written
    second wherever the other scene overlaps it."""
    a = np.full((6, 6), np.nan, dtype="float32")
    a[:3, :3] = 0.8                      # fire A burns the top-left
    b = np.full((6, 6), np.nan, dtype="float32")
    b[3:, 3:] = 0.5                      # fire B burns the bottom-right
    pa = scene_tif(a, name="a.tif")
    pb = scene_tif(b, name="b.tif")      # identical footprint, disjoint data
    out, _, _ = burn.mosaic([pa, pb], _bounds_4326(pa), dst_crs="EPSG:3857")
    vals = out[np.isfinite(out)].astype("float64")
    assert np.isclose(vals, 0.8).any()               # A survived
    assert np.isclose(vals, 0.5).any()               # ... and so did B
    assert np.isclose(vals, 0.8).sum() >= 9
    assert np.isclose(vals, 0.5).sum() >= 9


def test_mosaic_takes_the_more_severe_value_where_scenes_overlap(scene_tif):
    a = np.full((4, 4), 0.2, dtype="float32")
    b = np.full((4, 4), 0.7, dtype="float32")
    pa, pb = scene_tif(a, name="a.tif"), scene_tif(b, name="b.tif")
    out, _, _ = burn.mosaic([pa, pb], _bounds_4326(pa), dst_crs="EPSG:3857")
    assert np.nanmax(out) == pytest.approx(0.7)
    # order must not matter — fmax is commutative, 'first'/'last' would not be
    out2, _, _ = burn.mosaic([pb, pa], _bounds_4326(pa), dst_crs="EPSG:3857")
    assert np.nanmax(out2) == pytest.approx(0.7)


def test_mosaic_honours_an_explicit_nodata_tag(scene_tif):
    arr = np.array([[1.0, -9999.0], [0.5, 0.25]], dtype="float32")
    p = scene_tif(arr, nodata=-9999.0)
    out, _, _ = burn.mosaic([p], _bounds_4326(p), dst_crs="EPSG:3857")
    assert np.nanmax(out) == pytest.approx(1.0)
    assert not (out[np.isfinite(out)] < 0).any()


def test_categorical_mosaic_drops_the_mask_classes(scene_tif):
    """The BAER rasters' own palette paints 0 and 5+ the same black as 'outside',
    so those are a mask; charting them as a severity class would invent area."""
    arr = np.array([[0, 1, 2], [3, 4, 5]], dtype="uint8")
    p = scene_tif(arr, dtype="uint8")
    out, _, _ = burn.mosaic([p], _bounds_4326(p), dst_crs="EPSG:3857",
                            categorical=True)
    kept = set(np.unique(out[np.isfinite(out)]).tolist())
    assert kept == {1.0, 2.0, 3.0, 4.0}


def test_categorical_mosaic_does_not_interpolate_classes(scene_tif):
    """Averaging class 1 and class 3 into 2 would be a fabricated severity."""
    arr = np.array([[1, 4], [4, 1]], dtype="uint8")
    p = scene_tif(arr, dtype="uint8")
    out, _, _ = burn.mosaic([p], _bounds_4326(p), dst_crs="EPSG:3857",
                            categorical=True)
    assert set(np.unique(out[np.isfinite(out)]).tolist()) <= {1.0, 4.0}


# --------------------------------------------------------------------------- #
# catalog + index caching
# --------------------------------------------------------------------------- #
_LISTING = """<html><body>
<a href="?C=N;O=D">Name</a>
<a href="/pub/realearth/">Parent Directory</a>
<a href="Ward-NV-dNBR_20260813_235959.tif">Ward-NV-dNBR_20260813_235959.tif</a>
<a href="Ward-NV-dNBR_20260814_235959.tif">Ward-NV-dNBR_20260814_235959.tif</a>
<a href="Willow-CO-dNBR_20260814_235959.tif">Willow-CO-dNBR_20260814_235959.tif</a>
<a href="qgis_BRISK_dNBR_colorscale_v2.txt">colorscale</a>
</body></html>"""


def test_list_year_parses_the_directory_index(tmp_path, monkeypatch):
    monkeypatch.setattr(burn, "_list_url", lambda url, timeout=60.0:
                        burn._HREF.findall(_LISTING))
    df = burn.list_year(2026, cache_dir=str(tmp_path), verbose=False)
    assert len(df) == 3                      # the .txt is not a scene
    assert set(df.fire) == {"Ward", "Willow"}
    assert df.url.iloc[0].startswith(burn.BRISK_BASE + "/2026/")


def test_a_past_year_index_is_cached_forever(tmp_path, monkeypatch):
    """Closed years never gain scenes, so re-listing them is pure latency."""
    calls = []

    def fake(url, timeout=60.0):
        calls.append(url)
        return burn._HREF.findall(_LISTING)

    monkeypatch.setattr(burn, "_list_url", fake)
    past = dt.date.today().year - 1
    burn.list_year(past, cache_dir=str(tmp_path), verbose=False)
    burn.list_year(past, cache_dir=str(tmp_path), verbose=False)
    assert len(calls) == 1

    # and the cached frame round-trips its dates as dates, not strings
    df = burn.list_year(past, cache_dir=str(tmp_path), verbose=False)
    assert isinstance(df.date.iloc[0], dt.date)


def test_the_current_year_index_expires(tmp_path, monkeypatch):
    """'Near real-time' is only true if today's listing is allowed to go stale."""
    calls = []

    def fake(url, timeout=60.0):
        calls.append(url)
        return burn._HREF.findall(_LISTING)

    monkeypatch.setattr(burn, "_list_url", fake)
    now = dt.date.today().year
    burn.list_year(now, cache_dir=str(tmp_path), verbose=False)
    burn.list_year(now, cache_dir=str(tmp_path), ttl_h=6.0, verbose=False)
    assert len(calls) == 1                       # still fresh
    burn.list_year(now, cache_dir=str(tmp_path), ttl_h=-1.0, verbose=False)
    assert len(calls) == 2                       # forced stale -> re-listed


def test_list_year_falls_back_to_a_stale_cache_when_the_server_is_down(
        tmp_path, monkeypatch):
    """A dead archive should degrade to yesterday's catalog, not kill the run."""
    monkeypatch.setattr(burn, "_list_url", lambda url, timeout=60.0:
                        burn._HREF.findall(_LISTING))
    now = dt.date.today().year
    burn.list_year(now, cache_dir=str(tmp_path), verbose=False)

    def boom(url, timeout=60.0):
        raise OSError("connection refused")

    monkeypatch.setattr(burn, "_list_url", boom)
    df = burn.list_year(now, cache_dir=str(tmp_path), ttl_h=-1.0, verbose=False)
    assert len(df) == 3


def test_list_year_raises_when_there_is_no_cache_to_fall_back_on(
        tmp_path, monkeypatch):
    def boom(url, timeout=60.0):
        raise OSError("connection refused")

    monkeypatch.setattr(burn, "_list_url", boom)
    with pytest.raises(RuntimeError, match="failed to list"):
        burn.list_year(2026, cache_dir=str(tmp_path), verbose=False)


def test_unknown_product_is_rejected():
    with pytest.raises(ValueError, match="product must be one of"):
        burn._product("landsat")


# --------------------------------------------------------------------------- #
# AOI screening
# --------------------------------------------------------------------------- #
def _stub_catalog(monkeypatch, rows):
    import pandas as pd
    df = pd.DataFrame(rows)
    monkeypatch.setattr(burn, "catalog",
                        lambda *a, **k: df.copy())


def test_find_scenes_keeps_only_fires_that_touch_the_aoi(tmp_path, monkeypatch):
    _stub_catalog(monkeypatch, [
        dict(fire="Near", state="NV", date=dt.date(2026, 8, 14),
             filename="a.tif", url="u_near"),
        dict(fire="Far", state="ME", date=dt.date(2026, 8, 14),
             filename="b.tif", url="u_far"),
    ])
    monkeypatch.setattr(burn, "scene_bounds", lambda *a, **k: {
        "u_near": (-115.0, 39.0, -114.7, 39.3),
        "u_far": (-70.0, 44.0, -69.7, 44.3),
    })
    got = burn.find_scenes((-115.0, 39.05, -114.8, 39.22),
                           cache_dir=str(tmp_path), verbose=False)
    assert list(got.fire) == ["Near"]


def test_find_scenes_takes_the_latest_scene_per_fire(tmp_path, monkeypatch):
    """A fire has one file per day; mapping all of them would stack the same
    scar repeatedly and pick an arbitrary winner."""
    _stub_catalog(monkeypatch, [
        dict(fire="Ward", state="NV", date=dt.date(2026, 8, 12),
             filename="a", url="u_old"),
        dict(fire="Ward", state="NV", date=dt.date(2026, 8, 14),
             filename="b", url="u_new"),
    ])
    monkeypatch.setattr(burn, "scene_bounds", lambda urls, *a, **k:
                        {u: (-115.0, 39.0, -114.7, 39.3) for u in urls})
    got = burn.find_scenes((-115.0, 39.05, -114.8, 39.22),
                           cache_dir=str(tmp_path), verbose=False)
    assert len(got) == 1 and got.url.iloc[0] == "u_new"


def test_find_scenes_respects_an_as_of_date(tmp_path, monkeypatch):
    """Asking for the scar as it looked on a storm date must not return a scene
    published after the storm."""
    _stub_catalog(monkeypatch, [
        dict(fire="Ward", state="NV", date=dt.date(2026, 8, 10),
             filename="a", url="u_before"),
        dict(fire="Ward", state="NV", date=dt.date(2026, 8, 20),
             filename="b", url="u_after"),
    ])
    monkeypatch.setattr(burn, "scene_bounds", lambda urls, *a, **k:
                        {u: (-115.0, 39.0, -114.7, 39.3) for u in urls})
    got = burn.find_scenes((-115.0, 39.05, -114.8, 39.22), date="20260814",
                           cache_dir=str(tmp_path), verbose=False)
    assert list(got.url) == ["u_before"]


def test_find_scenes_can_drop_older_seasons(tmp_path, monkeypatch):
    _stub_catalog(monkeypatch, [
        dict(fire="LastYear", state="NV", date=dt.date(2025, 9, 1),
             filename="a", url="u_old"),
        dict(fire="ThisYear", state="NV", date=dt.date(2026, 8, 14),
             filename="b", url="u_new"),
    ])
    monkeypatch.setattr(burn, "scene_bounds", lambda urls, *a, **k:
                        {u: (-115.0, 39.0, -114.7, 39.3) for u in urls})
    got = burn.find_scenes((-115.0, 39.05, -114.8, 39.22), since="20260101",
                           cache_dir=str(tmp_path), verbose=False)
    assert list(got.fire) == ["ThisYear"]


def test_find_scenes_returns_an_empty_frame_for_an_empty_catalog(
        tmp_path, monkeypatch):
    """An AOI that never burned is the common case, not an error."""
    _stub_catalog(monkeypatch, [])
    got = burn.find_scenes((-115.0, 39.05, -114.8, 39.22),
                           cache_dir=str(tmp_path), verbose=False)
    assert len(got) == 0
    assert "geometry" in got.columns


def test_find_scenes_returns_an_empty_frame_when_no_fire_reaches_the_aoi(
        tmp_path, monkeypatch):
    """The catalog is full of fires and none of them are here — the usual case
    for an unburned AOI, and a different code path from an empty catalog."""
    _stub_catalog(monkeypatch, [
        dict(fire="Elsewhere", state="OR", date=dt.date(2026, 8, 14),
             filename="a", url="u_far"),
    ])
    monkeypatch.setattr(burn, "scene_bounds", lambda urls, *a, **k:
                        {u: (-120.0, 44.0, -119.7, 44.3) for u in urls})
    got = burn.find_scenes((-73.99, 40.70, -73.95, 40.75),
                           cache_dir=str(tmp_path), verbose=False)
    assert len(got) == 0
    assert "geometry" in got.columns and "fire" in got.columns


def test_burn_severity_returns_none_for_an_unburned_aoi(tmp_path, monkeypatch):
    _stub_catalog(monkeypatch, [
        dict(fire="Elsewhere", state="OR", date=dt.date(2026, 8, 14),
             filename="a", url="u_far"),
    ])
    monkeypatch.setattr(burn, "scene_bounds", lambda urls, *a, **k:
                        {u: (-120.0, 44.0, -119.7, 44.3) for u in urls})
    got = burn.burn_severity((-73.99, 40.70, -73.95, 40.75),
                             cache_dir=str(tmp_path), verbose=False)
    assert got is None


def test_scene_age_is_measured_from_the_fires_first_appearance(
        tmp_path, monkeypatch):
    """Age proxies how much post-fire imagery the composite has absorbed, so it
    must count from when the fire entered the archive."""
    _stub_catalog(monkeypatch, [
        dict(fire="Ward", state="NV", date=dt.date(2026, 8, 1),
             filename="a", url="u1"),
        dict(fire="Ward", state="NV", date=dt.date(2026, 8, 20),
             filename="b", url="u2"),
    ])
    monkeypatch.setattr(burn, "scene_bounds", lambda urls, *a, **k:
                        {u: (-115.0, 39.0, -114.7, 39.3) for u in urls})
    got = burn.find_scenes((-115.0, 39.05, -114.8, 39.22),
                           cache_dir=str(tmp_path), verbose=False)
    assert got.age_days.iloc[0] == 19


def test_scene_age_ignores_the_date_filter(tmp_path, monkeypatch):
    """Age is computed before --date/--since trim the catalog: trimming the early
    scenes would make a long-burning fire look brand new and wrongly immature."""
    _stub_catalog(monkeypatch, [
        dict(fire="Ward", state="NV", date=dt.date(2026, 8, 1),
             filename="a", url="u1"),
        dict(fire="Ward", state="NV", date=dt.date(2026, 8, 20),
             filename="b", url="u2"),
    ])
    monkeypatch.setattr(burn, "scene_bounds", lambda urls, *a, **k:
                        {u: (-115.0, 39.0, -114.7, 39.3) for u in urls})
    got = burn.find_scenes((-115.0, 39.05, -114.8, 39.22), since="20260815",
                           cache_dir=str(tmp_path), verbose=False)
    assert len(got) == 1 and got.age_days.iloc[0] == 19      # not 0


def test_min_age_drops_an_immature_composite(tmp_path, monkeypatch):
    """An immature composite reads LOW, which for post-fire hazard work is the
    dangerous direction, so --min-age must actually exclude it."""
    _stub_catalog(monkeypatch, [
        dict(fire="Fresh", state="NV", date=dt.date(2026, 8, 3),
             filename="a", url="u_new"),
        dict(fire="Old", state="NV", date=dt.date(2026, 6, 1),
             filename="b", url="u_old1"),
        dict(fire="Old", state="NV", date=dt.date(2026, 8, 3),
             filename="c", url="u_old2"),
    ])
    monkeypatch.setattr(burn, "scene_bounds", lambda urls, *a, **k:
                        {u: (-115.0, 39.0, -114.7, 39.3) for u in urls})
    aoi = (-115.0, 39.05, -114.8, 39.22)
    both = burn.find_scenes(aoi, cache_dir=str(tmp_path), verbose=False)
    assert set(both.fire) == {"Fresh", "Old"}
    mature = burn.find_scenes(aoi, min_age_days=14, cache_dir=str(tmp_path),
                              verbose=False)
    assert list(mature.fire) == ["Old"]


def test_min_age_reports_only_fires_in_this_aoi(tmp_path, monkeypatch, capsys):
    """Screening maturity across the whole national catalog would name hundreds
    of irrelevant fires; it must happen after the AOI intersection."""
    _stub_catalog(monkeypatch, [
        dict(fire="Here", state="NV", date=dt.date(2026, 8, 3),
             filename="a", url="u_here"),
        dict(fire="Elsewhere", state="ME", date=dt.date(2026, 8, 3),
             filename="b", url="u_far"),
    ])
    monkeypatch.setattr(burn, "scene_bounds", lambda urls, *a, **k: {
        "u_here": (-115.0, 39.0, -114.7, 39.3),
        "u_far": (-70.0, 44.0, -69.7, 44.3)})
    burn.find_scenes((-115.0, 39.05, -114.8, 39.22), min_age_days=14,
                     cache_dir=str(tmp_path), verbose=True)
    out = capsys.readouterr().out
    assert "Here" in out and "Elsewhere" not in out


def test_min_age_returns_a_usable_empty_frame(tmp_path, monkeypatch):
    _stub_catalog(monkeypatch, [
        dict(fire="Fresh", state="NV", date=dt.date(2026, 8, 3),
             filename="a", url="u_new"),
    ])
    monkeypatch.setattr(burn, "scene_bounds", lambda urls, *a, **k:
                        {u: (-115.0, 39.0, -114.7, 39.3) for u in urls})
    got = burn.find_scenes((-115.0, 39.05, -114.8, 39.22), min_age_days=99,
                           cache_dir=str(tmp_path), verbose=False)
    assert len(got) == 0 and "age_days" in got.columns


def test_maturity_threshold_is_the_validated_one():
    """14 d is where agreement with BAER's own dNBR settles (examples/
    brisk_vs_baer.py); it is a measured number, not a round guess."""
    assert burn.MATURITY_DAYS == 14


# --------------------------------------------------------------------------- #
# the BAER teams' own dNBR
# --------------------------------------------------------------------------- #
def test_baer_dnbr_is_a_known_product_scaled_by_1000():
    """BAER ships dNBR as int16 x1000 (the BARC convention). Forgetting the
    divide leaves every value 1000x high while the correlation still looks fine."""
    cfg = burn._product("baer_dnbr")
    assert cfg["base"] == burn.BAER_BASE
    assert cfg["scale"] == 1000.0
    assert burn._product("dnbr")["scale"] == 1.0


def test_parse_name_strips_the_baer_prelim_infix():
    """Left in, the state code is lost and the fire stops matching its own BRISK
    scenes -- which is exactly what the cross-product comparison needs."""
    r = burn.parse_name("Alder-Springs-OR-prelim-dNBR_20250625_235959.tif")
    assert r["fire"] == "Alder-Springs" and r["state"] == "OR"


def test_scene_bounds_memoises_header_reads(tmp_path, monkeypatch):
    """Scenes are immutable once published, so a cached bound is always valid —
    this is what makes repeat AOI screens instant instead of a network sweep."""
    calls = []

    def fake(url):
        calls.append(url)
        return (-115.0, 39.0, -114.7, 39.3)

    monkeypatch.setattr(burn, "_read_bounds", fake)
    urls = ["u1", "u2"]
    first = burn.scene_bounds(urls, cache_dir=str(tmp_path), verbose=False)
    second = burn.scene_bounds(urls, cache_dir=str(tmp_path), verbose=False)
    assert first == second and len(calls) == 2
    assert os.path.exists(tmp_path / "bounds_dnbr.json")


def test_scene_bounds_skips_a_scene_whose_header_will_not_read(
        tmp_path, monkeypatch):
    """One corrupt scene must not take the whole AOI screen down with it."""
    def fake(url):
        if url == "bad":
            raise OSError("truncated")
        return (-115.0, 39.0, -114.7, 39.3)

    monkeypatch.setattr(burn, "_read_bounds", fake)
    got = burn.scene_bounds(["good", "bad"], cache_dir=str(tmp_path),
                            verbose=False)
    assert set(got) == {"good"}


def test_scene_bounds_survives_a_corrupt_cache(tmp_path, monkeypatch):
    (tmp_path / "bounds_dnbr.json").write_text("{not json")
    monkeypatch.setattr(burn, "_read_bounds",
                        lambda url: (-115.0, 39.0, -114.7, 39.3))
    got = burn.scene_bounds(["u1"], cache_dir=str(tmp_path), verbose=False)
    assert set(got) == {"u1"}
    json.loads((tmp_path / "bounds_dnbr.json").read_text())   # rewritten valid


# --------------------------------------------------------------------------- #
# fetch caching
# --------------------------------------------------------------------------- #
def test_fetch_scene_reuses_a_cached_file(tmp_path, monkeypatch):
    calls = []

    def fake(url, path):
        calls.append(url)
        with open(path, "wb") as fh:
            fh.write(b"tif")

    monkeypatch.setattr(burn.urllib.request, "urlretrieve", fake)
    url = "https://example.invalid/Ward-NV-dNBR_20260814_235959.tif"
    p1 = burn.fetch_scene(url, str(tmp_path), verbose=False)
    p2 = burn.fetch_scene(url, str(tmp_path), verbose=False)
    assert p1 == p2 and len(calls) == 1


def test_fetch_scene_leaves_no_partial_file_behind(tmp_path, monkeypatch):
    """A half-downloaded scene cached under the real name would be reused
    forever and silently corrupt every later map."""
    def boom(url, path):
        with open(path, "wb") as fh:
            fh.write(b"half")
        raise OSError("connection reset")

    monkeypatch.setattr(burn.urllib.request, "urlretrieve", boom)
    monkeypatch.setattr(burn.time, "sleep", lambda s: None)
    with pytest.raises(RuntimeError, match="failed to download"):
        burn.fetch_scene("https://example.invalid/x-dNBR_20260814_235959.tif",
                         str(tmp_path), retries=1, verbose=False)
    assert list(tmp_path.iterdir()) == []


# --------------------------------------------------------------------------- #
# colour ramp
# --------------------------------------------------------------------------- #
def test_brisk_colormap_reproduces_the_published_anchors():
    """Users cross-check these maps against the CIMSS portal, so the ramp has to
    land on the portal's own colours at its own breaks."""
    cmap = burn.register_brisk_cmap("brisk_test_anchors")
    for value, rgb in ((0.40, (255, 232, 32)), (0.70, (168, 0, 0)),
                       (1.00, (114, 0, 0))):
        got = tuple(int(round(255 * c)) for c in cmap(value)[:3])
        assert all(abs(a - b) <= 2 for a, b in zip(got, rgb)), (value, got, rgb)


def test_brisk_colormap_keeps_its_step_at_the_unburned_break():
    """The portal ramp jumps colour at 0.10 rather than fading through it, which
    is what makes 'burned at all' readable on the map. A coarse lookup table
    would smear that edge across a visible band of dNBR."""
    cmap = burn.register_brisk_cmap("brisk_test_step")
    below = np.array(cmap(0.095)[:3])
    above = np.array(cmap(0.105)[:3])
    assert np.abs(above - below).max() > 0.1
    # and it really is a step: barely-unburned matches fully-unburned
    assert np.allclose(below, np.array(cmap(0.0)[:3]))


def test_register_brisk_cmap_is_idempotent():
    a = burn.register_brisk_cmap("brisk_test_twice")
    b = burn.register_brisk_cmap("brisk_test_twice")
    assert a is not None and b is not None


def test_the_ramp_anchors_are_the_baer_class_colours():
    """Verified against the products: all 77 of the 2025 BAER soil-burn-severity
    rasters embed this identical palette, and it is what BRISK publishes. If
    these drift apart, a stormscape map no longer matches a BAER deliverable."""
    ramp = {rgb for _, rgb in burn.BAER_ANCHORS}
    for rgb in burn.BAER_CLASS_COLORS:
        assert rgb in ramp
    assert burn.BRISK_ANCHORS is burn.BAER_ANCHORS      # one scheme, two names


def test_severity_colors_uses_the_official_table_for_four_classes():
    """A four-class scheme maps exactly onto the four published BAER colours —
    no interpolation, so the map is byte-comparable with a BAER product."""
    cmap, norm, ticks, labels = burn.severity_colors("brisk")
    assert len(labels) == 4 and len(ticks) == 4
    got = [tuple(int(round(255 * c)) for c in cmap(i)[:3]) for i in range(4)]
    assert got == list(burn.BAER_CLASS_COLORS)


def test_severity_colors_bands_at_the_scheme_breaks():
    """Each class must occupy its own band, so a dNBR either side of a break
    lands in different colours — that is what 'classed' means."""
    _, norm, _, _ = burn.severity_colors("brisk")
    breaks = burn.SEVERITY_SCHEMES["brisk"]["breaks"]
    for b in breaks:
        assert norm(b - 1e-4) != norm(b + 1e-4)
    assert norm(0.0) == 0                    # unburned
    assert norm(0.95) == len(breaks)         # top class


def test_severity_colors_covers_dnbr_beyond_the_unit_interval():
    """dNBR is legitimately negative (enhanced regrowth) and can exceed 1;
    neither may fall off the end of the norm and render as an out-of-range hole."""
    _, norm, _, _ = burn.severity_colors("usgs")
    n = len(burn.SEVERITY_SCHEMES["usgs"]["labels"])
    for v in (-0.35, -0.05, 0.0, 0.5, 1.0, 1.51):
        assert 0 <= int(norm(v)) <= n - 1


def test_severity_colors_supports_the_five_class_scheme():
    cmap, norm, ticks, labels = burn.severity_colors("usgs")
    assert len(labels) == 5 and len(ticks) == 5
    assert cmap.N == 5


def test_severity_colors_ticks_sit_inside_their_own_band():
    """A tick that lands outside its class would label the wrong colour."""
    for scheme in sorted(burn.SEVERITY_SCHEMES):
        _, norm, ticks, labels = burn.severity_colors(scheme)
        for i, t in enumerate(ticks):
            assert int(norm(t)) == i, (scheme, i, t)


def test_severity_colors_rejects_an_unknown_scheme():
    with pytest.raises(ValueError, match="scheme must be one of"):
        burn.severity_colors("nonsense")


def test_drape_accepts_a_norm_and_class_tick_labels(tmp_path, field_tif):
    """The classed burn map rides on drape_i15, so it must take a norm instead of
    vmin/vmax and put names on the colour bar."""
    import matplotlib.pyplot as plt

    from stormscape.plot import drape_i15
    arr = np.array([[0.05, 0.2], [0.5, 0.8]], dtype="float32")
    p = field_tif(arr, name="d.tif")
    cmap, norm, ticks, labels = burn.severity_colors("brisk")
    out = str(tmp_path / "m.png")
    fig, ax = drape_i15(None, p, out_path=out, work_crs="EPSG:4326",
                        wet_min=0.10, cmap=cmap, norm=norm,
                        cbar_ticks=ticks, cbar_ticklabels=labels,
                        scale_ticks=False, north_arrow=False, legend=None)
    plt.close(fig)
    assert os.path.exists(out)


# --------------------------------------------------------------------------- #
# downstream wiring
# --------------------------------------------------------------------------- #
def test_export_scales_burn_fields_on_their_own_range():
    """dNBR tops out near 1; the rainfall autoscale floor of 10 would render an
    entire scar in the bottom tenth of the colour ramp."""
    from stormscape import export
    arr = np.array([[0.0, 0.3, 0.9]])
    cmap, vmin, vmax, norm, mask_below = export._field_style("dnbr", arr)
    assert vmax == pytest.approx(1.0)
    assert mask_below == pytest.approx(0.10)
    assert norm is None


def test_export_treats_a_severity_class_field_as_categorical():
    """Interpolating between class 2 and class 4 would invent a severity."""
    from stormscape import export
    assert export._is_categorical("severity")
    assert export._is_categorical("Ward_severity")


def test_severity_class_field_gets_a_class_scale():
    from stormscape import export
    arr = np.array([[0.0, 1.0, 2.0, 3.0, 4.0]])
    _, _, vmax, _, mask_below = export._field_style("severity", arr)
    assert vmax == pytest.approx(4.0)
    assert mask_below == pytest.approx(0.5)      # class 0 = unburned, see-through


def test_burn_fields_have_panel_specs():
    """`panels` reads these by field name; a missing entry drops the panel."""
    from stormscape.plot import _PANEL_SPECS
    assert "dnbr" in _PANEL_SPECS and "severity" in _PANEL_SPECS


def test_burn_fields_are_not_in_the_rainfall_dry_mask():
    """_MASK_DRY cuts below 0.5, which on a 0-1 dNBR field would erase
    everything short of high severity."""
    from stormscape.plot import _MASK_DRY
    assert "dnbr" not in _MASK_DRY and "severity" not in _MASK_DRY


def test_brisk_cache_stays_at_the_event_root():
    """Caches are inputs, not products: sorting them into rasters/ would mix
    hundreds of raw scenes in with the event's own outputs."""
    from stormscape import layout
    assert "brisk_cache" in layout.RESERVED
    assert layout.subdir(".", "brisk_cache", make=False).endswith("brisk_cache")


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def test_burn_subcommand_parses():
    from stormscape.cli import main
    with pytest.raises(SystemExit):
        main(["burn", "--help"])


def _run_burn_cli(monkeypatch, argv):
    """Dispatch the `burn` parser into a stand-in, capturing the parsed args.

    ``main()`` looks ``_cmd_burn`` up as a module global while it builds the
    parser, so patching it first keeps the command entirely offline.
    """
    from stormscape import cli
    seen = {}
    monkeypatch.setattr(cli, "_cmd_burn", lambda args: seen.update(vars(args)))
    cli.main(argv)
    return seen


def test_burn_cli_defaults_suit_severity_not_rainfall(monkeypatch):
    """_add_overlays supplies YlGnBu and a 5 mm/h dry cut; both are wrong for a
    0-1 dNBR field, so the parser has to override them. The default palette is
    BAER's, so a stormscape map can be laid beside a BAER product."""
    got = _run_burn_cli(monkeypatch,
                        ["burn", "--bbox", "-115", "39", "-114", "40"])
    assert got["cmap"] == "baer"
    assert got["wet_min"] is None          # resolved per product, not 5 mm/h
    assert got["vmax"] is None
    assert got["continuous"] is False      # classed, the way BAER publishes
    assert got["product"] == "dnbr" and got["scheme"] == "usgs"


def test_burn_display_defaults_differ_by_product():
    """One shared scale would either clip the 1-4 class map to its lowest class
    or paint unburned ground on the dNBR map."""
    from stormscape.cli import _burn_display_defaults
    dnbr_cut, dnbr_vmax = _burn_display_defaults("dnbr")
    sbs_cut, sbs_vmax = _burn_display_defaults("sbs")
    assert (dnbr_cut, dnbr_vmax) == pytest.approx((0.10, 1.0))
    assert sbs_vmax >= 4.0                 # must reach class 4 (high)
    assert 1.0 < sbs_cut < 2.0             # cuts class 1, keeps class 2+


def test_burn_cli_accepts_the_portal_ramp_and_the_baer_product(monkeypatch):
    got = _run_burn_cli(monkeypatch,
                        ["burn", "--bbox", "-115", "39", "-114", "40",
                         "--product", "sbs", "--scheme", "brisk",
                         "--cmap", "brisk", "--date", "20260814",
                         "--since", "20260101", "--fire", "Ward", "Willow"])
    assert got["product"] == "sbs" and got["scheme"] == "brisk"
    assert got["cmap"] == "brisk" and got["fire"] == ["Ward", "Willow"]
    assert got["date"] == "20260814" and got["since"] == "20260101"
