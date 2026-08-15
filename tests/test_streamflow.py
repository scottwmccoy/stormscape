"""USGS stream gauges: WaterML parsing, unit and timezone handling, summaries.

Entirely offline -- ``requests.get`` is replaced by a double replaying canned
NWIS and OGC payloads, so nothing here touches the network.

The invariants pinned here are the ones that quietly corrupt a hydrograph:
NWIS stamps local time with an offset rather than UTC, missing values arrive as
-999999 rather than null, the site service returns one row per time series so a
gauge appears several times, and the two USGS backends must return the same
gauges or switching sources silently changes the answer.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from stormscape import streamflow as sf


# --------------------------------------------------------------------------- #
# doubles
# --------------------------------------------------------------------------- #
class _Resp:
    def __init__(self, payload=None, text="", status=200):
        self._payload, self.text, self.status_code = payload, text, status

    def json(self):
        if self._payload is None:
            raise ValueError("not json")
        return self._payload


def _install(monkeypatch, handler):
    """Replace requests.get with ``handler(url, params) -> _Resp``."""
    monkeypatch.setattr(sf.requests, "get",
                        lambda url, params=None, timeout=None:
                        handler(url, dict(params or {})))


def _iv_payload(site="10348000", rows=None, code="00060"):
    rows = rows or [("2026-08-12T17:00:00.000-07:00", "330"),
                    ("2026-08-12T17:05:00.000-07:00", "335")]
    return {"value": {"timeSeries": [{
        "sourceInfo": {"siteCode": [{"value": site}], "siteName": "TEST RIVER"},
        "variable": {"variableCode": [{"value": code}],
                     "variableName": "Streamflow",
                     "unit": {"unitCode": "ft3/s"}, "noDataValue": -999999.0},
        "values": [{"value": [{"dateTime": t, "value": v} for t, v in rows]}],
    }]}}


SITE_RDB = "\n".join([
    "# comment line",
    "agency_cd\tsite_no\tstation_nm\tsite_tp_cd\tdec_lat_va\tdec_long_va\t"
    "alt_va\thuc_cd\tdrain_area_va",
    "5s\t15s\t50s\t7s\t16s\t16s\t8s\t16s\t8s",
    "USGS\t10348000\tTRUCKEE RV AT RENO, NV\tST\t39.5301\t-119.7954\t4448\t"
    "16050102\t1067",
    # the site service repeats a gauge once per available time series
    "USGS\t10348000\tTRUCKEE RV AT RENO, NV\tST\t39.5301\t-119.7954\t4448\t"
    "16050102\t1067",
    "USGS\t10311750\tCARSON R BLW DAYTON, NV\tST\t39.2809\t-119.5251\t4295\t"
    "16050202\t1121",
])


# --------------------------------------------------------------------------- #
# RDB parsing + sites
# --------------------------------------------------------------------------- #
def test_rdb_drops_comments_and_the_format_row():
    df = sf._rdb_to_frame(SITE_RDB)
    assert len(df) == 3                       # 2 Truckee rows + 1 Carson
    assert "5s" not in df.site_no.tolist()    # the format row is not data


def test_rdb_on_empty_input_is_an_empty_frame():
    assert not len(sf._rdb_to_frame(""))
    assert not len(sf._rdb_to_frame("# only a comment"))


def test_sites_dedupe_one_row_per_time_series(monkeypatch):
    """A gauge measuring discharge, stage and temperature appears three times;
    counting the rows would treble the gauge count."""
    _install(monkeypatch, lambda url, p: _Resp(text=SITE_RDB))
    g = sf.stream_sites((-120, 39, -119, 40))
    assert len(g) == 2
    assert sorted(g.site_no) == ["10311750", "10348000"]


def test_sites_convert_drainage_area_to_km2(monkeypatch):
    _install(monkeypatch, lambda url, p: _Resp(text=SITE_RDB))
    g = sf.stream_sites((-120, 39, -119, 40)).set_index("site_no")
    assert g.loc["10348000", "drain_area_km2"] == pytest.approx(
        1067 * 2.589988110336)


def test_sites_are_clipped_to_the_aoi(monkeypatch):
    """The service bbox filter is generous at the edges; the Carson gauge sits
    outside this tighter box and must not come back."""
    _install(monkeypatch, lambda url, p: _Resp(text=SITE_RDB))
    g = sf.stream_sites((-120, 39.4, -119, 40))
    assert sorted(g.site_no) == ["10348000"]


def test_sites_never_send_the_mutually_exclusive_pair(monkeypatch):
    """`outputDataTypeCd` together with `siteOutput=expanded` is an HTTP 400."""
    seen = {}

    def h(url, p):
        seen.update(p)
        return _Resp(text=SITE_RDB)

    _install(monkeypatch, h)
    sf.stream_sites((-120, 39, -119, 40))
    assert not ("outputDataTypeCd" in seen and "siteOutput" in seen)


def test_unknown_source_is_rejected():
    with pytest.raises(ValueError, match="unknown source"):
        sf.stream_sites((-120, 39, -119, 40), source="nope")


def test_empty_sites_carry_the_full_schema(monkeypatch):
    _install(monkeypatch, lambda url, p: _Resp(text=""))
    g = sf.stream_sites((-120, 39, -119, 40))
    assert not len(g)
    for c in ("site_no", "name", "drain_area_km2", "source"):
        assert c in g.columns
    assert g.crs is not None


def test_a_failed_site_query_warns_rather_than_returning_a_quiet_empty(
        monkeypatch, capsys):
    _install(monkeypatch, lambda url, p: _Resp(text="", status=500))
    g = sf.stream_sites((-120, 39, -119, 40))
    assert not len(g)
    assert "nwis sites" in capsys.readouterr().err


def test_404_is_no_sites_not_a_failure(monkeypatch, capsys):
    """NWIS answers 404 when nothing matches, which is a real answer."""
    _install(monkeypatch, lambda url, p: _Resp(text="", status=404))
    assert not len(sf.stream_sites((-120, 39, -119, 40)))
    assert "warning" not in capsys.readouterr().err


# --------------------------------------------------------------------------- #
# series: timezone, sentinels, units
# --------------------------------------------------------------------------- #
def test_local_timestamps_are_normalised_to_utc(monkeypatch):
    """NWIS stamps the gauge's local time WITH an offset. Read naively, a
    Pacific gauge's hydrograph lands 7 hours from the radar that caused it."""
    _install(monkeypatch, lambda url, p: _Resp(_iv_payload()))
    s = sf.stream_series(["10348000"], "2026-08-13T00:00Z", "2026-08-13T01:00Z")
    idx = s["10348000"].index
    assert str(idx.tz) == "UTC"
    assert idx[0] == pd.Timestamp("2026-08-13T00:00:00Z")


def test_nodata_sentinel_becomes_missing(monkeypatch):
    rows = [("2026-08-13T00:00:00.000+00:00", "330"),
            ("2026-08-13T00:05:00.000+00:00", "-999999"),
            ("2026-08-13T00:10:00.000+00:00", "340")]
    _install(monkeypatch, lambda url, p: _Resp(_iv_payload(rows=rows)))
    df = sf.stream_series(["10348000"], "2026-08-13T00:00Z",
                          "2026-08-13T01:00Z")["10348000"]
    assert len(df) == 2                       # the sentinel row is dropped
    assert df.discharge_cfs.max() == 340


def test_both_unit_systems_are_present_and_exact(monkeypatch):
    _install(monkeypatch, lambda url, p: _Resp(_iv_payload()))
    df = sf.stream_series(["10348000"], "2026-08-13T00:00Z",
                          "2026-08-13T01:00Z")["10348000"]
    assert {"discharge_cms", "discharge_cfs"} <= set(df.columns)
    assert df.discharge_cms.iloc[0] == pytest.approx(330 * 0.028316846592)


def test_stage_is_converted_to_metres(monkeypatch):
    _install(monkeypatch, lambda url, p: _Resp(_iv_payload(
        rows=[("2026-08-13T00:00:00.000+00:00", "4.18")], code="00065")))
    df = sf.stream_series(["10348000"], "2026-08-13T00:00Z",
                          "2026-08-13T01:00Z")["10348000"]
    assert df.stage_m.iloc[0] == pytest.approx(4.18 * 0.3048)


def test_series_are_sorted_and_deduplicated(monkeypatch):
    rows = [("2026-08-13T00:10:00.000+00:00", "3"),
            ("2026-08-13T00:00:00.000+00:00", "1"),
            ("2026-08-13T00:10:00.000+00:00", "9")]     # duplicate stamp
    _install(monkeypatch, lambda url, p: _Resp(_iv_payload(rows=rows)))
    df = sf.stream_series(["10348000"], "2026-08-13T00:00Z",
                          "2026-08-13T01:00Z")["10348000"]
    assert df.index.is_monotonic_increasing
    assert len(df) == 2 and df.discharge_cfs.iloc[-1] == 9   # last wins


def test_empty_site_list_short_circuits(monkeypatch):
    def boom(url, p):
        raise AssertionError("should not have been called")

    _install(monkeypatch, boom)
    assert sf.stream_series([], "2026-08-13T00:00Z", "2026-08-13T01:00Z") == {}


def test_series_window_is_sent_as_utc(monkeypatch):
    seen = {}

    def h(url, p):
        seen.update(p)
        return _Resp(_iv_payload())

    _install(monkeypatch, h)
    sf.stream_series(["10348000"], "2026-08-13T00:00Z", "2026-08-13T06:00Z")
    assert seen["startDT"].endswith("Z") and seen["endDT"].endswith("Z")
    assert seen["parameterCd"] == "00060,00065"


# --------------------------------------------------------------------------- #
# summaries
# --------------------------------------------------------------------------- #
def _series(vals, freq="5min", start="2026-08-13T00:00Z"):
    idx = pd.date_range(start, periods=len(vals), freq=freq, tz="UTC")
    q = pd.Series(vals, index=idx, dtype="float64")
    return pd.DataFrame({"discharge_cfs": q, "discharge_cms": q * 0.028316846592,
                         "stage_ft": q / 100.0, "stage_m": q / 100.0 * 0.3048})


def _sites_gdf():
    import geopandas as gpd
    return gpd.GeoDataFrame(
        dict(site_no=["A"], name=["TEST"], drain_area_km2=[100.0]),
        geometry=gpd.points_from_xy([-119.5], [39.5]), crs=4326)


def test_summary_reports_peak_and_time_of_peak():
    s = {"A": _series([10, 50, 30])}
    out = sf.flow_summary(_sites_gdf(), s, units="cfs")
    assert out.peak_discharge.iloc[0] == 50
    assert out.peak_time.iloc[0] == pd.Timestamp("2026-08-13T00:05:00Z")


def test_summary_rise_is_relative_to_the_window_start():
    """On a big regulated river the absolute peak says little; the rise above
    where the window started is the storm response."""
    out = sf.flow_summary(_sites_gdf(), {"A": _series([10, 50, 30])}, units="cfs")
    assert out.rise_discharge.iloc[0] == 40
    assert out.rise_ratio.iloc[0] == pytest.approx(5.0)


def test_summary_rise_ratio_is_nan_when_the_river_started_dry():
    out = sf.flow_summary(_sites_gdf(), {"A": _series([0, 5])}, units="cfs")
    assert np.isnan(out.rise_ratio.iloc[0])


def test_summary_units_switch_changes_the_column_read():
    s = {"A": _series([10, 50, 30])}
    si = sf.flow_summary(_sites_gdf(), s, units="si").peak_discharge.iloc[0]
    cfs = sf.flow_summary(_sites_gdf(), s, units="cfs").peak_discharge.iloc[0]
    assert si == pytest.approx(cfs * 0.028316846592)


def test_summary_measures_reporting_cadence():
    out = sf.flow_summary(_sites_gdf(), {"A": _series([1, 2, 3], freq="15min")})
    assert out.report_min.iloc[0] == pytest.approx(15.0)


def test_report_min_needs_two_observations():
    assert np.isnan(sf._report_min(pd.DatetimeIndex([])))
    assert np.isnan(sf._report_min(pd.DatetimeIndex(["2026-08-13"])))


def test_summary_computes_unit_discharge_per_km2():
    out = sf.flow_summary(_sites_gdf(), {"A": _series([10, 50])}, units="cfs")
    assert out.unit_discharge.iloc[0] == pytest.approx(50 / 100.0)


def test_peak_on_the_last_observation_is_flagged():
    """A peak sitting on the final sample is not a peak -- the window cut the
    hydrograph while it was still rising, and everything derived reads low."""
    out = sf.flow_summary(_sites_gdf(), {"A": _series([1, 5, 20])})
    assert bool(out.peak_at_edge.iloc[0]) is True


def test_peak_in_the_middle_is_not_flagged():
    out = sf.flow_summary(_sites_gdf(), {"A": _series([1, 20, 5, 4, 3])})
    assert bool(out.peak_at_edge.iloc[0]) is False


def test_summary_keeps_gauges_that_reported_nothing():
    """A gauge with no data in the window is still a gauge on the map."""
    out = sf.flow_summary(_sites_gdf(), {})
    assert len(out) == 1 and np.isnan(out.peak_discharge.iloc[0])


# --------------------------------------------------------------------------- #
# store round trip
# --------------------------------------------------------------------------- #
def test_event_store_round_trips(tmp_path, monkeypatch):
    calls = {"n": 0}

    def h(url, p):
        calls["n"] += 1
        return _Resp(text=SITE_RDB) if "site" in url else _Resp(_iv_payload())

    _install(monkeypatch, h)
    summary, series = sf.fetch_stream_event(
        (-120, 39, -119, 40), "2026-08-13T00:00Z", "2026-08-13T01:00Z",
        str(tmp_path), "ev", layout="flat")
    assert len(summary) == 2 and "10348000" in series
    back = sf.load_event_series(str(tmp_path), "ev")
    assert "10348000" in back
    assert str(back["10348000"].index.tz) == "UTC"
    pd.testing.assert_series_equal(
        back["10348000"].discharge_cfs, series["10348000"].discharge_cfs,
        check_freq=False)


def test_store_writes_the_geojson_and_per_gauge_csvs(tmp_path, monkeypatch):
    _install(monkeypatch, lambda url, p:
             _Resp(text=SITE_RDB) if "site" in url else _Resp(_iv_payload()))
    sf.fetch_stream_event((-120, 39, -119, 40), "2026-08-13T00:00Z",
                          "2026-08-13T01:00Z", str(tmp_path), "ev",
                          layout="flat")
    assert (tmp_path / "ev_streamgauges.geojson").exists()
    assert (tmp_path / "StreamGaugeData" / "ev_stream_10348000.csv").exists()


def test_stream_store_dir_is_reserved_at_the_event_root():
    """Like RainGaugeData: a data store, not a figure directory."""
    from stormscape import layout
    assert "StreamGaugeData" in layout.RESERVED


def test_load_event_series_on_a_missing_store_is_empty(tmp_path):
    assert sf.load_event_series(str(tmp_path), "nope") == {}


# --------------------------------------------------------------------------- #
# the OGC backend
# --------------------------------------------------------------------------- #
def _ogc_tsm(sites_ends):
    return {"features": [
        {"properties": {"monitoring_location_id": f"USGS-{s}",
                        "parameter_code": "00060", "begin": "1990-01-01",
                        "end": end}} for s, end in sites_ends], "links": []}


def _ogc_ml(sites):
    return {"features": [
        {"properties": {"monitoring_location_number": s,
                        "monitoring_location_name": f"SITE {s}",
                        "drainage_area": 100, "altitude": 4000,
                        "site_type_code": "ST"},
         "geometry": {"coordinates": [-119.5, 39.5]}} for s in sites],
        "links": []}


def test_ogc_sites_require_an_instantaneous_discharge_series(monkeypatch):
    """monitoring-locations alone answers 'is this a stream site', not 'does it
    measure discharge every 15 minutes' -- filtering on it returns a set several
    times too large, so the two sources would disagree."""
    now = pd.Timestamp.now(tz="UTC")
    recent = (now - pd.Timedelta(days=1)).isoformat()

    def h(url, p):
        if "time-series-metadata" in url:
            return _Resp(_ogc_tsm([("10348000", recent)]))
        return _Resp(_ogc_ml(["10348000", "99999999"]))   # extra, no discharge

    _install(monkeypatch, h)
    g = sf.stream_sites((-120, 39, -119, 40), source="ogc")
    assert sorted(g.site_no) == ["10348000"]


def test_ogc_active_filter_drops_discontinued_records(monkeypatch):
    now = pd.Timestamp.now(tz="UTC")
    recent = (now - pd.Timedelta(days=1)).isoformat()

    def h(url, p):
        if "time-series-metadata" in url:
            return _Resp(_ogc_tsm([("10348000", recent),
                                   ("10350400", "1997-01-02T00:00:00Z")]))
        return _Resp(_ogc_ml(["10348000", "10350400"]))

    _install(monkeypatch, h)
    assert sorted(sf.stream_sites((-120, 39, -119, 40),
                                  source="ogc").site_no) == ["10348000"]
    both = sf.stream_sites((-120, 39, -119, 40), source="ogc",
                           active_only=False)
    assert sorted(both.site_no) == ["10348000", "10350400"]


def test_ogc_api_key_is_sent_when_present(monkeypatch):
    seen = []

    def h(url, p):
        seen.append(p)
        return _Resp({"features": [], "links": []})

    _install(monkeypatch, h)
    sf.stream_sites((-120, 39, -119, 40), source="ogc", api_key="SEKRIT")
    assert seen and seen[0].get("api_key") == "SEKRIT"


def test_ogc_reads_the_key_from_the_environment(monkeypatch):
    monkeypatch.setenv("STORMSCAPE_USGS_API_KEY", "envkey")
    assert sf._api_key("ogc", None) == "envkey"
    assert sf._api_key("nwis", None) is None      # legacy needs no key


def test_ogc_pagination_follows_next_links(monkeypatch):
    pages = [
        {"features": [{"properties": {"time": "2026-08-13T00:00:00+00:00",
                                      "value": 1}}],
         "links": [{"rel": "next", "href": "http://x/page2"}]},
        {"features": [{"properties": {"time": "2026-08-13T00:05:00+00:00",
                                      "value": 2}}], "links": []},
    ]
    seq = iter(pages)
    _install(monkeypatch, lambda url, p: _Resp(next(seq)))
    got = list(sf._ogc_pages("http://x", {}, 10, "t", None))
    assert len(got) == 2


# --------------------------------------------------------------------------- #
# figures
# --------------------------------------------------------------------------- #
def test_hydrograph_draws_discharge_and_stage():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from stormscape.plot import hydrograph
    fig, ax = hydrograph(_series([1, 5, 3]), name="TEST")
    assert ax.lines, "no discharge line drawn"
    assert len(fig.axes) == 2, "stage twin axis missing"
    plt.close(fig)


def test_hydrograph_survives_an_all_missing_window():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from stormscape.plot import hydrograph
    df = _series([1, 2])
    df[:] = np.nan
    fig, ax = hydrograph(df, name="TEST")
    assert any("no discharge" in t.get_text() for t in ax.texts)
    plt.close(fig)


def test_hydrograph_adds_a_hyetograph_panel_when_rain_is_given():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from stormscape.plot import hydrograph
    idx = pd.date_range("2026-08-13T00:00Z", periods=3, freq="5min", tz="UTC")
    rain = pd.DataFrame({"i15_mmph": [0.0, 12.0, 3.0]}, index=idx)
    fig, ax = hydrograph(_series([1, 5, 3]), name="TEST", rain=rain)
    # hyetograph + hydrograph + stage twin
    assert len(fig.axes) == 3
    rax = fig.axes[0]
    assert rax.get_ylim()[0] > rax.get_ylim()[1], "hyetograph not inverted"
    plt.close(fig)


def test_atlas_orders_panels_by_peak_discharge():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from stormscape.plot import hydrograph_atlas
    series = {"small": _series([1, 2]), "big": _series([1, 90]),
              "mid": _series([1, 20])}
    fig, axes = hydrograph_atlas(series, ncols=3)
    titles = [a.get_title() for a in axes[0]]
    assert titles[:3] == ["big", "mid", "small"]
    plt.close(fig)


def test_atlas_min_peak_drops_the_irrigation_drains():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from stormscape.plot import hydrograph_atlas
    series = {"drain": _series([0.01, 0.02]), "river": _series([1, 90])}
    fig, axes = hydrograph_atlas(series, min_peak=1.0, units="cfs")
    assert [a.get_title() for a in axes[0] if a.get_title()] == ["river"]
    plt.close(fig)


def test_atlas_on_empty_input_returns_none():
    from stormscape.plot import hydrograph_atlas
    assert hydrograph_atlas({}) == (None, None)


def test_stream_gauge_marker_differs_from_the_rain_gauge():
    """A map can carry both; they must not be mistaken for each other."""
    from stormscape import plot
    assert plot.STREAM_GAUGE_STYLE["marker"] != "o"
    assert plot.STREAM_GAUGE_STYLE["color"].lower() != "yellow"


def test_add_stream_gauges_on_empty_input_is_a_no_op():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    import geopandas as gpd
    from stormscape.plot import add_stream_gauges
    fig, ax = plt.subplots()
    empty = gpd.GeoDataFrame(geometry=gpd.GeoSeries([], crs=4326), crs=4326)
    assert add_stream_gauges(ax, "EPSG:4326", empty) == ([], None)
    plt.close(fig)


def test_add_stream_gauges_colours_by_a_value_column():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from stormscape.plot import add_stream_gauges
    g = _sites_gdf()
    g["peak_discharge"] = [12.0]
    fig, ax = plt.subplots()
    handles, sc = add_stream_gauges(ax, "EPSG:4326", g, value="peak_discharge")
    assert sc is not None and handles
    plt.close(fig)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _run_flow_cli(monkeypatch, argv):
    from stormscape import cli
    seen = {}
    monkeypatch.setattr(cli, "_cmd_flow", lambda args: seen.update(vars(args)))
    cli.main(argv)
    return seen


def test_flow_cli_defaults(monkeypatch):
    got = _run_flow_cli(monkeypatch, ["flow", "--bbox", "-120", "39", "-119",
                                      "40", "--date", "20260813"])
    assert got["source"] == "nwis"          # legacy default, no key needed
    assert got["units"] == "si"             # matches the rest of stormscape
    assert got["include_inactive"] is False
    assert got["no_atlas"] is False


def test_flow_cli_accepts_the_ogc_source_and_units(monkeypatch):
    got = _run_flow_cli(monkeypatch,
                        ["flow", "--bbox", "-120", "39", "-119", "40",
                         "--date", "20260813", "--source", "ogc",
                         "--units", "cfs", "--min-peak", "5", "--detail"])
    assert got["source"] == "ogc" and got["units"] == "cfs"
    assert got["min_peak"] == 5.0 and got["detail"] is True


def test_cmd_flow_only_reads_args_the_parser_defines():
    """Every ``args.<name>`` the command touches must exist in the namespace.

    The same guard that caught `mines` reading an undefined ``args.vmax``:
    patched-out CLI tests cannot see it and the real path needs the network.
    """
    import ast
    import inspect

    from stormscape import cli

    tree = ast.parse(inspect.getsource(cli._cmd_flow))
    used = {n.attr for n in ast.walk(tree)
            if isinstance(n, ast.Attribute) and isinstance(n.value, ast.Name)
            and n.value.id == "args"}
    seen = {}
    real, cli._cmd_flow = cli._cmd_flow, lambda a: seen.update(vars(a))
    try:
        cli.main(["flow", "--bbox", "-120", "39", "-119", "40",
                  "--date", "20260813"])
    finally:
        cli._cmd_flow = real
    missing = sorted(used - set(seen))
    assert not missing, f"_cmd_flow reads undefined args: {missing}"


@pytest.mark.parametrize("cmd", ["map", "run", "nexrad", "zoom", "burn",
                                 "export"])
def test_map_commands_expose_the_stream_gauge_overlay(cmd, monkeypatch):
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
    cli.main([cmd] + base + ["--stream-gauges"])
    assert seen["stream_gauges"] is True


def test_flow_does_not_advertise_an_inert_overlay_flag():
    """`flow` always draws its own gauges, so --stream-gauges would do nothing
    there; a flag that silently does nothing is worse than no flag."""
    from stormscape import cli
    import argparse
    import contextlib
    import io
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), pytest.raises(SystemExit):
        cli.main(["flow", "--help"])
    assert "--stream-gauges" not in buf.getvalue()


def test_stream_kwargs_empty_when_off():
    from argparse import Namespace

    from stormscape import cli
    assert cli._stream_kwargs(Namespace(stream_gauges=False)) == {}
    assert cli._stream_kwargs(Namespace()) == {}
