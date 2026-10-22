"""Synoptic access routing: the academic tier toggles, stormscape must not care.

Synoptic reports most failures inside an HTTP-200 body. A token cut back to the
Basic Precipitation Service gets ``RESPONSE_CODE 403`` from stations/metadata,
and stormscape used to read that as zero stations. It also answers an
overlapping call (one at a time per account) with the same 403 text.
"""
import json
from datetime import datetime

import pytest

from stormscape import gauges

DENIED = {"SUMMARY": {"RESPONSE_CODE": 403, "RESPONSE_MESSAGE": "Account associated "
                      "with this token does not have access to the requested service."}}
TOO_MANY = {"SUMMARY": {"RESPONSE_CODE": -1, "RESPONSE_MESSAGE": "Querying too many "
                        "station hours. Please reduce the request size."}}
OLD = {"SUMMARY": {"RESPONSE_CODE": -1, "RESPONSE_MESSAGE": "Data outside the "
                   "allowed period of record for this account."}}


def _station(stid, lon=-120.04417, lat=39.41417):
    return {"STID": stid, "NAME": stid, "STATE": "CA", "ELEVATION": "6800",
            "MNET_ID": "2", "STATUS": "ACTIVE", "TIMEZONE": "America/Los_Angeles",
            "LATITUDE": str(lat), "LONGITUDE": str(lon), "OBSERVATIONS": {}}


def _ok(*stids):
    return {"SUMMARY": {"RESPONSE_CODE": 1, "RESPONSE_MESSAGE": "OK"},
            "STATION": [_station(s) for s in stids]}


class _R:
    def __init__(self, js):
        self.status_code, self._js = 200, js

    def json(self):
        return self._js


@pytest.fixture(autouse=True)
def _fresh(monkeypatch):
    gauges.reset_access()
    monkeypatch.setattr(gauges.time, "sleep", lambda s: None)
    yield
    gauges.reset_access()


def _serve(monkeypatch, responses):
    calls = []

    def fake_get(url, params=None, timeout=None):
        calls.append((url, dict(params or {})))
        return _R(responses[len(calls) - 1])
    monkeypatch.setattr(gauges.requests, "get", fake_get)
    return calls


BOX = (-120.1, 39.3, -119.9, 39.5)
T0, T1 = datetime(2026, 10, 21, 20), datetime(2026, 10, 22, 8)


def test_precip_only_token_still_finds_stations(monkeypatch):
    calls = _serve(monkeypatch, [DENIED, DENIED, _ok("XOUC1")])
    with pytest.warns(UserWarning, match="precipitation service"):
        st = gauges.get_stations(BOX, token="t", starttime=T0, endtime=T1)
    assert list(st.station_id) == ["XOUC1"]
    assert st.geometry.iloc[0].x == pytest.approx(-120.04417)
    url, params = calls[2]
    assert url == gauges.TS_URL and params["precip"] == "1"
    assert (params["start"], params["end"]) == ("202610212000", "202610220800")
    rep = gauges.access_report()
    assert rep["metadata"]["status"] == "denied"
    assert rep["timeseries_precip"]["status"] == "ok"


def test_one_403_is_a_busy_account_not_a_denial(monkeypatch):
    calls = _serve(monkeypatch, [DENIED, _ok("KRNO")])
    st = gauges.get_stations(BOX, token="t")
    assert list(st.station_id) == ["KRNO"] and len(calls) == 2
    assert all(u == gauges.META_URL for u, _ in calls)
    assert gauges.access_report()["metadata"]["status"] == "ok"


def test_denied_route_is_skipped_for_the_rest_of_the_run(monkeypatch):
    calls = _serve(monkeypatch, [DENIED, DENIED, _ok("A"), _ok("B")])
    with pytest.warns(UserWarning):
        gauges.get_stations(BOX, token="t", starttime=T0, endtime=T1)
    with pytest.warns(UserWarning):
        st = gauges.get_stations(BOX, token="t", starttime=T0, endtime=T1)
    assert list(st.station_id) == ["B"]
    assert [u for u, _ in calls].count(gauges.META_URL) == 2     # never asked again


def test_both_routes_denied_raises_access_error(monkeypatch):
    _serve(monkeypatch, [DENIED] * 4)
    with pytest.raises(gauges.SynopticAccessError, match="timeseries_precip"):
        gauges.get_stations(BOX, token="t", starttime=T0, endtime=T1)


def test_discovery_splits_an_oversized_window(monkeypatch):
    calls = _serve(monkeypatch, [DENIED, DENIED, TOO_MANY, _ok("A", "B"), _ok("B", "C")])
    with pytest.warns(UserWarning):
        st = gauges.get_stations(BOX, token="t", starttime=T0, endtime=T1)
    assert sorted(st.station_id) == ["A", "B", "C"]
    assert len(calls) == 5


def test_refused_window_is_named_not_swallowed(monkeypatch):
    _serve(monkeypatch, [OLD])
    with pytest.raises(gauges.SynopticError, match="2026-10-21 20:00 to 2026-10-22 08:00"):
        gauges.get_rainfall(["XOUC1"], T0, T1, token="t")
    assert gauges.access_report() == {}          # not mistaken for a service denial


def test_access_record_is_written(monkeypatch, tmp_path):
    _serve(monkeypatch, [DENIED, DENIED, _ok("XOUC1")])
    with pytest.warns(UserWarning):
        gauges.get_stations(BOX, token="t", starttime=T0, endtime=T1)
    path = gauges.write_access_report(str(tmp_path), "ev", layout="flat")
    rec = json.loads(open(path).read())
    assert rec["metadata"]["status"] == "denied" and "checked_utc" in rec["metadata"]
