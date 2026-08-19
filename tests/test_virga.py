"""Virga-risk mask: classification truth table + raster front-end."""
from __future__ import annotations

import numpy as np
import pytest

from stormscape import virga
from stormscape.virga import SUPPORTED, VIRGA_RISK, UNDERREAD, NODATA, classify


def test_classify_truth_table():
    mosaic  = np.array([[34.0,  3.3, 20.0, 1.0,  0.0, np.nan]])
    support = np.array([[ 1.2, 22.1, 18.0, 1.0, 30.0, 15.0]])
    out = classify(mosaic, support)
    #        Paiute      WValley    agree      drizzle   zero-vs-30  nan
    assert out[0, 0] == VIRGA_RISK          # 34 vs 1.2: aloft only
    assert out[0, 1] == UNDERREAD           # 3.3 vs 22.1: low tilt saw more
    assert out[0, 2] == SUPPORTED           # within 3x
    assert out[0, 3] == NODATA              # both below min_mmph
    assert out[0, 4] == UNDERREAD           # hard zero vs heavy rain
    assert out[0, 5] == NODATA              # missing data


def test_drizzle_agreement_is_not_evidence():
    """Both fields at 2 mm/h agree perfectly -- and must stay NODATA."""
    out = classify(np.full((3, 3), 2.0), np.full((3, 3), 2.0))
    assert (out == NODATA).all()


def test_one_sided_intensity_is_assessed():
    """min_mmph is an OR: 34 mm/h vs 0.5 must be flagged, not skipped."""
    out = classify(np.array([[34.0]]), np.array([[0.5]]))
    assert out[0, 0] == VIRGA_RISK


def test_ratio_threshold_is_inclusive_and_tunable():
    m, s = np.array([[30.0]]), np.array([[10.0]])
    assert classify(m, s, ratio=3.0)[0, 0] == VIRGA_RISK      # exactly 3x
    assert classify(m, s, ratio=3.1)[0, 0] == SUPPORTED


def test_shape_mismatch_raises():
    with pytest.raises(ValueError, match="regrid"):
        classify(np.zeros((2, 2)), np.zeros((3, 3)))


def _write(path, arr, transform, crs="EPSG:4326"):
    import rasterio
    with rasterio.open(path, "w", driver="GTiff", height=arr.shape[0],
                       width=arr.shape[1], count=1, dtype="float32",
                       crs=crs, transform=transform, nodata=np.nan) as ds:
        ds.write(arr.astype("float32"), 1)


def test_virga_mask_regrids_and_writes(tmp_path):
    """Different grids in, classified uint8 + ratio tifs out, layout-sorted."""
    import rasterio
    tr_s = rasterio.transform.from_origin(-119.5, 39.8, 0.005, 0.005)   # fine
    tr_m = rasterio.transform.from_origin(-119.5, 39.8, 0.01, 0.01)     # coarse
    sup = np.full((20, 20), 1.0);  sup[5:10, 5:10] = 40.0               # real core
    mos = np.full((10, 10), 30.0)                                       # aloft blanket
    _write(tmp_path / "sup.tif", sup, tr_s)
    _write(tmp_path / "mos.tif", mos, tr_m)
    res = virga.virga_mask(str(tmp_path / "mos.tif"), str(tmp_path / "sup.tif"),
                           str(tmp_path), "ev")
    with rasterio.open(res["virgarisk_tif"]) as ds:
        cls = ds.read(1)
        assert ds.tags()["CLASSES"].startswith("0=supported")
    assert cls.shape == sup.shape                       # support grid wins
    assert cls[7, 7] == SUPPORTED                       # 30 vs 40: agree
    assert cls[15, 15] == VIRGA_RISK                    # 30 vs 1: aloft only
    assert "rasters" in res["virgarisk_tif"]            # sorted layout
    assert res["counts"]["virga_risk"] > 0
    assert res["percent"]["virga_risk"] + res["percent"]["supported"] \
        + res["percent"]["underread"] == pytest.approx(100.0, abs=0.2)


def test_near_radar_exclusion(tmp_path):
    """Cells inside exclude_km of the radar become NODATA even if extreme."""
    import rasterio
    tr = rasterio.transform.from_origin(-119.5, 39.8, 0.005, 0.005)
    sup = np.full((20, 20), 1.0)
    mos = np.full((20, 20), 50.0)
    _write(tmp_path / "sup.tif", sup, tr)
    _write(tmp_path / "mos.tif", mos, tr)
    centre = (-119.5 + 10 * 0.005, 39.8 - 10 * 0.005)
    res = virga.virga_mask(str(tmp_path / "mos.tif"), str(tmp_path / "sup.tif"),
                           str(tmp_path), "ev", radar_lonlat=centre,
                           exclude_km=1.0)
    with rasterio.open(res["virgarisk_tif"]) as ds:
        cls = ds.read(1)
    assert cls[10, 10] == NODATA                        # at the radar
    assert cls[0, 0] == VIRGA_RISK                      # far corner still flagged


def test_cli_smoke(tmp_path, capsys):
    import rasterio
    from stormscape.cli import main
    tr = rasterio.transform.from_origin(-119.5, 39.8, 0.005, 0.005)
    _write(tmp_path / "sup.tif", np.full((8, 8), 1.0), tr)
    _write(tmp_path / "mos.tif", np.full((8, 8), 30.0), tr)
    main(["virga", "--mosaic", str(tmp_path / "mos.tif"),
          "--support", str(tmp_path / "sup.tif"),
          "--out-dir", str(tmp_path), "--key", "t"])
    out = capsys.readouterr().out
    assert "virga_risk" in out and "wrote" in out
