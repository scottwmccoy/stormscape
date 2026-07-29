"""Vector reference layers (streams, roads, place names) for an AOI.

Fetches authoritative, AOI-scoped vectors from public ArcGIS REST services so
they can be drawn as crisp labelled overlays on a hillshade -- the same
underlying data the USGS pfdf tutorials use, but queried per bounding box
(small, fast) instead of downloading whole hydrologic-unit bundles:

  * streams  -- USGS **NHDPlus HR** network flowlines (``gnis_name``); this is
    the National Hydrography Dataset, the source of named creeks/rivers;
  * roads    -- US **Census TIGER/Line** primary & secondary roads (``NAME``);
  * places   -- USGS **GNIS** populated places (``gaz_name``).

Each function returns a GeoDataFrame in EPSG:4326 with a tidy ``name`` column
(plus a ``kind`` for roads/places). Service hiccups degrade gracefully to an
empty GeoDataFrame so a flaky endpoint never kills a figure.
"""

from __future__ import annotations

import warnings

import geopandas as gpd
import requests

from .aoi import load_aoi

warnings.filterwarnings("ignore")

NHDPLUS_HR = ("https://hydro.nationalmap.gov/arcgis/rest/services/"
              "NHDPlus_HR/MapServer")
TIGER_TRANS = ("https://tigerweb.geo.census.gov/arcgis/rest/services/"
               "TIGERweb/Transportation/MapServer")
GNIS = ("https://carto.nationalmap.gov/arcgis/rest/services/"
        "geonames/MapServer")

# NHD FCodes: StreamRiver 460**, ArtificialPath 55800, Connector 33400,
# CanalDitch 336**, Pipeline 428**. "Watercourse" excludes ditches/pipelines.
_WATERCOURSE = {46000, 46003, 46006, 46007, 55800, 33400}


def _query(service, layer, bounds, out_fields="*", where="1=1", page=2000,
           timeout=60):
    """Page an ArcGIS REST layer over a (W,S,E,N) bbox -> GeoDataFrame(4326)."""
    geom = f"{bounds[0]},{bounds[1]},{bounds[2]},{bounds[3]}"
    feats, offset = [], 0
    while True:
        params = dict(geometry=geom, geometryType="esriGeometryEnvelope",
                      inSR=4326, outSR=4326, spatialRel="esriSpatialRelIntersects",
                      where=where, outFields=out_fields, returnGeometry="true",
                      f="geojson", resultRecordCount=page, resultOffset=offset)
        try:
            r = requests.get(f"{service}/{layer}/query", params=params,
                             timeout=timeout)
            fc = r.json()
        except Exception as e:                         # noqa: BLE001
            warnings.warn(f"reference query failed ({service}/{layer}): "
                          f"{repr(e)[:100]}")
            break
        batch = fc.get("features", [])
        feats.extend(batch)
        if len(batch) < page or not fc.get("exceededTransferLimit"):
            break
        offset += page
    if not feats:
        return gpd.GeoDataFrame(geometry=[], crs=4326)
    return gpd.GeoDataFrame.from_features(feats, crs=4326)


def streams(aoi, named_only=False, watercourse_only=True):
    """USGS NHDPlus-HR flowlines for the AOI (``name`` = gnis_name).

    ``watercourse_only`` drops canals/ditches/pipelines (keeps streams,
    rivers, connectors). ``named_only`` keeps only flowlines that carry a
    GNIS name.
    """
    bounds, _ = load_aoi(aoi)
    g = _query(NHDPLUS_HR, 3, bounds, out_fields="gnis_name,fcode")
    if not len(g):
        return g
    g["name"] = g.get("gnis_name")
    if watercourse_only and "fcode" in g:
        g = g[g.fcode.fillna(0).astype(int).isin(_WATERCOURSE)]
    if named_only:
        g = g[g.name.notna() & (g.name.astype(str).str.len() > 0)]
    return g.reset_index(drop=True)


def roads(aoi, local=False):
    """TIGER primary + secondary (+ optional local) roads (``name`` = NAME)."""
    bounds, _ = load_aoi(aoi)
    layers = [(2, "primary"), (6, "secondary")]
    if local:
        layers.append((8, "local"))
    parts = []
    for lid, kind in layers:
        g = _query(TIGER_TRANS, lid, bounds, out_fields="NAME")
        if len(g):
            g["kind"] = kind
            parts.append(g)
    if not parts:
        return gpd.GeoDataFrame(geometry=[], crs=4326)
    out = gpd.GeoDataFrame(__import__("pandas").concat(parts,
                                                       ignore_index=True),
                           crs=4326)
    out["name"] = (out.get("NAME").astype(str)
                   .str.replace(r"\s+", " ", regex=True).str.strip())
    return out


def places(aoi, populated=True, incorporated=True):
    """GNIS place points for the AOI (``name`` = gaz_name, ``kind``)."""
    bounds, _ = load_aoi(aoi)
    layers = []
    if incorporated:
        layers.append((1, "incorporated"))
    if populated:
        layers.append((3, "populated"))
    parts = []
    for lid, kind in layers:
        g = _query(GNIS, lid, bounds, out_fields="gaz_name,gaz_featureclass")
        if len(g):
            g["kind"] = kind
            g["name"] = g.get("gaz_name")
            parts.append(g)
    if not parts:
        return gpd.GeoDataFrame(geometry=[], crs=4326)
    out = gpd.GeoDataFrame(__import__("pandas").concat(parts,
                                                       ignore_index=True),
                           crs=4326)
    return out.drop_duplicates("name").reset_index(drop=True)


def reference_layers(aoi, local_roads=False):
    """Fetch all three reference layers for an AOI -> dict of GeoDataFrames."""
    return dict(streams=streams(aoi), roads=roads(aoi, local=local_roads),
                places=places(aoi))
