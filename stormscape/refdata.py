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

import sys
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


def _more_pages(fc) -> bool:
    """Whether an ArcGIS response says it truncated the result set.

    The flag lives in **two different places** depending on the server: a
    MapServer (NHD, TIGER, GNIS) puts ``exceededTransferLimit`` at the top
    level of the GeoJSON, while a hosted FeatureServer (USMIN) puts it only
    under ``properties``. Checking one location silently truncates the other's
    results at the page size -- USMIN returns exactly 2,000 of 3,134 features
    over the Hidden Valley AOI, with no error to notice.
    """
    return bool(fc.get("exceededTransferLimit")
                or (fc.get("properties") or {}).get("exceededTransferLimit"))


def arcgis_query(service, layer, bounds, out_fields="*", where="1=1", page=2000,
                 timeout=60, token=None, paginates=True, what=None):
    """Page an ArcGIS REST layer over a (W,S,E,N) bbox -> GeoDataFrame(4326).

    ``token`` authenticates against a gated service (never logged).
    ``paginates=False`` suits a server that rejects ``resultOffset`` outright
    (the NBMG USMIN mirror answers "Pagination is not supported"), in which
    case a single page is fetched. ``what`` labels this layer in any failure
    message.
    """
    geom = f"{bounds[0]},{bounds[1]},{bounds[2]},{bounds[3]}"
    label = what or f"{service}/{layer}"
    feats, offset = [], 0
    while True:
        params = dict(geometry=geom, geometryType="esriGeometryEnvelope",
                      inSR=4326, outSR=4326, spatialRel="esriSpatialRelIntersects",
                      where=where, outFields=out_fields, returnGeometry="true",
                      f="geojson")
        if paginates:
            params.update(resultRecordCount=page, resultOffset=offset)
        if token:
            params["token"] = token
        try:
            r = requests.get(f"{service}/{layer}/query", params=params,
                             timeout=timeout)
            fc = r.json()
        except Exception as e:                         # noqa: BLE001
            _fail(label, repr(e)[:120])
            break
        # An ArcGIS error is HTTP 200 with an {"error": ...} body, so it would
        # otherwise read as "no features here" -- indistinguishable from an AOI
        # that genuinely has none.
        if isinstance(fc, dict) and "error" in fc:
            err = fc["error"]
            _fail(label, f"{err.get('code')} {err.get('message')}")
            break
        batch = fc.get("features", [])
        feats.extend(batch)
        if not paginates or not batch or not _more_pages(fc):
            break
        offset += page
    if not feats:
        return gpd.GeoDataFrame(geometry=[], crs=4326)
    return gpd.GeoDataFrame.from_features(feats, crs=4326)


def _fail(label, detail):
    """Report a query failure loudly.

    Printed as well as warned because this module installs a global
    ``warnings.filterwarnings("ignore")``, which would otherwise swallow the
    warning and leave a dead service looking like an empty AOI.
    """
    msg = f"query failed ({label}): {detail}"
    warnings.warn(msg)
    print(f"warning: {msg}", file=sys.stderr)


#: back-compat alias -- this was private before mines.py needed it.
_query = arcgis_query


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
