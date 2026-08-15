"""Abandoned mine features for an AOI -- USGS **USMIN** topographic mine symbols.

Post-fire debris flows run out of steep, freshly burned catchments, and in the
Great Basin those catchments are full of historic mining. Mine dumps and
tailings are loose, often fine-grained, sometimes contaminated material sitting
on or near the channel network -- exactly the stuff a convective cell mobilises.
This module puts those features on the same map as the rain.

Why this source and not an AML hazard database
----------------------------------------------
There is **no public point-level abandoned-mine hazard inventory**, and that is
policy rather than oversight. USGS Fact Sheet 2025-3003 (2025) describes the
national abandoned-mine-feature database being assembled under USMIN and states
that it "will not publish specific location information of any abandoned mine
workings, and the detailed national abandoned mine feature database will not be
publicly available" -- because the locations could be used to enter hazardous
workings or vandalise historic structures. Only aggregated derivatives (counts
per county or watershed) are planned for release. Nevada matches the pattern:
the Division of Minerals' operational AML layers (``NVPoints``, ``NVSites``,
``InternFieldDataCaptureLayers``) exist on their ArcGIS Online organisation but
every one answers ``499 Token Required`` to an anonymous request.

What *is* public is USMIN's **Prospect- and Mine-Related Features from USGS 7.5-
and 15-minute topographic quadrangle maps** -- mine symbols digitised from
already-published historical topo sheets, so it reveals nothing the printed maps
did not. It is national, uniform, and served live. Nevada alone carries 121,193
point features; every western state is covered.

Read it for what it is: **a historical map compilation, not a hazard
inventory.** A feature is where a topographer drew a symbol on a sheet dated
1950-1994 (``topo_date``). There is no hazard ranking, no securing status, no
confirmation the feature still exists, and no guarantee that a mine opened after
the last map revision appears at all. For the authoritative soil-and-hazard
picture you need the state programme's own database -- see ``SOURCES`` below,
which is a registry precisely so a credentialed NDOM feed can be dropped in
without touching call sites.

Feature groups
--------------
USMIN's vocabulary is 55 feature types nationally, most of which are noise for
this purpose -- two-thirds of the features over a typical Nevada AOI are
prospect pits, shallow exploration diggings that would bury a rainfall field in
dots. :data:`GROUPS` collapses the vocabulary into six classes, and
:data:`DEFAULT_KINDS` selects the two that matter for storm work: ``waste``
(dumps, tailings, slag, process ponds) and ``openings`` (adits and shafts).
Pass ``kinds="all"`` for the raw firehose, or name individual feature types.

**The waste is in the polygon layer.** USMIN splits features between a point
layer and a polygon layer according to how the topographer drew the symbol, and
dumps and tailings were nearly always drawn as an extent -- 14,815 polygons
against 413 points nationally. ``geometry="both"`` is the default for that
reason; ask for points only and a query for mine waste comes back all but
empty, which looks exactly like an AOI that has none.

Density instead of points
-------------------------
Plotting every feature over a mining district is unreadable. :func:`density_grid`
bins features onto an equal-area grid (EPSG:5070) so a symbol can be sized by
the count per square kilometre -- the district reads as intensity rather than a
smear of markers. See :func:`stormscape.plot.add_mines`.
"""

from __future__ import annotations

import os
import sys

import geopandas as gpd
import pandas as pd

from .aoi import load_aoi
from .refdata import arcgis_query

USMIN = ("https://energy.usgs.gov/arcgis/rest/services/Hosted/"
         "USMin_Prospect_and_mine_related_map_features/FeatureServer")

#: NBMG's mirror of USMIN (version 4.0, Nevada). Older ArcGIS Server that
#: rejects pagination outright, so it can only return one page -- the national
#: service above is newer and pages properly. Kept for provenance.
USMIN_NBMG = ("https://gisweb.unr.edu/nbmg/rest/services/MineralsAndEnergy/"
              "USMIN/MapServer")

#: Nevada Division of Minerals AML layers. Token-gated: anonymous requests get
#: ``499 Token Required``. Registered so that access, once granted, is a
#: ``--source ndom --token ...`` away rather than a code change. The field
#: candidates are unverified guesses -- ``_pick`` tolerates absent columns, so a
#: wrong guess yields an empty column rather than a silently mislabelled one.
NDOM = "https://services.arcgis.com/CXYUMoYknZtf5Qr3/arcgis/rest/services"

#: logical field -> candidate source column names, first match wins. USMIN's
#: national service is lower-case; the NBMG mirror is ``Ftr_Type`` -- carrying
#: both means one spec serves either.
_USMIN_FIELDS = dict(
    ftr_type=("ftr_type", "Ftr_Type"),
    name=("ftr_name", "Ftr_Name"),
    state=("state", "State"),
    county=("county", "County"),
    topo_name=("topo_name", "Topo_Name"),
    topo_date=("topo_date", "Topo_Date"),
)

SOURCES = {
    "usmin": dict(
        label="USGS USMIN topographic mine symbols",
        service=USMIN,
        layers={"points": 17, "areas": 18},
        fields=_USMIN_FIELDS,
        token_env=None,
        public=True,
        note="digitised from published USGS topo sheets; historical, not a "
             "hazard inventory",
    ),
    "usmin_nbmg": dict(
        label="USMIN v4.0 mirror (NBMG)",
        service=USMIN_NBMG,
        layers={"points": 0},
        fields=_USMIN_FIELDS,
        token_env=None,
        public=True,
        paginates=False,          # server rejects resultOffset entirely
        note="Nevada only; single page per query -- prefer 'usmin'",
    ),
    "ndom": dict(
        label="Nevada Division of Minerals AML (token required)",
        service=NDOM,
        layers={"points": "NVPoints/FeatureServer/0",
                "areas": "NVSites/FeatureServer/0"},
        fields=dict(ftr_type=("Type", "TYPE", "FeatureType", "HazardType"),
                    name=("Name", "NAME", "SiteName"),
                    state=("State", "STATE"),
                    county=("County", "COUNTY"),
                    topo_name=("Quad", "QUAD", "Topo_Name"),
                    topo_date=("Date", "DATE", "SurveyDate")),
        token_env="STORMSCAPE_NDOM_TOKEN",
        public=False,
        note="operational AML database; requires credentials from NDOM",
    ),
}

#: group -> one-line description, in the order a legend should list them.
GROUPS = {
    "waste": "mine waste (dumps, tailings, slag, process ponds)",
    "openings": "underground access (adits, shafts)",
    "surface": "open workings (pits, quarries, placer/strip mines)",
    "aggregate": "sand, gravel, borrow and industrial-mineral pits",
    "prospect": "exploration (prospect pits, diggings, trenches)",
    "other": "mill sites, tipples, unrecognised types",
}

#: the two groups that matter for storm/debris-flow context: erodible,
#: potentially contaminated material and open ground-access hazards.
DEFAULT_KINDS = ("waste", "openings")

# Exact feature-type -> group. Built from the full national vocabulary of both
# USMIN layers (47 point types + 48 polygon types, 55 distinct).
_EXACT = {
    # waste -- loose material and impoundments, the storm-relevant class
    "mine dump": "waste", "slag pile": "waste",
    "ore stockpile/storage": "waste",
    "leach pond": "waste", "settling pond": "waste",
    # openings -- ways into the ground
    "adit": "openings", "mine shaft": "openings", "air shaft": "openings",
    "glory hole": "openings",
    # surface -- extraction workings and disturbed ground
    "open pit mine": "surface", "open pit mine or quarry": "surface",
    "strip mine": "surface", "placer mine": "surface",
    "hydraulic mine": "surface", "mine": "surface", "coal mine": "surface",
    "uranium mine": "surface", "silica mine": "surface",
    "iron pit": "surface", "lignite pit": "surface",
    "evaporation pond": "surface", "salt evaporator": "surface",
    # aggregate -- construction and industrial minerals, usually shallow
    "borrow pit": "aggregate", "gravel pit": "aggregate",
    "sand pit": "aggregate", "sand and gravel pit": "aggregate",
    "gravel/borrow pit - undifferentiated": "aggregate",
    "cinder pit": "aggregate", "scoria pit": "aggregate",
    "clay pit": "aggregate", "marl pit": "aggregate",
    "caliche pit": "aggregate", "chert pit": "aggregate",
    "bentonite pit": "aggregate", "shale pit": "aggregate",
    "shell pit": "aggregate",
    # prospect -- exploration scratches; the bulk of the record
    "prospect pit": "prospect", "diggings": "prospect", "trench": "prospect",
    # other -- structures
    "mill site": "other", "tipple": "other",
}

# Prefix families, checked after the exact map. These exist so that a subtype
# USMIN adds later ("Tailings - Something") lands in the right group instead of
# silently falling through to "other".
_PREFIX = (
    ("tailings", "waste"),
    ("quarry", "surface"),
    ("disturbed surface", "surface"),
)


def register_source(name, **spec):
    """Add or replace a source in :data:`SOURCES`.

    Lets a credentialed feed (e.g. NDOM once access is granted) be plugged in
    at runtime without editing this module::

        mines.register_source("ndom", service=..., layers={...},
                              fields={...}, token_env="STORMSCAPE_NDOM_TOKEN")
    """
    SOURCES[name] = dict(spec)
    return SOURCES[name]


def group_of(ftr_type) -> str:
    """Group name for a USMIN feature type (``"other"`` if unrecognised)."""
    t = str(ftr_type or "").strip().lower()
    if not t:
        return "other"
    if t in _EXACT:
        return _EXACT[t]
    for prefix, grp in _PREFIX:
        if t.startswith(prefix):
            return grp
    return "other"


def _resolve_kinds(kinds):
    """Normalise ``kinds`` -> (set of group names, set of explicit types).

    Accepts group names (``"waste"``), explicit feature types
    (``"Mine Dump"``), ``None``/``"all"`` for everything, and any iterable or
    comma-separated string mixing them.
    """
    if kinds is None:
        return set(GROUPS), set()
    if isinstance(kinds, str):
        kinds = [k for k in kinds.replace(",", " ").split() if k]
    kinds = [str(k).strip() for k in kinds if str(k).strip()]
    if not kinds or any(k.lower() == "all" for k in kinds):
        return set(GROUPS), set()
    groups, types = set(), set()
    for k in kinds:
        if k.lower() in GROUPS:
            groups.add(k.lower())
        else:
            types.add(k.lower())
    return groups, types


def _types_for(groups):
    """Known feature types belonging to ``groups`` (exact matches only)."""
    return {t for t, g in _EXACT.items() if g in groups}


def _where_for(groups, types, field):
    """Server-side ``where`` clause for the wanted kinds, or ``None``.

    Pushing the type filter into the query keeps a mining district from coming
    down the wire in full -- over the Hidden Valley AOI ``waste`` is ~90 of
    3,134 features. Returns ``None`` when the selection cannot be expressed
    server-side (``other`` is a catch-all for *unrecognised* types, which by
    definition cannot be enumerated), leaving the caller to filter locally.
    """
    if "other" in groups or groups == set(GROUPS):
        return None
    wanted = _types_for(groups) | set(types)
    if not wanted:
        return None
    quoted = ", ".join("'" + t.replace("'", "''") + "'" for t in sorted(wanted))
    clause = f"LOWER({field}) IN ({quoted})"
    # prefix families: keep future subtypes without another round trip
    for prefix, grp in _PREFIX:
        if grp in groups:
            clause += f" OR LOWER({field}) LIKE '{prefix}%'"
    return clause


def _pick(df, candidates):
    """First present column among ``candidates`` -> Series, else all-NA."""
    for c in candidates:
        if c in df.columns:
            return df[c]
    return pd.Series([pd.NA] * len(df), index=df.index)


def _blank_to_na(s):
    """USMIN writes unnamed features as ``''``, not null -- 2,990 of 3,134 over
    the Hidden Valley AOI. Left alone, ``.notna()`` calls every feature named."""
    out = s.astype("object").where(s.notna(), None)
    return out.map(lambda v: None if v is None or str(v).strip() == ""
                   else str(v).strip())


def _resolve_token(spec, token):
    """Explicit token wins, else the source's env var. Never logged."""
    if token:
        return token
    env = spec.get("token_env")
    return os.environ.get(env) if env else None


def _empty():
    return gpd.GeoDataFrame(
        {c: pd.Series(dtype="object") for c in
         ("name", "ftr_type", "group", "state", "county", "topo_name",
          "topo_date", "geom_kind", "source")},
        geometry=gpd.GeoSeries([], crs=4326), crs=4326)


def mine_features(aoi, kinds=DEFAULT_KINDS, geometry="both", source="usmin",
                  named_only=False, token=None, pad_deg=0.0, timeout=60):
    """Mine features intersecting ``aoi`` -> GeoDataFrame (EPSG:4326).

    Parameters
    ----------
    aoi
        Anything :func:`stormscape.aoi.load_aoi` accepts -- a ``(W,S,E,N)``
        bbox, a vector path, or a shapely geometry.
    kinds
        Which features to keep: group names from :data:`GROUPS`, explicit
        feature types (``"Mine Dump"``), or ``"all"``. Defaults to
        :data:`DEFAULT_KINDS` (``waste`` + ``openings``). The filter is pushed
        into the service query where it can be expressed, so narrow selections
        are cheap.
    geometry
        ``"both"`` (default), ``"points"``, or ``"areas"``. **Default to both
        or you will miss the waste.** USMIN splits its record across a point
        layer and a polygon layer by how the symbol was drawn on the sheet, and
        dumps and tailings were almost always drawn as mapped extents: 14,815
        of them are polygons against 413 points nationally, 778 against 4 in
        Nevada. Asking for ``kinds="waste", geometry="points"`` is therefore a
        near-empty query that looks like a quiet AOI. Areas carry real
        footprints to intersect with a catchment; points are cheaper to count.
    source
        Key into :data:`SOURCES`. ``"usmin"`` is the public national service.
    named_only
        Keep only features carrying a name on the topo sheet.
    token
        ArcGIS token for a gated source; falls back to the source's
        ``token_env``. Passed straight to the service and never logged.

    Returns a frame with ``name``, ``ftr_type``, ``group``, ``state``,
    ``county``, ``topo_name``, ``topo_date``, ``geom_kind`` and ``source``.
    A dead service degrades to an empty frame with a message on stderr.
    """
    if source not in SOURCES:
        raise ValueError(f"unknown mine source {source!r}; "
                         f"have {sorted(SOURCES)}")
    spec = SOURCES[source]
    if geometry not in ("points", "areas", "both"):
        raise ValueError("geometry must be 'points', 'areas' or 'both'")
    if not spec.get("public", True) and not _resolve_token(spec, token):
        raise PermissionError(
            f"source {source!r} ({spec['label']}) needs a token: pass "
            f"token=... or set ${spec.get('token_env')}. "
            f"{spec.get('note', '')}".strip())

    bounds, _ = load_aoi(aoi, pad_deg=pad_deg)
    groups, types = _resolve_kinds(kinds)
    fields = spec["fields"]
    tok = _resolve_token(spec, token)

    wants = (["points", "areas"] if geometry == "both" else [geometry])
    parts = []
    for kind in wants:
        layer = spec["layers"].get(kind)
        if layer is None:
            continue
        where = _where_for(groups, types, fields["ftr_type"][0]) or "1=1"
        g = arcgis_query(spec["service"], layer, bounds, where=where,
                         token=tok, timeout=timeout,
                         paginates=spec.get("paginates", True),
                         what=f"{source} {kind}")
        if not len(g):
            continue
        out = gpd.GeoDataFrame(
            dict(name=_blank_to_na(_pick(g, fields["name"])),
                 ftr_type=_pick(g, fields["ftr_type"]),
                 state=_pick(g, fields["state"]),
                 county=_pick(g, fields["county"]),
                 topo_name=_pick(g, fields["topo_name"]),
                 topo_date=_pick(g, fields["topo_date"])),
            geometry=g.geometry, crs=4326)
        out["geom_kind"] = "point" if kind == "points" else "area"
        parts.append(out)

    if not parts:
        return _empty()
    out = gpd.GeoDataFrame(pd.concat(parts, ignore_index=True), crs=4326)
    out["group"] = out.ftr_type.map(group_of)
    out["source"] = source
    # the server-side clause is an optimisation, not the contract -- re-apply
    # locally so an un-pushable selection ('other') and a source whose `where`
    # was ignored both still honour `kinds`
    if groups != set(GROUPS) or types:
        keep = out.group.isin(groups)
        if types:
            keep |= out.ftr_type.astype(str).str.lower().isin(types)
        out = out[keep]
    if named_only:
        out = out[out.name.notna()]
    cols = ["name", "ftr_type", "group", "state", "county", "topo_name",
            "topo_date", "geom_kind", "source", "geometry"]
    return out[cols].reset_index(drop=True)


def _as_points(gdf):
    """Representative points for any geometry type.

    Polygon layers (mine dumps, tailings ponds) need a point to be counted or
    plotted; ``representative_point`` is guaranteed inside the polygon, unlike
    a centroid of a crescent-shaped dump.
    """
    geom = gdf.geometry
    is_pt = geom.geom_type == "Point"
    if is_pt.all():
        return geom
    return geom.where(is_pt, geom.representative_point())


def density_grid(mines, cell_km=1.0, crs="EPSG:5070", groups=None):
    """Count mine features per equal-area grid cell.

    Returns a GeoDataFrame (EPSG:4326) of **cell centres** carrying ``count``
    (features in the cell), ``per_km2`` (``count`` / cell area) and ``cell_km``,
    ready to drive a graduated symbol.

    Binning happens in ``crs`` -- EPSG:5070 CONUS Albers by default, which is
    equal-area, so "per square kilometre" means the same thing at the top and
    bottom of the AOI. (Doing this in 4326 or Web Mercator would inflate cell
    area with latitude, exactly the error that makes BRISK's "60 m" pixels 46.5 m
    of ground at 39 degrees N.) The grid is anchored on the projection origin
    rather than on the data, so the same cell boundaries fall in the same place
    for every AOI and the counts are reproducible.

    ``groups`` optionally restricts the count to some of :data:`GROUPS` -- a
    density map is only meaningful for one class at a time, since a dump and a
    prospect pit are not interchangeable units.
    """
    import numpy as np

    if groups is not None:
        if isinstance(groups, str):
            groups = [groups]
        groups = {str(g).lower() for g in groups}
        mines = mines[mines.group.isin(groups)] if "group" in mines else mines
    if mines is None or not len(mines):
        return gpd.GeoDataFrame({"count": pd.Series(dtype="int64"),
                                 "per_km2": pd.Series(dtype="float64"),
                                 "cell_km": pd.Series(dtype="float64")},
                                geometry=gpd.GeoSeries([], crs=4326), crs=4326)
    if cell_km <= 0:
        raise ValueError("cell_km must be > 0")

    pts = gpd.GeoDataFrame(geometry=_as_points(mines), crs=mines.crs)
    pts = pts.to_crs(crs)
    cell = float(cell_km) * 1000.0
    ix = np.floor(pts.geometry.x.to_numpy() / cell).astype("int64")
    iy = np.floor(pts.geometry.y.to_numpy() / cell).astype("int64")

    counts = (pd.DataFrame({"ix": ix, "iy": iy})
              .value_counts(["ix", "iy"]).rename("count").reset_index())
    cx = (counts.ix.to_numpy() + 0.5) * cell
    cy = (counts.iy.to_numpy() + 0.5) * cell
    area = float(cell_km) ** 2
    out = gpd.GeoDataFrame(
        dict(count=counts["count"].astype("int64"),
             per_km2=counts["count"].to_numpy() / area,
             cell_km=float(cell_km)),
        geometry=gpd.points_from_xy(cx, cy), crs=crs)
    return out.to_crs(4326).sort_values("count", ascending=False).reset_index(
        drop=True)


def density_raster(mines, bounds, cell_km=1.0, groups=None, crs="EPSG:5070"):
    """Mine-feature density as a raster (features per km^2).

    Returns an ``mrms``-style result dict -- ``fields``, ``profile``, ``meta``
    -- so :func:`stormscape.mrms.save_fields` writes it and
    :func:`stormscape.plot.drape_i15` drapes it like any other field.

    Built in an equal-area projection (EPSG:5070) on a grid anchored to the
    projection origin, so a cell is the same number of square kilometres
    everywhere and the same ground is binned identically from run to run.
    ``bounds`` is the ``(W,S,E,N)`` extent in degrees to cover.
    """
    import numpy as np
    import rasterio
    from rasterio.transform import from_origin

    if groups is not None:
        if isinstance(groups, str):
            groups = [groups]
        groups = {str(g).lower() for g in groups}
        if mines is not None and len(mines) and "group" in mines:
            mines = mines[mines.group.isin(groups)]
    if cell_km <= 0:
        raise ValueError("cell_km must be > 0")

    cell = float(cell_km) * 1000.0
    box = gpd.GeoDataFrame(
        geometry=gpd.GeoSeries.from_wkt(
            [f"POLYGON(({bounds[0]} {bounds[1]}, {bounds[2]} {bounds[1]}, "
             f"{bounds[2]} {bounds[3]}, {bounds[0]} {bounds[3]}, "
             f"{bounds[0]} {bounds[1]}))"]), crs=4326).to_crs(crs)
    x0, y0, x1, y1 = box.total_bounds
    # snap outward to whole cells so the grid lines up with density_grid()
    x0 = np.floor(x0 / cell) * cell
    y0 = np.floor(y0 / cell) * cell
    x1 = np.ceil(x1 / cell) * cell
    y1 = np.ceil(y1 / cell) * cell
    nx = max(1, int(round((x1 - x0) / cell)))
    ny = max(1, int(round((y1 - y0) / cell)))

    counts = np.zeros((ny, nx), dtype="float32")
    if mines is not None and len(mines):
        pts = gpd.GeoDataFrame(geometry=_as_points(mines),
                               crs=mines.crs).to_crs(crs)
        ix = np.floor((pts.geometry.x.to_numpy() - x0) / cell).astype("int64")
        iy = np.floor((pts.geometry.y.to_numpy() - y0) / cell).astype("int64")
        ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
        # rasters count rows downward from the top, the grid indexes upward
        np.add.at(counts, (ny - 1 - iy[ok], ix[ok]), 1.0)
    counts /= float(cell_km) ** 2                      # -> features per km^2

    profile = dict(driver="GTiff", height=ny, width=nx, count=1,
                   dtype="float32", crs=rasterio.crs.CRS.from_string(crs),
                   transform=from_origin(x0, y1, cell, cell), nodata=None,
                   compress="deflate")
    return dict(fields={"mine_density": counts}, profile=profile,
                meta=dict(cell_km=float(cell_km), crs=crs,
                          groups=sorted(groups) if groups else "all",
                          n_features=int(len(mines)) if mines is not None else 0,
                          peak_per_km2=float(counts.max())))


def group_counts(mines):
    """Tidy per-group tally of a mine frame -> DataFrame.

    Columns ``group``, ``description``, ``n``, ``types`` (the distinct feature
    types present, so the reader can see what a group actually contained here
    rather than trusting the grouping blind).
    """
    if mines is None or not len(mines):
        return pd.DataFrame(columns=["group", "description", "n", "types"])
    rows = []
    for grp, g in mines.groupby("group"):
        rows.append(dict(group=grp, description=GROUPS.get(grp, ""), n=len(g),
                         types="; ".join(sorted(
                             {str(t) for t in g.ftr_type.dropna()}))))
    out = pd.DataFrame(rows)
    order = {g: i for i, g in enumerate(GROUPS)}
    return (out.sort_values("group", key=lambda s: s.map(order))
            .reset_index(drop=True))


def describe_sources():
    """Print the source registry -- what is public, what needs a token."""
    for name, spec in SOURCES.items():
        gate = "public" if spec.get("public", True) else "TOKEN REQUIRED"
        print(f"{name:<12} {gate:<14} {spec['label']}")
        if spec.get("note"):
            print(f"{'':<27}{spec['note']}")
        if not spec.get("public", True) and spec.get("token_env"):
            have = "set" if os.environ.get(spec["token_env"]) else "not set"
            print(f"{'':<27}${spec['token_env']} is {have}")
    print(f"\ngroups: {', '.join(GROUPS)}", file=sys.stdout)
