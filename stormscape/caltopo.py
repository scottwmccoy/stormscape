"""stormscape.caltopo -- CalTopo-native GeoJSON, grouped into toggleable folders.

CalTopo imports plain GeoJSON, but it imports it *flat*: every LineString in the
file becomes its own independent map object with its own checkbox. Hand it a
debris-flow network of 1,400 segments and the sidebar gets 1,400 rows, with no
way to switch the layer off as a unit -- the single most common complaint about
exporting a modelled network into the field.

The fix is CalTopo's own container object, which it round-trips through GeoJSON
and reads on import. A **Folder** is a geometry-less feature::

    {"type": "Feature", "id": "<uuid>", "geometry": null,
     "properties": {"class": "Folder", "title": "Predicted DF network",
                    "visible": true, "labelVisible": false}}

and every shape that names it in ``folderId`` lands inside it::

    {"type": "Feature", "id": "<uuid>",
     "geometry": {"type": "LineString", "coordinates": [...]},
     "properties": {"class": "Shape", "title": "seg 41", "folderId": "<uuid>",
                    "description": "...", "stroke": "#D7191C",
                    "stroke-width": 3, "stroke-opacity": 1, "pattern": "solid"}}

One checkbox per folder, and the objects stay individually clickable for their
attributes. This module writes that structure from GeoDataFrames.

What CalTopo reads
------------------
* ``class`` -- ``Folder``, ``Shape`` (lines/polygons) or ``Marker`` (points).
* ``title`` -- the object label. **Not** ``name``: a shapefile's ``name``
  column imports as an unlabelled object, which is why :class:`Layer` takes a
  ``label`` column and copies it into ``title``.
* ``description`` -- free text in the object's info panel. The only place a
  modelled attribute can travel into the field, so :class:`Layer` folds chosen
  columns into it.
* Styling follows the Mapbox **simplestyle** spec: ``stroke``,
  ``stroke-width``, ``stroke-opacity``, ``fill``, ``fill-opacity``,
  ``pattern`` for shapes; ``marker-color``, ``marker-symbol``, ``marker-size``
  for markers. Colours are ``#RRGGBB``.

Vectors go in **EPSG:4326** -- GeoJSON is WGS84 by specification, and CalTopo
assumes it. Rasters are the opposite: they must be **EPSG:3857** and no larger
than 40 MB, which is :func:`stormscape.export.export_geotiffs`' job.
"""
from __future__ import annotations

import json
import os
import uuid
from dataclasses import dataclass, field
from typing import Any, Optional, Sequence, Union

WGS84 = "EPSG:4326"

# CalTopo writes 4-element positions -- [lon, lat, elevation, timestamp] -- and
# fills the last two with 0 when it has neither. RFC 7946 says SHOULD NOT go
# past three, so matching CalTopo costs something on the way back: GDAL warns
# ("too many members in array ... At most 3 are handled"), keeps the third, and
# reopens the file as 3-D geometry with z=0. Harmless in CalTopo, which is the
# consumer this file is for; set ``coord_len=2`` for a QGIS-clean strict file.
COORD_LEN = 4

# Colorblind-safe five-class ramp (ColorBrewer RdYlBu, reversed): blue = low
# through red = high. Class breaks are the caller's; these are the colours.
CLASS_COLORS = ("#2C7BB6", "#ABD9E9", "#FFFFBF", "#FDAE61", "#D7191C")


def _hex(color: str) -> str:
    """Normalize ``RRGGBB`` / ``#rrggbb`` to CalTopo's ``#RRGGBB``."""
    s = str(color).strip().lstrip("#").upper()
    return f"#{s}"


def _new_id() -> str:
    return str(uuid.uuid4())


def _round_coords(coords, ndigits: int, coord_len: int):
    """Round a nested coordinate structure and pad positions to ``coord_len``.

    Recurses on the nesting rather than switching on the geometry type, so one
    routine covers Point through MultiPolygon.
    """
    if coords and isinstance(coords[0], (int, float)):
        pos = [round(float(c), ndigits) for c in coords[:2]]
        return pos + [0] * (coord_len - 2)
    return [_round_coords(c, ndigits, coord_len) for c in coords]


def folder(title: str, *, visible: bool = True, label_visible: bool = False,
           id: Optional[str] = None) -> dict:
    """A CalTopo Folder feature -- the container that gives a layer one checkbox.

    ``label_visible`` defaults **off**: a folder with labels on paints its
    title beside every object inside it, which is unreadable at network scale.
    """
    return {"type": "Feature", "id": id or _new_id(), "geometry": None,
            "properties": {"class": "Folder", "title": str(title),
                           "visible": bool(visible),
                           "labelVisible": bool(label_visible)}}


def classify(values, breaks: Sequence[float],
             colors: Sequence[str] = CLASS_COLORS) -> list:
    """Per-value colours from ``breaks``, low class first.

    ``len(colors)`` must be ``len(breaks) + 1``; a value lands in class *i* when
    ``breaks[i-1] <= v < breaks[i]``. Non-finite values return ``None`` so the
    caller can fall back to the layer colour rather than paint a NaN.
    """
    if len(colors) != len(breaks) + 1:
        raise ValueError(f"{len(breaks)} breaks need {len(breaks) + 1} colours, "
                         f"got {len(colors)}")
    out = []
    for v in values:
        try:
            x = float(v)
        except (TypeError, ValueError):
            out.append(None)
            continue
        if x != x:                                   # NaN
            out.append(None)
            continue
        i = 0
        while i < len(breaks) and x >= breaks[i]:
            i += 1
        out.append(_hex(colors[i]))
    return out


def class_labels(breaks: Sequence[float], *, fmt: str = "{:g}",
                 unit: str = "") -> list:
    """Human-readable names for the classes :func:`classify` assigns.

    Used for folder titles when a layer is split by class, so the folder list
    doubles as the legend CalTopo does not draw.
    """
    u = f" {unit}" if unit else ""
    names = [f"< {fmt.format(breaks[0])}{u}"]
    names += [f"{fmt.format(a)}-{fmt.format(b)}{u}"
              for a, b in zip(breaks[:-1], breaks[1:])]
    names.append(f"≥ {fmt.format(breaks[-1])}{u}")
    return names


@dataclass
class Layer:
    """One folder's worth of features.

    ``data`` is a GeoDataFrame; it is reprojected to WGS84 if it is not already.
    ``color`` is a single ``#RRGGBB`` or one colour per row (from
    :func:`classify`); a ``None`` in a per-row sequence falls back to
    ``fallback_color``. ``label`` names a column to copy into CalTopo's
    ``title``, and ``fields`` names columns to fold into the ``description`` --
    the only route a modelled attribute takes into the field.
    """
    title: str
    data: Any                                        # GeoDataFrame
    color: Union[str, Sequence[Optional[str]]] = "#D7191C"
    fallback_color: str = "#808080"
    width: float = 3.0
    opacity: float = 1.0
    fill: Optional[str] = None                       # None -> no fill written
    fill_opacity: float = 0.15
    pattern: str = "solid"
    label: Optional[str] = None
    fields: Sequence[str] = ()
    description: Optional[str] = None                # constant, prepended
    marker_symbol: str = "point"
    marker_size: str = "1"
    visible: bool = True
    label_visible: bool = False
    simplify_m: float = 0.0                          # 0 = keep every vertex


def _describe(row, fields: Sequence[str], const: Optional[str]) -> str:
    """``key: value`` lines for the object's info panel, blanks dropped."""
    parts = [const] if const else []
    for f in fields:
        if f not in row:
            continue
        v = row[f]
        if v is None or v != v or (isinstance(v, str) and not v.strip()):
            continue
        parts.append(f"{f}: {round(v, 4) if isinstance(v, float) else v}")
    return "\n".join(str(p) for p in parts)


def _layer_features(layer: Layer, folder_id: str, *, ndigits: int,
                    coord_len: int) -> list:
    """Every row of ``layer`` as a CalTopo feature inside ``folder_id``."""
    import geopandas as gpd  # noqa: F401  (documents the expected input type)

    g = layer.data
    if g is None or not len(g):
        return []
    if g.crs is not None and str(g.crs).upper() not in (WGS84, "EPSG:4326"):
        g = g.to_crs(WGS84)
    if layer.simplify_m:
        # Simplify in a metric CRS -- a tolerance in degrees is a different
        # distance at every latitude, and the whole point is a metre budget.
        metric = g.to_crs("EPSG:5070")
        g = g.set_geometry(metric.geometry.simplify(layer.simplify_m)
                           .to_crs(WGS84))

    n = len(g)
    colors = ([layer.color] * n if isinstance(layer.color, str)
              else list(layer.color))
    if len(colors) != n:
        raise ValueError(f"{layer.title}: {len(colors)} colours for {n} rows")

    out = []
    for (_, row), color in zip(g.iterrows(), colors):
        geom = row.geometry
        if geom is None or geom.is_empty:
            continue
        gj = geom.__geo_interface__
        gtype = gj["type"]
        col = _hex(color) if color else _hex(layer.fallback_color)
        props = {
            "class": "Marker" if gtype in ("Point", "MultiPoint") else "Shape",
            "title": (str(row[layer.label]) if layer.label and layer.label in row
                      and row[layer.label] == row[layer.label] else ""),
            "folderId": folder_id,
            "description": _describe(row, layer.fields, layer.description),
            "visible": True,
        }
        if props["class"] == "Marker":
            props.update({"marker-color": col,
                          "marker-symbol": layer.marker_symbol,
                          "marker-size": str(layer.marker_size)})
        else:
            props.update({"stroke": col, "stroke-width": layer.width,
                          "stroke-opacity": layer.opacity,
                          "pattern": layer.pattern})
            if layer.fill is not None:
                props.update({"fill": _hex(layer.fill),
                              "fill-opacity": layer.fill_opacity})
        out.append({
            "type": "Feature", "id": _new_id(),
            "geometry": {"type": gtype,
                         "coordinates": _round_coords(gj["coordinates"],
                                                      ndigits, coord_len)},
            "properties": props})
    return out


def build(layers: Sequence[Layer], *, ndigits: int = 6,
          coord_len: int = COORD_LEN) -> dict:
    """The CalTopo FeatureCollection for ``layers`` -- one folder each.

    Folders are emitted before their contents so a reader that streams the
    collection has the container before anything references it. ``ndigits``
    rounds coordinates (6 = ~0.1 m, and a large fraction of the file size).
    """
    feats, folders = [], []
    for layer in layers:
        fid = _new_id()
        shapes = _layer_features(layer, fid, ndigits=ndigits,
                                 coord_len=coord_len)
        if not shapes:
            continue
        folders.append(folder(layer.title, visible=layer.visible,
                              label_visible=layer.label_visible, id=fid))
        feats.extend(shapes)
    return {"type": "FeatureCollection", "features": folders + feats}


def write(layers: Sequence[Layer], out_path: str, *, ndigits: int = 6,
          coord_len: int = COORD_LEN) -> str:
    """Write :func:`build`'s collection to ``out_path``. Returns the path."""
    fc = build(layers, ndigits=ndigits, coord_len=coord_len)
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(fc, fh, separators=(",", ":"))
    return out_path


def summary(path: str) -> dict:
    """``{folder title: feature count}`` for a written file, plus its size.

    The check worth running before a field trip: it reads back what CalTopo
    will see, so an empty folder shows up here rather than on a phone with no
    signal.
    """
    with open(path, encoding="utf-8") as fh:
        fc = json.load(fh)
    names, counts = {}, {}
    for f in fc["features"]:
        p = f.get("properties", {})
        if p.get("class") == "Folder":
            names[f["id"]] = p.get("title", "")
        else:
            counts[p.get("folderId")] = counts.get(p.get("folderId"), 0) + 1
    return {"path": path, "bytes": os.path.getsize(path),
            "folders": {names[k]: counts.get(k, 0) for k in names},
            "ungrouped": counts.get(None, 0)}
