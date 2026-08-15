"""Where a storm's products live inside its ``--out-dir``.

A finished event is a few dozen files of four different kinds -- GeoTIFF fields,
PNG figures, CSV tables, vector AOIs -- and dropping them all in one directory
makes the figures hard to find and the whole folder hard to hand to anyone else.
So by default stormscape sorts them::

    <out-dir>/
      figures/    *.png *.pdf *.svg      (+ VirtualGaugeFigures/ per-gauge detail)
      rasters/    *.tif *.tiff *.vrt
      tables/     *.csv *.md *.txt *.json
      vectors/    *.geojson *.gpkg *.shp
      RainGaugeData/  nexrad_cache/  atlas14_cache/     (unchanged)

Two rules make this safe to turn on by default:

**Readers auto-detect.** :func:`find` looks in the sorted subdirectory first and
then falls back to the flat path, so every ``--from-dir`` / ``--radar-dir``
chain keeps working against event folders written by older versions -- and
against hand-made folders. Nothing needs migrating.

**It is one flag to turn off.** ``--flat`` on any writing subcommand (or
``layout="flat"`` in the library, or ``STORMSCAPE_LAYOUT=flat``) restores the
old single-directory behaviour exactly.

Caches and ``RainGaugeData/`` are deliberately left at the top level: they are
inputs and intermediate stores rather than products, and ``nexrad_cache/`` in
particular is large and wants to stay obvious to delete.
"""
from __future__ import annotations

import os

#: extension -> subdirectory. Anything not listed stays at the top level.
SUBDIR = {
    ".tif": "rasters", ".tiff": "rasters", ".vrt": "rasters",
    ".png": "figures", ".pdf": "figures", ".svg": "figures",
    ".csv": "tables", ".md": "tables", ".txt": "tables", ".json": "tables",
    ".geojson": "vectors", ".gpkg": "vectors", ".shp": "vectors",
}

#: directories that keep their own name under the event root, unsorted.
RESERVED = ("RainGaugeData", "StreamGaugeData", "nexrad_cache", "atlas14_cache",
            "brisk_cache", "figures", "rasters", "tables", "vectors")

#: files that describe the whole event folder and belong at its root.
ROOT_FILES = {"README.md", "readme.md", "NOTES.md", "notes.md"}

_ENV = "STORMSCAPE_LAYOUT"


def sorted_layout(layout=None) -> bool:
    """Whether to sort products into subdirectories.

    ``layout`` is ``None`` (use ``$STORMSCAPE_LAYOUT``, default sorted),
    ``"sorted"``/``True``, or ``"flat"``/``False``.
    """
    if layout is None:
        layout = os.environ.get(_ENV, "sorted")
    if isinstance(layout, bool):
        return layout
    return str(layout).strip().lower() != "flat"


def subdir_for(filename: str) -> str | None:
    """Subdirectory a product with this name belongs in, or ``None`` for root."""
    base = os.path.basename(filename)
    if base in ROOT_FILES:
        return None
    return SUBDIR.get(os.path.splitext(base)[1].lower())


def out_path(out_dir: str, filename: str, layout=None, make=True) -> str:
    """Absolute path to write ``filename`` under ``out_dir``.

    Creates the subdirectory unless ``make=False``. With a flat layout, or an
    extension we do not sort, this is just ``out_dir/filename``.
    """
    sub = subdir_for(filename) if sorted_layout(layout) else None
    d = os.path.join(out_dir, sub) if sub else out_dir
    if make:
        os.makedirs(d, exist_ok=True)
    return os.path.join(d, os.path.basename(filename))


def _candidates(in_dir: str, base: str):
    """Where ``base`` might live, best first, given ``in_dir``.

    Also searches the *parent* when ``in_dir`` is itself one of the sorted
    subdirectories. Pointing ``--from-dir`` at ``<event>/rasters`` is a natural
    thing to do -- it is where the GeoTIFFs are -- but the event AOI lives in
    ``<event>/vectors`` and the gauge store beside it, so a lookup anchored at
    ``rasters/`` would miss them and silently degrade (the climate figures fall
    back to the i15 footprint instead of the event AOI). Treating a layout
    subdirectory as "inside the event" makes that call behave exactly like
    passing the event root.

    A directory that merely happens to be named ``rasters`` only ever gains a
    fallback: local hits still win.
    """
    sub = subdir_for(base)
    out = []
    if sub:
        out.append(os.path.join(in_dir, sub, base))
    out.append(os.path.join(in_dir, base))
    if os.path.basename(os.path.normpath(in_dir)) in RESERVED:
        parent = os.path.dirname(os.path.normpath(in_dir))
        if sub:
            out.append(os.path.join(parent, sub, base))
        out.append(os.path.join(parent, base))
    return out


def find(in_dir: str, filename: str) -> str:
    """Path to read ``filename`` from ``in_dir``, sorted layout or flat.

    Prefers the sorted subdirectory, falls back to the flat path, then to the
    same pair on the parent when ``in_dir`` is a layout subdirectory (see
    :func:`_candidates`). When nothing exists the *sorted* path under ``in_dir``
    is returned, so a "file not found" message points at where a fresh run would
    have put it.
    """
    cand = _candidates(in_dir, os.path.basename(filename))
    for p in cand:
        if os.path.exists(p):
            return p
    return cand[0]


def subdir(out_dir: str, name: str, layout=None, make=True) -> str:
    """Path to a named product *directory* (e.g. ``VirtualGaugeFigures``).

    Figure directories nest under ``figures/`` in a sorted layout;
    :data:`RESERVED` names stay at the event root either way.
    """
    if sorted_layout(layout) and name not in RESERVED:
        d = os.path.join(out_dir, "figures", name)
    else:
        d = os.path.join(out_dir, name)
    if make:
        os.makedirs(d, exist_ok=True)
    return d


def find_subdir(in_dir: str, name: str) -> str:
    """Locate a product directory written by either layout (see :func:`find`).

    Climbs to the parent on a layout subdirectory for the same reason
    :func:`find` does -- ``--from-dir <event>/rasters`` must still find
    ``RainGaugeData/`` at the event root.
    """
    cand = [os.path.join(in_dir, "figures", name), os.path.join(in_dir, name)]
    if os.path.basename(os.path.normpath(in_dir)) in RESERVED:
        parent = os.path.dirname(os.path.normpath(in_dir))
        cand += [os.path.join(parent, "figures", name), os.path.join(parent, name)]
    for p in cand:
        if os.path.isdir(p):
            return p
    return cand[1]
