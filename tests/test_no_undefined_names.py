"""Every global a module reads must be bound somewhere in that module.

This exists because of a real escape. Routing the readers through
``stormscape.layout`` replaced ``os.path.join(...)`` with ``find(...)`` in
`plot.py`, but the accompanying ``from .layout import find`` never landed --
and the full suite still passed, because the three call sites are inside figure
functions (`diagnostic_panels`, `climatology_comparison`, the field-map helper)
that the offline tests do not execute. It surfaced only when a live `climate`
run raised ``NameError: name 'find' is not defined``.

Import-time checks cannot catch that: an undefined *global* is only resolved
when the line runs. So walk the symbol table instead and assert that every
global a function reads is either bound at module level or a builtin. That is
the pyflakes F821 check, done with the stdlib so it needs no new dependency.

Lazy imports inside a function are bound in that function's own scope, so they
are correctly ignored.
"""
from __future__ import annotations

import builtins
import pathlib
import symtable

import pytest

PKG = pathlib.Path(__file__).resolve().parent.parent / "stormscape"
MODULES = sorted(PKG.glob("*.py"))
BUILTINS = set(dir(builtins))


def _module_bindings(table: symtable.SymbolTable) -> set[str]:
    """Names bound at module scope: imports, defs, classes, assignments."""
    return {s.get_name() for s in table.get_symbols()
            if s.is_assigned() or s.is_imported() or s.is_namespace()}


def _free_globals(table: symtable.SymbolTable):
    """(name, scope) for every global a nested scope *reads* without binding."""
    for child in table.get_children():
        for sym in child.get_symbols():
            if sym.is_global() and not sym.is_assigned():
                yield sym.get_name(), child.get_name()
        yield from _free_globals(child)


@pytest.mark.parametrize("path", MODULES, ids=lambda p: p.name)
def test_module_has_no_undefined_globals(path):
    src = path.read_text()
    top = symtable.symtable(src, str(path), "exec")
    bound = _module_bindings(top) | BUILTINS

    missing = sorted({f"{name}  (used in {scope})"
                      for name, scope in _free_globals(top)
                      if name not in bound})
    assert not missing, (
        f"{path.name} reads globals it never binds -- a NameError waiting for "
        f"the right code path:\n  " + "\n  ".join(missing))
