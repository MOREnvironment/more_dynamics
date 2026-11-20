"""Mechanics shared by every vehicle part function: cutting a block to the
outputs a slot needs, and building a named CasADi function from symbols.

Author:    Enio Krizman
Date:      2026-10-08
"""

import casadi as ca


def function_from(name, inputs, outputs):
    """``ca.Function`` from ``{name: SX}`` inputs and outputs.
    ``allow_duplicate_io_names``: a primitive handed on unchanged as a
    coupling (e.g. ``hull_mass``) keeps its name on both sides."""
    return ca.Function(name, list(inputs.values()), list(outputs.values()), list(inputs), list(outputs),
                       {"allow_duplicate_io_names": True})


def restrict(function, outputs, name, rename=None):
    """``function`` cut to ``outputs``: the inputs those outputs do not
    depend on are dropped (checked with ``ca.depends_on``, so a dropped
    input cannot change a kept output). ``rename`` maps an output name to
    the coupling name a slot expects."""
    rename = rename or {}
    ins = {n: function.sx_in(i) for i, n in enumerate(function.name_in())}
    out = function.call(ins)
    kept = {rename.get(o, o): out[o] for o in outputs}
    stacked = ca.vertcat(*[ca.vec(e) for e in kept.values()])
    needed = [n for n in function.name_in() if ca.depends_on(stacked, ins[n])]
    return ca.Function(name, [ins[n] for n in needed], list(kept.values()), needed, list(kept.keys()))


def with_passthrough(function, name, extra):
    """``function`` plus some of its own inputs handed on unchanged as
    couplings (e.g. a rigid-body part passing its ``hull_mass`` on to the
    added-mass slot)."""
    ins = {n: function.sx_in(i) for i, n in enumerate(function.name_in())}
    out = function.call(ins)
    return function_from(name, ins, {**out, **{e: ins[e] for e in extra}})
