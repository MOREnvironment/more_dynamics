"""Declaration helpers shared by the parts of the propulsor: the ranges that recur
(a strictly positive number, a fraction of one, a free scalar, a body vector) and
the one way a part becomes a named CasADi function (its couplings first, then
its own declared numbers, each input and output named).

A convention of this library, no source equation.

Author:    Enio Krizman
Date:      2026-10-23
"""
import casadi as ca
from more_transformations.more_casadi_transformations import Parameter


def positive(name, unit, meaning):
    return Parameter(name, (1, 1), unit, meaning, 0.0, minimum_exclusive=True)


def scalar(name, unit, meaning):
    return Parameter(name, (1, 1), unit, meaning)


def fraction(name, meaning):
    return Parameter(name, (1, 1), "1", meaning, 0.0, 1.0)


def vector(name, unit, meaning, rows=3):
    return Parameter(name, (rows, 1), unit, meaning)


def part_function(name, inputs, outputs, symbols, declared):
    """``inputs`` / ``outputs``: ordered ``{coupling name: SX}``; ``symbols``: the part's
    own symbols by declared name. Inputs are the couplings, then the declared numbers."""
    names = [d.name for d in declared]
    return ca.Function(name, [*inputs.values(), *[symbols[n] for n in names]], list(outputs.values()),
                       [*inputs, *names], list(outputs))
