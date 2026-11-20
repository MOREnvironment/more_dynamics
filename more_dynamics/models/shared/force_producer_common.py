"""Shared CasADi helpers of the force-producer blocks.

Rotations and skew matrices come from
``more_transformations.more_casadi_transformations``; none is defined here.

Equations (keys in References):

* Point force ``f`` at ``r`` (BODY, from the CO): ``tau = [f; r x f] =
  [f; S(r) f]`` (Fossen 2011, eq. 12.226, p. 400; ``S(r) a = r x a``,
  Fossen 2011, eqs. 2.9-2.10, p. 20).
* Saturation ``min(max(x, lower), upper)``, the order of ``numpy.clip``.
* Safe norm ``|v| = sqrt(v . v)`` with a zero derivative at ``v = 0``.

References
----------
[Fossen 2011] Fossen, T. I. (2011). *Handbook of Marine Craft Hydrodynamics
    and Motion Control*, 1st ed. John Wiley & Sons, Chichester. Ch. 2,
    eqs. 2.9-2.10, p. 20; Ch. 12, Sec. 12.3.1, eq. 12.226, p. 400.

Author:    Enio Krizman
Date:      2026-10-05
"""

import casadi as ca
from more_transformations.more_casadi_transformations import MatrixTransforms


NAN_COMMAND_REFUSED = "saturate refuses a NaN command: the input is not a number (a NaN would be clipped to a limit)"


def saturate(value, lower, upper):
    """Elementwise ``min(max(x, lower), upper)`` (the order of
    ``numpy.clip``). Unit: that of ``value``.

    A NaN command is refused, never clipped (``fmin``/``fmax`` would turn it
    into a limit). By type of ``value``:

    * ``MX``: an assertion in the graph; evaluating the Function with a NaN
      input raises ``RuntimeError`` carrying ``NAN_COMMAND_REFUSED``.
    * ``SX``: an SX graph cannot assert, so a NaN input gives a NaN output
      that the caller sees; finite values are unchanged.
    * numbers (``float``, ``list``, ``DM``): raises ``ValueError``.
    """
    if isinstance(value, ca.MX):
        clipped = ca.fmin(ca.fmax(value, lower), upper)
        return ca.MX.attachAssert(clipped, ca.mmin(value >= -ca.inf), NAN_COMMAND_REFUSED)  # x >= -inf is false only for NaN
    if isinstance(value, ca.SX):
        clipped = ca.fmin(ca.fmax(value, lower), upper)
        return ca.if_else(value >= -ca.inf, clipped, float("nan"))
    if any(v != v for v in ca.DM(value).nonzeros()):  # v != v holds only for NaN
        raise ValueError(NAN_COMMAND_REFUSED)
    return ca.fmin(ca.fmax(value, lower), upper)


def safe_norm(vector):
    """``|v|`` with the same value as ``sqrt(v . v)`` and a zero (not NaN)
    derivative at ``v = 0`` (a guard of this module: the branch not taken
    gets a harmless argument)."""
    squared = ca.sumsqr(vector)  # v . v (definition of the Euclidean norm)
    nonzero = squared > 0.0
    return ca.if_else(nonzero, ca.sqrt(ca.if_else(nonzero, squared, 1.0)), 0.0)  # guard against sqrt'(0)


def wrench(force, position):
    """``[f; S(r) f]`` for a point force at ``position`` (BODY, about the CO;
    Fossen 2011, eq. 12.226, p. 400). ``position`` may be numbers or a
    CasADi expression."""
    return ca.vertcat(force, MatrixTransforms.skew(position) @ force)  # (Fossen 2011, eq. 12.226, p. 400; eq. 2.9, p. 20)
