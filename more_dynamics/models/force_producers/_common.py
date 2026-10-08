"""Shared CasADi helpers of the force-producer blocks.

Rotations and skew matrices come from
``more_transformations.more_casadi_transformations``; none is defined here.

Equations (keys in References):

* Point force ``f`` at ``r`` (BODY, from the CO): ``tau = [f; r x f] =
  [f; S(r) f]`` (Fossen 2011, eq. 12.226, p. 400; ``S(r) a = r x a``,
  Fossen 2011, eqs. 2.9-2.10, p. 20; MSS ``Smtrx.m`` 11-13).
* Saturation ``min(max(x, lower), upper)`` (MSS ``sat.m`` / ``satlim.m``
  behaviour, e.g. ``remus100.m`` 113-115, ``otter.m`` 220).

References
----------
[Fossen 2011] Fossen, T. I. (2011). *Handbook of Marine Craft Hydrodynamics
    and Motion Control*, 1st ed. John Wiley & Sons, Chichester. Ch. 2,
    eqs. 2.9-2.10, p. 20; Ch. 12, §12.3.1, eq. 12.226, p. 400.
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2
    with the fixes of 2026-10-07. https://github.com/cybergalactic/MSS, MIT
    licence, revision ``cc07579``: ``LIBRARY/kinematics/Smtrx.m`` 11-13;
    ``CRAFT/AUV/models/remus100.m`` 113-115; ``CRAFT/USV/models/otter.m``
    220.

Author:    Enio Krizman
Date:      2026-10-05
"""

import casadi as ca
from more_transformations.more_casadi_transformations import MatrixTransforms


def saturate(value, lower, upper):
    """Elementwise ``min(max(x, lower), upper)`` (the order of ``np.clip``;
    MSS ``satlim.m`` behaviour, ``otter.m`` 220)."""
    return ca.fmin(ca.fmax(value, lower), upper)  # (otter.m 220; remus100.m 113-115)


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
