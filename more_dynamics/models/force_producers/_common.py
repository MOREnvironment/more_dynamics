"""Shared numpy checks and CasADi helpers of the force-producer blocks.

Rotations and skew matrices come from ``more_transformations`` (numpy, in
the pre-processing) and ``more_transformations.more_casadi_transformations``
(in the graphs); none is defined here.

Equations (keys in References):

* Point force ``f`` at ``r`` (BODY, from the CO): ``tau = [f; r x f] =
  [f; S(r) f]`` (Fossen 2011, eq. 12.226, p. 400; MSS ``Smtrx.m`` 11-13).
* Saturation ``min(max(x, lower), upper)`` (MSS ``sat.m`` / ``satlim.m``
  behaviour, e.g. ``remus100.m`` 113-115, ``otter.m`` 219).

References
----------
[Fossen 2011] Fossen, T. I. (2011). *Handbook of Marine Craft Hydrodynamics
    and Motion Control*, 1st ed. John Wiley & Sons, Chichester. Ch. 12, §12.3.1, eq. 12.226, p. 400.
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2.
    https://github.com/cybergalactic/MSS, MIT licence, revision ``72656d1``:
    ``LIBRARY/kinematics/Smtrx.m`` 11-13; ``CRAFT/AUV/models/remus100.m``
    113-115; ``CRAFT/USV/models/otter.m`` 219.
"""

from typing import Sequence

import casadi as ca
import numpy as np
from more_transformations.more_casadi_transformations.matrix_transforms import (
    MatrixTransforms as CasadiMatrixTransforms,
)


def positive(name: str, value: float) -> float:
    value = float(value)
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be a positive finite value")
    return value


def finite(name: str, value: float) -> float:
    value = float(value)
    if not np.isfinite(value):
        raise ValueError(f"{name} must be a finite value")
    return value


def array(values: Sequence, shape: tuple, name: str) -> np.ndarray:
    result = np.asarray(values, dtype=float)
    if result.shape != shape:
        raise ValueError(f"{name} must have shape {shape}, got {result.shape}")
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain only finite values")
    return result


def saturate(value, lower, upper):
    """Elementwise ``clip`` (``np.clip`` order: ``min(max(x, lower), upper)``)."""
    return ca.fmin(ca.fmax(value, lower), upper)


def safe_norm(vector: ca.SX) -> ca.SX:
    """``|v|`` with the same value as ``sqrt(v . v)`` and a zero (not NaN)
    derivative at ``v = 0``."""
    squared = ca.sumsqr(vector)
    nonzero = squared > 0.0
    return ca.if_else(nonzero, ca.sqrt(ca.if_else(nonzero, squared, 1.0)), 0.0)


def wrench(force: ca.SX, position: np.ndarray) -> ca.SX:
    """``[f; S(r) f]`` for a point force at ``position`` (BODY, about the CO;
    Fossen 2011, eq. 12.226, p. 400)."""
    return ca.vertcat(force, ca.mtimes(CasadiMatrixTransforms.skew(position), force))  # (Fossen 2011, eq. 12.226, p. 400)
