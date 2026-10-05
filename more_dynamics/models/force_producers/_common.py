"""Shared numpy checks and CasADi helpers of the force-producer blocks."""

from typing import Sequence

import casadi as ca
import numpy as np


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


def rotation_zyx(angles: np.ndarray) -> np.ndarray:
    """``R_zyx(phi, theta, psi)`` (Fossen 2021, eq. 2.52)."""
    phi, theta, psi = angles
    cphi, sphi = np.cos(phi), np.sin(phi)
    cth, sth = np.cos(theta), np.sin(theta)
    cpsi, spsi = np.cos(psi), np.sin(psi)
    return np.array(
        [
            [cpsi * cth, -spsi * cphi + cpsi * sth * sphi, spsi * sphi + cpsi * cphi * sth],
            [spsi * cth, cpsi * cphi + sphi * sth * spsi, -cpsi * sphi + sth * spsi * cphi],
            [-sth, cth * sphi, cth * cphi],
        ]
    )


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
    """``[f; r x f]`` for a point force at ``position`` (BODY, about the CO)."""
    r = ca.DM(position)
    return ca.vertcat(force, ca.cross(r, force))
