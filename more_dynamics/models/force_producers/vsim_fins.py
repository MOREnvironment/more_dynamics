"""DUNE VSIM-style fin assembly: force ``K (delta / delta_max) U^2`` per fin.

Ported from ``more_generic_models``
``dynamics/propulsion/fins/fins_auv_actuation_model/fins_auv_actuation_vsim.py``
(``FinsActuationVSIM``, lines 21-93), including its zero-speed branch
(``|nu_r[0:3]| <= 1e-6`` gives zero) and its moment mapping
(``tau[4] = T_z``, ``tau[5] = T_y``).
"""

from dataclasses import dataclass
from typing import Sequence

import casadi as ca
import numpy as np

from ._common import positive, saturate

ZERO_SPEED = 1e-6  # fins_auv_actuation_vsim.py line 58


@dataclass(frozen=True)
class VsimFinsConstants:
    """``force_gain`` (3,N) and ``moment_gain`` (3,N) per unit normalised deflection."""

    max_deflection: float
    force_gain: np.ndarray
    moment_gain: np.ndarray


def preprocess_vsim_fins(
    max_forces: Sequence[Sequence[float]],
    positions: Sequence[Sequence[float]],
    max_deflection: float,
) -> VsimFinsConstants:
    """Source names ``max_force`` (N,3), ``position`` (N,3), ``max_act`` (rad)."""
    forces = np.asarray(max_forces, dtype=float)
    r = np.asarray(positions, dtype=float)
    if forces.ndim != 2 or forces.shape[1] != 3 or forces.shape != r.shape:
        raise ValueError("max_forces and positions must have the same shape (N, 3)")
    if not (np.all(np.isfinite(forces)) and np.all(np.isfinite(r))):
        raise ValueError("max_forces and positions must contain only finite values")
    moment = np.vstack(
        [
            np.sqrt(r[:, 1] ** 2 + r[:, 2] ** 2) * np.abs(forces[:, 0]) / 2,
            r[:, 0] * np.abs(forces[:, 1]) / 2,
            r[:, 0] * np.abs(forces[:, 2]) / 2,
        ]
    )
    return VsimFinsConstants(
        max_deflection=positive("max_deflection", max_deflection),
        force_gain=(forces / 2.0).T,
        moment_gain=moment,
    )


def vsim_fins_casadi(constants: VsimFinsConstants) -> ca.Function:
    """``(delta, nu_r) -> tau``; ``delta`` (Nx1) rad, saturated."""
    n_fins = constants.force_gain.shape[1]
    delta = ca.SX.sym("delta", n_fins)
    nu_r = ca.SX.sym("nu_r", 6)
    limit = constants.max_deflection
    normalized = saturate(delta, -limit, limit) / limit
    speed_squared = ca.sumsqr(nu_r[0:3])
    force = speed_squared * ca.mtimes(ca.DM(constants.force_gain), normalized)
    torque = speed_squared * ca.mtimes(ca.DM(constants.moment_gain), normalized)
    tau = ca.vertcat(force, torque[0], torque[2], torque[1])
    tau = ca.if_else(ca.sqrt(speed_squared) <= ZERO_SPEED, ca.DM.zeros(6), tau)
    return ca.Function("vsim_fins", [delta, nu_r], [tau], ["delta", "nu_r"], ["tau"])
