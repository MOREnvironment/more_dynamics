"""Two fixed propellers with a quadratic thrust law (differential thrust).

Ported from ``more_generic_models``
``dynamics/propulsion/thruster/thruster_differential/thruster_differential.py``
(``DifferentialThruster``): limits (lines 84-91, the ``n_max_override`` /
``n_min_override`` branch), ``saturation`` (152-154), ``thrust_from_speed``
(156-171), ``tau_from_thrust`` (194-196) and ``_compute_B`` (272-279). The
shaft-speed lag ``update_state`` (136-143) is a state, not part of this map.
Same law as MSS ``otter.m`` 132-137, 220-232.
"""

from dataclasses import dataclass
from typing import Sequence

import casadi as ca
import numpy as np

from ._common import array, saturate


@dataclass(frozen=True)
class DifferentialThrusterConstants:
    """Per thruster ``[left, right]``; ``allocation_matrix`` maps thrust to tau."""

    positive_thrust_coefficients: np.ndarray
    negative_thrust_coefficients: np.ndarray
    max_speed: np.ndarray
    min_speed: np.ndarray
    allocation_matrix: np.ndarray


def preprocess_differential_thruster(
    positive_thrust_coefficients: Sequence[float],
    negative_thrust_coefficients: Sequence[float],
    thruster_positions: Sequence[Sequence[float]],
    thruster_directions: Sequence[Sequence[float]],
    max_forward_thrust: Sequence[float],
    max_reverse_thrust: Sequence[float],
) -> DifferentialThrusterConstants:
    """Source names: ``k_pos``, ``k_neg`` (N s^2/rad^2), ``r_thruster_*``,
    ``d_thruster_*`` (normalised here), ``F_fwd_ref`` / ``F_rev_ref`` as forces
    in N. ``n_max = sqrt(F_fwd / k_pos)``, ``n_min = -sqrt(F_rev / k_neg)``."""
    k_pos = array(positive_thrust_coefficients, (2,), "positive_thrust_coefficients")
    k_neg = array(negative_thrust_coefficients, (2,), "negative_thrust_coefficients")
    positions = array(thruster_positions, (2, 3), "thruster_positions")
    directions = array(thruster_directions, (2, 3), "thruster_directions")
    forward = array(max_forward_thrust, (2,), "max_forward_thrust")
    reverse = array(max_reverse_thrust, (2,), "max_reverse_thrust")
    for name, values in (("positive_thrust_coefficients", k_pos),
                         ("negative_thrust_coefficients", k_neg),
                         ("max_forward_thrust", forward), ("max_reverse_thrust", reverse)):
        if np.any(values <= 0.0):
            raise ValueError(f"{name} must be positive")

    norms = np.linalg.norm(directions, axis=1)
    if np.any(norms <= 0.0):
        raise ValueError("thruster_directions must be non-zero")
    directions = directions / norms[:, None]

    allocation = np.zeros((6, 2))
    for i in range(2):
        allocation[0:3, i] = directions[i]
        allocation[3:6, i] = np.cross(positions[i], directions[i])

    return DifferentialThrusterConstants(
        positive_thrust_coefficients=k_pos,
        negative_thrust_coefficients=k_neg,
        max_speed=np.sqrt(forward / k_pos),
        min_speed=-np.sqrt(reverse / k_neg),
        allocation_matrix=allocation,
    )


def differential_thruster_casadi(constants: DifferentialThrusterConstants) -> ca.Function:
    """``(n, nu_r) -> tau``; ``n`` (2x1) shaft speed in rad/s, saturated;
    ``T = k_pos n|n|`` for ``n > 0``, else ``k_neg n|n|``; ``tau = B T``."""
    n = ca.SX.sym("n", 2)
    nu_r = ca.SX.sym("nu_r", 6)
    n_sat = saturate(n, ca.DM(constants.min_speed), ca.DM(constants.max_speed))
    quadratic = n_sat * ca.fabs(n_sat)
    thrust = ca.if_else(
        n_sat > 0.0,
        ca.DM(constants.positive_thrust_coefficients) * quadratic,
        ca.DM(constants.negative_thrust_coefficients) * quadratic,
    )
    tau = ca.mtimes(ca.DM(constants.allocation_matrix), thrust)
    return ca.Function("differential_thruster", [n, nu_r], [tau], ["n", "nu_r"], ["tau"])
