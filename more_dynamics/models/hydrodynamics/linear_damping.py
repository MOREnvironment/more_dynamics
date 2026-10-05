"""Linear damping: ``D`` and ``tau = -D nu_r`` (numpy constants, CasADi algebra).

Ported from ``more_generic_models`` ``dynamics/plant/matrices/linear_damping.py``:
``DampingMatrix.get_damping_coeff`` (lines 16-84), ``D_linear`` (88-112),
``D_linear_catamaran`` (165-197) and ``get_tau_damping`` (200-226); translated
against MSS (MIT, T. I. Fossen) ``LIBRARY/modeling/Dmtrx.m`` lines 30, 44-49,
``remus100.m`` line 216 and ``otter.m`` lines 197-207, 235-240.

Two sign conventions meet here. The submerged form follows ``Dmtrx.m``: a
positive ``D`` that is subtracted. The surface form follows ``otter.m``:
negative derivatives ``[Xu Yv Zw Kp Mq Nr]`` whose product with ``nu_r`` is
added. Both functions return a positive ``D`` and the force ``tau`` on the
vehicle, so the caller adds ``tau`` to the right-hand side in either case.
"""

from dataclasses import dataclass
from typing import Sequence

import casadi as ca
import numpy as np


@dataclass(frozen=True)
class SubmergedLinearDampingConstants:
    """``Dmtrx.m`` diagonal before any speed fade."""

    damping_coefficients: np.ndarray
    sway_damping_fade: bool


@dataclass(frozen=True)
class SurfaceLinearDampingConstants:
    """``[Xu Yv Zw Kp Mq Nr]`` (negative) and the quadratic yaw factor."""

    damping_derivatives: np.ndarray
    yaw_damping_nonlinearity: float


def _matrix(values: np.ndarray, name: str) -> np.ndarray:
    matrix = np.asarray(values, dtype=float)
    if matrix.shape != (6, 6):
        raise ValueError(f"{name} must be a 6x6 matrix")
    if not np.all(np.isfinite(matrix)):
        raise ValueError(f"{name} must contain only finite values")
    return matrix


def _vector(values: Sequence[float], size: int, name: str) -> np.ndarray:
    vector = np.asarray(values, dtype=float)
    if vector.shape != (size,):
        raise ValueError(f"{name} must contain exactly {size} values")
    if not np.all(np.isfinite(vector)):
        raise ValueError(f"{name} must contain only finite values")
    return vector


def _positive(values: np.ndarray, name: str) -> None:
    if np.any(np.asarray(values) <= 0.0):
        raise ValueError(f"{name} must be positive")


def _relative_speed(nu_r: ca.SX) -> ca.SX:
    """``U_r = sqrt(u_r^2 + v_r^2 + w_r^2)`` (``remus100.m`` line 128)."""
    return ca.sqrt(nu_r[0] ** 2 + nu_r[1] ** 2 + nu_r[2] ** 2)


def preprocess_submerged_linear_damping(
    rigid_body_mass_matrix: np.ndarray,
    added_mass_matrix: np.ndarray,
    weight: float,
    center_of_gravity: Sequence[float],
    center_of_buoyancy: Sequence[float],
    time_constants: Sequence[float],
    damping_ratios: Sequence[float],
    sway_damping_fade: bool = False,
) -> SubmergedLinearDampingConstants:
    """``Dmtrx.m`` diagonal (``get_damping_coeff``).

    Source names: ``MRB``, ``MA``, ``W``, ``r_bg``, ``r_bb``, ``T_126``,
    ``zeta_45``. ``T3 = T2``; ``w4``, ``w5`` from the metacentric height
    ``z_g - z_b``. ``sway_damping_fade=True`` is the template departure
    (``D_linear`` fades sway as well as surge); the default fades surge only.
    """
    mass = _matrix(rigid_body_mass_matrix, "rigid_body_mass_matrix") + _matrix(
        added_mass_matrix, "added_mass_matrix"
    )
    r_g = _vector(center_of_gravity, 3, "center_of_gravity")
    r_b = _vector(center_of_buoyancy, 3, "center_of_buoyancy")
    t1, t2, t6 = _vector(time_constants, 3, "time_constants")
    zeta4, zeta5 = _vector(damping_ratios, 2, "damping_ratios")
    _positive(np.array([t1, t2, t6]), "time_constants")
    _positive(np.diag(mass), "diagonal of the total mass matrix")
    if not np.isfinite(weight) or weight <= 0.0:
        raise ValueError("weight must be a positive finite value")
    metacentric_height = r_g[2] - r_b[2]
    if metacentric_height <= 0.0:
        raise ValueError("center_of_gravity must lie below center_of_buoyancy (z down)")

    t3 = t2
    w4 = np.sqrt(weight * metacentric_height / mass[3, 3])
    w5 = np.sqrt(weight * metacentric_height / mass[4, 4])
    coefficients = np.array(
        [
            mass[0, 0] / t1,
            mass[1, 1] / t2,
            mass[2, 2] / t3,
            mass[3, 3] * 2 * zeta4 * w4,
            mass[4, 4] * 2 * zeta5 * w5,
            mass[5, 5] / t6,
        ]
    )
    return SubmergedLinearDampingConstants(
        damping_coefficients=coefficients,
        sway_damping_fade=bool(sway_damping_fade),
    )


def submerged_linear_damping_casadi(
    constants: SubmergedLinearDampingConstants,
) -> ca.Function:
    """``nu_r -> (D, tau)``; ``D(1,1) *= exp(-3 U_r)`` (``remus100.m`` 216).

    ``D(2,2)`` fades the same way only with ``sway_damping_fade`` (``D_linear``
    lines 109-110). ``tau = -D nu_r``.
    """
    nu_r = ca.SX.sym("nu_r", 6)
    fade = ca.exp(-3 * _relative_speed(nu_r))
    diagonal = [ca.SX(c) for c in constants.damping_coefficients]
    diagonal[0] = diagonal[0] * fade
    if constants.sway_damping_fade:
        diagonal[1] = diagonal[1] * fade
    damping = ca.diag(ca.vertcat(*diagonal))
    return ca.Function(
        "submerged_linear_damping",
        [nu_r],
        [damping, -damping @ nu_r],
        ["nu_r"],
        ["D", "tau"],
    )


def preprocess_surface_linear_damping(
    mass_matrix: np.ndarray,
    restoring_matrix: np.ndarray,
    max_forward_thrust: float,
    max_speed: float,
    time_constants: Sequence[float],
    damping_ratios: Sequence[float],
    yaw_damping_nonlinearity: float,
) -> SurfaceLinearDampingConstants:
    """``otter.m`` 197-207 (``D_linear_catamaran``).

    ``mass_matrix`` = ``M_RB + M_A``; ``restoring_matrix`` about the centre of
    flotation (``G_CF``; only ``G33``, ``G44``, ``G55`` used);
    ``Xu = -max_forward_thrust / max_speed`` (the source's ``coeff[0] g / Umax``
    with the thrust in newtons); ``time_constants = (T_sway, T_yaw)``;
    ``damping_ratios = (zeta3, zeta4, zeta5)``.
    """
    mass = _matrix(mass_matrix, "mass_matrix")
    restoring = _matrix(restoring_matrix, "restoring_matrix")
    t_sway, t_yaw = _vector(time_constants, 2, "time_constants")
    zeta3, zeta4, zeta5 = _vector(damping_ratios, 3, "damping_ratios")
    _positive(np.array([t_sway, t_yaw]), "time_constants")
    _positive(np.diag(mass), "diagonal of the mass matrix")
    if not np.isfinite(max_speed) or max_speed <= 0.0:
        raise ValueError("max_speed must be a positive finite value")
    if not np.isfinite(max_forward_thrust) or max_forward_thrust < 0.0:
        raise ValueError("max_forward_thrust must be a non-negative finite value")
    if not np.isfinite(yaw_damping_nonlinearity) or yaw_damping_nonlinearity < 0.0:
        raise ValueError("yaw_damping_nonlinearity must be a non-negative finite value")
    stiffness = np.array([restoring[2, 2], restoring[3, 3], restoring[4, 4]])
    if np.any(stiffness < 0.0):
        raise ValueError("G33, G44 and G55 of restoring_matrix must be non-negative")

    w3, w4, w5 = np.sqrt(stiffness / np.diag(mass)[2:5])
    derivatives = np.array(
        [
            -max_forward_thrust / max_speed,
            -mass[1, 1] / t_sway,
            -2 * zeta3 * w3 * mass[2, 2],
            -2 * zeta4 * w4 * mass[3, 3],
            -2 * zeta5 * w5 * mass[4, 4],
            -mass[5, 5] / t_yaw,
        ]
    )
    return SurfaceLinearDampingConstants(
        damping_derivatives=derivatives,
        yaw_damping_nonlinearity=float(yaw_damping_nonlinearity),
    )


def surface_linear_damping_casadi(
    constants: SurfaceLinearDampingConstants,
) -> ca.Function:
    """``nu_r -> (D, tau)``; ``D = -diag([Xu Yv Zw Kp Mq Nr])`` (positive).

    ``tau`` = ``otter.m`` 235-240 (``get_tau_damping``): the derivatives times
    ``nu_r``, yaw times ``1 + k |r|`` with ``k = yaw_damping_nonlinearity``.
    """
    nu_r = ca.SX.sym("nu_r", 6)
    derivatives = constants.damping_derivatives
    forces = [derivatives[i] * nu_r[i] for i in range(5)]
    forces.append(
        derivatives[5]
        * (1 + constants.yaw_damping_nonlinearity * ca.fabs(nu_r[5]))
        * nu_r[5]
    )
    damping = ca.SX(ca.DM(-np.diag(derivatives)))
    return ca.Function(
        "surface_linear_damping",
        [nu_r],
        [damping, ca.vertcat(*forces)],
        ["nu_r"],
        ["D", "tau"],
    )
