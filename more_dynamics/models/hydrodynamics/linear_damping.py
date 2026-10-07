"""Linear damping: ``D`` and ``tau = -D nu_r`` (numpy constants, CasADi algebra).

Three forms (keys in References):

* submerged: MSS ``Dmtrx.m``'s submerged branch with the surge fade of
  ``remus100.m``;
* floating: ``Dmtrx.m``'s surface-craft branch, all six terms from time
  constants and damping ratios (the numpy source's ``coeff`` multipliers are
  not kept: ``coeff_i`` scales ``D_ii``, which ``T_i / coeff_i`` or
  ``coeff_i zeta_i`` gives);
* surface: the form of MSS ``otter.m``, surge from a thrust and a top speed,
  yaw with a quadratic term.

Equations:

* ``D = -diag(Xu, Yv, Zw, Kp, Mq, Nr)`` positive, ``tau = -D nu_r``
  (Fossen 2011, eq. 6.62, p. 123; dissipative, Property 6.3, eqs.
  6.58-6.60, p. 123).
* Surge, sway, yaw from time constants ``D_ii = M_ii / T_i``; heave, roll,
  pitch from damping ratios ``D_ii = 2 zeta_i w_i M_ii`` (Fossen 2011,
  eqs. 6.71, p. 124, and 6.76-6.81, p. 125; MSS ``Dmtrx.m`` 30, 48-49,
  62-63).
* Natural frequencies: submerged ``w4 = sqrt(W (z_g - z_b) / M44)``,
  ``w5 = sqrt(W (z_g - z_b) / M55)``, ``T3 = T2`` (MSS ``Dmtrx.m`` 44-46);
  floating ``w_i = sqrt(G_ii / M_ii)``, i = 3, 4, 5 (Fossen 2011,
  eqs. 4.51-4.53, p. 68; MSS ``Dmtrx.m`` 58-60, ``otter.m`` 196-198).
* Submerged surge fade ``D11 = D11 exp(-3 U_r)``, ``U_r = sqrt(u_r^2 + v_r^2
  + w_r^2)`` (MSS ``remus100.m`` 127, 218).
* Surface: ``Xu = -thrust_max / U_max``, ``Yv = -M22 / T_sway``,
  ``Zw = -2 zeta3 w3 M33``, ``Kp = -2 zeta4 w4 M44``, ``Mq = -2 zeta5 w5 M55``,
  ``Nr = -M66 / T_yaw`` (MSS ``otter.m`` 201-206);
  ``tau = [Xu u_r, Yv v_r, Zw w_r, Kp p_r, Mq q_r, Nr (1 + k |r_r|) r_r]``
  (MSS ``otter.m`` 234-239).

Two sign conventions meet here. ``Dmtrx.m`` gives a positive ``D`` that is
subtracted. ``otter.m`` gives negative derivatives ``[Xu Yv Zw Kp Mq Nr]``
whose product with ``nu_r`` is added. Every function returns a positive ``D``
and the force ``tau`` on the vehicle, so the caller adds ``tau`` to the
right-hand side in each case.

Deviations from MSS:

* Inputs outside their domain raise ``ValueError`` naming the input, where
  MSS returns a complex, infinite or NaN number (a centre of gravity not below
  the centre of buoyancy, a negative restoring stiffness, a non-positive mass
  diagonal, time constant, weight or top speed). Negative ``damping_ratios``
  raise too (zero is allowed: an undamped mode); MSS does not reject them,
  and a negative ratio there gives a negative damping diagonal, energy fed
  into the mode against a positive restoring stiffness (Fossen 2011,
  Property 6.3, p. 123).
* MSS fixes ``zeta3 = 0.2`` inside the surface-craft branch (``Dmtrx.m`` 57);
  here it is a parameter.
* ``smooth_speed_epsilon`` (submerged form, default 0 = MSS exactly). The
  speed ``U_r`` of ``remus100.m`` 127 has the derivative 0/0 at ``u_r = v_r =
  w_r = 0``, so the CasADi Jacobian of ``tau`` has non-finite entries there. A
  positive value ``eps`` (m/s) replaces it by ``sqrt(u_r^2 + v_r^2 + w_r^2 +
  eps^2)``, which is smooth everywhere and departs from MSS by at most
  ``eps`` in the speed.

Ported from the numpy source [MGM] ``dynamics/plant/matrices/linear_damping.py``:
``DampingMatrix.get_damping_coeff`` (16-84), ``D_linear`` (88-112),
``get_D_linear_surface_vessel`` (115-160), ``D_linear_catamaran`` (165-197)
and ``get_tau_damping`` (200-226).

References
----------
[Fossen 2011] Fossen, T. I. (2011). *Handbook of Marine Craft Hydrodynamics
    and Motion Control*, 1st ed. John Wiley & Sons, Chichester. Ch. 4,
    eqs. 4.51-4.53, p. 68; Ch. 6, §6.4.1, eqs. 6.58-6.81, pp. 123-125.
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2.
    https://github.com/cybergalactic/MSS, MIT licence, revision ``72656d1``:
    ``LIBRARY/modeling/Dmtrx.m`` 30-63; ``CRAFT/AUV/models/remus100.m`` 127,
    217-218; ``CRAFT/USV/models/otter.m`` 196-206, 234-239.
[MGM] Krizman, E. *more_generic_models*.
    https://github.com/MOREnvironment/more_generic_models (no licence file),
    revision ``524e336``: ``more_generic_models/dynamics/plant/matrices/
    linear_damping.py``, lines listed above.
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
    smooth_speed_epsilon: float


@dataclass(frozen=True)
class FloatingLinearDampingConstants:
    """``Dmtrx.m`` diagonal of the surface-craft branch."""

    damping_coefficients: np.ndarray


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


def _non_negative(values: np.ndarray, name: str) -> None:
    if np.any(np.asarray(values) < 0.0):
        raise ValueError(f"{name} must be non-negative")


def _relative_speed(nu_r: ca.SX, smooth_speed_epsilon: float) -> ca.SX:
    """``U_r = sqrt(u_r^2 + v_r^2 + w_r^2)`` (MSS ``remus100.m`` 127).

    With ``smooth_speed_epsilon > 0`` its square is added under the root.
    """
    speed_squared = nu_r[0] ** 2 + nu_r[1] ** 2 + nu_r[2] ** 2  # (remus100.m 127)
    if smooth_speed_epsilon > 0.0:  # departure from remus100.m 127 (module docstring)
        speed_squared = speed_squared + smooth_speed_epsilon**2
    return ca.sqrt(speed_squared)


def preprocess_submerged_linear_damping(
    rigid_body_mass_matrix: np.ndarray,
    added_mass_matrix: np.ndarray,
    weight: float,
    center_of_gravity: Sequence[float],
    center_of_buoyancy: Sequence[float],
    time_constants: Sequence[float],
    damping_ratios: Sequence[float],
    sway_damping_fade: bool = False,
    smooth_speed_epsilon: float = 0.0,
) -> SubmergedLinearDampingConstants:
    """MSS ``Dmtrx.m`` 30-49 diagonal, submerged branch (Fossen 2011,
    eqs. 6.76-6.80, p. 125; source ``get_damping_coeff``).

    Source names: ``MRB``, ``MA``, ``W``, ``r_bg``, ``r_bb``, ``T_126``,
    ``zeta_45``. ``T3 = T2``; ``w4``, ``w5`` from the metacentric height
    ``z_g - z_b``. ``sway_damping_fade=True`` is the template departure
    (``D_linear`` fades sway as well as surge); the default fades surge only.
    ``smooth_speed_epsilon`` (m/s, default 0 = MSS exactly) is added in
    quadrature to the speed of the fade (module docstring).
    """
    mass = _matrix(rigid_body_mass_matrix, "rigid_body_mass_matrix") + _matrix(
        added_mass_matrix, "added_mass_matrix"
    )
    r_g = _vector(center_of_gravity, 3, "center_of_gravity")
    r_b = _vector(center_of_buoyancy, 3, "center_of_buoyancy")
    t1, t2, t6 = _vector(time_constants, 3, "time_constants")
    zeta4, zeta5 = _vector(damping_ratios, 2, "damping_ratios")
    _non_negative(np.array([zeta4, zeta5]), "damping_ratios")
    _positive(np.array([t1, t2, t6]), "time_constants")
    _positive(np.diag(mass), "diagonal of the total mass matrix")
    if not np.isfinite(weight) or weight <= 0.0:
        raise ValueError("weight must be a positive finite value")
    metacentric_height = r_g[2] - r_b[2]
    if metacentric_height <= 0.0:
        raise ValueError("center_of_gravity must lie below center_of_buoyancy (z down)")
    if not np.isfinite(smooth_speed_epsilon) or smooth_speed_epsilon < 0.0:
        raise ValueError("smooth_speed_epsilon must be a non-negative finite value")

    t3 = t2  # (Dmtrx.m 44)
    w4 = np.sqrt(weight * metacentric_height / mass[3, 3])  # (Dmtrx.m 45)
    w5 = np.sqrt(weight * metacentric_height / mass[4, 4])  # (Dmtrx.m 46)
    # M_ii / T_i and 2 zeta_i w_i M_ii (Fossen 2011, eqs. 6.76-6.81, p. 125; Dmtrx.m 48-49)
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
        smooth_speed_epsilon=float(smooth_speed_epsilon),
    )


def submerged_linear_damping_casadi(
    constants: SubmergedLinearDampingConstants,
) -> ca.Function:
    """``nu_r -> (D, tau)``; ``D(1,1) *= exp(-3 U_r)`` (MSS ``remus100.m`` 218).

    ``D(2,2)`` fades the same way only with ``sway_damping_fade`` (``D_linear``
    lines 109-110). ``tau = -D nu_r``.
    """
    nu_r = ca.SX.sym("nu_r", 6)
    fade = ca.exp(-3 * _relative_speed(nu_r, constants.smooth_speed_epsilon))  # (remus100.m 218)
    diagonal = [ca.SX(c) for c in constants.damping_coefficients]
    diagonal[0] = diagonal[0] * fade  # (remus100.m 218)
    if constants.sway_damping_fade:  # (source D_linear 109-110)
        diagonal[1] = diagonal[1] * fade
    damping = ca.diag(ca.vertcat(*diagonal))
    return ca.Function(
        "submerged_linear_damping",
        [nu_r],
        [damping, -damping @ nu_r],  # tau = -D nu_r (Fossen 2011, eq. 6.62, p. 123; remus100.m 258)
        ["nu_r"],
        ["D", "tau"],
    )


def preprocess_floating_linear_damping(
    rigid_body_mass_matrix: np.ndarray,
    added_mass_matrix: np.ndarray,
    restoring_matrix: np.ndarray,
    time_constants: Sequence[float],
    damping_ratios: Sequence[float],
) -> FloatingLinearDampingConstants:
    """MSS ``Dmtrx.m`` 51-63, the surface-craft branch (Fossen 2011,
    eqs. 6.76-6.81, p. 125; source ``get_D_linear_surface_vessel``).

    Source names: ``MRB``, ``MA``, ``G`` (only ``G33``, ``G44``, ``G55`` used),
    ``T_126``, ``zeta_345``. ``time_constants = (T1, T2, T6)``;
    ``damping_ratios = (zeta3, zeta4, zeta5)``. MSS fixes ``zeta3`` inside the
    function (``Dmtrx.m`` 57); here it is a parameter.
    """
    mass = _matrix(rigid_body_mass_matrix, "rigid_body_mass_matrix") + _matrix(
        added_mass_matrix, "added_mass_matrix"
    )
    restoring = _matrix(restoring_matrix, "restoring_matrix")
    t1, t2, t6 = _vector(time_constants, 3, "time_constants")
    zeta3, zeta4, zeta5 = _vector(damping_ratios, 3, "damping_ratios")
    _non_negative(np.array([zeta3, zeta4, zeta5]), "damping_ratios")
    _positive(np.array([t1, t2, t6]), "time_constants")
    _positive(np.diag(mass), "diagonal of the total mass matrix")
    stiffness = np.array([restoring[2, 2], restoring[3, 3], restoring[4, 4]])
    if np.any(stiffness < 0.0):
        raise ValueError("G33, G44 and G55 of restoring_matrix must be non-negative")

    # w_i = sqrt(G_ii / M_ii) (Fossen 2011, eqs. 4.51-4.53, p. 68; Dmtrx.m 58-60)
    w3 = np.sqrt(stiffness[0] / mass[2, 2])
    w4 = np.sqrt(stiffness[1] / mass[3, 3])
    w5 = np.sqrt(stiffness[2] / mass[4, 4])
    # (Fossen 2011, eqs. 6.76-6.81, p. 125; Dmtrx.m 62-63)
    coefficients = np.array(
        [
            mass[0, 0] / t1,
            mass[1, 1] / t2,
            mass[2, 2] * 2 * zeta3 * w3,
            mass[3, 3] * 2 * zeta4 * w4,
            mass[4, 4] * 2 * zeta5 * w5,
            mass[5, 5] / t6,
        ]
    )
    return FloatingLinearDampingConstants(damping_coefficients=coefficients)


def floating_linear_damping_casadi(
    constants: FloatingLinearDampingConstants,
) -> ca.Function:
    """``nu_r -> (D, tau)``; ``D`` constant (MSS ``Dmtrx.m`` 62-63), ``tau = -D
    nu_r`` (Fossen 2011, eq. 6.62, p. 123)."""
    nu_r = ca.SX.sym("nu_r", 6)
    damping = ca.SX(ca.DM(np.diag(constants.damping_coefficients)))
    return ca.Function(
        "floating_linear_damping",
        [nu_r],
        [damping, -damping @ nu_r],  # (Fossen 2011, eq. 6.62, p. 123)
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
    """MSS ``otter.m`` 196-206 (source ``D_linear_catamaran``).

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
    _non_negative(np.array([zeta3, zeta4, zeta5]), "damping_ratios")
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

    w3, w4, w5 = np.sqrt(stiffness / np.diag(mass)[2:5])  # (Fossen 2011, eqs. 4.51-4.53, p. 68; otter.m 196-198)
    derivatives = np.array(
        [
            -max_forward_thrust / max_speed,  # Xu (otter.m 201)
            -mass[1, 1] / t_sway,  # Yv (Fossen 2011, eq. 6.77, p. 125; otter.m 202)
            -2 * zeta3 * w3 * mass[2, 2],  # Zw (Fossen 2011, eq. 6.78, p. 125; otter.m 203)
            -2 * zeta4 * w4 * mass[3, 3],  # Kp (Fossen 2011, eq. 6.79, p. 125; otter.m 204)
            -2 * zeta5 * w5 * mass[4, 4],  # Mq (Fossen 2011, eq. 6.80, p. 125; otter.m 205)
            -mass[5, 5] / t_yaw,  # Nr (Fossen 2011, eq. 6.81, p. 125; otter.m 206)
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

    ``tau`` = MSS ``otter.m`` 234-239 (source ``get_tau_damping``): the
    derivatives times ``nu_r``, yaw times ``1 + k |r|`` with ``k =
    yaw_damping_nonlinearity``.
    """
    nu_r = ca.SX.sym("nu_r", 6)
    derivatives = constants.damping_derivatives
    forces = [derivatives[i] * nu_r[i] for i in range(5)]  # (otter.m 234-238)
    forces.append(  # (otter.m 239)
        derivatives[5]
        * (1 + constants.yaw_damping_nonlinearity * ca.fabs(nu_r[5]))
        * nu_r[5]
    )
    damping = ca.SX(ca.DM(-np.diag(derivatives)))  # (Fossen 2011, eq. 6.62, p. 123)
    return ca.Function(
        "surface_linear_damping",
        [nu_r],
        [damping, ca.vertcat(*forces)],
        ["nu_r"],
        ["D", "tau"],
    )
