"""Single-screw propeller with Wageningen B-series coefficients.

Equations (keys in References), ``n`` in rpm saturated at ``n_max``,
``n_p = n / 60`` (MSS ``remus100.m`` 115):

* ``Va = (1 - w) U_r`` (``remus100.m`` 150, there ``Va = 0.944 U_r``);
* ``"bollard"``: ``X = rho D^4 K_T0 n_p |n_p|``, ``K = rho D^5 K_Q0 n_p |n_p|``
  (``remus100.m`` 175-176; Fossen 2011, eq. 12.265, p. 411);
* ``"linearized"``: for ``n_p > 0``, ``X = rho D^4 (K_T0 n_p |n_p| +
  (K_Tmax - K_T0) / Ja_max (Va / D) |n_p|)`` and the same for ``K`` with
  ``D^5`` and ``K_Q`` (``remus100.m`` 163-172); for ``n_p <= 0`` the bollard
  form (``remus100.m`` 173-176);
* ``"full"``: ``J = Va / (|n_p| D)``, ``X = rho D^4 K_T(J) n_p |n_p|``,
  ``K = rho D^5 K_Q(J) n_p |n_p|`` with the Wageningen polynomial
  (``wageningen_kt_kq.py``);
* body wrench: thrust ``(1 - t) X`` along the shaft axis at ``position``
  (moment ``r x f``, Fossen 2011, eq. 12.226, p. 400), roll torque
  ``scale K`` about the shaft axis (``remus100.m`` 249, 252, where the scale
  is ``1/10``).

The ``"linearized"`` model equals MSS ``remus100.m`` 148-176, 249 and 252.
Deviations from MSS: the polynomial coefficients ``K_T0, K_Q0, K_Tmax,
K_Qmax`` are computed from the Wageningen polynomial at ``J = 0`` and
``J = Ja_max`` unless given (``remus100.m`` 156-161 types them, rounded to
four digits); the shaft axis and position are parameters (``remus100.m``
puts the thrust on the x axis at the CO).

Frames: the shaft axis is ``R_zyx(orientation)[:, 0]`` in BODY
(``more_transformations`` ``MatrixTransforms.Rzyx_explicit``; Fossen 2011,
eq. 2.18, p. 22).

Ported from the numpy source [MGM]
``dynamics/propulsion/thruster/thruster_wagenigen/thruster_wagenigen.py``
(``ThrusterWageningen``): coefficients (62-63), ``_compute_kt_kq``
(104-115), ``_compute_linearized_open_water`` (117-137), ``_compute_bollard``
(139-142), ``compute_open_water`` (144-175) and ``compute_force`` (217-268),
with the saturation of ``thruster_base.py``. Motor lag and the thrust-command
input type are states / interfaces, not part of this map.

References
----------
[Fossen 2011] Fossen, T. I. (2011). *Handbook of Marine Craft Hydrodynamics
    and Motion Control*, 1st ed. John Wiley & Sons, Chichester. Ch. 2, eq. 2.18, p. 22; Ch. 12, eq. 12.226, p. 400; eq. 12.265,
    p. 411.
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2.
    https://github.com/cybergalactic/MSS, MIT licence, revision ``72656d1``:
    ``CRAFT/AUV/models/remus100.m`` 110, 115, 127, 148-176, 249, 252.
[MGM] Krizman, E. *more_generic_models*.
    https://github.com/MOREnvironment/more_generic_models (no licence file),
    revision ``524e336``:
    the files and lines listed above.
"""

from dataclasses import dataclass
from typing import Optional, Sequence

import casadi as ca
import numpy as np
from more_transformations.matrix_transforms import MatrixTransforms

from ._common import array, positive, safe_norm, saturate, wrench
from .wageningen_kt_kq import (
    WageningenConstants,
    preprocess_wageningen,
    wageningen_expression,
    wageningen_values,
)

OPEN_WATER_MODELS = ("linearized", "full", "bollard")
MIN_ADVANCE_DENOMINATOR = 1e-6  # thruster_wagenigen.py lines 110, 161


@dataclass(frozen=True)
class PropellerConstants:
    """``thrust_torque_coefficients`` = ``(KT_0, KQ_0, KT_max, KQ_max)`` in use."""

    diameter: float
    max_speed: float
    thrust_deduction: float
    wake_fraction: float
    max_advance_number: float
    roll_moment_scale: float
    water_density: float
    position: np.ndarray
    shaft_axis: np.ndarray
    open_water_model: str
    wageningen: WageningenConstants
    thrust_torque_coefficients: np.ndarray


def preprocess_propeller(
    diameter: float,
    max_speed: float,
    thrust_deduction: float,
    wake_fraction: float,
    pitch_diameter_ratio: float,
    blade_area_ratio: float,
    blade_count: float,
    max_advance_number: float,
    roll_moment_scale: float,
    water_density: float,
    position: Sequence[float],
    orientation: Sequence[float],
    open_water_model: str,
    thrust_torque_coefficients: Optional[Sequence[float]] = None,
    clip_advance_ratio: bool = False,
) -> PropellerConstants:
    """Source names: ``D_prop``, ``n_max`` (rpm), ``t_prop``, ``w_factor``
    (used as ``Va = (1 - w_factor) U``), ``PD``, ``AEAO``, ``z``, ``Ja_max``,
    ``scale_roll``, ``rho``, ``position``, ``orientation_rpy`` (rad), ``mode``.

    ``thrust_torque_coefficients`` overrides the polynomial values at ``J = 0``
    and ``J = max_advance_number`` with exact ones (a coefficient is
    overridden only when an exact value exists).
    ``clip_advance_ratio`` is passed to the Wageningen polynomial (default:
    MSS, no clip; ``True``: the numpy source's ``min(|J|, 1.3)``).
    """
    if open_water_model not in OPEN_WATER_MODELS:
        raise ValueError(f"open_water_model must be one of {OPEN_WATER_MODELS}")
    diameter = positive("diameter", diameter)
    max_speed = positive("max_speed", max_speed)
    max_advance_number = positive("max_advance_number", max_advance_number)
    water_density = positive("water_density", water_density)
    thrust_deduction = float(thrust_deduction)
    wake_fraction = float(wake_fraction)
    roll_moment_scale = float(roll_moment_scale)
    for name, value in (("thrust_deduction", thrust_deduction), ("wake_fraction", wake_fraction),
                        ("roll_moment_scale", roll_moment_scale)):
        if not np.isfinite(value):
            raise ValueError(f"{name} must be a finite value")
    r = array(position, (3,), "position")
    rotation = MatrixTransforms.Rzyx_explicit(array(orientation, (3,), "orientation"))

    wageningen = preprocess_wageningen(
        pitch_diameter_ratio, blade_area_ratio, blade_count, clip_advance_ratio
    )
    if thrust_torque_coefficients is None:
        coefficients = np.array(
            [*wageningen_values(wageningen, 0.0), *wageningen_values(wageningen, max_advance_number)]
        )
    else:
        coefficients = array(thrust_torque_coefficients, (4,), "thrust_torque_coefficients")

    return PropellerConstants(
        diameter=diameter,
        max_speed=max_speed,
        thrust_deduction=thrust_deduction,
        wake_fraction=wake_fraction,
        max_advance_number=max_advance_number,
        roll_moment_scale=roll_moment_scale,
        water_density=water_density,
        position=r,
        shaft_axis=rotation[:, 0].copy(),
        open_water_model=open_water_model,
        wageningen=wageningen,
        thrust_torque_coefficients=coefficients,
    )


def _open_water(c: PropellerConstants, n_rps, advance_speed):
    """Shaft thrust ``X_prop`` and torque ``K_prop`` of the chosen model."""
    kt0, kq0, kt_max, kq_max = (float(v) for v in c.thrust_torque_coefficients)
    rho, d = c.water_density, c.diameter
    quadratic = ca.fabs(n_rps) * n_rps

    if c.open_water_model == "bollard":  # (remus100.m 175-176)
        return rho * d**4 * kt0 * quadratic, rho * d**5 * kq0 * quadratic

    if c.open_water_model == "linearized":  # (remus100.m 163-176)
        kt_slope = (kt_max - kt0) / c.max_advance_number
        kq_slope = (kq_max - kq0) / c.max_advance_number
        inflow = (advance_speed / d) * ca.fabs(n_rps)
        thrust = ca.if_else(  # (remus100.m 169-170, 175)
            n_rps > 0.0,
            rho * d**4 * (kt0 * quadratic + kt_slope * inflow),
            rho * d**4 * kt0 * quadratic,
        )
        torque = ca.if_else(  # (remus100.m 171-172, 176)
            n_rps > 0.0,
            rho * d**5 * (kq0 * quadratic + kq_slope * inflow),
            rho * d**5 * kq0 * quadratic,
        )
        return thrust, torque

    # "full": J clipped to [0, 2 Ja_max] (both paths), then the polynomial's
    # own clip only with ``clip_advance_ratio``;
    # K_T(0), K_Q(0) when n <= 0; zero force at n = 0 through |n| n.
    advance_number = advance_speed / ca.fmax(ca.fabs(n_rps) * d, MIN_ADVANCE_DENOMINATOR)
    advance_number = saturate(advance_number, 0.0, 2.0 * c.max_advance_number)
    kt_j, kq_j = wageningen_expression(c.wageningen, advance_number)
    kt = ca.if_else(n_rps <= 0.0, kt0, kt_j)
    kq = ca.if_else(n_rps <= 0.0, kq0, kq_j)
    return rho * d**4 * kt * quadratic, rho * d**5 * kq * quadratic  # (Fossen 2011, eq. 12.265, p. 411)


def propeller_casadi(constants: PropellerConstants) -> ca.Function:
    """``(n, nu_r) -> tau``; ``n`` (1x1) rpm, saturated at ``max_speed``;
    advance speed ``(1 - w) |nu_r[0:3]|``."""
    n = ca.SX.sym("n", 1)
    nu_r = ca.SX.sym("nu_r", 6)
    n_rps = saturate(n, -constants.max_speed, constants.max_speed) / 60.0  # (remus100.m 115)
    advance_speed = (1.0 - constants.wake_fraction) * safe_norm(nu_r[0:3])  # (remus100.m 127, 150)
    thrust, torque = _open_water(constants, n_rps, advance_speed)

    axis = ca.DM(constants.shaft_axis)
    force = axis * ((1.0 - constants.thrust_deduction) * thrust)  # (remus100.m 249)
    shaft_torque = axis * (constants.roll_moment_scale * torque)  # (remus100.m 252)
    tau = wrench(force, constants.position) + ca.vertcat(ca.DM.zeros(3), shaft_torque)  # (Fossen 2011, eq. 12.226, p. 400)
    return ca.Function("propeller", [n, nu_r], [tau], ["n", "nu_r"], ["tau"])
