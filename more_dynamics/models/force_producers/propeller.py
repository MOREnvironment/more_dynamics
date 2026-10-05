"""Single-screw propeller with Wageningen B-series coefficients.

Ported from ``more_generic_models``
``dynamics/propulsion/thruster/thruster_wagenigen/thruster_wagenigen.py``
(``ThrusterWageningen``): coefficients (lines 62-63), ``_compute_kt_kq``
(104-115), ``_compute_linearized_open_water`` (117-137), ``_compute_bollard``
(139-142), ``compute_open_water`` (144-175) and ``compute_force`` (217-268),
with the saturation of ``thruster_base.py``. The ``"linearized"`` model equals
MSS ``remus100.m`` 149-178, 247 and 250. Motor lag and the thrust-command
input type are states / interfaces, not part of this map.

Frames: the shaft axis is ``R_zyx(orientation)[:, 0]`` in BODY; thrust acts
along it at ``position``; the roll torque acts about it.
"""

from dataclasses import dataclass
from typing import Optional, Sequence

import casadi as ca
import numpy as np

from ._common import array, positive, rotation_zyx, safe_norm, saturate
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
    and ``J = max_advance_number`` with exact ones (E-16 parameter style).
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
    rotation = rotation_zyx(array(orientation, (3,), "orientation"))

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

    if c.open_water_model == "bollard":
        return rho * d**4 * kt0 * quadratic, rho * d**5 * kq0 * quadratic

    if c.open_water_model == "linearized":
        kt_slope = (kt_max - kt0) / c.max_advance_number
        kq_slope = (kq_max - kq0) / c.max_advance_number
        inflow = (advance_speed / d) * ca.fabs(n_rps)
        thrust = ca.if_else(
            n_rps > 0.0,
            rho * d**4 * (kt0 * quadratic + kt_slope * inflow),
            rho * d**4 * kt0 * quadratic,
        )
        torque = ca.if_else(
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
    return rho * d**4 * kt * quadratic, rho * d**5 * kq * quadratic


def propeller_casadi(constants: PropellerConstants) -> ca.Function:
    """``(n, nu_r) -> tau``; ``n`` (1x1) rpm, saturated at ``max_speed``;
    advance speed ``(1 - w) |nu_r[0:3]|``."""
    n = ca.SX.sym("n", 1)
    nu_r = ca.SX.sym("nu_r", 6)
    n_rps = saturate(n, -constants.max_speed, constants.max_speed) / 60.0
    advance_speed = (1.0 - constants.wake_fraction) * safe_norm(nu_r[0:3])
    thrust, torque = _open_water(constants, n_rps, advance_speed)

    axis = ca.DM(constants.shaft_axis)
    force = axis * ((1.0 - constants.thrust_deduction) * thrust)
    moment = ca.cross(ca.DM(constants.position), force) + axis * (constants.roll_moment_scale * torque)
    tau = ca.vertcat(force, moment)
    return ca.Function("propeller", [n, nu_r], [tau], ["n", "nu_r"], ["tau"])
