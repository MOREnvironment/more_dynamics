"""Steerable electric outboard motor: force map of the rpm and throttle models.

Equations (keys in References):

* rpm model: ``n_p = n / 60``, ``Va = c_Va (1 - w) u_r``, ``J = Va / (n_p D)``,
  ``K_T(J)``, ``K_Q(J)`` from the Wageningen polynomial; ``T = rho D^4 K_T n_p
  |n_p|``, ``Q = rho D^5 K_Q n_p |n_p|`` (Fossen 2011, eq. 12.265, p. 411);
  ``T = (1 - t) T``; shaft power ``P = |2 pi n_p Q|`` capped at ``eta P_max``
  by scaling ``T``; ``T`` saturated at ``+-T_max``;
* throttle model: ``T`` saturated at ``min(T_max, P_prop / max(U_h, 0.1))``;
* steered wrench ``f = T [cos delta, sin delta, 0]`` (an azimuth thruster,
  Fossen 2011, Table 12.3, p. 398), ``tau = [f; r x f]`` (Fossen 2011,
  eq. 12.226, p. 400).

No MSS file holds this model; it is the numpy source's own (lines below).

Only the algebraic force map is here. The rpm lag, the thrust lag, the
steering rate limit and the battery are states of the assembly: the inputs
are the *actual* propeller speed (rpm) or thrust (N) and the *actual*
steering angle ``delta`` (rad, port-positive, the source's internal
convention). A flat battery (zero force in the source) is not modelled.

Reverse thrust: for ``n < 0`` the reverse factors scale ``K_T`` and ``K_Q``
only, ``|n| n`` carries the sign, so full reverse gives a negative surge
force (source lines 317-321 at ``524e336``; an earlier revision of the source
had an extra ``sign(n)`` that pushed forward at full reverse, corrected by
the owner's decision of 2026-10-05).

Ported from the numpy source [MGM]
``dynamics/propulsion/thruster/electrical_outboard_motor/electrical_outboard_motor.py``:
``ElectricOutboardMotor.step`` (throttle, 105-184) and
``ElectricOutboardMotorRPM._prop_thrust_torque`` / ``step`` (rpm, 311-433);
parameters in ``config/dataclass/thruster/electrical_outboard_motor_params.py``.

References
----------
[Fossen 2011] Fossen, T. I. (2011). *Handbook of Marine Craft Hydrodynamics
    and Motion Control*, 1st ed. John Wiley & Sons, Chichester. Ch. 12, Table 12.3, p. 398; eq. 12.226, p. 400; eq. 12.265,
    p. 411.
[MGM] Krizman, E. *more_generic_models*.
    https://github.com/MOREnvironment/more_generic_models (no licence file),
    revision ``524e336``:
    the files and lines listed above.
"""

from dataclasses import dataclass
from typing import Sequence

import casadi as ca
import numpy as np

from ._common import array, finite, positive, safe_norm, saturate, wrench
from .wageningen_kt_kq import WageningenConstants, preprocess_wageningen, wageningen_expression

ZERO_SPEED_RPS = 1e-9         # electrical_outboard_motor.py line 312
MIN_EFFICIENCY = 1e-3         # lines 391, 176
MIN_HORIZONTAL_SPEED = 0.1    # line 157, m/s


@dataclass(frozen=True)
class OutboardMotorRpmConstants:
    max_thrust: float
    shaft_power_limit: float
    position: np.ndarray
    max_speed: float
    propeller_diameter: float
    thrust_deduction: float
    water_density: float
    advance_speed_gain: float
    reverse_thrust_factor: float
    reverse_torque_factor: float
    power_limit: bool
    wageningen: WageningenConstants


@dataclass(frozen=True)
class OutboardMotorThrottleConstants:
    max_thrust: float
    max_propulsive_power: float
    position: np.ndarray


def _steered_wrench(thrust, delta, position: np.ndarray) -> ca.SX:
    """``f = T [cos delta, sin delta, 0]``, ``tau = [f; r x f]`` (Fossen 2011,
    eq. 12.226, p. 400; source lines 162-171, 403-414)."""
    force = ca.vertcat(thrust * ca.cos(delta), thrust * ca.sin(delta), 0.0)
    return wrench(force, position)


def preprocess_outboard_motor_rpm(
    max_thrust: float,
    max_power: float,
    efficiency: float,
    position: Sequence[float],
    max_speed: float,
    propeller_diameter: float,
    pitch_diameter_ratio: float,
    blade_area_ratio: float,
    blade_count: float,
    thrust_deduction: float,
    water_density: float,
    advance_speed_factor: float,
    propwash_factor: float,
    reverse_thrust_factor: float,
    reverse_torque_factor: float,
    power_limit: bool = True,
    clip_advance_ratio: bool = False,
) -> OutboardMotorRpmConstants:
    """Source names: ``max_thrust`` (N), ``max_power`` (W, electrical input),
    ``efficiency``, ``r_body``, ``n_max_rpm``, ``D_prop``, ``PD``, ``AEAO``,
    ``z``, ``t_prop``, ``rho``, ``Va_factor``, ``w_prop``,
    ``reverse_KT_factor``, ``reverse_KQ_factor``; ``power_limit`` is
    ``not disable_power_limit``. Shaft power cap ``max_power max(eta, 1e-3)``;
    advance speed ``Va_factor (1 - w_prop) u_r``. ``clip_advance_ratio`` is
    passed to the Wageningen polynomial (default MSS, no clip)."""
    return OutboardMotorRpmConstants(
        max_thrust=positive("max_thrust", max_thrust),
        shaft_power_limit=positive("max_power", max_power)
        * max(finite("efficiency", efficiency), MIN_EFFICIENCY),
        position=array(position, (3,), "position"),
        max_speed=positive("max_speed", max_speed),
        propeller_diameter=positive("propeller_diameter", propeller_diameter),
        thrust_deduction=finite("thrust_deduction", thrust_deduction),
        water_density=positive("water_density", water_density),
        advance_speed_gain=finite("advance_speed_factor", advance_speed_factor)
        * (1.0 - finite("propwash_factor", propwash_factor)),
        reverse_thrust_factor=finite("reverse_thrust_factor", reverse_thrust_factor),
        reverse_torque_factor=finite("reverse_torque_factor", reverse_torque_factor),
        power_limit=bool(power_limit),
        wageningen=preprocess_wageningen(
            pitch_diameter_ratio, blade_area_ratio, blade_count, clip_advance_ratio
        ),
    )


def outboard_motor_rpm_casadi(constants: OutboardMotorRpmConstants) -> ca.Function:
    """``(n, delta, nu_r) -> tau``; ``n`` (1x1) rpm, saturated at ``max_speed``."""
    n = ca.SX.sym("n", 1)
    delta = ca.SX.sym("delta", 1)
    nu_r = ca.SX.sym("nu_r", 6)
    c = constants

    n_rps = saturate(n, -c.max_speed, c.max_speed) / 60.0
    advance_speed = c.advance_speed_gain * nu_r[0]  # Va (source line 382)
    stopped = ca.fabs(n_rps) < ZERO_SPEED_RPS
    advance_number = advance_speed / (ca.if_else(stopped, 1.0, n_rps) * c.propeller_diameter)
    kt, kq = wageningen_expression(c.wageningen, advance_number)
    reverse = n_rps < 0.0
    kt = ca.if_else(reverse, c.reverse_thrust_factor * kt, kt)
    kq = ca.if_else(reverse, c.reverse_torque_factor * kq, kq)

    quadratic = ca.fabs(n_rps) * n_rps
    # (Fossen 2011, eq. 12.265, p. 411; source _prop_thrust_torque 311-325)
    thrust = ca.if_else(stopped, 0.0, c.water_density * c.propeller_diameter**4 * kt * quadratic)
    torque = ca.if_else(stopped, 0.0, c.water_density * c.propeller_diameter**5 * kq * quadratic)
    thrust = (1.0 - c.thrust_deduction) * thrust  # (source line 385)

    if c.power_limit:
        shaft_power = ca.fabs(2.0 * np.pi * n_rps * torque)  # P = 2 pi n Q (source lines 388-395)
        limit = c.shaft_power_limit
        thrust = ca.if_else(shaft_power > limit, thrust * (limit / ca.fmax(shaft_power, limit)), thrust)

    thrust = saturate(thrust, -c.max_thrust, c.max_thrust)
    tau = _steered_wrench(thrust, delta, c.position)
    return ca.Function("outboard_motor_rpm", [n, delta, nu_r], [tau], ["n", "delta", "nu_r"], ["tau"])


def preprocess_outboard_motor_throttle(
    max_thrust: float,
    max_propulsive_power: float,
    position: Sequence[float],
) -> OutboardMotorThrottleConstants:
    """Source names ``max_thrust`` (N), ``max_propulsive_power`` (W delivered
    to the water), ``r_body``. Thrust cap ``min(max_thrust, P / max(U_h, 0.1))``."""
    return OutboardMotorThrottleConstants(
        max_thrust=positive("max_thrust", max_thrust),
        max_propulsive_power=positive("max_propulsive_power", max_propulsive_power),
        position=array(position, (3,), "position"),
    )


def outboard_motor_throttle_casadi(constants: OutboardMotorThrottleConstants) -> ca.Function:
    """``(thrust, delta, nu_r) -> tau``; ``thrust`` (1x1) the actual thrust state in N."""
    thrust = ca.SX.sym("thrust", 1)
    delta = ca.SX.sym("delta", 1)
    nu_r = ca.SX.sym("nu_r", 6)
    c = constants
    horizontal_speed = ca.fmax(safe_norm(nu_r[0:2]), MIN_HORIZONTAL_SPEED)
    limit = ca.fmin(c.max_thrust, c.max_propulsive_power / horizontal_speed)
    tau = _steered_wrench(saturate(thrust, -limit, limit), delta, c.position)
    return ca.Function(
        "outboard_motor_throttle", [thrust, delta, nu_r], [tau], ["thrust", "delta", "nu_r"], ["tau"]
    )
