"""Steerable electric outboard motor: force map of the rpm and throttle models (CasADi).

Two blocks: ``outboard_motor_rpm_parameters()`` / ``outboard_motor_rpm_casadi
(...)`` build ``(n, delta, nu_r, <parameters by name>) -> tau``;
``outboard_motor_throttle_parameters()`` / ``outboard_motor_throttle_casadi()``
build ``(thrust, delta, nu_r, <parameters by name>) -> tau``.

Equations (keys in References):

* rpm model: ``n_p = n / 60``, ``Va = c_Va (1 - w) u_r``, ``J = Va / (n_p D)``
  (MSS ``wageningen.m`` 6), ``K_T(J)``, ``K_Q(J)`` from the Wageningen
  polynomial; ``T = rho D^4 K_T n_p |n_p|``, ``Q = rho D^5 K_Q n_p |n_p|``
  (Fossen 2011, eq. 12.265, p. 411); ``T = (1 - t) T`` (the thrust
  deduction of MSS ``remus100.m`` 250); shaft power ``P = |2 pi n_p Q|``
  capped at ``eta P_max`` by scaling ``T``; ``T`` saturated at ``+-T_max``;
* throttle model: ``T`` saturated at ``min(T_max, P_prop / max(U_h, 0.1))``
  (the propulsive power ``P = T U`` at the horizontal speed ``U_h``, bounded
  below at 0.1 m/s);
* steered wrench ``f = T [cos delta, sin delta, 0]`` (an azimuth thruster,
  Fossen 2011, Table 12.3, p. 398, and p. 400 below eq. 12.226),
  ``tau = [f; r x f]`` (Fossen 2011, eq. 12.226, p. 400).

No MSS file holds this motor; the power cap, the throttle cap, the reverse
factors, the efficiency floor and the stop guard are constructions of this
module (stated beside each line), the propeller laws above are published.

Only the algebraic force map is here. The rpm lag, the thrust lag, the
steering rate limit and the battery are states of the assembly: the inputs
are the *actual* propeller speed (rpm) or thrust (N) and the *actual*
steering angle ``delta`` (rad, port-positive: the force is ``T [cos delta,
sin delta, 0]`` with +y to starboard, so a positive angle turns the bow to port
for a motor aft of the CO). A flat battery (zero force)
is not modelled.

Admissible ranges of the dimensionless numbers of the rpm model (refused
outside them by ``check_values``, naming the parameter): ``efficiency`` in
``(0, 1]`` (a drive cannot deliver more shaft power than it draws; zero
would make the cap zero); ``thrust_deduction`` ``t`` and ``propwash_factor``
``w`` in ``[0, 1]`` (fractions of the thrust and of the inflow speed; the
equations need no exclusive endpoint: ``t = 1`` gives zero thrust, ``w = 1``
zero advance speed, ``J = 0``); ``reverse_thrust_factor`` and
``reverse_torque_factor`` ``> 0`` (a factor that keeps the sign of ``K_T``
and ``K_Q``; a negative one would turn reverse thrust into forward thrust).
The efficiency floor ``max(eta, 1e-3)`` in the power cap stays: the
declared range ``eta > 0`` keeps the cap positive without it, but values
``0 < eta < 1e-3`` still reach it, and removing it would change the force
there.

Reverse thrust: for ``n < 0`` the reverse factors scale ``K_T`` and ``K_Q``
only, ``|n| n`` carries the sign, so full reverse gives a negative surge
force (owner's decision of 2026-10-05: an earlier form with an extra
``sign(n)`` pushed forward at full reverse).

References
----------
[Fossen 2011] Fossen, T. I. (2011). *Handbook of Marine Craft Hydrodynamics
    and Motion Control*, 1st ed. John Wiley & Sons, Chichester. Ch. 12,
    Table 12.3, p. 398; eq. 12.226, p. 400; eq. 12.265, p. 411.
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2
    with the fixes of 2026-10-07. https://github.com/cybergalactic/MSS, MIT
    licence, revision ``6e78e98``: ``LIBRARY/modeling/wageningen.m`` 6;
    ``CRAFT/AUV/models/remus100.m`` 250.
[Krizman] Krizman, E. (2026). Model file ``grethe.m @ 679a6e8 (retired 2026-10-30; file deleted, kept as tag archive/grethe-vehicle)``, MSS fork branch ``grethe-vehicle``, revision
    ``679a6e8``: ``CRAFT/USV/models/grethe.m @ 679a6e8 (retired 2026-10-30; file deleted, kept as tag archive/grethe-vehicle)`` 269-271, 400-408 (the leg normal force, ``outboard_leg_normal_force``;
    the file may still be revised). Whicker, L. F. and Fehlner, L. F. (1958), the rudder lift slope the file adapts
    (page not read).

Author:    Enio Krizman
Date:      2026-10-05
"""

import math

import casadi as ca
from more_transformations.more_casadi_transformations import Parameter, symbols

from more_dynamics.models.shared.force_producer_common import safe_norm, saturate, wrench
from more_dynamics.models.propeller.wageningen_kt_kq import wageningen_expression, wageningen_parameters

ZERO_SPEED_RPS = 1e-9         # rev/s: below it the propeller is stopped (guard of this module)
MIN_EFFICIENCY = 1e-3         # floor of the efficiency in the power cap (guard of this module)
MIN_HORIZONTAL_SPEED = 0.1    # m/s: floor of U_h in the throttle cap (guard of this module)
LEG_SPEED_GUARD = 1e-6        # m^2/s^2: added under the root of the flow speed at the leg (guard of the model file)


def _positive(name, unit, meaning):
    return Parameter(name, (1, 1), unit, meaning, 0.0, minimum_exclusive=True)


def _scalar(name, unit, meaning):
    return Parameter(name, (1, 1), unit, meaning)


def _fraction(name, meaning):
    return Parameter(name, (1, 1), "1", meaning, 0.0, 1.0)


_POSITION = Parameter("position", (3, 1), "m", "CO -> motor, body axes (FRD)")

_RPM_DECLARED = (
    _positive("max_thrust", "N", "thrust limit T_max"),
    _positive("max_power", "W", "electrical input power limit P_max"),
    Parameter("efficiency", (1, 1), "1", "drive efficiency eta (shaft power cap eta P_max)", 0.0, 1.0,
              minimum_exclusive=True),
    _POSITION,
    _positive("max_speed", "rpm", "propeller speed limit"),
    _positive("propeller_diameter", "m", "propeller diameter D"),
    *wageningen_parameters(),
    _fraction("thrust_deduction", "thrust deduction t of (1 - t) T"),
    _positive("water_density", "kg/m^3", "water density rho"),
    _scalar("advance_speed_factor", "1", "c_Va of Va = c_Va (1 - w) u_r"),
    _fraction("propwash_factor", "w of Va = c_Va (1 - w) u_r"),
    _positive("reverse_thrust_factor", "1", "factor on K_T for n < 0"),
    _positive("reverse_torque_factor", "1", "factor on K_Q for n < 0"),
)

_THROTTLE_DECLARED = (
    _positive("max_thrust", "N", "thrust limit T_max"),
    _positive("max_propulsive_power", "W", "power delivered to the water P_prop"),
    _POSITION,
)


def outboard_motor_rpm_parameters():
    """The declared parameters of the rpm model, in the order of its inputs."""
    return _RPM_DECLARED


def outboard_motor_throttle_parameters():
    """The declared parameters of the throttle model, in the order of its inputs."""
    return _THROTTLE_DECLARED


def _steered_wrench(thrust, delta, position):
    """``f = T [cos delta, sin delta, 0]``, ``tau = [f; r x f]`` (Fossen 2011,
    Table 12.3, p. 398; eq. 12.226, p. 400)."""
    force = ca.vertcat(thrust * ca.cos(delta), thrust * ca.sin(delta), 0.0)  # azimuth thruster (Fossen 2011, p. 400)
    return wrench(force, position)  # (Fossen 2011, eq. 12.226, p. 400)


def signed_advance_open_water(density, diameter, pitch_diameter_ratio, blade_area_ratio, blade_count,
                              reverse_thrust_factor, reverse_torque_factor, n_rps, advance_speed,
                              clip_advance_ratio=False):
    """Shaft thrust and torque ``(T, Q)`` of the open-water law with a signed advance number: ``J = Va / (n D)``
    (MSS ``wageningen.m`` 6), ``K_T`` and ``K_Q`` from the Wageningen polynomial at that ``J``, the reverse
    factors on both for ``n < 0`` and zero force when the shaft is stopped (constructions of this module);
    ``T = rho D^4 K_T n |n|``, ``Q = rho D^5 K_Q n |n|`` (Fossen 2011, eq. 12.265, p. 411); ``n_rps`` in rev/s."""
    d, rho = diameter, density
    stopped = ca.fabs(n_rps) < ZERO_SPEED_RPS  # guard of this module
    advance_number = advance_speed / (ca.if_else(stopped, 1.0, n_rps) * d)  # J = Va / (n D) (wageningen.m 6), signed
    kt, kq = wageningen_expression(advance_number, pitch_diameter_ratio, blade_area_ratio, blade_count,
                                   clip_advance_ratio)
    reverse = n_rps < 0.0
    kt = ca.if_else(reverse, reverse_thrust_factor * kt, kt)  # reverse factor on K_T: construction
    kq = ca.if_else(reverse, reverse_torque_factor * kq, kq)  # reverse factor on K_Q: construction
    quadratic = ca.fabs(n_rps) * n_rps  # n_p |n_p| (Fossen 2011, eq. 12.265, p. 411)
    thrust = ca.if_else(stopped, 0.0, rho * d**4 * kt * quadratic)  # T = rho D^4 K_T n|n| (Fossen 2011, eq. 12.265, p. 411)
    torque = ca.if_else(stopped, 0.0, rho * d**5 * kq * quadratic)  # Q = rho D^5 K_Q n|n| (Fossen 2011, eq. 12.265 form)
    return thrust, torque


def outboard_leg_normal_force(nu_r, delta, position, density, leg_area, leg_aspect_ratio):
    """Wrench ``tau`` (6x1, BODY, about the CO) of the normal force of the outboard leg: flat-plate-like, the
    force on the leg proportional to the flow component normal to the motor axis and to the flow speed at the leg.

    Flow at the leg ``u_l = u_r - y_P r_r``, ``v_l = v_r + x_P r_r``; unit normal to the motor axis ``n = [-sin d,
    cos d]``; ``v_n = u_l n_x + v_l n_y``; ``U_l = sqrt(u_l^2 + v_l^2 + 1e-6)``; ``F_n = -1/2 rho A C_L' U_l v_n``
    with the slope ``C_L' = 6.13 AR / (AR + 2.25)``; the force ``F_n n`` acts at the motor position (the
    model file ``grethe.m @ 679a6e8 (retired 2026-10-30; file deleted, kept as tag archive/grethe-vehicle)`` 400-408 with the slope of 269-271, in the author's MSS fork at 679a6e8: the slope is
    the rudder lift slope of Whicker and Fehlner adapted to the leg, a construction; the ``1e-6`` under the root
    is the model file's guard against a zero flow). Wrench of a point force ``[f; r x f]`` (Fossen 2011, eq.
    12.226, p. 400)."""
    slope = 6.13 * leg_aspect_ratio / (leg_aspect_ratio + 2.25)  # normal-force slope per rad (grethe.m @ 679a6e8 (retired 2026-10-30; file deleted, kept as tag archive/grethe-vehicle) 271)
    flow_x = nu_r[0] - position[1] * nu_r[5]  # flow at the leg (grethe.m @ 679a6e8 (retired 2026-10-30; file deleted, kept as tag archive/grethe-vehicle) 401)
    flow_y = nu_r[1] + position[0] * nu_r[5]  # (grethe.m @ 679a6e8 (retired 2026-10-30; file deleted, kept as tag archive/grethe-vehicle) 402)
    normal_x, normal_y = -ca.sin(delta), ca.cos(delta)  # normal to the motor axis (grethe.m @ 679a6e8 (retired 2026-10-30; file deleted, kept as tag archive/grethe-vehicle) 403)
    normal_speed = flow_x * normal_x + flow_y * normal_y  # (grethe.m @ 679a6e8 (retired 2026-10-30; file deleted, kept as tag archive/grethe-vehicle) 404)
    speed = ca.sqrt(flow_x**2 + flow_y**2 + LEG_SPEED_GUARD)  # (grethe.m @ 679a6e8 (retired 2026-10-30; file deleted, kept as tag archive/grethe-vehicle) 405)
    normal_force = -0.5 * density * leg_area * slope * speed * normal_speed  # (grethe.m @ 679a6e8 (retired 2026-10-30; file deleted, kept as tag archive/grethe-vehicle) 406)
    return wrench(ca.vertcat(normal_force * normal_x, normal_force * normal_y, 0.0), position)  # (grethe.m @ 679a6e8 (retired 2026-10-30; file deleted, kept as tag archive/grethe-vehicle) 407-408; Fossen 2011, eq. 12.226, p. 400)


def outboard_motor_rpm_casadi(*, power_limit=True, clip_advance_ratio=False):
    """The rpm block: ``(n, delta, nu_r, <parameters by name>) -> tau``.

    Contract
    --------
    ``power_limit`` and ``clip_advance_ratio`` are bools (the shaft-power cap
    on or off; the Wageningen clip). Inputs: ``n`` (1x1, actual propeller
    speed in rpm, saturated at ``max_speed``), ``delta`` (1x1, actual
    steering angle, rad, port-positive: positive turns the bow to port for a motor
    aft of the CO), ``nu_r`` (6x1) and the parameters
    of ``outboard_motor_rpm_parameters()`` by name. Output: ``tau`` (6x1).
    """
    if not isinstance(power_limit, bool) or not isinstance(clip_advance_ratio, bool):
        raise ValueError("power_limit and clip_advance_ratio must be True or False")
    p = symbols(_RPM_DECLARED)
    n = ca.SX.sym("n", 1)
    delta = ca.SX.sym("delta", 1)
    nu_r = ca.SX.sym("nu_r", 6)
    d, rho = p["propeller_diameter"], p["water_density"]

    n_rps = saturate(n, -p["max_speed"], p["max_speed"]) / 60.0  # rpm -> rev/s, saturated
    advance_speed = p["advance_speed_factor"] * (1.0 - p["propwash_factor"]) * nu_r[0]  # Va = c_Va (1 - w) u_r: construction
    thrust, torque = signed_advance_open_water(
        rho, d, p["pitch_diameter_ratio"], p["blade_area_ratio"], p["blade_count"], p["reverse_thrust_factor"],
        p["reverse_torque_factor"], n_rps, advance_speed, clip_advance_ratio)
    thrust = (1.0 - p["thrust_deduction"]) * thrust  # (1 - t) T (remus100.m 250)

    if power_limit:
        shaft_power = ca.fabs(2.0 * math.pi * n_rps * torque)  # P = 2 pi n Q (definition of shaft power)
        limit = p["max_power"] * ca.fmax(p["efficiency"], MIN_EFFICIENCY)  # eta P_max, eta floored (reached for 0 < eta < 1e-3): construction
        # T scaled by limit / P when P exceeds the limit: construction of this module
        thrust = ca.if_else(shaft_power > limit, thrust * (limit / ca.fmax(shaft_power, limit)), thrust)

    thrust = saturate(thrust, -p["max_thrust"], p["max_thrust"])  # thrust cap +-T_max: construction of this module
    tau = _steered_wrench(thrust, delta, p["position"])
    names = [d_.name for d_ in _RPM_DECLARED]
    return ca.Function("outboard_motor_rpm", [n, delta, nu_r, *[p[name] for name in names]], [tau],
                       ["n", "delta", "nu_r", *names], ["tau"])


def outboard_motor_throttle_casadi():
    """The throttle block: ``(thrust, delta, nu_r, max_thrust,
    max_propulsive_power, position) -> tau``; ``thrust`` (1x1) the actual
    thrust state in N, saturated at ``min(T_max, P_prop / max(U_h, 0.1))``."""
    p = symbols(_THROTTLE_DECLARED)
    thrust = ca.SX.sym("thrust", 1)
    delta = ca.SX.sym("delta", 1)
    nu_r = ca.SX.sym("nu_r", 6)
    horizontal_speed = ca.fmax(safe_norm(nu_r[0:2]), MIN_HORIZONTAL_SPEED)  # U_h = |[u, v]|, floored (guard)
    limit = ca.fmin(p["max_thrust"], p["max_propulsive_power"] / horizontal_speed)  # T = P / U (power = thrust x speed): construction
    tau = _steered_wrench(saturate(thrust, -limit, limit), delta, p["position"])
    names = [d.name for d in _THROTTLE_DECLARED]
    return ca.Function("outboard_motor_throttle", [thrust, delta, nu_r, *[p[name] for name in names]], [tau],
                       ["thrust", "delta", "nu_r", *names], ["tau"])
