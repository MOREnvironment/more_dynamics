"""The ``current`` slot: a uniform current fixed in NED, rotated into BODY,
and its body-frame acceleration (needed for the relative-velocity equation
of motion).  Two MSS forms, both kept: neither is the physics default (the form with
the full attitude and rotation rate awaits a derivation from a read source
before either becomes a default).

References
----------
[MSS] Fossen, T. I. MSS, MIT, @ cc07579: CRAFT/AUV/models/remus100.m 118-125
    (current rotated by yaw only; body acceleration keeps the yaw-rate terms
    only); CRAFT/USV/models/otter.m 113-118 (body acceleration from the full
    body-rate skew matrix).

Author:    Enio Krizman
Date:      2026-10-08
"""

import casadi as ca

from more_transformations.more_casadi_transformations import MatrixTransforms, Parameter

from more_dynamics.models.wiring import function_from

NO_CURRENT_PARAMETERS = ()
HORIZONTAL_CURRENT_PARAMETERS = (
    Parameter("current_speed", (1, 1), "m/s", "current speed Vc (NED, horizontal)", 0.0),
    Parameter("current_direction", (1, 1), "rad", "current set (NED, from-north, clockwise)"),
    Parameter("current_vertical_speed", (1, 1), "m/s", "current w_c (NED, down positive)"),
)


def no_current():
    """No current: zero velocity and acceleration."""
    return function_from("no_current", {}, {
        "current_velocity": ca.SX.zeros(6), "current_acceleration": ca.SX.zeros(6)})


def _horizontal_current(name, full_rotation_rate):
    eta = ca.SX.sym("eta", 6)
    nu = ca.SX.sym("nu", 6)
    current_speed = ca.SX.sym("current_speed")
    current_direction = ca.SX.sym("current_direction")
    current_vertical_speed = ca.SX.sym("current_vertical_speed")
    psi = eta[5]
    u_c = current_speed * ca.cos(current_direction - psi)  # (remus100.m 118; otter.m 113)
    v_c = current_speed * ca.sin(current_direction - psi)  # (remus100.m 119; otter.m 114)
    nu_c = ca.vertcat(u_c, v_c, current_vertical_speed, 0, 0, 0)  # (remus100.m 121; otter.m 115 with w_c = 0)
    if full_rotation_rate:
        # nu_c_dot = [-S(omega) nu_c(1:3); 0] (otter.m 117-118)
        nu_c_dot = ca.vertcat(-MatrixTransforms.skew(nu[3:6]) @ nu_c[0:3], ca.SX.zeros(3))
    else:
        # yaw-rate terms only: [r v_c, -r u_c, 0, 0, 0, 0] (remus100.m 122)
        nu_c_dot = ca.vertcat(nu[5] * v_c, -nu[5] * u_c, 0, 0, 0, 0)
    return function_from(name, {"eta": eta, "nu": nu, "current_speed": current_speed,
                                "current_direction": current_direction,
                                "current_vertical_speed": current_vertical_speed},
                         {"current_velocity": nu_c, "current_acceleration": nu_c_dot})


def horizontal_current_yaw_rate_terms():
    """MSS form: the current's body acceleration keeps only the yaw-rate
    terms (remus100.m 118-122)."""
    return _horizontal_current("horizontal_current_yaw_rate_terms", False)


def horizontal_current_full_rotation_rate():
    """MSS form: the current's body acceleration uses the full body-rate
    skew matrix (otter.m 113-118)."""
    return _horizontal_current("horizontal_current_full_rotation_rate", True)
