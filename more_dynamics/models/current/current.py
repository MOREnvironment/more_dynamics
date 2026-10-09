"""The ``current`` slot: a uniform current fixed in NED, rotated into BODY,
and its body-frame acceleration (needed for the relative-velocity equation
of motion). Three forms:

* ``uniform_current_full_attitude`` (the physics form, the default of the
  vehicle types): ``v_c^b = R(Theta)^T v_c^n`` with the full attitude
  ``R = R_z,psi R_y,theta R_x,phi`` and ``nu_c_dot = [-S(omega) v_c^b; 0]``
  for a current constant in NED (Fossen 2011, eqs. 8.138-8.141 and 8.157,
  pp. 221-225; eq. 2.33, Sec. 2.2, pp. 20-27). Exact for every attitude and rotation rate.
  Raise: a measured current profile (a current that varies in time adds
  ``R^T dv_c^n/dt``, not modelled here).
* ``horizontal_current_yaw_rate_terms`` (MSS ``remus100.m``): the current is
  rotated by the heading only and its acceleration keeps the yaw-rate terms
  only; exact for a level craft without roll and pitch rate.
* ``horizontal_current_full_rotation_rate`` (MSS ``otter.m``): the current is
  rotated by the heading only and its acceleration is the full body-rate
  skew matrix on that vector; it adds a heave term its velocity does not
  have when roll or pitch is non-zero. Both MSS forms are shortcuts kept to
  reproduce MSS; the first is exact only at ``phi = theta = 0`` and
  ``p = q = 0``, the second only at ``phi = theta = 0``.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics
    and Motion Control. Wiley. Sec. 8.3, eqs. 8.137-8.141, 8.155-8.157,
    pp. 221-225; Sec. 2.2, eq. 2.18, p. 22 and Theorem 2.2, eq. 2.33, pp. 20-27
    (R-dot = R S(omega)); the text is Sec. 8.3 pp. 221-225 and Sec. 2.2
    pp. 20-27 (printed pages).
[MSS] Fossen, T. I. MSS, MIT, @ cc07579: CRAFT/AUV/models/remus100.m 118-125
    (current rotated by yaw only; body acceleration keeps the yaw-rate terms
    only); CRAFT/USV/models/otter.m 113-118 (body acceleration from the full
    body-rate skew matrix).

Author:    Enio Krizman
Date:      2026-10-08
"""

import casadi as ca

from more_transformations.more_casadi_transformations import MatrixTransforms, Parameter

from more_dynamics.models.shared.wiring import function_from

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


def uniform_current_full_attitude():
    """``(eta, nu, current_speed, current_direction, current_vertical_speed) ->
    (current_velocity, current_acceleration)``: a current constant in NED,
    ``v_c^n = [V_c cos(beta_c), V_c sin(beta_c), w_c]`` (``remus100.m`` 118-119
    with the heading removed; Fossen 2011, eq. 8.156, pp. 221-225), taken into BODY
    with the full attitude, ``v_c^b = R(Theta)^T v_c^n`` (eq. 8.157), and
    ``nu_c = [v_c^b; 0]`` (eqs. 8.137, 8.155), ``nu_c_dot = [-S(omega) v_c^b; 0]``
    (eq. 8.141, from ``R-dot = R S(omega)``, eq. 2.33, Sec. 2.2, pp. 20-27)."""
    eta = ca.SX.sym("eta", 6)
    nu = ca.SX.sym("nu", 6)
    current_speed = ca.SX.sym("current_speed")
    current_direction = ca.SX.sym("current_direction")
    current_vertical_speed = ca.SX.sym("current_vertical_speed")
    v_ned = ca.vertcat(current_speed * ca.cos(current_direction),  # (remus100.m 118: V_c cos(beta_c))
                       current_speed * ca.sin(current_direction),  # (remus100.m 119: V_c sin(beta_c))
                       current_vertical_speed)  # (otter.m 115: w_c, positive down)
    v_body = MatrixTransforms.Rzyx(eta[3:6]).T @ v_ned  # v_c^b = R^T v_c^n (Fossen 2011, eq. 8.157, pp. 221-225)
    nu_c = ca.vertcat(v_body, ca.SX.zeros(3))  # (Fossen 2011, eqs. 8.137, 8.155)
    nu_c_dot = ca.vertcat(-MatrixTransforms.skew(nu[3:6]) @ v_body, ca.SX.zeros(3))  # (Fossen 2011, eq. 8.141, pp. 221-225)
    return function_from("uniform_current_full_attitude", {
        "eta": eta, "nu": nu, "current_speed": current_speed, "current_direction": current_direction,
        "current_vertical_speed": current_vertical_speed},
        {"current_velocity": nu_c, "current_acceleration": nu_c_dot})
