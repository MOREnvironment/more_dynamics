"""Lift and drag of a slender body (CasADi).

``lift_drag_parameters(smooth_speed=...)`` declares the numbers of the block
(name, shape, SI unit, meaning, admissible range): five, and
``smooth_speed_epsilon`` with ``smooth_speed=True``; ``lift_drag_casadi(...)``
builds ``(nu_r, <parameters by name>) -> tau``. No numbers live in the module.
The keyword ``smooth_speed`` is a selector (it decides the structure of the
graph); every number in an equation, the regularisation ``eps`` included, is a
declared parameter and a named input.

Equations (keys in References):

* ``AR = b^2 / S``, ``C_L_alpha = pi AR / (1 + sqrt(1 + (AR/2)^2))``,
  ``C_L = C_L_alpha alpha`` (linear lift, no stall blending: ``sigma = 0``),
  ``C_D = C_D0 + C_L^2 / (pi e AR)`` (MSS ``coeffLiftDrag.m`` 54-68, which
  attributes them to Beard and McLain (2012), not read here).
* ``F_drag = 1/2 rho U_r^2 S C_D``, ``F_lift = 1/2 rho U_r^2 S C_L``; from
  flow to body axes ``X = -F_drag cos(alpha) + F_lift sin(alpha)``,
  ``Z = -F_drag sin(alpha) - F_lift cos(alpha)`` (MSS ``forceLiftDrag.m``
  28-40).
* ``alpha = atan2(w_r, u_r)``, ``U_r = |nu_r[0:3]|`` (MSS ``remus100.m``
  126-127).

Lift and drag act in the x-z plane. Deviations from MSS: MSS hard-codes
``rho = 1026`` (``forceLiftDrag.m`` 26) and ``e = 0.3`` (``coeffLiftDrag.m``
54); here both are parameters. Inputs outside their domain (a non-positive
span, area, efficiency or density, a negative ``C_D0``) are refused by the
declared ranges, naming the input.

``smooth_speed`` (a selector of the graph, default ``False`` = MSS exactly).
``atan2(w_r, u_r)`` has 0/0 derivatives at ``u_r = w_r = 0`` (at rest, and in
a pure sideways flow), so the CasADi Jacobian of ``tau`` has non-finite
entries there. ``smooth_speed=True`` adds the declared parameter
``smooth_speed_epsilon`` = ``eps`` (m/s, ``> 0``; a numerical regularisation
chosen by the composition author, not a property of the vehicle) and changes
two things (a construction of this module, a departure from ``remus100.m``
126 and ``forceLiftDrag.m`` 34-40):

* the angle that sets the coefficients is the half-angle form

      alpha = 2 atan( w_r / (sqrt(u_r^2 + w_r^2 + eps^2) + u_r) ),

  which equals ``atan2(w_r, u_r)`` for ``eps = 0``, is smooth for every
  ``(u_r, w_r)`` (the denominator is positive) and is 0 at ``u_r = w_r = 0``;
* the force is resolved on the flow direction itself,
  ``(u_r, w_r) / sqrt(u_r^2 + w_r^2 + eps^2)``, in place of ``(cos alpha,
  sin alpha)``. Then ``X u_r + Z w_r = -q C_D (u_r^2 + w_r^2) / sqrt(u_r^2 +
  w_r^2 + eps^2) <= 0`` for every state: drag opposes the flow, lift does no
  work, whatever the smooth angle does (dissipative, Fossen 2011,
  Property 6.3, p. 123). (Resolving on ``cos alpha, sin alpha`` of the
  smooth angle, an earlier form, pushed a craft going straight astern further
  astern: at ``u_r < 0, w_r = 0`` the smooth angle is 0, not pi.)

The force is built on ``U_r^2``, which needs no root, so it stays exactly
zero at rest. ``smooth_speed=False`` keeps the MSS lines unchanged
(``atan2``, ``cos``, ``sin``); the smooth form tends to them as ``eps -> 0``.

Band where ``eps > 0`` departs from MSS. With ``K >= 10`` the force equals
the MSS form to within about ``0.6 / K`` of its magnitude outside

    B(eps, K) = { u_r^2 + w_r^2 < (K eps)^2 }  (near rest, any v_r)
              U { u_r < 0 and |u_r w_r| < K eps^2 }  (the astern cut),

(measured 2026-10-07 on 200,000 log-spaced states for ``eps`` = 1e-6 and
1e-3 m/s, ``K`` = 10, 100, 1000, with span 0.2 m, area 0.36 m^2, ``C_D0``
0.1, ``e`` 0.3; the factor 0.6 grows with the lift slope; a test holds
``K = 1000`` on the spheroid parameter set; a construction of this module,
no published source). The
sources: the direction differs by at most ``eps^2 / (2 (u_r^2 + w_r^2))``
relative; the angle by order ``eps^2 / (u_r^2 + w_r^2)`` ahead and
``eps^2 / |u_r w_r|`` astern. Inside the astern strip the angle passes from
``-pi`` to ``pi`` through 0 instead of jumping (MSS's own force jumps there:
its lift flips sign), so lift and induced drag fade and the force is the
parasitic drag opposing the motion. A smooth angle cannot avoid such a band:
``atan2`` winds once around the origin. The linear-lift model is a
forward-speed model in either form.

References
----------
[Fossen 2011] Fossen, T. I. (2011). *Handbook of Marine Craft Hydrodynamics
    and Motion Control*, 1st ed. John Wiley & Sons, Chichester. Ch. 6,
    Property 6.3, p. 123.
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2
    with the fixes of 2026-10-07. https://github.com/cybergalactic/MSS, MIT
    licence, revision ``cc07579``:
    ``LIBRARY/modeling/forceLiftDrag.m`` 26-40, ``coeffLiftDrag.m`` 54-68;
    ``CRAFT/AUV/models/remus100.m`` 126-127.

Author:    Enio Krizman
Date:      2026-10-05
"""

import math

import casadi as ca
from more_transformations.more_casadi_transformations import Parameter, symbols

_DECLARED = (
    Parameter("span", (1, 1), "m", "span b (for a body of revolution: its diameter)", 0.0,
              minimum_exclusive=True),
    Parameter("planform_area", (1, 1), "m^2", "planform area S", 0.0, minimum_exclusive=True),
    Parameter("parasitic_drag_coefficient", (1, 1), "1", "parasitic drag coefficient C_D0", 0.0),
    Parameter("oswald_efficiency", (1, 1), "1", "Oswald efficiency e", 0.0, minimum_exclusive=True),
    Parameter("water_density", (1, 1), "kg/m^3", "water density rho", 0.0, minimum_exclusive=True),
)

_SMOOTH_SPEED_EPSILON = Parameter(
    "smooth_speed_epsilon", (1, 1), "m/s",
    "regularisation eps of the smooth angle of attack and flow direction (smooth_speed=True); "
    "numerical, not a vehicle quantity",
    0.0, minimum_exclusive=True,
)


def lift_drag_parameters(*, smooth_speed=False):
    """The declared parameter set: a tuple of ``Parameter`` (name, shape, SI
    unit, meaning, admissible range), in the order of the block's inputs;
    with ``smooth_speed=True`` ``smooth_speed_epsilon`` is the last one."""
    _check_selectors(smooth_speed)
    return _DECLARED + ((_SMOOTH_SPEED_EPSILON,) if smooth_speed else ())


def _angle_of_attack(nu_r, smooth_speed_epsilon):
    """``alpha = atan2(w_r, u_r)`` (MSS ``remus100.m`` 126), or its smooth form
    (module docstring) when ``smooth_speed_epsilon`` is a symbol, not ``None``."""
    u_r, w_r = nu_r[0], nu_r[2]
    if smooth_speed_epsilon is not None:
        offset_squared = w_r**2 + smooth_speed_epsilon**2
        plane_speed = ca.sqrt(u_r**2 + offset_squared)  # sqrt(u_r^2 + w_r^2 + eps^2): construction of this module
        # plane_speed + u_r, written without cancellation for u_r < 0 (identity (a + u)(a - u) = a^2 - u^2)
        denominator = ca.if_else(
            u_r >= 0, plane_speed + u_r, offset_squared / (plane_speed - u_r)
        )
        return 2 * ca.atan(w_r / denominator)  # half-angle form, construction of this module (module docstring)
    return ca.atan2(w_r, u_r)  # (remus100.m 126)


def _check_selectors(smooth_speed):
    if not isinstance(smooth_speed, bool):
        raise ValueError("smooth_speed must be True or False")


def lift_drag_casadi(*, smooth_speed=False):
    """The block: ``(nu_r, span, planform_area, parasitic_drag_coefficient,
    oswald_efficiency, water_density[, smooth_speed_epsilon]) -> tau = [X 0 Z
    0 0 0]`` (MSS ``forceLiftDrag.m`` 28-40).

    Contract
    --------
    ``smooth_speed`` (bool, default ``False`` = MSS) is a selector: it chooses
    the graph (module docstring); anything but a bool raises ``ValueError``.
    Inputs: ``nu_r`` (6x1, relative to the water, body axes, SI) and the
    parameters of ``lift_drag_parameters(smooth_speed=smooth_speed)`` by
    name. Output: ``tau`` (6x1, the force on the vehicle). No value is
    checked inside the graph (``check_values`` with the declaration does
    that, ``smooth_speed_epsilon > 0`` included).
    """
    declared = lift_drag_parameters(smooth_speed=smooth_speed)
    p = symbols(declared)
    epsilon = p.get("smooth_speed_epsilon")
    nu_r = ca.SX.sym("nu_r", 6)
    aspect_ratio = p["span"] * p["span"] / p["planform_area"]  # AR = b^2 / S (coeffLiftDrag.m 57)
    lift_slope = math.pi * aspect_ratio / (1.0 + ca.sqrt(1.0 + (0.5 * aspect_ratio) ** 2))  # (coeffLiftDrag.m 60)
    induced_drag_factor = 1.0 / (math.pi * p["oswald_efficiency"] * aspect_ratio)  # 1 / (pi e AR) (coeffLiftDrag.m 68)
    pressure_area = 0.5 * p["water_density"] * p["planform_area"]  # 1/2 rho S (forceLiftDrag.m 30-31)

    alpha = _angle_of_attack(nu_r, epsilon)
    speed_squared = nu_r[0] ** 2 + nu_r[1] ** 2 + nu_r[2] ** 2  # U_r^2 (remus100.m 127)
    lift_coefficient = lift_slope * alpha  # (coeffLiftDrag.m 61, sigma = 0)
    drag_coefficient = p["parasitic_drag_coefficient"] + lift_coefficient**2 * induced_drag_factor  # (coeffLiftDrag.m 68)
    pressure = pressure_area * speed_squared  # 1/2 rho U_r^2 S (forceLiftDrag.m 30-31)
    if epsilon is not None:
        # flow direction in the x-z plane, construction of this module (module docstring): dissipative
        plane_speed = ca.sqrt(nu_r[0] ** 2 + nu_r[2] ** 2 + epsilon**2)
        sin_alpha, cos_alpha = nu_r[2] / plane_speed, nu_r[0] / plane_speed  # direction cosines of (u_r, w_r): construction
    else:
        sin_alpha, cos_alpha = ca.sin(alpha), ca.cos(alpha)
    surge = -pressure * (drag_coefficient * cos_alpha - lift_coefficient * sin_alpha)  # (forceLiftDrag.m 35)
    heave = -pressure * (drag_coefficient * sin_alpha + lift_coefficient * cos_alpha)  # (forceLiftDrag.m 37)
    tau = ca.vertcat(surge, 0.0, heave, 0.0, 0.0, 0.0)
    names = [d.name for d in declared]
    return ca.Function("lift_drag", [nu_r, *[p[name] for name in names]], [tau], ["nu_r", *names], ["tau"])
