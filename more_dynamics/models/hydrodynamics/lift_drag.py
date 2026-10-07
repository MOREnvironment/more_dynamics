"""Lift and drag of a slender body (numpy constants, CasADi algebra).

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
span, area, efficiency or density, a negative ``C_D0``) raise ``ValueError``
naming the input.

``smooth_speed_epsilon`` (default 0 = MSS exactly). ``atan2(w_r, u_r)`` has
0/0 derivatives at ``u_r = w_r = 0`` (at rest, and in a pure sideways flow),
so the CasADi Jacobian of ``tau`` has non-finite entries there. A positive
value ``eps`` (m/s) changes two things (a construction of this module, a
departure from ``remus100.m`` 126 and ``forceLiftDrag.m`` 34-40):

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
zero at rest. ``eps = 0`` keeps the MSS lines unchanged (``atan2``, ``cos``,
``sin``).

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

Ported from the numpy source [MGM] ``dynamics/plant/matrices/hydrodynamics.py``:
``HydroForces.get_lift_drag_coeff`` (62-75, ``sigma = 0``) and
``force_lift_drag`` (78-88).

References
----------
[Fossen 2011] Fossen, T. I. (2011). *Handbook of Marine Craft Hydrodynamics
    and Motion Control*, 1st ed. John Wiley & Sons, Chichester. Ch. 6,
    Property 6.3, p. 123.
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2.
    https://github.com/cybergalactic/MSS, MIT licence, revision ``72656d1``:
    ``LIBRARY/modeling/forceLiftDrag.m`` 26-40, ``coeffLiftDrag.m`` 54-68;
    ``CRAFT/AUV/models/remus100.m`` 126-127.
[MGM] Krizman, E. *more_generic_models*.
    https://github.com/MOREnvironment/more_generic_models (no licence file),
    revision ``524e336``: ``more_generic_models/dynamics/plant/matrices/
    hydrodynamics.py``, lines listed above.
"""

from dataclasses import dataclass

import casadi as ca
import numpy as np


@dataclass(frozen=True)
class LiftDragConstants:
    """Lift slope, induced-drag factor and dynamic-pressure area."""

    lift_slope: float
    parasitic_drag_coefficient: float
    induced_drag_factor: float
    pressure_area: float
    smooth_speed_epsilon: float


def preprocess_lift_drag(
    span: float,
    planform_area: float,
    parasitic_drag_coefficient: float,
    oswald_efficiency: float,
    water_density: float,
    smooth_speed_epsilon: float = 0.0,
) -> LiftDragConstants:
    """Coefficients of ``forceLiftDrag(b, S, CD_0, alpha, U_r)``.

    Source names: ``b``, ``S``, ``CD_0``, ``e``, ``rho``.
    ``AR = b^2 / S``; ``C_L_alpha = pi AR / (1 + sqrt(1 + (AR/2)^2))``.
    ``smooth_speed_epsilon`` (m/s, default 0 = MSS exactly) smooths the angle
    of attack and resolves the force on the flow direction (module
    docstring, with the band where it departs from MSS).
    """
    positive = {
        "span": span,
        "planform_area": planform_area,
        "oswald_efficiency": oswald_efficiency,
        "water_density": water_density,
    }
    for name, value in positive.items():
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be a positive finite value")
    if not np.isfinite(parasitic_drag_coefficient) or parasitic_drag_coefficient < 0.0:
        raise ValueError("parasitic_drag_coefficient must be a non-negative finite value")
    if not np.isfinite(smooth_speed_epsilon) or smooth_speed_epsilon < 0.0:
        raise ValueError("smooth_speed_epsilon must be a non-negative finite value")

    aspect_ratio = span * span / planform_area  # (coeffLiftDrag.m 57)
    return LiftDragConstants(
        lift_slope=np.pi * aspect_ratio / (1.0 + np.sqrt(1.0 + (0.5 * aspect_ratio) ** 2)),  # (coeffLiftDrag.m 60)
        parasitic_drag_coefficient=float(parasitic_drag_coefficient),
        induced_drag_factor=1.0 / (np.pi * oswald_efficiency * aspect_ratio),  # (coeffLiftDrag.m 68)
        pressure_area=0.5 * water_density * planform_area,  # (forceLiftDrag.m 30-31)
        smooth_speed_epsilon=float(smooth_speed_epsilon),
    )


def _angle_of_attack(nu_r: ca.SX, smooth_speed_epsilon: float) -> ca.SX:
    """``alpha = atan2(w_r, u_r)`` (MSS ``remus100.m`` 126), or its smooth form
    (module docstring)."""
    u_r, w_r = nu_r[0], nu_r[2]
    if smooth_speed_epsilon > 0.0:
        offset_squared = w_r**2 + smooth_speed_epsilon**2
        plane_speed = ca.sqrt(u_r**2 + offset_squared)
        # plane_speed + u_r, written without cancellation for u_r < 0
        denominator = ca.if_else(
            u_r >= 0, plane_speed + u_r, offset_squared / (plane_speed - u_r)
        )
        return 2 * ca.atan(w_r / denominator)  # half-angle form (module docstring)
    return ca.atan2(w_r, u_r)  # (remus100.m 126)


def lift_drag_casadi(constants: LiftDragConstants) -> ca.Function:
    """``nu_r -> tau = [X 0 Z 0 0 0]`` (MSS ``forceLiftDrag.m`` 28-40)."""
    nu_r = ca.SX.sym("nu_r", 6)
    epsilon = constants.smooth_speed_epsilon
    alpha = _angle_of_attack(nu_r, epsilon)
    speed_squared = nu_r[0] ** 2 + nu_r[1] ** 2 + nu_r[2] ** 2  # U_r^2 (remus100.m 127)
    lift_coefficient = constants.lift_slope * alpha  # (coeffLiftDrag.m 61, sigma = 0)
    drag_coefficient = (  # (coeffLiftDrag.m 68)
        constants.parasitic_drag_coefficient
        + lift_coefficient**2 * constants.induced_drag_factor
    )
    pressure = constants.pressure_area * speed_squared  # 1/2 rho U_r^2 S (forceLiftDrag.m 30-31)
    if epsilon > 0.0:
        # flow direction in the x-z plane (module docstring): dissipative
        plane_speed = ca.sqrt(nu_r[0] ** 2 + nu_r[2] ** 2 + epsilon**2)
        sin_alpha, cos_alpha = nu_r[2] / plane_speed, nu_r[0] / plane_speed
    else:
        sin_alpha, cos_alpha = ca.sin(alpha), ca.cos(alpha)
    surge = -pressure * (drag_coefficient * cos_alpha - lift_coefficient * sin_alpha)  # (forceLiftDrag.m 35)
    heave = -pressure * (drag_coefficient * sin_alpha + lift_coefficient * cos_alpha)  # (forceLiftDrag.m 37)
    tau = ca.vertcat(surge, 0.0, heave, 0.0, 0.0, 0.0)
    return ca.Function("lift_drag", [nu_r], [tau], ["nu_r"], ["tau"])
