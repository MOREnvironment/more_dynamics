"""Quadratic axial drag and quadratic roll damping of a slender body (CasADi).

Frame: BODY FRD about the origin CO (x forward, y starboard, z down).
``nu_r = [u, v, w, p, q, r]`` is the velocity relative to the water, SI units; the
result is a force and moment ``tau`` in N and N m about the CO.

Equations (keys in References):

* Axial drag ``X = -1/2 rho c_d A_f u |u|`` with the frontal area
  ``A_f = pi d^2 / 4`` and the axial drag coefficient ``c_d`` (Prestero 2001,
  eqs. 4.5-4.6, p. 25; the viscous drag always opposes the motion, which is why
  ``u |u|`` and not ``u^2``). ``d`` is the maximum diameter of the body of
  revolution (the ``span`` of the vehicle).
* Roll damping ``K = K_p|p| p |p|`` (Prestero 2001, eqs. 4.10-4.11, p. 27, from
  the cross-flow drag of the fins). Dissipative for ``K_p|p| <= 0``, which the
  declared range enforces. The source calls its own estimate "at best a rough
  approximation" and prefers experimental data, so the coefficient is a value
  the composition gives, with no default.

References
----------
[Prestero 2001] Prestero, T. (2001). Verification of a six-degree of freedom
    simulation model for the REMUS autonomous underwater vehicle. MIT/WHOI MSc
    thesis. Eqs. 4.5-4.6, p. 25; eqs. 4.10-4.11, p. 27 (read as text).

Author:    Enio Krizman
"""

import math

import casadi as ca
from more_transformations.more_casadi_transformations import Parameter

from more_dynamics.models.shared.wiring import function_from

QUADRATIC_AXIAL_DRAG_PARAMETERS = (
    Parameter("hull_frontal_drag_coefficient", (1, 1), "1", "axial drag coefficient c_d of the hull on its frontal area",
              0.0, minimum_exclusive=True),
)
QUADRATIC_ROLL_DAMPING_PARAMETERS = (
    Parameter("roll_quadratic_damping", (1, 1), "N*m*s^2/rad^2",
              "coefficient K_p|p| of the quadratic roll damping K = K_p|p| p |p|; not positive (dissipative)",
              None, 0.0),
)


def quadratic_axial_drag():
    """``(nu_r, span, water_density, hull_frontal_drag_coefficient) -> tau = [X 0 0 0 0 0]`` with
    ``X = -1/2 rho c_d (pi d^2/4) u |u|`` (Prestero 2001, eqs. 4.5-4.6, p. 25)."""
    nu_r = ca.SX.sym("nu_r", 6)
    span, density = ca.SX.sym("span"), ca.SX.sym("water_density")
    drag = ca.SX.sym("hull_frontal_drag_coefficient")
    frontal_area = math.pi * span**2 / 4.0  # A_f of the circle of diameter d
    surge = -0.5 * density * drag * frontal_area * nu_r[0] * ca.fabs(nu_r[0])  # (Prestero 2001, eq. 4.5, p. 25)
    tau = ca.vertcat(surge, 0.0, 0.0, 0.0, 0.0, 0.0)
    return function_from("quadratic_axial_drag", {"nu_r": nu_r, "span": span, "water_density": density,
                                                  "hull_frontal_drag_coefficient": drag}, {"tau": tau})


def quadratic_roll_damping():
    """``(nu_r, roll_quadratic_damping) -> tau = [0 0 0 K 0 0]`` with ``K = K_p|p| p |p|`` (Prestero 2001, eq. 4.10,
    p. 27)."""
    nu_r = ca.SX.sym("nu_r", 6)
    coefficient = ca.SX.sym("roll_quadratic_damping")
    roll = coefficient * nu_r[3] * ca.fabs(nu_r[3])  # (Prestero 2001, eq. 4.10, p. 27)
    tau = ca.vertcat(0.0, 0.0, 0.0, roll, 0.0, 0.0)
    return function_from("quadratic_roll_damping", {"nu_r": nu_r, "roll_quadratic_damping": coefficient}, {"tau": tau})
