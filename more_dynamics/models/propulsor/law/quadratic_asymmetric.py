"""Quadratic thrust law with separate coefficients ahead and astern: ``T = k omega |omega|`` with
``omega = 2 pi n`` the shaft speed in rad/s, ``k_pos`` for ``omega > 0`` and ``k_neg`` otherwise;
no torque.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics and Motion Control.
    Wiley. Ch. 12, eq. 12.265, p. 411 (quadratic thrust, separate forward and reverse coefficients).
[MSS] Fossen, T. I. (2026). Marine Systems Simulator, MIT, CRAFT/USV/models/otter_legacy.m 238-244
    @ 6e78e98.

Author:    Enio Krizman
Date:      2026-10-23
"""
import math

import casadi as ca
from more_transformations.more_casadi_transformations import symbols

from more_dynamics.models.propulsor.declarations import part_function, positive

LAW_INPUT_UNIT = "rev/s"


def quadratic_asymmetric_parameters():
    return (positive("positive_thrust_coefficient", "N*s^2/rad^2", "k_pos: T = k_pos omega |omega| for omega > 0"),
            positive("negative_thrust_coefficient", "N*s^2/rad^2", "k_neg: T = k_neg omega |omega| for omega <= 0"))


def quadratic_asymmetric_casadi():
    declared = quadratic_asymmetric_parameters()
    p = symbols(declared)
    n, advance, density = ca.SX.sym("drive_output"), ca.SX.sym("advance_speed"), ca.SX.sym("water_density")
    omega = 2.0 * math.pi * n  # rev/s -> rad/s
    quadratic = omega * ca.fabs(omega)  # omega |omega| (Fossen 2011, eq. 12.265, p. 411)
    thrust = ca.if_else(omega > 0.0, p["positive_thrust_coefficient"] * quadratic,
                        p["negative_thrust_coefficient"] * quadratic)  # (otter_legacy.m 238-244)
    return part_function("quadratic_asymmetric", {"drive_output": n, "advance_speed": advance, "water_density": density},
                         {"thrust": thrust, "torque": ca.SX(0.0)}, p, declared)
