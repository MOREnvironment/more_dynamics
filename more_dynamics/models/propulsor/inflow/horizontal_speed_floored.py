"""Inflow of a unit whose power limit is set by the horizontal speed through the water: the advance speed is the
horizontal speed ``sqrt(u_r^2 + v_r^2)`` bounded below by ``speed_floor`` (the power limit ``P / U`` stays finite at
rest), thrust factor 1, torque factor 1. The speed has a zero (not NaN) derivative at standstill.

References
----------
[Krizman] Krizman, E. Model file ``grethe.m``, MSS fork branch ``grethe-vehicle`` @ 679a6e8,
    CRAFT/USV/models/grethe.m 391 (``U_h = max(sqrt(nu_r(1)^2 + nu_r(2)^2), 0.1)``; read, the file may still be
    revised).

Author:    Enio Krizman
Date:      2026-10-26
"""
import casadi as ca
from more_transformations.more_casadi_transformations import symbols

from more_dynamics.models.propulsor.declarations import part_function, positive
from more_dynamics.models.shared.force_producer_common import safe_norm


def horizontal_speed_floored_parameters():
    return (positive("speed_floor", "m/s", "lowest horizontal speed the advance speed takes (numerical floor)"),)


def horizontal_speed_floored_casadi():
    declared = horizontal_speed_floored_parameters()
    p = symbols(declared)
    inputs = {"nu_r": ca.SX.sym("nu_r", 6), "thrust_axis": ca.SX.sym("thrust_axis", 3),
              "thrust_position": ca.SX.sym("thrust_position", 3)}
    advance = ca.fmax(safe_norm(inputs["nu_r"][0:2]), p["speed_floor"])  # U_h = max(|nu_r(1:2)|, floor) (grethe.m 391)
    return part_function("horizontal_speed_floored", inputs,
                         {"advance_speed": advance, "thrust_factor": ca.SX(1.0), "torque_factor": ca.SX(1.0)},
                         p, declared)
