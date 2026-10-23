"""Housing of a unit that adds no load of its own: zero wrench.

References
----------
No equation (a convention of this library).

Author:    Enio Krizman
Date:      2026-10-23
"""
import casadi as ca

from more_dynamics.models.propulsor.declarations import part_function


def none_parameters():
    return ()


def none_casadi():
    inputs = {"nu_r": ca.SX.sym("nu_r", 6), "thrust_axis": ca.SX.sym("thrust_axis", 3),
              "thrust_position": ca.SX.sym("thrust_position", 3), "steering_angle": ca.SX.sym("steering_angle"),
              "water_density": ca.SX.sym("water_density")}
    return part_function("none", inputs, {"tau": ca.SX.zeros(6)}, {}, ())
