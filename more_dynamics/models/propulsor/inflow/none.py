"""Inflow of a unit that sees no flow: advance speed zero, thrust and torque factors one.

References
----------
No equation: the unit takes nothing from the water (a convention of this library).

Author:    Enio Krizman
Date:      2026-10-23
"""
import casadi as ca

from more_dynamics.models.propulsor.declarations import part_function


def none_parameters():
    return ()


def none_casadi():
    inputs = {"nu_r": ca.SX.sym("nu_r", 6), "thrust_axis": ca.SX.sym("thrust_axis", 3),
              "thrust_position": ca.SX.sym("thrust_position", 3)}
    return part_function("none", inputs, {"advance_speed": ca.SX(0.0), "thrust_factor": ca.SX(1.0),
                                          "torque_factor": ca.SX(1.0)}, {}, ())
