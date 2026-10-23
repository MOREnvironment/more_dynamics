"""Steering of a unit that does not steer: no command, no state, deflection zero.

References
----------
A fixed unit has no steering equation; the form exists so every unit has the same
couplings (a convention of this library).

Author:    Enio Krizman
Date:      2026-10-23
"""
import casadi as ca

from more_dynamics.models.propulsor.declarations import part_function


def fixed_parameters():
    return ()


def fixed_casadi():
    return part_function("fixed", {"command": ca.SX.sym("command", 0), "servo_state": ca.SX.sym("servo_state", 0)},
                         {"deflection": ca.SX(0.0), "servo_state_dot": ca.SX.zeros(0, 1)}, {}, ())
