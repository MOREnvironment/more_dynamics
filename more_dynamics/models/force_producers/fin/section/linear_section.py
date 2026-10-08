"""Linear fin lift with constant zero-lift drag.

References
----------
[Prestero 2001] Prestero, T. (2001). Verification of a six-degree of freedom simulation model for the REMUS autonomous underwater vehicle. MIT/WHOI MSc thesis, eq. 4.37, p. 31.

Author:    Enio Krizman
Date:      2026-10-08
"""
import casadi as ca
from more_transformations.more_casadi_transformations import Parameter, symbols


def linear_section_parameters():
    return (Parameter("lift_slope", (1, 1), "1/rad", "lift coefficient per angle", 0.0, minimum_exclusive=True),
            Parameter("zero_lift_drag", (1, 1), "1", "measured drag coefficient at zero lift", 0.0))


def linear_section_casadi():
    p = symbols(linear_section_parameters())
    alpha = ca.SX.sym("angle_of_attack")
    lift = p["lift_slope"] * alpha  # (Prestero 2001, eq. 4.37, p. 31)
    drag = p["zero_lift_drag"]  # measured option; Prestero 2001, eq. 4.37, p. 31 omits fin drag
    return ca.Function("linear_section", [alpha, *p.values()], [lift, drag],
                       ["angle_of_attack", *p], ["lift_coefficient", "drag_coefficient"])
