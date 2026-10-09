"""Linear lift and quadratic deflection drag, the MSS comparison form.

References
----------
[MSS] Fossen, T. I. (2026). Marine Systems Simulator, MIT, CRAFT/AUV/models/remus100.m 238-245 @ cc07579.

Author:    Enio Krizman
Date:      2026-10-08
"""
import casadi as ca
from more_transformations.more_casadi_transformations import Parameter, symbols


def quadratic_drag_parameters():
    return (Parameter("lift_slope", (1, 1), "1/rad", "lift coefficient per angle", 0.0, minimum_exclusive=True),)


def quadratic_drag_casadi():
    p = symbols(quadratic_drag_parameters())
    alpha = ca.SX.sym("angle_of_attack")
    lift = p["lift_slope"] * alpha  # (MSS remus100.m 242, 245)
    drag = p["lift_slope"] * alpha**2  # (MSS remus100.m 238-239; comparison shortcut)
    return ca.Function("quadratic_drag", [alpha, p["lift_slope"]], [lift, drag],
                       ["angle_of_attack", "lift_slope"], ["lift_coefficient", "drag_coefficient"])
