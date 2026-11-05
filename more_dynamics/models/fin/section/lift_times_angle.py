"""Linear lift and the drag ``C_D = C_L alpha`` (``= a alpha^2``), the drag law ``lift_times_angle``.

``C_L = a alpha`` and ``C_D = C_L alpha = a alpha^2``, ``a`` the lift slope, ``alpha`` the effective angle of the fin: the
drag of the angle-free rudder and stern-plane form ``X = -1/2 rho U^2 A C_L delta^2``. A shortcut with no zero-lift and no
induced drag.

References
----------
[MSS] Fossen, T. I. (2026). Marine Systems Simulator, MIT, CRAFT/AUV/models/remus100.m 239-246 @ 6e78e98.

Author:    Enio Krizman
Date:      2026-10-08
"""
import casadi as ca
from more_transformations.more_casadi_transformations import Parameter, symbols


def lift_times_angle_parameters():
    return (Parameter("lift_slope", (1, 1), "1/rad", "lift coefficient per angle", 0.0, minimum_exclusive=True),)


def lift_times_angle_casadi():
    p = symbols(lift_times_angle_parameters())
    alpha = ca.SX.sym("angle_of_attack")
    lift = p["lift_slope"] * alpha  # (MSS remus100.m 243, 246)
    drag = p["lift_slope"] * alpha**2  # (MSS remus100.m 239-240; comparison shortcut)
    return ca.Function("lift_times_angle", [alpha, p["lift_slope"]], [lift, drag],
                       ["angle_of_attack", "lift_slope"], ["lift_coefficient", "drag_coefficient"])
