"""Section-plane speed without a flow-angle correction.

References
----------
[MSS] Fossen, T. I. (2026). Marine Systems Simulator, MIT, CRAFT/AUV/models/remus100.m 234-235 @ cc07579.

Author:    Enio Krizman
Date:      2026-10-08
"""
import casadi as ca


def none_parameters():
    return ()


def none_casadi():
    v, chord, lift = ca.SX.sym("fin_velocity", 3), ca.SX.sym("chord_axis", 3), ca.SX.sym("lift_axis", 3)
    speed_squared = ca.dot(v, chord)**2 + ca.dot(v, lift)**2  # (MSS remus100.m 234-235, in fin section axes)
    return ca.Function("no_flow_angle", [v, chord, lift], [ca.SX(0), speed_squared],  # beta = 0 (MSS remus100.m 234-245)
                       ["fin_velocity", "chord_axis", "lift_axis"], ["flow_angle", "speed_squared"])
