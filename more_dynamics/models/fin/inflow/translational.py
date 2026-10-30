"""Fin inflow from vehicle translation alone.

References
----------
[MSS] Fossen, T. I. (2026). Marine Systems Simulator, MIT, CRAFT/AUV/models/remus100.m 235-236 @ 6e78e98.

Author:    Enio Krizman
Date:      2026-10-08
"""
import casadi as ca


def translational_parameters():
    return ()


def translational_casadi():
    nu_r, r = ca.SX.sym("nu_r", 6), ca.SX.sym("fin_position", 3)
    velocity = nu_r[:3]  # (MSS remus100.m 235-236)
    return ca.Function("translational", [nu_r, r], [velocity],
                       ["nu_r", "fin_position"], ["fin_velocity"])
