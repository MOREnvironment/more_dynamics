"""Unit fin-interference factors.

References
----------
[Prestero 2001] Prestero, T. (2001). Verification of a six-degree of freedom simulation model for the REMUS autonomous underwater vehicle. MIT/WHOI MSc thesis, eqs. 4.41-4.43, pp. 32-33.
[MSS] Fossen, T. I. (2026). Marine Systems Simulator, MIT, CRAFT/AUV/models/remus100.m 238-245 @ cc07579.

Author:    Enio Krizman
Date:      2026-10-08
"""
import casadi as ca


def none_parameters():
    return ()


def none_casadi():
    # Unit factors leave deflection and flow angle unchanged (Prestero 2001, eq. 4.41, p. 32; MSS remus100.m 238-245).
    return ca.Function("no_interference", [], [ca.DM(1), ca.DM(1)], [],
                       ["deflection_factor", "flow_angle_factor"])
