"""Small fin flow angle from section-plane velocity.

The source assumes nonzero chordwise flow. At exactly zero chordwise speed,
the finite guard sets the angle to zero; this guard is outside that model.

References
----------
[Prestero 2001] Prestero, T. (2001). Verification of a six-degree of freedom simulation model for the REMUS autonomous underwater vehicle. MIT/WHOI MSc thesis, eqs. 4.41-4.43, pp. 32-33.

Author:    Enio Krizman
Date:      2026-10-08
"""
import casadi as ca


def small_angle_parameters():
    return ()


def small_angle_casadi():
    v, chord, lift = ca.SX.sym("fin_velocity", 3), ca.SX.sym("chord_axis", 3), ca.SX.sym("lift_axis", 3)
    u_c, u_n = ca.dot(v, chord), ca.dot(v, lift)
    flow_angle = ca.if_else(u_c == 0, 0, u_n / ca.if_else(u_c == 0, 1, u_c))  # (Prestero 2001, eq. 4.42, p. 32; zero guard)
    speed_squared = u_c**2  # (Prestero 2001, eq. 4.43, p. 33)
    return ca.Function("small_angle", [v, chord, lift], [flow_angle, speed_squared],
                       ["fin_velocity", "chord_axis", "lift_axis"], ["flow_angle", "speed_squared"])
