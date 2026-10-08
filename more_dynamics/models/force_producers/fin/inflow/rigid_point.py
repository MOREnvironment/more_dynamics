"""Fin inflow at a rigidly attached point, retaining all rotational terms.

References
----------
[Prestero 2001] Prestero, T. (2001). Verification of a six-degree of freedom simulation model for the REMUS autonomous underwater vehicle. MIT/WHOI MSc thesis, eq. 4.40, p. 32.

Author:    Enio Krizman
Date:      2026-10-08
"""
import casadi as ca
from more_transformations.more_casadi_transformations import MatrixTransforms


def rigid_point_parameters():
    return ()


def rigid_point_casadi():
    nu_r, r = ca.SX.sym("nu_r", 6), ca.SX.sym("fin_position", 3)
    velocity = nu_r[:3] + MatrixTransforms.skew(nu_r[3:]) @ r  # (Prestero 2001, eq. 4.40, p. 32; y,z terms retained)
    return ca.Function("rigid_point", [nu_r, r], [velocity],
                       ["nu_r", "fin_position"], ["fin_velocity"])
