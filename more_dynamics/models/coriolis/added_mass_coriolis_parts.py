"""The ``added_mass_coriolis`` slot: ``C_A(nu_r)`` from the added-mass
matrix, physics form and the MSS shortcut that removes the Munk moment, as
two labelled parts instead of a flag inside one part.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics
    and Motion Control. Wiley. Eqs. 6.43, 6.52, p. 120-121, through
    ``more_dynamics.models.rigid_body.kinetics.added_mass_coriolis_casadi``,
    which cites them line by line.
[MSS] Fossen, T. I. MSS, MIT, @ cc07579: LIBRARY/modeling/m2c.m 33-48;
    CRAFT/AUV/models/remus100.m 207-210 (Munk moment removed).

Author:    Enio Krizman
Date:      2026-10-08
"""

import casadi as ca

from more_dynamics.models.rigid_body.kinetics import added_mass_coriolis_casadi

from more_dynamics.models.shared.wiring import function_from


def _part(name, stabilize):
    m_a = ca.SX.sym("M_A", 6, 6)
    nu_r = ca.SX.sym("nu_r", 6)
    return function_from(name, {"M_A": m_a, "nu_r": nu_r},
                         {"C_A": added_mass_coriolis_casadi(m_a, nu_r, stabilize)})


KIRCHHOFF_FULL_PARAMETERS = ()        # M_A, nu_r are both couplings/kinematic; no own parameter
MUNK_COUPLINGS_REMOVED_PARAMETERS = ()


def kirchhoff_full():
    """The physics form: ``C_A = m2c(M_A, nu_r)`` (m2c.m 33-48), every term
    kept."""
    return _part("kirchhoff_full", False)


def munk_couplings_removed():
    """MSS shortcut: the Munk moment of a slender body removed
    (remus100.m 207-210)."""
    return _part("munk_couplings_removed", True)
