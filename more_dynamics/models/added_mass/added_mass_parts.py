"""The ``added_mass`` slot: the fluid's inertia, apart from the rigid body
so a Capytaine or identified matrix can be
plugged into any hull without a new rigid-body form.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics
    and Motion Control. Wiley. Eq. 6.53, p. 121 (``added_mass_matrix``, cited
    in ``more_dynamics.models.added_mass.added_mass``).
[MSS] Fossen, T. I. MSS, MIT, @ cc07579: LIBRARY/modeling/imlay61.m 31-59
    (Lamb's spheroid derivatives, through the rigid-body block's ``spheroid``
    form); CRAFT/USV/models/otter.m 152-159 (scaled derivatives).

Author:    Enio Krizman
Date:      2026-10-08
"""

import casadi as ca

from more_transformations.more_casadi_transformations import Parameter

from more_dynamics.models.added_mass.added_mass import added_mass_matrix, scaled_added_mass_derivatives
from more_dynamics.models.rigid_body.kinetics import rigid_body_casadi
from more_dynamics.models.rigid_body.mass_properties import rigid_body_parameters

from more_dynamics.models.shared.wiring import function_from, restrict

SCALED_DERIVATIVES_PARAMETERS = (
    Parameter("added_mass_coefficients", (6, 1), "1",
             "c of M_A = -diag(c [A11, m_hull, m_hull, I11, I22, I33]) (otter.m 152-159)"),
)


def lamb_spheroid():
    """``M_A`` of the ideal homogeneous spheroid (imlay61.m 31-59), read off
    the rigid-body block's own ``spheroid`` form."""
    return restrict(rigid_body_casadi(mass_properties="spheroid"), ["M_A"], "lamb_spheroid")


def lamb_spheroid_parameters():
    """The part's own declared parameter: ``roll_added_inertia_ratio``.
    ``semi_major_axis``/``semi_minor_axis`` (hull_form) and ``water_density``
    (site) remain as inputs but are vehicle couplings."""
    by_name = {d.name: d for d in rigid_body_parameters("spheroid")}
    return (by_name["roll_added_inertia_ratio"],)


def scaled_derivatives():
    """``M_A = -diag(c [A11, m_hull, m_hull, I11, I22, I33])`` from the hull's
    own couplings (otter.m 152-159)."""
    hull_mass = ca.SX.sym("hull_mass")
    length = ca.SX.sym("length")
    water_density = ca.SX.sym("water_density")
    inertia = ca.SX.sym("inertia", 3, 3)
    coefficients = ca.SX.sym("added_mass_coefficients", 6)
    m_a = added_mass_matrix(scaled_added_mass_derivatives(hull_mass, length, water_density, inertia, coefficients))
    return function_from("scaled_derivatives",
                         {"hull_mass": hull_mass, "length": length, "water_density": water_density,
                          "inertia": inertia, "added_mass_coefficients": coefficients},
                         {"M_A": m_a})
