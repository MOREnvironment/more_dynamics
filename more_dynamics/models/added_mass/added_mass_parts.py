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


GIVEN_DERIVATIVE_ADDED_MASS_NAMES = ("Xudot", "Yvdot", "Ypdot", "Yrdot", "Zwdot", "Zqdot", "Kvdot", "Kpdot", "Krdot",
                                     "Mwdot", "Mqdot", "Nvdot", "Npdot", "Nrdot")


def given_derivative_added_mass_parameters():
    """The 14 published nondimensional acceleration derivatives, each its own
    declared parameter (rule 16: primitives, not a re-derivation)."""
    return tuple(Parameter(f"nondim_{n}", (1, 1), "1", f"nondimensional {n} (npsauv.m 166-171)")
                for n in GIVEN_DERIVATIVE_ADDED_MASS_NAMES)


def given_derivative_added_mass():
    """``M_A`` of a body whose added mass is given as published nondimensional
    acceleration derivatives (not computed from a hull form): the sparse
    symmetric matrix of ``npsauv.m`` 166-172, de-scaled by Fossen's
    prime-scaling (Appendix D.2) with ``r3 = 1/2 rho L^3`` and
    ``Tinv = diag(1,1,1,L,L,L)``. Couplings ``length``, ``water_density``."""
    names = [f"nondim_{n}" for n in GIVEN_DERIVATIVE_ADDED_MASS_NAMES]
    s = {n: ca.SX.sym(n) for n in names}
    d = {n: s[f"nondim_{n}"] for n in GIVEN_DERIVATIVE_ADDED_MASS_NAMES}
    length = ca.SX.sym("length")
    water_density = ca.SX.sym("water_density")
    zero = ca.SX(0)
    nondim_matrix = -ca.vertcat(  # (npsauv.m 166-171)
        ca.horzcat(d["Xudot"], zero, zero, zero, zero, zero),
        ca.horzcat(zero, d["Yvdot"], zero, d["Ypdot"], zero, d["Yrdot"]),
        ca.horzcat(zero, zero, d["Zwdot"], zero, d["Zqdot"], zero),
        ca.horzcat(zero, d["Kvdot"], zero, d["Kpdot"], zero, d["Krdot"]),
        ca.horzcat(zero, zero, d["Mwdot"], zero, d["Mqdot"], zero),
        ca.horzcat(zero, d["Nvdot"], zero, d["Npdot"], zero, d["Nrdot"]))
    r3 = 0.5 * water_density * length ** 3  # (Fossen's prime-scaling, Appendix D.2; npsauv.m 119)
    t_inv = ca.diag(ca.vertcat(1.0, 1.0, 1.0, length, length, length))  # (npsauv.m 117)
    m_a = r3 * t_inv @ nondim_matrix @ t_inv  # (npsauv.m 172)
    return function_from("given_derivative_added_mass", {**s, "length": length, "water_density": water_density},
                         {"M_A": m_a})


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
