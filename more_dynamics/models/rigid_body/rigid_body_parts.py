"""The ``rigid_body`` slot, cut from the rigid-body block (``M_RB``, ``C_RB``,
``mass``, ``center_of_gravity``, ``inertia``; added mass is its own slot).

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics
    and Motion Control. Wiley. Through ``more_dynamics.models.rigid_body``,
    which cites it line by line.
[MSS] Fossen, T. I. MSS, MIT, @ cc07579: LIBRARY/modeling/spheroid.m 35-42
    (homogeneous spheroid); CRAFT/USV/models/otter.m 94-98, 122-128 (hull
    with a point payload).

Author:    Enio Krizman
Date:      2026-10-08
"""

import math

import casadi as ca

from more_transformations.more_casadi_transformations import Parameter

from more_dynamics.models.rigid_body.kinetics import rigid_body_casadi
from more_dynamics.models.rigid_body.mass_properties import rigid_body_parameters

from more_dynamics.models.shared.wiring import function_from, restrict, with_passthrough

RIGID_BODY_OUTPUTS = ["M_RB", "C_RB", "mass", "center_of_gravity", "inertia"]

# Declared (own) names, by hand against each part's actual ``name_in()``
# after ``restrict``: ``semi_major_axis``/``semi_minor_axis`` (hull_form) and
# ``water_density`` (site) also remain as inputs but are vehicle couplings,
# resolved by name at the composite plugin, never frozen here.
_HOMOGENEOUS_SPHEROID_OWN = ("body_density", "body_center_of_gravity")
_HULL_WITH_POINT_PAYLOAD_OWN = ("hull_mass", "payload_mass", "hull_center_of_gravity",
                               "payload_position", "radii_of_gyration")


SPHEROID_BODY_MASS_PARAMETERS = (
    Parameter("body_mass", (1, 1), "kg", "body mass m of the vehicle", 0.0, minimum_exclusive=True),
)


def spheroid_body_density():
    """``(body_mass, semi_major_axis, semi_minor_axis) -> body_density``: the mean density of the homogeneous
    spheroid that has the given mass, ``rho_b = m / (4/3 pi a b^2)`` (the volume of a prolate spheroid, MSS
    ``spheroid.m`` 36, ``m = rho_b 4/3 pi a b^2`` solved for ``rho_b``; the mass of ``remus100.m`` 3 is a given)."""
    body_mass = ca.SX.sym("body_mass")
    semi_major_axis = ca.SX.sym("semi_major_axis")
    semi_minor_axis = ca.SX.sym("semi_minor_axis")
    volume = 4.0 / 3.0 * math.pi * semi_major_axis * semi_minor_axis**2  # (spheroid.m 36)
    return function_from("spheroid_body_density", {
        "body_mass": body_mass, "semi_major_axis": semi_major_axis, "semi_minor_axis": semi_minor_axis},
        {"body_density": body_mass / volume})


def homogeneous_spheroid():
    """Mass and inertia from a homogeneous prolate spheroid (spheroid.m
    35-42)."""
    return restrict(rigid_body_casadi(mass_properties="spheroid", coriolis="co"), RIGID_BODY_OUTPUTS,
                    "homogeneous_spheroid")


def homogeneous_spheroid_parameters(*, density_given=True):
    """The part's own declared parameters (``_HOMOGENEOUS_SPHEROID_OWN``); its
    remaining inputs are vehicle couplings, not plugin parameters. With
    ``density_given=False`` ``body_density`` is not declared: the vehicle
    computes it (``spheroid_body_density``) and the part reads it by name."""
    by_name = {d.name: d for d in rigid_body_parameters("spheroid")}
    names = _HOMOGENEOUS_SPHEROID_OWN if density_given else tuple(n for n in _HOMOGENEOUS_SPHEROID_OWN if n != "body_density")
    return tuple(by_name[n] for n in names)


def hull_with_point_payload():
    """Hull mass plus a point payload (otter.m 94-98, 122-128); ``hull_mass``
    is handed on as a coupling for the added-mass slot's scaled
    derivatives."""
    cut = restrict(rigid_body_casadi(mass_properties="hull_with_payload", coriolis="co"), RIGID_BODY_OUTPUTS,
                   "hull_with_point_payload_cut")
    return with_passthrough(cut, "hull_with_point_payload", ["hull_mass"])


def hull_with_point_payload_parameters():
    """The part's own declared parameters (``_HULL_WITH_POINT_PAYLOAD_OWN``);
    ``length``/``beam`` remain as inputs but are ``hull_form`` couplings."""
    by_name = {d.name: d for d in rigid_body_parameters("hull_with_payload")}
    return tuple(by_name[n] for n in _HULL_WITH_POINT_PAYLOAD_OWN)
