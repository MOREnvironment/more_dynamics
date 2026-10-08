"""The ``rigid_body`` slot, cut from the rigid-body block (``M_RB``, ``C_RB``,
``mass``, ``center_of_gravity``, ``inertia``; added mass is its own slot).

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics
    and Motion Control. Wiley. Through ``more_dynamics.models.vehicles.hull_parts.rigid_body``,
    which cites it line by line.
[MSS] Fossen, T. I. MSS, MIT, @ cc07579: LIBRARY/modeling/spheroid.m 35-42
    (homogeneous spheroid); CRAFT/USV/models/otter.m 94-98, 122-128 (hull
    with a point payload).

Author:    Enio Krizman
Date:      2026-10-08
"""

from more_dynamics.models.vehicles.hull_parts.rigid_body.kinetics import rigid_body_casadi
from more_dynamics.models.vehicles.hull_parts.rigid_body.mass_properties import rigid_body_parameters

from more_dynamics.models.wiring import restrict, with_passthrough

RIGID_BODY_OUTPUTS = ["M_RB", "C_RB", "mass", "center_of_gravity", "inertia"]

# Declared (own) names, by hand against each part's actual ``name_in()``
# after ``restrict``: ``semi_major_axis``/``semi_minor_axis`` (hull_form) and
# ``water_density`` (site) also remain as inputs but are vehicle couplings,
# resolved by name at the composite plugin, never frozen here.
_HOMOGENEOUS_SPHEROID_OWN = ("body_density", "body_center_of_gravity")
_HULL_WITH_POINT_PAYLOAD_OWN = ("hull_mass", "payload_mass", "hull_center_of_gravity",
                               "payload_position", "radii_of_gyration")


def homogeneous_spheroid():
    """Mass and inertia from a homogeneous prolate spheroid (spheroid.m
    35-42)."""
    return restrict(rigid_body_casadi(mass_properties="spheroid", coriolis="co"), RIGID_BODY_OUTPUTS,
                    "homogeneous_spheroid")


def homogeneous_spheroid_parameters():
    """The part's own declared parameters (``_HOMOGENEOUS_SPHEROID_OWN``); its
    remaining inputs are vehicle couplings, not plugin parameters."""
    by_name = {d.name: d for d in rigid_body_parameters("spheroid")}
    return tuple(by_name[n] for n in _HOMOGENEOUS_SPHEROID_OWN)


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
