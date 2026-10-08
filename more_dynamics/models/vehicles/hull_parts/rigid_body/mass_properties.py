"""Mass properties of the rigid-body + added-mass block (CasADi): the declared
parameters of each form and the symbolic ``m``, ``r_g``, ``I`` (about the CG)
and added-mass derivative set built from them.

A form is chosen by the keyword ``mass_properties`` of the builders in
``kinetics.py``, never by a vehicle name. ``rigid_body_parameters(form)``
is the declaration of that form (name, shape, SI unit, meaning, admissible
range of every parameter); ``mass_properties_casadi(form, p)`` takes one
symbol or number per declared name and returns the derived quantities by
name.

Equations per form (keys in References):

* ``"hull_with_payload"``: ``m = m_h + m_p``, ``r_g = (m_h r_h + m_p r_p) /
  m``, ``I = I_h - m_h S(r_h - r_g)^2 - m_p S(r_p - r_g)^2`` about the
  combined CG, ``I_h = m_h diag((R_s [B, L, L])^2)`` (parallel axes, Fossen
  2011, Theorem 3.1, eq. 3.34, p. 50; MSS ``otter.m`` 122-128).
* ``"displacement_hull"``: ``nabla = Cb L B T``, ``m = rho nabla``,
  ``I = m diag((R_s [B, L, L])^2)`` about the CG, wetted surface
  ``L B + 2 T B`` (definitions of this module, see below; the radii-of-
  gyration form as MSS ``rbody.m`` 35). Rotational added mass scales
  ``I`` about the CG (MSS ``otter.m`` 155-157).
* ``"spheroid"``: ``m = rho 4/3 pi a b^2``, ``I = diag(2/5 m b^2,
  1/5 m (a^2 + b^2), 1/5 m (a^2 + b^2))`` (MSS ``spheroid.m`` 35-42,
  ``imlay61.m`` 31-34); Lamb added mass of the displaced fluid
  (``added_mass.py``). Two densities: MSS uses 1025 for the body
  (``spheroid.m`` 35) and 1026 for the displaced fluid (``imlay61.m`` 31);
  the two parameters may be given the same value.
* ``"explicit"``: given ``m``, the principal moments of inertia about the CG
  (the inertia is ``diag`` of them), ``r_g = r_cg - origin`` (the CO -> CG
  vector of Fossen 2011, eq. 3.24, p. 49), the derivative set
  ``[X_du, Y_dv, Z_dw, K_dp, M_dq, N_dr]`` (Fossen 2011, eq. 6.53, p. 121).

Two choices of this module, stated against the references:

* ``"hull_with_payload"``: the inertia is taken about the combined CG, as
  MSS ``otter.m`` 122-128 does since MSS commit ``880b2ef`` (2026-04-20);
  shifting it to the CO first and then applying ``H`` would move it twice.
* ``"displacement_hull"``: the inertia stays about the CG and is moved to
  the CO once, by ``H`` (as ``rbody.m`` 35-45); the Newton-Euler wrench
  (Fossen 2011, eqs. 3.33 and 3.40, pp. 50-51) agrees with that. The hull
  mass ``rho Cb L B T`` and the wetted surface ``L B + 2 T B`` (a box of
  the hull's length, beam and draft, bottom and two sides) are definitions
  of this module; no source equation.
* ``"explicit"``: the inertia is declared by its three principal moments, so
  a non-diagonal inertia cannot be given.

References
----------
[Fossen 2011] Fossen, T. I. (2011). *Handbook of Marine Craft Hydrodynamics
    and Motion Control*, 1st ed. John Wiley & Sons, Chichester. Ch. 3,
    eq. 3.24, p. 49; eqs. 3.33-3.40, pp. 50-51; Ch. 6, eq. 6.53, p. 121.
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2.
    https://github.com/cybergalactic/MSS, MIT licence, revision ``72656d1``:
    ``CRAFT/USV/models/otter.m`` 122-128, 152-157;
    ``LIBRARY/modeling/rbody.m`` 35, ``spheroid.m`` 35-42, ``imlay61.m``
    31-34.

Author:    Enio Krizman
Date:      2026-10-05
"""

import math

import casadi as ca
from more_transformations.more_casadi_transformations import MatrixTransforms, Parameter

from more_dynamics.models.vehicles.hull_parts.added_mass.added_mass import (
    lamb_k_factors,
    scaled_added_mass_derivatives,
    spheroid_added_mass_derivatives,
)

MASS_PROPERTIES_FORMS = ("hull_with_payload", "displacement_hull", "spheroid", "explicit")


def _positive(name, unit, meaning):
    return Parameter(name, (1, 1), unit, meaning, 0.0, minimum_exclusive=True)


_LENGTH = _positive("length", "m", "hull length L")
_BEAM = _positive("beam", "m", "hull beam B")
_WATER_DENSITY = _positive("water_density", "kg/m^3", "water density rho")
_RADII_OF_GYRATION = Parameter(
    "radii_of_gyration", (3, 1), "1",
    "radii of gyration about the CG divided by [B, L, L] (roll, pitch, yaw)",
)
_ADDED_MASS_COEFFICIENTS = Parameter(
    "added_mass_coefficients", (6, 1), "1",
    "c of [X_du, Y_dv, Z_dw, K_dp, M_dq, N_dr] = c * [A11, m, m, I11, I22, I33]; negative c adds inertia",
)
_HULL_CENTER_OF_GRAVITY = Parameter(
    "hull_center_of_gravity", (3, 1), "m", "CO -> hull CG r_h, body axes (FRD)"
)

_DECLARATIONS = {
    "hull_with_payload": (
        _LENGTH,
        _BEAM,
        _WATER_DENSITY,
        _positive("hull_mass", "kg", "hull mass m_h, payload excluded"),
        Parameter("payload_mass", (1, 1), "kg", "point payload mass m_p", 0.0),
        _HULL_CENTER_OF_GRAVITY,
        Parameter("payload_position", (3, 1), "m", "CO -> payload r_p, body axes (FRD)"),
        _ADDED_MASS_COEFFICIENTS,
        _RADII_OF_GYRATION,
    ),
    "displacement_hull": (
        _WATER_DENSITY,
        _LENGTH,
        _BEAM,
        _positive("draft", "m", "hull draft T"),
        _positive("block_coefficient", "1", "block coefficient Cb = nabla / (L B T)"),
        _RADII_OF_GYRATION,
        _HULL_CENTER_OF_GRAVITY,
        _ADDED_MASS_COEFFICIENTS,
    ),
    "spheroid": (
        _positive("semi_major_axis", "m", "semi-axis a along x_b, a > b"),
        _positive("semi_minor_axis", "m", "semi-axis b across x_b"),
        _positive("body_density", "kg/m^3", "mean density of the body (m = rho_b 4/3 pi a b^2)"),
        _WATER_DENSITY,
        Parameter("roll_added_inertia_ratio", (1, 1), "1",
                  "r44: added roll inertia / roll inertia of the displaced fluid", 0.0),
        Parameter("body_center_of_gravity", (3, 1), "m", "CO -> CG r_g, body axes (FRD)"),
    ),
    "explicit": (
        _positive("body_mass", "kg", "body mass m"),
        Parameter("principal_inertia", (3, 1), "kg*m^2",
                  "principal moments of inertia [Ix, Iy, Iz] about the CG, axes along the body axes",
                  0.0, minimum_exclusive=True),
        Parameter("body_center_of_gravity", (3, 1), "m",
                  "CG position, axes parallel to the body axes (FRD); CO -> CG when origin is zero"),
        Parameter("origin", (3, 1), "m",
                  "CO position in the axes of body_center_of_gravity; r_g = body_center_of_gravity - origin"),
        Parameter("added_mass_derivatives", (6, 1), "kg (1-3), kg*m^2 (4-6)",
                  "[X_du, Y_dv, Z_dw, K_dp, M_dq, N_dr]; M_A = -diag"),
    ),
}

# Outputs every form gives (after M, C_RB, C_A, M_RB, M_A in the block) and
# the extra ones of one form. No output carries the name of a parameter.
COMMON_QUANTITIES = ("mass", "center_of_gravity", "inertia")
FORM_QUANTITIES = {
    "hull_with_payload": (),
    "displacement_hull": ("displaced_volume", "wetted_surface"),
    "spheroid": ("lamb_k_factors",),
    "explicit": (),
}


def _check_form(mass_properties):
    if mass_properties not in _DECLARATIONS:
        raise ValueError(
            f"mass_properties must be one of {MASS_PROPERTIES_FORMS}, got {mass_properties!r}"
        )


def rigid_body_parameters(mass_properties="hull_with_payload"):
    """The declared parameter set of one form: a tuple of ``Parameter``
    (name, shape, SI unit, meaning, admissible range), in the order of the
    block's inputs. An unknown form raises ``ValueError``."""
    _check_form(mass_properties)
    return _DECLARATIONS[mass_properties]


def _hull_with_payload(p):
    """Hull with a point payload. Radii of gyration are scaled by ``[beam, length,
    length]`` and apply to the hull only; added mass uses the hull mass."""
    hull_mass, payload_mass = p["hull_mass"], p["payload_mass"]
    r_hull, r_payload = p["hull_center_of_gravity"], p["payload_position"]
    mass = hull_mass + payload_mass  # m + mp (otter.m 124, 142)
    radii = p["radii_of_gyration"] * ca.vertcat(p["beam"], p["length"], p["length"])  # R44, R55, R66 (otter.m 95-97)
    inertia_hull_cg = ca.diag(hull_mass * radii**2)  # (otter.m 122)
    r_g = (hull_mass * r_hull + payload_mass * r_payload) / mass  # (otter.m 124)
    s_hull = MatrixTransforms.skew(r_hull - r_g)  # (otter.m 125)
    s_payload = MatrixTransforms.skew(r_payload - r_g)  # (otter.m 126)
    # parallel axes to the combined CG (Fossen 2011, eq. 3.34, p. 50; otter.m 127)
    inertia = inertia_hull_cg - hull_mass * s_hull @ s_hull - payload_mass * s_payload @ s_payload
    derivatives = scaled_added_mass_derivatives(  # (pattern of otter.m 152-157)
        hull_mass, p["length"], p["water_density"], inertia, p["added_mass_coefficients"]
    )
    return {
        "mass": mass,
        "center_of_gravity": r_g,
        "inertia": inertia,
        "added_mass_derivatives": derivatives,
    }


def _displacement_hull(p):
    """Neutrally buoyant hull: ``nabla = Cb L B T``, ``m = rho nabla``,
    ``I = m diag((R_s [B, L, L])^2)`` about the CG, wetted surface
    ``L B + 2 T B`` (definitions of this module; the wetted surface is used
    by the surge damping)."""
    length, beam, draft = p["length"], p["beam"], p["draft"]
    displaced_volume = p["block_coefficient"] * length * beam * draft  # nabla = Cb L B T: definition of the block coefficient; no source equation
    mass = p["water_density"] * displaced_volume  # m = rho nabla: the hull floats at its draft (Archimedes); no source equation
    radii = p["radii_of_gyration"] * ca.vertcat(beam, length, length)  # R44, R55, R66 = R_s [B, L, L] (the scaling of otter.m 95-97)
    inertia = ca.diag(mass * radii**2)  # about the CG (rbody.m 35)
    derivatives = scaled_added_mass_derivatives(  # (pattern of otter.m 152-157)
        mass, length, p["water_density"], inertia, p["added_mass_coefficients"]
    )
    return {
        "mass": mass,
        "center_of_gravity": p["hull_center_of_gravity"],
        "inertia": inertia,
        "added_mass_derivatives": derivatives,
        "displaced_volume": displaced_volume,
        "wetted_surface": length * beam + 2.0 * draft * beam,  # S = L B + 2 T B: bottom and two sides of a box L x B x T; approximation of this module, no source equation
    }


def _spheroid_mass_and_inertia(density, a, b):
    """``m = rho 4/3 pi a b^2``, ``I = diag(2/5 m b^2, 1/5 m (a^2+b^2),
    1/5 m (a^2+b^2))`` (MSS ``spheroid.m`` 35-42, ``imlay61.m`` 31-34)."""
    mass = density * 4.0 / 3.0 * math.pi * a * b**2  # (spheroid.m 36; imlay61.m 32)
    lateral = 0.2 * mass * (a**2 + b**2)  # (spheroid.m 40-41; imlay61.m 34)
    return mass, ca.diag(ca.vertcat(0.4 * mass * b**2, lateral, lateral))  # (spheroid.m 39, 42; imlay61.m 33)


def _spheroid(p):
    """Prolate spheroid ``a > b``: body from ``body_density``, Lamb added
    mass from the displaced fluid (``water_density``); ``M_A(4,4) = r44 I_x``
    of the fluid (MSS ``imlay61.m`` 47)."""
    a, b = p["semi_major_axis"], p["semi_minor_axis"]
    k_factors = lamb_k_factors(a, b)  # (Lamb 1932, Arts. 114-115, pp. 153-155; imlay61.m 50-56)
    mass, inertia = _spheroid_mass_and_inertia(p["body_density"], a, b)  # (spheroid.m 35-42)
    fluid_mass, fluid_inertia = _spheroid_mass_and_inertia(p["water_density"], a, b)  # (imlay61.m 31-34)
    derivatives = spheroid_added_mass_derivatives(  # (imlay61.m 47, 59)
        k_factors, fluid_mass, fluid_inertia, p["roll_added_inertia_ratio"]
    )
    return {
        "mass": mass,
        "center_of_gravity": p["body_center_of_gravity"],
        "inertia": inertia,
        "added_mass_derivatives": derivatives,
        "lamb_k_factors": k_factors,
    }


def _explicit(p):
    """Given ``m``, the principal moments about the CG, the CG, the origin
    and the derivative set ``[X_du, Y_dv, Z_dw, K_dp, M_dq, N_dr]``;
    ``r_g = body_center_of_gravity - origin``."""
    return {
        "mass": p["body_mass"],
        # r_g = r_cg - r_co, the CO -> CG vector of (Fossen 2011, eq. 3.24, p. 49)
        "center_of_gravity": p["body_center_of_gravity"] - p["origin"],
        "inertia": ca.diag(p["principal_inertia"]),  # I_g of (Fossen 2011, eq. 3.21, p. 49) with principal axes along the body axes
        "added_mass_derivatives": p["added_mass_derivatives"],
    }


_FORMS = {
    "hull_with_payload": _hull_with_payload,
    "displacement_hull": _displacement_hull,
    "spheroid": _spheroid,
    "explicit": _explicit,
}


def mass_properties_casadi(mass_properties, parameters):
    """The derived quantities of one form, by name.

    ``parameters`` maps every declared name of the form to a CasADi symbol
    or a number. Returns a dict with ``mass``, ``center_of_gravity`` (CO ->
    CG, 3x1), ``inertia`` (3x3, about the CG), ``added_mass_derivatives``
    (6x1, ``M_A = -diag``; not an output of the block, which gives ``M_A``),
    and the quantities of ``FORM_QUANTITIES`` for that form
    (``displaced_volume``, ``wetted_surface``; ``lamb_k_factors``).

    A parameter and an output of the block never share a name: the inputs
    that locate the CG are ``hull_center_of_gravity`` or
    ``body_center_of_gravity``, the output is ``center_of_gravity``.
    """
    _check_form(mass_properties)
    return _FORMS[mass_properties](parameters)
