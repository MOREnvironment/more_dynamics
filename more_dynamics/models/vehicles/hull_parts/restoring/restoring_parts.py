"""The ``restoring`` slot: ``g`` (subtracted in the equation of motion), and
whatever quantities its form derives (weight, buoyancy, the restoring
matrix, displaced volume, draft) as couplings for the hydrodynamic-load
parts that need them.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics
    and Motion Control. Wiley. Through ``more_dynamics.models.vehicles.hull_parts.restoring``,
    which cites it line by line (ch. 4, pp. 59-80).
[MSS] Fossen, T. I. MSS, MIT, @ cc07579: CRAFT/AUV/models/remus100.m 214
    (neutral buoyancy, W = m g = B); CRAFT/USV/models/otter.m 121-122, 172-193
    (equilibrium draft of the twin pontoons).

Author:    Enio Krizman
Date:      2026-10-08
"""

import casadi as ca

from more_transformations.more_casadi_transformations import Parameter

from more_dynamics.models.vehicles.hull_parts.restoring.submerged import submerged_hydrostatics_casadi
from more_dynamics.models.vehicles.hull_parts.restoring.surface import surface_hydrostatics_casadi

from more_dynamics.models.shared.wiring import function_from

SUBMERGED_NEUTRAL_PARAMETERS = (
    Parameter("center_of_buoyancy", (3, 1), "m", "CO -> CB r_bb, body axes (FRD); above the CG"),
)
TWIN_HULL_METACENTRIC_EQUILIBRIUM_DRAFT_PARAMETERS = (
    Parameter("longitudinal_inertia_factor", (1, 1), "1", "c_L: I_L = n c_L (1/12) B L^3",
             0.0, minimum_exclusive=True),
    Parameter("longitudinal_center_of_flotation", (1, 1), "m", "LCF x_F from the CO, body x"),
    Parameter("reference_point", (3, 1), "m", "CO -> point P where G is given, body axes (FRD)"),
)


def submerged_neutral():
    """``W = m g``, ``B = W`` (remus100.m 214), fed into the submerged
    restoring block."""
    block = submerged_hydrostatics_casadi()
    eta = ca.SX.sym("eta", 6)
    mass = ca.SX.sym("mass")
    gravity = ca.SX.sym("gravity")
    center_of_gravity = ca.SX.sym("center_of_gravity", 3)
    center_of_buoyancy = ca.SX.sym("center_of_buoyancy", 3)
    weight = mass * gravity  # W = m g (remus100.m 214)
    buoyancy = weight  # B = W, neutral by assertion (remus100.m 214)
    g = block(eta=eta, weight=weight, buoyancy=buoyancy,
              center_of_gravity=center_of_gravity, center_of_buoyancy=center_of_buoyancy)["g"]
    return function_from(
        "submerged_neutral",
        {"eta": eta, "mass": mass, "gravity": gravity,
         "center_of_gravity": center_of_gravity, "center_of_buoyancy": center_of_buoyancy},
        {"g": g, "weight": weight, "buoyancy": buoyancy, "center_of_buoyancy": center_of_buoyancy},
    )


def twin_hull_metacentric_equilibrium_draft():
    """``nabla = m / rho``, ``T = nabla / (2 Cb L B_pont)`` (otter.m 121-122)
    fed into the surface restoring block (two hulls, gravity given by the
    site)."""
    block = surface_hydrostatics_casadi(hull_count=2, gravity_source="value")
    names = ("mass", "water_density", "gravity", "length", "pontoon_beam", "pontoon_block_coefficient",
             "pontoon_waterplane_coefficient", "pontoon_lateral_offset", "longitudinal_inertia_factor",
             "longitudinal_center_of_flotation")
    s = {n: ca.SX.sym(n) for n in names}
    s["center_of_gravity"] = ca.SX.sym("center_of_gravity", 3)
    s["reference_point"] = ca.SX.sym("reference_point", 3)
    s["eta"] = ca.SX.sym("eta", 6)
    displaced_volume = s["mass"] / s["water_density"]  # nabla = m / rho (otter.m 121)
    draft = displaced_volume / (2 * s["pontoon_block_coefficient"] * s["pontoon_beam"] * s["length"])  # (otter.m 122)
    out = block(eta=s["eta"], length=s["length"], beam=s["pontoon_beam"], draft=draft,
                displacement_volume=displaced_volume,
                waterplane_area_coefficient=s["pontoon_waterplane_coefficient"],
                center_of_gravity=s["center_of_gravity"], hull_lateral_offset=s["pontoon_lateral_offset"],
                longitudinal_inertia_factor=s["longitudinal_inertia_factor"],
                longitudinal_center_of_flotation=s["longitudinal_center_of_flotation"],
                reference_point=s["reference_point"], water_density=s["water_density"], gravity=s["gravity"])
    return function_from("twin_hull_metacentric_equilibrium_draft", s, {
        "g": out["g"], "G": out["G"], "G_CF": out["G_CF"],
        "displaced_volume": displaced_volume, "draft": draft})
