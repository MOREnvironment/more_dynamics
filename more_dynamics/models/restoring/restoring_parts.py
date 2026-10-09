"""Restoring of a vehicle as two functions of the vehicle's own quantities: a
submerged body, whose buoyancy is set by one of three methods, and a surface
craft of one or two hulls at its equilibrium draft.

Submerged: ``W = m g`` (Fossen 2011, eq. 4.1, p. 59) and ``B`` by
``buoyancy_method``: ``"neutral"`` (``B = W``, Fossen 2011, eq. 4.7, p. 61;
MSS ``remus100.m`` 214), ``"from_volume"`` (``B = rho g nabla`` with the
displaced volume given, Fossen 2011, eq. 4.1, p. 59) or ``"given"`` (``B`` as
a value, for a measured reserve of buoyancy). Fidelity: neutral < from_volume
(needs the volume of the hull) < given (needs a weighing of the reserve
buoyancy at the quay); every method keeps the centre of buoyancy and the
centre of gravity apart (``r_bb``, ``r_bg``).

Surface: ``nabla = m / rho``, ``T = nabla / (n C_b L B_hull)`` (MSS
``otter.m`` 121-122 for two hulls) fed into the surface restoring block
(``gravity`` given by the site). The wetted surface of the hulls at that draft
is an output too: the Mumford approximation ``S = n 1.025 L (C_b B_hull +
1.7 T)`` of a displacement hull (MSS ``XuuITTC.m`` 38, which attributes it to
Mumford; the original not read), summed over the ``n`` hulls without
interference between them, or given as a value (``wetted_surface_given``).
Raise: the wetted surface from the hull lines.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics
    and Motion Control. Wiley. Eqs. 4.1, 4.5, 4.7, pp. 59-61 (through
    ``more_dynamics.models.restoring.submerged``, which cites them line by
    line).
[MSS] Fossen, T. I. MSS, MIT, @ cc07579: CRAFT/AUV/models/remus100.m 214
    (neutral buoyancy, W = m g = B); CRAFT/USV/models/otter.m 121-122,
    172-193 (equilibrium draft of the twin pontoons);
    LIBRARY/modeling/XuuITTC.m 38 (Mumford wetted-area approximation).

Author:    Enio Krizman
Date:      2026-10-08
"""

import casadi as ca

from more_transformations.more_casadi_transformations import Parameter

from more_dynamics.models.restoring.submerged import submerged_hydrostatics_casadi
from more_dynamics.models.restoring.surface import surface_hydrostatics_casadi

from more_dynamics.models.shared.wiring import function_from

BUOYANCY_METHODS = ("neutral", "from_volume", "given")


def submerged_restoring_parameters(buoyancy_method):
    """The parameters of the chosen method: the centre of buoyancy always,
    the displaced volume for ``"from_volume"``, the buoyancy for ``"given"``."""
    if buoyancy_method not in BUOYANCY_METHODS:
        raise ValueError(f"buoyancy_method must be one of {BUOYANCY_METHODS}, got {buoyancy_method!r}")
    declared = [Parameter("center_of_buoyancy", (3, 1), "m", "CO -> CB r_bb, body axes (FRD)")]
    if buoyancy_method == "from_volume":
        declared.append(Parameter("displaced_volume", (1, 1), "m^3", "displaced volume nabla of the submerged hull",
                                  0.0, minimum_exclusive=True))
    if buoyancy_method == "given":
        declared.append(Parameter("buoyancy", (1, 1), "N", "buoyancy B", 0.0, minimum_exclusive=True))
    return tuple(declared)


def submerged_restoring(buoyancy_method):
    """``(eta, mass, gravity, center_of_gravity, center_of_buoyancy[, water_density, displaced_volume | buoyancy]) ->
    (g, weight, buoyancy, center_of_buoyancy)``; ``W = m g`` (Fossen 2011, eq. 4.1, p. 59), ``B`` by the method."""
    submerged_restoring_parameters(buoyancy_method)
    block = submerged_hydrostatics_casadi()
    eta = ca.SX.sym("eta", 6)
    mass = ca.SX.sym("mass")
    gravity = ca.SX.sym("gravity")
    center_of_gravity = ca.SX.sym("center_of_gravity", 3)
    center_of_buoyancy = ca.SX.sym("center_of_buoyancy", 3)
    inputs = {"eta": eta, "mass": mass, "gravity": gravity, "center_of_gravity": center_of_gravity}
    weight = mass * gravity  # W = m g (Fossen 2011, eq. 4.1, p. 59; remus100.m 214)
    if buoyancy_method == "neutral":
        buoyancy = weight  # B = W (Fossen 2011, eq. 4.7, p. 61; remus100.m 214)
    elif buoyancy_method == "from_volume":
        water_density = ca.SX.sym("water_density")
        displaced_volume = ca.SX.sym("displaced_volume")
        inputs.update({"water_density": water_density, "displaced_volume": displaced_volume})
        buoyancy = water_density * gravity * displaced_volume  # B = rho g nabla (Fossen 2011, eq. 4.1, p. 59)
    else:
        buoyancy = ca.SX.sym("buoyancy")
        inputs["buoyancy"] = buoyancy
    inputs["center_of_buoyancy"] = center_of_buoyancy
    g = block(eta=eta, weight=weight, buoyancy=buoyancy,
              center_of_gravity=center_of_gravity, center_of_buoyancy=center_of_buoyancy)["g"]
    return function_from("submerged_restoring", inputs,
                         {"g": g, "weight": weight, "buoyancy": buoyancy, "center_of_buoyancy": center_of_buoyancy})


SURFACE_RESTORING_PARAMETERS = (
    Parameter("longitudinal_inertia_factor", (1, 1), "1", "c_L: I_L = n c_L (1/12) B L^3",
              0.0, minimum_exclusive=True),
    Parameter("longitudinal_center_of_flotation", (1, 1), "m", "LCF x_F from the CO, body x"),
    Parameter("reference_point", (3, 1), "m", "CO -> point P where G is given, body axes (FRD)"),
)


WETTED_SURFACE_PARAMETER = Parameter("wetted_surface", (1, 1), "m^2", "wetted surface S of the hulls",
                                     0.0, minimum_exclusive=True)
MUMFORD_AREA_FACTOR = 1.025  # (XuuITTC.m 38)
MUMFORD_DRAFT_FACTOR = 1.7  # (XuuITTC.m 38)


def surface_restoring(hull_count, *, wetted_surface_given=False):
    """``nabla = m / rho``, ``T = nabla / (n C_b L B_hull)`` (otter.m 121-122) fed into the surface restoring block
    (``hull_count`` 1 or 2, gravity given by the site). The hull geometry (``length``, ``hull_beam``,
    ``hull_block_coefficient``, ``hull_waterplane_coefficient``, ``hull_lateral_offset`` for two hulls) and the
    mass, density, gravity and centre of gravity are the vehicle's quantities. The output ``wetted_surface`` is
    the Mumford approximation of the ``hull_count`` hulls at the equilibrium draft, or the given input of that name
    when ``wetted_surface_given``."""
    block = surface_hydrostatics_casadi(hull_count=hull_count, gravity_source="value")
    names = ("mass", "water_density", "gravity", "length", "hull_beam", "hull_block_coefficient",
             "hull_waterplane_coefficient")
    s = {n: ca.SX.sym(n) for n in names}
    if hull_count == 2:
        s["hull_lateral_offset"] = ca.SX.sym("hull_lateral_offset")
    for p in SURFACE_RESTORING_PARAMETERS:
        s[p.name] = ca.SX.sym(p.name, *p.shape)
    s["center_of_gravity"] = ca.SX.sym("center_of_gravity", 3)
    s["eta"] = ca.SX.sym("eta", 6)
    displaced_volume = s["mass"] / s["water_density"]  # nabla = m / rho (otter.m 121)
    draft = displaced_volume / (hull_count * s["hull_block_coefficient"] * s["hull_beam"] * s["length"])  # (otter.m 122)
    args = dict(eta=s["eta"], length=s["length"], beam=s["hull_beam"], draft=draft,
                displacement_volume=displaced_volume, waterplane_area_coefficient=s["hull_waterplane_coefficient"],
                center_of_gravity=s["center_of_gravity"],
                longitudinal_inertia_factor=s["longitudinal_inertia_factor"],
                longitudinal_center_of_flotation=s["longitudinal_center_of_flotation"],
                reference_point=s["reference_point"], water_density=s["water_density"], gravity=s["gravity"])
    if hull_count == 2:
        args["hull_lateral_offset"] = s["hull_lateral_offset"]
    out = block(**args)
    if wetted_surface_given:
        s["wetted_surface"] = ca.SX.sym("wetted_surface")
        wetted_surface = s["wetted_surface"]
    else:  # S = n 1.025 L (C_b B + 1.7 T) per hull (XuuITTC.m 38, Mumford)
        wetted_surface = hull_count * MUMFORD_AREA_FACTOR * s["length"] * (
            s["hull_block_coefficient"] * s["hull_beam"] + MUMFORD_DRAFT_FACTOR * draft)
    return function_from(f"surface_restoring_{hull_count}_hull", s, {
        "g": out["g"], "G": out["G"], "G_CF": out["G_CF"],
        "displaced_volume": displaced_volume, "draft": draft, "wetted_surface": wetted_surface})
