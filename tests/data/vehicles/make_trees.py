"""Writes the vehicle compositions as ``.rppws`` trees, in the layout of
Luka's ``.rppws/parts/more_dynamics__hull_vessel/11111111-.../`` tree:
``parts/<plugin id>/<uuid>/{description.json, params/parameters.py[,
callbacks.py], SOURCE.md}`` with children under ``subcomponents/<uuid>/``,
and one script description per configuration
(``script_descriptions/<name>.json``). rpp's own ``ComponentContextBuilder``
builds them. Every number is read from the frozen parameter files in this
folder (``remus100/remus100_parameters_consistent.json``,
``otter/otter_mss_cc07579.json``); none is typed here.

The gate trees under ``tests/data/vehicles/.rppws`` and the user vehicles
under ``scripts/vehicles/.rppws`` are the output of

    python make_trees.py

(run it with the environment that has rpp installed: it takes the plugin
ids from rpp's naming). The tests import this module to write perturbed or
refused variants of a tree into a temporary folder.

Author:    Enio Krizman
Date:      2026-10-08
"""

import json
import math
import re
import shutil
import uuid
from pathlib import Path

HERE = Path(__file__).resolve().parent
LIBRARY = HERE.parents[2]
GATE_TREES = HERE / ".rppws"
USER_TREES = LIBRARY / "scripts" / "vehicles" / ".rppws"
NAMESPACE = uuid.UUID("6d0e1a00-0000-4000-8000-000000000001")

REMUS_FILE = json.loads((HERE / "remus100" / "remus100_parameters_consistent.json").read_text())
REMUS = {name: entry["value"] for name, entry in REMUS_FILE["parameters"].items()}
OTTER = json.loads((HERE / "otter" / "otter_mss_cc07579.json").read_text())
_RIGID_TEST = (HERE.parents[1] / "vehicles" / "hull_parts" / "rigid_body" / "test_rigid_body_block.py").read_text()
CURRENT_SPEED = float(re.search(r"^CURRENT_SPEED = ([0-9.]+)$", _RIGID_TEST, re.M).group(1))
CURRENT_DIRECTION = math.radians(float(
    re.search(r"^CURRENT_DIRECTION = np.deg2rad\(([0-9.]+)\)$", _RIGID_TEST, re.M).group(1)))

# the one kinematic viscosity of both vehicles (cylinderDrag.m 78-80, nu_water = 1e-6 m^2/s)
KINEMATIC_VISCOSITY = 1e-6

TYPES = {  # plugin -> plugin type
    "MarineCraft6DOF": "VehicleModel3D", "SiteAtLatitude": "SiteModel", "SiteGiven": "SiteModel",
    "NoCurrent": "CurrentModel", "HorizontalCurrentYawRateTerms": "CurrentModel",
    "HorizontalCurrentFullRotationRate": "CurrentModel", "ProlateSpheroidMainDimensions": "HullForm",
    "TwinPontoons": "HullForm", "HomogeneousSpheroid": "RigidBodyModel", "HullWithPointPayload": "RigidBodyModel",
    "LambSpheroid": "AddedMassModel", "ScaledDerivatives": "AddedMassModel",
    "KirchhoffFull": "AddedMassCoriolisModel", "MunkCouplingsRemoved": "AddedMassCoriolisModel",
    "SubmergedNeutral": "HydrostaticsModel", "TwinHullMetacentricEquilibriumDraft": "HydrostaticsModel",
    "LinearSurfaceHydrostatics": "HydrostaticsModel",
    "TimeConstantDampingSubmerged": "HydrodynamicsModel", "TimeConstantDampingSurface": "HydrodynamicsModel",
    "HullLiftDrag": "HydrodynamicsModel", "CrossFlowStrip": "HydrodynamicsModel",
    "SurgeResistanceIttc": "HydrodynamicsModel", "LinearSurfaceHydrodynamics": "HydrodynamicsModel",
    "CircularCylinderReynolds": "SectionDragModel", "RectangularSectionHoerner": "SectionDragModel",
    "PrescribedWrench": "ForceProducer", "FinPairsDeflectionOnly": "ForceProducer",
    "PropellerLinearizedOpenWater": "ForceProducer", "JetNozzle": "ForceProducer",
    "LiftingFin": "ForceProducer", "ForceProducerSet": "ForceProducer", "Servo": "ActuatorServo",
    "FinInflowRigidPoint": "FinInflow", "FinInflowTranslational": "FinInflow",
    "FinFlowAngleSmallAngle": "FinFlowAngle", "FinFlowAngleNone": "FinFlowAngle",
    "FinInterferenceNone": "FinInterference", "FinSectionQuadraticDrag": "FinSection",
    "FinSectionLinearSection": "FinSection", "FinInterferenceSlenderBody": "FinInterference",
    "FinSectionLiftingLine": "FinSection",
}
SPEC = {"MarineCraft6DOF": {
    "site": "more_dynamics::SiteModel", "current": "more_dynamics::CurrentModel",
    "hull_form": "more_dynamics::HullForm", "rigid_body": "more_dynamics::RigidBodyModel",
    "added_mass": "more_dynamics::AddedMassModel", "added_mass_coriolis": "more_dynamics::AddedMassCoriolisModel",
    "restoring": "more_dynamics::HydrostaticsModel",
    "hydrodynamic_loads": "List[more_dynamics::HydrodynamicsModel]",
    "force_producers": "List[more_dynamics::ForceProducer]"},
    "CrossFlowStrip": {"section": "more_dynamics::SectionDragModel"},
    "LiftingFin": {"servo": "more_dynamics::ActuatorServo", "inflow": "more_dynamics::FinInflow",
                   "flow_angle": "more_dynamics::FinFlowAngle", "interference": "more_dynamics::FinInterference",
                   "section": "more_dynamics::FinSection"},
    "ForceProducerSet": {"producers": "List[more_dynamics::ForceProducer]"}}
LIST_SLOTS = {"hydrodynamic_loads", "force_producers", "producers"}

# What each part is: fidelity level, assumptions, and the measurement or identification that raises it to the
# next level (one line per part, written into its SOURCE.md).
LABELS = {
    "MarineCraft6DOF": "vehicle | relative-velocity equation of motion with the current's acceleration "
        "(remus100.m 257-259, otter.m 261-263) | assumes a uniform current and rigid-body kinematics | "
        "raise: identification from logs of the parts' own parameters",
    "SiteAtLatitude": "physics form | WGS-84 normal gravity at the latitude (gravity.m 11-12) | assumes sea level "
        "| raise: a gravimeter value at the site (SiteGiven)",
    "SiteGiven": "given value | gravity, density, viscosity as stated by the vehicle's source (otter.m 90-91) | "
        "assumes constant water properties | raise: a CTD cast for density and viscosity on the day",
    "NoCurrent": "physics form | still water | assumes no current | raise: a measured current profile",
    "HorizontalCurrentYawRateTerms": "MSS form | uniform horizontal current, body acceleration with the yaw-rate "
        "terms only (remus100.m 118-122) | assumes the current set is constant in NED | raise: the form with "
        "the full attitude and rotation rate once derived from a read source",
    "HorizontalCurrentFullRotationRate": "MSS form | uniform horizontal current, body acceleration from the full "
        "body-rate skew matrix (otter.m 113-118) | assumes the current set is constant in NED | raise: a "
        "measured current profile",
    "ProlateSpheroidMainDimensions": "physics form | beam, draft, span, section beam = diameter, semi-axes from "
        "length and diameter (remus100.m 131-132, 220-221) | assumes a body of revolution | raise: the "
        "measured outline of the hull",
    "TwinPontoons": "physics form | twin-pontoon primitives (otter.m 92-93, 104-107) | assumes identical "
        "pontoons | raise: the measured hull lines",
    "HomogeneousSpheroid": "physics form | mass and inertia of a homogeneous prolate spheroid (spheroid.m 35-42) "
        "| assumes uniform density | raise: weighing and a swing test for the inertia",
    "HullWithPointPayload": "physics form | hull mass plus a point payload (otter.m 94-98, 122-128) | assumes "
        "given radii of gyration | raise: a swing test for the inertia",
    "LambSpheroid": "physics form | added mass of a prolate spheroid from the Lamb factors (imlay61.m 31-59) | "
        "assumes an ideal fluid | raise: a Capytaine or identified matrix (the part slot takes any)",
    "ScaledDerivatives": "MSS form | added mass from scaled derivatives (otter.m 152-159) | assumes the printed "
        "coefficients | raise: a Capytaine matrix or identification from logs",
    "KirchhoffFull": "physics form | every Coriolis term of the added mass (m2c.m 33-48) | assumes a constant "
        "added-mass matrix | raise: a frequency-dependent matrix",
    "MunkCouplingsRemoved": "MSS shortcut | the Munk couplings removed (remus100.m 207-210) | comparison form, "
        "not the default | raise: KirchhoffFull",
    "SubmergedNeutral": "MSS form | weight equals buoyancy (remus100.m 214) | assumes a neutrally buoyant body | "
        "raise: a measured weight and buoyancy from a quay-side trim test",
    "TwinHullMetacentricEquilibriumDraft": "physics form | draft from the displaced volume at equilibrium and "
        "the metacentric restoring (otter.m 121-122, 177, 192) | assumes box-like pontoons | raise: the "
        "measured draft",
    "TimeConstantDampingSubmerged": "MSS calibration | linear damping from time constants and damping ratios "
        "(Dmtrx.m; remus100.m 217-218) | assumes the printed time constants | raise: identification from "
        "decay tests",
    "TimeConstantDampingSurface": "MSS calibration | linear damping from top speed and thrust (otter.m 195-240) "
        "| assumes the top-speed calibration | raise: identification from logs",
    "HullLiftDrag": "physics form | slender-body lift and induced drag (forceLiftDrag.m) | assumes the printed "
        "planform and Oswald factor | raise: towing-tank or log identification",
    "CrossFlowStrip": "physics form | cross-flow drag by the strip integral (crossFlowDrag.m 54-69) with the "
        "section law of its child | assumes strip theory | raise: a section law from a tank test",
    "SurgeResistanceIttc": "physics form | ITTC-1957 friction line with a form factor and a blend to linear "
        "damping (Fossen 2011, eqs. 6.82-6.85, p. 125) | assumes the printed form factor | raise: a measured "
        "resistance curve",
    "CircularCylinderReynolds": "physics form | drag coefficient of a circular cylinder against Reynolds number "
        "with the aspect-ratio correction (cylinderDrag.m 78-110) | assumes 1e-6 m^2/s viscosity | raise: a "
        "tank-measured section drag",
    "RectangularSectionHoerner": "published table | drag of a rectangular section against beam over twice the "
        "draft (Hoerner.m 47-51) | assumes Hoerner's table | raise: a tank-measured section drag",
    "PrescribedWrench": "boundary | the commanded wrench acts on the vehicle unchanged (hydroVessel.m) | no "
        "actuator model | raise: a fin and propeller model",
    "FinPairsDeflectionOnly": "MSS form | lift of each fin pair proportional to its deflection (remus100.m "
        "228-245) | assumes the printed lift slopes | raise: fins with the local flow (effective angle of attack)",
    "Servo": "one servo with switches | first-order lag, rate limit and angle limit on the command, all at the "
        "plugin's defaults (Murray-Smith 2016 Fig. 1 p. 246; Sarhadi 2026 Fig. 4 p. 4: 0.1 s, 20 deg, 30 deg/s) | "
        "assumes the study's values, a demanding choice, not a measurement | raise: identify the time constant "
        "and rate limit of the real servo from a logged step in the fin command",
    "LiftingFin": "composite | one fin from five parts, force along fixed body axes (Prestero 2001 eqs. 4.37, "
        "4.41-4.43 pp. 31-33; Fossen 2011 eq. 12.226 p. 400) | small-angle form, no fin-body interference | "
        "raise: the local flow and the interference factors once their sources are read",
    "ForceProducerSet": "composite | the children's wrenches summed behind a command map (Fossen 2011 eq. 12.226 "
        "p. 400) | one density for every child | raise: a measured allocation of the real fins",
    "FinInflowTranslational": "MSS form | fin velocity = vehicle translation (remus100.m 234-235) | ignores the "
        "rotation at the fin | raise: FinInflowRigidPoint",
    "FinInflowRigidPoint": "physics form | fin velocity at a point fixed to the body (Prestero 2001 eq. 4.40 p. 32) "
        "| assumes a rigid hull | raise: the flow measured at the fin",
    "FinFlowAngleNone": "MSS form | no flow angle, section speed from the two axis components (remus100.m "
        "234-245) | assumes the flow along the chord | raise: FinFlowAngleSmallAngle",
    "FinFlowAngleSmallAngle": "published form | flow angle u_n / u_c (Prestero 2001 eqs. 4.42-4.43 pp. 32-33) | "
        "small angles | raise: a large-angle flow-axis model",
    "FinInterferenceNone": "MSS form | both interference factors one (Prestero 2001 eq. 4.41 p. 32) | assumes "
        "no fin-body interference | raise: the slender-body factors once Pitts 1957 is read",
    "FinSectionQuadraticDrag": "MSS comparison form | linear lift, drag proportional to the squared angle "
        "(remus100.m 238-245) | assumes the printed lift slope | raise: a measured section polar",
    "FinSectionLinearSection": "published form | linear lift, constant zero-lift drag (Prestero 2001 eq. 4.37 p. 31) "
        "| assumes no stall | raise: a measured section polar or the lifting-line form",
    "PropellerLinearizedOpenWater": "MSS form | linearised open-water characteristic with given coefficients "
        "(remus100.m 148-176) | assumes the printed coefficients | raise: a measured open-water curve or the "
        "Wageningen series",
}


def plugin_id(name):
    from rpp_plugin_registrator.plugin_descriptors.core import plugin_id_from_name
    return plugin_id_from_name(f"more_dynamics::{name}")


def node(plugin, params=None, source=None, children=None, linked=None):
    """A component: plugin name, its parameters (name -> value), SOURCE.md lines, children by slot. ``linked``
    (the name of a top-level part of the same workspace) makes it a pointer to that part, as ``IsLinked`` in
    ``rpp_control``'s cascade tree."""
    return {"plugin": plugin, "params": params or {}, "source": source or [], "children": children or {},
            "linked": linked}


def _entry(name, value, kind, place):
    return f"{name} | {json.dumps(value)} | {kind} | {place}"


def remus_lines(values):
    """SOURCE.md value lines of a REMUS 100 part: value | kind | place, from the frozen parameter file."""
    lines = []
    for name, value in values.items():
        key = {"diameter": "beam"}.get(name, name)
        if name == "latitude":
            lines.append(_entry(name, value, "published", REMUS_FILE["gravity"]["place"]))
        elif name == "kinematic_viscosity":
            lines.append(_entry(name, value, "published", "cylinderDrag.m 78-80, nu_water = 1e-6 m^2/s"))
        else:
            entry = REMUS_FILE["parameters"][key]
            lines.append(_entry(name, value, entry["kind"], entry["place"]))
    return lines


def otter_lines(values, place):
    return [_entry(name, value, "published", place[name]) for name, value in values.items()]


OTTER_PLACE = {
    "gravity": "otter.m 90-91", "water_density": "otter.m 90-91",
    "kinematic_viscosity": "cylinderDrag.m 78-80, nu_water = 1e-6 m^2/s",
    "length": "otter.m 92-93, 104-107", "beam": "otter.m 92-93, 104-107", "pontoon_beam": "otter.m 104-107",
    "pontoon_lateral_offset": "otter.m 104-107", "pontoon_block_coefficient": "otter.m 104-107",
    "pontoon_waterplane_coefficient": "otter.m 104-107",
    "hull_mass": "otter.m 94-98", "hull_center_of_gravity": "otter.m 94-98", "radii_of_gyration": "otter.m 94-98",
    "payload_mass": "generate_catamaran_mss.m 61-62", "payload_position": "generate_catamaran_mss.m 61-62",
    "added_mass_coefficients": "otter.m 152-157", "longitudinal_inertia_factor": "otter.m 179",
    "longitudinal_center_of_flotation": "otter.m 177, 192", "reference_point": "the body-frame origin (CO)",
    "max_forward_thrust": "otter.m, 24.4 kgf = 24.4 * 9.81 N", "max_speed": "otter.m, 6 knots = 6 * 0.5144 m/s",
    "time_constants": "otter.m 99-101, 202-240", "damping_ratios": "otter.m 99-101, 202-240",
    "yaw_damping_nonlinearity": "otter.m 99-101, 202-240",
}


def _pick(source, names):
    return {n: source[n] for n in names}


def remus_part(plugin, values, children=None):
    return node(plugin, values, remus_lines(values), children)


def torpedo(coriolis="KirchhoffFull", full=False, open_current=True,
            diagnostics=("mass_matrix", "relative_velocity"), values=None, producers=None):
    """REMUS 100: the hull alone (stage 1, a prescribed wrench) or with fin pairs and propeller (``full``).
    ``values`` replaces the frozen parameter set (a perturbed set in a test); ``producers`` replaces the
    force-producer list."""
    r = REMUS if values is None else values
    current_names = ["current_speed", "current_direction", "current_vertical_speed"]
    open_inputs = [f"current.{n}" for n in current_names] if open_current else []
    if full:
        producers = [
            remus_part("FinPairsDeflectionOnly", _pick(r, (
                "rudder_area", "stern_plane_area", "rudder_lift_coefficient", "stern_plane_lift_coefficient",
                "rudder_position", "stern_plane_position", "max_deflection"))),
            remus_part("PropellerLinearizedOpenWater", _pick(r, (
                "propeller_diameter", "max_shaft_speed", "thrust_deduction", "wake_fraction",
                "pitch_diameter_ratio", "blade_area_ratio", "blade_count", "max_advance_number",
                "roll_moment_scale", "position", "orientation", "thrust_torque_coefficients")))]
    elif producers is None:
        producers = [node("PrescribedWrench")]
    current = node("HorizontalCurrentYawRateTerms",
                   {"open_parameters": current_names} if open_current else {"current_speed": 0.0,
                                                                            "current_direction": 0.0,
                                                                            "current_vertical_speed": 0.0},
                   [])
    return node("MarineCraft6DOF",
                {"integration_max_step": 0.05, "open_inputs": open_inputs, "diagnostic_outputs": list(diagnostics)},
                ["integration_max_step | 0.05 | published | SIMremus100.m 60 (h = 0.05 s)"],
                {"site": remus_part("SiteAtLatitude", {"latitude": REMUS_FILE["gravity"]["latitude_rad"],
                                                       "water_density": r["water_density"],
                                                       "kinematic_viscosity": KINEMATIC_VISCOSITY}),
                 "current": current,
                 "hull_form": remus_part("ProlateSpheroidMainDimensions", {"length": r["length"],
                                                                           "diameter": r["beam"]}),
                 "rigid_body": remus_part("HomogeneousSpheroid", _pick(r, ("body_density",
                                                                           "body_center_of_gravity"))),
                 "added_mass": remus_part("LambSpheroid", _pick(r, ("roll_added_inertia_ratio",))),
                 "added_mass_coriolis": node(coriolis),
                 "restoring": remus_part("SubmergedNeutral", _pick(r, ("center_of_buoyancy",))),
                 "hydrodynamic_loads": [
                     remus_part("TimeConstantDampingSubmerged", _pick(r, ("time_constants", "damping_ratios"))),
                     remus_part("HullLiftDrag", _pick(r, ("planform_area", "parasitic_drag_coefficient",
                                                          "oswald_efficiency"))),
                     node("CrossFlowStrip", {}, [], {"section": node("CircularCylinderReynolds")})],
                 "force_producers": producers})


def catamaran(current=False, diagnostics=("mass_matrix", "rigid_body_mass_matrix", "added_mass_matrix",
                                          "rigid_body_coriolis_matrix", "added_mass_coriolis_matrix",
                                          "restoring_matrix", "hydrodynamic_loads.0.damping_matrix",
                                          "hydrodynamic_loads.0.hydrodynamic_force",
                                          "hydrodynamic_loads.1.hydrodynamic_force")):
    """Otter: no force producers (its hull blocks are the gate), with or without the rigid-body test's current."""
    o = OTTER

    def part(plugin, names, extra=None):
        values = {**_pick(o, names), **(extra or {})}
        return node(plugin, values, otter_lines(values, OTTER_PLACE))

    if current:
        current_node = node("HorizontalCurrentFullRotationRate", {
            "current_speed": CURRENT_SPEED, "current_direction": CURRENT_DIRECTION, "current_vertical_speed": 0.0},
            ["0.3 m/s at 30 deg | test fixture | tests/vehicles/hull_parts/rigid_body/test_rigid_body_block.py CURRENT_SPEED, "
             "CURRENT_DIRECTION (matlab_reference_mss_current.csv)"])
    else:
        current_node = node("NoCurrent")
    return node("MarineCraft6DOF",
                {"integration_max_step": 0.05, "open_inputs": [], "diagnostic_outputs": list(diagnostics)},
                ["integration_max_step | 0.05 | published | SIMremus100.m 60 (h = 0.05 s), the step of the "
                 "MSS vehicle simulations"],
                {"site": part("SiteGiven", ("gravity", "water_density"), {"kinematic_viscosity": KINEMATIC_VISCOSITY}),
                 "current": current_node,
                 "hull_form": part("TwinPontoons", ("length", "beam", "pontoon_beam", "pontoon_lateral_offset",
                                                    "pontoon_block_coefficient", "pontoon_waterplane_coefficient")),
                 "rigid_body": part("HullWithPointPayload", ("hull_mass", "payload_mass", "hull_center_of_gravity",
                                                             "payload_position", "radii_of_gyration")),
                 "added_mass": part("ScaledDerivatives", ("added_mass_coefficients",)),
                 "added_mass_coriolis": node("KirchhoffFull"),
                 "restoring": part("TwinHullMetacentricEquilibriumDraft",
                                   ("longitudinal_inertia_factor", "longitudinal_center_of_flotation"),
                                   {"reference_point": [0.0, 0.0, 0.0]}),
                 "hydrodynamic_loads": [
                     part("TimeConstantDampingSurface", ("max_forward_thrust", "max_speed", "time_constants",
                                                         "damping_ratios", "yaw_damping_nonlinearity")),
                     node("CrossFlowStrip", {}, [], {"section": node("RectangularSectionHoerner")})],
                 "force_producers": []})


def lukas_parts_in_catamaran():
    """Luka's committed LinearSurfaceHydrostatics and LinearSurfaceHydrodynamics (his defaults) in the restoring
    slot and the hull-load list: his input names ``pose`` and ``velocity`` resolve by name."""
    tree = catamaran(False)
    tree["params"]["diagnostic_outputs"] = ["restoring_force", "displaced_volume",
                                            "hydrodynamic_loads.0.hydrodynamic_force"]
    tree["children"]["restoring"] = node("LinearSurfaceHydrostatics")
    tree["children"]["hydrodynamic_loads"] = [node("LinearSurfaceHydrodynamics")]
    return tree


def jet_nozzle_in_catamaran():
    """Luka's JetNozzle (two states, two command inputs; his defaults) in the force-producer list."""
    tree = catamaran(False)
    tree["params"]["diagnostic_outputs"] = ["force_producers.0.generated_thrust"]
    tree["children"]["force_producers"] = [node("JetNozzle")]
    return tree


def _parameters_py(params):
    lines = ["from __future__ import annotations", "", "", "class ComponentParameters:"]
    if not params:
        lines.append("    pass")
    for name, value in params.items():
        lines.append(f"    {name} = {json.dumps(value)}".replace("true", "True").replace("false", "False"))
    return "\n".join(lines) + "\n"


def _write(folder, tree_name, path, n, parent):
    folder.mkdir(parents=True)
    cid = str(uuid.uuid5(NAMESPACE, f"{tree_name}/{path}"))
    subcomponents = {}
    for slot, child in n["children"].items():
        items = child if isinstance(child, list) else [child]
        infos = []
        for k, c in enumerate(items):
            cpath = f"{path}/{slot}.{k}"
            ccid = str(uuid.uuid5(NAMESPACE, f"{tree_name}/{cpath}"))
            infos.append({"Id": ccid, "PluginType": f"more_dynamics::{TYPES[c['plugin']]}",
                          "PluginName": f"more_dynamics::{c['plugin']}", "SlotName": slot,
                          "Library": "more_dynamics", "IsLinked": c["linked"] is not None})
            parent_info = {"Id": cid, "PluginType": f"more_dynamics::{TYPES[n['plugin']]}",
                           "PluginName": f"more_dynamics::{n['plugin']}", "SlotName": slot,
                           "Library": "more_dynamics", "IsLinked": False}
            sub = folder / "subcomponents" / ccid
            if c["linked"] is not None:  # a pointer to a top-level part (rpp resolves it one level below the top)
                sub.mkdir(parents=True)
                (sub / "description.json").write_text(json.dumps({
                    "Id": ccid, "Name": f"{c['plugin']}_comp", "LinkedComponentId": top_id(c["linked"]),
                    "LinkedComponentWorkspace": "more_dynamics", "ParentComponentInfo": parent_info},
                    indent=4) + "\n")
            else:
                _write(sub, tree_name, cpath, c, parent_info)
        subcomponents[slot] = infos if slot in LIST_SLOTS else infos[0]
    description = {"Id": cid, "Name": f"{n['plugin']}_comp" if parent else tree_name,
                   "PluginType": f"more_dynamics::{TYPES[n['plugin']]}", "PluginName": f"more_dynamics::{n['plugin']}",
                   "Library": "more_dynamics", "SubcomponentSpec": SPEC.get(n["plugin"], {}),
                   "Subcomponents": subcomponents, "ParentComponentInfo": parent}
    (folder / "description.json").write_text(json.dumps(description, indent=4) + "\n")
    (folder / "params").mkdir()
    (folder / "params" / "parameters.py").write_text(_parameters_py(n["params"]))
    if parent is None:
        (folder / "callbacks.py").write_text("from __future__ import annotations\n")
    label = LABELS.get(n["plugin"])
    lines = ([f"{n['plugin']}: {label}"] if label else [f"{n['plugin']}: Luka's committed plugin, his defaults"])
    lines += ["", "value | kind | place" if n["source"] else ""] + n["source"]
    (folder / "SOURCE.md").write_text("\n".join(line for line in lines if line is not None).rstrip() + "\n")
    return cid


def top_id(tree_name):
    return str(uuid.uuid5(NAMESPACE, f"{tree_name}/"))


def top(workspace, tree_name, n):
    cid = top_id(tree_name)
    _write(workspace / "parts" / plugin_id(n["plugin"]) / cid, tree_name, "", n, None)
    return cid


def write_workspace(workspace, trees, script_path, descriptions, parts=None):
    """Write ``trees`` ({name: node}) under ``workspace`` and one script description per entry of ``descriptions``
    ({file name: [tree names]}); the first tree of a description is its active configuration. ``parts`` ({name:
    node}) are top-level parts that no configuration builds by itself: the targets of linked children."""
    workspace = Path(workspace)
    if workspace.exists():
        shutil.rmtree(workspace)
    (workspace / "script_descriptions").mkdir(parents=True)
    ids = {name: top(workspace, name, n) for name, n in {**trees, **(parts or {})}.items()}
    for file_name, names in descriptions.items():
        configurations = {name: {"Description": name, "Components": {"vessels": [
            {"Id": ids[name], "PluginName": "more_dynamics::MarineCraft6DOF"}]}} for name in names}
        (workspace / "script_descriptions" / f"{file_name}.json").write_text(json.dumps({
            "ScriptPath": script_path, "Language": "python", "Configurations": configurations,
            "ActiveConfiguration": names[0], "Spec": {"vessels": "List[more_dynamics::VehicleModel3D]"}},
            indent=4) + "\n")
    return ids


FIN_DIAGNOSTICS = ("force_producers.1.generated_force", "force_producers.1.deflection",
                   "force_producers.1.angle_of_attack")


def _remus_place(name):
    return REMUS_FILE["parameters"][name]["place"]


def fin_pair_set(r=None):
    """The rudder and the stern plane of REMUS 100 as two ``LiftingFin`` compositions (servo ``ideal``,
    translational inflow, no flow angle, no interference, quadratic-drag section) under one ``ForceProducerSet``
    with an identity command map: the deflection-only fin pairs of ``remus100.m`` 228-245 in the fin skeleton."""
    r = REMUS if r is None else r

    def fin(prefix, lift_axis):
        servo = node("Servo", {"dynamics": "none", "rate_limit": False, "angle_limit": "on_command",
                               "max_deflection": r["max_deflection"]},
                     [_entry("max_deflection", r["max_deflection"], "published", _remus_place("max_deflection"))])
        slope = r[f"{prefix}_lift_coefficient"]
        section = node("FinSectionQuadraticDrag", {"lift_slope": slope},
                       [_entry("lift_slope", slope, "published", _remus_place(f"{prefix}_lift_coefficient"))])
        position = [r[f"{prefix}_position"], 0.0, 0.0]
        area = r[f"{prefix}_area"]
        return node("LiftingFin", {"fin_position": position, "chord_axis": [1.0, 0.0, 0.0],
                                   "lift_axis": lift_axis, "fin_area": area},
                    [_entry("fin_position", position, "derived", _remus_place(f"{prefix}_position")),
                     "chord_axis, lift_axis | convention | body axes: chord along x, lift along the stated axis "
                     "(a rudder lifts along -y, a stern plane along -z in remus100.m 238-254)",
                     _entry("fin_area", area, "published", _remus_place(f"{prefix}_area"))],
                    {"servo": servo, "inflow": node("FinInflowTranslational"),
                     "flow_angle": node("FinFlowAngleNone"), "interference": node("FinInterferenceNone"),
                     "section": section})

    return node("ForceProducerSet", {"command_count": 2, "command_map": [[1.0, 0.0], [0.0, 1.0]]},
                ["command_count, command_map | convention | one command per fin, identity map"],
                {"producers": [fin("rudder", [0.0, -1.0, 0.0]), fin("stern_plane", [0.0, 0.0, -1.0])]})


def prestero_fin():
    """One REMUS 100 fin at the fin post with the servo of the plugin's defaults, the rigid-point inflow, the
    small-angle flow angle, no interference and the quadratic-drag section: the fin of the linked-part gate."""
    prestero = json.loads((HERE.parent / "force_producers" / "fin_parts" / "prestero_2001_remus_fins.json")
                          .read_text())["parameters"]

    def line(name, key):
        return f"{name} | {prestero[key]['value']} | published | Prestero 2001, {prestero[key]['where']}"

    position = [prestero["fin_position_x"]["value"], 0.0, 0.0]
    return node("LiftingFin", {"fin_position": position, "chord_axis": [1.0, 0.0, 0.0], "lift_axis": [0.0, 1.0, 0.0],
                               "fin_area": prestero["fin_area"]["value"]},
                [line("fin_position x", "fin_position_x"), line("fin_area", "fin_area"),
                 "chord_axis, lift_axis | convention | chord along x, a rudder lifting along +y"],
                {"servo": node("Servo"), "inflow": node("FinInflowRigidPoint"),
                 "flow_angle": node("FinFlowAngleSmallAngle"), "interference": node("FinInterferenceNone"),
                 "section": node("FinSectionQuadraticDrag", {"lift_slope": prestero["lift_slope"]["value"]},
                                 [line("lift_slope", "lift_slope")])})


def fin_vehicle_trees():
    """The REMUS hull with its fins as ``LiftingFin`` compositions: the deflection-only pairs in a set, next to
    the same pairs as ``FinPairsDeflectionOnly``; and one lag-servo fin written inline and as a linked part."""
    trees = {
        "remus100_fin_pairs": torpedo("MunkCouplingsRemoved", diagnostics=("force_producers.0.generated_force",),
                                      producers=[remus_part(
            "FinPairsDeflectionOnly", _pick(REMUS, (
                "rudder_area", "stern_plane_area", "rudder_lift_coefficient", "stern_plane_lift_coefficient",
                "rudder_position", "stern_plane_position", "max_deflection")))]),
        "remus100_fin_set": torpedo("MunkCouplingsRemoved", diagnostics=("force_producers.0.generated_force",),
                                    producers=[fin_pair_set()]),
        "remus100_hull_fin_inline": torpedo("MunkCouplingsRemoved", diagnostics=FIN_DIAGNOSTICS,
                                            producers=[node("PrescribedWrench"), prestero_fin()]),
        "remus100_hull_fin_linked": torpedo("MunkCouplingsRemoved", diagnostics=FIN_DIAGNOSTICS, producers=[
            node("PrescribedWrench"), node("LiftingFin", linked="lifting_fin")]),
    }
    return trees, {"lifting_fin": prestero_fin()}


def gate_trees():
    refused_open = torpedo("MunkCouplingsRemoved")
    refused_open["params"]["open_inputs"] = refused_open["params"]["open_inputs"][:2] + [
        "current.current_vertical_sped"]
    twin_pontoons = node("TwinPontoons", _pick(OTTER, (
        "length", "beam", "pontoon_beam", "pontoon_lateral_offset", "pontoon_block_coefficient",
        "pontoon_waterplane_coefficient")))
    refused_coupling = torpedo("MunkCouplingsRemoved")
    refused_coupling["children"]["hull_form"] = twin_pontoons  # no spheroid semi-axes for the Lamb added mass
    return {
        "remus100_hull": torpedo("MunkCouplingsRemoved"),
        "remus100_hull_kirchhoff": torpedo("KirchhoffFull"),
        "remus100_full": torpedo("MunkCouplingsRemoved", full=True),
        "remus100_full_kirchhoff": torpedo("KirchhoffFull", full=True),
        "otter": catamaran(False),
        "otter_current": catamaran(True),
        "otter_lukas_parts": lukas_parts_in_catamaran(),
        "otter_jet_nozzle": jet_nozzle_in_catamaran(),
        "refused_missing_coupling": refused_coupling,
        "refused_open_input": refused_open,
        **fin_vehicle_trees()[0],
    }


def gate_parts():
    """The top-level parts the gate trees link to."""
    return fin_vehicle_trees()[1]


def user_trees():
    return {
        "remus100": torpedo("KirchhoffFull", full=True, open_current=False, diagnostics=()),
        "otter": catamaran(False, diagnostics=()),
    }


def main():
    gate = gate_trees()
    write_workspace(GATE_TREES, gate, "tests/vehicles/vehicle_contract.py", {"vehicles": list(gate)},
                    parts=gate_parts())
    user = user_trees()
    write_workspace(USER_TREES, user, "scripts/vehicles/simulate_vehicle.py",
                    {name: [name] for name in user})
    print(f"wrote {len(gate)} gate trees under {GATE_TREES} and {len(user)} user vehicles under {USER_TREES}")


if __name__ == "__main__":
    main()
