"""The parameter sets of the fidelity-option gates, read from the files beside
this one and from the one-value REMUS 100 base set (``../remus100/``): the
named part "REMUS 100 (estimated)" and the comparison sets built from the same
base, each a flat set by the plugins' parameter names
(``reference_torpedo_auv.TorpedoAuv``'s input).

Cases (one mechanism changed at a time against the published reference forms;
the pressure centre at 0.60, 0.65 and 0.70 of the length from the nose, Prestero
2001, p. 31):

* ``published_forms`` -- every form of the published reference model: shortcut
  Coriolis (the Munk-moment couplings removed), wing lift and drag on the
  vertical-plane angle, linear roll damping, the fin pairs deflection-only with
  the pair's published lift coefficients;
* ``munk_kept_alone`` -- the same with every added-mass Coriolis term kept;
* ``geometric_fin_slope_effective_angle`` -- published hull forms, four
  effective-angle fins with the lift slope from the fin geometry;
* ``paired_published_fin_slope`` -- Munk kept, point-force body lift,
  effective-angle fins with the published lift coefficients (expected
  straight-line unstable: the gate must be able to fail);
* ``paired_geometric_fin_slope`` -- the named part's forms.

References
----------
Prestero, T. (2001). Verification of a six-degree of freedom simulation
    model for the REMUS autonomous underwater vehicle. MSc thesis, MIT/WHOI.
    Table 2.2, p. 17; eqs. 4.35, 4.38-4.39, p. 31.

Author:    Enio Krizman
Date:      2026-10-13
"""

import json
from pathlib import Path

from reference_torpedo_auv import low_aspect_ratio_lift_slope

HERE = Path(__file__).resolve().parent
PART_FILE = HERE / "remus100_estimated_parameters.json"
PRESSURE_CENTRE_FRACTIONS = (0.60, 0.65, 0.70)   # Prestero 2001, p. 31 (Hoerner's range)
KINEMATIC_VISCOSITY = 1e-6                       # m^2/s, fresh water


def part_file():
    return json.loads(PART_FILE.read_text())


def base_values():
    """The one-value base set by the plugins' names (the part file's ``base``)."""
    record = json.loads((HERE / part_file()["base"]["file"]).read_text())
    v = {name: entry["value"] for name, entry in record["parameters"].items()}
    v["diameter"] = v["beam"]
    v["body_mass"] = record["one_value"]["mass"]
    v["latitude"] = record["gravity"]["latitude_rad"]
    v["kinematic_viscosity"] = KINEMATIC_VISCOSITY
    v.update({name: entry["value"] for name, entry in part_file()["base_additions"].items()})
    return v


def cruise_shaft_speed():
    """The cruise shaft speed of the part (its ``cruise`` entry, rpm): the one that holds the cruise speed straight in the
    reference model with the published forms (``reference_torpedo_auv.shaft_speed_holding``, tested)."""
    return part_file()["cruise"]["shaft_speed"]["value"]


def _propeller(v):
    names = ("propeller_diameter", "max_shaft_speed", "thrust_deduction", "wake_fraction",
             "max_advance_number", "roll_moment_scale", "thrust_torque_coefficients")
    return {n: v[n] for n in names}


def _published_forms(v):
    return {"coriolis_form": "munk_couplings_removed", "body_lift_form": "vertical_plane_wing",
            "hull_drag_form": "wing_drag", "roll_damping_form": "linear"}


def published_fin_pairs(v, inflow="translational", flow_angle="none"):
    """The published reference model's two fin pairs, each as one fin of the pair's area (two fins of the part's fin area)."""
    common = {"chord_axis": [1, 0, 0], "max_deflection": v["max_deflection"], "inflow": inflow, "flow_angle": flow_angle}
    return [
        {**common, "fin_position": [v["rudder_position"], 0, 0], "lift_axis": [0, -1, 0],
         "fin_area": v["rudder_area"], "lift_slope": v["rudder_lift_coefficient"]},
        {**common, "fin_position": [v["stern_plane_position"], 0, 0], "lift_axis": [0, 0, -1],
         "fin_area": v["stern_plane_area"], "lift_slope": v["stern_plane_lift_coefficient"]},
    ]


def part_fins(v, lift_slopes=None):
    """The named part's four fins (its ``fin`` and ``fins`` entries), the lift slope computed by
    ``low_aspect_ratio`` unless ``lift_slopes`` (rudder, stern plane) replaces it."""
    record = part_file()
    fin = record["fin"]
    slope = low_aspect_ratio_lift_slope(fin["fin_span"]["value"], fin["fin_area"]["value"],
                                        fin["section_efficiency"]["value"], fin["aspect_ratio_factor"]["value"])
    fins = []
    for entry in record["fins"]:
        rudder = entry["role"] == "rudder"
        fins.append({
            "fin_position": [-v["length"] / 2.0, 0, 0],                  # at the tail, the hull centred on the origin
            "lift_axis": entry["lift_axis"], "chord_axis": fin["chord_axis"]["value"],
            "fin_area": fin["fin_area"]["value"], "max_deflection": v[fin["max_deflection"]["base"]],
            "inflow": fin["inflow"]["value"], "flow_angle": fin["flow_angle"]["value"],
            "lift_slope": slope if lift_slopes is None else lift_slopes[0 if rudder else 1]})
    return fins


def estimated(fraction=None, lift_slopes=None):
    """The named part (``fraction`` replaces its pressure-centre fraction; ``lift_slopes`` its fin slope)."""
    v = base_values()
    record = part_file()
    part = {**v, **{k: e["value"] for k, e in record["vehicle"].items()},
            **{k: e["value"] for k, e in record["hydrodynamics"].items()}}
    if fraction is not None:
        part["pressure_centre_fraction"] = fraction
    part["fins"] = part_fins(v, lift_slopes)
    part.update({k: e["value"] for k, e in record["propeller"].items()})   # the part's own propeller values, by the plugin's name
    part["propeller"] = _propeller(part)
    return part


def published_forms():
    v = base_values()
    return {**v, **_published_forms(v), "fins": published_fin_pairs(v), "propeller": _propeller(v)}


def munk_kept_alone():
    return {**published_forms(), "coriolis_form": "kirchhoff_full"}


def geometric_fin_slope_effective_angle():
    v = base_values()
    return {**v, **_published_forms(v), "fins": part_fins(v), "propeller": _propeller(v)}


def paired_published_fin_slope(fraction):
    """The part's forms with the published reference model's lift coefficients on each fin (the pair's area is two fins)."""
    v = base_values()
    return estimated(fraction, lift_slopes=(v["rudder_lift_coefficient"], v["stern_plane_lift_coefficient"]))


def cases():
    """``(name, fraction or None, parameter set)`` of every reference case."""
    out = [("published_forms", None, published_forms()), ("munk_kept_alone", None, munk_kept_alone()),
           ("geometric_fin_slope_effective_angle", None, geometric_fin_slope_effective_angle())]
    out += [("paired_published_fin_slope", f, paired_published_fin_slope(f)) for f in PRESSURE_CENTRE_FRACTIONS]
    out += [("paired_geometric_fin_slope", f, estimated(f)) for f in PRESSURE_CENTRE_FRACTIONS]
    return out


def cases_by_name():
    """The cases without a pressure-centre fraction, by name."""
    return {name: part for name, fraction, part in cases() if fraction is None}


def rudder_command(part, rudder, stern_plane=0.0, shaft_speed=None):
    """The vehicle input for a part: one command per fin by its role, then the shaft speed."""
    commands = []
    for fin in part["fins"]:
        commands.append(rudder if abs(fin["lift_axis"][1]) > 0 else stern_plane)
    return commands + [cruise_shaft_speed() if shaft_speed is None else shaft_speed]
