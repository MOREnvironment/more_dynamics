"""The parameter sets of the fidelity-option gates, read from the files beside
this one and from the one-value REMUS 100 reference (``../remus100/``): the
named part "REMUS 100 (Munk + body lift)" and the comparison sets built from
the same base, each a flat set by the plugins' parameter names
(``reference_torpedo_auv.TorpedoAuv``'s input).

Cases (one mechanism changed at a time against the MSS forms; the pressure
centre at 0.60, 0.65 and 0.70 of the length from the nose, Prestero 2001,
p. 31):

* ``mss_forms`` -- every MSS form, the fin pairs deflection-only with MSS's
  calibrated lift coefficients (remus100.m 179-189, 230-245);
* ``munk_kept_alone`` -- the same with every added-mass Coriolis term kept;
* ``geometric_fin_slope_effective_angle`` -- MSS forms, four effective-angle
  fins with the lift slope from the fin geometry;
* ``paired_mss_fin_slope`` -- Munk kept, point-force body lift,
  effective-angle fins with MSS's calibrated coefficients (expected
  straight-line unstable: the gate must be able to fail);
* ``paired_geometric_fin_slope`` -- the named part's forms.

References
----------
Prestero, T. (2001). Verification of a six-degree of freedom simulation
    model for the REMUS autonomous underwater vehicle. MSc thesis, MIT/WHOI.
    Table 2.2, p. 17; eqs. 4.35, 4.38-4.39, p. 31.
Fossen, T. I. MSS, MIT, @ cc07579: CRAFT/AUV/models/remus100.m 109,
    179-189, 207-210, 230-245.

Author:    Enio Krizman
Date:      2026-10-13
"""

import csv
import json
from pathlib import Path

from reference_torpedo_auv import low_aspect_ratio_lift_slope

HERE = Path(__file__).resolve().parent
PART_FILE = HERE / "remus100_munk_body_lift_parameters.json"
PRESSURE_CENTRE_FRACTIONS = (0.60, 0.65, 0.70)   # Prestero 2001, p. 31 (Hoerner's range)
KINEMATIC_VISCOSITY = 1e-6                       # cylinderDrag.m 81


def part_file():
    return json.loads(PART_FILE.read_text())


def base_values():
    """The one-value MSS set by the plugins' names (the part file's ``base``)."""
    record = json.loads((HERE / part_file()["base"]["file"]).read_text())
    v = {name: entry["value"] for name, entry in record["parameters"].items()}
    v["diameter"] = v["beam"]
    v["body_mass"] = record["one_value"]["mass"]
    v["latitude"] = record["gravity"]["latitude_rad"]
    v["kinematic_viscosity"] = KINEMATIC_VISCOSITY
    return v


def cruise_shaft_speed():
    """The cruise shaft speed of the part (its ``cruise`` entry: a frozen MSS reference cell)."""
    cell = part_file()["cruise"]["shaft_speed"]
    with open(HERE / cell["file"], newline="") as handle:
        rows = list(csv.DictReader(handle))
    return float(rows[cell["row"]][cell["column"]])


def _propeller(v):
    names = ("propeller_diameter", "max_shaft_speed", "thrust_deduction", "wake_fraction",
             "max_advance_number", "roll_moment_scale", "thrust_torque_coefficients")
    return {n: v[n] for n in names}


def _mss_forms(v):
    return {"coriolis_form": "munk_couplings_removed", "body_lift_form": "vertical_plane_wing",
            "hull_drag_form": "wing_drag", "roll_damping_form": "linear"}


def mss_fin_pairs(v, inflow="translational", flow_angle="none"):
    """MSS's two fin pairs as two fins of the pair area (remus100.m 179-189, 230-245)."""
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
                                        fin["section_lift_efficiency"]["value"])
    fins = []
    for entry in record["fins"]:
        rudder = entry["role"] == "rudder"
        fins.append({
            "fin_position": [-v["length"] / 2.0, 0, 0],                  # remus100.m 184, 189
            "lift_axis": entry["lift_axis"], "chord_axis": fin["chord_axis"]["value"],
            "fin_area": fin["fin_area"]["value"], "max_deflection": v[fin["max_deflection"]["base"]],
            "inflow": fin["inflow"]["value"], "flow_angle": fin["flow_angle"]["value"],
            "lift_slope": slope if lift_slopes is None else lift_slopes[0 if rudder else 1]})
    return fins


def munk_body_lift(fraction=None, lift_slopes=None):
    """The named part (``fraction`` replaces its pressure-centre fraction; ``lift_slopes`` its fin slope)."""
    v = base_values()
    record = part_file()
    part = {**v, **{k: e["value"] for k, e in record["vehicle"].items()},
            **{k: e["value"] for k, e in record["hydrodynamics"].items()}}
    if fraction is not None:
        part["pressure_centre_fraction"] = fraction
    part["fins"] = part_fins(v, lift_slopes)
    part["propeller"] = _propeller(v)
    return part


def mss_forms():
    v = base_values()
    return {**v, **_mss_forms(v), "fins": mss_fin_pairs(v), "propeller": _propeller(v)}


def munk_kept_alone():
    return {**mss_forms(), "coriolis_form": "kirchhoff_full"}


def geometric_fin_slope_effective_angle():
    v = base_values()
    return {**v, **_mss_forms(v), "fins": part_fins(v), "propeller": _propeller(v)}


def paired_mss_fin_slope(fraction):
    """The part's forms with MSS's calibrated fin coefficients on each fin (the pair's area is two fins)."""
    v = base_values()
    return munk_body_lift(fraction, lift_slopes=(v["rudder_lift_coefficient"], v["stern_plane_lift_coefficient"]))


def cases():
    """``(name, fraction or None, parameter set)`` of every reference case."""
    out = [("mss_forms", None, mss_forms()), ("munk_kept_alone", None, munk_kept_alone()),
           ("geometric_fin_slope_effective_angle", None, geometric_fin_slope_effective_angle())]
    out += [("paired_mss_fin_slope", f, paired_mss_fin_slope(f)) for f in PRESSURE_CENTRE_FRACTIONS]
    out += [("paired_geometric_fin_slope", f, munk_body_lift(f)) for f in PRESSURE_CENTRE_FRACTIONS]
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
