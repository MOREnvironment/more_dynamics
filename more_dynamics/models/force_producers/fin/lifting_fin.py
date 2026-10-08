"""One lifting fin assembled from five independently selected parts.

The force is resolved along fixed BODY FRD fin axes; this is the small-angle
form used by the cited sources, not a large-angle flow-axis model.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics and Motion Control, 1st ed., Wiley, eq. 12.226, p. 400.
[Prestero 2001] Prestero, T. (2001). Verification of a six-degree of freedom simulation model for the REMUS autonomous underwater vehicle. MIT/WHOI MSc thesis, eqs. 4.37, 4.41-4.43, pp. 31-33.
[MSS] Fossen, T. I. (2026). Marine Systems Simulator, MIT, CRAFT/AUV/models/remus100.m 238-254 @ cc07579.

Author:    Enio Krizman
Date:      2026-10-08
"""
import casadi as ca
from more_transformations.more_casadi_transformations import Parameter, symbols, check_values
from .._common import wrench


_COUPLINGS = {
    "servo": ({"command": 1, "servo_state": None}, {"deflection": 1, "servo_state_dot": None}),
    "inflow": ({"nu_r": 6, "fin_position": 3}, {"fin_velocity": 3}),
    "flow_angle": ({"fin_velocity": 3, "chord_axis": 3, "lift_axis": 3}, {"flow_angle": 1, "speed_squared": 1}),
    "interference": ({}, {"deflection_factor": 1, "flow_angle_factor": 1}),
    "section": ({"angle_of_attack": 1}, {"lift_coefficient": 1, "drag_coefficient": 1}),
}


def lifting_fin_parameters():
    return (Parameter("fin_position", (3, 1), "m", "CO to fin in BODY FRD"),
            Parameter("chord_axis", (3, 1), "1", "unit chord direction in BODY FRD"),
            Parameter("lift_axis", (3, 1), "1", "unit positive-lift direction in BODY FRD"),
            Parameter("fin_area", (1, 1), "m^2", "planform area", 0.0, minimum_exclusive=True))


def check_lifting_fin_values(values):
    checked = check_values(lifting_fin_parameters(), values)
    chord, lift = checked["chord_axis"], checked["lift_axis"]
    for name, axis in (("chord_axis", chord), ("lift_axis", lift)):
        if abs(float(ca.norm_2(axis)) - 1.0) > 1e-9:
            raise ValueError(f"{name} must be a unit vector")
    if abs(float(ca.dot(chord, lift))) > 1e-9:
        raise ValueError("lift_axis must be perpendicular to chord_axis")
    return checked


def _validate(parts):
    for slot in parts:
        if slot not in _COUPLINGS:
            raise ValueError(f"unknown slot {slot}")
    for slot, (inputs, outputs) in _COUPLINGS.items():
        if slot not in parts:
            raise ValueError(f"missing slot {slot}")
        part = parts[slot]
        if not isinstance(part, ca.Function):
            raise ValueError(f"slot {slot} requires a CasADi function")
        if "water_density" in part.name_in():
            raise ValueError(f"slot {slot}: water_density belongs to the vehicle feed")
        for names, expected, direction in ((part.name_in(), inputs, "in"), (part.name_out(), outputs, "out")):
            for name, rows in expected.items():
                if name not in names:
                    raise ValueError(f"slot {slot}: missing {name}")
                size = part.size_in(name) if direction == "in" else part.size_out(name)
                if rows is not None and size != (rows, 1):
                    raise ValueError(f"slot {slot}: {name} must have {rows} rows")
        if slot == "servo" and part.size1_in("servo_state") != part.size1_out("servo_state_dot"):
            raise ValueError("slot servo: state and state_dot sizes differ")


def lifting_fin_casadi(parts):
    _validate(parts)
    p = symbols(lifting_fin_parameters())
    command = ca.SX.sym("command")
    state = ca.SX.sym("state", parts["servo"].size1_in("servo_state"))
    nu_r, density = ca.SX.sym("nu_r", 6), ca.SX.sym("water_density")
    couplings = {"command": command, "servo_state": state, "nu_r": nu_r,
                 "fin_position": p["fin_position"], "chord_axis": p["chord_axis"],
                 "lift_axis": p["lift_axis"]}
    open_symbols = []
    open_names = []

    def invoke(slot):
        part = parts[slot]
        args = {}
        for name in part.name_in():
            if name in _COUPLINGS[slot][0]:
                args[name] = couplings[name]
            else:
                key = f"{slot}.{name}"
                sym = ca.SX.sym(key, *part.size_in(name))
                args[name] = sym
                open_names.append(key)
                open_symbols.append(sym)
        return part(**args)

    servo = invoke("servo")
    couplings["fin_velocity"] = invoke("inflow")["fin_velocity"]
    angle = invoke("flow_angle")
    factors = invoke("interference")
    alpha = factors["deflection_factor"] * servo["deflection"] - factors["flow_angle_factor"] * angle["flow_angle"]  # (Prestero 2001, eq. 4.41, p. 32; interference factors equal one in cited form)
    couplings["angle_of_attack"] = alpha
    section = invoke("section")
    force = 0.5 * density * p["fin_area"] * angle["speed_squared"] * (section["lift_coefficient"] * p["lift_axis"] - section["drag_coefficient"] * p["chord_axis"])  # (Prestero 2001, eqs. 4.37, 4.43, pp. 31-33; MSS remus100.m 238-245)
    tau = wrench(force, p["fin_position"])  # (Fossen 2011, eq. 12.226, p. 400)
    geometry = [d.name for d in lifting_fin_parameters()]
    return ca.Function("lifting_fin", [command, state, nu_r, density, *[p[n] for n in geometry], *open_symbols],
                       [tau, servo["servo_state_dot"], servo["deflection"], alpha],
                       ["command", "state", "nu_r", "water_density", *geometry, *open_names],
                       ["tau", "state_dot", "deflection", "angle_of_attack"])
