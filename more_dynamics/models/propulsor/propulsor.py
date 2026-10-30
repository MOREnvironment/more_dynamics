"""One propulsor from six typed parts: drive, steering, mounting, inflow, law and housing.

Each part is one CasADi function whose inputs are the couplings of its type, then its own declared numbers,
and whose outputs are the couplings of its type. The skeleton calls the parts in this order and never
knows which form a part has:

1. steering: ``command`` and ``servo_state`` -> the deflection ``d`` and the servo's derivative;
2. mounting: ``d`` -> the thrust axis ``a`` (BODY, unit or not as the part says) and the thrust point ``r``;
3. inflow: ``nu_r``, ``a``, ``r`` -> the advance speed ``Va``, the thrust factor ``k_T`` and the torque
   factor ``k_Q``;
4. drive: ``command``, ``drive_state``, load torque -> the shaft output ``n``, the state derivative, the thrust
   scale ``s`` and the thrust limit ``T_lim``; called first to get ``n`` (which must not depend on the load: a
   speed that depends on the load is a state of the drive), then again with the load ``k_Q Q``;
5. law: ``n``, ``Va``, density -> thrust ``T`` and shaft torque ``Q``;
6. housing: ``nu_r``, ``a``, ``r``, ``d``, density -> an extra wrench.

``thrust = sat(s k_T T, +-T_lim)`` (the thrust along the axis after deduction, power cap and limit) and
``tau = [thrust a; r x thrust a] + [0; c k_Q Q a] + housing`` with the shaft-torque scale ``c``: the
wrench of a point force, ``[f; r x f]`` (Fossen 2011, eq. 12.226, p. 400), plus the shaft torque as a
couple about the axis (MSS ``remus100.m`` 250-253, where the scale is 1/10; the physical reaction has
scale 1). Several propulsors on one vehicle are several force producers of the vehicle's own list.

Inputs of the function: ``command`` (the drive's command, then the steering's), ``state`` (the drive's,
then the steering's), ``nu_r`` (6, relative to the water, BODY), ``water_density`` (fed by the vehicle),
``shaft_torque_scale`` and every still-open part number as ``"<slot>.<name>"``. Outputs: ``tau``,
``state_dot``, ``thrust`` and ``steering_angle``.

A drive and a law each name the unit of what they pass (``units``): the shaft output of a speed drive is
in rev/s, a normalised effort has unit 1. A pair in different units is refused.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics and Motion Control. Wiley.
    Ch. 12, eq. 12.226, p. 400.
[MSS] Fossen, T. I. (2026). Marine Systems Simulator, MIT, CRAFT/AUV/models/remus100.m 250-253
    @ 6e78e98.
[Mandic] Mandic, L. rpp_control, ``plugins/controllers/cascade_controller.hpp`` (the form: a skeleton
    with typed slots that holds parts it does not look into).

Author:    Enio Krizman
Date:      2026-10-23
"""
import casadi as ca
from more_transformations.more_casadi_transformations import Parameter, check_values

from more_dynamics.models.shared.force_producer_common import saturate, wrench

SLOTS = ("drive", "steering", "mounting", "inflow", "law", "housing")
DENSITY_SLOTS = ("law", "housing")

# slot -> (coupling inputs {name: rows, None = any}, coupling outputs {name: rows, None = any}), in this order
COUPLINGS = {
    "drive": ({"command": None, "drive_state": None, "shaft_torque": 1},
              {"drive_output": 1, "drive_state_dot": None, "thrust_scale": 1, "thrust_limit": 1}),
    "steering": ({"command": None, "servo_state": None}, {"deflection": 1, "servo_state_dot": None}),
    "mounting": ({"steering_angle": 1}, {"thrust_axis": 3, "thrust_position": 3}),
    "inflow": ({"nu_r": 6, "thrust_axis": 3, "thrust_position": 3},
               {"advance_speed": 1, "thrust_factor": 1, "torque_factor": 1}),
    "law": ({"drive_output": 1, "advance_speed": 1, "water_density": 1}, {"thrust": 1, "torque": 1}),
    "housing": ({"nu_r": 6, "thrust_axis": 3, "thrust_position": 3, "steering_angle": 1, "water_density": 1},
                {"tau": 6}),
}
STATE_PAIRS = {"drive": ("drive_state", "drive_state_dot"), "steering": ("servo_state", "servo_state_dot")}


def propulsor_parameters():
    return (Parameter("shaft_torque_scale", (1, 1), "1",
                      "share of the shaft torque that reacts on the body (physical: 1; a comparison value is below it)",
                      0.0),)


def check_propulsor_values(values):
    """The skeleton's number checked (``check_values``): a negative scale is refused, naming it."""
    return check_values(propulsor_parameters(), values)


def _check_part(slot, part):
    if not isinstance(part, ca.Function):
        raise ValueError(f"slot {slot!r}: a part must be a CasADi function, got {type(part).__name__}")
    inputs, outputs = COUPLINGS[slot]
    for kind, expected, names, size in (("input", inputs, part.name_in(), part.size_in),
                                        ("output", outputs, part.name_out(), part.size_out)):
        for name, rows in expected.items():
            if name not in names:
                raise ValueError(f"slot {slot!r}: the part has no {kind} {name!r} (its type's couplings: "
                                 f"{list(inputs)} -> {list(outputs)})")
            if rows is not None and size(name) != (rows, 1):
                raise ValueError(f"slot {slot!r}: {kind} {name!r} must have {rows} rows, got {size(name)[0]}")
    if slot not in DENSITY_SLOTS and "water_density" in part.name_in():
        raise ValueError(f"slot {slot!r}: a part of this slot cannot take 'water_density' (only law and housing do)")
    if slot in STATE_PAIRS:
        state, derivative = STATE_PAIRS[slot]
        if part.size1_in(state) != part.size1_out(derivative):
            raise ValueError(f"slot {slot!r}: {state!r} has {part.size1_in(state)} rows but {derivative!r} has "
                             f"{part.size1_out(derivative)}")


def check_units(units):
    """``units`` is ``{"drive": unit of the drive's output, "law": unit the law accepts}``; a pair in different
    units is refused, naming both."""
    for slot in ("drive", "law"):
        if slot not in units:
            raise ValueError(f"units: no unit given for the {slot}")
    if units["drive"] != units["law"]:
        raise ValueError(f"drive and law units differ: the drive's output is {units['drive']!r}, "
                         f"the law accepts {units['law']!r}")


def propulsor_casadi(parts, units):
    """The propulsor function of ``parts`` (slot -> part function, frozen or with numbers open) and ``units``
    (``{"drive": unit of the drive's output, "law": unit the law accepts}``); see the module docstring."""
    unknown = [slot for slot in parts if slot not in SLOTS]
    if unknown:
        raise ValueError(f"unknown slot {unknown[0]!r}: the slots are {list(SLOTS)}")
    for slot in SLOTS:
        if slot not in parts:
            raise ValueError(f"missing slot {slot!r}: the slots are {list(SLOTS)}")
        _check_part(slot, parts[slot])
    check_units(units)

    p = {d.name: ca.SX.sym(d.name, *d.shape) for d in propulsor_parameters()}
    opened = {}
    for slot in SLOTS:
        for name in parts[slot].name_in():
            if name not in COUPLINGS[slot][0]:
                key = f"{slot}.{name}"
                opened[key] = ca.SX.sym(key, *parts[slot].size_in(name))

    def call(slot, **couplings):
        args = {name: couplings[name] if name in couplings else opened[f"{slot}.{name}"]
                for name in parts[slot].name_in()}
        return parts[slot](**args)

    drive_size, steering_size = parts["drive"].size1_in("command"), parts["steering"].size1_in("command")
    drive_state_size, steering_state_size = parts["drive"].size1_in("drive_state"), parts["steering"].size1_in("servo_state")
    command = ca.SX.sym("command", drive_size + steering_size)
    state = ca.SX.sym("state", drive_state_size + steering_state_size)
    nu_r, density, load = ca.SX.sym("nu_r", 6), ca.SX.sym("water_density"), ca.SX.sym("shaft_torque")

    steering = call("steering", command=command[drive_size:], servo_state=state[drive_state_size:])
    mounting = call("mounting", steering_angle=steering["deflection"])
    axis, point = mounting["thrust_axis"], mounting["thrust_position"]
    inflow = call("inflow", nu_r=nu_r, thrust_axis=axis, thrust_position=point)
    drive = call("drive", command=command[:drive_size], drive_state=state[:drive_state_size], shaft_torque=load)
    if ca.depends_on(drive["drive_output"], load):
        raise ValueError("slot 'drive': the shaft output must not depend on the load torque (a load-dependent "
                         "speed is a state of the drive, not an algebraic loop)")
    law = call("law", drive_output=drive["drive_output"], advance_speed=inflow["advance_speed"], water_density=density)
    shaft_torque = inflow["torque_factor"] * law["torque"]  # load on the drive: k_Q Q
    drive_dot, scale, limit = (ca.substitute(drive[name], load, shaft_torque)
                               for name in ("drive_state_dot", "thrust_scale", "thrust_limit"))
    thrust = saturate(scale * (inflow["thrust_factor"] * law["thrust"]), -limit, limit)
    housing = call("housing", nu_r=nu_r, thrust_axis=axis, thrust_position=point,
                   steering_angle=steering["deflection"], water_density=density)
    tau = (wrench(axis * thrust, point)  # [f; r x f] (Fossen 2011, eq. 12.226, p. 400)
           + ca.vertcat(ca.DM.zeros(3), axis * (p["shaft_torque_scale"] * shaft_torque))  # shaft torque as a couple (remus100.m 253)
           + housing["tau"])
    names = ["command", "state", "nu_r", "water_density", *p, *opened]
    return ca.Function("propulsor", [command, state, nu_r, density, *p.values(), *opened.values()],
                       [tau, ca.vertcat(drive_dot, steering["servo_state_dot"]), thrust, steering["deflection"]],
                       names, ["tau", "state_dot", "thrust", "steering_angle"])
