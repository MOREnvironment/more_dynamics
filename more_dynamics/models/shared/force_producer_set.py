"""Sum a list of force producers behind an explicit command map.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics and Motion Control, 1st ed., Wiley, eq. 12.226, p. 400 (superposition of applied wrenches).

Author:    Enio Krizman
Date:      2026-10-08
"""
import casadi as ca
from more_transformations.more_casadi_transformations import Parameter


def _layout(producers):
    rows = 0
    states = 0
    for i, producer in enumerate(producers):
        for name, size in (("command", None), ("state", None), ("nu_r", 6), ("water_density", 1)):
            if name not in producer.name_in():
                raise ValueError(f"producers[{i}]: missing {name}")
            if size is not None and producer.size_in(name) != (size, 1):
                raise ValueError(f"producers[{i}]: invalid {name}")
        for name, size in (("tau", 6), ("state_dot", producer.size1_in("state"))):
            if name not in producer.name_out() or producer.size_out(name) != (size, 1):
                raise ValueError(f"producers[{i}]: invalid {name}")
        rows += producer.size1_in("command")
        states += producer.size1_in("state")
    return rows, states


def force_producer_set_parameters(producers, command_count):
    rows, _ = _layout(producers)
    if not isinstance(command_count, int) or command_count < 1:
        raise ValueError("command_count must be a positive integer")
    return (Parameter("command_map", (rows, command_count), "1", "producer commands per vehicle command"),)


def force_producer_set_casadi(producers, command_count):
    declared = force_producer_set_parameters(producers, command_count)
    rows, states = _layout(producers)
    command, state = ca.SX.sym("command", command_count), ca.SX.sym("state", states)
    nu_r, density = ca.SX.sym("nu_r", 6), ca.SX.sym("water_density")
    command_map = ca.SX.sym("command_map", rows, command_count)
    mapped = command_map @ command  # command allocation by the composition's map
    free_inputs, free_names = [], []
    forces, rates = [], []
    command_offset = state_offset = 0
    for i, producer in enumerate(producers):
        nc, ns = producer.size1_in("command"), producer.size1_in("state")
        args = {"command": mapped[command_offset:command_offset + nc], "state": state[state_offset:state_offset + ns],
                "nu_r": nu_r, "water_density": density}
        for name in producer.name_in()[4:]:
            prefixed = f"producers[{i}].{name}"
            sym = ca.SX.sym(prefixed, *producer.size_in(name))
            free_inputs.append(sym)
            free_names.append(prefixed)
            args[name] = sym
        out = producer(**args)
        forces.append(out["tau"])
        rates.append(out["state_dot"])
        command_offset += nc
        state_offset += ns
    tau = sum(forces, ca.SX.zeros(6, 1))  # wrench superposition (Fossen 2011, eq. 12.226, p. 400)
    state_dot = ca.vertcat(*rates) if rates else ca.SX.zeros(0, 1)
    return ca.Function("force_producer_set", [command, state, nu_r, density, command_map, *free_inputs],
                       [tau, state_dot], ["command", "state", "nu_r", "water_density", "command_map", *free_names],
                       ["tau", "state_dot"])
