"""A list of force producers behind one command map, as Luka's
``ForceProducer``: the set's wrench is the superposition of its children's
wrenches (Fossen 2011, eq. 12.226, p. 400), the children's states are stacked
in list order, and the set's commands reach the children through
``command_map`` (rows: the children's commands in list order, columns: the
set's commands).

Each child is a ``ForceProducer`` read through its payload: ``relative_velocity``
and ``water_density`` are fed to every child (one density for all), every
other input of a child is a command (the vehicle's rule for a force producer),
and its first output is its wrench. Inputs of the set:
``force_producer_command`` (``command_count`` values), ``relative_velocity``,
``water_density``. Output: ``generated_force``; state: ``producer_state``.

Defaults: one set command driving a single one-command child (a map of one);
every other set declares its own ``command_count`` and ``command_map``.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics
    and Motion Control. Wiley. Eq. 12.226, p. 400 (superposition of applied
    wrenches).

Author:    Enio Krizman
Date:      2026-10-08
"""
import casadi as ca
from rpp_plugin_types.more_dynamics import ForceProducer
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_common.casadi_graph import RppCasadiGraph
from more_transformations.more_casadi_transformations import freeze
from more_dynamics.models.force_producers.force_producer_set import (force_producer_set_casadi,
                                                                    force_producer_set_parameters)
from more_dynamics.plugins.shared.payload_io import PayloadBuilder

FED = {"relative_velocity": 6, "water_density": 1}


def producer_function(child, index):
    """A child force producer as the model-layer set's producer function:
    inputs ``command, state, nu_r, water_density``, outputs ``tau, state_dot``."""
    label = f"producers[{index}] ({type(child).__name__})"
    payload = child.graph()
    graph = RppCasadiGraph(payload)
    sizes = {d.name: d.size for d in payload.inputDescription}
    for name, size in FED.items():
        if sizes.get(name) != size:
            raise ValueError(f"{label}: input {name!r} ({size}) is missing or has size {sizes.get(name)}")
    if not payload.outputDescription or payload.outputDescription[0].size != 6:
        raise ValueError(f"{label}: the first output must be the 6-value generalized force")
    commands = sum(size for name, size in sizes.items() if name not in FED)
    command, state = ca.SX.sym("command", commands), ca.SX.sym("state", graph.num_states)
    nu_r, density = ca.SX.sym("nu_r", 6), ca.SX.sym("water_density")
    args, k = [], 0
    for d in payload.inputDescription:
        if d.name == "relative_velocity":
            args.append(nu_r)
        elif d.name == "water_density":
            args.append(density)
        else:
            args.append(command[k:k + d.size])
            k += d.size
    stacked = ca.vertcat(*args) if args else ca.SX(0, 1)
    tau = graph.output(state, stacked)[:6]
    state_dot = graph.step(state, stacked) if graph.step is not None else ca.SX(0, 1)
    return ca.Function(f"producer_{index}", [command, state, nu_r, density], [tau, state_dot],
                       ["command", "state", "nu_r", "water_density"], ["tau", "state_dot"])


class ForceProducerSet(ForceProducer):
    COMPONENTS = {"producers": "List[more_dynamics::ForceProducer]"}
    PARAMETERS = [
        ParameterDescription("command_count", 1),  # one set command (convention)
        ParameterDescription("command_map", [[1.0]]),  # one child command per set command (convention)
    ]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        children = context.get_component("producers")
        if not children:
            raise ValueError("ForceProducerSet: slot 'producers' is empty")
        producers = [producer_function(child, i) for i, child in enumerate(children)]
        count = context.get_parameter("command_count")
        declared = force_producer_set_parameters(producers, count)
        self._model = freeze(force_producer_set_casadi(producers, count), declared,
                             {"command_map": context.get_parameter("command_map")})

    def graph(self) -> ForceProducer.CasadyPayload:
        if self._model is None:
            raise RuntimeError("ForceProducerSet must be initialized before graph()")
        io = PayloadBuilder(ForceProducer.CasadyPayload())
        n = self._model.size1_in("state")
        state = io.state("producer_state", n, [0.0] * n, "stacked states of the producers, list order") \
            if n else ca.SX.sym("state", 0)
        out = io.call(self._model, rename={"command": "force_producer_command"}, given={"state": state})
        io.output("generated_force", out["tau"], "summed wrench about the CO, BODY, N and N m")
        if n:
            io.state_dot(out["state_dot"])
        return io.payload()
