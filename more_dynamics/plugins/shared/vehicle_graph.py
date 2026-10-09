"""The graph of a marine-craft vehicle type, in the form of ``HullVessel``
(``hull_vessel.py``): the children's payloads are wrapped in
``RppCasadiGraph``, the actuators' states are stacked before the 12 vessel
states, the signed generalized forces of the hydrostatics, the hydrodynamics
and every actuator are summed, and the equation of motion is the model
layer's (``models/equation_of_motion``).

What it adds to ``HullVessel``: a child's input is connected **by the name
and size in its ``inputDescription``** to a quantity the vehicle type
computed from its own parameters (``pose``, ``velocity``,
``relative_velocity``, ``mass_matrix``, ``water_density`` ...), or to an
output of an earlier child (hydrostatics before hydrodynamics before the
actuators); a missing or mis-sized one is refused, naming the slot and the
child's plugin. An actuator input that is no vehicle quantity is a command:
the vehicle's inputs, in list order, after the current inputs the type asks
for. No ``sensors`` slot: the output is the 12 vessel states; sensors attach
outside the vehicle.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics
    and Motion Control. Wiley. Eq. 6.48, p. 120 (M = M_RB + M_A).
[MSS] Fossen, T. I. MSS, MIT, @ cc07579: LIBRARY/numericalMethods/rk4.m
    26-38 (the step's sub-stepping); CRAFT/AUV/models/SIMremus100.m 60
    (integration step h = 0.05 s).
[HullVessel] Mandic, L. more_dynamics ``HullVessel`` (hull_vessel.py 132-182):
    the composite form (children by slot, graphs combined with
    ``RppCasadiGraph``, output layout) this module follows.

Author:    Enio Krizman
Date:      2026-10-09
"""
import math
from dataclasses import dataclass, field

import casadi as ca
import numpy as np

from rpp_plugin_types.more_dynamics import VehicleModel3D
from rpp_schema.more_dynamics.IODescription import IODescription
from rpp_schema.more_dynamics.StateDescription import StateDescription
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_common.casadi_graph import RppCasadiGraph, graph_to_bytes
from more_transformations.more_casadi_transformations import EulerQuaternionTransforms
from more_dynamics.models.equation_of_motion import marine_craft_6dof_casadi
from more_dynamics.plugins.shared.payload_io import frozen_block, payload_name

INTEGRATION_STEP = ParameterDescription("integration_max_step", 0.05)  # s, SIMremus100.m 60 (h = 0.05)


class CompositionError(ValueError):
    """A composition that cannot be assembled; the message names the slot and the child."""


@dataclass
class Craft:
    """What a vehicle type computes from its own parameters, as expressions of the pose, the velocity and the
    current inputs: the inertia, the Coriolis matrix, the current, and every named quantity its children may read
    (payload names)."""
    mass_matrix: ca.SX
    coriolis_matrix: ca.SX
    current_velocity: ca.SX
    current_acceleration: ca.SX
    known: dict = field(default_factory=dict)


def call_model(function, known, label):
    """Call a model-layer function whose inputs are all vehicle quantities: each input is read from ``known`` under
    its payload name. Returns ``{model output name: expression}``."""
    args = {}
    for name in function.name_in():
        key = payload_name(name)
        if key not in known:
            raise CompositionError(f"{label}: needs {key!r}, which no quantity of this vehicle gives "
                                   f"(known: {sorted(known)})")
        value = known[key]
        if value.numel() != function.numel_in(name):
            raise CompositionError(f"{label}: {key!r} has {value.numel()} values, the model takes "
                                   f"{function.numel_in(name)}")
        args[name] = ca.reshape(value, *function.size_in(name))
    return function.call(args)


class Parts:
    """The model-layer functions a vehicle type calls in order, each frozen to the type's own parameters and
    read from / written to one dictionary of named quantities (payload names)."""

    def __init__(self, context: ComponentContext, plugin: str, known: dict):
        self.context, self.plugin, self.known = context, plugin, known

    def run(self, label, function, declared=()):
        """Freeze ``function`` to the parameters ``declared`` (their values from the vehicle's parameters), call it
        on the quantities known so far and add its outputs under their payload names. Returns the outputs by their
        model names."""
        frozen = frozen_block(self.context, function, declared)
        out = call_model(frozen, self.known, f"{self.plugin}: {label}")
        for name in frozen.name_out():
            key = payload_name(name)
            if key in self.known:
                raise CompositionError(f"{self.plugin}: {label} gives {key!r}, produced already (one quantity, "
                                       "one value)")
            self.known[key] = out[name]
        return out

    def select(self, parameter, forms):
        """The form named by the vehicle's parameter ``parameter``, refusing an unknown one."""
        form = self.context.get_parameter(parameter)
        if form not in forms:
            raise CompositionError(f"{self.plugin}: {parameter} must be one of {sorted(forms)}, got {form!r}")
        return forms[form]


def _split(descriptions, vector, label, plugin):
    out, k = {}, 0
    for d in descriptions:
        if d.name in out:
            raise CompositionError(f"slot {label} ({plugin}): output {d.name!r} declared twice")
        out[d.name] = vector[k:k + d.size]
        k += d.size
    if k != vector.numel():
        raise CompositionError(f"slot {label} ({plugin}): descriptions cover {k} outputs, the graph gives "
                               f"{vector.numel()}")
    return out


class VehicleAssembly:
    """The assembled vehicle: its payload pieces, a numeric RK4 step and ``signals`` for reading any named
    quantity back (tests, diagnostics); built once in ``initialize``."""

    def __init__(self, context: ComponentContext, plugin: str, compute, current_inputs=()):
        self.plugin = plugin
        hydrostatics = context.get_component("hydrostatics")
        hydrodynamics = context.get_component("hydrodynamics")
        for slot, child in (("hydrostatics", hydrostatics), ("hydrodynamics", hydrodynamics)):
            if child == []:
                raise CompositionError(f"{plugin}: slot {slot!r} is empty: a vehicle needs a {slot} model")
        actuators = list(context.get_component("actuators"))

        self._max_step = float(context.get_parameter("integration_max_step"))
        if not self._max_step > 0.0:
            raise CompositionError(f"{plugin}: integration_max_step must be positive")

        children, state_index = [], 0
        for slot, child in (("hydrostatics", hydrostatics), ("hydrodynamics", hydrodynamics)):
            graph = RppCasadiGraph(child.graph())
            if graph.num_states or graph.step is not None:
                raise CompositionError(f"slot {slot} ({type(child).__name__}): a {slot} model declares no states "
                                       "and no dynamics")
            children.append((slot, slot, type(child).__name__, graph))
        for k, child in enumerate(actuators):
            graph = RppCasadiGraph(child.graph(), start_index_state=state_index)
            if graph.step is None and graph.num_states:
                raise CompositionError(f"slot actuators.{k} ({type(child).__name__}): states without dynamics")
            state_index += graph.num_states
            children.append(("actuators", f"actuators.{k}", type(child).__name__, graph))
        self.child_states = state_index

        states = ca.SX.sym("state", self.child_states + 12)
        vessel = states[-12:]
        pose, velocity = vessel[0:6], vessel[6:12]
        inputs, input_descriptions = [], []
        current = {}
        for name, size, description in current_inputs:
            current[name] = ca.SX.sym(name, size)
            inputs.append(current[name])
            input_descriptions.append(IODescription(size, name=name, description=description))

        craft = compute(pose, velocity, current)
        known = {"pose": pose, "velocity": velocity, **craft.known}
        known.setdefault("relative_velocity", velocity - craft.current_velocity)  # nu_r = nu - nu_c
        known["mass_matrix"] = craft.mass_matrix  # M = M_RB + M_A (Fossen 2011, eq. 6.48, p. 120)
        origin = {name: f"{plugin}" for name in known}
        visible = {}
        force = ca.SX.zeros(6)
        state_dots = []

        for slot, label, child_plugin, graph in children:
            args = []
            for d in graph.payload.inputDescription:
                where = f"slot {label} ({child_plugin}): input {d.name!r} ({d.size})"
                if d.name in known:
                    if known[d.name].numel() != d.size:
                        raise CompositionError(f"{where} but {origin[d.name]} gives {known[d.name].numel()} values")
                    args.append(ca.vec(known[d.name]))
                elif slot == "actuators":
                    value = ca.SX.sym(f"{label}.{d.name}", d.size)
                    inputs.append(value)
                    input_descriptions.append(IODescription(d.size, list(d.min), list(d.max), d.name, d.description))
                    args.append(value)
                else:
                    raise CompositionError(f"{where}: no vehicle quantity and no earlier slot produces it "
                                           f"(known: {sorted(known)})")
            stacked = ca.vertcat(*args) if args else ca.SX(0, 1)
            child_state = states[graph.slice_state()] if slot == "actuators" else ca.SX.zeros(0, 1)
            outputs = _split(graph.payload.outputDescription, graph.output(child_state, stacked), label,
                             child_plugin)
            if slot == "actuators" and graph.step is not None:
                state_dots.append(graph.step(child_state, stacked))
            first = graph.payload.outputDescription[0] if graph.payload.outputDescription else None
            if first is None or first.size != 6:
                raise CompositionError(f"slot {label} ({child_plugin}): the first output must be the 6-value "
                                       "signed generalized force")
            force += outputs[first.name]
            for name, value in outputs.items():
                visible[f"{label}.{name}"] = value
                if slot != "actuators" and name != first.name:
                    if name in known:
                        raise CompositionError(f"slot {label} ({child_plugin}): {name!r} is produced already by "
                                               f"{origin[name]} (one quantity, one value)")
                    known[name] = value
                    origin[name] = f"slot {label} ({child_plugin})"

        vessel_dot = marine_craft_6dof_casadi()(vessel, craft.current_velocity, craft.current_acceleration,
                                                ca.reshape(craft.mass_matrix, 6, 6),
                                                ca.reshape(craft.coriolis_matrix, 6, 6), force)
        u = ca.vertcat(*inputs) if inputs else ca.SX(0, 1)
        self.dynamics = ca.Function("dynamics", [states, u], [ca.vertcat(*state_dots, vessel_dot)],
                                    ["state", "input"], ["state_dot"])
        self.output = ca.Function("output", [states, u], [vessel], ["state", "input"], ["output"])
        self.input_descriptions = input_descriptions
        self.output_descriptions = [IODescription(12, name="output")]
        self.state_descriptions = [d for _, _, _, g in children for d in g.payload.stateDescription] + \
            [StateDescription(12, ic=[0.0] * 12, name="state")]

        everything = {**{n: ca.vec(v) for n, v in known.items()}, **{n: ca.vec(v) for n, v in visible.items()}}
        self.signal_names = list(everything)
        self.signal_sizes = [v.numel() for v in everything.values()]
        self._signals = ca.Function("signals", [states, u], [ca.vertcat(*everything.values())])
        self.wiring = {n: origin.get(n, "child output") for n in everything}

        h = ca.SX.sym("h")
        k1 = self.dynamics(states, u)
        k2 = self.dynamics(states + h / 2 * k1, u)
        k3 = self.dynamics(states + h / 2 * k2, u)
        k4 = self.dynamics(states + h * k3, u)
        self._rk4 = ca.Function("rk4", [states, u, h], [states + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)])  # (rk4.m 26-38)

    def signals(self, x, u=()):
        """``{name: values}`` of every quantity of the vehicle and of its children's outputs (``<slot>.<name>``) at
        state ``x`` (children's states first, then the 12 vessel states) and input ``u``."""
        vector = np.asarray(self._signals(np.asarray(x, float).ravel(), np.asarray(u, float).ravel())).ravel()
        out, k = {}, 0
        for name, size in zip(self.signal_names, self.signal_sizes):
            out[name] = vector[k:k + size]
            k += size
        return out

    def payload(self):
        payload = VehicleModel3D.CasadyPayload()
        payload.inputDescription.extend(self.input_descriptions)
        payload.outputDescription.extend(self.output_descriptions)
        payload.stateDescription.extend(self.state_descriptions)
        payload.dynamics = graph_to_bytes(self.dynamics)
        payload.output = graph_to_bytes(self.output)
        return payload

    def step_vector(self, x, u, dt):
        """``ceil(dt / integration_max_step)`` RK4 sub-steps of the vehicle's own dynamics, zero-order hold on ``u``
        (rk4.m 26-38)."""
        x = np.asarray(x, float).ravel()
        if dt == 0:
            return x
        n = math.ceil(dt / self._max_step - 1e-12)
        for _ in range(n):
            x = np.asarray(self._rk4(x, u, dt / n)).ravel()
        return x

    def step(self, state, command, t, dt):
        if self.child_states:
            raise NotImplementedError("step() with actuator states waits on an owner decision "
                                      "(where those states live between calls is not decided)")
        q = state.pose.orientation
        euler = [float(a) for a in EulerQuaternionTransforms.quaternion_to_euler([q.w, q.x, q.y, q.z])]
        p, v, w = state.pose.position, state.twist.linear, state.twist.angular
        x = [p.x, p.y, p.z, *euler, v.x, v.y, v.z, w.x, w.y, w.z]
        u = [value for c in command for value in c.data]
        x = self.step_vector(x, u, dt)
        out = VehicleModel3D.Odometry3D()
        out.pose.position.x, out.pose.position.y, out.pose.position.z = (float(a) for a in x[0:3])
        quaternion = np.asarray(EulerQuaternionTransforms.euler_to_quaternion(*x[3:6])).ravel()
        out.pose.orientation.w, out.pose.orientation.x, out.pose.orientation.y, out.pose.orientation.z = \
            (float(a) for a in quaternion)
        out.twist.linear.x, out.twist.linear.y, out.twist.linear.z = (float(a) for a in x[6:9])
        out.twist.angular.x, out.twist.angular.y, out.twist.angular.z = (float(a) for a in x[9:12])
        return out
