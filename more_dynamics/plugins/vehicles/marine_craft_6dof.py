"""A 6-DOF marine craft assembled from typed child slots, in the form of
``HullVessel`` (``hull_vessel.py``): children are taken by slot, each
child's ``graph()`` payload is wrapped in ``RppCasadiGraph``, children with
states are stacked before the 12 vessel states, the signed forces are
summed, and the equation of motion is the model layer's
(``models/vehicles/marine_craft_6dof``).

What it adds to ``HullVessel``: a child's input is connected **by the name
and size in its ``inputDescription``** — to a vehicle quantity (``pose``,
``velocity``, ``relative_velocity``, ``mass_matrix``) or to an output of an
earlier single-slot child of the same name; a missing or mis-sized one is
refused, naming the slot and the child's plugin. The slot order is the
physics (water -> current -> hull form -> rigid body -> added mass -> its
Coriolis -> restoring -> hull loads -> force producers), so a child can only
consume what is already known.

* single slots: every output becomes a vehicle quantity under its name (one
  name, one producer); required outputs per slot are checked.
* ``restoring`` / ``hydrodynamic_loads`` / ``force_producers``: the first
  output entry is a 6-vector signed generalized force the vehicle adds
  (``vehicle_model.capnp`` 11-21; ``hull_vessel.py`` 122, 131, 138).
* list members' further outputs are visible as ``<slot>.<k>.<name>``
  (``diagnostic_outputs``), never consumed.
* a force producer's inputs that are no vehicle quantity are commands: the
  vehicle's inputs, in list order (``hull_vessel.py`` 159-171).
* ``open_inputs`` (``"<slot>.<name>"``): a part input left open on purpose
  (the part's ``open_parameters``) becomes a vehicle input of that name.
* no ``sensors`` slot: the output is the 12 vessel states (plus the declared
  ``diagnostic_outputs``); sensors attach outside the vehicle.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics
    and Motion Control. Wiley. Eq. 6.48, p. 120 (M = M_RB + M_A).
[MSS] Fossen, T. I. MSS, MIT, @ cc07579: LIBRARY/numericalMethods/rk4.m
    26-38 (the step's sub-stepping); CRAFT/AUV/models/SIMremus100.m 60
    (integration step h = 0.05 s).
[HullVessel] Mandic, L. more_dynamics ``HullVessel`` (hull_vessel.py): the
    composite form (children by slot, graphs combined with
    ``RppCasadiGraph``, output layout) this class follows.

Author:    Enio Krizman
Date:      2026-10-08
"""
import math
from typing import List

import casadi as ca
import numpy as np

from rpp_plugin_types.more_dynamics import VehicleModel3D
from rpp_schema.more_dynamics.IODescription import IODescription
from rpp_schema.more_dynamics.StateDescription import StateDescription
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_common.casadi_graph import RppCasadiGraph, graph_to_bytes
from more_transformations.more_casadi_transformations import EulerQuaternionTransforms
from more_dynamics.models.vehicles.marine_craft_6dof import marine_craft_6dof_casadi

# slot -> is a list; walked in this order (the physics)
SLOTS = (("site", False), ("current", False), ("hull_form", False), ("rigid_body", False),
         ("added_mass", False), ("added_mass_coriolis", False), ("restoring", False),
         ("hydrodynamic_loads", True), ("force_producers", True))
REQUIRED = {"site": {"water_density": 1, "gravity": 1},
            "current": {"current_velocity": 6, "current_acceleration": 6},
            "hull_form": {"length": 1},
            "rigid_body": {"rigid_body_mass_matrix": 36, "rigid_body_coriolis_matrix": 36, "mass": 1},
            "added_mass": {"added_mass_matrix": 36},
            "added_mass_coriolis": {"added_mass_coriolis_matrix": 36}}
SIGNED_FORCE = ("restoring", "hydrodynamic_loads", "force_producers")


class CompositionError(ValueError):
    """A tree that cannot be assembled; the message names the slot and the child."""


class MarineCraft6DOF(VehicleModel3D):

    COMPONENTS = {
        "site": "more_dynamics::SiteModel",
        "current": "more_dynamics::CurrentModel",
        "hull_form": "more_dynamics::HullForm",
        "rigid_body": "more_dynamics::RigidBodyModel",
        "added_mass": "more_dynamics::AddedMassModel",
        "added_mass_coriolis": "more_dynamics::AddedMassCoriolisModel",
        "restoring": "more_dynamics::HydrostaticsModel",
        "hydrodynamic_loads": "List[more_dynamics::HydrodynamicsModel]",
        "force_producers": "List[more_dynamics::ForceProducer]",
    }

    PARAMETERS = [
        ParameterDescription("integration_max_step", 0.05),  # s, SIMremus100.m 60 (h = 0.05)
        ParameterDescription("open_inputs", []),
        ParameterDescription("diagnostic_outputs", []),
    ]

    def __init__(self):
        self._children = []          # (label, plugin class name, RppCasadiGraph)
        self._input_descriptions = []
        self._state_descriptions = []
        self._output_descriptions = []
        self._dynamics = None
        self._output = None
        self._rk4 = None
        self._max_step = None
        self._child_states = 0
        self.wiring = {}             # "<slot label>.<input>" -> where it came from (read by tests and reports)

    # -- build ------------------------------------------------------------------------------
    def initialize(self, context: ComponentContext):
        state_index = 0
        for slot, is_list in SLOTS:
            children = context.get_component(slot)
            if not is_list:
                if children == []:
                    raise CompositionError(f"slot {slot!r} is empty: every single slot needs a part")
                children = [children]
            for k, child in enumerate(children):
                label = f"{slot}.{k}" if is_list else slot
                graph = RppCasadiGraph(child.graph(), start_index_state=state_index)
                if graph.num_states and slot != "force_producers":
                    raise CompositionError(f"slot {label} ({type(child).__name__}): only force producers may "
                                           "declare states")
                if graph.step is None and graph.num_states:
                    raise CompositionError(f"slot {label} ({type(child).__name__}): states without dynamics")
                state_index += graph.num_states
                self._children.append((slot, label, type(child).__name__, graph))
        self._child_states = state_index
        self._max_step = float(context.get_parameter("integration_max_step"))
        if not self._max_step > 0.0:
            raise CompositionError("integration_max_step must be positive")
        self._assemble(list(context.get_parameter("open_inputs") or []),
                       list(context.get_parameter("diagnostic_outputs") or []))

    def _assemble(self, open_inputs, diagnostic_outputs):
        states = ca.SX.sym("state", self._child_states + 12)
        vessel = states[-12:]
        known = {"pose": (vessel[0:6], "vehicle"), "velocity": (vessel[6:12], "vehicle")}
        visible = {}
        inputs, input_descriptions, state_dots = [], [], []
        force = ca.SX.zeros(6)
        unused_open = set(open_inputs)

        for slot, label, plugin, graph in self._children:
            args = []
            for d in graph.payload.inputDescription:
                where = f"slot {label} ({plugin}): input {d.name!r} ({d.size})"
                if d.name in known:
                    value, origin = known[d.name]
                    if value.numel() != d.size:
                        raise CompositionError(f"{where} but {origin} gives {value.numel()} values")
                    self.wiring[f"{label}.{d.name}"] = origin
                elif f"{label}.{d.name}" in open_inputs:
                    value = ca.SX.sym(f"{label}.{d.name}", d.size)
                    inputs.append(value)
                    input_descriptions.append(IODescription(d.size, name=f"{label}.{d.name}",
                                                            description=f"open parameter of {plugin}"))
                    unused_open.discard(f"{label}.{d.name}")
                    self.wiring[f"{label}.{d.name}"] = "open input"
                elif slot == "force_producers":
                    value = ca.SX.sym(f"{label}.{d.name}", d.size)
                    inputs.append(value)
                    input_descriptions.append(IODescription(d.size, list(d.min), list(d.max), d.name, d.description))
                    self.wiring[f"{label}.{d.name}"] = "command"
                else:
                    listed = [n for n in open_inputs if n.startswith(f"{label}.")]
                    hint = f"; open_inputs for this slot: {listed}" if listed else f" (known: {sorted(known)})"
                    raise CompositionError(f"{where}: no vehicle quantity and no earlier slot produces it{hint}")
                args.append(value)
            stacked = ca.vertcat(*args) if args else ca.SX(0, 1)
            child_state = states[graph.slice_state()]
            outputs = _split(graph.payload.outputDescription, graph.output(child_state, stacked), label, plugin)
            if graph.step is not None:
                state_dots.append(graph.step(child_state, stacked))
            for name, size in REQUIRED.get(slot, {}).items():
                if name not in outputs or outputs[name].numel() != size:
                    raise CompositionError(f"slot {label} ({plugin}): a {slot} part must produce {name!r} ({size})")
            if slot in SIGNED_FORCE:
                first = graph.payload.outputDescription[0] if graph.payload.outputDescription else None
                if first is None or first.size != 6:
                    raise CompositionError(f"slot {label} ({plugin}): the first output must be the 6-value "
                                           "signed generalized force")
                force += outputs[first.name]
            if slot in ("hydrodynamic_loads", "force_producers"):
                visible.update({f"{label}.{n}": v for n, v in outputs.items()})
            else:
                for name, value in outputs.items():
                    if name in known:
                        raise CompositionError(f"slot {label} ({plugin}): {name!r} is produced already by "
                                               f"{known[name][1]} (one quantity, one value)")
                    known[name] = (value, f"slot {label} ({plugin})")
            if slot == "current":
                known["relative_velocity"] = (vessel[6:12] - known["current_velocity"][0], "vehicle")  # nu_r = nu - nu_c
            if slot == "added_mass":
                known["mass_matrix"] = (known["rigid_body_mass_matrix"][0] + known["added_mass_matrix"][0],
                                        "vehicle")  # M = M_RB + M_A (Fossen 2011, eq. 6.48, p. 120)
        if unused_open:
            raise CompositionError(f"open_inputs {sorted(unused_open)}: no such open part input")

        mass_matrix = ca.reshape(known["mass_matrix"][0], 6, 6)
        coriolis = ca.reshape(known["rigid_body_coriolis_matrix"][0] + known["added_mass_coriolis_matrix"][0], 6, 6)
        vessel_dot = marine_craft_6dof_casadi()(vessel, known["current_velocity"][0],
                                                known["current_acceleration"][0], mass_matrix, coriolis, force)
        u = ca.vertcat(*inputs) if inputs else ca.SX(0, 1)
        self._dynamics = ca.Function("dynamics", [states, u], [ca.vertcat(*state_dots, vessel_dot)],
                                     ["state", "input"], ["state_dot"])

        outputs = [vessel]
        self._output_descriptions = [IODescription(12, name="output")]  # (hull_vessel.py 175-177)
        everything = {**{n: v for n, (v, _) in known.items()}, **visible}
        for name in diagnostic_outputs:
            if name not in everything:
                raise CompositionError(f"diagnostic output {name!r}: no such quantity (known: {sorted(everything)})")
            outputs.append(everything[name])
            self._output_descriptions.append(IODescription(everything[name].numel(), name=name))
        self._output = ca.Function("output", [states, u], [ca.vertcat(*outputs)], ["state", "input"], ["output"])

        self._input_descriptions = input_descriptions
        self._state_descriptions = [d for _, _, _, g in self._children for d in g.payload.stateDescription] + \
            [StateDescription(12, ic=[0.0] * 12, name="state")]  # (hull_vessel.py 181-183)
        h = ca.SX.sym("h")
        k1 = self._dynamics(states, u)
        k2 = self._dynamics(states + h / 2 * k1, u)
        k3 = self._dynamics(states + h / 2 * k2, u)
        k4 = self._dynamics(states + h * k3, u)
        self._rk4 = ca.Function("rk4", [states, u, h], [states + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)])  # (rk4.m 26-38)

    # -- the plugin type's methods ---------------------------------------------------------------
    def graph(self):
        if self._dynamics is None:
            raise RuntimeError("MarineCraft6DOF must be initialized before graph()")
        payload = VehicleModel3D.CasadyPayload()
        payload.inputDescription.extend(self._input_descriptions)
        payload.outputDescription.extend(self._output_descriptions)
        payload.stateDescription.extend(self._state_descriptions)
        payload.dynamics = graph_to_bytes(self._dynamics)
        payload.output = graph_to_bytes(self._output)
        return payload

    def getInputDescriptions(self):
        return list(self._input_descriptions)

    def step_vector(self, x, u, dt):
        """``ceil(dt / integration_max_step)`` RK4 sub-steps of the vehicle's own
        dynamics, zero-order hold on ``u`` (rk4.m 26-38)."""
        x = np.asarray(x, float).ravel()
        if dt == 0:
            return x
        n = math.ceil(dt / self._max_step - 1e-12)
        for _ in range(n):
            x = np.asarray(self._rk4(x, u, dt / n)).ravel()
        return x

    def step(self, state: VehicleModel3D.Odometry3D, command: List[VehicleModel3D.Command], t: float, dt: float,
             **kwargs) -> VehicleModel3D.Odometry3D:
        if self._child_states:
            raise NotImplementedError("step() with force-producer states waits on an owner decision "
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
