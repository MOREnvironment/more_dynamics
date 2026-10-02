from typing import List
import casadi as ca
from rpp_plugin_types.more_dynamics import VehicleModel3D
from rpp_schema.more_dynamics.IODescription import IODescription
from rpp_schema.more_dynamics.StateDescription import StateDescription
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_common.casadi_graph import RppCasadiGraph, graph_to_bytes
from more_dynamics.models.hull_vessel import (
    preprocess_hull_mass_properties,
    vessel_model_casadi,
)

class HullVessel(VehicleModel3D):

    COMPONENTS = {
        "sensors" : "List[more_sensors::Sensor]",
        "actuators" : "List[more_dynamics::ForceProducer]",
        "hydrostatics": "more_dynamics::HydrostaticsModel",
        "hydrodynamics": "more_dynamics::HydrodynamicsModel",
    }

    PARAMETERS = [
        ParameterDescription("mass", 801.01905),
        ParameterDescription(
            "inertia",
            [453.5820434315624, 1353.7221945000001, 1353.7221945000001],
        ),
        ParameterDescription(
            "added_mass",
            [
                67.85979986031616,
                1201.528575,
                801.01905,
                90.81653606756248,
                1083.378265125,
                2301.32773065,
            ],
        ),
        ParameterDescription("center_of_gravity", [0.0, 0.0, 0.025]),
    ]

    def __init__(self):
        self.actuators = []
        self.sensors = []
        self._model = None
        self._mass_properties = None
        self._actuator_graphs = []
        self._sensor_graphs = []
        self._hydrostatics_graph = None
        self._hydrodynamics_graph = None
        self.current_input_idx = 0
        self.current_state_idx = 0
        self.current_output_idx = 0

    def initialize(self, context: ComponentContext):
        hydrostatics = context.get_component("hydrostatics")
        self._hydrostatics_graph = RppCasadiGraph(hydrostatics.graph())
        if self._hydrostatics_graph.num_inputs != 6:
            raise ValueError("hydrostatics graph must accept a 6-value pose")
        if self._hydrostatics_graph.num_outputs < 6:
            raise ValueError(
                "hydrostatics graph must return a 6-value generalized force"
            )
        if self._hydrostatics_graph.num_states != 0:
            raise ValueError("hydrostatics graph must not declare states")
        if self._hydrostatics_graph.step is not None:
            raise ValueError("hydrostatics graph must not define dynamics")
        hydrodynamics = context.get_component("hydrodynamics")
        self._hydrodynamics_graph = RppCasadiGraph(hydrodynamics.graph())
        if self._hydrodynamics_graph.num_inputs != 6:
            raise ValueError("hydrodynamics graph must accept a 6-value velocity")
        if self._hydrodynamics_graph.num_outputs < 6:
            raise ValueError(
                "hydrodynamics graph must return a 6-value generalized force"
            )
        if self._hydrodynamics_graph.num_states != 0:
            raise ValueError("hydrodynamics graph must not declare states")
        if self._hydrodynamics_graph.step is not None:
            raise ValueError("hydrodynamics graph must not define dynamics")
        for actuator in context.get_component("actuators"):
            self.actuators.append(actuator)
            self._actuator_graphs.append(self._get_actuator_graph(actuator))
        for sensor in context.get_component("sensors"):
            self.sensors.append(sensor)
            self._sensor_graphs.append(self._get_sensor_graph(sensor))

        self._mass_properties = preprocess_hull_mass_properties(
            mass=context.get_parameter("mass"),
            inertia=context.get_parameter("inertia"),
            center_of_gravity=context.get_parameter("center_of_gravity"),
            added_mass=context.get_parameter("added_mass"),
        )
        self._model = vessel_model_casadi(
            mass_properties=self._mass_properties,
        )

    def step(self, state: VehicleModel3D.Odometry3D,
            command: List[VehicleModel3D.Command], t: float, dt: float, **kwargs) -> VehicleModel3D.Odometry3D:
        forces = []
        for actuator in self.actuators:
            forces.append(actuator.getForce(command))
        return state


    def graph(self):
        payload = VehicleModel3D.CasadyPayload()
        payload = self._graph_io(payload)

        sum_inputs = sum(actuator.num_inputs for actuator in self._actuator_graphs)
        sum_states = sum(actuator.num_states for actuator in self._actuator_graphs) + 12

        inputs = ca.SX.sym("input", sum_inputs)
        states = ca.SX.sym("state", sum_states)
        force = ca.SX.zeros(6)

        actuators_dot = []
        for actuator in self._actuator_graphs:
            slice_inputs = inputs[actuator.slice_input()]
            slice_state = states[actuator.slice_state()]
            force += actuator.output(slice_state, slice_inputs)
            if actuator.step is not None:
                actuators_dot.append(actuator.step(slice_state, slice_inputs))

        pose = states[-12:-6]
        hydrostatics_output = self._hydrostatics_graph.output(
            ca.SX.zeros(0, 1),
            pose,
        )
        force += hydrostatics_output[:6]

        velocity = states[-6:]
        hydrodynamics_output = self._hydrodynamics_graph.output(
            ca.SX.zeros(0, 1),
            velocity,
        )
        force += hydrodynamics_output[:6]

        state_dot = self._model(states[-12:], force)

        dyn_function = ca.Function("dynamics", [states, inputs], [ca.vertcat(*actuators_dot, state_dot)])
        vessel_state = states[-12:]
        sensor_outputs = [
            sensor.output(ca.SX.zeros(0, 1), vessel_state)
            for sensor in self._sensor_graphs
        ]
        output_function = ca.Function(
            "output",
            [states, inputs],
            [ca.vertcat(vessel_state, *sensor_outputs)],
        )

        payload.dynamics = graph_to_bytes(dyn_function)
        payload.output = graph_to_bytes(output_function)

        return payload

    def getInputDescription(self):
        payload = []
        for actuator in self._actuator_graphs:
            for act_input in actuator.payload.inputDescription:
                payload.append(act_input)
        return payload

    def _graph_io(self, payload : VehicleModel3D.CasadyPayload):
        for actuator in self._actuator_graphs:
            for act_input in actuator.payload.inputDescription:
                payload.inputDescription.append(act_input)
            for act_state in actuator.payload.stateDescription:
                payload.stateDescription.append(act_state)
            # for act_output in actuator.payload.outputDescription:
            #     payload.outputDescription.append(act_output)

        payload.outputDescription.append(
            IODescription(12, name="output")
        )
        for sensor in self._sensor_graphs:
            for sensor_output in sensor.payload.outputDescription:
                payload.outputDescription.append(sensor_output)
        payload.stateDescription.append(
            StateDescription(12, name="state")
        )

        return payload

    def _get_actuator_graph(self, actuator):
        payload = actuator.graph()
        graph =  RppCasadiGraph(payload,
            start_index_input=self.current_input_idx,
            start_index_state=self.current_state_idx,
            start_index_output=self.current_output_idx)
        if graph.step is None and graph.num_states != 0:
            raise ValueError(
                "force producer without dynamics must not declare states"
            )
        self.current_input_idx += graph.num_inputs
        self.current_state_idx += graph.num_states
        self.current_output_idx += graph.num_outputs
        return graph

    @staticmethod
    def _get_sensor_graph(sensor):
        graph = RppCasadiGraph(sensor.graph())
        if graph.num_inputs != 12:
            raise ValueError("sensor graph must accept the 12-value vessel state")
        if graph.num_states != 0:
            raise ValueError("sensor graph must not declare states")
        if graph.step is not None:
            raise ValueError("sensor graph must not define dynamics")
        return graph
