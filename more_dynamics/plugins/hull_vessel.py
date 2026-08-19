from typing import List
import casadi as ca
from rpp_plugin_types.more_dynamics import VehicleModel3D
from rpp_schema.more_dynamics.IODescription import IODescription
from rpp_schema.more_dynamics.StateDescription import StateDescription
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from .casadi_helpers import vessel_model_casadi, graph_to_bytes, RppCasadiGraph

class HullVessel(VehicleModel3D):

    COMPONENTS = {
        "sensors" : "List[more_sensors::Sensor]",
        "actuators" : "List[more_dynamics::ForceProducer]"
    }

    PARAMETERS = [
        ParameterDescription("mass", 1000.0),
        ParameterDescription("inertia", [100.0, 100.0, 100.0]),
        ParameterDescription("drag_coefficient", 0.5),
        ParameterDescription("added_mass", [50.0, 50.0, 50.0]),
        ParameterDescription("center_of_gravity", [0.0, 0.0, 0.0]),
        ParameterDescription("center_of_buoyancy", [0.0, 0.0, 0.0]),
        ParameterDescription("water_density", 1000.0),
        ParameterDescription("gravity", 9.81),
    ]

    def __init__(self):
        self.actuators = []
        self.sensors = []
        self._model = None
        self._actuator_graphs = []
        self.current_input_idx = 0
        self.current_state_idx = 0
        self.current_output_idx = 0


    def initialize(self, context: ComponentContext):
        for actuator in context.get_component("actuators"):
            self.actuators.append(actuator)
            self._actuator_graphs.append(self._get_actuator_graph(actuator))
        for sensor in context.get_component("sensors"):
            self.sensors.append(sensor)

        self._model = vessel_model_casadi(
            mass=context.get_parameter("mass"),
            inertia=context.get_parameter("inertia"),
            drag_coefficient=context.get_parameter("drag_coefficient"),
            added_mass=context.get_parameter("added_mass"),
            center_of_buoyancy=context.get_parameter("center_of_buoyancy"),
            center_of_gravity=context.get_parameter("center_of_gravity"),
            gravity=context.get_parameter("gravity"),
            water_density=context.get_parameter("water_density")
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
            actuators_dot.append(actuator.step(slice_state, slice_inputs))

        state_dot = self._model(states[-12:], force)

        dyn_function = ca.Function("dynamics", [states, inputs], [ca.vertcat(*actuators_dot, state_dot)])
        output_function = ca.Function("output", [states, inputs], [states[-12:]])

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
        self.current_input_idx += graph.num_inputs
        self.current_state_idx += graph.num_states
        self.current_output_idx += graph.num_outputs
        return graph