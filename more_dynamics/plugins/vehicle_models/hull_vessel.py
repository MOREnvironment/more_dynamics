from math import asin, atan2, cos, sin
from typing import List

import casadi as ca
import numpy as np
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
        ParameterDescription("include_coriolis", True),
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
        self._vessel_step = None
        self._include_coriolis = True
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

        self._include_coriolis = context.get_parameter("include_coriolis")
        if not isinstance(self._include_coriolis, bool):
            raise ValueError("include_coriolis must be a boolean")

        self._mass_properties = preprocess_hull_mass_properties(
            mass=context.get_parameter("mass"),
            inertia=context.get_parameter("inertia"),
            center_of_gravity=context.get_parameter("center_of_gravity"),
            added_mass=context.get_parameter("added_mass"),
        )
        resolved_inertia = (
            self._mass_properties.inertia_at_center_of_gravity.diagonal()
            .tolist()
        )
        resolved_added_mass = (
            self._mass_properties.added_mass_matrix.diagonal().tolist()
        )
        effective_surge_mass = self._mass_properties.total_mass_matrix[0, 0]
        # resolved_parameters = (
        #     "HullVessel resolved parameters: "
        #     f"mass={self._mass_properties.mass}, "
        #     f"inertia={resolved_inertia}, "
        #     f"center_of_gravity={context.get_parameter('center_of_gravity')}, "
        #     f"added_mass={resolved_added_mass}, "
        #     f"effective_surge_mass={effective_surge_mass}, "
        #     f"include_coriolis={self._include_coriolis}"
        # )
        # context.get_logger().info(resolved_parameters)
        self._model = vessel_model_casadi(
            mass_properties=self._mass_properties,
            include_coriolis=self._include_coriolis,
        )
        self._vessel_step = self._create_vessel_step()

    def step(self, state: VehicleModel3D.Odometry3D,
            command: List[VehicleModel3D.Command], t: float, dt: float, **kwargs) -> VehicleModel3D.Odometry3D:
        """Advance the vessel pose and twist by dt.

        The command values are read in order as the actuator inputs of
        graph(). Actuators are taken at the steady state of their command,
        because an Odometry3D carries no actuator state.
        """
        if self._vessel_step is None:
            raise RuntimeError("HullVessel must be initialized before step()")

        values = [float(value) for item in command for value in item.data]
        wrench = np.zeros(6)
        for actuator, graph in zip(self.actuators, self._actuator_graphs):
            inputs = values[graph.slice_input()]
            actuator_command = VehicleModel3D.Command()
            actuator_command.data.extend(
                inputs + [0.0] * (graph.num_inputs - len(inputs))
            )
            actuator_wrench = actuator.step(state, actuator_command, t, dt)
            wrench += [
                actuator_wrench.force.x,
                actuator_wrench.force.y,
                actuator_wrench.force.z,
                actuator_wrench.torque.x,
                actuator_wrench.torque.y,
                actuator_wrench.torque.z,
            ]

        next_state = np.asarray(
            self._vessel_step(self._state_from_odometry(state), wrench, dt)
        ).reshape(-1)
        return self._odometry_from_state(next_state)

    def _create_vessel_step(self) -> ca.Function:
        """Build one RK4 step of the vessel under a constant actuator wrench."""
        state = ca.SX.sym("state", 12)
        actuator_wrench = ca.SX.sym("actuator_wrench", 6)
        delta_t = ca.SX.sym("delta_t")

        def derivative(vessel_state):
            force = (
                actuator_wrench
                + self._hydrostatics_graph.output(
                    ca.SX.zeros(0, 1), vessel_state[:6]
                )[:6]
                + self._hydrodynamics_graph.output(
                    ca.SX.zeros(0, 1), vessel_state[6:]
                )[:6]
            )
            return self._model(vessel_state, force)

        k1 = derivative(state)
        k2 = derivative(state + delta_t * k1 / 2.0)
        k3 = derivative(state + delta_t * k2 / 2.0)
        k4 = derivative(state + delta_t * k3)
        next_state = state + delta_t * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
        return ca.Function(
            "vessel_step", [state, actuator_wrench, delta_t], [next_state]
        )

    @staticmethod
    def _state_from_odometry(odometry: VehicleModel3D.Odometry3D) -> np.ndarray:
        """Return position, roll-pitch-yaw, and body twist of an odometry."""
        position = odometry.pose.position
        orientation = odometry.pose.orientation
        norm = (
            orientation.x**2 + orientation.y**2
            + orientation.z**2 + orientation.w**2
        ) ** 0.5
        if norm < 1e-12:
            raise ValueError("odometry orientation must be a quaternion")
        x, y, z, w = (
            orientation.x / norm,
            orientation.y / norm,
            orientation.z / norm,
            orientation.w / norm,
        )
        linear = odometry.twist.linear
        angular = odometry.twist.angular
        return np.array(
            [
                position.x, position.y, position.z,
                atan2(2.0 * (w * x + y * z), 1.0 - 2.0 * (x * x + y * y)),
                asin(min(max(2.0 * (w * y - z * x), -1.0), 1.0)),
                atan2(2.0 * (w * z + x * y), 1.0 - 2.0 * (y * y + z * z)),
                linear.x, linear.y, linear.z,
                angular.x, angular.y, angular.z,
            ]
        )

    @staticmethod
    def _odometry_from_state(state: np.ndarray) -> VehicleModel3D.Odometry3D:
        """Build an odometry from position, roll-pitch-yaw, and body twist."""
        half_roll, half_pitch, half_yaw = (float(angle) / 2.0 for angle in state[3:6])
        odometry = VehicleModel3D.Odometry3D()
        position = odometry.pose.position
        position.x, position.y, position.z = (float(v) for v in state[0:3])
        orientation = odometry.pose.orientation
        orientation.x = (
            sin(half_roll) * cos(half_pitch) * cos(half_yaw)
            - cos(half_roll) * sin(half_pitch) * sin(half_yaw)
        )
        orientation.y = (
            cos(half_roll) * sin(half_pitch) * cos(half_yaw)
            + sin(half_roll) * cos(half_pitch) * sin(half_yaw)
        )
        orientation.z = (
            cos(half_roll) * cos(half_pitch) * sin(half_yaw)
            - sin(half_roll) * sin(half_pitch) * cos(half_yaw)
        )
        orientation.w = (
            cos(half_roll) * cos(half_pitch) * cos(half_yaw)
            + sin(half_roll) * sin(half_pitch) * sin(half_yaw)
        )
        linear = odometry.twist.linear
        linear.x, linear.y, linear.z = (float(v) for v in state[6:9])
        angular = odometry.twist.angular
        angular.x, angular.y, angular.z = (float(v) for v in state[9:12])
        return odometry


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
        velocity = states[-6:]
        hydrostatics_output = self._hydrostatics_graph.output(
            ca.SX.zeros(0, 1),
            pose,
        )
        force += hydrostatics_output[:6]

        hydrodynamics_output = self._hydrodynamics_graph.output(
            ca.SX.zeros(0, 1),
            velocity,
        )
        force += hydrodynamics_output[:6]

        state_dot = self._model(states[-12:], force)

        dyn_function = ca.Function("dynamics", [states, inputs], [ca.vertcat(*actuators_dot, state_dot)])
        vessel_state = states[-12:]
        vessel_acceleration = state_dot[6:12]
        sensor_outputs = [
            sensor.output(
                ca.SX.zeros(0, 1),
                ca.vertcat(vessel_state, vessel_acceleration)
                if sensor.num_inputs == 18
                else vessel_state,
            )
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
        if graph.num_inputs not in (12, 18):
            raise ValueError(
                "sensor graph must accept the 12-value vessel state, "
                "optionally followed by the 6-value body acceleration"
            )
        if graph.num_states != 0:
            raise ValueError("sensor graph must not declare states")
        if graph.step is not None:
            raise ValueError("sensor graph must not define dynamics")
        return graph
