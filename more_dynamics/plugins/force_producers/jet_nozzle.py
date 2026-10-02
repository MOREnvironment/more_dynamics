import casadi as ca

from rpp_plugin_types.more_dynamics import ForceProducer
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription
from rpp_schema.more_dynamics.IODescription import IODescription
from rpp_schema.more_dynamics.StateDescription import StateDescription

from more_common.casadi_graph import graph_to_bytes

class JetNozzle(ForceProducer):

    PARAMETERS = [
        ParameterDescription("location", [0.0, 0.0, 0.0]),
        ParameterDescription("trim_angle", 0.0),
        ParameterDescription("yaw_bias", 0.0),
        ParameterDescription("max_thrust", 100.0),
        ParameterDescription("thrust_coefficient", 1.0),
        ParameterDescription("thrust_rise_time", 0.5),
        ParameterDescription("thrust_fall_time", 0.1),
        ParameterDescription("nozzle_velocity_radpersec", 30.0),
        ParameterDescription("min_angle", -1),
        ParameterDescription("max_angle", 1),
    ]

    def __init__(self):
        self.location = 0
        self.max_thrust = 0
        self.yaw_bias = 0.0
        self.trim_angle = 0.0
        self.thrust_coefficient = 0
        self.thrust_rise_time = 0.0
        self.thrust_fall_time = 0.0
        self.nozzle_velocity_radpersec = 0.0
        self.min_angle = 0.0
        self.max_angle = 0.0

    def initialize(self, context: ComponentContext):
        self.location = context.get_parameter("location")
        self.yaw_bias = context.get_parameter("yaw_bias")
        self.trim_angle = context.get_parameter("trim_angle")
        self.max_thrust = context.get_parameter("max_thrust")
        self.thrust_coefficient = context.get_parameter("thrust_coefficient")
        self.thrust_rise_time = context.get_parameter("thrust_rise_time")
        self.thrust_fall_time = context.get_parameter("thrust_fall_time")
        self.nozzle_velocity_radpersec = context.get_parameter("nozzle_velocity_radpersec")
        self.min_angle = context.get_parameter("min_angle")
        self.max_angle = context.get_parameter("max_angle")

    def graph(self) -> ForceProducer.CasadyPayload:
        payload = ForceProducer.CasadyPayload()
        payload = self._graph_io(payload)

        location = ca.DM(self.location)
        states = ca.SX.sym("state", 2)
        inputs = ca.SX.sym("input", 2)

        thrust = ca.fmax(ca.fmin(states[0], 1.0), -1.0)
        nozzle_angle = states[1]

        nozzle_angle = ca.fmax(ca.fmin(nozzle_angle, self.max_angle), self.min_angle)
        thrust_ref = ca.fmax(ca.fmin(inputs[0], 1.0), -1.0)
        normalized_nozzle_angle_ref = ca.fmax(ca.fmin(inputs[1], 1.0), -1.0)
        nozzle_angle_ref = ca.if_else(
            normalized_nozzle_angle_ref >= 0.0,
            normalized_nozzle_angle_ref * self.max_angle,
            -normalized_nozzle_angle_ref * self.min_angle,
        )
        thrust_dot = ca.if_else(
            thrust_ref > thrust,
            (thrust_ref - thrust) / self.thrust_rise_time,
            (thrust_ref - thrust) / self.thrust_fall_time
        )

        nozzle_angle_dot = ca.if_else(
            nozzle_angle_ref > nozzle_angle,
            ca.fmin(nozzle_angle_ref - nozzle_angle, self.nozzle_velocity_radpersec),
            ca.fmax(nozzle_angle_ref - nozzle_angle, -self.nozzle_velocity_radpersec)
        )

        x_dot = ca.vertcat(thrust_dot, nozzle_angle_dot)

        dynamics_fn = ca.Function(
            "dynamics", [states, inputs], [x_dot], ["state", "input"], ["state_dot"]
        )

        payload.dynamics = graph_to_bytes(dynamics_fn)

        generated_thrust = (
            self.max_thrust * self.thrust_coefficient * thrust
        ) * ca.vertcat(
            ca.cos(nozzle_angle + self.yaw_bias) * ca.cos(self.trim_angle),
            ca.sin(nozzle_angle + self.yaw_bias),
            -ca.sin(self.trim_angle),
        )
        generated_moment = ca.cross(location, generated_thrust)

        payload.output = graph_to_bytes(ca.Function(
            "output",
            [states, inputs],
            [ca.vertcat(generated_thrust, generated_moment)],
            ["state", "input"],
            ["output"],
        ))


        return payload


    def _graph_io(self, payload: ForceProducer.CasadyPayload):
        payload.inputDescription.append(
            IODescription(
                name="desired_thrust",
                size=1,
                min=[-1.0],
                max=[1.0],
            )
        )
        payload.inputDescription.append(
            IODescription(
                name="desired_normalized_nozzle_angle",
                size=1,
                min=[-1.0],
                max=[1.0],
            )
        )
        payload.outputDescription.append(
            IODescription(
                name="generated_thrust",
                size=6,
                min=[],
                max=[],
            )
        )

        payload.stateDescription.append(
            StateDescription(
                name="thrust",
                size=1,
                min=[-1.0],
                max=[1.0],
                ic=[0.0]
            )
        )

        payload.stateDescription.append(
            StateDescription(
                name="nozzle_angle",
                size=1,
                min=[self.min_angle],
                max=[self.max_angle],
                ic=[0.0]
            )
        )

        return payload
