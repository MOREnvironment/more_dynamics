from typing import List
import casadi as ca
from rpp_schema.more_dynamics.CasadyPayload import CasadyPayload


class RppCasadiGraph:
    def __init__(self, payload: CasadyPayload,
                 start_index_input=0, start_index_state=0,
                 start_index_output=0):
        self.payload = payload
        self.output = graph_from_bytes(payload.output)
        self.step = graph_from_bytes(payload.dynamics)
        self.num_inputs = sum(desc.size for desc in payload.inputDescription)
        self.num_states = sum(desc.size for desc in payload.stateDescription)
        self.num_outputs = sum(desc.size for desc in payload.outputDescription)
        self.start_index_input = start_index_input
        self.start_index_state = start_index_state
        self.start_index_output = start_index_output
        self.end_index_input = start_index_input + self.num_inputs
        self.end_index_state = start_index_state + self.num_states
        self.end_index_output = start_index_output + self.num_outputs

    def slice_input(self):
        return slice(self.start_index_input, self.end_index_input)

    def slice_state(self):
        return slice(self.start_index_state, self.end_index_state)

    def slice_output(self):
        return slice(self.start_index_output, self.end_index_output)


def graph_from_bytes(graph_bytes: bytes) -> ca.Function:
    """
    Load a CasADi graph from a string representation.

    Args:
        graph_bytes (bytes): The byte representation of the CasADi graph.
    """
    return ca.Function.deserialize(graph_bytes.decode('utf-8'))

def graph_to_bytes(graph: ca.Function) -> bytes:
    """
    Serialize a CasADi graph to a byte representation.

    Args:
        graph (ca.Function): The CasADi graph to serialize.
    """
    return graph.serialize().encode('utf-8')

def vessel_model_casadi(
        mass: float,
        inertia: List[float],
        drag_coefficient: float,
        added_mass: List[float],
        center_of_gravity: List[float],
        center_of_buoyancy: List[float],
        water_density: float,
        gravity: float) -> ca.Function:

    state = ca.SX.sym('state', 12)
    tau = ca.SX.sym('input', 6)

    x, y, z = state[0], state[1], state[2]
    roll, pitch, yaw = state[3], state[4], state[5]
    vx, vy, vz = state[6], state[7], state[8]
    wx, wy, wz = state[9], state[10], state[11]

    s_roll, c_roll = ca.sin(roll), ca.cos(roll)
    s_pitch, c_pitch = ca.sin(pitch), ca.cos(pitch)
    s_yaw, c_yaw = ca.sin(yaw), ca.cos(yaw)

    body_velocity = ca.vertcat(vx, vy, vz)
    body_angular_rate = ca.vertcat(wx, wy, wz)

    # In the NED convention, the force vector is expressed in the world/NED frame,
    # so the translational state derivative is taken directly in that frame.
    # The body-rate kinematics still use the Euler-angle Jacobian for orientation.
    euler_rates = ca.vertcat(
        ca.horzcat(1.0, s_roll * ca.tan(pitch), c_roll * ca.tan(pitch)),
        ca.horzcat(0.0, c_roll, -s_roll),
        ca.horzcat(0.0, s_roll / ca.cos(pitch), c_roll / ca.cos(pitch)),
    ) @ body_angular_rate

    total_mass = ca.vertcat(
        mass + added_mass[0],
        mass + added_mass[1],
        mass + added_mass[2],
    )
    total_inertia = ca.vertcat(inertia[0], inertia[1], inertia[2])

    control_force = tau[0:3]
    control_moment = tau[3:6]

    # Neutral buoyancy approximation: the vessel has a displaced volume equivalent
    # to the mass divided by water density. In the NED frame, +z points downward,
    # so gravity acts in +z while buoyancy acts in -z.
    displaced_volume = mass / (water_density + 1e-6)
    buoyancy_force = ca.vertcat(0.0, 0.0, -water_density * gravity * displaced_volume)
    gravity_force = ca.vertcat(0.0, 0.0, mass * gravity)

    cg = ca.vertcat(center_of_gravity[0], center_of_gravity[1], center_of_gravity[2])
    cb = ca.vertcat(center_of_buoyancy[0], center_of_buoyancy[1], center_of_buoyancy[2])

    restoring_moment = ca.cross(cg, gravity_force) + ca.cross(cb, buoyancy_force)

    drag_force = -drag_coefficient * body_velocity
    drag_moment = -drag_coefficient * body_angular_rate

    force_total = control_force + gravity_force + buoyancy_force + drag_force
    moment_total = control_moment + restoring_moment + drag_moment

    linear_acceleration = force_total / total_mass
    angular_acceleration = moment_total / total_inertia

    position_dot = body_velocity
    orientation_dot = euler_rates

    state_dot = ca.vertcat(
        position_dot[0],
        position_dot[1],
        position_dot[2],
        orientation_dot[0],
        orientation_dot[1],
        orientation_dot[2],
        linear_acceleration[0],
        linear_acceleration[1],
        linear_acceleration[2],
        angular_acceleration[0],
        angular_acceleration[1],
        angular_acceleration[2],
    )

    return ca.Function(
        'vessel_model',
        [state, tau],
        [state_dot],
        ['state', 'input'],
        ['state_dot'],
    )
