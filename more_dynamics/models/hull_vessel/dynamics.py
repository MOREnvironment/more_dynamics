import casadi as ca

from .mass_properties import HullMassProperties


def vessel_model_casadi(
    mass_properties: HullMassProperties,
    include_coriolis: bool = True,
) -> ca.Function:
    """Build the current dynamics graph with a preprocessed mass matrix."""
    state = ca.SX.sym("state", 12)
    tau = ca.SX.sym("input", 6)

    roll, pitch, yaw = state[3], state[4], state[5]
    body_linear_velocity = state[6:9]
    body_angular_rate = state[9:12]

    s_roll, c_roll = ca.sin(roll), ca.cos(roll)
    s_pitch, c_pitch = ca.sin(pitch), ca.cos(pitch)
    s_yaw, c_yaw = ca.sin(yaw), ca.cos(yaw)
    body_to_global = ca.vertcat(
        ca.horzcat(
            c_yaw * c_pitch,
            c_yaw * s_pitch * s_roll - s_yaw * c_roll,
            c_yaw * s_pitch * c_roll + s_yaw * s_roll,
        ),
        ca.horzcat(
            s_yaw * c_pitch,
            s_yaw * s_pitch * s_roll + c_yaw * c_roll,
            s_yaw * s_pitch * c_roll - c_yaw * s_roll,
        ),
        ca.horzcat(
            -s_pitch,
            c_pitch * s_roll,
            c_pitch * c_roll,
        ),
    )
    global_velocity = body_to_global @ body_linear_velocity
    euler_rates = ca.vertcat(
        ca.horzcat(1.0, s_roll * ca.tan(pitch), c_roll * ca.tan(pitch)),
        ca.horzcat(0.0, c_roll, -s_roll),
        ca.horzcat(0.0, s_roll / ca.cos(pitch), c_roll / ca.cos(pitch)),
    ) @ body_angular_rate

    total_mass_matrix = ca.DM(mass_properties.total_mass_matrix)
    if include_coriolis:
        body_velocity = ca.vertcat(
            body_linear_velocity,
            body_angular_rate,
        )
        momentum = total_mass_matrix @ body_velocity
        linear_momentum = momentum[:3]
        angular_momentum = momentum[3:]

        # The body-frame momentum balance includes both rigid-body and
        # added-mass Coriolis/centripetal terms. Deriving this from the
        # complete mass matrix preserves their centre-of-gravity coupling.
        coriolis_wrench = ca.vertcat(
            ca.cross(body_angular_rate, linear_momentum),
            ca.cross(body_linear_velocity, linear_momentum)
            + ca.cross(body_angular_rate, angular_momentum),
        )
    else:
        coriolis_wrench = ca.SX.zeros(6)
    acceleration = ca.solve(total_mass_matrix, tau - coriolis_wrench)

    state_dot = ca.vertcat(global_velocity, euler_rates, acceleration)
    return ca.Function(
        "vessel_model",
        [state, tau],
        [state_dot],
        ["state", "input"],
        ["state_dot"],
    )
