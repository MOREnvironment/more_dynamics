import casadi as ca

from .mass_properties import HullMassProperties


def vessel_model_casadi(
    mass_properties: HullMassProperties,
) -> ca.Function:
    """Build the current dynamics graph with a preprocessed mass matrix."""
    state = ca.SX.sym("state", 12)
    tau = ca.SX.sym("input", 6)

    roll, pitch = state[3], state[4]
    body_velocity = state[6:9]
    body_angular_rate = state[9:12]

    s_roll, c_roll = ca.sin(roll), ca.cos(roll)
    euler_rates = ca.vertcat(
        ca.horzcat(1.0, s_roll * ca.tan(pitch), c_roll * ca.tan(pitch)),
        ca.horzcat(0.0, c_roll, -s_roll),
        ca.horzcat(0.0, s_roll / ca.cos(pitch), c_roll / ca.cos(pitch)),
    ) @ body_angular_rate

    acceleration = ca.solve(
        ca.DM(mass_properties.total_mass_matrix),
        tau,
    )

    state_dot = ca.vertcat(body_velocity, euler_rates, acceleration)
    return ca.Function(
        "vessel_model",
        [state, tau],
        [state_dot],
        ["state", "input"],
        ["state_dot"],
    )
