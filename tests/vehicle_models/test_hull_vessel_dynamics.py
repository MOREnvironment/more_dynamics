import casadi as ca
import numpy as np

from rpp_py.data_manager import DataManager
_DATA_MANAGER = DataManager()

from more_dynamics.models.hull_vessel import (
    coriolis_matrices_casadi,
    preprocess_hull_mass_properties,
    vessel_model_casadi,
)


DEFAULT_PARAMETERS = {
    "mass": 801.01905,
    "inertia": [
        453.5820434315624,
        1353.7221945000001,
        1353.7221945000001,
    ],
    "center_of_gravity": [0.0, 0.0, 0.025],
    "added_mass": [
        67.85979986031616,
        1201.528575,
        801.01905,
        90.81653606756248,
        1083.378265125,
        2301.32773065,
    ],
}


def test_preprocessing_uses_explicit_hull_mass_properties():
    properties = preprocess_hull_mass_properties(**DEFAULT_PARAMETERS)

    np.testing.assert_allclose(properties.mass, DEFAULT_PARAMETERS["mass"])
    np.testing.assert_allclose(
        properties.inertia_at_center_of_gravity,
        np.diag(DEFAULT_PARAMETERS["inertia"]),
    )
    expected_total_mass = np.array(
        [
            [868.8788498603161, 0.0, 0.0, 0.0, 20.02547625, 0.0],
            [0.0, 2002.547625, 0.0, -20.02547625, 0.0, 0.0],
            [0.0, 0.0, 1602.0381, 0.0, 0.0, 0.0],
            [0.0, -20.02547625, 0.0, 544.8992164053749, 0.0, 0.0],
            [20.02547625, 0.0, 0.0, 0.0, 2437.6010965312503, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 3655.04992515],
        ]
    )
    np.testing.assert_allclose(
        properties.total_mass_matrix,
        expected_total_mass,
    )
    np.testing.assert_allclose(
        properties.total_mass_matrix,
        properties.total_mass_matrix.T,
    )
    assert np.all(np.linalg.eigvalsh(properties.total_mass_matrix) > 0.0)


def test_casadi_graph_uses_preprocessed_total_mass_matrix():
    parameters = DEFAULT_PARAMETERS | {
        "center_of_gravity": [0.0, 0.0, 0.0]
    }
    properties = preprocess_hull_mass_properties(**parameters)
    model = vessel_model_casadi(
        mass_properties=properties,
    )
    force = np.array([10.0, -4.0, 2.0, 1.0, -3.0, 5.0])

    state_dot = np.asarray(model(ca.DM.zeros(12), force)).reshape(-1)

    np.testing.assert_allclose(state_dot[:6], np.zeros(6), atol=1e-12)
    np.testing.assert_allclose(
        state_dot[6:],
        np.linalg.solve(properties.total_mass_matrix, force),
    )


def test_casadi_graph_rotates_body_velocity_into_global_frame():
    properties = preprocess_hull_mass_properties(**DEFAULT_PARAMETERS)
    model = vessel_model_casadi(mass_properties=properties)
    state = np.zeros(12)
    state[5] = np.pi / 2.0
    state[6:9] = [4.0, 1.0, 0.0]

    state_dot = np.asarray(model(state, np.zeros(6))).reshape(-1)

    np.testing.assert_allclose(state_dot[:3], [-1.0, 4.0, 0.0])


def test_casadi_graph_applies_body_frame_centripetal_acceleration():
    parameters = DEFAULT_PARAMETERS | {
        "center_of_gravity": [0.0, 0.0, 0.0]
    }
    properties = preprocess_hull_mass_properties(**parameters)
    model = vessel_model_casadi(mass_properties=properties)
    state = np.zeros(12)
    surge_velocity = 8.0
    yaw_rate = 0.7
    state[6] = surge_velocity
    state[11] = yaw_rate

    state_dot = np.asarray(model(state, np.zeros(6))).reshape(-1)
    expected_sway_acceleration = -(
        properties.total_mass_matrix[0, 0]
        / properties.total_mass_matrix[1, 1]
        * surge_velocity
        * yaw_rate
    )

    np.testing.assert_allclose(
        state_dot[7],
        expected_sway_acceleration,
    )


def test_casadi_graph_can_disable_coriolis_acceleration():
    properties = preprocess_hull_mass_properties(**DEFAULT_PARAMETERS)
    model = vessel_model_casadi(
        mass_properties=properties,
        include_coriolis=False,
    )
    state = np.zeros(12)
    state[6] = 8.0
    state[11] = 0.7

    state_dot = np.asarray(model(state, np.zeros(6))).reshape(-1)

    np.testing.assert_allclose(state_dot[6:], np.zeros(6), atol=1e-12)


def test_casadi_graph_coriolis_force_does_not_create_kinetic_energy():
    properties = preprocess_hull_mass_properties(**DEFAULT_PARAMETERS)
    model = vessel_model_casadi(mass_properties=properties)
    state = np.zeros(12)
    state[6:] = [4.0, -1.0, 0.2, 0.3, -0.1, 0.7]

    state_dot = np.asarray(model(state, np.zeros(6))).reshape(-1)
    velocity = state[6:]
    kinetic_energy_rate = velocity @ properties.total_mass_matrix @ state_dot[6:]

    np.testing.assert_allclose(kinetic_energy_rate, 0.0, atol=1e-12)


def test_coriolis_decomposition_matches_fossen_rigid_and_added_mass_terms():
    properties = preprocess_hull_mass_properties(**DEFAULT_PARAMETERS)
    velocity = np.array([4.0, -1.0, 0.2, 0.3, -0.1, 0.7])
    velocity_symbol = ca.SX.sym("velocity", 6)
    rigid_body, added_mass = coriolis_matrices_casadi(
        properties,
        velocity_symbol,
    )
    model = ca.Function("coriolis", [velocity_symbol], [rigid_body, added_mass])

    rigid_body_value, added_mass_value = model(ca.DM(velocity))
    rigid_body_value = np.asarray(rigid_body_value)
    added_mass_value = np.asarray(added_mass_value)

    def skew(vector):
        x, y, z = vector
        return np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])

    center_of_gravity = properties.center_of_gravity
    transform = np.block(
        [
            [np.eye(3), skew(center_of_gravity).T],
            [np.zeros((3, 3)), np.eye(3)],
        ]
    )
    angular_velocity = velocity[3:]
    rigid_body_at_cg = np.block(
        [
            [
                properties.mass * skew(angular_velocity),
                np.zeros((3, 3)),
            ],
            [
                np.zeros((3, 3)),
                -skew(properties.inertia_at_center_of_gravity @ angular_velocity),
            ],
        ]
    )
    expected_rigid_body = transform.T @ rigid_body_at_cg @ transform
    added_linear_momentum = properties.added_mass_matrix[:3] @ velocity
    added_angular_momentum = properties.added_mass_matrix[3:] @ velocity
    expected_added_mass = np.block(
        [
            [np.zeros((3, 3)), -skew(added_linear_momentum)],
            [-skew(added_linear_momentum), -skew(added_angular_momentum)],
        ]
    )

    np.testing.assert_allclose(rigid_body_value, expected_rigid_body)
    np.testing.assert_allclose(added_mass_value, expected_added_mass)
    np.testing.assert_allclose(
        (rigid_body_value + added_mass_value)
        + (rigid_body_value + added_mass_value).T,
        np.zeros((6, 6)),
        atol=1e-12,
    )


def test_high_speed_turn_couples_into_roll_for_an_elevated_centre_of_gravity():
    centered_parameters = DEFAULT_PARAMETERS | {
        "center_of_gravity": [0.0, 0.0, 0.0]
    }
    elevated_properties = preprocess_hull_mass_properties(**DEFAULT_PARAMETERS)
    centered_properties = preprocess_hull_mass_properties(**centered_parameters)
    elevated_model = vessel_model_casadi(mass_properties=elevated_properties)
    centered_model = vessel_model_casadi(mass_properties=centered_properties)
    state = np.zeros(12)
    state[6] = 8.0
    state[11] = 0.7

    elevated_state_dot = np.asarray(
        elevated_model(state, np.zeros(6))
    ).reshape(-1)
    centered_state_dot = np.asarray(
        centered_model(state, np.zeros(6))
    ).reshape(-1)

    assert abs(elevated_state_dot[9]) > 1e-12
    np.testing.assert_allclose(centered_state_dot[9], 0.0, atol=1e-12)
