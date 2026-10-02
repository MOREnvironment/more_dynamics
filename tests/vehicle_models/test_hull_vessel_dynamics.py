import casadi as ca
import numpy as np

from rpp_py.data_manager import DataManager
_DATA_MANAGER = DataManager()

from more_dynamics.models.hull_vessel import (
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
