import casadi as ca
import numpy as np

from rpp_py.data_manager import DataManager
_DATA_MANAGER = DataManager()

from more_dynamics.models.hull_vessel import (
    preprocess_hull_mass_properties,
    vessel_model_casadi,
)


DEFAULT_PARAMETERS = {
    "length": 5.2,
    "beam": 2.15,
    "draft": 0.3,
    "block_coefficient": 0.233,
    "water_density": 1025.0,
    "radii_of_gyration": [0.35, 0.25, 0.25],
    "center_of_gravity": [0.0, 0.0, 0.025],
    "center_of_buoyancy": [0.0, 0.0, 0.0],
    "added_mass_coefficients": [-1.0, -1.5, -1.0, -0.2, -0.8, -1.7],
}


def test_preprocessing_matches_hull_usv_mass_formulas():
    properties = preprocess_hull_mass_properties(**DEFAULT_PARAMETERS)

    expected_volume = 0.233 * 5.2 * 2.15 * 0.3
    expected_mass = 1025.0 * expected_volume
    expected_radii = np.array([0.35 * 2.15, 0.25 * 5.2, 0.25 * 5.2])

    np.testing.assert_allclose(properties.displaced_volume, expected_volume)
    np.testing.assert_allclose(properties.mass, expected_mass)
    np.testing.assert_allclose(
        properties.inertia_at_center_of_gravity,
        np.diag(expected_mass * expected_radii**2),
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
