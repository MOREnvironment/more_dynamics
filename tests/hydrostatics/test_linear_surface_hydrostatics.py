import casadi as ca
import numpy as np

from more_dynamics.models.hydrostatics import (
    linear_surface_hydrostatics_casadi,
    preprocess_linear_surface_hydrostatics,
)


DEFAULT_PARAMETERS = {
    "length": 5.2,
    "beam": 2.15,
    "draft": 0.3,
    "block_coefficient": 0.233,
    "waterplane_coefficient": 0.8,
    "water_density": 1025.0,
    "gravity": 9.81,
    "center_of_gravity": [0.0, 0.0, 0.025],
    "longitudinal_center_of_flotation": -0.3,
    "coefficient_scales": [1.0, 0.7, 1.0],
    "reference_point": [0.0, 0.0, 0.0],
}


def test_preprocessing_matches_hull_usv_surface_restoring_formulas():
    properties = preprocess_linear_surface_hydrostatics(
        **DEFAULT_PARAMETERS
    )

    expected_volume = 0.233 * 5.2 * 2.15 * 0.3
    expected_area = 5.2 * 2.15 + 2.0 * 0.3 * 2.15

    np.testing.assert_allclose(properties.displaced_volume, expected_volume)
    np.testing.assert_allclose(properties.wetted_surface_area, expected_area)
    np.testing.assert_allclose(
        properties.stiffness_matrix,
        properties.stiffness_matrix.T,
    )
    np.testing.assert_allclose(
        properties.stiffness_matrix,
        np.array(
            [
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 89934.156, 0.0, 26980.2468, 0.0],
                [0.0, 0.0, 0.0, 28000.02855884294, 0.0, 0.0],
                [0.0, 0.0, 26980.2468, 0.0, 184988.9375388429, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            ]
        ),
        rtol=1e-10,
    )


def test_casadi_graph_returns_signed_force_volume_and_area():
    properties = preprocess_linear_surface_hydrostatics(
        **DEFAULT_PARAMETERS
    )
    model = linear_surface_hydrostatics_casadi(properties)
    pose = np.array([1.0, -2.0, 0.15, 0.04, -0.03, 0.2])

    output = np.asarray(model(ca.DM(pose))).reshape(-1)

    np.testing.assert_allclose(
        output[:6],
        -properties.stiffness_matrix @ pose,
    )
    np.testing.assert_allclose(output[6], properties.displaced_volume)
    np.testing.assert_allclose(output[7], properties.wetted_surface_area)
