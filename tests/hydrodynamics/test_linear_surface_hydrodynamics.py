import casadi as ca
import numpy as np
import pytest

from more_dynamics.models.hydrodynamics import (
    linear_surface_hydrodynamics_casadi,
    preprocess_linear_surface_hydrodynamics,
)


DEFAULT_COEFFICIENTS = np.array(
    [
        0.0,
        20.02547625,
        7201.948346468731,
        1562.418311261956,
        14535.734088287643,
        3655.04992515,
    ]
)


def test_preprocessing_builds_original_hull_usv_damping_matrix():
    properties = preprocess_linear_surface_hydrodynamics(
        DEFAULT_COEFFICIENTS
    )

    np.testing.assert_allclose(
        properties.damping_matrix,
        np.diag(DEFAULT_COEFFICIENTS),
    )


def test_casadi_graph_returns_signed_linear_damping_force():
    properties = preprocess_linear_surface_hydrodynamics(
        DEFAULT_COEFFICIENTS
    )
    model = linear_surface_hydrodynamics_casadi(properties)
    velocity = np.array([2.0, -1.0, 0.2, -0.1, 0.05, -0.3])

    force = np.asarray(model(ca.DM(velocity))).reshape(-1)

    np.testing.assert_allclose(
        force,
        -properties.damping_matrix @ velocity,
    )


def test_preprocessing_rejects_negative_damping():
    coefficients = DEFAULT_COEFFICIENTS.copy()
    coefficients[2] = -1.0

    with pytest.raises(ValueError, match="non-negative"):
        preprocess_linear_surface_hydrodynamics(coefficients)
