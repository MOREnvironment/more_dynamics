import casadi as ca
import numpy as np
import pytest

from more_dynamics.models.hydrodynamics import (
    crossflow_surface_hydrodynamics_casadi,
    preprocess_crossflow_surface_hydrodynamics,
)


PARAMETERS = {
    "damping_coefficients": [240.0, 1800.0, 7000.0, 900.0, 6500.0, 8000.0],
    "length": 3.5,
    "beam": 1.2,
    "draft": 0.16,
    "water_density": 1025.0,
    "crossflow_scale": 1.0,
    "strip_count": 20,
}


def test_crossflow_model_returns_zero_force_at_rest():
    properties = preprocess_crossflow_surface_hydrodynamics(**PARAMETERS)
    model = crossflow_surface_hydrodynamics_casadi(properties)

    force = np.asarray(model(ca.DM.zeros(6))).reshape(-1)

    np.testing.assert_allclose(force, np.zeros(6), atol=1e-12)


def test_crossflow_model_dissipates_sway_and_yaw_motion():
    properties = preprocess_crossflow_surface_hydrodynamics(**PARAMETERS)
    model = crossflow_surface_hydrodynamics_casadi(properties)
    velocity = np.array([20.0, 0.5, 0.0, 0.0, 0.0, 0.2])

    force = np.asarray(model(ca.DM(velocity))).reshape(-1)

    assert force[1] < 0.0
    assert force[5] < 0.0
    assert velocity @ force < 0.0


def test_crossflow_model_rejects_invalid_strip_count():
    with pytest.raises(ValueError, match="at least 2"):
        preprocess_crossflow_surface_hydrodynamics(
            **(PARAMETERS | {"strip_count": 1})
        )
