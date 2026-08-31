from dataclasses import dataclass
from typing import Sequence

import casadi as ca
import numpy as np


@dataclass(frozen=True)
class LinearSurfaceHydrodynamicsProperties:
    """Constant linear damping properties."""

    damping_matrix: np.ndarray


def preprocess_linear_surface_hydrodynamics(
    damping_coefficients: Sequence[float],
) -> LinearSurfaceHydrodynamicsProperties:
    coefficients = np.asarray(damping_coefficients, dtype=float)
    if coefficients.shape != (6,):
        raise ValueError("damping_coefficients must contain exactly 6 values")
    if not np.all(np.isfinite(coefficients)):
        raise ValueError("damping_coefficients must contain only finite values")
    if np.any(coefficients < 0.0):
        raise ValueError("damping_coefficients must be non-negative")

    return LinearSurfaceHydrodynamicsProperties(
        damping_matrix=np.diag(coefficients),
    )


def linear_surface_hydrodynamics_casadi(
    properties: LinearSurfaceHydrodynamicsProperties,
) -> ca.Function:
    """Build the signed linear damping force ``velocity -> -D @ velocity``."""
    velocity = ca.SX.sym("velocity", 6)
    hydrodynamic_force = -ca.DM(properties.damping_matrix) @ velocity
    return ca.Function(
        "linear_surface_hydrodynamics",
        [velocity],
        [hydrodynamic_force],
        ["velocity"],
        ["hydrodynamic_force"],
    )
