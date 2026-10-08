from dataclasses import dataclass
from typing import Sequence

import casadi as ca
import numpy as np


@dataclass(frozen=True)
class CrossflowSurfaceHydrodynamicsProperties:
    """Constants for linear damping plus surface-vessel cross-flow drag."""

    damping_matrix: np.ndarray
    strip_positions: np.ndarray
    strip_width: float
    crossflow_coefficient: float
    water_density: float
    draft: float
    crossflow_scale: float


_HOERNER_ASPECT_RATIOS = np.array(
    [
        0.0108623, 0.176606, 0.353025, 0.451863, 0.472838, 0.492877,
        0.493252, 0.558473, 0.646401, 0.833589, 0.988002, 1.30807,
        1.63918, 1.85998, 2.31288, 2.59998, 3.00877, 3.45075, 3.7379,
        4.00309,
    ],
    dtype=float,
)
_HOERNER_COEFFICIENTS = np.array(
    [
        1.96608, 1.96573, 1.89756, 1.78718, 1.58374, 1.27862, 1.21082,
        1.08356, 0.998631, 0.87959, 0.828415, 0.759941, 0.691442,
        0.657076, 0.630693, 0.596186, 0.586846, 0.585909, 0.559877,
        0.559315,
    ],
    dtype=float,
)


def _positive_finite(value: float, name: str) -> float:
    value = float(value)
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be a positive finite value")
    return value


def preprocess_crossflow_surface_hydrodynamics(
    damping_coefficients: Sequence[float],
    length: float,
    beam: float,
    draft: float,
    water_density: float,
    crossflow_scale: float,
    strip_count: int,
) -> CrossflowSurfaceHydrodynamicsProperties:
    """Preprocess the Fossen/Hoerner strip-theory cross-flow model.

    The formulation is the surface-vessel cross-flow term from
    ``more_generic_models.HydroForces.cross_flow_drag``. It supplements the
    configured linear damping matrix.
    """
    coefficients = np.asarray(damping_coefficients, dtype=float)
    if coefficients.shape != (6,):
        raise ValueError("damping_coefficients must contain exactly 6 values")
    if not np.all(np.isfinite(coefficients)):
        raise ValueError("damping_coefficients must contain only finite values")
    if np.any(coefficients < 0.0):
        raise ValueError("damping_coefficients must be non-negative")

    length = _positive_finite(length, "length")
    beam = _positive_finite(beam, "beam")
    draft = _positive_finite(draft, "draft")
    water_density = _positive_finite(water_density, "water_density")
    crossflow_scale = float(crossflow_scale)
    if not np.isfinite(crossflow_scale) or crossflow_scale < 0.0:
        raise ValueError("crossflow_scale must be a non-negative finite value")
    if isinstance(strip_count, bool) or int(strip_count) != strip_count:
        raise ValueError("strip_count must be an integer")
    strip_count = int(strip_count)
    if strip_count < 2:
        raise ValueError("strip_count must be at least 2")

    aspect_ratio = beam / (2.0 * draft)
    crossflow_coefficient = float(
        np.interp(
            aspect_ratio,
            _HOERNER_ASPECT_RATIOS,
            _HOERNER_COEFFICIENTS,
        )
    )
    strip_width = length / strip_count
    strip_positions = np.linspace(
        -length / 2.0,
        length / 2.0,
        strip_count + 1,
    )

    return CrossflowSurfaceHydrodynamicsProperties(
        damping_matrix=np.diag(coefficients),
        strip_positions=strip_positions,
        strip_width=strip_width,
        crossflow_coefficient=crossflow_coefficient,
        water_density=water_density,
        draft=draft,
        crossflow_scale=crossflow_scale,
    )


def crossflow_surface_hydrodynamics_casadi(
    properties: CrossflowSurfaceHydrodynamicsProperties,
) -> ca.Function:
    """Build ``nu -> -D nu + tau_crossflow(nu)`` for an FRD vessel body."""
    velocity = ca.SX.sym("velocity", 6)
    surge, sway, heave, roll_rate, pitch_rate, yaw_rate = (
        velocity[index] for index in range(6)
    )
    strip_positions = ca.DM(properties.strip_positions)
    lateral_strip_velocity = sway + strip_positions * yaw_rate
    vertical_strip_velocity = heave + strip_positions * pitch_rate
    lateral_quadratic_velocity = (
        lateral_strip_velocity * ca.fabs(lateral_strip_velocity)
    )
    vertical_quadratic_velocity = (
        vertical_strip_velocity * ca.fabs(vertical_strip_velocity)
    )
    strip_force_scale = (
        -0.5
        * properties.water_density
        * properties.draft
        * properties.crossflow_coefficient
        * properties.strip_width
        * properties.crossflow_scale
    )
    sway_force = strip_force_scale * ca.sum1(lateral_quadratic_velocity)
    heave_force = strip_force_scale * ca.sum1(vertical_quadratic_velocity)
    pitch_moment = strip_force_scale * ca.dot(
        strip_positions,
        vertical_quadratic_velocity,
    )
    yaw_moment = strip_force_scale * ca.dot(
        strip_positions,
        lateral_quadratic_velocity,
    )
    linear_force = -ca.DM(properties.damping_matrix) @ velocity
    force = linear_force + ca.vertcat(
        0.0,
        sway_force,
        heave_force,
        0.0,
        pitch_moment,
        yaw_moment,
    )
    return ca.Function(
        "crossflow_surface_hydrodynamics",
        [velocity],
        [force],
        ["velocity"],
        ["hydrodynamic_force"],
    )
