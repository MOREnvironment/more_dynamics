from dataclasses import dataclass
from typing import Sequence

import casadi as ca
import numpy as np


@dataclass(frozen=True)
class LinearSurfaceHydrostaticsProperties:
    """Constants used by the linear surface-vessel hydrostatics graph."""

    stiffness_matrix: np.ndarray
    displaced_volume: float
    wetted_surface_area: float


def _vector(values: Sequence[float], size: int, name: str) -> np.ndarray:
    vector = np.asarray(values, dtype=float)
    if vector.shape != (size,):
        raise ValueError(f"{name} must contain exactly {size} values")
    if not np.all(np.isfinite(vector)):
        raise ValueError(f"{name} must contain only finite values")
    return vector


def _skew(vector: np.ndarray) -> np.ndarray:
    x, y, z = vector
    return np.array(
        [[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]],
        dtype=float,
    )


def _transform(vector: np.ndarray) -> np.ndarray:
    return np.block(
        [
            [np.eye(3), _skew(vector).T],
            [np.zeros((3, 3)), np.eye(3)],
        ]
    )


def preprocess_linear_surface_hydrostatics(
    length: float,
    beam: float,
    draft: float,
    block_coefficient: float,
    waterplane_coefficient: float,
    water_density: float,
    gravity: float,
    center_of_gravity: Sequence[float],
    longitudinal_center_of_flotation: float,
    coefficient_scales: Sequence[float],
    reference_point: Sequence[float],
) -> LinearSurfaceHydrostaticsProperties:
    """Port the constant ``HullUSV`` surface-vessel restoring model."""
    positive_values = {
        "length": length,
        "beam": beam,
        "draft": draft,
        "block_coefficient": block_coefficient,
        "waterplane_coefficient": waterplane_coefficient,
        "water_density": water_density,
        "gravity": gravity,
    }
    for name, value in positive_values.items():
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be a positive finite value")

    if not np.isfinite(longitudinal_center_of_flotation):
        raise ValueError(
            "longitudinal_center_of_flotation must be finite"
        )

    center_of_gravity_vector = _vector(
        center_of_gravity,
        3,
        "center_of_gravity",
    )
    scales = _vector(coefficient_scales, 3, "coefficient_scales")
    reference_point_vector = _vector(reference_point, 3, "reference_point")

    displaced_volume = block_coefficient * length * beam * draft
    unscaled_waterplane_area = waterplane_coefficient * length * beam
    transverse_waterplane_inertia = (
        scales[0]
        * length
        * beam**3
        / 12.0
        * (
            6.0
            * waterplane_coefficient**3
            / (
                (1.0 + waterplane_coefficient)
                * (1.0 + 2.0 * waterplane_coefficient)
            )
        )
    )
    longitudinal_waterplane_inertia = (
        scales[1] * beam * length**3 / 12.0
    )
    center_of_buoyancy_above_keel = (
        5.0 * draft / 2.0
        - scales[2] * displaced_volume / unscaled_waterplane_area
    ) / 3.0
    center_of_gravity_above_keel = draft - center_of_gravity_vector[2]
    transverse_metacentric_height = (
        center_of_buoyancy_above_keel
        + transverse_waterplane_inertia / displaced_volume
        - center_of_gravity_above_keel
    )
    longitudinal_metacentric_height = (
        center_of_buoyancy_above_keel
        + longitudinal_waterplane_inertia / displaced_volume
        - center_of_gravity_above_keel
    )
    waterplane_area = scales[0] * unscaled_waterplane_area

    stiffness_at_center_of_flotation = np.diag(
        [
            0.0,
            0.0,
            water_density * gravity * waterplane_area,
            water_density
            * gravity
            * displaced_volume
            * transverse_metacentric_height,
            water_density
            * gravity
            * displaced_volume
            * longitudinal_metacentric_height,
            0.0,
        ]
    )
    center_of_flotation = np.array(
        [longitudinal_center_of_flotation, 0.0, 0.0]
    )
    transform_center_of_flotation = _transform(center_of_flotation)
    stiffness_at_origin = (
        transform_center_of_flotation.T
        @ stiffness_at_center_of_flotation
        @ transform_center_of_flotation
    )
    transform_reference_point = _transform(reference_point_vector)
    stiffness_matrix = (
        transform_reference_point.T
        @ stiffness_at_origin
        @ transform_reference_point
    )

    return LinearSurfaceHydrostaticsProperties(
        stiffness_matrix=stiffness_matrix,
        displaced_volume=displaced_volume,
        wetted_surface_area=length * beam + 2.0 * draft * beam,
    )


def linear_surface_hydrostatics_casadi(
    properties: LinearSurfaceHydrostaticsProperties,
) -> ca.Function:
    """Build ``eta -> [g(eta), displaced volume, wetted area]``."""
    pose = ca.SX.sym("pose", 6)
    restoring_force = -ca.DM(properties.stiffness_matrix) @ pose
    output = ca.vertcat(
        restoring_force,
        properties.displaced_volume,
        properties.wetted_surface_area,
    )
    return ca.Function(
        "linear_surface_hydrostatics",
        [pose],
        [output],
        ["pose"],
        ["output"],
    )
