from dataclasses import dataclass
from typing import Sequence

import numpy as np


@dataclass(frozen=True)
class HullMassProperties:
    """Constant hull properties computed before building the CasADi graph."""

    mass: float
    displaced_volume: float
    wetted_surface_area: float
    inertia_at_center_of_gravity: np.ndarray
    inertia_at_center_of_origin: np.ndarray
    rigid_body_mass_matrix: np.ndarray
    added_mass_matrix: np.ndarray
    total_mass_matrix: np.ndarray


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


def preprocess_hull_mass_properties(
    length: float,
    beam: float,
    draft: float,
    block_coefficient: float,
    water_density: float,
    radii_of_gyration: Sequence[float],
    center_of_gravity: Sequence[float],
    center_of_buoyancy: Sequence[float],
    added_mass_coefficients: Sequence[float],
) -> HullMassProperties:
    """Compute the constant mass properties used by ``HullUSV``."""
    positive_values = {
        "length": length,
        "beam": beam,
        "draft": draft,
        "block_coefficient": block_coefficient,
        "water_density": water_density,
    }
    for name, value in positive_values.items():
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be a positive finite value")

    radius_scale = _vector(radii_of_gyration, 3, "radii_of_gyration")
    r_cg = _vector(center_of_gravity, 3, "center_of_gravity")
    r_cb = _vector(center_of_buoyancy, 3, "center_of_buoyancy")
    added_mass_scale = _vector(
        added_mass_coefficients,
        6,
        "added_mass_coefficients",
    )

    displaced_volume = block_coefficient * length * beam * draft
    mass = water_density * displaced_volume
    radii = radius_scale * np.array([beam, length, length], dtype=float)
    inertia_cg = np.diag(mass * np.square(radii))

    r_bg = r_cg - r_cb
    skew_r_bg = _skew(r_bg)
    inertia_co = inertia_cg - mass * (skew_r_bg @ skew_r_bg)

    transform = np.block(
        [
            [np.eye(3), skew_r_bg.T],
            [np.zeros((3, 3)), np.eye(3)],
        ]
    )
    mass_at_cg = np.block(
        [
            [mass * np.eye(3), np.zeros((3, 3))],
            [np.zeros((3, 3)), inertia_cg],
        ]
    )
    rigid_body_mass = transform.T @ mass_at_cg @ transform

    surge_added_mass = (
        2.7
        * water_density
        * displaced_volume ** (5.0 / 3.0)
        / length**2
    )
    added_mass_derivatives = added_mass_scale * np.array(
        [
            surge_added_mass,
            mass,
            mass,
            inertia_co[0, 0],
            inertia_co[1, 1],
            inertia_co[2, 2],
        ]
    )
    added_mass = np.diag(-added_mass_derivatives)
    total_mass = rigid_body_mass + added_mass

    if not np.allclose(total_mass, total_mass.T):
        raise ValueError("total hull mass matrix must be symmetric")
    if not np.all(np.linalg.eigvalsh(total_mass) > 0.0):
        raise ValueError("total hull mass matrix must be positive definite")

    return HullMassProperties(
        mass=mass,
        displaced_volume=displaced_volume,
        wetted_surface_area=length * beam + 2.0 * draft * beam,
        inertia_at_center_of_gravity=inertia_cg,
        inertia_at_center_of_origin=inertia_co,
        rigid_body_mass_matrix=rigid_body_mass,
        added_mass_matrix=added_mass,
        total_mass_matrix=total_mass,
    )
