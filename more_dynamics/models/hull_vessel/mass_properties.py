from dataclasses import dataclass
from typing import Sequence

import numpy as np


@dataclass(frozen=True)
class HullMassProperties:
    """Constant rigid-body and added-mass properties for ``HullVessel``."""

    mass: float
    center_of_gravity: np.ndarray
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
    mass: float,
    inertia: Sequence[float],
    center_of_gravity: Sequence[float],
    added_mass: Sequence[float],
) -> HullMassProperties:
    """Build mass matrices from explicit rigid-body and added-mass values.

    Hull geometry belongs to hydrostatic models.  These values define only the
    vehicle inertia about the body-frame origin used by the dynamics state.
    """
    if not np.isfinite(mass) or mass <= 0.0:
        raise ValueError("mass must be a positive finite value")

    inertia_cg = _vector(inertia, 3, "inertia")
    if np.any(inertia_cg <= 0.0):
        raise ValueError("inertia must contain positive values")

    r_cg = _vector(center_of_gravity, 3, "center_of_gravity")
    added_mass_diagonal = _vector(added_mass, 6, "added_mass")
    if np.any(added_mass_diagonal < 0.0):
        raise ValueError("added_mass must contain non-negative values")

    inertia_cg_matrix = np.diag(inertia_cg)
    skew_r_cg = _skew(r_cg)
    inertia_co = inertia_cg_matrix - mass * (skew_r_cg @ skew_r_cg)

    transform = np.block(
        [
            [np.eye(3), skew_r_cg.T],
            [np.zeros((3, 3)), np.eye(3)],
        ]
    )
    mass_at_cg = np.block(
        [
            [mass * np.eye(3), np.zeros((3, 3))],
            [np.zeros((3, 3)), inertia_cg_matrix],
        ]
    )
    rigid_body_mass = transform.T @ mass_at_cg @ transform

    added_mass_matrix = np.diag(added_mass_diagonal)
    total_mass = rigid_body_mass + added_mass_matrix

    if not np.allclose(total_mass, total_mass.T):
        raise ValueError("total hull mass matrix must be symmetric")
    if not np.all(np.linalg.eigvalsh(total_mass) > 0.0):
        raise ValueError("total hull mass matrix must be positive definite")

    return HullMassProperties(
        mass=mass,
        center_of_gravity=r_cg,
        inertia_at_center_of_gravity=inertia_cg_matrix,
        inertia_at_center_of_origin=inertia_co,
        rigid_body_mass_matrix=rigid_body_mass,
        added_mass_matrix=added_mass_matrix,
        total_mass_matrix=total_mass,
    )
