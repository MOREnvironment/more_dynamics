"""Preprocessed constants of the rigid-body + added-mass block (numpy).

Ported from ``more_generic_models``
``dynamics/plant/matrices/rigid_body_kinetics_surface_vessel.py::get_params_catamaran``
(lines 17-51) and ``asv_catamaran.py::compute_constant_values`` (121-153).

Inertia (owner decision E-11 c, 2026-10-05): by default the inertia is taken
about the combined CG, as MSS ``otter.m`` lines 123-129 since its revision of
2026-04-20. ``legacy_otter_inertia=True`` reproduces the numpy source and the
frozen MATLAB reference: ``I_CG - m_hull S(r_g)^2 - m_payload S(r_p)^2``,
already about the CO and shifted a second time by ``H(r_g)``. The flag is
removed when the reference CSVs are regenerated.
"""

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .added_mass import added_mass_derivatives, added_mass_matrix
from .kinetics import rigid_body_mass_matrix, skew


@dataclass(frozen=True)
class RigidBodyConstants:
    """Constants computed before building the CasADi graph."""

    mass: float
    center_of_gravity: np.ndarray
    inertia_hull_at_cg: np.ndarray
    inertia: np.ndarray
    legacy_otter_inertia: bool
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


def preprocess_rigid_body(
    length: float,
    beam: float,
    water_density: float,
    hull_mass: float,
    payload_mass: float,
    hull_center_of_gravity: Sequence[float],
    payload_position: Sequence[float],
    added_mass_coefficients: Sequence[float],
    radii_of_gyration: Sequence[float],
    legacy_otter_inertia: bool = False,
) -> RigidBodyConstants:
    """Compute ``M_RB``, ``M_A`` and ``M`` for a hull with a point payload.

    Source names: ``l``, ``b``, ``rho``, ``m_hull``, ``m_payload``, ``r_hull``,
    ``r_payload``, ``coeff_added_mass``, ``R_456_scale``. Radii of gyration are
    scaled by ``[beam, length, length]`` and apply to the hull only.
    """
    positive_values = {
        "length": length,
        "beam": beam,
        "water_density": water_density,
        "hull_mass": hull_mass,
    }
    for name, value in positive_values.items():
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be a positive finite value")
    if not np.isfinite(payload_mass) or payload_mass < 0.0:
        raise ValueError("payload_mass must be a non-negative finite value")

    r_hull = _vector(hull_center_of_gravity, 3, "hull_center_of_gravity")
    r_payload = _vector(payload_position, 3, "payload_position")
    added_mass_scale = _vector(added_mass_coefficients, 6, "added_mass_coefficients")
    radius_scale = _vector(radii_of_gyration, 3, "radii_of_gyration")

    mass = hull_mass + payload_mass
    radii = radius_scale * np.array([beam, length, length], dtype=float)
    inertia_hull_cg = np.diag(hull_mass * np.square(radii))
    r_g = (hull_mass * r_hull + payload_mass * r_payload) / mass

    if legacy_otter_inertia:
        hull_offset, payload_offset = r_g, r_payload
    else:
        hull_offset, payload_offset = r_hull - r_g, r_payload - r_g
    inertia = (
        inertia_hull_cg
        - hull_mass * skew(hull_offset) @ skew(hull_offset)
        - payload_mass * skew(payload_offset) @ skew(payload_offset)
    )

    rigid_body_mass = rigid_body_mass_matrix(mass, inertia, r_g)
    added_mass = added_mass_matrix(
        added_mass_derivatives(
            hull_mass, length, water_density, inertia, added_mass_scale
        )
    )
    total_mass = rigid_body_mass + added_mass

    if not np.allclose(total_mass, total_mass.T):
        raise ValueError("total mass matrix must be symmetric")
    if not np.all(np.linalg.eigvalsh(total_mass) > 0.0):
        raise ValueError("total mass matrix must be positive definite")

    return RigidBodyConstants(
        mass=mass,
        center_of_gravity=r_g,
        inertia_hull_at_cg=inertia_hull_cg,
        inertia=inertia,
        legacy_otter_inertia=legacy_otter_inertia,
        rigid_body_mass_matrix=rigid_body_mass,
        added_mass_matrix=added_mass,
        total_mass_matrix=total_mass,
    )
