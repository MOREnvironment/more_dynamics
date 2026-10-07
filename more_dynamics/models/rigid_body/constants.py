"""Preprocessed constants of the rigid-body + added-mass block (numpy).

A form is chosen by explicit arguments of ``preprocess_rigid_body``, never by
a vehicle name: ``mass_properties`` says how ``m``, ``I`` (about the CG),
``r_g`` and ``M_A`` are obtained, ``coriolis`` which ``C_RB``
parametrisation the graph builds (``kinetics.py``),
``stabilize_added_mass_coriolis`` whether ``C_A`` drops its Munk entries.

Equations per form (keys in References):

* ``"hull_with_payload"``: ``m = m_h + m_p``, ``r_g = (m_h r_h + m_p r_p) /
  m``, ``I = I_h - m_h S(r_h - r_g)^2 - m_p S(r_p - r_g)^2`` about the
  combined CG, ``I_h = m_h diag((R_s [B, L, L])^2)`` (parallel axes, Fossen
  2011, Theorem 3.1, eq. 3.34, p. 50; MSS ``otter.m`` 122-128).
* ``"displacement_hull"``: ``nabla = Cb L B T``, ``m = rho nabla``,
  ``I = m diag((R_s [B, L, L])^2)`` about the CG, wetted surface
  ``L B + 2 T B`` (the numpy source's ``get_params_hull``; the radii-of-
  gyration form as MSS ``rbody.m`` 35). Rotational added mass scales ``I``
  about the CG (MSS ``otter.m`` 155-157).
* ``"spheroid"``: ``m = rho 4/3 pi a b^2``, ``I = diag(2/5 m b^2,
  1/5 m (a^2 + b^2), 1/5 m (a^2 + b^2))`` (MSS ``spheroid.m`` 35-42,
  ``imlay61.m`` 31-34); Lamb added mass of the displaced fluid
  (``added_mass.py``). Two densities: MSS uses 1025 for the body
  (``spheroid.m`` 35) and 1026 for the displaced fluid (``imlay61.m`` 31),
  the numpy source one ``rho`` for both.
* ``"explicit"``: given ``m``, diagonal ``I`` about the CG, ``r_g = r_cg -
  origin`` (the CO -> CG vector of Fossen 2011, eq. 3.24, p. 49),
  ``M_A = -diag(derivatives)`` (Fossen 2011, eq. 6.53, p. 121).
* All forms: ``M = M_RB + M_A`` (Fossen 2011, eq. 6.48, p. 120), symmetric
  and positive definite (Fossen 2011, Property 3.1, eq. 3.43, p. 52, and
  Property 6.1, p. 118).

Ported from the numpy source [MGM] ``dynamics/plant/matrices/``:

* ``"hull_with_payload"``: ``rigid_body_kinetics_surface_vessel.py::
  get_params_catamaran`` (17-51) and ``plant/asv_catamaran/asv_catamaran.py::
  compute_constant_values``. Deviation from that source: the inertia is taken
  about the combined CG as MSS ``otter.m`` 122-128 does since MSS commit
  ``880b2ef`` (2026-04-20); the source and its frozen reference keep the older inertia,
  shifted twice by the parallel-axis theorem.
* ``"displacement_hull"``: ``get_params_hull`` (54-94). Deviation from that
  source: it returns ``I_CG - m S(r)^2`` (about the CO) and its
  ``C_RB_runtime`` shifts it a second time with ``H``; that double shift is
  not ported (the Newton-Euler wrench, Fossen 2011, eqs. 3.33 and 3.40,
  pp. 50-51, agrees with ``rbody.m``, not with the source).
* ``"spheroid"``: ``rigid_body_kinetics_spheroid.py::get_I_gb`` (31-59) and
  the Lamb form of ``added_mass.py``.
* ``"explicit"``: ``rigid_body_kinetics.py::RigidBody6DOF`` (``__init__``
  23-57, ``get_r_bg`` 72-77) with ``M_A_6dof``. Deviations from that source:
  its "CG not below CO" warning and ``normalize_zero`` are not ported
  (owner's decision of 2026-10-06); the inertia must be exactly diagonal,
  where the source accepts off-diagonal entries below 1e-9.

References
----------
[Fossen 2011] Fossen, T. I. (2011). *Handbook of Marine Craft Hydrodynamics
    and Motion Control*, 1st ed. John Wiley & Sons, Chichester. Ch. 3,
    eq. 3.24, p. 49; eqs. 3.33-3.44, pp. 50-52; Ch. 6, Property 6.1 and
    eqs. 6.48, 6.53, pp. 118-121.
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2.
    https://github.com/cybergalactic/MSS, MIT licence, revision ``72656d1``:
    ``CRAFT/USV/models/otter.m`` 122-128, 155-157;
    ``LIBRARY/modeling/rbody.m`` 35, ``spheroid.m`` 35-42, ``imlay61.m``
    31-34.
[MGM] Krizman, E. *more_generic_models*.
    https://github.com/MOREnvironment/more_generic_models (no licence file),
    revision ``524e336``: ``more_generic_models/dynamics/plant/`` files and
    lines listed above.
"""

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
from more_transformations.matrix_transforms import MatrixTransforms

from .added_mass import (
    added_mass_matrix,
    lamb_k_factors,
    scaled_added_mass_derivatives,
    spheroid_added_mass_derivatives,
)
from .kinetics import CORIOLIS_FORMS, rigid_body_mass_matrix

MASS_PROPERTIES_FORMS = ("hull_with_payload", "displacement_hull", "spheroid", "explicit")


@dataclass(frozen=True)
class RigidBodyConstants:
    """Constants computed before building the CasADi graph.

    ``inertia`` is about the CG, ``center_of_gravity`` is CO -> CG.
    ``lamb_k_factors`` is set by ``"spheroid"`` only; ``displaced_volume``
    and ``wetted_surface`` by ``"displacement_hull"`` only (else ``None``).
    """

    mass_properties: str
    coriolis: str
    stabilize_added_mass_coriolis: bool
    mass: float
    center_of_gravity: np.ndarray
    inertia: np.ndarray
    added_mass_derivatives: np.ndarray
    rigid_body_mass_matrix: np.ndarray
    added_mass_matrix: np.ndarray
    total_mass_matrix: np.ndarray
    lamb_k_factors: Optional[np.ndarray]
    displaced_volume: Optional[float]
    wetted_surface: Optional[float]


@dataclass(frozen=True)
class _MassProperties:
    mass: float
    center_of_gravity: np.ndarray
    inertia: np.ndarray
    added_mass_derivatives: np.ndarray
    lamb_k_factors: Optional[np.ndarray] = None
    displaced_volume: Optional[float] = None
    wetted_surface: Optional[float] = None


def _vector(values: Sequence[float], size: int, name: str) -> np.ndarray:
    vector = np.asarray(values, dtype=float)
    if vector.shape != (size,):
        raise ValueError(f"{name} must contain exactly {size} values")
    if not np.all(np.isfinite(vector)):
        raise ValueError(f"{name} must contain only finite values")
    return vector


def _positive(**values: float) -> None:
    for name, value in values.items():
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be a positive finite value")


def _hull_with_payload(
    *,
    length: float,
    beam: float,
    water_density: float,
    hull_mass: float,
    payload_mass: float,
    hull_center_of_gravity: Sequence[float],
    payload_position: Sequence[float],
    added_mass_coefficients: Sequence[float],
    radii_of_gyration: Sequence[float],
) -> _MassProperties:
    """Hull with a point payload. Source names: ``l``, ``b``, ``rho``,
    ``m_hull``, ``m_payload``, ``r_hull``, ``r_payload``, ``coeff_added_mass``,
    ``R_456_scale``. Radii of gyration are scaled by ``[beam, length,
    length]`` and apply to the hull only; added mass uses the hull mass."""
    _positive(length=length, beam=beam, water_density=water_density, hull_mass=hull_mass)
    if not np.isfinite(payload_mass) or payload_mass < 0.0:
        raise ValueError("payload_mass must be a non-negative finite value")
    r_hull = _vector(hull_center_of_gravity, 3, "hull_center_of_gravity")
    r_payload = _vector(payload_position, 3, "payload_position")
    coefficients = _vector(added_mass_coefficients, 6, "added_mass_coefficients")
    radius_scale = _vector(radii_of_gyration, 3, "radii_of_gyration")

    mass = hull_mass + payload_mass
    radii = radius_scale * np.array([beam, length, length], dtype=float)
    inertia_hull_cg = np.diag(hull_mass * np.square(radii))  # (otter.m 122)
    r_g = (hull_mass * r_hull + payload_mass * r_payload) / mass  # (otter.m 124)
    s_hull = MatrixTransforms.skew(r_hull - r_g)  # (otter.m 125)
    s_payload = MatrixTransforms.skew(r_payload - r_g)  # (otter.m 126)
    # parallel axes to the combined CG (Fossen 2011, eq. 3.34, p. 50; otter.m 127)
    inertia = inertia_hull_cg - hull_mass * s_hull @ s_hull - payload_mass * s_payload @ s_payload
    derivatives = scaled_added_mass_derivatives(
        hull_mass, length, water_density, inertia, coefficients
    )
    return _MassProperties(mass, r_g, inertia, derivatives)


def _displacement_hull(
    *,
    water_density: float,
    length: float,
    beam: float,
    draft: float,
    block_coefficient: float,
    radii_of_gyration: Sequence[float],
    center_of_gravity: Sequence[float],
    added_mass_coefficients: Sequence[float],
) -> _MassProperties:
    """Neutrally buoyant hull: ``nabla = Cb L B T``, ``m = rho nabla``,
    ``I = m diag((R_s [B, L, L])^2)`` about the CG, wetted surface
    ``L B + 2 T B`` (source ``get_params_hull``; used by the surge damping)."""
    _positive(
        water_density=water_density,
        length=length,
        beam=beam,
        draft=draft,
        block_coefficient=block_coefficient,
    )
    radius_scale = _vector(radii_of_gyration, 3, "radii_of_gyration")
    r_g = _vector(center_of_gravity, 3, "center_of_gravity")
    coefficients = _vector(added_mass_coefficients, 6, "added_mass_coefficients")

    displaced_volume = block_coefficient * length * beam * draft  # (source get_params_hull)
    mass = water_density * displaced_volume  # neutrally buoyant
    radii = radius_scale * np.array([beam, length, length], dtype=float)
    inertia = np.diag(mass * np.square(radii))  # about the CG (rbody.m 35)
    derivatives = scaled_added_mass_derivatives(
        mass, length, water_density, inertia, coefficients
    )
    return _MassProperties(
        mass,
        r_g,
        inertia,
        derivatives,
        displaced_volume=displaced_volume,
        wetted_surface=length * beam + 2.0 * draft * beam,  # (source get_params_hull)
    )


def _spheroid_mass_and_inertia(density: float, a: float, b: float):
    """``m = rho 4/3 pi a b^2``, ``I = diag(2/5 m b^2, 1/5 m (a^2+b^2),
    1/5 m (a^2+b^2))`` (MSS ``spheroid.m`` 35-42, ``imlay61.m`` 31-34)."""
    mass = density * 4.0 / 3.0 * np.pi * a * b**2  # (spheroid.m 36; imlay61.m 32)
    lateral = 0.2 * mass * (a**2 + b**2)  # (spheroid.m 40-41; imlay61.m 34)
    return mass, np.diag([0.4 * mass * b**2, lateral, lateral])  # (spheroid.m 39, 42; imlay61.m 33)


def _spheroid(
    *,
    semi_major_axis: float,
    semi_minor_axis: float,
    body_density: float,
    water_density: float,
    roll_added_inertia_ratio: float,
    center_of_gravity: Sequence[float],
) -> _MassProperties:
    """Prolate spheroid ``a > b``: body from ``body_density``, Lamb added
    mass from the displaced fluid (``water_density``); ``M_A(4,4) = r44 I_x``
    of the fluid (MSS ``imlay61.m`` 47)."""
    _positive(body_density=body_density, water_density=water_density)
    if not np.isfinite(roll_added_inertia_ratio) or roll_added_inertia_ratio < 0.0:
        raise ValueError("roll_added_inertia_ratio must be a non-negative finite value")
    r_g = _vector(center_of_gravity, 3, "center_of_gravity")
    k_factors = lamb_k_factors(semi_major_axis, semi_minor_axis)
    a, b = float(semi_major_axis), float(semi_minor_axis)
    mass, inertia = _spheroid_mass_and_inertia(body_density, a, b)
    fluid_mass, fluid_inertia = _spheroid_mass_and_inertia(water_density, a, b)
    derivatives = spheroid_added_mass_derivatives(
        k_factors, fluid_mass, fluid_inertia, roll_added_inertia_ratio
    )
    return _MassProperties(mass, r_g, inertia, derivatives, lamb_k_factors=k_factors)


def _explicit(
    *,
    mass: float,
    inertia: Sequence[Sequence[float]],
    center_of_gravity: Sequence[float],
    added_mass_derivatives: Sequence[float],
    origin: Sequence[float] = (0.0, 0.0, 0.0),
) -> _MassProperties:
    """Given ``m``, diagonal ``I`` about the CG, the CG and the derivative set
    ``[X_du, Y_dv, Z_dw, K_dp, M_dq, N_dr]``; ``r_g = center_of_gravity -
    origin`` (source ``get_r_bg``)."""
    _positive(mass=mass)
    inertia = np.asarray(inertia, dtype=float)
    if inertia.shape != (3, 3) or not np.all(np.isfinite(inertia)):
        raise ValueError("inertia must be a finite 3x3 matrix")
    if np.any(inertia != np.diag(np.diag(inertia))):
        raise ValueError("inertia must be diagonal (principal axes along the body axes)")
    if np.any(np.diag(inertia) <= 0.0):
        raise ValueError("inertia must have positive diagonal entries")
    # r_g = r_cg - r_co (source get_r_bg 72-77)
    r_g = _vector(center_of_gravity, 3, "center_of_gravity") - _vector(origin, 3, "origin")
    derivatives = _vector(added_mass_derivatives, 6, "added_mass_derivatives")
    return _MassProperties(float(mass), r_g, inertia, derivatives)


_FORMS = {
    "hull_with_payload": _hull_with_payload,
    "displacement_hull": _displacement_hull,
    "spheroid": _spheroid,
    "explicit": _explicit,
}


def preprocess_rigid_body(
    *,
    mass_properties: str = "hull_with_payload",
    coriolis: str = "co",
    stabilize_added_mass_coriolis: bool = False,
    **form_keywords,
) -> RigidBodyConstants:
    """Compute ``M_RB``, ``M_A`` and ``M`` of the chosen form.

    ``form_keywords`` are the keywords of the chosen ``mass_properties``
    form (module docstring; the contract is the header of
    ``tests/rigid_body/test_rigid_body_forms.py``). A missing or foreign
    keyword raises ``TypeError``, an unknown form ``ValueError``.
    """
    if mass_properties not in _FORMS:
        raise ValueError(
            f"mass_properties must be one of {MASS_PROPERTIES_FORMS}, got {mass_properties!r}"
        )
    if coriolis not in CORIOLIS_FORMS:
        raise ValueError(f"coriolis must be one of {CORIOLIS_FORMS}, got {coriolis!r}")
    if not isinstance(stabilize_added_mass_coriolis, bool):
        raise ValueError("stabilize_added_mass_coriolis must be True or False")

    properties = _FORMS[mass_properties](**form_keywords)
    rigid_body_mass = rigid_body_mass_matrix(
        properties.mass, properties.inertia, properties.center_of_gravity
    )
    added_mass = added_mass_matrix(properties.added_mass_derivatives)
    total_mass = rigid_body_mass + added_mass  # (Fossen 2011, eq. 6.48, p. 120)

    # M = M^T > 0 (Fossen 2011, Property 3.1, eq. 3.43, p. 52; Property 6.1, p. 118)
    if not np.allclose(total_mass, total_mass.T):
        raise ValueError("total mass matrix must be symmetric")
    if not np.all(np.linalg.eigvalsh(total_mass) > 0.0):
        raise ValueError("total mass matrix must be positive definite")

    return RigidBodyConstants(
        mass_properties=mass_properties,
        coriolis=coriolis,
        stabilize_added_mass_coriolis=stabilize_added_mass_coriolis,
        mass=properties.mass,
        center_of_gravity=properties.center_of_gravity,
        inertia=properties.inertia,
        added_mass_derivatives=properties.added_mass_derivatives,
        rigid_body_mass_matrix=rigid_body_mass,
        added_mass_matrix=added_mass,
        total_mass_matrix=total_mass,
        lamb_k_factors=properties.lamb_k_factors,
        displaced_volume=properties.displaced_volume,
        wetted_surface=properties.wetted_surface,
    )
