"""Rigid-body kinetics: ``M_RB`` (numpy) and ``C_RB(nu)`` (CasADi).

Ported from ``more_generic_models``
``dynamics/plant/matrices/rigid_body_kinetics_surface_vessel.py``:
``RigidBody6DOFSurfaceVessel.M_RB`` (lines 96-121), ``C_RB`` (123-138) and
``C_RB_runtime`` (147-170). Body frame z down, CO is the origin, ``r_g`` points
CO -> CG, ``H(r) = [[I, S(r)^T], [0, I]]`` (Fossen 2021, eq. 3.31).
"""

from typing import TYPE_CHECKING

import casadi as ca
import numpy as np

from .added_mass import added_mass_coriolis_casadi

if TYPE_CHECKING:
    from .constants import RigidBodyConstants


def skew(vector: np.ndarray) -> np.ndarray:
    """``S(v)`` with ``S(v) a = v x a`` (Fossen 2021, eq. 2.47)."""
    x, y, z = vector
    return np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]], dtype=float)


def _skew_casadi(vector: ca.SX) -> ca.SX:
    return ca.vertcat(
        ca.horzcat(0.0, -vector[2], vector[1]),
        ca.horzcat(vector[2], 0.0, -vector[0]),
        ca.horzcat(-vector[1], vector[0], 0.0),
    )


def transform_matrix(r_g: np.ndarray) -> np.ndarray:
    """``H(r_g)``, CG -> CO (``MatrixTransforms.H_matrix``, Fossen 2021 eq. 3.31)."""
    return np.block([[np.eye(3), skew(r_g).T], [np.zeros((3, 3)), np.eye(3)]])


def rigid_body_mass_matrix(
    mass: float, inertia: np.ndarray, r_g: np.ndarray
) -> np.ndarray:
    """``M_RB = H(r_g)^T diag(m I3, I) H(r_g)`` (source ``M_RB``, lines 96-121)."""
    mass_at_cg = np.block(
        [[mass * np.eye(3), np.zeros((3, 3))], [np.zeros((3, 3)), inertia]]
    )
    transform = transform_matrix(r_g)
    return transform.T @ mass_at_cg @ transform


def rigid_body_coriolis_casadi(
    mass: float, inertia: np.ndarray, r_g: np.ndarray, nu: ca.SX
) -> ca.SX:
    """``C_RB(nu) = H^T diag(m S(w), -S(I w)) H``, ``w = nu[3:6]``.

    Source ``C_RB_runtime`` (lines 147-170); the CG form moved to the CO, not
    the ``m2c(M_RB, nu)`` form (same ``C nu``, different matrix).
    """
    omega = nu[3:6]
    coriolis_at_cg = ca.SX.zeros(6, 6)
    coriolis_at_cg[0:3, 0:3] = mass * _skew_casadi(omega)
    coriolis_at_cg[3:6, 3:6] = -_skew_casadi(ca.DM(inertia) @ omega)
    transform = ca.DM(transform_matrix(r_g))
    return transform.T @ coriolis_at_cg @ transform


def rigid_body_casadi(constants: "RigidBodyConstants") -> ca.Function:
    """``(nu, nu_r) -> (M, C_RB(nu), C_A(nu_r))`` from preprocessed constants."""
    nu = ca.SX.sym("nu", 6)
    nu_r = ca.SX.sym("nu_r", 6)
    total_mass = ca.SX(ca.DM(constants.total_mass_matrix))
    c_rb = rigid_body_coriolis_casadi(
        constants.mass, constants.inertia, constants.center_of_gravity, nu
    )
    c_a = added_mass_coriolis_casadi(constants.added_mass_matrix, nu_r)
    return ca.Function(
        "rigid_body",
        [nu, nu_r],
        [total_mass, c_rb, c_a],
        ["nu", "nu_r"],
        ["M", "C_RB", "C_A"],
    )
