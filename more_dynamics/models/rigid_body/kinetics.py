"""Rigid-body + added-mass kinetics: ``M_RB`` (numpy), ``C_RB(nu)`` and
``C_A(nu_r)`` (CasADi).

Body frame z down, CO is the origin, ``r_g`` points CO -> CG, ``I`` is the
inertia about the CG. ``S`` and ``H`` come from ``more_transformations``:
numpy in the pre-processing, its CasADi mirror in the graph.

Equations (keys in References):

* ``H(r) = [[I3, S(r)^T], [0, I3]]`` (Fossen 2011, eq. 3.24, p. 49;
  MSS ``Hmtrx.m`` 16-18).
* ``M_RB = H(r_g)^T diag(m I3, I) H(r_g)`` (Fossen 2011, eqs. 3.26, 3.28,
  p. 50, and 3.44, p. 52; MSS ``rbody.m`` 36-44).
* Three ``C_RB`` parametrisations, one force ``C_RB(nu) nu``:

  - ``"co"``: ``C_RB = H^T diag(m S(w), -S(I w)) H`` (Fossen 2011, eqs. 3.21,
    p. 49, and 3.27, p. 50; MSS ``rbody.m`` 39-45, ``spheroid.m`` 46-52),
    depends on ``w = nu[3:6]`` only.
  - ``"book"``: ``C_RB = [[m S(w), -m S(w) S(r)], [m S(r) S(w), -S(I_O w)]]``,
    ``I_O = I - m S(r)^2`` (Fossen 2011, eqs. 3.29, p. 50, and 3.57, p. 55;
    parallel axes eq. 3.34, p. 50). ``co - book = m |r|^2 S(w)`` in the moment
    block only; the two give the same ``C_RB(nu) nu`` (Fossen 2011, eq. 3.30,
    p. 50).
  - ``"lagrangian"``: ``C = m2c(M_RB, nu) = [[0, -S(M11 v1 + M12 v2)],
    [-S(M11 v1 + M12 v2), -S(M21 v1 + M22 v2)]]`` (Fossen 2011, Theorem 3.2,
    eq. 3.46, p. 53, and eq. 3.55, p. 55; MSS ``m2c.m`` 35-48), depends on
    the linear velocity too.

* ``C_A(nu_r) = m2c(M_A, nu_r)`` with ``M_A`` symmetrised,
  ``M = (M + M^T)/2`` (Fossen 2011, eq. 6.43, p. 120, and eq. 6.40, p. 118;
  MSS ``m2c.m`` 33). For a diagonal ``M_A`` it is Fossen 2011 eq. 6.54,
  p. 121.
* Stabilised ``C_A``: the pitch-heave, pitch-surge, yaw-surge and yaw-sway
  entries set to zero in both triangles (MSS ``remus100.m`` 207-210). This is
  a modelling choice of that file, not physics: it removes the Munk moment
  of a slender body (Fossen 2011, eq. 6.52, p. 121), which MSS states in
  ``remus100.m`` 202-206.

Ported from the numpy source [MGM] ``dynamics/plant/matrices/``:
``rigid_body_kinetics.py`` (``get_M_RB_book`` 143-154, ``get_M_RB_co``
156-176, ``get_M_RB`` 179-206, ``get_C_RB_book`` 210-233, ``get_C_RB_co``
270-300, ``get_C_RB_lagrangian`` 303-352, ``get_C_RB`` 357-402),
``rigid_body_kinetics_spheroid.py`` (``M_RB`` 62-88, ``C_RB`` 95-127,
``C_RB_`` 136-159, ``C_RB_explicit`` 162-197, ``get_RB_matrices`` 203-207),
``rigid_body_kinetics_surface_vessel.py`` (``M_RB`` 97-121, ``C_RB`` 124-138,
``C_RB_runtime`` 147-170) and ``added_mass.py`` (``C_A_6dof`` 97-113,
``C_A_6dof_lagrangian`` 116-122, ``stabilize_C_A`` 156-171,
``get_C_RB_lagrangian`` 363-388). Two deviations from that source: the
Lagrangian ``C_A`` symmetrises ``M_A`` as MSS ``m2c.m`` 33 does (the source's
``get_C_RB_lagrangian`` does not; no difference for the symmetric ``M_A``
every form builds); the source's ``get_C_RB_book_corrected`` equals ``"co"``
and is not a separate form here.

References
----------
[Fossen 2011] Fossen, T. I. (2011). *Handbook of Marine Craft Hydrodynamics
    and Motion Control*, 1st ed. John Wiley & Sons, Chichester. Ch. 3,
    eqs. 3.21-3.57, pp. 49-55; Ch. 6, eqs. 6.38-6.54, pp. 118-121.
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2.
    https://github.com/cybergalactic/MSS, MIT licence, revision ``72656d1``:
    ``LIBRARY/kinematics/Hmtrx.m`` 16-18; ``LIBRARY/modeling/rbody.m`` 35-45,
    ``m2c.m`` 33-48, ``spheroid.m`` 45-52; ``CRAFT/AUV/models/remus100.m``
    202-210.
[MGM] Krizman, E. *more_generic_models*.
    https://github.com/MOREnvironment/more_generic_models (no licence file),
    revision ``524e336``: ``more_generic_models/dynamics/plant/matrices/``
    files and lines listed above.
"""

from typing import TYPE_CHECKING

import casadi as ca
import numpy as np
from more_transformations.matrix_transforms import MatrixTransforms
from more_transformations.more_casadi_transformations.matrix_transforms import (
    MatrixTransforms as CasadiMatrixTransforms,
)

from .added_mass import STABILIZED_ADDED_MASS_CORIOLIS_PAIRS

if TYPE_CHECKING:
    from .constants import RigidBodyConstants

CORIOLIS_FORMS = ("co", "book", "lagrangian")


def rigid_body_mass_matrix(
    mass: float, inertia: np.ndarray, r_g: np.ndarray
) -> np.ndarray:
    """``M_RB = H(r_g)^T diag(m I3, I) H(r_g)`` (Fossen 2011, eqs. 3.26, 3.28,
    p. 50; MSS ``rbody.m`` 36-44).

    Equals the source's ``get_M_RB_book`` ``[[m I3, -m S(r)], [m S(r),
    I - m S(r)^2]]`` (Fossen 2011, eq. 3.44, p. 52) to rounding; one function
    serves every form.
    """
    # M_RB^CG = diag(m I3, I) (Fossen 2011, eq. 3.21, p. 49; rbody.m 36-37)
    mass_at_cg = np.block(
        [[mass * np.eye(3), np.zeros((3, 3))], [np.zeros((3, 3)), inertia]]
    )
    transform = MatrixTransforms.H_matrix(r_g)  # (Fossen 2011, eq. 3.24, p. 49; Hmtrx.m 16-18)
    return transform.T @ mass_at_cg @ transform  # (Fossen 2011, eq. 3.26, p. 50; rbody.m 43-44)


def lagrangian_coriolis_casadi(mass_matrix: np.ndarray, nu: ca.SX) -> ca.SX:
    """``m2c(M, nu)`` 6-DOF: ``[[0, -S(M11 v1 + M12 v2)], [-S(M11 v1 + M12 v2),
    -S(M21 v1 + M22 v2)]]`` with ``M`` symmetrised and ``M21 = M12^T``
    (Fossen 2011, Theorem 3.2, eq. 3.46, p. 53; MSS ``m2c.m`` 33-48)."""
    # M = (M + M^T)/2 (Fossen 2011, eq. 6.40, p. 118; m2c.m 33)
    matrix = ca.DM(MatrixTransforms.symmetrize(mass_matrix))
    v1, v2 = nu[0:3], nu[3:6]
    linear = matrix[0:3, 0:3] @ v1 + matrix[0:3, 3:6] @ v2  # (Fossen 2011, eq. 3.49, p. 54; m2c.m 44)
    angular = matrix[0:3, 3:6].T @ v1 + matrix[3:6, 3:6] @ v2  # (Fossen 2011, eq. 3.50, p. 54; m2c.m 39, 45)
    s_linear = CasadiMatrixTransforms.skew(linear)
    # (Fossen 2011, eq. 3.46, p. 53; m2c.m 47-48)
    return ca.vertcat(
        ca.horzcat(ca.SX.zeros(3, 3), -s_linear),
        ca.horzcat(-s_linear, -CasadiMatrixTransforms.skew(angular)),
    )


def rigid_body_coriolis_casadi(
    mass: float, inertia: np.ndarray, r_g: np.ndarray, nu: ca.SX, coriolis: str
) -> ca.SX:
    """``C_RB(nu)`` in the chosen parametrisation (module docstring)."""
    omega = nu[3:6]
    s_omega = CasadiMatrixTransforms.skew(omega)
    if coriolis == "co":
        zeros = ca.SX.zeros(3, 3)
        # C_RB^CG = diag(m S(w), -S(I w)) (Fossen 2011, eq. 3.21, p. 49; rbody.m 39-40)
        coriolis_at_cg = ca.vertcat(
            ca.horzcat(mass * s_omega, zeros),
            ca.horzcat(zeros, -CasadiMatrixTransforms.skew(ca.DM(inertia) @ omega)),
        )
        transform = CasadiMatrixTransforms.H_matrix(r_g)
        return transform.T @ coriolis_at_cg @ transform  # (Fossen 2011, eq. 3.27, p. 50; rbody.m 43-45)
    if coriolis == "book":
        s_r = ca.DM(MatrixTransforms.skew(r_g))
        # I_O = I - m S(r)^2 (Fossen 2011, Theorem 3.1, eq. 3.34, p. 50)
        inertia_at_co = inertia - mass * MatrixTransforms.skew(r_g) @ MatrixTransforms.skew(r_g)
        # (Fossen 2011, eqs. 3.29, p. 50, and 3.57, p. 55)
        return ca.vertcat(
            ca.horzcat(mass * s_omega, -mass * s_omega @ s_r),
            ca.horzcat(
                mass * s_r @ s_omega,
                -CasadiMatrixTransforms.skew(ca.DM(inertia_at_co) @ omega),
            ),
        )
    if coriolis == "lagrangian":
        # m2c(M_RB, nu) (Fossen 2011, eq. 3.55, p. 55; m2c.m 35-48)
        return lagrangian_coriolis_casadi(rigid_body_mass_matrix(mass, inertia, r_g), nu)
    raise ValueError(f"coriolis must be one of {CORIOLIS_FORMS}, got {coriolis!r}")


def added_mass_coriolis_casadi(
    added_mass: np.ndarray, nu_r: ca.SX, stabilize: bool
) -> ca.SX:
    """``C_A(nu_r) = m2c(M_A, nu_r)`` (Fossen 2011, eq. 6.43, p. 120; MSS
    ``m2c.m`` 33-48); ``stabilize`` zeroes the pitch-heave, pitch-surge,
    yaw-surge and yaw-sway entries (MSS ``remus100.m`` 207-210).

    The stabilised form is a modelling choice, not physics: it removes the
    Munk moment of a slender body (Fossen 2011, eq. 6.52, p. 121), as MSS
    states in ``remus100.m`` 202-206.
    """
    coriolis = lagrangian_coriolis_casadi(added_mass, nu_r)  # (Fossen 2011, eq. 6.43, p. 120)
    if stabilize:  # (remus100.m 207-210)
        for row, column in STABILIZED_ADDED_MASS_CORIOLIS_PAIRS:
            coriolis[row, column] = 0.0
            coriolis[column, row] = 0.0
    return coriolis


def rigid_body_casadi(constants: "RigidBodyConstants") -> ca.Function:
    """``(nu, nu_r) -> (M, C_RB(nu), C_A(nu_r))`` from preprocessed constants;
    the forms are the ones stored on ``constants``."""
    nu = ca.SX.sym("nu", 6)
    nu_r = ca.SX.sym("nu_r", 6)
    total_mass = ca.SX(ca.DM(constants.total_mass_matrix))
    c_rb = rigid_body_coriolis_casadi(
        constants.mass,
        constants.inertia,
        constants.center_of_gravity,
        nu,
        constants.coriolis,
    )
    c_a = added_mass_coriolis_casadi(
        constants.added_mass_matrix, nu_r, constants.stabilize_added_mass_coriolis
    )
    return ca.Function(
        "rigid_body",
        [nu, nu_r],
        [total_mass, c_rb, c_a],
        ["nu", "nu_r"],
        ["M", "C_RB", "C_A"],
    )
