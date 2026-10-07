"""3-DOF (surge, sway, yaw) forms as a reduction of the 6-DOF block.

The owner's decision of 2026-10-06: models are 6-DOF and the 3-DOF form is
obtained from them (*"we should have both, but usually our models should be
6 dof high fidelity and from 6 dof we should somehow scale to 3"*). One
source of truth: every 3-DOF matrix here is the surge, sway, yaw rows and
columns (``PLANAR_DOFS``) of the 6-DOF matrix built by ``kinetics.py`` with
heave, roll and pitch velocities set to zero. ``C(nu) nu`` of the reduction
is therefore exactly the surge, sway, yaw part of the 6-DOF force for planar
motion; heave, roll and pitch are neglected, as in every 3-DOF model.

Equations (keys in References):

* Reduction: ``X3 = X[PLANAR_DOFS, PLANAR_DOFS]`` of a 6x6 matrix, at
  ``nu = [u, v, 0, 0, 0, r]`` (a construction of this module; the selection
  of surge, sway and yaw is that of Fossen 2011, Example 6.1, p. 121).
* ``M_A3 = -[[X_du, 0, 0], [0, Y_dv, Y_dr], [0, Y_dr, N_dr]]`` (Fossen 2011,
  eq. 6.50, p. 121).
* ``C3(nu) = m2c(M, [u, v, 0, 0, 0, r])[PLANAR_DOFS, PLANAR_DOFS]``; with
  ``p = M3 [u, v, r]`` this is ``[[0, 0, -p_y], [0, 0, p_x], [p_y, -p_x, 0]]``
  (Fossen 2011, Theorem 3.2, eq. 3.46, p. 53, reduced; for ``M_A3`` it is
  eq. 6.51, p. 121).

Deviation from MSS ``m2c.m`` 52-54 (the 3-DOF branch): that branch keeps only
``M(1,1) u`` in ``p_x`` and ``M(2,2) v + M(2,3) r`` in ``p_y``. For ``M =
[[m11, 0, 0], [0, m22, m23], [0, m23, m33]]`` (the structure
``planar_added_mass_matrix`` builds, and the source's ``M_A_3dof``) both are
equal entry by entry; for a mass matrix that couples surge with sway or yaw
(a CG off the centre line, ``M_RB(1,6) = -m y_g``) the reduction keeps the
terms the 3-DOF branch drops, as the 6-DOF branch ``m2c.m`` 44-48 does.

Ported from the numpy source [MGM] ``dynamics/plant/matrices/``:
``added_mass.py::M_A_3dof`` (81-88), ``C_A_3dof`` (126-152) and the 3-DOF
branch of ``rigid_body_kinetics.py::get_C_RB_lagrangian`` (342-348).

References
----------
[Fossen 2011] Fossen, T. I. (2011). *Handbook of Marine Craft Hydrodynamics
    and Motion Control*, 1st ed. John Wiley & Sons, Chichester. Ch. 3,
    Theorem 3.2, eq. 3.46, p. 53; Ch. 6, Example 6.1, eqs. 6.50-6.51,
    p. 121.
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2.
    https://github.com/cybergalactic/MSS, MIT licence, revision ``72656d1``:
    ``LIBRARY/modeling/m2c.m`` 44-54.
[MGM] Krizman, E. *more_generic_models*.
    https://github.com/MOREnvironment/more_generic_models (no licence file),
    revision ``524e336``: ``more_generic_models/dynamics/plant/matrices/``
    files and lines listed above.
"""

from typing import TYPE_CHECKING

import casadi as ca
import numpy as np

from .kinetics import lagrangian_coriolis_casadi, rigid_body_casadi

PLANAR_DOFS = (0, 1, 5)  # surge, sway, yaw in nu = [u, v, w, p, q, r]

if TYPE_CHECKING:
    from .constants import RigidBodyConstants


def planar_reduction(matrix: np.ndarray) -> np.ndarray:
    """Surge, sway, yaw rows and columns of a 6x6 matrix."""
    matrix = np.asarray(matrix, dtype=float)
    if matrix.shape != (6, 6):
        raise ValueError("matrix must be 6x6")
    return matrix[np.ix_(PLANAR_DOFS, PLANAR_DOFS)]  # reduction rule (module docstring)


def _embed_matrix(matrix: np.ndarray) -> np.ndarray:
    """A 3x3 surge-sway-yaw matrix placed in a 6x6 zero matrix."""
    matrix = np.asarray(matrix, dtype=float)
    if matrix.shape != (3, 3) or not np.all(np.isfinite(matrix)):
        raise ValueError("mass_matrix must be a finite 3x3 matrix")
    full = np.zeros((6, 6))
    full[np.ix_(PLANAR_DOFS, PLANAR_DOFS)] = matrix
    return full


def _embed_velocity(nu: ca.SX) -> ca.SX:
    """``[u, v, r] -> [u, v, 0, 0, 0, r]``."""
    return ca.vertcat(nu[0], nu[1], 0.0, 0.0, 0.0, nu[2])


def planar_added_mass_matrix(X_du: float, Y_dv: float, Y_dr: float, N_dr: float) -> np.ndarray:
    """``-[[X_du, 0, 0], [0, Y_dv, Y_dr], [0, Y_dr, N_dr]]``: the reduction of
    the 6-DOF ``M_A`` with diagonal ``-[X_du, Y_dv, ., ., ., N_dr]`` and the
    sway-yaw coupling ``M_A(2,6) = M_A(6,2) = -Y_dr`` (Fossen 2011, eq. 6.50,
    p. 121; source ``M_A_3dof``)."""
    full = np.zeros((6, 6))  # 6-DOF M_A placed, then reduced (Fossen 2011, eq. 6.50, p. 121)
    full[0, 0] = -X_du
    full[1, 1] = -Y_dv
    full[5, 5] = -N_dr
    full[1, 5] = full[5, 1] = -Y_dr
    return planar_reduction(full)


def planar_coriolis_casadi(mass_matrix: np.ndarray) -> ca.Function:
    """``nu = [u, v, r] -> C = m2c(M, nu)`` 3-DOF, as the reduction of the
    6-DOF ``m2c`` of ``mass_matrix`` placed in surge, sway, yaw (Fossen 2011,
    eq. 3.46, p. 53, reduced; equals ``m2c.m`` 52-54 only without surge
    coupling, module docstring)."""
    full = _embed_matrix(mass_matrix)
    nu = ca.SX.sym("nu", 3)
    coriolis = lagrangian_coriolis_casadi(full, _embed_velocity(nu))  # (Fossen 2011, eq. 3.46, p. 53; m2c.m 44-48)
    index = list(PLANAR_DOFS)
    return ca.Function("planar_coriolis", [nu], [coriolis[index, index]], ["nu"], ["C"])


def planar_casadi(constants: "RigidBodyConstants") -> ca.Function:
    """``(nu, nu_r) -> (M, C_RB(nu), C_A(nu_r))``, 3x1 in, 3x3 out: the
    surge, sway, yaw reduction of ``rigid_body_casadi(constants)`` (the same
    constants object, the forms stored on it)."""
    full = rigid_body_casadi(constants)
    nu = ca.SX.sym("nu", 3)
    nu_r = ca.SX.sym("nu_r", 3)
    out = full(nu=_embed_velocity(nu), nu_r=_embed_velocity(nu_r))
    index = list(PLANAR_DOFS)
    return ca.Function(
        "planar_rigid_body",
        [nu, nu_r],
        [out["M"][index, index], out["C_RB"][index, index], out["C_A"][index, index]],
        ["nu", "nu_r"],
        ["M", "C_RB", "C_A"],
    )
