"""Added mass: ``M_A`` (numpy) and ``C_A(nu_r)`` (CasADi).

Ported from ``more_generic_models`` ``dynamics/plant/matrices/added_mass.py``:
``added_mass_surge`` (lines 313-335), ``get_added_mass_derivates_catamaran``
(269-289), ``M_A_6dof`` (24-41), ``C_A_6dof_lagrangian`` (116-122) and
``get_C_RB_lagrangian`` (362-388). Diagonal ``M_A``, no Munk-moment
cancellation (``stabilize=False``).
"""

import casadi as ca
import numpy as np


def surge_added_mass(mass: float, length: float, water_density: float) -> float:
    """``A11 = 2.7 rho nabla^(5/3) / L^2``, ``nabla = m / rho`` (``added_mass_surge``)."""
    displaced_volume = mass / water_density
    return 2.7 * water_density * displaced_volume ** (5.0 / 3.0) / length**2


def added_mass_derivatives(
    mass: float,
    length: float,
    water_density: float,
    inertia: np.ndarray,
    coefficients: np.ndarray,
) -> np.ndarray:
    """``[X_du, Y_dv, Z_dw, K_dp, M_dq, N_dr]`` (``get_added_mass_derivates_catamaran``).

    ``coefficients * [A11, m, m, I[0,0], I[1,1], I[2,2]]``; ``mass`` is the hull
    mass and ``inertia`` the one the rigid body uses.
    """
    return coefficients * np.array(
        [
            surge_added_mass(mass, length, water_density),
            mass,
            mass,
            inertia[0, 0],
            inertia[1, 1],
            inertia[2, 2],
        ]
    )


def added_mass_matrix(derivatives: np.ndarray) -> np.ndarray:
    """``M_A = -diag(derivatives)`` (``M_A_6dof``)."""
    return np.diag(-np.asarray(derivatives, dtype=float))


def _skew_casadi(vector: ca.SX) -> ca.SX:
    return ca.vertcat(
        ca.horzcat(0.0, -vector[2], vector[1]),
        ca.horzcat(vector[2], 0.0, -vector[0]),
        ca.horzcat(-vector[1], vector[0], 0.0),
    )


def added_mass_coriolis_casadi(added_mass: np.ndarray, nu_r: ca.SX) -> ca.SX:
    """``C_A(nu_r)`` in the ``m2c`` form (Fossen 2021, Theorem 3.2).

    Source ``get_C_RB_lagrangian``: ``[[0, -S(M11 v1 + M12 v2)],
    [-S(M11 v1 + M12 v2), -S(M21 v1 + M22 v2)]]``; ``M`` is not symmetrised.
    """
    matrix = ca.DM(added_mass)
    v1, v2 = nu_r[0:3], nu_r[3:6]
    linear = matrix[0:3, 0:3] @ v1 + matrix[0:3, 3:6] @ v2
    angular = matrix[3:6, 0:3] @ v1 + matrix[3:6, 3:6] @ v2
    coriolis = ca.SX.zeros(6, 6)
    coriolis[0:3, 3:6] = -_skew_casadi(linear)
    coriolis[3:6, 0:3] = -_skew_casadi(linear)
    coriolis[3:6, 3:6] = -_skew_casadi(angular)
    return coriolis
