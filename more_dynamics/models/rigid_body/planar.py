"""3-DOF (surge, sway, yaw) forms as a reduction of the 6-DOF block (CasADi).

Models are 6-DOF and the 3-DOF form is obtained from them. One source of
truth: every 3-DOF matrix here is the surge, sway, yaw rows and
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

Historical deviation at MSS ``72656d1``; ``m2c.m`` agrees from ``cc07579``.
At ``72656d1`` the 3-DOF branch of ``m2c.m`` (lines 52-54) kept only
``M(1,1) u`` in ``p_x`` and ``M(2,2) v + M(2,3) r`` in ``p_y``. For ``M =
[[m11, 0, 0], [0, m22, m23], [0, m23, m33]]`` (the structure
``planar_added_mass_matrix`` builds) both are
equal entry by entry; for a mass matrix that couples surge with sway or yaw
(a CG off the centre line, ``M_RB(1,6) = -m y_g``) the reduction keeps the
terms that branch dropped, as the 6-DOF branch ``m2c.m`` 44-48 does. At
``cc07579`` the 3-DOF branch is ``p = M nu`` (``m2c.m`` 54-57), the same
``C3`` as this reduction.

References
----------
[Fossen 2011] Fossen, T. I. (2011). *Handbook of Marine Craft Hydrodynamics
    and Motion Control*, 1st ed. John Wiley & Sons, Chichester. Ch. 3,
    Theorem 3.2, eq. 3.46, p. 53; Ch. 6, Example 6.1, eqs. 6.50-6.51,
    p. 121.
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2.
    https://github.com/cybergalactic/MSS, MIT licence, revision ``72656d1``:
    ``LIBRARY/modeling/m2c.m`` 44-54; revision ``cc07579`` (release 2.0.2
    with the fixes of 2026-10-07): ``LIBRARY/modeling/m2c.m`` 44-57.

Author:    Enio Krizman
Date:      2026-10-07
"""

import casadi as ca
from more_transformations.more_casadi_transformations import Parameter

from .kinetics import lagrangian_coriolis_casadi, rigid_body_casadi
from .mass_properties import rigid_body_parameters

PLANAR_DOFS = (0, 1, 5)  # surge, sway, yaw in nu = [u, v, w, p, q, r]

_PLANAR_MASS_MATRIX = Parameter(
    "mass_matrix", (3, 3), "kg, kg*m, kg*m^2",
    "3x3 inertia matrix in surge, sway, yaw (rigid body, added mass or their sum)",
)


def _embed_velocity(nu):
    """``[u, v, r] -> [u, v, 0, 0, 0, r]``."""
    return ca.vertcat(nu[0], nu[1], 0.0, 0.0, 0.0, nu[2])


def planar_added_mass_matrix(X_du, Y_dv, Y_dr, N_dr):
    """``-[[X_du, 0, 0], [0, Y_dv, Y_dr], [0, Y_dr, N_dr]]``: the surge, sway,
    yaw rows and columns of the 6-DOF ``M_A`` with diagonal ``-[X_du, Y_dv,
    ., ., ., N_dr]`` and the sway-yaw coupling ``M_A(2,6) = M_A(6,2) =
    -Y_dr`` (Fossen 2011, eq. 6.50, p. 121). Symbols or
    numbers in, a 3x3 CasADi matrix out."""
    return -ca.vertcat(  # (Fossen 2011, eq. 6.50, p. 121)
        ca.horzcat(X_du, 0.0, 0.0),
        ca.horzcat(0.0, Y_dv, Y_dr),
        ca.horzcat(0.0, Y_dr, N_dr),
    )


def planar_coriolis_parameters():
    """The declared parameter of ``planar_coriolis_casadi``: ``mass_matrix``
    (3x3)."""
    return (_PLANAR_MASS_MATRIX,)


def planar_coriolis_casadi():
    """``(nu = [u, v, r], mass_matrix) -> C = m2c(M, nu)`` 3-DOF, as the
    reduction of the 6-DOF ``m2c`` of ``mass_matrix`` placed in surge, sway,
    yaw (Fossen 2011, eq. 3.46, p. 53, reduced; equals ``m2c.m`` 52-54 only
    without surge coupling, module docstring)."""
    nu = ca.SX.sym("nu", 3)
    mass_matrix = ca.SX.sym("mass_matrix", 3, 3)
    index = list(PLANAR_DOFS)
    full = ca.SX.zeros(6, 6)
    full[index, index] = mass_matrix  # the 3x3 matrix placed in surge, sway, yaw (reduction rule, module docstring)
    coriolis = lagrangian_coriolis_casadi(full, _embed_velocity(nu))  # (Fossen 2011, eq. 3.46, p. 53; m2c.m 44-48)
    return ca.Function(
        "planar_coriolis", [nu, mass_matrix], [coriolis[index, index]], ["nu", "mass_matrix"], ["C"]
    )


def planar_casadi(
    *,
    mass_properties="hull_with_payload",
    coriolis="co",
    stabilize_added_mass_coriolis=False,
):
    """``(nu, nu_r, <each declared parameter by name>) -> (M, C_RB(nu),
    C_A(nu_r))``, velocities 3x1 ``[u, v, r]``, matrices 3x3: the surge,
    sway, yaw reduction of ``rigid_body_casadi`` built with the same
    keywords, evaluated at ``[u, v, 0, 0, 0, r]``. The parameters are those
    of ``rigid_body_parameters(mass_properties)``."""
    full = rigid_body_casadi(
        mass_properties=mass_properties,
        coriolis=coriolis,
        stabilize_added_mass_coriolis=stabilize_added_mass_coriolis,
    )
    names = [d.name for d in rigid_body_parameters(mass_properties)]
    p = {name: full.sx_in(name) for name in names}
    nu = ca.SX.sym("nu", 3)
    nu_r = ca.SX.sym("nu_r", 3)
    out = full(nu=_embed_velocity(nu), nu_r=_embed_velocity(nu_r), **p)
    index = list(PLANAR_DOFS)
    return ca.Function(
        "planar_rigid_body",
        [nu, nu_r, *[p[name] for name in names]],
        # reduction rule X3 = X[PLANAR_DOFS, PLANAR_DOFS] (module docstring; Fossen 2011, Example 6.1, p. 121)
        [out["M"][index, index], out["C_RB"][index, index], out["C_A"][index, index]],
        ["nu", "nu_r", *names],
        ["M", "C_RB", "C_A"],
    )
