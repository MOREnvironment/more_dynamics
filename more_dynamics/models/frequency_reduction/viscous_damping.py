"""Viscous damping that completes the reduced potential damping of a fully submerged hull, and the move of a 6 x 6
matrix between the body origin and the table reference point, as CasADi blocks.

The potential damping of a deeply submerged hull is small; the total damping of the unrestrained modes is set by a time
constant ``T_i`` (heave included: a submerged hull has no heave stiffness), the restored modes (roll, pitch) by a
damping-ratio increment ``delta_zeta_i``::

    Bv_ii = M_ii / T_i - B_eq,ii           (i = surge, sway, heave, yaw)
    Bv_ii = 2 delta_zeta_i sqrt(M_ii G_ii)  (i = roll, pitch)                          (Fossen 2025, the power-based method)

``M = M_RB + A_eq`` and ``G`` are about the table reference point, as the tables. For the scalar mode ``M_ii nu_dot +
D_ii nu = 0`` the decay time is ``M_ii / D_ii = T_i``; for ``M x_ddot + 2 zeta sqrt(M G) x_dot + G x = 0`` the added
damping ratio is ``delta_zeta`` (the mass-damper-spring form, Fossen 2011, eqs. 12.57-12.59, p. 365: d = 2 zeta sqrt(m k)). A time constant that
asks for a negative viscous damping is refused at the plugin boundary (``negative_time_constant_margin`` is the value
the plugin checks); the block itself returns ``max(0, M_ii/T_i - B_eq,ii)``.

Reference point (Fossen 2011, eq. 3.24, p. 49, and the congruence of eq. 3.27, p. 50)::

    X_CO = H(r)^T X_R H(r),        X_R = H(-r)^T X_CO H(-r)        (H(a) H(b) = H(a + b))

``r`` from the body origin CO to the point R, BODY FRD, metres; it holds for a mass matrix and a damping matrix alike.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics and Motion Control. Wiley. Eq. 3.24, p. 49;
    eq. 3.27, p. 50; eqs. 12.57-12.59, p. 365.
[Fossen 2025] Fossen, T. I. (2025). Maneuvering coefficient estimation from frequency-dependent added mass and damping:
    a power-based approach. Ocean Engineering 341:122494. (Not read: not on disk.)

Author:    Enio Krizman
Date:      2026-11-05
"""
import casadi as ca

from more_transformations.more_casadi_transformations import MatrixTransforms, Parameter

from more_dynamics.models.shared.wiring import function_from

TIME_CONSTANT_DOFS = (0, 1, 2, 5)    # surge, sway, heave, yaw
RATIO_DOFS = (3, 4)                  # roll, pitch
DAMPING_TOLERANCE = 1e-10            # relative: a margin below -tol max(1, |M/T|, |B_eq|) is refused

SUBMERGED_VISCOUS_DAMPING_PARAMETERS = (
    Parameter("T_1236", (4, 1), "s", "time constants of the total damping [surge, sway, heave, yaw]", 0.0,
              minimum_exclusive=True),
    Parameter("delta_zeta_45", (2, 1), "1", "damping-ratio increments of roll and pitch", 0.0),
)


def move_to_body_origin(source="matrix_at_table_reference", target="matrix", *, symmetric=False):
    """``(<source>, table_reference_point) -> <target>``: ``H(r)^T X_R H(r)`` (Fossen 2011, eq. 3.24, p. 49), with
    ``symmetric=True`` its symmetric part ``(X + X^T)/2`` (a mass matrix)."""
    x_r, r = ca.SX.sym(source, 6, 6), ca.SX.sym("table_reference_point", 3)
    h = MatrixTransforms.H_matrix(r)  # (Fossen 2011, eq. 3.24, p. 49)
    moved = h.T @ x_r @ h  # X_CO = H(r)^T X_R H(r) (Fossen 2011, eq. 3.27, p. 50)
    return function_from(f"move_{source}_to_body_origin", {source: x_r, "table_reference_point": r},
                         {target: MatrixTransforms.symmetrize(moved) if symmetric else moved})


def move_to_table_reference(source="matrix", target="matrix_at_table_reference"):
    """``(<source>, table_reference_point) -> <target>``: ``H(-r)^T X_CO H(-r)``, the inverse move."""
    x_o, r = ca.SX.sym(source, 6, 6), ca.SX.sym("table_reference_point", 3)
    h = MatrixTransforms.H_matrix(-r)  # H(r)^-1 = H(-r) (Fossen 2011, eq. 3.24, p. 49)
    return function_from(f"move_{source}_to_table_reference", {source: x_o, "table_reference_point": r},
                         {target: h.T @ x_o @ h})


def restoring_stiffness_345():
    """``hydrostatic_stiffness_345 -> restoring_stiffness``: the 6 x 6 restoring stiffness ``G`` about the table
    reference point with the given heave-roll-pitch block (row-major 3 x 3), zero elsewhere."""
    block = ca.SX.sym("hydrostatic_stiffness_345", 9)
    g = ca.SX.zeros(6, 6)
    g[2:5, 2:5] = ca.reshape(block, 3, 3).T
    return function_from("restoring_stiffness_345", {"hydrostatic_stiffness_345": block}, {"restoring_stiffness": g})


RESTORING_STIFFNESS_PARAMETERS = (
    Parameter("hydrostatic_stiffness_345", (9, 1), "N/m, N, N m/rad",
              "restoring stiffness of heave, roll and pitch about the table reference point, 3 x 3 row-major"),
)
TABLE_REFERENCE_PARAMETERS = (
    Parameter("database_reference", (3, 1), "m", "table reference point from the CO, BODY FRD"),
)


def table_reference_point():
    """``database_reference -> table_reference_point``: the point the tables are about, from the CO, BODY FRD."""
    r = ca.SX.sym("database_reference", 3)
    return function_from("table_reference_point", {"database_reference": r}, {"table_reference_point": r})


def submerged_viscous_damping():
    """``(mass_matrix_at_table_reference, potential_damping_at_table_reference, restoring_stiffness, T_1236,
    delta_zeta_45) -> (viscous_damping, negative_time_constant_margin)``: the diagonal ``Bv`` (module docstring) and, per
    time-constant degree of freedom, ``M_ii/T_i - B_eq,ii + tol max(1, |M_ii/T_i|, |B_eq,ii|)`` (negative = refused)."""
    m, b = ca.SX.sym("mass_matrix_at_table_reference", 6, 6), ca.SX.sym("potential_damping_at_table_reference", 6, 6)
    g = ca.SX.sym("restoring_stiffness", 6, 6)
    t, dz = ca.SX.sym("T_1236", 4), ca.SX.sym("delta_zeta_45", 2)
    bv = ca.SX.zeros(6, 6)
    margin = []
    for j, i in enumerate(TIME_CONSTANT_DOFS):
        target = m[i, i] / t[j]
        viscous = target - b[i, i]  # Bv_ii = M_ii / T_i - B_eq,ii (Fossen 2025)
        margin.append(viscous + DAMPING_TOLERANCE * ca.fmax(1.0, ca.fmax(ca.fabs(target), ca.fabs(b[i, i]))))
        bv[i, i] = ca.fmax(0.0, viscous)
    for j, i in enumerate(RATIO_DOFS):
        bv[i, i] = 2.0 * dz[j] * ca.sqrt(m[i, i] * g[i, i])  # Bv_ii = 2 delta_zeta sqrt(M_ii G_ii) (Fossen 2025)
    return function_from("submerged_viscous_damping",
                         {"mass_matrix_at_table_reference": m, "potential_damping_at_table_reference": b,
                          "restoring_stiffness": g, "T_1236": t, "delta_zeta_45": dz},
                         {"viscous_damping": bv, "negative_time_constant_margin": ca.vertcat(*margin)})


def scaled_table_damping():
    """``(damping_at_table_reference, relative_velocity, surge_damping_blend, kappa_456, table_reference_point) ->
    damping_matrix``: ``D_11 exp(-k_u |u_r|)`` and ``D_ii (1 + kappa_i |nu_r,i|)`` for roll, pitch and yaw at the table
    reference point, then moved to the body origin, ``H(r)^T D H(r)`` (the angular part of ``nu_r`` is the same at both
    points)."""
    d = ca.SX.sym("damping_at_table_reference", 6, 6)
    nu_r = ca.SX.sym("nu_r", 6)
    k_u, kappa, r = ca.SX.sym("surge_damping_blend"), ca.SX.sym("kappa_456", 3), ca.SX.sym("table_reference_point", 3)
    scale = [ca.exp(-k_u * ca.fabs(nu_r[0])), 1.0, 1.0] + [1.0 + kappa[j] * ca.fabs(nu_r[3 + j]) for j in range(3)]
    scaled = d + ca.diag(ca.diag(d) * (ca.vertcat(*scale) - 1.0))  # only the diagonal entries are scaled
    h = MatrixTransforms.H_matrix(r)  # (Fossen 2011, eq. 3.24, p. 49)
    return function_from("scaled_table_damping",
                         {"damping_at_table_reference": d, "nu_r": nu_r, "surge_damping_blend": k_u,
                          "kappa_456": kappa, "table_reference_point": r},
                         {"D": h.T @ scaled @ h})
