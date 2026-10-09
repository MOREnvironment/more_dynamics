"""The hull hydrodynamic force of a body whose damping, lift and the rigid
body's own Coriolis-centripetal terms are given together as one set of
published nondimensional derivatives (Healey and Lienard 1993, via MSS
``npsauv.m`` 241-274), rather than separate closed-form blocks
(``AuvHullLoads``'s cylinder damping + lift/drag).

``npsauv.m`` does not apply a separate ``-C(nu) nu_r`` term in its equation of
motion (unlike ``remus100.m`` 257-258): the rigid-body Coriolis-centripetal
force, evaluated with the relative velocity, is written directly into
``X_h..N_h`` alongside the nondimensional hydrodynamic derivatives (MSS
``npsauv.m`` 241-274). This block reproduces that sum term by term, so the
vehicle applies **no** rigid-body Coriolis matrix of its own for this type
(``DerivativeAuv.coriolis_matrix`` is zero) — the mass and the full inertia
tensor are couplings of this block, not of a separate ``C_RB`` computation.

**Candidate MSS inconsistency, reproduced, not fixed (rule 15, not
registered here — A-65 owns the register).** The pitch-moment row carries
``- (Iz - Ix) p r`` (``npsauv.m`` 265), the opposite sign of the same
gyroscopic-coupling pattern in the roll row, ``+ (Iy - Iz) q r``
(``npsauv.m`` 259), and the yaw row, ``+ (Ix - Iy) p q`` (``npsauv.m`` 271) —
checked against the standard rigid-body Coriolis form
(``more_dynamics.models.rigid_body.kinetics.rigid_body_coriolis_casadi``,
"co" form, which gives ``+ (Iz - Ix) p r`` for the same body): this block
reproduces MSS's literal sign (verified byte-for-byte against the frozen
``npsauv_hull_derivative.csv``, max |diff| 8.9e-16 over 540 cases), since the
gate compares this vehicle against MSS's own numbers, not against physics.

Inputs by name, from the vehicle: ``relative_velocity``, ``mass``,
``inertia_diagonal``, ``inertia_products``, ``center_of_gravity``,
``length``, ``water_density``. Own parameters: the 52 hull-only
nondimensional derivatives (the remaining 41 of the published 93 are the 14
``*dot`` acceleration derivatives, read by the added-mass block, and the 27
that couple to a fin deflection or the propeller, read by the actuators).

References
----------
[MSS] Fossen, T. I. MSS, MIT, CRAFT/AUV/models/npsauv.m 118-121, 241-274
    @ cc07579.

Author:    Enio Krizman
Date:      2026-10-09
"""
import casadi as ca

from more_transformations.more_casadi_transformations import Parameter

from more_dynamics.models.shared.wiring import function_from

HULL_ONLY_DERIVATIVE_NAMES = {
    "X": ("Xpp", "Xqq", "Xrr", "Xpr", "Xwq", "Xvp", "Xvr", "Xvv", "Xww"),
    "Y": ("Ypq", "Yqr", "Yp", "Yr", "Yvq", "Ywp", "Ywr", "Yv", "Yvw"),
    "Z": ("Zpp", "Zpr", "Zrr", "Zq", "Zvp", "Zvr", "Zw", "Zvv"),
    "K": ("Kpq", "Kqr", "Kp", "Kr", "Kvq", "Kwp", "Kwr", "Kv", "Kvw"),
    "M": ("Mpp", "Mpr", "Mrr", "Muq", "Mvp", "Mvr", "Muw", "Mvv"),
    "N": ("Npq", "Nqr", "Np", "Nr", "Nvq", "Nwp", "Nwr", "Nv", "Nvw"),
}


def coefficient_hull_loads_parameters():
    """The 52 hull-only derivatives, each its own declared parameter (rule
    16: primitives, not a re-derivation)."""
    return tuple(Parameter(f"nondim_{n}", (1, 1), "1", f"nondimensional {n} (npsauv.m 126-156)")
                for names in HULL_ONLY_DERIVATIVE_NAMES.values() for n in names)


def coefficient_hull_loads():
    """``(relative_velocity, mass, inertia_diagonal, inertia_products,
    center_of_gravity, length, water_density, <52 nondim_*>) ->
    hydrodynamic_force`` (6x1, N and N m, BODY, about the CO): the
    nondimensional-derivative hull force and the rigid-body
    Coriolis-centripetal terms of ``npsauv.m`` 241-274, together (module
    docstring)."""
    nu_r = ca.SX.sym("relative_velocity", 6)
    u_r, v_r, w_r, p, q, r = (nu_r[i] for i in range(6))
    mass = ca.SX.sym("mass")
    inertia_diagonal = ca.SX.sym("inertia_diagonal", 3)
    inertia_products = ca.SX.sym("inertia_products", 3)  # [Ixy, Iyz, Ixz]
    r_g = ca.SX.sym("center_of_gravity", 3)
    length = ca.SX.sym("length")
    water_density = ca.SX.sym("water_density")

    names = [f"nondim_{n}" for names in HULL_ONLY_DERIVATIVE_NAMES.values() for n in names]
    s = {n: ca.SX.sym(n) for n in names}
    d = {n: s[f"nondim_{n}"] for n in (n for names in HULL_ONLY_DERIVATIVE_NAMES.values() for n in names)}

    r2 = 0.5 * water_density * length ** 2  # (npsauv.m 118)
    r3 = 0.5 * water_density * length ** 3  # (npsauv.m 119)
    r4 = 0.5 * water_density * length ** 4  # (npsauv.m 120)
    r5 = 0.5 * water_density * length ** 5  # (npsauv.m 121)
    xG, yG, zG = r_g[0], r_g[1], r_g[2]
    Ix, Iy, Iz = inertia_diagonal[0], inertia_diagonal[1], inertia_diagonal[2]
    Ixy, Iyz, Ixz = inertia_products[0], inertia_products[1], inertia_products[2]

    # (npsauv.m 241-244, the m*(...) rigid-body terms of X_h kept with it, module docstring)
    X_h = r2 * (d["Xvv"] * v_r ** 2 + d["Xww"] * w_r ** 2) + \
        r3 * (d["Xvr"] * v_r * r + d["Xwq"] * w_r * q + d["Xvp"] * v_r * p) + \
        r4 * (d["Xqq"] * q ** 2 + d["Xrr"] * r ** 2 + d["Xpr"] * p * r + d["Xpp"] * p ** 2) + \
        mass * (v_r * r - w_r * q + xG * (q ** 2 + r ** 2) - yG * p * q - zG * p * r)

    # (npsauv.m 246-249)
    Y_h = r2 * (d["Yv"] * u_r * v_r + d["Yvw"] * v_r * w_r) + \
        r3 * (d["Yp"] * u_r * p + d["Yr"] * u_r * r + d["Yvq"] * v_r * q + d["Ywp"] * w_r * p +
              d["Ywr"] * w_r * r) + \
        r4 * (d["Ypq"] * p * q + d["Yqr"] * q * r) - \
        mass * (u_r * r - w_r * p + xG * p * q - yG * (p ** 2 + r ** 2) + zG * q * r)

    # (npsauv.m 251-254)
    Z_h = r2 * (d["Zw"] * w_r * u_r + d["Zvv"] * v_r ** 2) + \
        r3 * (d["Zq"] * u_r * q + d["Zvp"] * v_r * p + d["Zvr"] * v_r * r) + \
        r4 * (d["Zpp"] * p ** 2 + d["Zpr"] * p * r + d["Zrr"] * r ** 2) + \
        mass * (v_r * p - u_r * q + xG * p * r + yG * q * r - zG * (p ** 2 + q ** 2))

    # (npsauv.m 256-260; (Iy - Iz) q r positive, the roll-plane counterpart of the module docstring's flag)
    K_h = r3 * (d["Kv"] * u_r * v_r + d["Kvw"] * v_r * w_r) + \
        r4 * (d["Kp"] * u_r * p + d["Kr"] * u_r * r + d["Kvq"] * v_r * q + d["Kwp"] * w_r * p +
              d["Kwr"] * w_r * r) + \
        r5 * (d["Kpq"] * p * q + d["Kqr"] * q * r) + \
        (Iy - Iz) * q * r - Ixy * p * r - (r ** 2 - q ** 2) * Iyz + Ixz * p * q - \
        mass * (yG * (v_r * p - u_r * q) - zG * (u_r * r - w_r * p))

    # (npsauv.m 262-266; "- (Iz - Ix) p r" literal, the flagged candidate MSS sign, module docstring)
    M_h = r3 * (d["Muw"] * u_r * w_r + d["Mvv"] * v_r ** 2) + \
        r4 * (d["Muq"] * u_r * q + d["Mvp"] * v_r * p + d["Mvr"] * v_r * r) + \
        r5 * (d["Mpp"] * p ** 2 + d["Mpr"] * p * r + d["Mrr"] * r ** 2) - \
        (Iz - Ix) * p * r + Ixy * q * r - Iyz * p * q - (p ** 2 - r ** 2) * Ixz + \
        mass * (xG * (v_r * p - u_r * q) - zG * (w_r * q - v_r * r))

    # (npsauv.m 268-272)
    N_h = r3 * (d["Nv"] * u_r * v_r + d["Nvw"] * v_r * w_r) + \
        r4 * (d["Np"] * u_r * p + d["Nr"] * u_r * r + d["Nvq"] * v_r * q + d["Nwp"] * w_r * p +
              d["Nwr"] * w_r * r) + \
        r5 * (d["Npq"] * p * q + d["Nqr"] * q * r) + \
        (Ix - Iy) * p * q + (p ** 2 - q ** 2) * Ixy + Iyz * p * r - Ixz * q * r - \
        mass * (xG * (u_r * r - w_r * p) - yG * (w_r * q - v_r * r))

    tau = ca.vertcat(X_h, Y_h, Z_h, K_h, M_h, N_h)  # (npsauv.m 274)
    return function_from(
        "coefficient_hull_loads",
        {"relative_velocity": nu_r, "mass": mass, "inertia_diagonal": inertia_diagonal,
         "inertia_products": inertia_products, "center_of_gravity": r_g, "length": length,
         "water_density": water_density, **s},
        {"hydrodynamic_force": tau})
