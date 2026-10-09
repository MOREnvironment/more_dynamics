"""Added mass: the derivative sets and ``M_A`` (CasADi). ``C_A(nu_r)`` is
built in ``kinetics.py``.

Every function takes CasADi symbols or numbers and returns a CasADi
expression; nothing here checks a value (the admissible ranges are declared
in ``mass_properties.py`` and checked where numbers enter).

Equations (keys in References):

* ``M_A = -diag(X_du, Y_dv, Z_dw, K_dp, M_dq, N_dr)`` (Fossen 2011, eqs. 6.38,
  p. 118, and 6.53, p. 121; MSS ``otter.m`` 159).
* Surge: ``A11 = 2.7 rho nabla^(5/3) / L^2``, ``nabla = m / rho`` (MSS
  ``addedMassSurge.m`` 32-33, which attributes it to Söding (1982), not read
  here).
* Scaled derivatives: ``[X_du, Y_dv, Z_dw, K_dp, M_dq, N_dr] = c * [A11, m, m,
  I11, I22, I33]``, the rotational terms scaling the inertia about the CG
  (the pattern of MSS ``otter.m`` 152-157).
* Prolate spheroid, ``a > b``, eccentricity ``e = sqrt(1 - (b/a)^2)``:
  ``alpha_0 = 2 (1 - e^2) / e^3 (1/2 ln((1+e)/(1-e)) - e)``,
  ``beta_0 = 1/e^2 - (1 - e^2) / (2 e^3) ln((1+e)/(1-e))``,
  ``k1 = alpha_0 / (2 - alpha_0)``, ``k2 = beta_0 / (2 - beta_0)`` (Lamb 1932,
  Art. 114, eqs. 12-15, pp. 153-154; MSS ``imlay61.m`` 50-55);
  ``k' = e^4 (beta_0 - alpha_0) / ((2 - e^2)(2 e^2 - (2 - e^2)(beta_0 -
  alpha_0)))`` (Lamb 1932, Art. 115, eq. 8, p. 155; MSS ``imlay61.m`` 56);
  ``M_A = diag(k1 m_f, k2 m_f, k2 m_f, r44 I_x, k' I_y, k' I_y)`` of the
  displaced fluid (MSS ``imlay61.m`` 47, 59).
* Stabilised ``C_A`` entries: MSS ``remus100.m`` 207-210.

References
----------
[Fossen 2011] Fossen, T. I. (2011). *Handbook of Marine Craft Hydrodynamics
    and Motion Control*, 1st ed. John Wiley & Sons, Chichester. Ch. 6,
    eqs. 6.38, 6.53, pp. 118-121.
[Lamb 1932] Lamb, H. (1932). *Hydrodynamics*, 6th ed. Cambridge University
    Press (Dover reprint, 1945). Ch. V, Arts. 114-115, pp. 152-155.
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2.
    https://github.com/cybergalactic/MSS, MIT licence, revision ``72656d1``:
    ``LIBRARY/modeling/addedMassSurge.m`` 32-33, ``imlay61.m`` 47-59;
    ``CRAFT/USV/models/otter.m`` 152-159; ``CRAFT/AUV/models/remus100.m``
    207-210.

Author:    Enio Krizman
Date:      2026-10-05
"""

import casadi as ca

# C_A entries zeroed by the stabilised form, 0-based (row, column); each
# pair is zeroed in both triangles: pitch-heave, pitch-surge, yaw-surge,
# yaw-sway (MSS remus100.m 207-210).
STABILIZED_ADDED_MASS_CORIOLIS_PAIRS = ((4, 2), (4, 0), (5, 0), (5, 1))


def surge_added_mass(mass, length, water_density):
    """``A11 = 2.7 rho nabla^(5/3) / L^2``, ``nabla = m / rho`` (MSS
    ``addedMassSurge.m`` 32-33, after Söding 1982, not read)."""
    displaced_volume = mass / water_density  # (addedMassSurge.m 32)
    return 2.7 * water_density * displaced_volume ** (5.0 / 3.0) / length**2  # (addedMassSurge.m 33)


def scaled_added_mass_derivatives(mass, length, water_density, inertia, coefficients):
    """``[X_du, Y_dv, Z_dw, K_dp, M_dq, N_dr] = c * [A11, m, m, I11, I22, I33]``
    as a 6x1 column.

    ``mass`` is the hull mass, ``inertia`` (3x3) the one the rigid body uses (about the
    CG, as MSS ``otter.m`` 152-157). The coefficients carry the sign:
    negative coefficients give a positive ``M_A``.
    """
    base = ca.vertcat(
        surge_added_mass(mass, length, water_density),  # A11 (otter.m 152; addedMassSurge.m 32-33)
        mass,  # (otter.m 153)
        mass,  # (otter.m 154)
        inertia[0, 0],  # (otter.m 155)
        inertia[1, 1],  # (otter.m 156)
        inertia[2, 2],  # (otter.m 157)
    )
    return coefficients * base  # c * [A11, m, m, I11, I22, I33] (pattern of otter.m 152-157)


def lamb_k_factors(semi_major_axis, semi_minor_axis):
    """Lamb's ``[k1, k2, k']`` of a prolate spheroid ``a > b > 0``, a 3x1
    column.

    ``e = sqrt(1 - (b/a)^2)``, ``alpha_0``, ``beta_0``, ``k1 = alpha_0 /
    (2 - alpha_0)``, ``k2 = beta_0 / (2 - beta_0)`` (Lamb 1932, Art. 114,
    eqs. 12-15, pp. 153-154), ``k' = e^4 (beta_0 - alpha_0) / ((2 - e^2)(2 e^2
    - (2 - e^2)(beta_0 - alpha_0)))`` (Lamb 1932, Art. 115, eq. 8, p. 155);
    MSS ``imlay61.m`` 50-56.
    The expressions are real only for ``a > b > 0`` (eccentricity ``0 < e <
    1``); that condition is checked on numbers by
    ``kinetics.check_rigid_body_values``, stricter than MSS ``imlay61.m``
    40-42, which accepts ``b = 0``.
    """
    a, b = semi_major_axis, semi_minor_axis
    e = ca.sqrt(1.0 - (b / a) ** 2)  # (imlay61.m 50)
    log_term = ca.log((1.0 + e) / (1.0 - e))  # ln((1+e)/(1-e)) of (Lamb 1932, eq. 14, p. 154; imlay61.m 51-52)
    alpha_0 = (2.0 * (1.0 - e**2) / e**3) * (0.5 * log_term - e)  # (Lamb 1932, eq. 14, p. 154; imlay61.m 51)
    beta_0 = 1.0 / e**2 - (1.0 - e**2) / (2.0 * e**3) * log_term  # (Lamb 1932, eq. 14, p. 154; imlay61.m 52)
    k1 = alpha_0 / (2.0 - alpha_0)  # (Lamb 1932, eqs. 12, 15, pp. 153-154; imlay61.m 54)
    k2 = beta_0 / (2.0 - beta_0)  # (Lamb 1932, eq. 15, p. 154; imlay61.m 55)
    # (Lamb 1932, Art. 115, eq. 8, p. 155; imlay61.m 56)
    k_prime = e**4 * (beta_0 - alpha_0) / (
        (2.0 - e**2) * (2.0 * e**2 - (2.0 - e**2) * (beta_0 - alpha_0))
    )
    return ca.vertcat(k1, k2, k_prime)


def spheroid_added_mass_derivatives(k_factors, fluid_mass, fluid_inertia, roll_added_inertia_ratio):
    """``-[k1 m_f, k2 m_f, k2 m_f, r44 I_x, k' I_y, k' I_y]`` as a 6x1 column
    (MSS ``imlay61.m`` 47, 59).

    ``m_f`` and ``I_f`` (3x3) are the mass and inertia of the displaced
    fluid. With ``r44 = 0`` the roll term vanishes: Imlay (1961) gives a zero
    added moment in roll for a spheroid (MSS ``imlay61.m`` 36-38).
    """
    k1, k2, k_prime = k_factors[0], k_factors[1], k_factors[2]
    # -diag of MA = diag([m k1, m k2, m k2, r44 Ix, k' Iy, k' Iy]) (imlay61.m 47, 59)
    return -ca.vertcat(
        k1 * fluid_mass,  # (imlay61.m 59)
        k2 * fluid_mass,  # (imlay61.m 59)
        k2 * fluid_mass,  # (imlay61.m 59)
        roll_added_inertia_ratio * fluid_inertia[0, 0],  # MA_44 = r44 Ix (imlay61.m 47)
        k_prime * fluid_inertia[1, 1],  # (imlay61.m 59)
        k_prime * fluid_inertia[1, 1],  # (imlay61.m 59)
    )


def added_mass_matrix(derivatives):
    """``M_A = -diag(derivatives)`` (Fossen 2011, eq. 6.53, p. 121; MSS
    ``otter.m`` 159)."""
    return ca.diag(-derivatives)  # (Fossen 2011, eq. 6.53, p. 121; otter.m 159)
