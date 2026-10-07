"""Added mass: the derivative sets and ``M_A`` (numpy). ``C_A(nu_r)`` is
built in ``kinetics.py``.

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

Ported from the numpy source [MGM] ``dynamics/plant/matrices/added_mass.py``:
``M_A_6dof`` (24-41), ``M_A_lamb_6dof`` (45-78), ``get_added_mass_derivates``
(174-267), ``get_added_mass_derivates_catamaran`` (269-289) and
``get_added_mass_derivates_hull`` (291-310, the same formula: one function
here), ``added_mass_surge`` / ``added_mass_surge_static`` (313-360), and the
entries ``stabilize_C_A`` (156-171) zeroes. Deviation from that source: its
``normalize_zero`` (17-18), which turned values below 1e-8 into 0, is not
ported (dropped by the owner's decision of 2026-10-06).

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
[MGM] Krizman, E. *more_generic_models*.
    https://github.com/MOREnvironment/more_generic_models (no licence file),
    revision ``524e336``: ``more_generic_models/dynamics/plant/matrices/
    added_mass.py``, lines listed above.
"""

import numpy as np

# C_A entries zeroed by the stabilised form, 0-based (row, column); each
# pair is zeroed in both triangles: pitch-heave, pitch-surge, yaw-surge,
# yaw-sway (MSS remus100.m 207-210; source stabilize_C_A 160-169).
STABILIZED_ADDED_MASS_CORIOLIS_PAIRS = ((4, 2), (4, 0), (5, 0), (5, 1))


def surge_added_mass(mass: float, length: float, water_density: float) -> float:
    """``A11 = 2.7 rho nabla^(5/3) / L^2``, ``nabla = m / rho`` (MSS
    ``addedMassSurge.m`` 32-33, after Söding 1982, not read; source
    ``added_mass_surge``)."""
    displaced_volume = mass / water_density  # (addedMassSurge.m 32)
    return 2.7 * water_density * displaced_volume ** (5.0 / 3.0) / length**2  # (addedMassSurge.m 33)


def scaled_added_mass_derivatives(
    mass: float,
    length: float,
    water_density: float,
    inertia: np.ndarray,
    coefficients: np.ndarray,
) -> np.ndarray:
    """``[X_du, Y_dv, Z_dw, K_dp, M_dq, N_dr] = c * [A11, m, m, I11, I22, I33]``.

    Source ``get_added_mass_derivates_catamaran`` / ``_hull``; ``mass`` is
    the hull mass, ``inertia`` the one the rigid body uses (about the CG, as
    MSS ``otter.m`` 152-157). The coefficients carry the sign: negative
    coefficients give a positive ``M_A``.
    """
    # c * [A11, m, m, I11, I22, I33] (pattern of otter.m 152-157)
    return np.asarray(coefficients, dtype=float) * np.array(
        [
            surge_added_mass(mass, length, water_density),
            mass,
            mass,
            inertia[0, 0],
            inertia[1, 1],
            inertia[2, 2],
        ]
    )


def lamb_k_factors(semi_major_axis: float, semi_minor_axis: float) -> np.ndarray:
    """Lamb's ``[k1, k2, k']`` of a prolate spheroid, ``a > b > 0``.

    ``e = sqrt(1 - (b/a)^2)``, ``alpha_0``, ``beta_0``, ``k1 = alpha_0 /
    (2 - alpha_0)``, ``k2 = beta_0 / (2 - beta_0)`` (Lamb 1932, Art. 114,
    eqs. 12-15, pp. 153-154), ``k' = e^4 (beta_0 - alpha_0) / ((2 - e^2)(2 e^2
    - (2 - e^2)(beta_0 - alpha_0)))`` (Lamb 1932, Art. 115, eq. 8, p. 155);
    MSS ``imlay61.m`` 50-56; source ``get_added_mass_derivates`` 176-261.
    The input check is stricter than MSS ``imlay61.m`` 40-42, which accepts
    ``b = 0`` (e = 1, infinite logarithm): here ``b > 0`` is required.
    """
    a, b = float(semi_major_axis), float(semi_minor_axis)
    if not (np.isfinite(a) and np.isfinite(b)) or b <= 0.0 or a <= b:
        raise ValueError(
            "semi_major_axis must be larger than semi_minor_axis, both positive "
            "(prolate spheroid, eccentricity 0 < e < 1)"
        )
    e = np.sqrt(1.0 - (b / a) ** 2)  # (imlay61.m 50)
    log_term = np.log((1.0 + e) / (1.0 - e))
    alpha_0 = (2.0 * (1.0 - e**2) / e**3) * (0.5 * log_term - e)  # (Lamb 1932, eq. 14, p. 154; imlay61.m 51)
    beta_0 = 1.0 / e**2 - (1.0 - e**2) / (2.0 * e**3) * log_term  # (Lamb 1932, eq. 14, p. 154; imlay61.m 52)
    k1 = alpha_0 / (2.0 - alpha_0)  # (Lamb 1932, eqs. 12, 15, pp. 153-154; imlay61.m 54)
    k2 = beta_0 / (2.0 - beta_0)  # (Lamb 1932, eq. 15, p. 154; imlay61.m 55)
    # (Lamb 1932, Art. 115, eq. 8, p. 155; imlay61.m 56)
    k_prime = e**4 * (beta_0 - alpha_0) / (
        (2.0 - e**2) * (2.0 * e**2 - (2.0 - e**2) * (beta_0 - alpha_0))
    )
    return np.array([k1, k2, k_prime])


def spheroid_added_mass_derivatives(
    k_factors: np.ndarray,
    fluid_mass: float,
    fluid_inertia: np.ndarray,
    roll_added_inertia_ratio: float,
) -> np.ndarray:
    """``-[k1 m_f, k2 m_f, k2 m_f, r44 I_x, k' I_y, k' I_y]`` (MSS ``imlay61.m``
    47, 59; source ``M_A_lamb_6dof`` 68-75).

    ``m_f`` and ``I_f`` are the mass and inertia of the displaced fluid.
    With ``r44 = 0`` this equals the source's other derivative set
    (``get_added_mass_derivates`` 207-241, which cites Fossen 2021, 2nd ed.,
    eqs. 8.74-8.77, not on disk): ``N_dr = -k' I_y`` to rounding.
    """
    k1, k2, k_prime = k_factors
    # -diag of MA = diag([m k1, m k2, m k2, r44 Ix, k' Iy, k' Iy]) (imlay61.m 47, 59)
    return -np.array(
        [
            k1 * fluid_mass,
            k2 * fluid_mass,
            k2 * fluid_mass,
            roll_added_inertia_ratio * fluid_inertia[0, 0],
            k_prime * fluid_inertia[1, 1],
            k_prime * fluid_inertia[1, 1],
        ]
    )


def added_mass_matrix(derivatives: np.ndarray) -> np.ndarray:
    """``M_A = -diag(derivatives)`` (Fossen 2011, eq. 6.53, p. 121; MSS
    ``otter.m`` 159; source ``M_A_6dof``)."""
    return np.diag(-np.asarray(derivatives, dtype=float))  # (Fossen 2011, eq. 6.53, p. 121)
