"""U8a: independent references for Grethe's ``Monohull`` (numpy, no CasADi, no
``more_dynamics`` import — a second implementation of the cited equations, G5
against the published anchors, G2 against today's library formulas where they
already exist). There is no MSS model of Grethe (she is not in MSS), so these
are not a second-implementation check against a frozen MATLAB run: they are
written directly from the printed equations and checked for internal
consistency (symmetry, skew-symmetry, passivity) the way ``a64_math_check.py``
checked the library itself (`agents-more/20_sources/monohull_structure_for_grethe.md`
Section 7).

Run:  python3 -I generate_grethe_references.py
      writes every CSV/JSON beside this file; a re-run is checked with ``cmp``
      against the committed files (no external dependency: deterministic
      arithmetic, no MSS_DIR, no MATLAB_BIN).

Sources (numbers and forms, rule 4 — none typed from memory without a citation)
---------------------------------------------------------------------------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics and
    Motion Control. Wiley. eq. 2.10 p. 20 (skew); eq. 3.24 p. 49 (H); eq. 3.26,
    3.34 p. 50 (M_RB); eq. 3.46 p. 53 / eq. 6.43 p. 120 (m2c); eq. 6.53 p. 121
    (M_A); eq. 6.52 p. 121 (Munk moment); eq. 4.24 p. 65, 4.32-4.35 p. 65-66,
    7.250 p. 181 (metacentric restoring); eqs. 6.82-6.85 p. 125 (ITTC + C_T =
    (1+k) C_F + C_R); eqs. 6.91-6.92 p. 127 (cross-flow strips).
[MSS] Fossen, T. I. MSS, MIT, @ cc07579: otter.m 94-98, 104-107, 121-128,
    152-159, 172-193; addedMassSurge.m 32-33; Gmtrx.m 25-37; Hoerner.m 25-51;
    crossFlowDrag.m 36-69; XuuITTC.m 32-39; exShipHydrostatics.m 32-56 (the
    same lines the library's own modules cite; reproduced here independently
    in numpy, not imported).
[RA14] Radojcic, D., Zgradic, A., Kalajdzic, M., Simic, A. (2014). Resistance
    prediction for hard chine hulls in the pre-planing regime. Polish Maritime
    Research 21(2):9-26. Appendix 1 p. 24 (the simple model: R/Delta, S/V^(2/3)
    cubic in Fn_V, coefficients quartic/cubic in L/V^(1/3)); basis p. 11
    (Delta = 100000 lb, rho = 1026, nu = 1.1907e-6 m^2/s, ITTC-1957, C_A = 0).
    Coefficients copied from `agents-more/60_working/A64_scripts/a64_hull_regime_resistance.py`
    (A-64, transcribed and verified there from the cited page); not retyped
    from memory here (rule 4).
[E-101] `agents-more/80_owner/inbox/E-101_a64_grethe_monohull_choices_remus_density_at_repin.md`:
    Q1 a (ITTC + RA14 residual default), Q2 a (Y_r = N_v = 0, unidentified),
    Q4 a (CO at the waterline midpoint, x_P = -2.24 m), Q5 a (added mass on
    the displaced mass).
[A-43] `agents-more/60_working/2026-10-07_A43_vehicle_parameter_dossiers.md`
    Part 2 (M1-M50): every Grethe value's source and kind.
[A-64] `agents-more/20_sources/monohull_structure_for_grethe.md` Section 1.1
    (drawing measurements), Section 2 (fidelity ladder), Section 3.3 (the
    option names and defaults this file freezes).

Author:    Enio Krizman
Date:      2026-10-09
"""
import csv
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
G = 9.81  # estimate: grethe_mariner5.yaml hydrostatics.gravity ("standard"), A-43 Part 2 header
RHO = 1025.0  # published: standard sea water (A-43 M5)
NU_WATER = 1.0e-6  # Fossen 2011, p. 125, below eq. 6.85 (20 degC); cylinderDrag.m 78-80


# ---------------------------------------------------------------------------
# Rigid body + added mass (Fossen 2011, Ch. 3 and Sec. 6.3; MSS otter.m)
# ---------------------------------------------------------------------------

def skew(v):
    """S(v) a = v x a (Fossen 2011, eq. 2.10, p. 20; Smtrx.m 11-13)."""
    x, y, z = v
    return np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])


def h_matrix(r):
    """H(r) = [[I3, S(r)^T], [0, I3]] (Fossen 2011, eq. 3.24, p. 49; Hmtrx.m 16-18)."""
    s = skew(r)
    top = np.hstack([np.eye(3), s.T])
    bottom = np.hstack([np.zeros((3, 3)), np.eye(3)])
    return np.vstack([top, bottom])


def rigid_body_hull_with_payload(hull_mass, payload_mass, hull_cg, payload_position, radii_of_gyration,
                                 beam, length):
    """``m = m_h + m_p``, ``r_g = (m_h r_h + m_p r_p) / m``, ``I = I_h - m_h
    S(r_h - r_g)^2 - m_p S(r_p - r_g)^2`` (Fossen 2011, Theorem 3.1, eq. 3.34,
    p. 50; MSS otter.m 94-98, 122-128)."""
    hull_cg, payload_position = np.asarray(hull_cg, dtype=float), np.asarray(payload_position, dtype=float)
    mass = hull_mass + payload_mass  # (otter.m 124)
    radii = np.asarray(radii_of_gyration, dtype=float) * np.array([beam, length, length])  # (otter.m 95-97)
    inertia_hull_cg = np.diag(hull_mass * radii**2)  # (otter.m 122)
    r_g = (hull_mass * hull_cg + payload_mass * payload_position) / mass  # (otter.m 124)
    s_hull, s_payload = skew(hull_cg - r_g), skew(payload_position - r_g)  # (otter.m 125-126)
    inertia = inertia_hull_cg - hull_mass * s_hull @ s_hull - payload_mass * s_payload @ s_payload  # (otter.m 127)
    return mass, r_g, inertia


def m_rb(mass, inertia, r_g):
    """M_RB = H(r_g)^T diag(m I3, I) H(r_g) (Fossen 2011, eq. 3.26, p. 50; rbody.m 43-44)."""
    mass_at_cg = np.block([[mass * np.eye(3), np.zeros((3, 3))], [np.zeros((3, 3)), inertia]])
    h = h_matrix(r_g)
    return h.T @ mass_at_cg @ h


def c_rb_co(mass, inertia, r_g, nu):
    """C_RB(nu), ``"co"`` form (Fossen 2011, eq. 3.27, p. 50; rbody.m 39-45): depends on
    the angular velocity ``w = nu[3:6]`` only."""
    w = nu[3:6]
    s_w = skew(w)
    coriolis_at_cg = np.block([[mass * s_w, np.zeros((3, 3))], [np.zeros((3, 3)), -skew(inertia @ w)]])
    h = h_matrix(r_g)
    return h.T @ coriolis_at_cg @ h


def m2c(mass_matrix, nu):
    """``m2c(M, nu)`` (Fossen 2011, Theorem 3.2, eq. 3.46, p. 53, eq. 6.43,
    p. 120; MSS m2c.m 33-48), ``M`` symmetrised first."""
    matrix = 0.5 * (mass_matrix + mass_matrix.T)
    v1, v2 = nu[0:3], nu[3:6]
    linear = matrix[0:3, 0:3] @ v1 + matrix[0:3, 3:6] @ v2
    angular = matrix[0:3, 3:6].T @ v1 + matrix[3:6, 3:6] @ v2
    s_linear, s_angular = skew(linear), skew(angular)
    top = np.hstack([np.zeros((3, 3)), -s_linear])
    bottom = np.hstack([-s_linear, -s_angular])
    return np.vstack([top, bottom])


# pitch-heave, pitch-surge, yaw-surge, yaw-sway, both triangles (remus100.m 207-210)
STABILIZE_PAIRS = ((4, 2), (4, 0), (5, 0), (5, 1))


def c_a(added_mass_matrix, nu_r, stabilize):
    """C_A(nu_r) = m2c(M_A, nu_r) (Fossen 2011, eq. 6.43, p. 120); ``stabilize``
    zeroes the Munk-moment couplings (MSS remus100.m 207-210, eq. 6.52, p. 121)."""
    c = m2c(added_mass_matrix, nu_r)
    if stabilize:
        for i, j in STABILIZE_PAIRS:
            c[i, j] = 0.0
            c[j, i] = 0.0
    return c


def scaled_added_mass_matrix(basis_mass, length, water_density, inertia, coefficients):
    """``M_A = -diag(c * [A11, m, m, I11, I22, I33])``, ``A11 = 2.7 rho (m/rho)^(5/3)
    / L^2`` (MSS addedMassSurge.m 32-33, otter.m 152-159). ``basis_mass`` is the
    hull mass (today's library, deviation flagged A-64 Section 3.1 gap 6) or the
    displaced mass hull+payload (owner E-101 Q5 a)."""
    a11 = 2.7 * water_density * (basis_mass / water_density) ** (5.0 / 3.0) / length**2
    base = np.array([a11, basis_mass, basis_mass, inertia[0, 0], inertia[1, 1], inertia[2, 2]])
    derivatives = np.asarray(coefficients, dtype=float) * base
    return np.diag(-derivatives)


# ---------------------------------------------------------------------------
# Surface hydrostatics, one hull (Fossen 2011, Ch. 4; MSS Gmtrx.m, otter.m,
# exShipHydrostatics.m, through the library's own surface_restoring.py form,
# reproduced independently here)
# ---------------------------------------------------------------------------

def surface_restoring(*, length, beam, block_coefficient, waterplane_coefficient, mass, water_density, gravity,
                      center_of_gravity, longitudinal_inertia_factor, longitudinal_center_of_flotation,
                      reference_point, hull_count=1, hull_lateral_offset=0.0):
    """One or two identical hulls (``hull_count``): ``nabla = m / rho``,
    ``T = nabla / (n Cb L B)`` (otter.m 121-122); ``A_hull = Cw L B``
    (exShipHydrostatics.m 32); Munro-Smith ``I_T`` (exShipHydrostatics.m 38,
    47; otter.m 175-177); ``KB`` (Morrish, exShipHydrostatics.m 35; otter.m
    178); ``GM = KB + BM - KG`` (Fossen 2011, eq. 4.32, p. 65); ``G`` at the
    centre of flotation, the CO and ``reference_point`` (Fossen 2011,
    eq. 4.24, p. 65, eq. 7.250, p. 181; Gmtrx.m 29-37)."""
    n = float(hull_count)
    nabla = mass / water_density
    draft = nabla / (n * block_coefficient * beam * length)
    hull_area = waterplane_coefficient * length * beam
    waterplane_area = n * hull_area
    cw = waterplane_coefficient  # Cw, not Cb (exShipHydrostatics.m 38)
    munro_smith = 6.0 * cw**3 / ((1.0 + cw) * (1.0 + 2.0 * cw))
    transverse_inertia = n * length * beam**3 / 12.0 * munro_smith + n * hull_area * hull_lateral_offset**2
    longitudinal_inertia = n * longitudinal_inertia_factor * beam * length**3 / 12.0
    center_of_buoyancy_above_keel = (2.5 * draft - nabla / n / hull_area) / 3.0  # area_for_kb = hull_area (waterplane default)
    bm_t, bm_l = transverse_inertia / nabla, longitudinal_inertia / nabla
    kg = draft - center_of_gravity[2]
    km_t, km_l = center_of_buoyancy_above_keel + bm_t, center_of_buoyancy_above_keel + bm_l
    gm_t, gm_l = km_t - kg, km_l - kg
    rho_g = water_density * gravity
    g_cf = np.diag([0.0, 0.0, rho_g * waterplane_area, rho_g * nabla * gm_t, rho_g * nabla * gm_l, 0.0])
    h_f = h_matrix(np.array([longitudinal_center_of_flotation, 0.0, 0.0]))
    g_co = h_f.T @ g_cf @ h_f
    h_p = h_matrix(np.asarray(reference_point, dtype=float))
    g = h_p.T @ g_co @ h_p
    wetted_surface_mumford = n * 1.025 * length * (block_coefficient * beam + 1.7 * draft)  # (XuuITTC.m 38, Mumford)
    return {"draft": draft, "displaced_volume": nabla, "waterplane_area": waterplane_area,
            "transverse_metacentric_height": gm_t, "longitudinal_metacentric_height": gm_l,
            "center_of_buoyancy_above_keel": center_of_buoyancy_above_keel,
            "G": g, "G_CO": g_co, "G_CF": g_cf, "wetted_surface_mumford": wetted_surface_mumford}


# ---------------------------------------------------------------------------
# Hoerner cross-flow strips (Fossen 2011, eqs. 6.91-6.92, p. 127; MSS
# crossFlowDrag.m 54-69, Hoerner.m 25-51)
# ---------------------------------------------------------------------------

NUMBER_OF_STRIPS = 20  # (crossFlowDrag.m 37)

HOERNER_DRAG_DATA = (  # [B/(2T), Cd_2D] (Fossen 2011, Fig. 6.5, p. 128; Hoerner.m 25-45)
    (0.0108623, 1.96608), (0.176606, 1.96573), (0.353025, 1.89756), (0.451863, 1.78718),
    (0.472838, 1.58374), (0.492877, 1.27862), (0.493252, 1.21082), (0.558473, 1.08356),
    (0.646401, 0.998631), (0.833589, 0.87959), (0.988002, 0.828415), (1.30807, 0.759941),
    (1.63918, 0.691442), (1.85998, 0.657076), (2.31288, 0.630693), (2.59998, 0.596186),
    (3.00877, 0.586846), (3.45075, 0.585909), (3.7379, 0.559877), (4.00309, 0.559315),
)


def _clamped_interp(x, table):
    """Piecewise-linear f(x), end values held outside the table (cylinderDrag.m
    83-108; Hoerner.m 47-51)."""
    xs, fs = [row[0] for row in table], [row[1] for row in table]
    if x < xs[0]:
        return fs[0]
    if x >= xs[-1]:
        return fs[-1]
    for j in range(len(xs) - 1):
        if xs[j] <= x < xs[j + 1]:
            slope = (fs[j + 1] - fs[j]) / (xs[j + 1] - xs[j])
            return fs[j] + slope * (x - xs[j])
    return fs[-1]


def cross_flow_hoerner(nu_r, length, beam, draft, water_density):
    """Strip integral, 20 midpoint strips, Hoerner's rectangular section
    (Fossen 2011, eqs. 6.91-6.92, p. 127; crossFlowDrag.m 54-69, Hoerner.m 47-51)."""
    cd = _clamped_interp(beam / (2.0 * draft), HOERNER_DRAG_DATA)
    dx = length / NUMBER_OF_STRIPS
    positions = [-length / 2 + (i - 0.5) * dx for i in range(1, NUMBER_OF_STRIPS + 1)]
    factor = -0.5 * water_density * draft * cd * dx
    v_r, w_r, q, r = nu_r[1], nu_r[2], nu_r[4], nu_r[5]
    y_sum = z_sum = m_sum = n_sum = 0.0
    for x in positions:
        horizontal, vertical = v_r + x * r, w_r + x * q
        y_sum += abs(horizontal) * horizontal
        z_sum += abs(vertical) * vertical
        m_sum += x * abs(vertical) * vertical
        n_sum += x * abs(horizontal) * horizontal
    return np.array([0.0, factor * y_sum, factor * z_sum, 0.0, factor * m_sum, factor * n_sum])


# ---------------------------------------------------------------------------
# ITTC friction (Fossen 2011, eqs. 6.82-6.85, p. 125) and the RA14 residual
# prior (Appendix 1, p. 24)
# ---------------------------------------------------------------------------

ITTC_FRICTION_FACTOR = 0.075
ITTC_LOG10_OFFSET = 2.0


def ittc_friction_coefficient(u_r, length, kinematic_viscosity, reynolds_floor):
    """C_F = 0.075 / (log10 Rn - 2)^2, Rn = max(L |u_r| / nu, Rn_min)
    (Fossen 2011, eq. 6.83, p. 125; XuuITTC.m 33-36)."""
    reynolds = max(length * abs(u_r) / kinematic_viscosity, reynolds_floor)
    return ITTC_FRICTION_FACTOR / (math.log10(reynolds) - ITTC_LOG10_OFFSET) ** 2


def ittc_surge_force(u_r, *, mass, length, water_density, wetted_surface, time_constant, form_factor,
                     crossover_speed, kinematic_viscosity, surge_added_mass_factor, reynolds_floor,
                     residual_coefficient=0.0):
    """``X = sigma Xu u_r + (1 - sigma) (Xuu_friction + Xuu_residual) |u_r| u_r``,
    ``Xuu_friction = -1/2 rho S (1+k) C_F`` (Fossen 2011, eqs. 6.82-6.85, p. 125;
    forceSurgeDamping.m 79-82 @ ac77394, XuuITTC.m 38); ``residual_coefficient``
    is ``C_R(Fn)`` (Fossen 2011, eq. 6.82's decomposition ``C_T = (1+k) C_F + C_R``),
    0 reproduces today's ITTC-only form exactly."""
    displaced_volume = mass / water_density
    added_mass = surge_added_mass_factor * water_density * displaced_volume ** (5.0 / 3.0) / length**2
    linear_coefficient = -(mass + added_mass) / time_constant
    friction = ittc_friction_coefficient(u_r, length, kinematic_viscosity, reynolds_floor)
    quadratic_friction = -0.5 * water_density * wetted_surface * (1.0 + form_factor) * friction
    quadratic_residual = -0.5 * water_density * wetted_surface * residual_coefficient
    sigma = 1.0 - math.tanh(abs(u_r) / crossover_speed)
    return (sigma * linear_coefficient * u_r
            + (1.0 - sigma) * (quadratic_friction + quadratic_residual) * abs(u_r) * u_r)


# RA14 Appendix 1, p. 24: R/Delta, S/V^(2/3), LK/L as cubics in Fn_V, each
# coefficient quartic/cubic in s = L/V^(1/3). Copied with citation from
# `agents-more/60_working/A64_scripts/a64_hull_regime_resistance.py` (A-64,
# transcribed and verified there); not retyped from the paper in this job
# (rule 4).
R_COEF = {
    "A": (-0.0048694, 0.1057838, -0.8432151, 2.8994541, -3.5683179),
    "B": (0.0301221, -0.6562651, 5.2443776, -18.0560600, 22.1656778),
    "C": (-0.0508810, 1.1119496, -8.9006694, 30.6066779, -37.2112557),
    "D": (0.0209140, -0.4573261, 3.6580101, -12.5431467, 15.14353),
}
S_COEF = {
    "A": (0.0197989, -0.2876721, 1.3944044, -2.1175915),
    "B": (-0.1612880, 2.3295698, -11.3511782, 18.0264798),
    "C": (0.4226973, -6.1030996, 29.8350864, -49.7163001),
    "D": (-0.4432887, 6.7794188, -33.1423777, 58.8177152),
}


def _poly(coefficients, x):
    return sum(ci * x ** (len(coefficients) - 1 - i) for i, ci in enumerate(coefficients))


def _cubic_in_fn(coefficients, s, fn):
    a, b, c, d = (_poly(coefficients[k], s) for k in "ABCD")
    return a * fn**3 + b * fn**2 + c * fn + d


RA14_VALIDITY_FN_V = 0.6  # below it C_R is set to 0 (A-64 hidden assumption 3: outside the regression, C_R "not shown small")


def ra14_residual_coefficient_and_area(mass, length, speed):
    """``C_R(Fn)`` and the RA14 wetted-surface option ``S = (S/V^(2/3)) V^(2/3)``,
    Froude-scaled from the 100000 lb basis (RA14 p. 11: ITTC-1957, C_A = 0),
    ``R_COEF``/``S_COEF`` from Appendix 1 p. 24. The scaling procedure (same s,
    same Fn_V between Grethe and the standard craft) is the construction of
    A-64 (`20_sources/monohull_structure_for_grethe.md` hidden assumption 3),
    not an equation printed in the paper. Below ``Fn_V = 0.6`` the regression
    does not cover the hull and ``C_R`` is clamped to 0 (A-64 hidden
    assumption 3, same place)."""
    vol = mass / RHO
    s = length / vol ** (1.0 / 3.0)
    fn = speed / math.sqrt(G * vol ** (1.0 / 3.0))
    r_over_delta = _cubic_in_fn(R_COEF, s, fn)
    s_coefficient = _cubic_in_fn(S_COEF, s, fn)
    rho_s, nu_s, m_s = 1026.0, 1.1907e-6, 100000 * 0.45359237
    vol_s = m_s / rho_s
    lam = (vol_s / vol) ** (1.0 / 3.0)
    u_s, l_s = speed * math.sqrt(lam), length * lam
    area_s = s_coefficient * vol_s ** (2.0 / 3.0)
    rt_s = r_over_delta * m_s * G
    ct_s = rt_s / (0.5 * rho_s * area_s * u_s**2)
    cf_s = ITTC_FRICTION_FACTOR / (math.log10(max(u_s * l_s / nu_s, 1.0)) - ITTC_LOG10_OFFSET) ** 2
    c_r = ct_s - cf_s  # the residual, after the standard craft's own friction (RA14 p. 11 basis)
    if fn < RA14_VALIDITY_FN_V:
        c_r = 0.0
    area = s_coefficient * vol ** (2.0 / 3.0)
    return c_r, area, fn


# ---------------------------------------------------------------------------
# Linear coupled sway-yaw damper, zero cross terms (E-101 Q2 a)
# ---------------------------------------------------------------------------

def linear_coupled_sway_yaw(mass_matrix, time_sway, time_yaw, nu_r):
    """``Y_v = -M22 / T_sway``, ``N_r = -M66 / T_yaw`` (as the library's own
    surface damping form, otter.m 203, 207); ``Y_r = N_v = 0``, marked
    unidentified (E-101 Q2 a)."""
    y_v, n_r = -mass_matrix[1, 1] / time_sway, -mass_matrix[5, 5] / time_yaw
    y_r = n_v = 0.0
    tau_y = y_v * nu_r[1] + y_r * nu_r[5]
    tau_n = n_v * nu_r[1] + n_r * nu_r[5]
    return y_v, y_r, n_v, n_r, tau_y, tau_n


# ---------------------------------------------------------------------------
# Grethe's parameter set (A-43 Part 2, A-64 Section 1.1, E-101)
# ---------------------------------------------------------------------------

def parameter_set():
    return {
        "gravity": {"value": G, "unit": "m/s^2", "kind": "estimate",
                   "place": "grethe_mariner5.yaml hydrostatics.gravity (\"standard\"), A-43 header"},
        "water_density": {"value": RHO, "unit": "kg/m^3", "kind": "published",
                          "place": "A-43 M5: standard sea water"},
        "kinematic_viscosity": {"value": NU_WATER, "unit": "m^2/s", "kind": "published",
                                "place": "Fossen 2011 p. 125, below eq. 6.85 (20 degC)"},
        "length": {"value": 5.2, "unit": "m", "kind": "document",
                  "place": "A-43 M1: F003 p2, M001 p7 (overall length)"},
        "beam": {"value": 2.15, "unit": "m", "kind": "document",
                "place": "A-43 M2: F003 p2, M001 p7, D002 p1"},
        "waterline_length": {"value": 4.63, "unit": "m", "kind": "estimate",
                             "place": "A-64 Section 1.1: D002 p1 drawing measurement at T = 0.30 m "
                                      "(sensitivity case for the unmeasured projected chine length L_P, "
                                      "A-64 Section 1.2 footnote 7)"},
        "waterline_beam": {"value": 1.97, "unit": "m", "kind": "estimate",
                           "place": "A-64 Section 1.1: D002 p1 drawing measurement at T = 0.30 m"},
        "hull_block_coefficient": {"value": 0.234, "unit": "1", "kind": "derived",
                                   "place": "A-43 M4: 806 / (1025 * 5.2 * 2.15 * 0.3), a mass knob not hull fullness"},
        "hull_waterplane_coefficient": {"value": 0.8, "unit": "1", "kind": "estimate", "place": "A-43 M14"},
        "hull_mass": {"value": 806.0, "unit": "kg", "kind": "document", "place": "A-43 M7: F003 p2 (FAT, +/- 20 kg)"},
        "payload_mass": {"value": 0.0, "unit": "kg", "kind": "gap",
                         "place": "A-43 M8: the 2026 field payload (USBL mount, battery packs) is not weighed"},
        "hull_center_of_gravity": {"value": [0.0, 0.0, 0.025], "unit": "m", "kind": "estimate", "place": "A-43 M9"},
        "payload_position": {"value": [0.0, 0.0, 0.0], "unit": "m", "kind": "gap", "place": "A-43 M10 (unused, payload_mass 0)"},
        "radii_of_gyration": {"value": [0.35, 0.25, 0.25], "unit": "1", "kind": "estimate", "place": "A-43 M12"},
        "added_mass_coefficients": {"value": [-1.0, -1.5, -1.0, -0.2, -0.8, -1.7], "unit": "1", "kind": "estimate",
                                    "place": "A-43 M11 (the workspace .rppws value, open: library yaw is -1.2)"},
        "added_mass_form": {"value": "scaled_derivatives", "unit": "-", "kind": "selector", "place": "monohull.py 99"},
        "added_mass_mass_basis": {"value": "displaced", "unit": "-", "kind": "owner decision",
                                  "place": "E-101 Q5 a (not yet a library option, A-64 Section 3.1 gap 6)"},
        "coriolis_form": {"value": "kirchhoff_full", "unit": "-", "kind": "selector", "place": "monohull.py 100"},
        "current_form": {"value": "none", "unit": "-", "kind": "test convention",
                         "place": "still water: the G1 references below hold current fixed, not the form itself"},
        "site_form": {"value": "given", "unit": "-", "kind": "selector", "place": "monohull.py 84"},
        "hull_count": {"value": 1, "unit": "-", "kind": "convention", "place": "A-43 M15: monohull"},
        "reference_point": {"value": [-2.24, 0.0, 0.0], "unit": "m", "kind": "owner decision",
                            "place": "E-101 Q4 a: CO at the midpoint of the waterline length, on the waterline"},
        "longitudinal_inertia_factor": {"value": 1.0, "unit": "1", "kind": "estimate", "place": "A-43 M17"},
        "longitudinal_center_of_flotation": {"value": -0.3, "unit": "m", "kind": "estimate", "place": "A-43 M18"},
        "wetted_surface_method": {"value": "regression_table", "unit": "-", "kind": "owner decision",
                                  "place": "A-64 Section 3.3 R11: not yet a library option (today: computed/given only)"},
        "time_constants": {"value": [8.0, 16.0], "unit": "s", "kind": "estimate",
                           "place": "A-43 M24: [T_sway, T_yaw] (surface damping form)"},
        "damping_ratios": {"value": [0.3, 0.2, 0.5], "unit": "1", "kind": "estimate", "place": "A-43 M21"},
        "yaw_damping_nonlinearity": {"value": 10.0, "unit": "s/rad", "kind": "published (other vehicle)",
                                    "place": "A-43 M25: the Otter's value, not Grethe's (OTTER:240)"},
        "surge_resistance": {"value": "ittc_residual", "unit": "-", "kind": "owner decision",
                             "place": "E-101 Q1 a (not yet a library option: today's SurfaceHullLoads has none/ittc only)"},
        "time_constant": {"value": 15.0, "unit": "s", "kind": "estimate", "place": "A-43 M26 (= M20 T1)"},
        "form_factor": {"value": 0.1, "unit": "1", "kind": "published", "place": "Fossen 2011 p. 125; A-43 M27"},
        "crossover_speed": {"value": 2.0, "unit": "m/s", "kind": "published (MSS constant)", "place": "A-43 M28"},
        "surge_added_mass_factor": {"value": 2.7, "unit": "1", "kind": "published", "place": "addedMassSurge.m 33; A-43"},
        "ittc_reynolds_floor": {"value": 1.0e5, "unit": "1", "kind": "published + owner default", "place": "A-43 M29"},
        "manoeuvring_damping": {"value": "linear_coupled", "unit": "-", "kind": "owner decision",
                                "place": "A-64 Section 3.4 G0 set (not yet a library option)"},
        "sway_yaw_derivatives": {"value": None, "unit": "N s/m, N m s", "kind": "derived + unidentified",
                                 "place": "Y_v, N_r from time_constants (library's own formula, D_ii = M_ii / T_i); "
                                          "Y_r = N_v = 0 (E-101 Q2 a) — numeric values in grethe_damping.csv"},
        "modulus_coefficients": {"value": 0.0, "unit": "-", "kind": "unidentified", "place": "A-64 Section 3.3 (placeholder, not identified)"},
        "cross_flow": {"value": "hoerner_strips", "unit": "-", "kind": "selector", "place": "monohull.py / A-64 Section 3.4"},
    }


def main():
    out = HERE
    values = parameter_set()
    (out / "grethe_parameters.json").write_text(json.dumps({
        "vehicle": "Grethe (NTNU Mariner 5 USV, Maritime Robotics build 40402, Pioner 17 ft hull)",
        "source": "agents-more/60_working/2026-10-07_A43_vehicle_parameter_dossiers.md Part 2; "
                 "agents-more/20_sources/monohull_structure_for_grethe.md Section 1.1, 3.3; "
                 "agents-more/80_owner/inbox/E-101_a64_grethe_monohull_choices_remus_density_at_repin.md",
        "no_august_data": True,
        "parameters": values,
    }, indent=2) + "\n")

    # -- rigid body + added mass: two payload masses, two added-mass bases ---------------------
    hull_cg = np.array(values["hull_center_of_gravity"]["value"])
    payload_position = np.array(values["payload_position"]["value"])
    radii = values["radii_of_gyration"]["value"]
    length, beam = values["length"]["value"], values["beam"]["value"]
    coefficients = values["added_mass_coefficients"]["value"]
    rows = []
    for payload_mass in (0.0, 300.0):  # the bare hull and the owner's "+300 kg" check case (E-64 row 7)
        mass, r_g, inertia = rigid_body_hull_with_payload(
            values["hull_mass"]["value"], payload_mass, hull_cg, payload_position, radii, beam, length)
        for basis_name, basis_mass in (("hull", values["hull_mass"]["value"]), ("displaced", mass)):
            m_a = scaled_added_mass_matrix(basis_mass, length, RHO, inertia, coefficients)
            rows.append({
                "payload_mass": payload_mass, "added_mass_mass_basis": basis_name,
                "mass": mass, "x_g": r_g[0], "y_g": r_g[1], "z_g": r_g[2],
                "I11": inertia[0, 0], "I22": inertia[1, 1], "I33": inertia[2, 2],
                **{f"M_A_{i+1}{i+1}": m_a[i, i] for i in range(6)},  # the added_mass_matrix diagonal (M_A = -diag(derivatives), positive)
            })
    with open(out / "grethe_rigid_body_added_mass.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    # -- Coriolis: seeded 6-DOF states, bare hull, hull-mass added-mass basis ------------------
    mass0, r_g0, inertia0 = rigid_body_hull_with_payload(
        values["hull_mass"]["value"], 0.0, hull_cg, payload_position, radii, beam, length)
    m_a_hull = scaled_added_mass_matrix(values["hull_mass"]["value"], length, RHO, inertia0, coefficients)
    rng = np.random.default_rng(20261009)
    rows = []
    for _ in range(30):
        nu = np.r_[rng.uniform(-3.0, 3.0, 1), rng.uniform(-1.0, 1.0, 1), rng.uniform(-0.3, 0.3, 1),
                   rng.uniform(-0.5, 0.5, 2), rng.uniform(-1.0, 1.0, 1)]
        nu_r = nu  # still water (current_form "none")
        c_rb = c_rb_co(mass0, inertia0, r_g0, nu)
        c_a_full = c_a(m_a_hull, nu_r, stabilize=False)
        c_a_munk_removed = c_a(m_a_hull, nu_r, stabilize=True)
        row = {f"nu_{i+1:02d}": nu[i] for i in range(6)}
        row.update({f"C_RB_{i+1:02d}": c_rb.flatten(order="F")[i] for i in range(36)})
        row.update({f"C_A_full_{i+1:02d}": c_a_full.flatten(order="F")[i] for i in range(36)})
        row.update({f"C_A_munk_removed_{i+1:02d}": c_a_munk_removed.flatten(order="F")[i] for i in range(36)})
        rows.append(row)
    with open(out / "grethe_coriolis.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    # -- restoring: bare hull and +300 kg, CO at the waterline midpoint ------------------------
    rows = []
    for payload_mass in (0.0, 300.0):
        mass, r_g, _ = rigid_body_hull_with_payload(
            values["hull_mass"]["value"], payload_mass, hull_cg, payload_position, radii, beam, length)
        restoring = surface_restoring(  # center_of_gravity is the rigid body's own combined r_g (a vehicle
            # coupling by name, Monohull wires "center_of_gravity" from the rigid-body block, not hull_cg alone)
            length=length, beam=beam, block_coefficient=values["hull_block_coefficient"]["value"],
            waterplane_coefficient=values["hull_waterplane_coefficient"]["value"], mass=mass, water_density=RHO,
            gravity=G, center_of_gravity=r_g,
            longitudinal_inertia_factor=values["longitudinal_inertia_factor"]["value"],
            longitudinal_center_of_flotation=values["longitudinal_center_of_flotation"]["value"],
            reference_point=values["reference_point"]["value"], hull_count=1)
        rows.append({
            "payload_mass": payload_mass, "draft": restoring["draft"],
            "displaced_volume": restoring["displaced_volume"],
            "wetted_surface_mumford": restoring["wetted_surface_mumford"],
            "transverse_metacentric_height": restoring["transverse_metacentric_height"],
            "longitudinal_metacentric_height": restoring["longitudinal_metacentric_height"],
            "G_33": restoring["G"][2, 2], "G_44": restoring["G"][3, 3], "G_55": restoring["G"][4, 4],
        })
    with open(out / "grethe_restoring.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    # -- cross-flow strips: seeded nu_r, bare-hull draft; both today's basis (overall
    # length/beam, what single_hull() actually feeds the strip today) and the future
    # basis (waterline length/beam, once R11's selector reads them instead) ------------------
    draft0 = rows[0]["draft"]
    rng = np.random.default_rng(20261010)
    cf_rows = []
    for _ in range(30):
        nu_r = np.r_[rng.uniform(-3.0, 3.0, 1), rng.uniform(-1.0, 1.0, 1), rng.uniform(-0.3, 0.3, 1),
                    rng.uniform(-0.5, 0.5, 2), rng.uniform(-1.0, 1.0, 1)]
        tau_overall = cross_flow_hoerner(nu_r, length, beam, draft0, RHO)
        tau_waterline = cross_flow_hoerner(nu_r, values["waterline_length"]["value"],
                                           values["waterline_beam"]["value"], draft0, RHO)
        row = {f"nu_r_{i+1:02d}": nu_r[i] for i in range(6)}
        row.update({f"tau_overall_{i+1:02d}": tau_overall[i] for i in range(6)})
        row.update({f"tau_waterline_{i+1:02d}": tau_waterline[i] for i in range(6)})
        cf_rows.append(row)
    with open(out / "grethe_cross_flow.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(cf_rows[0]))
        writer.writeheader()
        writer.writerows(cf_rows)

    # -- surge resistance: 0..3.2 m/s, both L sensitivity cases, ITTC-only and ITTC+RA14 -------
    mass0_only_hull = values["hull_mass"]["value"]
    sr_rows = []
    for u in np.round(np.arange(0.0, 3.21, 0.1), 2):
        row = {"u": u}
        for label, length_case in (("L463", 4.63), ("L52", 5.2)):
            wetted_mumford = surface_restoring(
                length=length_case, beam=values["beam"]["value"] if label == "L52" else values["waterline_beam"]["value"],
                block_coefficient=values["hull_block_coefficient"]["value"],
                waterplane_coefficient=values["hull_waterplane_coefficient"]["value"], mass=mass0_only_hull,
                water_density=RHO, gravity=G, center_of_gravity=hull_cg,
                longitudinal_inertia_factor=values["longitudinal_inertia_factor"]["value"],
                longitudinal_center_of_flotation=values["longitudinal_center_of_flotation"]["value"],
                reference_point=values["reference_point"]["value"], hull_count=1)["wetted_surface_mumford"]
            c_r, s_ra14, fn_v = ra14_residual_coefficient_and_area(mass0_only_hull, length_case, max(u, 1e-6))
            ittc_only = ittc_surge_force(
                u, mass=mass0_only_hull, length=length_case, water_density=RHO, wetted_surface=wetted_mumford,
                time_constant=values["time_constant"]["value"], form_factor=values["form_factor"]["value"],
                crossover_speed=values["crossover_speed"]["value"], kinematic_viscosity=NU_WATER,
                surge_added_mass_factor=values["surge_added_mass_factor"]["value"],
                reynolds_floor=values["ittc_reynolds_floor"]["value"], residual_coefficient=0.0)
            ittc_residual = ittc_surge_force(
                u, mass=mass0_only_hull, length=length_case, water_density=RHO, wetted_surface=wetted_mumford,
                time_constant=values["time_constant"]["value"], form_factor=values["form_factor"]["value"],
                crossover_speed=values["crossover_speed"]["value"], kinematic_viscosity=NU_WATER,
                surge_added_mass_factor=values["surge_added_mass_factor"]["value"],
                reynolds_floor=values["ittc_reynolds_floor"]["value"], residual_coefficient=c_r if u > 0 else 0.0)
            row.update({f"wetted_surface_mumford_{label}": wetted_mumford, f"C_R_RA14_{label}": c_r,
                       f"S_RA14_{label}": s_ra14, f"Fn_V_{label}": fn_v,
                       f"X_ittc_only_{label}": ittc_only, f"X_ittc_residual_{label}": ittc_residual})
        sr_rows.append(row)
    with open(out / "grethe_surge_resistance.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(sr_rows[0]))
        writer.writeheader()
        writer.writerows(sr_rows)

    # -- linear coupled sway-yaw damper: bare hull, hull-mass added-mass basis -----------------
    mass_matrix = m_rb(mass0, inertia0, r_g0) + m_a_hull
    rng = np.random.default_rng(20261011)
    d_rows = []
    y_v, y_r, n_v, n_r, _, _ = linear_coupled_sway_yaw(
        mass_matrix, values["time_constants"]["value"][0], values["time_constants"]["value"][1], np.zeros(6))
    for _ in range(20):
        nu_r = np.r_[0.0, rng.uniform(-1.0, 1.0, 1), 0.0, 0.0, 0.0, rng.uniform(-1.0, 1.0, 1)]
        _, _, _, _, tau_y, tau_n = linear_coupled_sway_yaw(
            mass_matrix, values["time_constants"]["value"][0], values["time_constants"]["value"][1], nu_r)
        d_rows.append({"v_r": nu_r[1], "r_r": nu_r[5], "Y_v": y_v, "Y_r": y_r, "N_v": n_v, "N_r": n_r,
                       "tau_Y": tau_y, "tau_N": tau_n})
    with open(out / "grethe_damping.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(d_rows[0]))
        writer.writeheader()
        writer.writerows(d_rows)

    print("wrote: grethe_parameters.json, grethe_rigid_body_added_mass.csv, grethe_coriolis.csv, "
         "grethe_restoring.csv, grethe_cross_flow.csv, grethe_surge_resistance.csv, grethe_damping.csv")
    print(f"Y_v = {y_v:.6f} N s/m (from T_sway=8s), N_r = {n_r:.6f} N m s (from T_yaw=16s), Y_r = N_v = 0 (unidentified)")


if __name__ == "__main__":
    main()
