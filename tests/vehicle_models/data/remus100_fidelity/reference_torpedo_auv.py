"""Reference functions for the higher-fidelity forms of ``SpheroidAuv``,
``AuvHullLoads`` and ``Fin``, and a whole torpedo-AUV built from them, in
numpy -- the second implementation the fidelity-option gates compare the
plugins' CasADi graphs against, and the generator of the frozen turn and
stability references beside this file.

Every function takes a parameter set (geometry and coefficients in, forces
or derivatives out); every number comes from the caller's set (a named
part's parameter file), none is written here. The parameter names are the
plugins' own. Frames: BODY FRD about the body origin CO (x forward, y
starboard, z down); ``nu_r = [u, v, w, p, q, r]`` the velocity relative to
the water, SI units.

Forms (named by their physics; the published reference model's choice of each
mechanism is one form):

* body lift: ``vertical_plane_wing`` -- lift on the vertical-plane angle of
  attack only, no sway force, no moment; ``point_force_both_planes`` -- the
  slender-body lift ``F_y = -1/2 rho d^2 c_ydb u v``, ``F_z = -1/2 rho d^2 c_ydb
  u w`` (Prestero 2001, eqs. 4.27-4.29, 4.33-4.34, pp. 30-31; the sway half by
  the hull's symmetry of revolution) applied as one point force at the pressure
  centre ``r_cp = (x_cp, 0, 0)``: ``tau = [F; r_cp x F]`` (Fossen 2011,
  eq. 12.226, p. 400). Hence ``M_uwl = +1/2 rho d^2 c_ydb x_cp = -N_uvl``;
  Prestero's eq. 4.36 prints ``M_uwl`` with the opposite sign, his own
  linearised Table D.1 (p. 112) uses ``M_wl / Z_wl = -x_cp``, the sign used here.
* body-lift slope: ``slender_body_hoerner`` ``c_ydb = (l/d) c_yb 180/pi``
  with ``c_yb`` per degree (Prestero 2001, eqs. 4.30-4.32, p. 30, after
  Hoerner 1965, p. 13-3, not read; stated valid for 6.7 <= l/d <= 10);
  ``given``.
* pressure centre: ``fraction_of_length_from_nose`` ``x_cp = x_nose - f l``
  (Prestero 2001, eq. 4.35, p. 31, f = 0.65; Hoerner's range 0.6-0.7 l as
  Prestero p. 31 gives it), ``x_nose = l/2`` for a hull centred on the
  body origin; ``given``.
* hull drag: ``wing_drag`` -- parasitic plus induced drag of a wing on the
  vertical-plane angle of attack (the drag half of the wing block, Fossen 2025
  after Beard and McLain 2012, not read); ``quadratic_axial``
  ``X = -1/2 rho c_d A_f u|u|``, ``A_f = pi d^2 / 4`` (Prestero 2001,
  eqs. 4.5-4.6, p. 25).
* roll damping: ``linear`` (the time-constant damping alone);
  ``linear_plus_quadratic`` adds ``K = K_p|p| p|p|`` (Prestero 2001,
  eqs. 4.10-4.11, p. 27), the coefficient given only (his eq. 4.11 estimate
  is not reproduced from his own inputs).
* added-mass Coriolis: ``kirchhoff_full`` (every term, Fossen 2011,
  eq. 6.43, p. 120) or ``munk_couplings_removed`` (the couplings that carry the
  Munk moment of a slender body, Fossen 2011, eq. 6.52, p. 121, set to zero:
  a shortcut of the published reference model).
* fin lift slope: ``given``; ``low_aspect_ratio``
  ``c_La = [1/(2 a pi) + 1/(pi AR_e)]^-1``, ``AR_e = k_AR b^2 / S`` (Prestero
  2001, eqs. 4.38-4.39, p. 31; Hoerner 1965 p. 3-2 and Whicker and
  Fehlner 1958 behind it, not read), per fin, ``S`` the fin's own area, ``k_AR`` the
  ``aspect_ratio_factor`` (2).
* fins: inflow ``translational`` or ``rigid_point`` (Prestero 2001,
  eq. 4.40, p. 32), flow angle ``none`` (beta = 0, the speed in the section
  plane) or ``small_angle`` (Prestero 2001, eqs. 4.41-4.43, pp. 32-33), section
  ``lift_times_angle`` (``C_L = a alpha``, ``C_D = C_L alpha``); a positive
  deflection lifts along the fin's ``lift_axis``.

References
----------
Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics and Motion
    Control. Wiley. Ch. 2, eq. 2.18, p. 22; ch. 3, eq. 3.24, p. 49 and
    eq. 3.34, p. 50; ch. 4, eq. 4.6, p. 60; ch. 6, eq. 6.43, p. 120,
    eq. 6.52, p. 121, eqs. 6.76-6.81, p. 125; ch. 12, eq. 12.226, p. 400.
Prestero, T. (2001). Verification of a six-degree of freedom simulation
    model for the REMUS autonomous underwater vehicle. MSc thesis, MIT/WHOI.
    Ch. 4, eqs. 4.5-4.6, 4.10-4.11, 4.27-4.43, pp. 25-33.
Lamb, H. (1932). Hydrodynamics. Cambridge University Press (the ideal-fluid
    inertia coefficients of a prolate spheroid; page not located).
DNV (2010). Recommended Practice DNV-RP-C205, Fig. 6-6 (the drag coefficient
    of a circular cylinder against Reynolds number; the table is typed below).

Author:    Enio Krizman
Date:      2026-10-13
"""

import numpy as np

from more_transformations.ecef_ned_transforms import ECEFNEDtransform
from more_transformations.euler_ned_body_transforms import EulerNEDBodyTransforms
from more_transformations.matrix_transforms import MatrixTransforms

skew = MatrixTransforms.skew


def wrench(force, arm):
    """``tau = [f; r x f]`` (Fossen 2011, eq. 12.226, p. 400)."""
    force = np.asarray(force, float)
    return np.concatenate([force, skew(np.asarray(arm, float)) @ force])


# ---------------------------------------------------------------------------------------------
# Derived quantities (each is one option of a ``*_method``)
# ---------------------------------------------------------------------------------------------
def slender_body_lift_slope(length, diameter, lift_slope_per_degree):
    """c_ydb = (l/d) c_yb 180/pi, 1/rad (Prestero 2001, eqs. 4.30-4.32, p. 30)."""
    return (length / diameter) * lift_slope_per_degree * 180.0 / np.pi


def pressure_centre_from_nose(nose_x, length, fraction):
    """x_cp = x_nose - f l, BODY (Prestero 2001, eq. 4.35, p. 31)."""
    return nose_x - fraction * length


def low_aspect_ratio_lift_slope(span, area, efficiency, factor):
    """c_La per fin from its span, area, empirical factor and ``aspect_ratio_factor`` (Prestero 2001, eqs. 4.38-4.39, p. 31)."""
    effective_aspect_ratio = factor * span ** 2 / area                       # (Prestero 2001, eq. 4.39, p. 31)
    return 1.0 / (1.0 / (2.0 * efficiency * np.pi) + 1.0 / (np.pi * effective_aspect_ratio))   # (eq. 4.38)


# ---------------------------------------------------------------------------------------------
# Rigid body, added mass, Coriolis
# ---------------------------------------------------------------------------------------------
def spheroid_inertia_about_cg(mass, a, b):
    """Principal inertia of a homogeneous prolate spheroid about its centre, ``diag(2/5 m b^2, 1/5 m (a^2 + b^2), 1/5 m (a^2 + b^2))``."""
    return np.diag([2.0 / 5.0 * mass * b ** 2, 1.0 / 5.0 * mass * (a ** 2 + b ** 2), 1.0 / 5.0 * mass * (a ** 2 + b ** 2)])


def rigid_body(mass, inertia_cg, r_bg, omega):
    """M_RB, C_RB about the CO from the mass and the inertia about the CG
    (``H^T diag(m I3, I_g) H``, ``H^T diag(m S(w), -S(I_g w)) H``; H matrix Fossen 2011, eq. 3.24, p. 49)."""
    m_cg = np.zeros((6, 6))
    m_cg[:3, :3] = mass * np.eye(3)
    m_cg[3:, 3:] = inertia_cg
    c_cg = np.zeros((6, 6))
    c_cg[:3, :3] = mass * skew(omega)
    c_cg[3:, 3:] = -skew(inertia_cg @ omega)
    h = MatrixTransforms.H_matrix(np.asarray(r_bg, float))                 # (Fossen 2011, eq. 3.24, p. 49)
    return h.T @ m_cg @ h, h.T @ c_cg @ h


def inertia_about_co(mass, inertia_cg, r_bg):
    """I_o = I_g - m S(r_g)^2, the parallel-axis form of M_RB's lower block (Fossen 2011, eq. 3.24, p. 49)."""
    s = skew(np.asarray(r_bg, float))
    return inertia_cg - mass * s @ s


def lamb_spheroid_added_mass(a, b, rho, r44):
    """Lamb's added mass of a prolate spheroid (Lamb 1932: the coefficients k1, k2 and k' of the ideal fluid), with the roll added inertia ``r44`` times that of the displaced fluid."""
    m = 4.0 / 3.0 * np.pi * rho * a * b ** 2
    ix = 2.0 / 5.0 * m * b ** 2
    iy = 1.0 / 5.0 * m * (a ** 2 + b ** 2)
    e = np.sqrt(1.0 - (b / a) ** 2)
    alpha0 = (2.0 * (1.0 - e ** 2) / e ** 3) * (0.5 * np.log((1 + e) / (1 - e)) - e)
    beta0 = 1.0 / e ** 2 - (1.0 - e ** 2) / (2.0 * e ** 3) * np.log((1 + e) / (1 - e))
    k1 = alpha0 / (2.0 - alpha0)
    k2 = beta0 / (2.0 - beta0)
    kp = e ** 4 * (beta0 - alpha0) / ((2.0 - e ** 2) * (2.0 * e ** 2 - (2.0 - e ** 2) * (beta0 - alpha0)))
    return np.diag([m * k1, m * k2, m * k2, r44 * ix, kp * iy, kp * iy])


MUNK_ENTRIES = ((4, 2), (2, 4), (4, 0), (0, 4), (5, 0), (0, 5), (5, 1), (1, 5))   # 0-based: pitch-heave, pitch-surge, yaw-surge, yaw-sway and transposes


def added_mass_coriolis(form, m_a, nu_r):
    """C_A(nu_r) with every term kept (Fossen 2011, eq. 6.43, p. 120), or with the
    Munk-moment couplings removed (the pitch-heave, pitch-surge, yaw-surge and yaw-sway entries and their transposes)."""
    m = 0.5 * (m_a + m_a.T)
    p1 = m[:3, :3] @ nu_r[:3] + m[:3, 3:] @ nu_r[3:]
    p2 = m[3:, :3] @ nu_r[:3] + m[3:, 3:] @ nu_r[3:]
    c = np.zeros((6, 6))
    c[:3, 3:] = -skew(p1)
    c[3:, :3] = -skew(p1)
    c[3:, 3:] = -skew(p2)
    if form == "kirchhoff_full":
        return c
    if form == "munk_couplings_removed":
        for i, j in MUNK_ENTRIES:
            c[i, j] = 0.0
        return c
    raise ValueError(f"coriolis_form {form!r}")


# ---------------------------------------------------------------------------------------------
# Hull loads
# ---------------------------------------------------------------------------------------------
def time_constant_damping(m_total, weight, r_bg, r_bb, time_constants, damping_ratios, speed):
    """Time-constant damping, submerged form (Fossen 2011, eqs. 6.76-6.80, p. 125), with the surge term faded by ``exp(-3 U_r)``."""
    t1, t2, t6 = time_constants
    z4, z5 = damping_ratios
    w4 = np.sqrt(weight * (r_bg[2] - r_bb[2]) / m_total[3, 3])
    w5 = np.sqrt(weight * (r_bg[2] - r_bb[2]) / m_total[4, 4])
    d = np.diag([m_total[0, 0] / t1, m_total[1, 1] / t2, m_total[2, 2] / t2,
                 m_total[3, 3] * 2 * z4 * w4, m_total[4, 4] * 2 * z5 * w5, m_total[5, 5] / t6])
    d[0, 0] *= np.exp(-3.0 * speed)
    return d


def quadratic_roll_damping(k_pp, p):
    """K = K_p|p| p|p| (Prestero 2001, eqs. 4.10-4.11, p. 27)."""
    tau = np.zeros(6)
    tau[3] = k_pp * p * abs(p)
    return tau


def _wing_coefficients(span, area, cd0, oswald, alpha):
    aspect_ratio = span ** 2 / area
    cl_alpha = np.pi * aspect_ratio / (1 + np.sqrt(1 + (aspect_ratio / 2) ** 2))
    cl = cl_alpha * alpha
    cd = cd0 + cl ** 2 / (np.pi * oswald * aspect_ratio)
    return cl, cd


def vertical_plane_wing_lift(nu_r, rho, span, area, cd0, oswald):
    """The lift half of the wing block: on alpha = atan2(w, u) only; no sway, no moment."""
    alpha = np.arctan2(nu_r[2], nu_r[0])
    speed2 = nu_r[0] ** 2 + nu_r[1] ** 2 + nu_r[2] ** 2
    cl, _ = _wing_coefficients(span, area, cd0, oswald, alpha)
    lift = 0.5 * rho * speed2 * area * cl
    return np.array([np.sin(alpha) * lift, 0, -np.cos(alpha) * lift, 0, 0, 0])


def wing_drag(nu_r, rho, span, area, cd0, oswald):
    """The drag half of the wing block (parasitic plus induced on the vertical-plane alpha)."""
    alpha = np.arctan2(nu_r[2], nu_r[0])
    speed2 = nu_r[0] ** 2 + nu_r[1] ** 2 + nu_r[2] ** 2
    _, cd = _wing_coefficients(span, area, cd0, oswald, alpha)
    drag = 0.5 * rho * speed2 * area * cd
    return np.array([-np.cos(alpha) * drag, 0, -np.sin(alpha) * drag, 0, 0, 0])


def quadratic_axial_drag(nu_r, rho, diameter, drag_coefficient):
    """X = -1/2 rho c_d A_f u|u|, A_f = pi d^2/4 (Prestero 2001, eqs. 4.5-4.6, p. 25)."""
    frontal_area = np.pi * diameter ** 2 / 4.0
    return np.array([-0.5 * rho * drag_coefficient * frontal_area * nu_r[0] * abs(nu_r[0]), 0, 0, 0, 0, 0])


def point_force_body_lift(nu_r, rho, diameter, lift_slope, pressure_centre):
    """Body lift in both planes as one point force at (x_cp, 0, 0) (Prestero 2001, eqs. 4.27-4.29,
    4.33-4.34, pp. 30-31; Fossen 2011, eq. 12.226, p. 400)."""
    k = 0.5 * rho * diameter ** 2 * lift_slope
    u, v, w = nu_r[0], nu_r[1], nu_r[2]
    return wrench([0.0, -k * u * v, -k * u * w], [pressure_centre, 0.0, 0.0])


_CD_DATA = np.array([
    [10211.0405297256, 1.20769], [15543.5423211490, 1.20369], [24434.8681540938, 1.20823],
    [35719.8807222993, 1.21060], [60880.7508493947, 1.21093], [86174.7436190610, 1.20901],
    [118086.431856986, 1.21134], [149262.255516420, 1.21148], [214708.876786124, 1.21171],
    [230940.445891927, 1.20109], [271581.385758833, 1.16918], [297098.185911873, 1.11163],
    [325252.018159382, 1.00926], [361942.848517458, 0.89411], [409798.156149477, 0.70855],
    [475302.392082697, 0.53155], [524734.960427919, 0.40785], [578921.928291835, 0.32683],
    [638292.879987135, 0.28422], [744046.366191920, 0.29711], [853237.628592033, 0.32280],
    [1068685.70450627, 0.38909], [1393381.37138896, 0.47033], [1831850.38455424, 0.53878],
    [2258395.50114555, 0.58372], [2899011.91287642, 0.62868], [3663279.03082383, 0.64803],
    [4937305.49552761, 0.67808]])
_KAPPA_SUB = np.array([[2, 0.58], [5, 0.62], [10, 0.68], [20, 0.74], [40, 0.82], [50, 0.87], [100, 0.98]])
_KAPPA_SUP = np.array([[2, 0.80], [5, 0.80], [10, 0.82], [20, 0.90], [40, 0.98], [50, 0.99], [100, 1.00]])


def cylinder_strip_cross_flow(nu_r, rho, length, beam, draft, viscosity, n_strips=20):
    """Strip-theory cross-flow drag (Fossen 2011, eqs. 6.91-6.92, p. 127) with the circular-cylinder section drag against Reynolds number and the slenderness factor kappa."""
    re = np.sqrt(nu_r[1] ** 2 + nu_r[2] ** 2) * beam / viscosity
    cd = np.interp(re, _CD_DATA[:, 0], _CD_DATA[:, 1])                      # clamped at both ends, 84-90
    table = _KAPPA_SUB if re < 2e5 else _KAPPA_SUP
    cd2 = cd * np.interp(length / beam, table[:, 0], table[:, 1])
    dx = length / n_strips
    yh = zh = mh = nh = 0.0
    for i in range(n_strips):
        x = -length / 2 + (i + 0.5) * dx
        uh = abs(nu_r[1] + x * nu_r[5]) * (nu_r[1] + x * nu_r[5])
        uv = abs(nu_r[2] + x * nu_r[4]) * (nu_r[2] + x * nu_r[4])
        yh -= 0.5 * rho * draft * cd2 * uh * dx
        zh -= 0.5 * rho * draft * cd2 * uv * dx
        mh -= 0.5 * rho * draft * cd2 * x * uv * dx
        nh -= 0.5 * rho * draft * cd2 * x * uh * dx
    return np.array([0, yh, zh, 0, mh, nh])


# ---------------------------------------------------------------------------------------------
# Actuators
# ---------------------------------------------------------------------------------------------
def lifting_fin(nu_r, rho, deflection, fin):
    """One fin, BODY FRD: ``fin`` holds the plugin's names (``fin_position``, ``chord_axis``,
    ``lift_axis``, ``fin_area``, ``lift_slope``, ``inflow``, ``flow_angle``).

    v_fin = nu_r[:3] (translational) or nu_r[:3] + omega x r (rigid_point; Prestero 2001, eq. 4.40,
    p. 32); u_c = v_fin . chord, u_n = v_fin . lift_axis; flow angle none: beta = 0, s^2 = u_c^2 + u_n^2
    (the speed in the section plane); small_angle: beta = u_n / u_c, s^2 = u_c^2 (Prestero 2001, eqs. 4.41-4.43,
    pp. 32-33); alpha = delta - beta; C_L = c_La alpha, C_D = c_La alpha^2;
    F = 1/2 rho S s^2 (C_L lift_axis - C_D chord_axis), tau = [F; r x F]."""
    r = np.asarray(fin["fin_position"], float)
    n = np.asarray(fin["lift_axis"], float)
    c = np.asarray(fin["chord_axis"], float)
    v = nu_r[:3] + (skew(nu_r[3:]) @ r if fin["inflow"] == "rigid_point" else 0.0)   # (Prestero 2001, eq. 4.40)
    u_c, u_n = v @ c, v @ n
    if fin["flow_angle"] == "small_angle":
        beta, s2 = u_n / u_c, u_c ** 2                                       # (Prestero 2001, eqs. 4.41-4.43)
    elif fin["flow_angle"] == "none":
        beta, s2 = 0.0, u_c ** 2 + u_n ** 2
    else:
        raise ValueError(f"flow_angle {fin['flow_angle']!r}")
    alpha = deflection - beta
    cl = fin["lift_slope"] * alpha
    cd = fin["lift_slope"] * alpha ** 2
    return wrench(0.5 * rho * fin["fin_area"] * s2 * (cl * n - cd * c), r)


def thrust_torque_propeller(speed, shaft_speed_rpm, rho, p):
    """Thrust and torque from the printed open-water coefficients (linearised in the advance speed); the hull feels the reaction -scale K to the shaft torque."""
    n_p = np.clip(shaft_speed_rpm, -p["max_shaft_speed"], p["max_shaft_speed"]) / 60.0
    d = p["propeller_diameter"]
    kt0, kq0, ktm, kqm = p["thrust_torque_coefficients"]
    ja_max = p["max_advance_number"]
    va = (1 - p["wake_fraction"]) * speed
    if n_p > 0:
        thrust = rho * d ** 4 * (kt0 * abs(n_p) * n_p + (ktm - kt0) / ja_max * (va / d) * abs(n_p))
        torque = rho * d ** 5 * (kq0 * abs(n_p) * n_p + (kqm - kq0) / ja_max * (va / d) * abs(n_p))
    else:
        thrust = rho * d ** 4 * kt0 * abs(n_p) * n_p
        torque = rho * d ** 5 * kq0 * abs(n_p) * n_p
    tau = np.zeros(6)
    tau[0] = (1 - p["thrust_deduction"]) * thrust
    tau[3] = -p["roll_moment_scale"] * torque                               # reaction to the shaft torque (Prestero 2001, Sec. 4.6.2)
    return tau


def restoring(weight, buoyancy, r, r_bg, r_bb):
    """g(eta) of a submerged body (Fossen 2011, eq. 4.6, p. 60)."""
    w, b = weight, buoyancy
    return np.array([
        -(w - b) * r[2, 0], -(w - b) * r[2, 1], -(w - b) * r[2, 2],
        -(r_bg[1] * w - r_bb[1] * b) * r[2, 2] + (r_bg[2] * w - r_bb[2] * b) * r[2, 1],
        -(r_bg[2] * w - r_bb[2] * b) * r[2, 0] + (r_bg[0] * w - r_bb[0] * b) * r[2, 2],
        -(r_bg[0] * w - r_bb[0] * b) * r[2, 1] + (r_bg[1] * w - r_bb[1] * b) * r[2, 0]])


# ---------------------------------------------------------------------------------------------
# The vehicle
# ---------------------------------------------------------------------------------------------
class TorpedoAuv:
    """A torpedo AUV in still water; state ``x = [N E D phi theta psi u v w p q r]`` (the plugins'
    order), input ``[fin commands in fin order..., shaft speed rpm]``.

    ``part``: a flat set by the plugins' parameter names -- ``SpheroidAuv`` (``length``, ``diameter``,
    ``body_mass``, ``body_center_of_gravity``, ``roll_added_inertia_ratio``, ``coriolis_form``), the
    site (``water_density``, ``latitude``, ``kinematic_viscosity``), ``SubmergedRestoring``
    (``center_of_buoyancy``; buoyancy neutral), ``AuvHullLoads`` (``time_constants``, ``damping_ratios``,
    ``planform_area``, ``parasitic_drag_coefficient``, ``oswald_efficiency``, ``body_lift_form`` and its
    methods, ``hull_drag_form``, ``hull_frontal_drag_coefficient``, ``roll_damping_form``,
    ``roll_quadratic_damping``), ``fins`` (a list of ``Fin`` parameter sets with their ``lift_slope``
    resolved) and ``propeller`` (a ``Propeller`` set)."""

    def __init__(self, part):
        self.p = part
        a, b = part["length"] / 2.0, part["diameter"] / 2.0
        self.mass = part["body_mass"]
        self.inertia_cg = spheroid_inertia_about_cg(self.mass, a, b)
        self.g = ECEFNEDtransform.gravity(part["latitude"], unit="rad")
        self.weight = self.mass * self.g
        self.buoyancy = self.weight
        self.m_a = lamb_spheroid_added_mass(a, b, part["water_density"], part["roll_added_inertia_ratio"])
        if part["body_lift_form"] == "point_force_both_planes":
            self.lift_slope = (slender_body_lift_slope(part["length"], part["diameter"], part["body_lift_slope_per_degree"])
                               if part["body_lift_slope_method"] == "slender_body_hoerner" else part["body_lift_slope"])
            self.pressure_centre = (pressure_centre_from_nose(a, part["length"], part["pressure_centre_fraction"])
                                    if part["pressure_centre_method"] == "fraction_of_length_from_nose"
                                    else part["pressure_centre"])

    def xdot(self, x, u):
        p = self.p
        x = np.asarray(x, float)
        nu_r = x[6:12]                                                     # still water
        speed = np.sqrt(nu_r[0] ** 2 + nu_r[1] ** 2 + nu_r[2] ** 2)
        rho = p["water_density"]
        r_bg, r_bb = np.asarray(p["body_center_of_gravity"], float), np.asarray(p["center_of_buoyancy"], float)
        m_rb, c_rb = rigid_body(self.mass, self.inertia_cg, r_bg, x[9:12])
        m_total = m_rb + self.m_a
        c_a = added_mass_coriolis(p["coriolis_form"], self.m_a, nu_r)
        d = time_constant_damping(m_total, self.weight, r_bg, r_bb, p["time_constants"], p["damping_ratios"], speed)
        wing = (p["span"], p["planform_area"], p["parasitic_drag_coefficient"], p["oswald_efficiency"])
        if p["body_lift_form"] == "vertical_plane_wing":
            lift = vertical_plane_wing_lift(nu_r, rho, *wing)
        else:
            lift = point_force_body_lift(nu_r, rho, p["diameter"], self.lift_slope, self.pressure_centre)
        drag = (wing_drag(nu_r, rho, *wing) if p["hull_drag_form"] == "wing_drag"
                else quadratic_axial_drag(nu_r, rho, p["diameter"], p["hull_frontal_drag_coefficient"]))
        roll = (quadratic_roll_damping(p["roll_quadratic_damping"], nu_r[3])
                if p["roll_damping_form"] == "linear_plus_quadratic" else np.zeros(6))
        cross = cylinder_strip_cross_flow(nu_r, rho, p["length"], p["diameter"], p["diameter"], p["kinematic_viscosity"])
        fins = np.zeros(6)
        for k, fin in enumerate(p["fins"]):
            delta = float(np.clip(u[k], -fin["max_deflection"], fin["max_deflection"]))
            fins += lifting_fin(nu_r, rho, delta, fin)
        propeller = thrust_torque_propeller(speed, u[len(p["fins"])], rho, p["propeller"])
        rot = EulerNEDBodyTransforms.R_fossen(*x[3:6])                     # (Fossen 2011, eq. 2.18, p. 22)
        tau = (fins + propeller + lift + drag + cross + roll - (c_rb + c_a) @ nu_r - d @ nu_r
               - restoring(self.weight, self.buoyancy, rot, r_bg, r_bb))
        nu_dot = np.linalg.solve(m_total, tau)
        eta_dot = np.concatenate([rot @ x[6:9], EulerNEDBodyTransforms.T_euler_fossen(*x[3:5]) @ x[9:12]])
        return np.concatenate([eta_dot, nu_dot])


def rk4_step(f, x, u, h):
    """One fourth-order Runge-Kutta step, the input held over the step."""
    k1 = f(x, u)
    k2 = f(x + h / 2 * k1, u)
    k3 = f(x + h / 2 * k2, u)
    k4 = f(x + h * k3, u)
    return x + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


MOTION = [3, 4, 6, 7, 8, 9, 10, 11]      # phi theta u v w p q r: position and heading are neutral
SWAY_YAW = [7, 11]                       # v, r


def jacobian(f, x, u, eps=1e-7):
    """Central-difference Jacobian of f over the motion states (MOTION order)."""
    a = np.zeros((len(MOTION), len(MOTION)))
    for j, k in enumerate(MOTION):
        dx = np.zeros(len(x))
        dx[k] = eps
        a[:, j] = ((f(x + dx, u) - f(x - dx, u)) / (2 * eps))[MOTION]
    return a


def straight_running(f, u, x0, tolerance=1e-12, iterations=100):
    """The straight-running equilibrium at fixed inputs: xdot = 0 over MOTION, by Newton from x0
    (stable or not)."""
    x = np.asarray(x0, float).copy()
    for _ in range(iterations):
        residual = f(x, u)[MOTION]
        if np.abs(residual).max() < tolerance:
            break
        x[MOTION] -= np.linalg.solve(jacobian(f, x, u), residual)
    return x, np.abs(f(x, u)[MOTION]).max()


def sway_yaw_roots(a):
    """Roots of the sway-yaw block of the linearisation (Fossen 2011, eqs. 12.17-12.21, p. 351)."""
    i = [MOTION.index(k) for k in SWAY_YAW]
    return np.linalg.eigvals(a[np.ix_(i, i)])


def shaft_speed_holding(part, speed, bracket=(300.0, 1500.0), iterations=80):
    """The shaft speed (rpm) at which the straight-running equilibrium of ``part`` at zero fin angles has the surge speed
    ``speed``: bisection on the equilibrium surge speed, which increases with the shaft speed (the part's fins are
    ``part["fins"]``, the shaft speed the last command)."""
    model = TorpedoAuv(part)
    zero = np.zeros(len(part["fins"]) + 1)
    start = np.zeros(12)
    start[6] = speed

    def surge(n):
        command = zero.copy()
        command[-1] = n
        x, _ = straight_running(model.xdot, command, start)
        return x[6] - speed

    low, high = bracket
    if surge(low) > 0 or surge(high) < 0:
        raise ValueError(f"no shaft speed in {bracket} rpm holds {speed} m/s")
    for _ in range(iterations):
        middle = 0.5 * (low + high)
        low, high = (middle, high) if surge(middle) < 0 else (low, middle)
    return 0.5 * (low + high)
