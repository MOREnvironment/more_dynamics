"""Independent references and the parameter set of the research AUV Marie's
one named part, "Marie (estimated)" (a variant named by where its data come
from; the hull form is the floating one, exact at the surface and fully
submerged), and of a second torpedo geometry used only by the generality
check: Marie's body scaled geometrically by ``SCALE`` (every value derived by
similitude), so the same test file runs with only the data swapped. Numpy and
scipy only; the one library import is ``more_transformations`` (the 6x6
transformation H, as every generator takes it). No ``more_dynamics`` import:
this file is a second implementation of the cited equations, written from the
printed formulas, not a copy of the code under test.

Run:  python3 generate_marie_references.py
      reads the hydrodynamic file ``LAUV_marie.mat`` frozen beside the
      vessel-database tests (``tests/vessel_database/data/``, byte-identical
      copy of the MSS file, sha256 pinned in the test); writes
      ``marie_parameters.json``, ``marie_floating_reference.csv``,
      ``scaled_torpedo_parameters.json`` and ``scaled_torpedo_floating_reference.csv``
      beside this file. A re-run is checked with ``cmp`` against the committed
      files (deterministic arithmetic; no environment variable, no MATLAB).

Every value of ``marie_parameters.json`` has an entry in ``PRIMITIVES`` (typed once,
with its source) or is computed here from the primitives and the file (``kind``
"derived" or "bem"); nothing is typed twice.

Kinds, in the order of what they can be trusted for: ``published`` (a number
printed in a cited document), ``published_second_hand`` (printed in a document
that itself cites a page not read), ``bem`` (read from the boundary-element
file), ``deployed_configuration`` (the vehicle's own configuration file),
``derived`` (computed from other entries by a stated equation), ``convention``
(a modelling choice, not a measurement), ``estimate`` (an assumption of the
cited source), ``placeholder`` (a value that fills a slot because the model
needs one and nothing measured exists; it is a reported gap, not a claim).

Equations
---------
* Wet section of a circular station of radius R cut by the flat free surface:
  with the waterline at height c above the section centre, c clipped to
  [-R, R], wet area ``A = R^2 (asin(c/R) + pi/2) + c sqrt(R^2 - c^2)`` and
  first moment about the centre ``Q = -(2/3)(R^2 - c^2)^(3/2)`` (the area and
  the first moment of a circular segment; Weisstein, MathWorld, "Circular
  Segment").
* Body of revolution with a cylindrical midbody of length Lc and two
  half-ellipsoid caps of lengths Lf (forward) and La (aft): station radius
  ``R(x) = R`` on the midbody and ``R sqrt(1 - ((|x| - Lc/2)/L_cap)^2)`` on a
  cap; full volume ``V = pi R^2 (Lc + 2 (Lf + La) / 3)``.
* Waterline at station x for the body origin at height z above the calm surface
  and bow-up pitch theta: ``c(x) = (-z - x sin(theta)) / cos(theta)``.
  Submerged volume and centre ``V_s = int A dx``, ``x_B = int x A dx / V_s``,
  ``q_B = int Q dx / V_s`` (q up, body axes).
* Static floating equilibrium: ``rho g V_s = m g`` and the lines of action of
  weight and buoyancy on one vertical, ``x_B cos(theta) - q_B sin(theta) =
  x_G cos(theta) - q_G sin(theta)`` (Archimedes' principle and the moment
  balance of a free floating body; Fossen 2011, Ch. 4, eq. 4.5, p. 60 gives
  the same balance as the restoring vector g(eta) = 0).
* Steady surge: ``T(n) = 1/2 rho Cd pi (D/2)^2 U^2 + (M11 / T1) exp(-3 U) U``
  (axial drag of the hull at zero angle of attack, ``C_D = C_D0`` with
  ``S C_D0 = Cd pi (D/2)^2``; time-constant surge damping with the speed fade;
  Beard and McLain 2012 as cited by MSS ``coeffLiftDrag.m`` 54-68;
  ``forceLiftDrag.m`` 28-40; ``Dmtrx.m`` 30, 48-49; ``remus100.m`` 127, 218).
* Rigid body and added mass moved from the file's reference point (the centre
  of gravity) to the body origin: ``M_CO = H(r)^T M_CG H(r)`` with r the vector
  from the origin to the centre of gravity (Fossen 2011, eq. 3.24, p. 49).
* Geometric similitude at scale s (the second geometry): lengths s, areas
  s^2, volumes and masses s^3 (same density), inertias s^5; a 6x6 mass or
  added-mass matrix ``D M D`` with ``D = diag(s^1.5 I3, s^2.5 I3)`` (the
  added-mass integral rho int phi_i n_j dS of a similar body scales with s^3,
  s^4, s^5 for translation, coupling, rotation); times, speeds and shaft speeds
  at equal Froude number: s^0.5, s^0.5, s^-0.5. The floating equilibrium of
  the scaled body has the same pitch and s times the draft for offsets scaled
  by s: the generator asserts it (a check of its own quadrature).

References
----------
Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics and Motion
    Control. Wiley. Eq. 3.24, p. 49; Ch. 4, eq. 4.5, p. 60.
Fossen, T. I. MSS, MIT, ``HYDRO/vessels_capytaine/LAUV_marie/LAUV_marie.mat`` and
    ``README.md`` (hydrodynamic example of a 2.15 m, 34 kg body with 0.15 m
    diameter; the file's own statement: "a hydrodynamic example rather than a
    validated digital twin"), and MSS-Capytaine ``vessels_capytaine/LAUV_marie/
    {README.md, config.json}`` (cylindrical midbody 1.45 m, equal 0.35 m
    ellipsoidal tapers, radii of gyration of a uniform cylinder, centre of
    gravity 0.01 m below the centre of buoyancy); ``LIBRARY/modeling/
    coeffLiftDrag.m`` 51-53.
Silva, J. E., Terra, B., Martins, R., Borges de Sousa, J. (2007). Modeling and
    Simulation of the LAUV Autonomous Underwater Vehicle. MMAR 2007, pp. 3-4,
    eqs. 5-8 (fin face area 64 cm^2, fin lift coefficient C_LF = 3 in
    ``F = rho C_LF S u^2 delta`` for a pair of fins).
de Moraes, Faltinsen, Esperanca, Sphaier, Lugni (2025). Free-surface
    interaction of a fully appendaged AUV. Ocean Engineering 324:120670, p. 9
    (drag coefficient 0.63 on the bare-hull frontal area, 0.38 on the 0.03 m^2
    frontal area with all appendages, 3 hull diameters deep, Fn 0.42).
NTNU AUR-Lab, "LAUV Marie" (https://www.ntnu.edu/web/aur-lab/lauv-marie, read
    2026-10-10), "Technical specifications": length 215 cm, weight in air
    34 kg, speed 0.5-2.5 m/s.
Marie's deployed DUNE configuration (session 2026-05-21 092659_pf-sss,
    ``Config.ini``): lines 76-77 (fin limits), 235 (maximum motor speed),
    240-242 (speed table), 1779 (propeller position), 1785-1788 (fin positions).
Beard, R. W., McLain, T. W. (2012). Small Unmanned Aircraft. Princeton
    University Press (cited by MSS, not read here).

Author:    Enio Krizman
Date:      2026-10-10
"""

import hashlib
import json
import math
from pathlib import Path

import numpy as np
import scipy.io
from scipy import optimize

from more_transformations.matrix_transforms import MatrixTransforms

HERE = Path(__file__).resolve().parent
MAT_FILE = HERE.parents[2] / "vessel_database" / "data" / "LAUV_marie.mat"
DEG = math.pi / 180.0

# -- primitives: typed once, each with the place it comes from -----------------------------------------
# (value, kind, source, how it becomes a measurement)
PAGE = "NTNU AUR-Lab, LAUV Marie page (https://www.ntnu.edu/web/aur-lab/lauv-marie), read 2026-10-10"
PRIMITIVES = {
    "length": (2.15, "published", PAGE + ", technical specifications: length 215 cm",
               "tape the hull as deployed (quay sheet)"),
    "body_mass": (34.0, "published", PAGE + ", technical specifications: weight in air 34 kg",
                  "weigh the vehicle in air as deployed (quay sheet)"),
    "max_speed": (2.5, "published", PAGE + ", technical specifications: speed 0.5-2.5 m/s (the upper bound)",
                  "logged speed over ground at the maximum shaft speed in still water"),
    "diameter": (0.15, "published_second_hand",
                 "OceanScan LAUV-family hull diameter as printed in the MSS-Capytaine example README",
                 "calipers on the hull"),
    "cylinder_length": (1.45, "estimate",
                        "MSS-Capytaine example README: generic body, 1.45 m cylindrical midbody",
                        "tape the hull and photograph its profile"),
    "forward_cap_length": (0.35, "estimate", "MSS-Capytaine example README: 0.35 m ellipsoidal taper",
                           "tape the nose profile"),
    "aft_cap_length": (0.35, "estimate", "MSS-Capytaine example README: 0.35 m ellipsoidal taper",
                       "tape the tail profile"),
    "centre_of_gravity_below_buoyancy": (0.01, "estimate",
                                         "MSS-Capytaine example README: centre of gravity 0.01 m below the "
                                         "centre of buoyancy (read from the file as main.CG)",
                                         "in-water roll and pitch period; the manufacturer's mass sheet"),
    "hull_frontal_drag_coefficient": (0.63, "published",
                                      "de Moraes et al. (2025), p. 9: 0.63 on the bare-hull frontal area "
                                      "(0.38 on 0.03 m^2 with all appendages), a different, 1.8 m, "
                                      "fully appendaged vehicle of the same family, 3 hull diameters deep",
                                      "steady speed-versus-thrust pairs with a measured thrust law"),
    "oswald_efficiency": (0.3, "published",
                          "MSS coeffLiftDrag.m 51-53: 0.3-0.4 for cylinder-shaped AUVs (generic)",
                          "body-lift slope from logged turns and depth changes"),
    "planform_fraction": (0.7, "estimate",
                          "MSS rule S = 70 % of the rectangle length x diameter (remus100.m 133), "
                          "applied to this hull's own length and diameter",
                          "body-lift slope from logged turns and depth changes"),
    "time_constants": (None, "placeholder",
                       "MSS-Capytaine example: stored inputs of a damping recipe (surge 50 s, sway 5 s, "
                       "yaw 5 s), read from the file as powerBased.T_1236; not a measurement",
                       "speed steps (surge) and turns (sway, yaw) from the vehicle's logs"),
    "damping_ratios": (None, "placeholder",
                       "MSS-Capytaine example: stored inputs (roll 0.2, pitch 0.2), read from the file as "
                       "powerBased.delta_zeta_45; not a measurement",
                       "roll and pitch decay from the vehicle's logs"),
    "kinematic_viscosity": (1e-6, "published", "Fossen 2011, p. 125, below eq. 6.85 (20 degC)",
                            "water temperature from the vehicle's CTD"),
    "fin_area": (0.0064, "published",
                 "Silva et al. (2007), p. 4: fin face area S_fin = 64 cm^2 of the 2007 prototype (1.08 m, "
                 "18 kg), not Marie's fins",
                 "calipers on Marie's fins"),
    "fin_lift_slope": (3.0, "placeholder",
                       "Silva et al. (2007), p. 4: C_LF = 3, read as the per-fin slope of 1/2 rho A C u^2 "
                       "delta for the pair formula F = rho C_LF S u^2 delta (two fins: rho C_LF S = "
                       "2 * 1/2 rho C S); the 2007 prototype's, not Marie's",
                       "turn and pitch rate versus fin angle from the vehicle's logs"),
    "fin_time_constant": (0.1, "placeholder",
                          "no measured servo time constant exists for the vehicle; 0.1 s is a published value "
                          "for an AUV fin servo (Sarhadi 2026, arXiv:2601.01302v2, Fig. 4, p. 4)",
                          "logged step: commanded versus measured fin position"),
    "fin_max_deflection": (25.0 * DEG, "deployed_configuration",
                           "Config.ini line 76: Maximum Fin Rotation = 25.0 (degrees)",
                           "maximum measured fin position in the logs"),
    "fin_max_rate": (300.0 * DEG, "deployed_configuration",
                     "Config.ini line 77: Maximum Rotation Rate = 300.0 (degrees per second)",
                     "fin position rate in the logs"),
    "fin_x": (-0.455, "deployed_configuration",
              "Config.ini lines 1785-1788 (simulator fin positions); the origin of that file is the "
              "simulator's, taken here as this vehicle's origin (an assumption)",
              "tape from the body origin to each fin"),
    "fin_y_z_offset": (0.08, "deployed_configuration",
                       "Config.ini lines 1785-1788 (simulator fin positions)", "tape from the body axis"),
    "propeller_position": ([-0.555, 0.0, 0.0], "deployed_configuration",
                           "Config.ini line 1779 (simulator engine position); origin as for the fins",
                           "tape from the body origin to the propeller"),
    "max_shaft_speed": (2340.0, "deployed_configuration", "Config.ini line 235: Maximum Motor RPMs = 2340",
                        "maximum measured motor speed in the logs"),
    "fit_rpm": (1160.0, "deployed_configuration", "Config.ini lines 240-242: speed table, row 50 %",
                "steady speed-versus-rpm pairs in the logs"),
}
SPEED_TABLE = {  # Config.ini lines 240-242: actuation %, rpm, speed m/s
    "actuation_percent": [0, 25, 33, 50, 66, 75],
    "rpm": [0, 580, 760, 1160, 1535, 1755],
    "speed": [0.0, 0.65, 0.80, 1.40, 1.90, 2.20],
}
FIT_SPEED = 1.40  # the table's speed at fit_rpm
#: the equivalent spheroid of the same length and diameter (Lamb's added mass), printed for this hull by an
#: independent computation (diagonal at the origin, kg and kg m^2, 2026-10-06), and its ratio to the
#: boundary-element diagonal (2026-10-09); index: surge, sway, heave, pitch, yaw
LAMB_PRINTED = {"diagonal": {0: 0.30355, 1: 25.369, 2: 25.369, 4: 5.6225, 5: 5.6225},
                "ratio_to_bem": {0: 0.466, 1: 0.748, 2: 0.749, 4: 0.576, 5: 0.576}}
#: the datasets this file writes: (name, scale, file stem, what it is)
DATASETS = (
    ("Marie (estimated)", 1.0, "marie",
     "the research AUV Marie's one named part: values named by where they come from (estimated, until the quay "
     "weighing and identification replace them)"),
    ("Scaled torpedo (generality check)", 1.6, "scaled_torpedo",
     "not a vehicle: Marie's body scaled geometrically by 1.6, every value derived by similitude, so the same "
     "test file runs on a second torpedo geometry with only the data swapped"),
)
#: the power of the scale each scaled primitive takes (geometric similitude at equal Froude number)
POWERS = {"length": 1, "diameter": 1, "cylinder_length": 1, "forward_cap_length": 1, "aft_cap_length": 1,
          "body_mass": 3, "centre_of_gravity_below_buoyancy": 1, "time_constants": 0.5, "fin_area": 2,
          "fin_x": 1, "fin_y_z_offset": 1, "propeller_position": 1, "max_shaft_speed": -0.5, "fit_rpm": -0.5,
          "max_speed": 0.5}

# -- the second implementation's printed values (independent station quadrature, 768 Gauss-Legendre
# nodes, blind to this file), frozen as the anchors this file's own quadrature must reproduce: the
# reserve, the weight it implies (N), the origin height z above the calm surface (m, up positive),
# the submerged volume (m^3) and the centre of buoyancy height q_B in body axes (m, up positive).
ANCHORS = {
    "full_volume": 0.033870295797,
    "full_buoyancy": 340.574291808,
    "reserve_at_34_kg": 0.0210898,
    "rows": [
        # reserve, weight N, z m, volume m^3, q_B m
        [0.005, 338.879892, -0.071482332, 0.033701786862, -0.000364369],
        [0.010, 337.202269, -0.069490755, 0.033534946333, -0.000716663],
        [0.020, 333.896365, -0.066407442, 0.033206172350, -0.001395917],
        [0.030, 330.654652, -0.063887094, 0.032883782327, -0.002048017],
        [0.050, 324.356468, -0.059699385, 0.032257424568, -0.003286508],
    ],
    "pitch_cells_deg": [  # reserve, x_G - x_B (m), pitch (deg): converged rows of the second implementation
        [0.01, -0.010, 7.4], [0.02, -0.010, 0.8], [0.02, -0.020, 11.3], [0.03, -0.020, 2.1],
        [0.05, -0.020, 0.9], [0.05, -0.005, 0.2],
    ],
}


# -- the file -------------------------------------------------------------------------------------------
def load_vessel():
    return scipy.io.loadmat(str(MAT_FILE), struct_as_record=False, squeeze_me=True)["vessel"]


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def h_matrix(r):
    return np.asarray(MatrixTransforms.H_matrix(np.asarray(r, float)), float)


def move_to_origin(matrix_at_cg, r_origin_to_cg):
    """``H(r)^T M H(r)`` (Fossen 2011, eq. 3.24, p. 49)."""
    h = h_matrix(r_origin_to_cg)
    return h.T @ np.asarray(matrix_at_cg, float) @ h


_LEGENDRE = {}


def _legendre(nodes):
    if nodes not in _LEGENDRE:
        _LEGENDRE[nodes] = np.polynomial.legendre.leggauss(nodes)
    return _LEGENDRE[nodes]


# -- body of revolution, station quadrature -------------------------------------------------------------
class Hull:
    def __init__(self, diameter, cylinder_length, forward_cap, aft_cap):
        self.R, self.Lc, self.Lf, self.La = diameter / 2.0, cylinder_length, forward_cap, aft_cap
        self.length = cylinder_length + forward_cap + aft_cap
        # origin at the midpoint of the cylindrical midbody (the caps are equal for this hull)
        self.x0 = -(cylinder_length / 2.0 + aft_cap)
        self.x1 = cylinder_length / 2.0 + forward_cap

    def full_volume(self):
        return math.pi * self.R**2 * (self.Lc + 2.0 * (self.Lf + self.La) / 3.0)

    def radius(self, x):
        h = self.Lc / 2.0
        if x > h:
            return self.R * math.sqrt(max(0.0, 1.0 - ((x - h) / self.Lf) ** 2))
        if x < -h:
            return self.R * math.sqrt(max(0.0, 1.0 - ((-h - x) / self.La) ** 2))
        return self.R

    def _stations(self, nodes):
        """Gauss-Legendre nodes and weights on the three segments (cap, midbody, cap)."""
        h = self.Lc / 2.0
        t, w = _legendre(nodes)
        xs, ws = [], []
        for a, b in ((self.x0, -h), (-h, h), (h, self.x1)):
            xs.append(0.5 * (b - a) * t + 0.5 * (b + a))
            ws.append(0.5 * (b - a) * w)
        return np.concatenate(xs), np.concatenate(ws)

    def radius_array(self, x):
        h = self.Lc / 2.0
        r = np.full_like(x, self.R)
        fwd, aft = x > h, x < -h
        r[fwd] = self.R * np.sqrt(np.maximum(0.0, 1.0 - ((x[fwd] - h) / self.Lf) ** 2))
        r[aft] = self.R * np.sqrt(np.maximum(0.0, 1.0 - ((-h - x[aft]) / self.La) ** 2))
        return r

    def submerged(self, z_up, theta, nodes=1500):
        """Submerged volume, x_B and q_B by Gauss-Legendre stations on the three segments (the integrand has
        square-root ends where the waterline touches a section; the node count is converged, see
        ``quadrature_error``)."""
        x, w = self._stations(nodes)
        r = self.radius_array(x)
        c = np.clip((-z_up - x * math.sin(theta)) / math.cos(theta), -r, r)
        safe = np.where(r > 0.0, r, 1.0)
        area = np.where(r > 0.0, r * r * (np.arcsin(c / safe) + math.pi / 2.0)
                        + c * np.sqrt(np.maximum(r * r - c * c, 0.0)), 0.0)
        moment = -(2.0 / 3.0) * np.maximum(r * r - c * c, 0.0) ** 1.5
        volume = float(np.sum(w * area))
        if volume <= 0.0:
            return 0.0, 0.0, 0.0
        return volume, float(np.sum(w * x * area)) / volume, float(np.sum(w * moment)) / volume

    def quadrature_error(self, z_up, theta):
        """Largest relative difference of volume, x_B and q_B between 1500 and 6000 nodes per segment."""
        a, b = self.submerged(z_up, theta, 1500), self.submerged(z_up, theta, 6000)
        return max(abs(u - v) / max(abs(v), 1e-3) for u, v in zip(a, b))


def equilibrium(hull, rho, mass, g, cg_x, cg_q):
    """``(z_up, theta, volume, x_B, q_B)`` of the floating body; ``cg_q`` is up positive in body axes. The root
    is the first bow-up solution with the horizontal lines of action of weight and buoyancy equal (the level
    solution when the centre of gravity is under the centre of buoyancy)."""
    target = mass / rho

    def draft_for(theta):
        f = lambda z: hull.submerged(z, theta)[0] - target
        return optimize.brentq(f, -hull.length, hull.length, xtol=1e-15, rtol=1e-14, maxiter=200)

    def residual(theta):
        z = draft_for(theta)
        volume, x_b, q_b = hull.submerged(z, theta)
        return (x_b * math.cos(theta) - q_b * math.sin(theta)) - (cg_x * math.cos(theta) - cg_q * math.sin(theta))

    thetas = np.linspace(0.0, 1.3, 131)
    values = [residual(t) for t in thetas]
    if abs(values[0]) < 1e-14:
        root = 0.0
    else:
        root = None
        for a, b, fa, fb in zip(thetas[:-1], thetas[1:], values[:-1], values[1:]):
            if fa * fb < 0.0:
                root = optimize.brentq(residual, a, b, xtol=1e-15, rtol=1e-14, maxiter=200)
                break
        if root is None:
            raise RuntimeError(f"no bow-up equilibrium for cg_x={cg_x}")
    z = draft_for(root)
    volume, x_b, q_b = hull.submerged(z, root)
    return z, root, volume, x_b, q_b


# -- steady surge ---------------------------------------------------------------------------------------
def surge_drag(speed, rho, cd, diameter, m11, t1):
    return (0.5 * rho * cd * math.pi * (diameter / 2.0) ** 2 * speed**2
            + (m11 / t1) * math.exp(-3.0 * speed) * speed)


def steady_speed(rpm, k, rho, cd, diameter, m11, t1):
    thrust = k * (rpm / 60.0) ** 2
    return optimize.brentq(lambda u: surge_drag(u, rho, cd, diameter, m11, t1) - thrust, 0.0, 10.0, xtol=1e-14)


# -- assemble --------------------------------------------------------------------------------------------
def entry(value, kind, source, route=""):
    out = {"value": value, "kind": kind, "source": source}
    if route:
        out["route_to_measured"] = route
    return out


def similitude_matrix(scale):
    """``D = diag(s^1.5 I3, s^2.5 I3)``: ``D M D`` scales a 6x6 mass matrix by s^3, s^4, s^5 (translation,
    coupling, rotation)."""
    return np.diag([scale**1.5] * 3 + [scale**2.5] * 3)


def build(scale=1.0):
    vessel = load_vessel()
    main = vessel.main
    rho, g = float(main.rho), float(main.g)
    # the file's numbers are the published ones (the page and the file agree), then scaled
    assert float(main.m) == PRIMITIVES["body_mass"][0] and abs(float(main.Lpp) - PRIMITIVES["length"][0]) < 1e-12
    file_cg = np.asarray(main.CG, float) + 0.0  # origin -> centre of gravity in body axes (down positive)
    assert abs(file_cg[2] - PRIMITIVES["centre_of_gravity_below_buoyancy"][0]) < 1e-12
    d = similitude_matrix(scale)
    mrb_cg = d @ np.asarray(vessel.MRB, float) @ d
    a_cg = d @ np.asarray(vessel.A, float)[:, :, -1] @ d  # the 10 rad/s sentinel entry is the infinite-frequency limit
    assert float(np.asarray(vessel.freqs)[-1]) == 10.0
    r_cg = scale * file_cg
    added_mass_co = move_to_origin(a_cg, r_cg)
    rigid_co = move_to_origin(mrb_cg, r_cg)
    unscaled = move_to_origin(np.asarray(vessel.A, float)[:, :, -1], file_cg)
    assert np.abs(added_mass_co - d @ unscaled @ d).max() <= 1e-12 * max(1.0, np.abs(added_mass_co).max())
    power = vessel.powerBased
    t1236 = np.asarray(power.T_1236, float)
    damping_ratios = np.asarray(power.delta_zeta_45, float).tolist()

    p = {k: v[0] for k, v in PRIMITIVES.items()}
    p["time_constants"] = [float(t1236[0]), float(t1236[1]), float(t1236[3])]  # surge, sway, yaw (T3 = T2)
    for name, k in POWERS.items():
        value = p[name]
        p[name] = [v * scale**k for v in value] if isinstance(value, list) else value * scale**k
    mass = p["body_mass"]
    assert abs(float(rigid_co[0, 0]) - mass) <= 1e-12 * mass
    assert abs(rigid_co[0, 4] - mass * r_cg[2]) < 1e-12 * max(1.0, mass) and abs(rigid_co[1, 3] + mass * r_cg[2]) < 1e-12 * max(1.0, mass)
    inertia = rigid_co[3:6, 3:6]
    inertia_diagonal = np.diag(inertia).tolist()
    inertia_products = [-inertia[0, 1] + 0.0, -inertia[1, 2] + 0.0, -inertia[0, 2] + 0.0]  # I = [[Ix,-Ixy,-Ixz],[-Ixy,Iy,-Iyz],...]
    time_constants = p["time_constants"]

    def prim(name, value=None):
        """The entry of a primitive: as typed at scale 1; at another scale, the scaled value labelled derived."""
        value = p[name] if value is None else value
        kind, source, route = PRIMITIVES[name][1:] if name in PRIMITIVES else ("placeholder", "", "")
        if scale != 1.0 and POWERS.get(name, 0) != 0:
            return entry(value, "derived", f"geometric similitude, scale {scale} to the power {POWERS[name]}, of "
                                           f"Marie's value ({kind}: {source})")
        return entry(value, kind, source, route)

    def scaled(e):
        """A computed entry: at another scale its kind is derived from Marie's, by similitude."""
        if scale == 1.0:
            return e
        return entry(e["value"], "derived", f"geometric similitude, scale {scale}, of Marie's value ({e['kind']}: "
                                            f"{e['source']})")

    diameter = p["diameter"]
    length = p["cylinder_length"] + p["forward_cap_length"] + p["aft_cap_length"]
    assert abs(length - p["length"]) < 1e-12, "the profile's length is the published length"
    hull = Hull(diameter, p["cylinder_length"], p["forward_cap_length"], p["aft_cap_length"])
    full_volume = hull.full_volume()
    full_buoyancy = rho * g * full_volume
    weight = mass * g
    reserve = (full_buoyancy - weight) / weight
    planform = p["planform_fraction"] * length * diameter
    mean_density = mass / full_volume
    fit_speed = FIT_SPEED * scale**0.5
    speed_table_rpm = [r * scale**-0.5 for r in SPEED_TABLE["rpm"]]
    speed_table_speed = [v * scale**0.5 for v in SPEED_TABLE["speed"]]

    # the thrust coefficient: one fit to one row of the configured speed table (placeholder)
    m11 = float(rigid_co[0, 0] + added_mass_co[0, 0])
    drag_fit = surge_drag(fit_speed, rho, p["hull_frontal_drag_coefficient"], diameter, m11, time_constants[0])
    thrust_coefficient = drag_fit / (p["fit_rpm"] / 60.0) ** 2  # N per (rev/s)^2

    def fin(label, position, lift_axis, note):
        return {"label": label, "fin_position": scaled(entry(position, "deployed_configuration", note)),
                "lift_axis": entry(lift_axis, "convention",
                                   "unit direction of positive lift in body axes (forward-right-down); the "
                                   "simulator's fin force directions, Config.ini lines 1781-1784")}

    fx, off = p["fin_x"], p["fin_y_z_offset"]
    o = f"{off:.3g}"
    fins = [
        fin(f"vertical fin, z = +{o} m (below the axis, FRD)", [fx, 0.0, off], [0.0, 1.0, 0.0], "Config.ini line 1785 (fin 0)"),
        fin(f"horizontal fin, y = -{o} m (port)", [fx, -off, 0.0], [0.0, 0.0, -1.0], "Config.ini line 1786 (fin 1)"),
        fin(f"horizontal fin, y = +{o} m (starboard)", [fx, off, 0.0], [0.0, 0.0, -1.0], "Config.ini line 1787 (fin 2)"),
        fin(f"vertical fin, z = -{o} m (above the axis, FRD)", [fx, 0.0, -off], [0.0, 1.0, 0.0], "Config.ini line 1788 (fin 3)"),
    ]
    fin_common = {
        "chord_axis": entry([1.0, 0.0, 0.0], "convention", "unit chord direction along the hull axis"),
        "fin_area": prim("fin_area"),
        "lift_slope": prim("fin_lift_slope"),
        "zero_lift_drag": entry(0.0, "placeholder", "no value for this fin; read only by the linear section form"),
        "servo_dynamics": entry("first_order_lag", "convention", "servo with a lag, a rate limit and an angle limit"),
        "servo_rate_limit": entry(True, "convention", "servo with a lag, a rate limit and an angle limit"),
        "servo_angle_limit": entry("on_command", "convention", "servo with a lag, a rate limit and an angle limit"),
        "time_constant": prim("fin_time_constant"),
        "max_deflection": prim("fin_max_deflection"),
        "max_rate": prim("fin_max_rate"),
        "inflow": entry("rigid_point", "convention", "velocity at the fin's own point of the rigid body"),
        "flow_angle": entry("small_angle", "convention", "small-angle flow angle"),
        "interference": entry("none", "convention", "no body interference on the fin"),
        "section": entry("quadratic_drag", "convention", "lift slope with drag proportional to the squared angle"),
    }

    parameters = {
        "vehicle": {
            "length": prim("length", length),
            "diameter": prim("diameter"),
            "body_mass": prim("body_mass"),
            "body_center_of_gravity": scaled(entry(r_cg.tolist(), *PRIMITIVES["centre_of_gravity_below_buoyancy"][1:])),
            "rigid_body_form": entry("given_inertia", "convention",
                                     "inertia of a uniform cylinder (the file's radii of gyration), the same "
                                     "body the added mass was computed for"),
            "inertia_diagonal": scaled(entry(inertia_diagonal, "derived",
                                             "file MRB (about the centre of gravity) moved to the origin by "
                                             "H(r)^T MRB H(r), Fossen 2011 eq. 3.24 p. 49; an estimate (a uniform "
                                             "cylinder's radii of gyration)",
                                             "the manufacturer's mass properties; a pendulum test")),
            "inertia_products": scaled(entry(inertia_products, "derived", "as inertia_diagonal")),
            "body_density_method": entry("computed", "convention", "unread with rigid_body_form 'given_inertia'"),
            "body_density": entry(mean_density, "derived",
                                  "mean density of the idealised hull, mass / full volume; unread with "
                                  "rigid_body_form 'given_inertia'"),
            "roll_added_inertia_ratio": entry(0.0, "convention", "unread with added_mass_form 'given'"),
            "added_mass_form": entry("given", "convention", "boundary-element added mass of the idealised hull"),
            "added_mass_matrix": scaled(entry(added_mass_co.tolist(), "bem",
                                              "LAUV_marie.mat A(infinity) (frequency entry 36, the 10 rad/s sentinel; "
                                              "zero speed, 5 m submergence), shifted from the centre of gravity to the "
                                              "body origin by H(r)^T A H(r), r = main.CG, Fossen 2011 eq. 3.24 p. 49",
                                              "a boundary-element run of the real hull with fins, sonar head and "
                                              "antenna")),
            "coriolis_form": entry("kirchhoff_full", "convention", "every added-mass Coriolis term"),
            "current_form": entry("none", "convention", "still water"),
            "integration_max_step": entry(0.05, "convention", "numerical step of the vehicle's own integrator"),
        },
        "hydrostatics": {
            "buoyancy_method": entry("floating", "convention",
                                     "partly submerged hull form, integrated over the hull profile: exact at the "
                                     "surface and fully submerged (a regime is a state, not a separate part)"),
            "center_of_buoyancy": entry([0.0, 0.0, 0.0], "convention",
                                        "the body origin is the centre of buoyancy of the idealised hull "
                                        "(the file's CB differs by 2 micrometres)"),
            "displaced_volume": entry(full_volume, "derived",
                                      "pi R^2 (Lc + 2 (Lf + La)/3) of the idealised hull profile; read by "
                                      "'from_volume' only (the floating method integrates the profile)",
                                      "the hull lines; a weighing in water (reserve buoyancy)"),
            "buoyancy": entry(full_buoyancy, "derived", "rho g V of the full hull; read by 'given' only"),
            "cylinder_length": prim("cylinder_length"),
            "forward_cap_length": prim("forward_cap_length"),
            "aft_cap_length": prim("aft_cap_length"),
        },
        "hydrodynamics": {
            "time_constants": prim("time_constants"),
            "damping_ratios": entry(damping_ratios, *PRIMITIVES["damping_ratios"][1:]),
            "oswald_efficiency": prim("oswald_efficiency"),
            "planform_area_method": entry("computed", "convention", "S = f L D"),
            "planform_fraction": prim("planform_fraction"),
            "planform_area": entry(planform, "derived", "f L D; unread with planform_area_method 'computed'"),
            "parasitic_drag_coefficient_method": entry("computed", "convention", "C_D0 = Cd pi (D/2)^2 / S"),
            "hull_frontal_drag_coefficient": prim("hull_frontal_drag_coefficient"),
            "parasitic_drag_coefficient": entry(p["hull_frontal_drag_coefficient"] * math.pi * (diameter / 2) ** 2
                                                / planform, "derived",
                                                "Cd pi (D/2)^2 / S; unread with 'computed'"),
            "smooth_speed": entry(True, "convention", "numerical guard"),
            "smooth_speed_epsilon": entry(1e-9, "convention", "numerical guard"),
        },
        "fin_common": fin_common,
        "fins": fins,
        "propulsor": {
            "selectors": {
                "drive": entry("speed_command_rpm", "convention", "shaft speed in revolutions per minute"),
                "steering": entry("fixed", "convention", "fixed propeller"),
                "mounting": entry("fixed_euler", "convention", "shaft along the hull axis"),
                "inflow": entry("none", "convention", "the lumped law reads no advance speed"),
                "law": entry("quadratic_asymmetric", "placeholder",
                             "thrust = k w |w|, w = 2 pi n: a lumped law, because no propeller geometry of the "
                             "vehicle exists in any document"),
                "housing": entry("none", "convention", "no housing force"),
            },
            "min_command": entry(0.0, "placeholder",
                                 "no document says whether the vehicle runs its motor astern; none modelled",
                                 "the logged shaft speed range"),
            "max_command": prim("max_shaft_speed"),
            "position": prim("propeller_position"),
            "orientation": entry([0.0, 0.0, 0.0], "convention", "shaft along the hull axis"),
            "positive_thrust_coefficient": entry(
                thrust_coefficient / (2.0 * math.pi) ** 2, "placeholder",
                f"one fit: thrust equal to the modelled hull drag at {fit_speed:.4g} m/s and {p['fit_rpm']:.0f} rpm, "
                "the speed table's 50 % row (a configured table, not a measurement), written for w in rad/s",
                "thrust against shaft speed from a bollard pull and the logs"),
            "negative_thrust_coefficient": entry(
                thrust_coefficient / (2.0 * math.pi) ** 2, "placeholder",
                "equal to the forward coefficient; unread with min_command 0", "bollard pull astern"),
            "shaft_torque_scale": entry(0.0, "placeholder",
                                        "no value for the share of shaft torque that rolls the vehicle; not modelled",
                                        "steady roll offset against shaft speed in straight runs"),
        },
        "site": {
            "site_form": entry("given", "convention", "gravity and density as in the hydrodynamic file"),
            "latitude": entry(0.0, "placeholder", "unread with site_form 'given'"),
            "gravity": entry(g, "bem", "LAUV_marie.mat main.g"),
            "water_density": entry(rho, "bem", "LAUV_marie.mat main.rho"),
            "kinematic_viscosity": prim("kinematic_viscosity"),
        },
    }

    def model_speed(rpm):
        return steady_speed(rpm, thrust_coefficient, rho, p["hull_frontal_drag_coefficient"], diameter, m11,
                            time_constants[0])

    lamb_diagonal = {k: v * scale**(3 if k < 3 else 5) for k, v in LAMB_PRINTED["diagonal"].items()}
    derived = {
        "mat_sha256": sha256(MAT_FILE),
        "scale": scale,
        "rho": rho, "g": g, "mass": mass, "weight": weight,
        "full_volume": full_volume, "full_buoyancy": full_buoyancy, "reserve": reserve,
        "mesh_volume": float(main.nabla) * scale**3,
        "mesh_reserve": (rho * g * float(main.nabla) * scale**3 - weight) / weight,
        "rigid_body_mass_matrix_at_origin": rigid_co.tolist(),
        "added_mass_at_origin": added_mass_co.tolist(),
        "added_mass_at_centre_of_gravity": a_cg.tolist(),
        "lamb_spheroid": {
            "diagonal": {str(k): v for k, v in lamb_diagonal.items()},
            "ratio_to_bem": {str(k): v for k, v in LAMB_PRINTED["ratio_to_bem"].items()},
            "source": "the equivalent spheroid of the same length and diameter (Lamb's added mass), printed for "
                      "Marie's hull by an independent computation (2026-10-06), and its ratio to the "
                      "boundary-element diagonal (2026-10-09); at another scale, by similitude (s^3, s^5; the "
                      "ratios do not change)",
        },
        "m11": m11,
        "thrust_coefficient": thrust_coefficient,
        "speed_table": {
            "rpm": speed_table_rpm[1:], "table_speed": speed_table_speed[1:],
            "model_speed": [model_speed(r) for r in speed_table_rpm[1:]],
            "max_rpm": p["max_shaft_speed"], "max_model_speed": model_speed(p["max_shaft_speed"]),
            "published_max_speed": p["max_speed"],
            "published_max_speed_source": prim("max_speed")["source"],
        },
    }
    return hull, parameters, derived, {"rho": rho, "g": g, "mass": mass, "cg": r_cg}


NOMINAL_OFFSETS = (0.0, -0.005, -0.010, -0.020)  # m, x_G - x_B at scale 1


def floating_cases(mass, g, full_buoyancy, scale):
    """The grid of the floating-equilibrium reference: the declared part, then reserve x longitudinal offset
    cells (the reserve is set through the mass: m = B_full / (g (1 + reserve)); the offset x_G - x_B is the
    centre of gravity's forward coordinate with the centre of buoyancy at the origin, scaled with the body)."""
    cases = [("declared", mass, 0.0)]
    cases += [(f"declared_dx_{dx * scale * 1000:+.0f}mm", mass, dx * scale) for dx in (-0.010, -0.020)]
    for reserve in (0.01, 0.02, 0.03, 0.04):
        for dx in NOMINAL_OFFSETS:
            if reserve <= 0.011 and dx <= -0.0199:
                continue  # a 45 degree trim: outside the range a vehicle floats in
            cases.append((f"reserve_{reserve:.2f}_dx_{dx * scale * 1000:+.0f}mm",
                          full_buoyancy / (g * (1.0 + reserve)), dx * scale))
    return cases


def write_dataset(name, scale, stem, description):
    hull, parameters, derived, fixed = build(scale)
    rho, g, mass, cg = fixed["rho"], fixed["g"], fixed["mass"], fixed["cg"]
    full_buoyancy = derived["full_buoyancy"]
    q_g = -float(cg[2])  # body axes, up positive
    rows, solved = [], []
    for label, m, dx in floating_cases(mass, g, full_buoyancy, scale):
        z_up, theta, volume, x_b, q_b = equilibrium(hull, rho, m, g, dx, q_g)
        solved.append((z_up, theta))
        reserve = (full_buoyancy - m * g) / (m * g)
        rows.append([label, f"{m:.17g}", f"{dx:.17g}", f"{reserve:.17g}", f"{-z_up:.17g}", f"{theta:.17g}",
                     f"{volume:.17g}", f"{x_b:.17g}", f"{-q_b:.17g}"])
    header = ["case", "mass_kg", "cg_x_m", "reserve", "origin_depth_m", "pitch_rad", "submerged_volume_m3",
              "buoyancy_x_m", "buoyancy_z_down_m"]
    lines = [",".join(header)] + [",".join(r) for r in rows]
    (HERE / f"{stem}_floating_reference.csv").write_text("\n".join(lines) + "\n")

    if scale == 1.0:
        anchors = dict(ANCHORS)
        anchors["quadrature_check"] = []
        for reserve, weight, z, volume, q_b in ANCHORS["rows"]:
            m = full_buoyancy / (g * (1.0 + reserve))
            z_up, theta, vol, x_b, qb = equilibrium(hull, rho, m, g, 0.0, q_g)
            anchors["quadrature_check"].append({"reserve": reserve, "weight_error": abs(m * g - weight),
                                                "z_error": abs(z_up - z), "volume_error": abs(vol - volume),
                                                "q_b_error": abs(qb - q_b)})
        derived["anchor_check"] = anchors["quadrature_check"]
        derived["full_volume_anchor_error"] = abs(derived["full_volume"] - ANCHORS["full_volume"])
        derived["full_buoyancy_anchor_error"] = abs(full_buoyancy - ANCHORS["full_buoyancy"])
        derived["pitch_cells"] = []
        for reserve, dx, deg in ANCHORS["pitch_cells_deg"]:
            m = full_buoyancy / (g * (1.0 + reserve))
            _, theta, *_ = equilibrium(hull, rho, m, g, dx, q_g)
            derived["pitch_cells"].append({"reserve": reserve, "dx": dx, "pitch_deg": theta / DEG, "printed_deg": deg})

    out = {
        "description": f"Parameter set of {description}. Values are grouped by the component that declares them; "
                       "every entry names its kind and source; see SOURCE.md.",
        "part": name,
        "parameters": parameters,
        "derived": derived,
    }
    return out, solved


def main():
    results = {}
    for name, scale, stem, description in DATASETS:
        results[stem] = write_dataset(name, scale, stem, description)
    # similitude check of the quadrature: the scaled body floats at the same pitch, s times deeper
    (marie, marie_solved), (other, other_solved) = results["marie"], results["scaled_torpedo"]
    scale = other["derived"]["scale"]
    pitch_error = max(abs(a[1] - b[1]) for a, b in zip(marie_solved, other_solved))
    draft_error = max(abs(b[0] - scale * a[0]) for a, b in zip(marie_solved, other_solved))
    assert pitch_error < 1e-9 and draft_error < 1e-9, (pitch_error, draft_error)
    other["derived"]["similitude_check"] = {"against": "marie_floating_reference.csv", "pitch_error_rad": pitch_error,
                                            "draft_error_m": draft_error}
    for stem, (out, _) in results.items():
        (HERE / f"{stem}_parameters.json").write_text(json.dumps(out, indent=1, sort_keys=False) + "\n")
        d = out["derived"]
        print(stem, "full volume", d["full_volume"], "reserve", d["reserve"])
        if "anchor_check" in d:
            print("anchor errors", max(max(c["z_error"], c["volume_error"], c["q_b_error"]) for c in d["anchor_check"]))
            for cell in d["pitch_cells"]:
                print(cell)
        print("speed table", d["speed_table"])
    print("similitude check", other["derived"]["similitude_check"])


if __name__ == "__main__":
    main()
