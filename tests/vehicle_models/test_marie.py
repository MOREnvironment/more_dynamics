"""Gates for the research AUV Marie as one named part of the table (boundary-element) hull, "Marie (estimated)":
``SubmergedBemAuv`` with ``SubmergedRestoring`` (``buoyancy_method = "immersed_body_of_revolution"``), ``BemHullLoads``
and four ``Fin``. A variant is named by where its data come from (later "Marie (measured)", "Marie (identified)"),
never by a physics regime: the floating restoring is exact at the surface and fully submerged, so the surface phase and
the dive are states of one part. The numbers are ``data/marie/marie_parameters.json`` (each value with its kind).

Every gate runs on a second geometry too, with only the data swapped (``data/marie/scaled_torpedo_parameters.json``:
the body scaled by 1.6 with Froude similitude, written by ``generate_marie_references.py``; a check of the generality
of the vehicle type, not a vehicle).

The references are independent second implementations (numpy and scipy, written from the printed equations), closed
forms and identities, never the code under test: the power-based reduction with scipy's shape-preserving interpolant,
the viscous damping from its defining decay time, the surge force in closed form, eq. 4.5 for the submerged restoring,
the floating equilibria by adaptive quadrature (the generator), the fin lift from the low-aspect-ratio slope, an RK4 of
the graph's own dynamics and the Froude similitude of the two sets.

Tolerances: ``G1`` equals a second implementation of the same equation (1e-9 relative to the scale of the quantity);
``REDUCTION`` (1e-12 relative) the reduction against its second implementation (the two sum in another order);
``QUADRATURE_*`` allow the station sums of the floating form their numerical error against the converged quadrature of
the generator (the tolerances of the earlier form of this file, unchanged).

References
----------
Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics and Motion Control. Wiley. Eq. 3.24, p. 49; eq. 3.34,
    p. 50; eq. 3.44, p. 52; eq. 4.5, p. 60; eqs. 6.43, 6.52, pp. 120-121; eqs. 12.17-12.21, p. 351; eqs. 12.57-12.59,
    p. 365; eq. 12.226, p. 400.
Fossen, T. I. (2025). Maneuvering coefficient estimation from frequency-dependent added mass and damping: a
    power-based approach. Ocean Engineering 341:122494 (not read; the method as stated in the library module).
DNV GL (2017). DNVGL-RP-C205, sec. 3.5.5, pp. 64-65 (JONSWAP spectrum).
Prestero, T. (2001). Verification of a six-degree of freedom simulation model for the REMUS autonomous underwater
    vehicle. MIT/WHOI MSc thesis. Eqs. 4.38-4.39, p. 31 (fin lift slope).
Weisstein, E. W. Circular Segment. MathWorld.

Author:    Enio Krizman
Date:      2026-10-10
"""

import json
import math
import re
from pathlib import Path

import numpy as np
import pytest

import vehicle_contract as vc

DATA_DIR = vc.DATA / "marie"
LIBRARY = Path(__file__).resolve().parents[2]

DATASETS = {"Marie (estimated)": "marie", "Scaled torpedo (generality check)": "scaled_torpedo"}
MARIE = "Marie (estimated)"
SCALED = "Scaled torpedo (generality check)"
SCALE = 1.6

G1 = 1e-9
REDUCTION = 1e-12
QUADRATURE_DEPTH = 5e-5        # m, floating form against the converged second implementation (per metre of scale)
QUADRATURE_FORCE = 1e-2        # N on the heave row at the reference equilibrium (per unit of scale^3)
QUADRATURE_MOMENT = 1e-2       # N m on the pitch row at the reference equilibrium (per unit of scale^4)
KINDS = {"published", "fitted to field data", "estimate", "design choice"}
CURRENT_INPUTS = 3             # the current enters first (current_form full_attitude)


def _floating_rows(path):
    lines = path.read_text().splitlines()
    names = lines[0].split(",")
    return [{"case": cells[0], **{n: float(c) for n, c in zip(names[1:], cells[1:])}}
            for cells in (line.split(",") for line in lines[1:])]


class Dataset:
    def __init__(self, name, stem):
        self.name = name
        self.parameters = json.loads((DATA_DIR / f"{stem}_parameters.json").read_text())
        self.rows = _floating_rows(DATA_DIR / f"{stem}_floating_reference.csv")
        self.scale = SCALE if stem == "scaled_torpedo" else 1.0

    def values(self, group):
        return {name: entry["value"] for name, entry in self.parameters["parameters"][group].items()}

    def fin_values(self, index):
        fin = self.parameters["parameters"]["fins"][index]
        return {**self.values("fin_common"), "fin_position": fin["fin_position"]["value"],
                "lift_axis": fin["lift_axis"]["value"]}

    @property
    def weight(self):
        return self.values("vehicle")["body_mass"] * self.values("site")["gravity"]

    @property
    def full_buoyancy(self):
        site = self.values("site")
        return site["water_density"] * site["gravity"] * self.values("hydrostatics")["displaced_volume"]


DATA = {name: Dataset(name, stem) for name, stem in DATASETS.items()}
PARTS = tuple(DATA)
IDS = {MARIE: "marie", SCALED: "scaled"}
each_part = pytest.mark.parametrize("part", PARTS, ids=[IDS[p] for p in PARTS])

_CACHE = {}


def vehicle(part=MARIE, *, overrides=None, hydrostatics=None, hydrodynamics=None, fin=None):
    """The initialised vehicle of a dataset; the keyword dictionaries are test-side parameter overrides."""
    key = json.dumps([part, overrides, hydrostatics, hydrodynamics, fin], sort_keys=True)
    if key in _CACHE:
        return _CACHE[key]
    data, p = DATA[part], vc.plugins()
    actuators = [(p.Fin, {**data.fin_values(k), **(fin or {})}) for k in range(4)]
    vessel = vc.build_context(p.SubmergedBemAuv, {**data.values("vehicle"), **(overrides or {})},
                              site=(p.Site, data.values("site")),
                              hydrostatics=(p.SubmergedRestoring, {**data.values("hydrostatics"), **(hydrostatics or {})}),
                              hydrodynamics=(p.BemHullLoads, {**data.values("hydrodynamics"), **(hydrodynamics or {})}),
                              actuators=actuators, reference=False)
    _CACHE[key] = vessel
    return vessel


def state(eta=None, nu=None):
    return np.r_[np.zeros(6) if eta is None else eta, np.zeros(6) if nu is None else nu]


def inputs(vessel, fins=None):
    u = np.zeros(vc.input_size(vessel))
    if fins is not None:
        u[CURRENT_INPUTS:CURRENT_INPUTS + len(fins)] = fins
    return u


def signal(vessel, name, x, u=None):
    return vessel.signals(x, inputs(vessel) if u is None else u)[name]


def restoring_force(vessel, eta):
    return signal(vessel, "hydrostatics.restoring_force", state(eta))


def submerged_attitudes(data, seed, count=6):
    """Random poses deep enough that the whole hull is below the surface (depth >= one hull length)."""
    length = data.values("vehicle")["length"]
    rng = np.random.default_rng(seed)
    return [np.r_[rng.uniform(-1, 1, 2), rng.uniform(length, 3.0 * length), rng.uniform(-0.4, 0.4, 2),
                  rng.uniform(-1, 1)] for _ in range(count)]


def h_matrix(r):
    """H(r) = [[I, S(r)^T], [0, I]] (Fossen 2011, eq. 3.24, p. 49), second implementation."""
    s = np.array([[0, -r[2], r[1]], [r[2], 0, -r[0]], [-r[1], r[0], 0]])
    return np.block([[np.eye(3), s.T], [np.zeros((3, 3)), np.eye(3)]])


def _trapezoid(y, x, axis=0):
    y = np.moveaxis(np.asarray(y, float), axis, 0)
    return np.tensordot(np.diff(x), 0.5 * (y[1:] + y[:-1]), axes=(0, 0))


def reduction(data, table, omega_p=None, gamma=None):
    """The power-based reduction (Fossen 2025) with scipy's shape-preserving cubic Hermite interpolant (Fritsch and
    Carlson 1980) and numpy's trapezoid rule, the JONSWAP shape of DNVGL-RP-C205 2017, p. 64: a second implementation."""
    from scipy.interpolate import PchipInterpolator
    v = data.values("vehicle")
    w = np.asarray(v["frequencies"], float)
    values = np.asarray(v[table], float).reshape(w.size, 36)
    keep = w < 10.0
    w, values = w[keep], values[keep]
    omega_p = v["omega_p"] if omega_p is None else omega_p
    gamma = v["spectrum_gamma"] if gamma is None else gamma
    fine = np.linspace(1e-6 if w[0] == 0.0 else w[0], w[-1], 1000)
    sigma = np.where(fine <= omega_p, 0.07, 0.09)
    spectrum = fine ** -5 * np.exp(-1.25 * (omega_p / fine) ** 4) * gamma ** np.exp(-0.5 * ((fine - omega_p) / (sigma * omega_p)) ** 2)
    weight = spectrum / _trapezoid(spectrum, fine)
    return _trapezoid(PchipInterpolator(w, values)(fine) * weight[:, None], fine).reshape(6, 6)


def at_origin(data, matrix, symmetric=False):
    h = h_matrix(np.asarray(data.values("vehicle")["database_reference"], float))
    moved = h.T @ matrix @ h
    return 0.5 * (moved + moved.T) if symmetric else moved


def rigid_body(data):
    """``[[m I, -m S(r_g)], [m S(r_g), I_o]]`` with the tensor about the CO (Fossen 2011, eq. 3.44, p. 52)."""
    v = data.values("vehicle")
    m, r = v["body_mass"], np.asarray(v["body_center_of_gravity"], float)
    (ix, iy, iz), (ixy, iyz, ixz) = v["inertia_diagonal"], v["inertia_products"]
    s = np.array([[0, -r[2], r[1]], [r[2], 0, -r[0]], [-r[1], r[0], 0]])
    inertia = np.array([[ix, -ixy, -ixz], [-ixy, iy, -iyz], [-ixz, -iyz, iz]])
    return np.block([[m * np.eye(3), -m * s], [m * s, inertia]])


# --------------------------------------------------------------------------------------------------------
# 1. Builds, the data files, generality
# --------------------------------------------------------------------------------------------------------
@each_part
def test_the_named_part_builds_with_four_fins(part):
    vessel = vehicle(part)
    assert type(vessel).__name__ == "SubmergedBemAuv"
    assert vc.input_size(vessel) == CURRENT_INPUTS + 4            # the current, four fin commands
    assert vc.graph_of(vessel).step.size1_in(0) == 12             # the fins follow their command: no servo state


@each_part
def test_the_parameter_file_is_one_set_and_every_entry_has_a_kind(part):
    data = DATA[part]
    assert data.parameters["part"] == part and "parts" not in data.parameters
    assert data.values("hydrostatics")["buoyancy_method"] == "immersed_body_of_revolution"

    def walk(node, where):
        if isinstance(node, dict):
            if "value" in node:
                assert node["kind"] in KINDS, (where, node)
                assert set(node) == {"value", "kind"}, (where, sorted(node))          # a value and its kind, nothing else
                return 1
            return sum(walk(v, f"{where}.{k}") for k, v in node.items() if k != "command")
        if isinstance(node, list):
            return sum(walk(v, f"{where}[{i}]") for i, v in enumerate(node))
        return 0

    assert walk(data.parameters["parameters"], "parameters") > 50


def test_the_published_values_of_the_named_part():
    """Length, mass and section height are published values of the vehicle; the tables are those of an idealised hull
    of that length and diameter (an estimate), until an identification from the logs replaces them."""
    group = DATA[MARIE].parameters["parameters"]["vehicle"]
    assert {name for name, e in group.items() if e["kind"] == "published"} == {"length", "body_mass", "height"}
    for name in ("added_mass", "potential_damping", "frequencies"):
        assert group[name]["kind"] == "estimate"


def test_the_data_folder_carries_no_record_of_where_a_value_came_from():
    """Values by kind only: no address, date, configuration file, line number or source text in the data."""
    pattern = re.compile(r"https?://|\b20\d\d-\d\d-\d\d\b|\.ini\b|\bline \d|\blines \d", re.IGNORECASE)
    for name in ("marie_parameters.json", "scaled_torpedo_parameters.json", "SOURCE.md"):
        text = (DATA_DIR / name).read_text()
        assert not pattern.search(text), (name, pattern.search(text).group(0))


@each_part
def test_the_part_sets_every_value_its_plugins_read_unset(part):
    """The part gives every unset parameter of the four plugin types itself (a named part holds its numbers)."""
    from more_transformations.more_casadi_transformations import is_unset
    data, p = DATA[part], vc.plugins()
    unread = {"SubmergedRestoring": {"buoyancy"}, "Fin": {"lift_slope", "time_constant", "max_rate"}}
    for group, cls in (("vehicle", p.SubmergedBemAuv), ("hydrostatics", p.SubmergedRestoring),
                       ("hydrodynamics", p.BemHullLoads), ("fin_common", p.Fin)):
        given = set(data.values(group)) | ({"fin_position", "lift_axis"} if group == "fin_common" else set())
        missing = {d.name for d in cls.PARAMETERS if is_unset(d.default_value)} - given - unread.get(cls.__name__, set())
        assert not missing, (cls.__name__, sorted(missing))


def test_no_library_module_names_the_vehicle():
    """The plugins, the models and the plugin types know physics, not vehicles: ``grep -ri marie`` over the three
    folders finds nothing; the vehicle family's name appears only where a line cites a file or a paper."""
    hits = []
    family = re.compile(r"\blauv\w*\b", re.IGNORECASE)
    cite = re.compile(r"[\w./-]+\.(m|mat|py|json|yaml|csv|pdf)\b")
    for folder in ("models", "plugins", "plugin_types"):
        for path in (LIBRARY / "more_dynamics" / folder).rglob("*"):
            if not path.is_file() or "__pycache__" in path.parts:
                continue
            text = path.read_bytes().decode("utf-8", errors="ignore")
            hits += [f"{path.name}:{n}: {line}" for n, line in enumerate(text.splitlines(), 1) if "marie" in line.lower()]
            if path.suffix == ".py" and not any(s in str(path) for s in ("hull_vessel", "linear_surface", "jet_nozzle")):
                hits += [f"{path.name}:{n}: {line}" for n, line in enumerate(text.splitlines(), 1)
                         if family.search(line) and not cite.search(line) and "Silva" not in line and "Moraes" not in line]
    assert not hits, hits


def test_the_second_geometry_is_the_froude_similitude_of_the_first():
    """Generality: the scaled set shares no number with the part, yet at Froude-similar states (lengths s, speeds
    s^0.5, rates s^-0.5) its accelerations are the part's, linear ones equal and angular ones times s^-1. The strips
    use Hoerner's section here (the cylinder's Reynolds number is not Froude-similar)."""
    hoerner = {"cross_flow_model": "hoerner"}
    small, large = vehicle(MARIE, hydrodynamics=hoerner), vehicle(SCALED, hydrodynamics=hoerner)
    rng = np.random.default_rng(3)
    s = SCALE
    for _ in range(6):
        eta = np.r_[rng.uniform(-1, 1, 2), rng.uniform(-0.1, 3.0), rng.uniform(-0.5, 0.5, 2), rng.uniform(-3, 3)]
        nu = np.r_[rng.uniform(0.5, 2.0), rng.uniform(-0.2, 0.2, 2), rng.uniform(-0.3, 0.3, 3)]
        fins = rng.uniform(-0.3, 0.3, 4)
        a = np.asarray(vc.graph_of(small).step(state(eta, nu), inputs(small, fins))).ravel()[6:]
        b = np.asarray(vc.graph_of(large).step(state(np.r_[eta[:3] * s, eta[3:]], np.r_[nu[:3] * s ** 0.5, nu[3:] * s ** -0.5]),
                                               inputs(large, fins))).ravel()[6:]
        expected = np.r_[a[:3], a[3:] / s]
        assert np.abs(b - expected).max() <= 1e-9 * max(1.0, np.abs(expected).max())


# --------------------------------------------------------------------------------------------------------
# 2. Reduced tables, rigid body, damping
# --------------------------------------------------------------------------------------------------------
@each_part
def test_the_added_mass_is_the_spectrum_averaged_table_moved_to_the_origin(part):
    data = DATA[part]
    vessel = vehicle(part)
    resolved = vc.matrix(signal(vessel, "added_mass_matrix", state()))
    expected = at_origin(data, reduction(data, "added_mass"), symmetric=True)
    assert np.abs(resolved - expected).max() <= REDUCTION * np.abs(expected).max()


@each_part
def test_the_potential_damping_is_the_spectrum_averaged_table_moved_to_the_origin(part):
    data = DATA[part]
    resolved = vc.matrix(signal(vehicle(part), "potential_damping_matrix", state()))
    expected = at_origin(data, reduction(data, "potential_damping"))
    assert np.abs(resolved - expected).max() <= REDUCTION * np.abs(expected).max()


def test_the_reduction_of_a_table_constant_in_frequency_is_that_constant_and_ignores_the_sentinel():
    """Closed forms: the normalised spectrum integrates to one, so a table constant in frequency reduces to itself;
    the entry at 10 rad/s (infinite frequency of the data format) does not enter."""
    from more_dynamics.models.frequency_reduction.power_based import power_based_reduction
    w = DATA[MARIE].values("vehicle")["frequencies"]
    rng = np.random.default_rng(1)
    constant = rng.normal(size=36)
    table = np.tile(constant, len(w))
    table[-36:] = 1e6                                                       # the sentinel row
    out = power_based_reduction(w)(added_mass=table, potential_damping=table, omega_p=2.0, spectrum_gamma=3.3)
    assert np.abs(np.asarray(out["added_mass_at_table_reference"]) - constant.reshape(6, 6)).max() <= 1e-12


@pytest.mark.parametrize("omega_p,gamma", [(1.0, 1.0), (2.7, 3.3), (6.0, 7.0)])
def test_the_reduction_at_other_spectra_equals_the_second_implementation(omega_p, gamma):
    from more_dynamics.models.frequency_reduction.power_based import power_based_reduction
    data = DATA[MARIE]
    v = data.values("vehicle")
    out = power_based_reduction(v["frequencies"])(added_mass=v["added_mass"], potential_damping=v["potential_damping"],
                                                  omega_p=omega_p, spectrum_gamma=gamma)
    for name in ("added_mass", "potential_damping"):
        expected = reduction(data, name, omega_p, gamma)
        assert np.abs(np.asarray(out[f"{name}_at_table_reference"]) - expected).max() <= REDUCTION * np.abs(expected).max()


def test_the_reduction_preserves_the_shape_of_a_monotone_table_and_has_finite_derivatives():
    """Shape preservation (Fritsch and Carlson 1980): a table monotone in frequency reduces to a value between its
    smallest and largest entry; the Jacobian with respect to the tables is finite at a flat table (every if_else
    branch guarded)."""
    import casadi as ca
    from more_dynamics.models.frequency_reduction.power_based import power_based_reduction
    w = DATA[MARIE].values("vehicle")["frequencies"]
    block = power_based_reduction(w)
    steps = np.repeat(np.minimum(np.arange(len(w)), 10), 36).astype(float)       # rising, then flat
    out = np.asarray(block(added_mass=steps, potential_damping=steps, omega_p=2.0, spectrum_gamma=3.3)["added_mass_at_table_reference"])
    assert np.all(out >= 0.0) and np.all(out <= 10.0)
    x = ca.SX.sym("x", 36 * len(w))
    jacobian = ca.Function("j", [x], [ca.jacobian(block(x, x, 2.0, 3.3)[0], x)])
    for table in (np.zeros(36 * len(w)), steps):
        assert np.all(np.isfinite(np.asarray(jacobian(table))))


@each_part
def test_the_rigid_body_is_the_full_tensor_about_the_origin(part):
    data = DATA[part]
    resolved = vc.matrix(signal(vehicle(part), "rigid_body_mass_matrix", state()))
    assert np.abs(resolved - rigid_body(data)).max() <= G1 * np.abs(rigid_body(data)).max()
    v = data.values("vehicle")
    r = np.asarray(v["body_center_of_gravity"], float)
    s = np.array([[0, -r[2], r[1]], [r[2], 0, -r[0]], [-r[1], r[0], 0]])
    inertia_g = rigid_body(data)[3:, 3:] + v["body_mass"] * s @ s                  # I_g = I_o + m S^2 (eq. 3.34, p. 50)
    assert np.all(np.linalg.eigvalsh(inertia_g) > 0.0)


@each_part
def test_the_total_damping_decays_each_free_mode_in_its_time_constant_and_adds_the_damping_ratios(part):
    """``M_ii / D_ii = T_i`` for surge, sway, heave and yaw and ``Bv_ii = 2 delta_zeta sqrt(M_ii G_ii)`` for roll and
    pitch, about the table reference point (the mass-damper-spring form, Fossen 2011, eqs. 12.57-12.59, p. 365)."""
    data = DATA[part]
    v, d = data.values("vehicle"), data.values("hydrodynamics")
    mass = rigid_body(data)
    h_inverse = h_matrix(-np.asarray(v["database_reference"], float))
    m_r = h_inverse.T @ mass @ h_inverse + reduction(data, "added_mass")
    b_r = reduction(data, "potential_damping")
    g = np.zeros((6, 6))
    g[2:5, 2:5] = np.asarray(v["hydrostatic_stiffness_345"], float).reshape(3, 3)
    viscous = vc.matrix(signal(vehicle(part), "hydrodynamics.viscous_damping", state(eta=[0, 0, 20.0, 0, 0, 0])))
    for t, i in zip(d["T_1236"], (0, 1, 2, 5)):
        assert m_r[i, i] / (b_r[i, i] + viscous[i, i]) == pytest.approx(t, rel=1e-10)
    for dz, i in zip(d["delta_zeta_45"], (3, 4)):
        assert viscous[i, i] == pytest.approx(2.0 * dz * math.sqrt(m_r[i, i] * g[i, i]), rel=1e-10)


@each_part
@pytest.mark.parametrize("speed", [0.5, 1.5, 2.5])
def test_the_surge_force_at_pure_surge_is_the_damping_and_the_parasitic_drag_in_closed_form(part, speed):
    """``X = -(M_11 / T_1) exp(-k_u U) U - 1/2 rho c_d pi (H/2)^2 U^2``: the surge damping of its time constant with its
    speed fade and the hull's parasitic drag on its frontal area (no lift at zero angle, no cross-flow)."""
    data = DATA[part]
    v, d, site = data.values("vehicle"), data.values("hydrodynamics"), data.values("site")
    m11 = v["body_mass"] + reduction(data, "added_mass")[0, 0]
    expected = (-(m11 / d["T_1236"][0]) * math.exp(-d["surge_damping_blend"] * speed) * speed
                - 0.5 * site["water_density"] * d["c_d"] * math.pi * (v["height"] / 2) ** 2 * speed ** 2)
    force = signal(vehicle(part), "hydrodynamics.hydrodynamic_force", state(eta=[0, 0, 20.0, 0, 0, 0], nu=[speed, 0, 0, 0, 0, 0]))
    assert force[0] == pytest.approx(expected, rel=1e-12)


def _lagrangian_coriolis(a, nu):
    """``C_A(nu) = [[0, -S(A11 v1 + A12 v2)], [-S(A11 v1 + A12 v2), -S(A21 v1 + A22 v2)]]`` (Fossen 2011, eq. 6.43,
    p. 120), second implementation."""
    def skew(v):
        return np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    s1 = skew(a[:3, :3] @ nu[:3] + a[:3, 3:] @ nu[3:])
    s2 = skew(a[3:, :3] @ nu[:3] + a[3:, 3:] @ nu[3:])
    return np.block([[np.zeros((3, 3)), -s1], [-s1, -s2]])


@each_part
def test_the_part_selects_the_munk_shortcut_and_the_full_form_is_the_lagrangian_matrix(part):
    """The part's form is ``munk_couplings_removed`` until a straight-line stability check of the hull with its fins
    decides; the full form (the type's default) is every term of eq. 6.43, the Munk moment of eq. 6.52 (p. 121)
    included, and the shortcut the same matrix with the pitch-heave, pitch-surge, yaw-surge and yaw-sway entries and
    their transposes zero."""
    data = DATA[part]
    assert data.values("vehicle")["coriolis_form"] == "munk_couplings_removed"
    nu = np.r_[1.5, 0.1, -0.05, 0.02, 0.03, -0.04]
    x = state(eta=[0, 0, 20.0, 0, 0, 0], nu=nu)
    full = vehicle(part, overrides={"coriolis_form": "kirchhoff_full"})
    a = vc.matrix(signal(full, "added_mass_matrix", x))
    expected = _lagrangian_coriolis(a, nu)
    assert np.abs(vc.matrix(signal(full, "added_mass_coriolis_matrix", x)) - expected).max() <= G1 * np.abs(expected).max()
    for i, j in ((4, 2), (4, 0), (5, 0), (5, 1)):
        expected[i, j] = expected[j, i] = 0.0
    shortcut = vc.matrix(signal(vehicle(part), "added_mass_coriolis_matrix", x))
    assert np.abs(shortcut - expected).max() <= G1 * np.abs(expected).max()


# --------------------------------------------------------------------------------------------------------
# 3. Statics fully submerged
# --------------------------------------------------------------------------------------------------------
def _g_submerged(eta, weight, buoyancy, r_g, r_b):
    """Restoring vector of a submerged body (Fossen 2011, eq. 4.5, p. 60), second implementation."""
    phi, theta = eta[3], eta[4]
    s_t, c_t, s_p, c_p = math.sin(theta), math.cos(theta), math.sin(phi), math.cos(phi)
    x_g, y_g, z_g = r_g
    x_b, y_b, z_b = r_b
    return np.array([
        (weight - buoyancy) * s_t,
        -(weight - buoyancy) * c_t * s_p,
        -(weight - buoyancy) * c_t * c_p,
        -(y_g * weight - y_b * buoyancy) * c_t * c_p + (z_g * weight - z_b * buoyancy) * c_t * s_p,
        (z_g * weight - z_b * buoyancy) * s_t + (x_g * weight - x_b * buoyancy) * c_t * c_p,
        -(x_g * weight - x_b * buoyancy) * c_t * s_p - (y_g * weight - y_b * buoyancy) * s_t])


@each_part
@pytest.mark.parametrize("x_g", [0.0, -0.005, 0.01], ids=["cg_on_cb", "cg_aft", "cg_forward"])
def test_submerged_restoring_equals_the_second_implementation_at_attitudes(part, x_g):
    """Fully below the surface: net heave force B - W (positively buoyant, never neutral) and the moments of the offset
    of the centres, at several attitudes. ``x_g`` is the offset at the part's scale 1."""
    data = DATA[part]
    weight, buoyancy = data.weight, data.full_buoyancy
    r_b = data.values("hydrostatics")["center_of_buoyancy"]
    r_g = [r_b[0] + x_g * data.scale, 0.0, data.values("vehicle")["body_center_of_gravity"][2]]
    vessel = vehicle(part, overrides={"body_center_of_gravity": r_g})
    for eta in submerged_attitudes(data, 20261010):
        expected = -_g_submerged(eta, weight, buoyancy, r_g, r_b)
        assert np.abs(restoring_force(vessel, eta) - expected).max() <= G1 * weight
    level = np.array([0.0, 0.0, 2.0 * data.values("vehicle")["length"], 0.0, 0.0, 0.0])
    at_rest = restoring_force(vessel, level)
    assert at_rest[2] == pytest.approx(weight - buoyancy, abs=1e-9 * weight)    # upward: negative in FRD
    assert at_rest[2] < 0.0
    if x_g < 0.0:
        assert at_rest[4] > 0.0                                                 # the centre of gravity aft: bow up


@each_part
def test_released_at_depth_the_part_rises(part):
    vessel = vehicle(part)
    x = state(eta=[0, 0, 10.0, 0, 0, 0])
    for _ in range(5):
        x = vessel.step_vector(x, inputs(vessel), 1.0)
    assert x[2] < 10.0 - 1e-3                      # depth (NED z) decreased
    assert x[8] < 0.0                              # heave velocity upward


# --------------------------------------------------------------------------------------------------------
# 4. Floating equilibrium
# --------------------------------------------------------------------------------------------------------
ROW_CASES = [(part, row) for part in PARTS for row in DATA[part].rows]


def _floating(part, row):
    data = DATA[part]
    r_g = [row["cg_x_m"], 0.0, data.values("vehicle")["body_center_of_gravity"][2]]
    return vehicle(part, overrides={"body_mass": row["mass_kg"], "body_center_of_gravity": r_g})


@pytest.mark.parametrize("part,row", ROW_CASES, ids=[f"{IDS[p]}-{r['case']}" for p, r in ROW_CASES])
def test_floating_equilibrium_against_the_second_implementation(part, row):
    """At the second implementation's draft and pitch the part's restoring force and moment vanish: the station sums
    of the vehicle's form agree with the converged quadrature to the stated tolerances. ``cg_x_m`` is the longitudinal
    position of the centre of gravity; the mass sets the reserve buoyancy."""
    data = DATA[part]
    vessel = _floating(part, row)
    eta = np.array([0.0, 0.0, row["origin_depth_m"], 0.0, row["pitch_rad"], 0.0])
    force = restoring_force(vessel, eta)
    assert abs(force[2]) <= QUADRATURE_FORCE * data.scale ** 3, force
    assert abs(force[4]) <= QUADRATURE_MOMENT * data.scale ** 4, force
    assert abs(force[0]) <= QUADRATURE_FORCE * data.scale ** 3
    step = 1e-2 * data.scale
    above, below = (restoring_force(vessel, eta + np.array([0, 0, h, 0, 0, 0]))[2] for h in (-step, step))
    assert above > force[2] > below                # heave: deeper lowers the net downward force (stiffness > 0)
    if abs(row["pitch_rad"]) < 0.2:      # near level the waterline cuts the parallel body: a waterplane stiffness
        assert (above - below) / (2 * step) > 100.0 * data.scale ** 2     # N/m, not a flat curve


@each_part
def test_floating_depth_and_pitch_by_solving_the_vehicle_for_the_declared_part_and_the_offsets(part):
    """The same cells by solving the vehicle for equilibrium from a perturbed guess: depth and bow-up pitch."""
    from scipy import optimize
    data = DATA[part]
    for row in data.rows:
        if not row["case"].startswith(("declared", "reserve_0.02", "reserve_0.03", "reserve_0.04")):
            continue
        vessel = _floating(part, row)

        def residual(v):
            f = restoring_force(vessel, np.array([0.0, 0.0, v[0], 0.0, v[1], 0.0]))
            return [f[2], f[4]]

        guess = [row["origin_depth_m"] + 5e-3 * data.scale, row["pitch_rad"] + 2e-2]
        solution = optimize.fsolve(residual, guess, xtol=1e-12)
        assert solution[0] == pytest.approx(row["origin_depth_m"],
                                            abs=(QUADRATURE_DEPTH + 1e-3 * abs(row["pitch_rad"])) * data.scale)
        assert solution[1] == pytest.approx(row["pitch_rad"], abs=3e-4 + 1e-3 * abs(row["pitch_rad"]))
        if row["cg_x_m"] < data.values("hydrostatics")["center_of_buoyancy"][0]:
            assert solution[1] > 0.0               # the centre of gravity aft of the centre of buoyancy: nose up


@each_part
def test_the_floating_form_equals_the_volume_form_when_every_station_is_submerged(part):
    """One part for the dive and the surface: fully below the surface the floating form is the submerged form with the
    same displaced volume, to 1e-12 of the weight."""
    data = DATA[part]
    floating = vehicle(part)
    from_volume = vehicle(part, hydrostatics={"buoyancy_method": "from_volume"})
    for eta in submerged_attitudes(data, 7, count=8):
        difference = restoring_force(floating, eta) - restoring_force(from_volume, eta)
        assert np.abs(difference).max() <= 1e-12 * data.weight * max(1.0, data.scale)


@each_part
def test_the_floating_restoring_has_finite_derivatives_through_the_surface(part):
    """Every section dry, cut or whole: the Jacobian of the vehicle's dynamics in the state is finite (the guarded
    circle segments)."""
    import casadi as ca
    vessel = vehicle(part)
    dynamics = vc.graph_of(vessel).step
    x = ca.SX.sym("x", 12)
    jacobian = ca.Function("j", [x], [ca.jacobian(dynamics(x, inputs(vessel)), x)])
    radius = DATA[part].values("vehicle")["height"] / 2
    for depth in (-2.0 * radius, -radius, 0.0, 0.5 * radius, radius, 3.0 * radius):
        for pitch in (0.0, 0.3):
            assert np.all(np.isfinite(np.asarray(jacobian(state(eta=[0, 0, depth, 0.1, pitch, 0], nu=[1, 0, 0, 0, 0, 0])))))


@each_part
def test_reserve_buoyancy_sweep_reported_against_the_pitch_it_gives(part):
    """Reported, not gated beyond the second implementation's own monotonic result: at an aft offset of the centre of
    gravity the equilibrium pitch falls as the reserve grows from 1 % to 4 % (the declared reserve is an estimate until
    the vehicle is weighed in water)."""
    rows = DATA[part].rows
    offset = next(r["cg_x_m"] for r in rows if r["case"].startswith("declared_dx_"))
    cells = [r for r in rows if r["case"].startswith("reserve_") and abs(r["cg_x_m"] - offset) < 1e-12]
    pitches = [r["pitch_rad"] * 180.0 / math.pi for r in cells]
    assert len(cells) == 4 and pitches == sorted(pitches, reverse=True)
    print(f"{part}: reserve -> bow-up pitch (deg) at {offset * 1000:+.0f} mm:",
          [(round(r["reserve"], 3), round(p, 2)) for r, p in zip(cells, pitches)])


# --------------------------------------------------------------------------------------------------------
# 5. Pop-up: control released at depth
# --------------------------------------------------------------------------------------------------------
@each_part
def test_popup_from_ten_metres_with_no_control_input_surfaces_and_settles_at_the_floating_equilibrium(part):
    data = DATA[part]
    vessel = vehicle(part)
    equilibrium = next(r for r in data.rows if r["case"] == "declared")
    radius = data.values("vehicle")["height"] / 2.0
    x = state(eta=[0, 0, 10.0, 0, 0, 0])
    depths, time_to_surface = [], None
    for second in range(1, 1201):
        x = vessel.step_vector(x, inputs(vessel), 1.0)
        assert np.all(np.isfinite(x))
        depths.append(x[2])
        if time_to_surface is None and x[2] < radius:
            time_to_surface = second
    assert time_to_surface is not None, "the vehicle never broke the surface"
    print(f"{part}: pop-up time to surface (reported, not gated): {time_to_surface} s; settled depth {depths[-1]:.5f} m")
    assert depths[-1] == pytest.approx(equilibrium["origin_depth_m"], abs=1e-3 * data.scale)
    assert x[4] == pytest.approx(equilibrium["pitch_rad"], abs=1e-3)
    assert np.abs(x[6:]).max() < 1e-3                                            # at rest


# --------------------------------------------------------------------------------------------------------
# 6. Fins: layout, lift, limits
# --------------------------------------------------------------------------------------------------------
def _fin_wrench(vessel, k, deflection, speed):
    fins = np.zeros(4)
    fins[k] = deflection
    out = vessel.signals(state(eta=[0, 0, 20.0, 0, 0, 0], nu=[speed, 0, 0, 0, 0, 0]), inputs(vessel, fins))
    return out[f"actuators.{k}.generated_force"]


def _low_aspect_ratio_slope(values):
    """``C_L,alpha = [1/(2 pi a) + 1/(pi AR_e)]^-1``, ``AR_e = k b^2 / S`` (Prestero 2001, eqs. 4.38-4.39, p. 31)."""
    aspect = values["aspect_ratio_factor"] * values["fin_span"] ** 2 / values["fin_area"]
    return 1.0 / (1.0 / (2.0 * math.pi * values["section_efficiency"]) + 1.0 / (math.pi * aspect))


@each_part
@pytest.mark.parametrize("k", range(4))
@pytest.mark.parametrize("speed", [0.8, 1.9])
def test_each_fin_wrench_is_the_lift_of_the_low_aspect_ratio_slope(part, k, speed):
    """No drag law: the wrench is [f; r x f], f = 1/2 rho A U^2 C_L,alpha d along the fin's lift axis (Fossen 2011,
    eq. 12.226, p. 400; pure surge, so the flow angle is zero and the angle of attack the deflection)."""
    data = DATA[part]
    vessel = vehicle(part)
    values = data.fin_values(k)
    rho = data.values("site")["water_density"]
    r, axis = np.asarray(values["fin_position"], float), np.asarray(values["lift_axis"], float)
    for delta in (0.05, 0.2, 0.4):
        lift = 0.5 * rho * values["fin_area"] * speed ** 2 * _low_aspect_ratio_slope(values) * delta * axis
        expected = np.r_[lift, np.cross(r, lift)]
        got = _fin_wrench(vessel, k, delta, speed)
        assert np.abs(got - expected).max() <= 1e-9 * max(1.0, np.abs(expected).max())


@each_part
def test_fin_layout_rudder_pair_gives_sway_and_yaw_stern_pair_heave_and_pitch(part):
    """The two rudder fins lift sideways and turn the hull in yaw, the two stern planes lift vertically and pitch it;
    no fin gives roll (every fin is on the hull axis)."""
    vessel = vehicle(part)
    yaw = _fin_wrench(vessel, 0, 0.2, 1.5) + _fin_wrench(vessel, 1, 0.2, 1.5)
    pitch = _fin_wrench(vessel, 2, 0.2, 1.5) + _fin_wrench(vessel, 3, 0.2, 1.5)
    assert abs(yaw[1]) > 1e-3 and abs(yaw[5]) > 1e-3 and abs(yaw[2]) < 1e-12 and abs(yaw[4]) < 1e-12
    assert abs(pitch[2]) > 1e-3 and abs(pitch[4]) > 1e-3 and abs(pitch[1]) < 1e-12 and abs(pitch[5]) < 1e-12
    assert abs(yaw[3]) < 1e-12 and abs(pitch[3]) < 1e-12


@each_part
def test_fin_deflection_limit_is_enforced_on_the_command(part):
    data = DATA[part]
    vessel = vehicle(part)
    limit = data.values("fin_common")["max_deflection"]
    assert np.allclose(_fin_wrench(vessel, 0, 1.0, 1.5), _fin_wrench(vessel, 0, limit, 1.5), rtol=0, atol=1e-12)
    assert np.abs(_fin_wrench(vessel, 0, 0.5 * limit, 1.5)).max() < np.abs(_fin_wrench(vessel, 0, limit, 1.5)).max()


# --------------------------------------------------------------------------------------------------------
# 7. Straight-line stability (reported: the part's Coriolis form is the owner's choice after this check)
# --------------------------------------------------------------------------------------------------------
def _sway_yaw_roots(vessel, speed, h=1e-6):
    """Roots of the 2 x 2 sway-yaw linearisation at straight level flight (Fossen 2011, eqs. 12.17-12.21, p. 351), by
    central differences of the graph's dynamics."""
    dynamics = vc.graph_of(vessel).step
    a = np.zeros((2, 2))
    for j, index in enumerate((1, 5)):
        dv = np.zeros(6)
        dv[index] = h
        plus = np.asarray(dynamics(state(eta=[0, 0, 20.0, 0, 0, 0], nu=np.r_[speed, 0, 0, 0, 0, 0] + dv), inputs(vessel))).ravel()
        minus = np.asarray(dynamics(state(eta=[0, 0, 20.0, 0, 0, 0], nu=np.r_[speed, 0, 0, 0, 0, 0] - dv), inputs(vessel))).ravel()
        a[:, j] = (plus[[7, 11]] - minus[[7, 11]]) / (2 * h)
    return np.linalg.eigvals(a)


@each_part
def test_straight_line_sway_yaw_roots_with_both_coriolis_forms_are_reported(part):
    """Reported, not gated: the largest real part of the sway-yaw roots at the cruise speed with the fins fixed. The
    Munk moment destabilises (the full form's largest root lies to the right of the shortcut's)."""
    speed = 1.5 * DATA[part].scale ** 0.5
    full = _sway_yaw_roots(vehicle(part, overrides={"coriolis_form": "kirchhoff_full"}), speed)
    shortcut = _sway_yaw_roots(vehicle(part), speed)
    print(f"{part}: sway-yaw roots at {speed:.2f} m/s: full {np.round(full, 4)}, shortcut {np.round(shortcut, 4)}")
    assert np.all(np.isfinite(full)) and np.all(np.isfinite(shortcut))
    assert full.real.max() > shortcut.real.max()


# --------------------------------------------------------------------------------------------------------
# 8. Integration
# --------------------------------------------------------------------------------------------------------
@each_part
def test_the_vehicles_own_step_equals_a_test_side_rk4_of_the_graphs_dynamics(part):
    data = DATA[part]
    vessel = vehicle(part)
    dynamics = vc.graph_of(vessel).step
    step = data.values("vehicle").get("integration_max_step", 0.05)
    rng = np.random.default_rng(11)
    x = state(eta=np.r_[rng.uniform(-1, 1, 2), data.values("vehicle")["length"], rng.uniform(-0.2, 0.2, 2),
                        rng.uniform(-1, 1)], nu=np.r_[1.0, rng.uniform(-0.1, 0.1, 5)])
    u = inputs(vessel, [0.05, -0.05, 0.02, 0.0])
    reference = x.copy()
    for _ in range(10):                             # 0.5 s in sub-steps of the declared integration step
        k1 = np.asarray(dynamics(reference, u)).ravel()
        k2 = np.asarray(dynamics(reference + 0.5 * step * k1, u)).ravel()
        k3 = np.asarray(dynamics(reference + 0.5 * step * k2, u)).ravel()
        k4 = np.asarray(dynamics(reference + step * k3, u)).ravel()
        reference = reference + step / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
    assert np.abs(vessel.step_vector(x, u, 0.5) - reference).max() <= 1e-10
