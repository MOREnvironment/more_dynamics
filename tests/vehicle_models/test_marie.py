"""Gates for the research AUV Marie as one named part of ``SpheroidAuv``,
"Marie (estimated)": a variant named by where its data come from (later
"Marie (measured)" after the quay weighing, "Marie (identified)" from the logs),
never by a physics regime. Its hull form is the floating one, exact at the
surface and when fully submerged, so the surface phase and the dive are states
of the same part. The numbers are ``data/marie/marie_parameters.json`` (every
value with its kind and source, written by ``generate_marie_references.py``).

Every gate runs on a second torpedo geometry too, with only the data swapped
(``data/marie/scaled_torpedo_parameters.json``: Marie's body scaled
geometrically by 1.6, every value derived by similitude; a check of the
generality of the vehicle type, not a vehicle). Nothing in this file is
specific to either set beyond the two file names in ``DATASETS``.

The references are independent second implementations (numpy and scipy, written
from the printed equations) or published tables, never the code under test.

What the vehicle is built from, and what is not built yet (each test that waits
on a build names every build it waits on in its failure message):

* **added mass and inertia given** (``added_mass_form = "given"`` with a 6x6
  ``added_mass_matrix`` at the body origin, as the vessel-database reader writes
  it, and ``rigid_body_form = "given_inertia"`` with ``inertia_diagonal`` and
  ``inertia_products``, the full inertia tensor of the body the added mass was
  computed for): both options of the vehicle type, not built yet.
  ``VEHICLE_SELECTORS`` is the one place to change the two names.
* **floating restoring** (``buoyancy_method = "floating"`` with
  ``cylinder_length``, ``forward_cap_length`` and ``aft_cap_length`` declared
  on ``SubmergedRestoring``; the hull diameter is the vehicle's): not built yet.
  ``HULL_PROFILE_NAMES`` is the one place to change the three names.
* **the propulsor** (one plugin ``Propulsor`` whose six sockets are selectors
  named as the propulsor contract names its parts, with the parts' own
  parameters declared on it): not built yet. ``PROPULSOR_SELECTOR_NAMES`` and
  the parameter names of the parts are the contract's.

A test that does not need a missing build runs today (the fin gates and the
Lamb comparison on a stand-in with the hull forms the plugin has today, the
data files, the rule scans).

Tolerances: ``G1`` equals a second implementation of the same equation (1e-9
relative to the scale of the quantity); ``QUADRATURE`` allows the floating
form's station quadrature its declared numerical error (the second
implementation's own quadrature is converged to better than 1e-6 relative);
every other number a test reports and does not gate says so in its name.

References
----------
Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics and Motion
    Control. Wiley. Eq. 4.5, p. 60 (restoring vector of a submerged body),
    eq. 3.24, p. 49 (the transformation between two reference points),
    eq. 12.226, p. 400 (fin lift).
Fossen, T. I. MSS, MIT, ``HYDRO/vessels_capytaine/LAUV_marie/LAUV_marie.mat``
    (added mass, inertia and the vessel's hydrostatic quantities; Fossen's
    idealised hydrodynamic example, not a measurement of the vehicle).
NTNU AUR-Lab, "LAUV Marie" (https://www.ntnu.edu/web/aur-lab/lauv-marie, read
    2026-10-10): length 215 cm, weight in air 34 kg, speed 0.5-2.5 m/s.
Silva, J. E., Terra, B., Martins, R., Borges de Sousa, J. (2007). Modeling and
    Simulation of the LAUV Autonomous Underwater Vehicle. MMAR 2007, pp. 3-4.
Weisstein, E. W. Circular Segment. MathWorld (area and first moment of a
    circular segment, used by the floating-equilibrium reference).

Author:    Enio Krizman
Date:      2026-10-10
"""

import importlib
import json
import math
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

import vehicle_contract as vc

DATA_DIR = vc.DATA / "marie"
LIBRARY = Path(__file__).resolve().parents[2]
MAT_FILE = LIBRARY / "tests" / "vessel_database" / "data" / "LAUV_marie.mat"
MAT_SHA256 = "7eed3a0c4643ad446b6e45db9f74659052cc56866b3bc3a9c540bd5befcadb61"

#: the named part, and the second geometry of the generality check (file stem of its data)
DATASETS = {"Marie (estimated)": "marie", "Scaled torpedo (generality check)": "scaled_torpedo"}
MARIE = "Marie (estimated)"

G1 = 1e-9
QUADRATURE_DEPTH = 5e-5        # m, floating form against the converged second implementation (per metre of scale)
QUADRATURE_FORCE = 1e-2        # N on the heave row at the reference equilibrium (per unit of scale^3)
QUADRATURE_MOMENT = 1e-2       # N m on the pitch row at the reference equilibrium (per unit of scale^4)
STEP = 0.05                    # s, the vehicle's own integration step (declared parameter, set by the part)
DEG = math.pi / 180.0

# The one place the proposed names live (see the module docstring).
VEHICLE_SELECTORS = ("added_mass_form", "rigid_body_form")
GIVEN_RIGID_BODY_NAMES = ("rigid_body_form", "inertia_diagonal", "inertia_products", "added_mass_matrix")
HULL_PROFILE_NAMES = ("cylinder_length", "forward_cap_length", "aft_cap_length")
PROPULSOR_SELECTOR_NAMES = ("drive", "steering", "mounting", "inflow", "law", "housing")
WAITS = {
    "given": "the given added mass and given inertia options of the vehicle type",
    "floating": "the floating buoyancy method of the submerged restoring",
    "propulsor": "the propulsor plugin",
}


# --------------------------------------------------------------------------------------------------------
# The data files
# --------------------------------------------------------------------------------------------------------
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
        self.derived = self.parameters["derived"]
        self.scale = self.derived["scale"]

    def values(self, group):
        return {name: entry["value"] for name, entry in self.parameters["parameters"][group].items()
                if isinstance(entry, dict) and "value" in entry}

    def fin_values(self, index):
        fin = self.parameters["parameters"]["fins"][index]
        return {**self.values("fin_common"), "fin_position": fin["fin_position"]["value"],
                "lift_axis": fin["lift_axis"]["value"]}

    def propulsor_values(self):
        group = self.parameters["parameters"]["propulsor"]
        selectors = {name: entry["value"] for name, entry in group["selectors"].items()}
        own = {name: entry["value"] for name, entry in group.items() if name != "selectors"}
        return {**selectors, **own}

    def similitude(self):
        """``D`` with ``D M D`` the similitude of a 6x6 mass matrix (s^3, s^4, s^5)."""
        return np.diag([self.scale**1.5] * 3 + [self.scale**2.5] * 3)


DATA = {name: Dataset(name, stem) for name, stem in DATASETS.items()}
PARTS = tuple(DATA)
IDS = {MARIE: "marie", "Scaled torpedo (generality check)": "scaled"}
each_part = pytest.mark.parametrize("part", PARTS, ids=[IDS[p] for p in PARTS])


# --------------------------------------------------------------------------------------------------------
# The builder
# --------------------------------------------------------------------------------------------------------
def missing_builds(p, *, propulsor, method):
    """The builds a vehicle of the part needs and the library does not have yet."""
    missing = []
    if not set(GIVEN_RIGID_BODY_NAMES) <= {d.name for d in p.SpheroidAuv.PARAMETERS}:
        missing.append("given")
    if method == "floating" and not set(HULL_PROFILE_NAMES) <= {d.name for d in p.SubmergedRestoring.PARAMETERS}:
        missing.append("floating")
    if propulsor:
        try:
            importlib.import_module("more_dynamics.plugins.force_producers.propulsor")
        except ImportError:
            missing.append("propulsor")
    return missing


_CACHE = {}


def vehicle(part=MARIE, *, propulsor=False, stand_in=False, overrides=None, hydrostatics=None, fin=None, site=None,
            hydrostatics_method=None):
    """The initialised vehicle of a dataset: ``SpheroidAuv`` with a private ``Site``, ``SubmergedRestoring``,
    ``AuvHullLoads`` and four ``Fin`` (and a ``Propulsor`` with ``propulsor=True``). ``overrides``,
    ``hydrostatics``, ``fin`` and ``site`` are test-side parameter overrides (JSON-serialisable);
    ``hydrostatics_method`` replaces the part's buoyancy method. ``stand_in`` builds the same hull on what the
    plugin has today (Lamb added mass, homogeneous-spheroid inertia, the submerged form of the restoring):
    used only where the test is about that difference or does not depend on it."""
    key = json.dumps([part, propulsor, stand_in, overrides, hydrostatics, fin, site, hydrostatics_method],
                     sort_keys=True)
    if key in _CACHE:
        return _CACHE[key]
    data, p = DATA[part], vc.plugins()
    method = hydrostatics_method or ("from_volume" if stand_in else data.values("hydrostatics")["buoyancy_method"])
    if not stand_in:
        missing = missing_builds(p, propulsor=propulsor, method=method)
        if missing:
            pytest.fail(f"{part} does not build yet (waits on {'; '.join(WAITS[m] for m in missing)})")
    vehicle_values = {**data.values("vehicle"), **(overrides or {})}
    if stand_in:
        vehicle_values = {**vc.pick(p.SpheroidAuv, vehicle_values), "added_mass_form": "lamb_spheroid",
                          "roll_added_inertia_ratio": 0.0, **(overrides or {})}
    declared = {d.name for d in p.SubmergedRestoring.PARAMETERS}
    hydrostatics_values = {**{n: v for n, v in data.values("hydrostatics").items() if n in declared},
                           "buoyancy_method": method, **(hydrostatics or {})}
    actuators = [(p.Fin, {**data.fin_values(k), **(fin or {})}) for k in range(4)]
    if propulsor:
        from more_dynamics.plugins.force_producers.propulsor import Propulsor
        actuators.append((Propulsor, data.propulsor_values()))
    vessel = vc.build_context(p.SpheroidAuv, vehicle_values, site=(p.Site, {**data.values("site"), **(site or {})}),
                              hydrostatics=(p.SubmergedRestoring, hydrostatics_values),
                              hydrodynamics=(p.AuvHullLoads, data.values("hydrodynamics")), actuators=actuators)
    _CACHE[key] = vessel
    return vessel


def n_states(vessel):
    return vc.graph_of(vessel).step.size1_in(0)


def state(vessel, *, eta=None, nu=None, fin_deflections=None):
    """The stacked state: the children's states first (one servo lag per fin), then ``[eta; nu]``."""
    x = np.zeros(n_states(vessel))
    x[-12:-6] = np.zeros(6) if eta is None else eta
    x[-6:] = np.zeros(6) if nu is None else nu
    if fin_deflections is not None:
        x[:4] = fin_deflections
    return x


def inputs(vessel, command=None):
    u = np.zeros(vc.input_size(vessel))
    if command is not None:
        u[:len(command)] = command
    return u


def signal(vessel, name, x, u=None):
    return vessel.signals(x, inputs(vessel) if u is None else u)[name]


def restoring_force(vessel, eta):
    return signal(vessel, "hydrostatics.restoring_force", state(vessel, eta=eta))


def submerged_attitudes(data, seed, count=6):
    """Random poses deep enough that the whole hull is below the surface (depth >= one hull length)."""
    length = data.values("vehicle")["length"]
    rng = np.random.default_rng(seed)
    return [np.r_[rng.uniform(-1, 1, 2), rng.uniform(length, 3.0 * length), rng.uniform(-0.4, 0.4, 2),
                  rng.uniform(-1, 1)] for _ in range(count)]


# --------------------------------------------------------------------------------------------------------
# 1. Builds, the data files, generality
# --------------------------------------------------------------------------------------------------------
@each_part
def test_the_named_part_builds_with_four_fins_and_a_propulsor(part):
    vessel = vehicle(part, propulsor=True)
    assert type(vessel).__name__ == "SpheroidAuv"
    assert vc.input_size(vessel) == 4 + 1                     # four fin commands, the shaft speed
    assert n_states(vessel) == 4 + 12                         # four servo lags, the vehicle's twelve states


@each_part
def test_the_parameter_file_is_one_set_and_every_entry_has_a_kind_and_a_source(part):
    data = DATA[part]
    assert data.parameters["part"] == part and "parts" not in data.parameters      # one set, no variants inside
    assert data.derived["mat_sha256"] == MAT_SHA256
    assert data.values("hydrostatics")["buoyancy_method"] == "floating"
    kinds = {"published", "published_second_hand", "bem", "deployed_configuration", "derived", "convention",
             "estimate", "placeholder"}

    def walk(node, where):
        if isinstance(node, dict):
            if "value" in node:
                assert node["kind"] in kinds, (where, node)
                assert node["source"], (where, node)
                return 1
            return sum(walk(v, f"{where}.{k}") for k, v in node.items())
        if isinstance(node, list):
            return sum(walk(v, f"{where}[{i}]") for i, v in enumerate(node))
        return 0

    assert walk(data.parameters["parameters"], "parameters") > 60


def test_the_published_length_and_mass_of_the_named_part():
    """Length and mass are the vehicle's published numbers (NTNU AUR-Lab page, read 2026-10-10); weighing and
    taping as deployed is still the route to measured."""
    vehicle_group = DATA[MARIE].parameters["parameters"]["vehicle"]
    for name in ("length", "body_mass"):
        assert vehicle_group[name]["kind"] == "published"
        assert "aur-lab/lauv-marie" in vehicle_group[name]["source"] and "2026-10-10" in vehicle_group[name]["source"]
        assert vehicle_group[name]["route_to_measured"]


#: numeric values that a part legitimately shares with the reference torpedo's set: published or generic numbers
#: of the equation itself, a numerical guard, a water property
SHARED_WITH_THE_REFERENCE_SET = {
    "oswald_efficiency",        # generic published efficiency of a cylinder-shaped hull
    "planform_fraction",        # a rule applied to this hull's own length and diameter
    "kinematic_viscosity",      # a water property
    "integration_max_step", "smooth_speed", "smooth_speed_epsilon",   # numerical
    "gravity", "water_density",  # the site's: the same water as the reference set's still-water value
}


@each_part
def test_no_number_of_the_reference_torpedo_reaches_the_part(part):
    """Rule 20 of the library: the part gives every declared parameter of the four plugin types explicitly, and
    none of its numbers equals the reference torpedo's set (the frozen MSS set) unless it is a shared published or
    numerical value."""
    data, p = DATA[part], vc.plugins()
    reference = vc.remus_values()
    groups = {"vehicle": p.SpheroidAuv, "hydrostatics": p.SubmergedRestoring, "hydrodynamics": p.AuvHullLoads,
              "fin_common": p.Fin}
    for group, cls in groups.items():
        given = data.values(group)
        if group == "fin_common":
            given = {**given, "fin_position": None, "lift_axis": None}
        unset = sorted({d.name for d in cls.PARAMETERS} - set(given))
        assert not unset, (f"{cls.__name__}: the part must set {unset} itself (a named part holds its numbers and "
                           "forms), or the test's data file is incomplete")
        for name, value in data.values(group).items():
            if name in SHARED_WITH_THE_REFERENCE_SET or name not in reference:
                continue
            if isinstance(value, (str, bool)) or not np.any(np.asarray(value, dtype=float)):
                continue                                  # a selector, or an exact zero (a convention)
            assert not np.array_equal(np.asarray(value, float), np.asarray(reference[name], float)), (group, name, value)


def test_no_library_module_names_the_vehicle():
    """Rule 20: the plugins, the models and the plugin types know physics, not vehicles. ``grep -ri marie`` over the
    three folders finds nothing, in any file; the vehicle family's name appears only where a line cites a file
    (``LAUV_marie.mat``) or a paper."""
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


# --------------------------------------------------------------------------------------------------------
# 2. Added mass and rigid body
# --------------------------------------------------------------------------------------------------------
def _reader():
    if str(LIBRARY) not in sys.path:
        sys.path.insert(0, str(LIBRARY))
    return importlib.import_module("scripts.vessel_database.vessel_database_reader")


_READING = []


def _reading():
    if not _READING:
        reader = _reader()
        vessel = reader.load_vessel_mat(MAT_FILE)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _READING.append(reader.read_bem_vessel(MAT_FILE, added_mass_frequency="infinite", body_kind="submerged",
                                                   co_to_reference_point=np.asarray(vessel.main.CG, float),
                                                   generator_revision="cc07579"))
    return _READING[0]


def expected_added_mass(data):
    """The vessel-database reader's infinite-frequency matrix at the origin, under the dataset's similitude."""
    d = data.similitude()
    return d @ _reading().added_mass_matrix @ d


@each_part
def test_the_parts_added_mass_is_the_readers_matrix_under_its_similitude(part):
    data = DATA[part]
    given = np.asarray(data.values("vehicle")["added_mass_matrix"])
    assert np.abs(given - expected_added_mass(data)).max() <= 1e-12 * data.scale**5


@each_part
def test_added_mass_is_the_vessel_database_infinite_frequency_matrix_at_the_origin(part):
    data = DATA[part]
    vessel = vehicle(part)
    resolved = vc.matrix(signal(vessel, "added_mass_matrix", state(vessel)))
    assert np.abs(resolved - expected_added_mass(data)).max() <= 1e-12 * data.scale**5


@each_part
def test_rigid_body_mass_matrix_is_the_files_moved_to_the_origin(part):
    """The body the added mass was computed for: the uniform-cylinder inertia of the file, moved from the centre of
    gravity to the origin by H(r)^T M H(r) (Fossen 2011, eq. 3.24, p. 49)."""
    data = DATA[part]
    vessel = vehicle(part)
    resolved = vc.matrix(signal(vessel, "rigid_body_mass_matrix", state(vessel)))
    expected = np.asarray(data.derived["rigid_body_mass_matrix_at_origin"])
    assert np.abs(resolved - expected).max() <= 1e-12 * data.scale**5
    assert resolved[0, 0] == pytest.approx(data.derived["mass"], abs=1e-12 * data.scale**3)


@each_part
def test_the_lamb_spheroid_differs_from_the_bem_matrix_by_the_documented_ratios(part):
    """A documented difference, not an error: the equivalent spheroid of the same length and diameter is the wrong
    shape for this tapered cylinder. Runs today on the stand-in: the printed values are the spheroid's added mass
    at the origin and the ratios to the boundary-element matrix (the data file holds both, with their source)."""
    data = DATA[part]
    vessel = vehicle(part, stand_in=True)
    lamb = vc.matrix(signal(vessel, "added_mass_matrix", state(vessel)))
    bem = expected_added_mass(data)
    printed = data.derived["lamb_spheroid"]
    for index, value in printed["diagonal"].items():
        i = int(index)
        assert lamb[i, i] == pytest.approx(value, rel=1e-3)
        assert lamb[i, i] / bem[i, i] == pytest.approx(printed["ratio_to_bem"][index], abs=1.5e-3)
    assert lamb[3, 3] == pytest.approx(0.0, abs=1e-12)                              # Lamb: no roll added inertia


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
    """Fully below the surface: net heave force B - W (positively buoyant, never neutral) and the pitch moment of
    the offset of the centres, at several attitudes. ``x_g`` is the offset at Marie's scale."""
    data = DATA[part]
    weight, buoyancy = data.derived["weight"], data.derived["full_buoyancy"]
    r_g = [x_g * data.scale, 0.0, data.values("vehicle")["body_center_of_gravity"][2]]
    vessel = vehicle(part, overrides={"body_center_of_gravity": r_g})
    for eta in submerged_attitudes(data, 20261010):
        expected = -_g_submerged(eta, weight, buoyancy, r_g, [0.0, 0.0, 0.0])
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
    x = state(vessel, eta=[0, 0, 10.0, 0, 0, 0])
    u = inputs(vessel)
    for _ in range(5):
        x = vessel.step_vector(x, u, 1.0)
    assert x[-10] < 10.0 - 1e-3                    # depth (NED z) decreased
    assert x[-6 + 2] < 0.0                         # heave velocity upward


# --------------------------------------------------------------------------------------------------------
# 4. Floating equilibrium
# --------------------------------------------------------------------------------------------------------
ROW_CASES = [(part, row) for part in PARTS for row in DATA[part].rows]


@pytest.mark.parametrize("part,row", ROW_CASES, ids=[f"{IDS[p]}-{r['case']}" for p, r in ROW_CASES])
def test_floating_equilibrium_against_the_second_implementation(part, row):
    """At the second implementation's draft and pitch the part's restoring force and moment vanish: the quadrature
    of the vehicle's form agrees with the converged one to the stated tolerances. ``cg_x_m`` is the longitudinal
    offset of the centre of gravity from the centre of buoyancy; the mass sets the reserve buoyancy."""
    data = DATA[part]
    r_g = [row["cg_x_m"], 0.0, data.values("vehicle")["body_center_of_gravity"][2]]
    vessel = vehicle(part, overrides={"body_mass": row["mass_kg"], "body_center_of_gravity": r_g})
    eta = np.array([0.0, 0.0, row["origin_depth_m"], 0.0, row["pitch_rad"], 0.0])
    force = restoring_force(vessel, eta)
    assert abs(force[2]) <= QUADRATURE_FORCE * data.scale**3, force
    assert abs(force[4]) <= QUADRATURE_MOMENT * data.scale**4, force
    assert abs(force[0]) <= QUADRATURE_FORCE * data.scale**3
    step = 1e-2 * data.scale
    above, below = (restoring_force(vessel, eta + np.array([0, 0, h, 0, 0, 0]))[2] for h in (-step, step))
    assert above > force[2] > below                # heave: deeper lowers the net downward force (stiffness > 0)
    assert (above - below) / (2 * step) > 100.0 * data.scale**2     # N/m, a waterplane stiffness, not a flat curve


@each_part
def test_floating_depth_and_pitch_by_solving_the_vehicle_for_the_declared_part_and_the_offsets(part):
    """The same cells by solving the vehicle for equilibrium from a perturbed guess: depth and bow-up pitch."""
    from scipy import optimize
    data = DATA[part]
    for row in data.rows:
        if not row["case"].startswith(("declared", "reserve_0.02", "reserve_0.03", "reserve_0.04")):
            continue
        r_g = [row["cg_x_m"], 0.0, data.values("vehicle")["body_center_of_gravity"][2]]
        vessel = vehicle(part, overrides={"body_mass": row["mass_kg"], "body_center_of_gravity": r_g})

        def residual(v):
            f = restoring_force(vessel, np.array([0.0, 0.0, v[0], 0.0, v[1], 0.0]))
            return [f[2], f[4]]

        guess = [row["origin_depth_m"] + 5e-3 * data.scale, row["pitch_rad"] + 2e-2]
        solution = optimize.fsolve(residual, guess, xtol=1e-12)
        assert solution[0] == pytest.approx(row["origin_depth_m"],
                                            abs=(QUADRATURE_DEPTH + 1e-3 * abs(row["pitch_rad"])) * data.scale)
        assert solution[1] == pytest.approx(row["pitch_rad"], abs=3e-4 + 1e-3 * abs(row["pitch_rad"]))
        if row["cg_x_m"] < 0.0:
            assert solution[1] > 0.0               # the centre of gravity aft of the centre of buoyancy: nose up


@each_part
def test_the_floating_form_equals_the_volume_form_when_every_station_is_submerged(part):
    """Generality of the one part: fully below the surface the floating form is the submerged form with the volume
    of the same hull profile, to 1e-12 of the weight. This is why one part serves the dive and the surface; the
    volume form is a test-side override of the same part, not a second part."""
    data = DATA[part]
    floating, from_volume = vehicle(part), vehicle(part, hydrostatics_method="from_volume")
    for eta in submerged_attitudes(data, 7, count=8):
        difference = restoring_force(floating, eta) - restoring_force(from_volume, eta)
        assert np.abs(difference).max() <= 1e-12 * data.derived["weight"] * max(1.0, data.scale)


@each_part
def test_reserve_buoyancy_sweep_reported_against_the_pitch_it_gives(part):
    """Reported, not gated beyond the second implementation's own monotonic result: at the declared part's first
    aft offset of the centre of gravity the equilibrium pitch falls as the reserve grows from 1 % to 4 % (the
    sensitivity the owner was told about; the declared reserve is an estimate until the vehicle is weighed in
    water)."""
    rows = DATA[part].rows
    offset = next(r["cg_x_m"] for r in rows if r["case"].startswith("declared_dx_"))
    cells = [r for r in rows if r["case"].startswith("reserve_") and r["cg_x_m"] == offset]
    pitches = [r["pitch_rad"] / DEG for r in cells]
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
    radius = data.values("vehicle")["diameter"] / 2.0
    x = state(vessel, eta=[0, 0, 10.0, 0, 0, 0])
    u = inputs(vessel)                              # zero commands: fins neutral, no thrust (no propulsor fitted)
    depths, time_to_surface = [], None
    for second in range(1, 601):
        x = vessel.step_vector(x, u, 1.0)
        assert np.all(np.isfinite(x))
        depths.append(x[-10])
        if time_to_surface is None and x[-10] < radius:
            time_to_surface = second
    assert time_to_surface is not None, "the vehicle never broke the surface"
    assert min(depths[:time_to_surface]) < 10.0 - 1.0
    print(f"{part}: pop-up time to surface (reported, not gated): {time_to_surface} s; settled depth {depths[-1]:.5f} m")
    assert depths[-1] == pytest.approx(equilibrium["origin_depth_m"], abs=1e-3 * data.scale)
    assert x[-8] == pytest.approx(0.0, abs=1e-3)                                # pitch
    assert np.abs(x[-6:]).max() < 1e-3                                           # at rest


# --------------------------------------------------------------------------------------------------------
# 6. Speed table (reported against the configured table and the published top speed, gated against the
#    independent steady balance)
# --------------------------------------------------------------------------------------------------------
def _steady_speed(part, rpm):
    """Neutral buoyancy for this check only (the net buoyancy would tilt the flow and this is a test of the surge
    balance), solved as the speed at which the surge acceleration vanishes: open loop the hull is unstable in pitch
    and a long integration would only amplify the rounding."""
    from scipy import optimize
    vessel = vehicle(part, propulsor=True, hydrostatics_method="neutral")
    dynamics = vc.graph_of(vessel).step
    u = inputs(vessel, [0, 0, 0, 0, rpm])

    def surge_acceleration(speed):
        return float(np.asarray(dynamics(state(vessel, nu=[speed, 0, 0, 0, 0, 0]), u)).ravel()[-6])

    return optimize.brentq(surge_acceleration, 1e-3, 10.0, xtol=1e-12)


@each_part
@pytest.mark.parametrize("index", range(5))
def test_steady_surge_speed_equals_the_independent_balance_and_is_reported_against_the_configured_table(part, index):
    """Thrust of the lumped law against the hull's axial drag and the time-constant surge damping, fitted at one
    row of the table. The other rows are reported against the table: a configured table, not a measurement; the
    law is fitted later from the logs."""
    table = DATA[part].derived["speed_table"]
    rpm, table_speed, model_speed = table["rpm"][index], table["table_speed"][index], table["model_speed"][index]
    speed = _steady_speed(part, rpm)
    print(f"{part}: {rpm:.0f} rpm: vehicle {speed:.3f} m/s, configured table {table_speed:.2f} m/s (reported, not gated)")
    assert speed == pytest.approx(model_speed, abs=1e-4)


@each_part
def test_steady_surge_at_the_maximum_shaft_speed_is_reported_against_the_published_top_speed(part):
    """The published speed range's upper bound is the check's reference, reported, not gated (the thrust law is a
    placeholder); the gate is the independent balance."""
    table = DATA[part].derived["speed_table"]
    speed = _steady_speed(part, table["max_rpm"])
    print(f"{part}: {table['max_rpm']:.0f} rpm: vehicle {speed:.3f} m/s, published top speed "
          f"{table['published_max_speed']:.2f} m/s (reported, not gated)")
    assert speed == pytest.approx(table["max_model_speed"], abs=1e-4)


def _propeller_wrench(out):
    return next(v for k, v in out.items() if k.startswith("actuators.4.") and np.size(v) == 6)


@each_part
def test_zero_shaft_speed_gives_no_wrench(part):
    vessel = vehicle(part, propulsor=True)
    x = state(vessel, nu=[1.0, 0, 0, 0, 0, 0])
    assert np.abs(_propeller_wrench(vessel.signals(x, inputs(vessel, [0, 0, 0, 0, 0.0])))).max() == 0.0


@each_part
def test_propeller_wrench_is_the_lumped_law_along_the_shaft_at_the_declared_position(part):
    """tau = [T e_x; r x T e_x], T = k w |w|, w = 2 pi n (the lumped law, second implementation); the position is the
    declared propeller position; the command is capped at the maximum shaft speed."""
    values = DATA[part].propulsor_values()
    vessel = vehicle(part, propulsor=True)
    top = values["max_command"]
    for rpm in (0.13 * top, 0.5 * top, top):
        thrust = values["positive_thrust_coefficient"] * (2.0 * math.pi * rpm / 60.0) ** 2
        r = np.asarray(values["position"], float)
        expected = np.r_[[thrust, 0.0, 0.0], np.cross(r, [thrust, 0.0, 0.0])]
        wrench = _propeller_wrench(vessel.signals(state(vessel), inputs(vessel, [0, 0, 0, 0, rpm])))
        assert np.abs(wrench - expected).max() <= G1 * max(thrust, 1.0)
    capped = _propeller_wrench(vessel.signals(state(vessel), inputs(vessel, [0, 0, 0, 0, 2.0 * top])))
    assert capped[0] == pytest.approx(values["positive_thrust_coefficient"] * (2.0 * math.pi * top / 60.0) ** 2,
                                      rel=1e-9)


# --------------------------------------------------------------------------------------------------------
# 7. Fins: layout, lift, limits
# --------------------------------------------------------------------------------------------------------
def _fin_wrench(vessel, k, deflection, speed):
    deflections = np.zeros(4)
    deflections[k] = deflection
    out = vessel.signals(state(vessel, nu=[speed, 0, 0, 0, 0, 0], fin_deflections=deflections), inputs(vessel))
    return out[f"actuators.{k}.generated_force"]


@each_part
@pytest.mark.parametrize("k", range(4))
@pytest.mark.parametrize("speed", [0.8, 1.9])
def test_the_lift_part_of_each_fin_wrench_equals_the_second_implementation(part, k, speed):
    """Drag is even in the deflection and lift odd, so half the difference of the wrenches at +d and -d is the lift
    part, [f; r x f] with f = 1/2 rho A U^2 C_L,alpha d along the fin's lift axis (Fossen 2011, eq. 12.226, p. 400;
    pure surge, so the flow angle is the deflection)."""
    data = DATA[part]
    vessel = vehicle(part, stand_in=True)
    values = data.fin_values(k)
    rho = data.values("site")["water_density"]
    r = np.asarray(values["fin_position"], float)
    axis = np.asarray(values["lift_axis"], float)
    for delta in (0.05, 0.2, 0.4):
        odd = 0.5 * (_fin_wrench(vessel, k, delta, speed) - _fin_wrench(vessel, k, -delta, speed))
        lift = 0.5 * rho * values["fin_area"] * speed**2 * values["lift_slope"] * delta * axis
        expected = np.r_[lift, np.cross(r, lift)]
        assert np.abs(odd - expected).max() <= 1e-9 * max(1.0, np.abs(expected).max())


@each_part
def test_fin_layout_yaw_by_the_vertical_pair_pitch_by_the_horizontal_pair_roll_cancels(part):
    """The layout of the deployed simulator's fins (hypothesis: the vehicle's own layout is not documented): the two
    vertical fins give yaw and no roll or pitch when commanded alike, the two horizontal ones give pitch and no roll
    or yaw; differential horizontal commands give roll."""
    vessel = vehicle(part, stand_in=True)
    speed, delta = 1.5, 0.2

    def odd(ks, signs):
        total = np.zeros(6)
        for k, s in zip(ks, signs):
            total += s * 0.5 * (_fin_wrench(vessel, k, delta, speed) - _fin_wrench(vessel, k, -delta, speed))
        return total

    yaw = odd((0, 3), (1, 1))
    pitch = odd((1, 2), (1, 1))
    roll = odd((1, 2), (1, -1))
    assert abs(yaw[5]) > 1e-3 and abs(yaw[3]) < 1e-9 and abs(yaw[4]) < 1e-9
    assert abs(pitch[4]) > 1e-3 and abs(pitch[3]) < 1e-9 and abs(pitch[5]) < 1e-9
    assert abs(roll[3]) > 1e-3 and abs(roll[4]) < 1e-9


@each_part
def test_fin_deflection_limit_is_enforced_on_the_command(part):
    data = DATA[part]
    vessel = vehicle(part, stand_in=True)
    limit = data.values("fin_common")["max_deflection"]
    time_constant = data.values("fin_common")["time_constant"]
    graph = vc.graph_of(vessel).step
    x = state(vessel)
    over = np.asarray(graph(x, inputs(vessel, [1.0, 0, 0, 0]))).ravel()[0]
    at_limit = np.asarray(graph(x, inputs(vessel, [limit, 0, 0, 0]))).ravel()[0]
    assert over == pytest.approx(at_limit, abs=1e-12)
    assert over == pytest.approx(limit / time_constant, rel=1e-9)
    held = state(vessel)
    for _ in range(100):
        held = vessel.step_vector(held, inputs(vessel, [1.0, 0, 0, 0]), 0.1)
    assert held[0] <= limit + 1e-9


@each_part
def test_fin_rate_limit_is_enforced_through_the_servo(part):
    """With a time constant short enough that the lag alone would exceed it (a test-side override), the state's rate
    is the declared limit."""
    data = DATA[part]
    vessel = vehicle(part, stand_in=True, fin={"time_constant": 0.005})
    rate = data.values("fin_common")["max_rate"]
    graph = vc.graph_of(vessel).step
    got = np.asarray(graph(state(vessel), inputs(vessel, [0.3, 0, 0, 0]))).ravel()[0]
    assert got == pytest.approx(rate, rel=1e-9)
    assert rate < 0.3 / 0.005


# --------------------------------------------------------------------------------------------------------
# 8. Integration
# --------------------------------------------------------------------------------------------------------
@each_part
def test_the_vehicles_own_step_equals_a_test_side_rk4_of_the_graphs_dynamics(part):
    data = DATA[part]
    vessel = vehicle(part, propulsor=True)
    dynamics = vc.graph_of(vessel).step
    rng = np.random.default_rng(11)
    depth = data.values("vehicle")["length"]                  # fully submerged
    x = state(vessel, eta=np.r_[rng.uniform(-1, 1, 2), depth, rng.uniform(-0.2, 0.2, 2), rng.uniform(-1, 1)],
              nu=np.r_[1.0, rng.uniform(-0.1, 0.1, 5)], fin_deflections=rng.uniform(-0.1, 0.1, 4))
    u = inputs(vessel, [0.05, -0.05, 0.02, 0.0, 0.4 * data.propulsor_values()["max_command"]])
    reference = x.copy()
    for _ in range(10):                             # 0.5 s in sub-steps of the declared integration step
        k1 = np.asarray(dynamics(reference, u)).ravel()
        k2 = np.asarray(dynamics(reference + 0.5 * STEP * k1, u)).ravel()
        k3 = np.asarray(dynamics(reference + 0.5 * STEP * k2, u)).ravel()
        k4 = np.asarray(dynamics(reference + STEP * k3, u)).ravel()
        reference = reference + STEP / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
    assert np.abs(vessel.step_vector(x, u, 0.5) - reference).max() <= 1e-10


@each_part
def test_step_on_the_standard_odometry_with_actuator_states_is_not_contracted(part):
    """``step()`` of a vehicle whose actuators have states (the servo lags) waits on the vehicle-model type's
    maintainer (Luka): where those states live between calls is his type's question. Reported as an expected
    failure, not asserted; ``graph()`` and ``step_vector`` are gated above."""
    vessel = vehicle(part)
    with pytest.raises(NotImplementedError, match="actuator states"):
        vessel.step(None, [], 0.0, 0.1)
    pytest.xfail("step() with actuator states waits on Luka: where actuator states live between step() calls")
