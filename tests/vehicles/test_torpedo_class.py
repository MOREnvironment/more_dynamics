"""REMUS 100 gates for the torpedo-AUV tree of ``MarineCraft6DOF``.

The reference set is the marked one-density, one-geometry MSS copy.  States
in its CSVs are [nu; eta]; the vehicle graph and step use [eta; nu].  The
reordering is confined to these tests.  The tree's stage 1 producer is a
prescribed wrench; stage 2 uses fins and a propeller.  The vehicle's inputs
are the three open current parameters, then the commands of its force
producers in list order.

Author:    Enio Krizman
Date:      2026-10-07
"""

import csv
import inspect
import json
import os
from pathlib import Path
import importlib

import numpy as np
import pytest
from more_transformations.matrix_transforms import MatrixTransforms

import vehicle_contract as vc
from vehicle_contract import CURRENT_INPUTS, DATA, SLOTS, make_trees

DATA_DIR = DATA / "remus100"
MSS_DIR_VARIABLE = "MSS_DIR"
HULL, FULL = 1, 2
G1_TOLERANCE = 1e-9
STEP = 0.05  # SIMremus100.m:60
DEG = np.pi / 180.0
PARAMETER_FILE = "remus100_parameters_consistent.json"
DERIVATIVE_FILE = "remus100_derivative_consistent.csv"
HULL_DERIVATIVE_FILE = "remus100_hull_derivative_consistent.csv"
MASS_MATRIX_FILE = "remus100_mass_matrix_consistent.csv"
TOLERANCE_FILE = "remus100_trajectory_tolerances_consistent.csv"
TRAJECTORY_SUFFIX = "_consistent"
MSS_DERIVATIVE_FILE = "remus100_derivative_mss.csv"
MARKED_COPIES = {
    "remus100_hull.m": "CRAFT/AUV/models/remus100.m",
    "remus100_munk.m": "CRAFT/AUV/models/remus100.m",
    "remus100_hull_munk.m": "CRAFT/AUV/models/remus100.m",
    "remus100_consistent.m": "CRAFT/AUV/models/remus100.m",
    "spheroid_consistent.m": "LIBRARY/modeling/spheroid.m",
    "imlay61_consistent.m": "LIBRARY/modeling/imlay61.m",
    "forceLiftDrag_consistent.m": "LIBRARY/modeling/forceLiftDrag.m",
    "crossFlowDrag_consistent.m": "LIBRARY/modeling/crossFlowDrag.m",
}


def _read_csv(name):
    with open(DATA_DIR / name, newline="") as handle:
        rows = list(csv.reader(handle))
    return {name: np.asarray(rows[1:], float)[:, i] for i, name in enumerate(rows[0])}


def _columns(table, prefix, count):
    return np.column_stack([table[f"{prefix}_{i:02d}"] for i in range(1, count + 1)])


def _parameter_file():
    return json.loads((DATA_DIR / PARAMETER_FILE).read_text())


def _parameter_set(stage):
    entries = _parameter_file()["parameters"]
    keep = ("hull",) if stage == HULL else ("hull", "actuators")
    return {n: e["value"] for n, e in entries.items() if e["stage"] in keep}


def _trajectory_tolerances():
    table = {}
    with open(DATA_DIR / TOLERANCE_FILE, newline="") as handle:
        for row in csv.DictReader(handle):
            table.setdefault(row["scenario"], np.zeros(12))[int(row["state"]) - 1] = float(row["tolerance"])
    return table


_mss_to_ours = vc.mss_to_ours
_ours_to_mss = vc.ours_to_mss


def _tree(stage, munk_kept=False, values=None):
    """The generator's tree of the stage with the chosen added-mass Coriolis part and parameter set."""
    part = "KirchhoffFull" if munk_kept else "MunkCouplingsRemoved"
    p = _parameter_set(stage) if values is None else values
    return make_trees.torpedo(part, full=stage == FULL, values=p)


_BUILT = {}


def _builder(stage, munk_kept=False, values=None, tmp_path=None):
    """The vehicle of the stage: the committed gate tree when the parameters are the frozen set, otherwise the
    perturbed set written into ``tmp_path``."""
    if values is None:
        name = f"remus100_{'full' if stage == FULL else 'hull'}{'_kirchhoff' if munk_kept else ''}"
        if name not in _BUILT:
            _BUILT[name] = vc.build(name)
        return _BUILT[name]
    return vc.build_variant(tmp_path, _tree(stage, munk_kept, values))


def _inputs(vehicle, current, command):
    """The vehicle's stacked input rows: the open current parameters, then the commands."""
    return np.column_stack([current, command]).T


def _evaluate(vehicle, x_mss, current, tau_ext=None, u=None):
    n = len(x_mss)
    size = vc.input_size(vehicle) - 3
    command = np.zeros((n, size))
    if tau_ext is not None and size == 6:
        command = tau_ext
    elif u is not None and size == 3:
        command = u
    xdot = vc.dynamics_map(vehicle, n)(_mss_to_ours(x_mss).T, _inputs(vehicle, current, command))
    mass = vc.matrix(vc.output_entries(vehicle, np.zeros(12), np.zeros(3 + size))["mass_matrix"])
    return _ours_to_mss(np.asarray(xdot).T), mass


def _integrate(vehicle, table):
    times = table["t"]
    steps = int(round((times[1] - times[0]) / STEP))
    has_u = "ui_01" in table
    commands = _columns(table, "ui" if has_u else "tau_ext", 3 if has_u else 6)
    current = np.column_stack([table["Vc"], table["betaVc"], table["w_c"]])
    states = np.zeros((len(times), 12))
    x = _mss_to_ours(_columns(table, "x", 12)[0])
    for row in range(len(times)):
        states[row] = _ours_to_mss(x)
        if row < len(times) - 1:
            u = np.concatenate([current[row], commands[row]])
            for _ in range(steps):
                x = vehicle.step_vector(x, u, STEP)
    return states


def _scenario(name):
    if name.startswith("hull_munk"):
        return HULL, True
    if name.startswith("hull"):
        return HULL, False
    if name.startswith("munk"):
        return FULL, True
    return FULL, False


TRAJECTORIES = sorted(
    p.stem[len("remus100_trajectory_"):-len(TRAJECTORY_SUFFIX)]
    for p in DATA_DIR.glob(f"remus100_trajectory_*{TRAJECTORY_SUFFIX}.csv")
    if not p.stem.startswith("remus100_trajectory_tolerances"))


def _trajectory(name):
    return _read_csv(f"remus100_trajectory_{name}{TRAJECTORY_SUFFIX}.csv")


SOURCE_PARAMETER_HEADING = "## Parameter set — one value per quantity (the default)"


def _source_section(heading):
    source = (DATA_DIR / "SOURCE.md").read_text()
    start = source.index(heading + "\n")
    end = source.find("\n## ", start + len(heading))
    return source[start:] if end < 0 else source[start:end]

def test_parameter_set_names_are_block_declarations_and_have_a_source_line():
    """Every value of the parameter set: a block of this library declares it
    (the declared-parameter form, or the arguments of a block still in the
    earlier form), and ``SOURCE.md`` has its line with the same value."""
    entries = _parameter_file()["parameters"]
    source = _source_section(SOURCE_PARAMETER_HEADING)
    modules = {
        "rigid_body": "more_dynamics.models.vehicles.hull_parts.rigid_body",
        "submerged_hydrostatics": "more_dynamics.models.vehicles.hull_parts.restoring.submerged",
        "submerged_linear_damping": "more_dynamics.models.vehicles.hull_parts.hydrodynamic_loads.damping.linear_damping",
        "lift_drag": "more_dynamics.models.vehicles.hull_parts.hydrodynamic_loads.lift_drag.lift_drag",
        "cross_flow": "more_dynamics.models.vehicles.hull_parts.hydrodynamic_loads.cross_flow.cross_flow",
        "propeller": "more_dynamics.models.force_producers.propulsor.propeller",
        "fins": "more_dynamics.models.force_producers.fins",
    }
    selectors = _parameter_file()["selectors"]

    def declared(block):
        """Names a block declares: ``*_parameters()`` (called with the
        selectors it accepts) or the arguments of ``preprocess_*``."""
        module = importlib.import_module(modules[block])
        names = set()
        for attribute, function in vars(module).items():
            if not inspect.isfunction(function):
                continue
            if attribute.endswith("_parameters"):
                accepted = inspect.signature(function).parameters
                try:
                    result = function(**{k: v for k, v in selectors.items() if k in accepted})
                except TypeError:
                    continue
                names |= {p.name for p in result}
            elif attribute.startswith("preprocess_"):
                names |= set(inspect.signature(function).parameters)
        return names

    problems = []
    for name, entry in entries.items():
        blocks = [b.strip().split(" ")[0] for b in entry["block"].split(",")]
        if not any(entry["block_parameter"] in declared(b) for b in blocks):
            problems.append(f"{name}: no block of {blocks} declares {entry['block_parameter']!r}")
        line = [ln for ln in source.splitlines() if ln.startswith(f"| `{name}` |")]
        if len(line) != 1:
            problems.append(f"{name}: {len(line)} lines in SOURCE.md (one expected)")
            continue
        printed = json.loads(line[0].split("|")[2].strip().strip("`"))
        if not np.array_equal(np.asarray(printed, dtype=float), np.asarray(entry["value"], dtype=float)):
            problems.append(f"{name}: SOURCE.md prints {printed}, the parameter file holds {entry['value']}")
    assert not problems, "\n".join(problems)


def test_parameter_set_holds_one_value_per_quantity():
    """One water density, one geometry, the mass given: the body density
    gives the given mass on the hull's own semi-axes, ``W = m g = B``
    (``remus100.m`` 214), the fins sit at ``-a``, and the lift/drag and
    cross-flow use the same length and diameter as the body."""
    data = _parameter_file()
    p = {name: e["value"] for name, e in data["parameters"].items()}
    assert not any("density" in name and name not in ("water_density", "body_density") for name in p)
    a, b = p["semi_major_axis"], p["semi_minor_axis"]
    mass = data["one_value"]["mass"]
    assert p["body_density"] * 4.0 / 3.0 * np.pi * a * b**2 == pytest.approx(mass, rel=1e-15)
    assert p["weight"] == p["buoyancy"] == pytest.approx(mass * data["gravity"]["g"], rel=1e-15)
    assert p["length"] == 2 * a and p["span"] == p["beam"] == p["draft"] == 2 * b
    assert p["rudder_position"] == p["stern_plane_position"] == -a
    assert p["planform_area"] == pytest.approx(0.7 * p["length"] * p["span"], rel=1e-15)  # remus100.m 133


def test_history_unmodified_mss_differs_by_one_density_and_one_geometry():
    """The unmodified MSS reference, kept as history, is not this class's
    reference, and by how much (measured by the generator, ``SOURCE.md``):

    * density alone (MSS geometry, one density; column ``xdot_one_density``):
      ~0.19 % of the accelerations (``max |d nudot| / max |nudot|`` over the
      640 cases), and -0.18 deg of heading after the 60 s turn at 15 deg
      rudder and 925 rpm (two water densities in one MSS model);
    * geometry and the given mass (one density -> one value per quantity):
      ~6 % of the accelerations and +2.9 deg of heading in the same turn.

    Bounds are the measured numbers rounded outward to their order."""
    mss, one = _read_csv(MSS_DERIVATIVE_FILE), _read_csv(DERIVATIVE_FILE)
    for prefix, count in (("x", 12), ("ui", 3)):
        assert np.array_equal(_columns(mss, prefix, count), _columns(one, prefix, count))

    def relative(a, b):
        scale = np.abs(a[:, 0:6]).max(axis=1)
        moving = scale > 0
        return (np.abs(b[moving, 0:6] - a[moving, 0:6]).max(axis=1) / scale[moving]).max()

    reference = _columns(mss, "xdot", 12)
    one_density = _columns(one, "xdot_one_density", 12)
    consistent = _columns(one, "xdot", 12)
    assert 1e-3 < relative(reference, one_density) < 3e-3
    assert 3e-2 < relative(one_density, consistent) < 1e-1
    assert np.abs(consistent - reference).max() > 1e3 * G1_TOLERANCE

    heading = {model: _read_csv(f"remus100_trajectory_rudder_step_15deg_at_925rpm_{model}.csv")["x_12"][-1]
               for model in ("mss", "one_density", "consistent")}
    density_change = (heading["one_density"] - heading["mss"]) / DEG
    geometry_change = (heading["consistent"] - heading["one_density"]) / DEG
    assert -0.3 < density_change < -0.1, f"density: {density_change:.4f} deg"
    assert 2.0 < geometry_change < 4.0, f"geometry: {geometry_change:.4f} deg"


def test_marked_copies_rebuild_their_mss_files():
    """Every marked copy differs from its MSS file only in its marked lines
    (drop ``%<added>`` lines, un-comment ``%<removed> ``)."""
    mss = os.environ.get(MSS_DIR_VARIABLE)
    if not mss:
        pytest.skip(f"{MSS_DIR_VARIABLE} is not set (path to the MSS checkout); copy check skipped")
    for copy, original_file in MARKED_COPIES.items():
        original = (Path(mss) / original_file).read_text()
        kept = []
        for line in (DATA_DIR / copy).read_text().split("\n"):
            if line.endswith("%<added>"):
                continue
            kept.append(line[len("%<removed> "):] if line.startswith("%<removed> ") else line)
        assert "\n".join(kept) == original, f"{copy} does not rebuild {MSS_DIR_VARIABLE}/{original_file}"


# ---------------------------------------------------------------------------
# Contract
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("stage", [HULL, FULL], ids=["hull", "full"])
def test_tree_holds_only_part_choices_their_parameters_and_the_vehicle_settings(stage):
    tree = _tree(stage)
    assert tuple(tree["children"]) == SLOTS
    assert set(tree["params"]) == {"integration_max_step", "open_inputs", "diagnostic_outputs"}
    import more_dynamics
    assert not (Path(more_dynamics.__file__).parent / "models" / "vehicles" / "torpedo_auv.py").exists()


def test_committed_gate_trees_are_the_generator_output(tmp_path):
    """No hand edit: every file of ``tests/vehicles/data/.rppws`` is what ``make_trees.py`` writes."""
    import filecmp
    pytest.importorskip("rpp_plugin_registrator.plugin_descriptors.core", reason="rpp is not installed")
    trees = make_trees.gate_trees()
    make_trees.write_workspace(tmp_path / ".rppws", trees, "tests/vehicles/vehicle_contract.py",
                               {"vehicles": list(trees)}, parts=make_trees.gate_parts())
    comparison = filecmp.dircmp(tmp_path / ".rppws", vc.GATE_TREES)

    def differences(cmp):
        found = cmp.left_only + cmp.right_only + cmp.diff_files
        for sub in cmp.subdirs.values():
            found += differences(sub)
        return found

    assert not differences(comparison)


@pytest.mark.parametrize("stage", [HULL, FULL], ids=["hull", "full"])
def test_tree_has_one_owner_for_each_parameter(stage):
    tree = _tree(stage)
    owners = {}
    for slot, item in tree["children"].items():
        for entry in item if isinstance(item, list) else [item]:
            for name in entry["params"]:
                owners.setdefault(name, []).append(slot)
    assert owners["water_density"] == ["site"]
    assert owners["kinematic_viscosity"] == ["site"]
    assert owners["length"] == ["hull_form"]
    assert owners["body_density"] == ["rigid_body"]


@pytest.mark.parametrize("stage, file, input_kind", [
    (FULL, DERIVATIVE_FILE, "ui"),
    (HULL, DERIVATIVE_FILE, "tau"),
    (HULL, HULL_DERIVATIVE_FILE, "tau_ext"),
], ids=["stage2", "stage1_actuator_wrench", "stage1_hull"])
@pytest.mark.parametrize("munk_kept, column", [(False, "xdot"), (True, "xdot_munk")],
                         ids=["munk_removed", "kirchhoff_full"])
def test_G1_derivative_equals_remus100_one_value(stage, file, input_kind, munk_kept, column):
    table = _read_csv(file)
    current = np.column_stack([table["Vc"], table["betaVc"], table["w_c"]])
    inputs = _columns(table, input_kind, 3 if input_kind == "ui" else 6)
    kwargs = {"u" if input_kind == "ui" else "tau_ext": inputs}
    result, _ = _evaluate(_builder(stage, munk_kept), _columns(table, "x", 12), current, **kwargs)
    error = np.abs(result - _columns(table, column, 12))
    assert error.max() <= G1_TOLERANCE, f"max |diff| {error.max():.3g}"


def test_G1_mass_matrix_equals_remus100_one_value_empty_call():
    _, mass = _evaluate(_builder(HULL), np.zeros((1, 12)), np.zeros((1, 3)))
    table = _read_csv(MASS_MATRIX_FILE)
    reference = np.column_stack([table[f"M_{j}"] for j in range(1, 7)])
    assert np.abs(mass - reference).max() <= G1_TOLERANCE


def test_added_mass_coriolis_parts_differ_only_in_acceleration():
    table = _read_csv(DERIVATIVE_FILE)
    current = np.column_stack([table["Vc"], table["betaVc"], table["w_c"]])
    x = _columns(table, "x", 12)
    tau = _columns(table, "tau", 6)
    a = _evaluate(_builder(HULL, False), x, current, tau_ext=tau)[0]
    b = _evaluate(_builder(HULL, True), x, current, tau_ext=tau)[0]
    assert np.abs(a[:, 0:6] - b[:, 0:6]).max() > 10 * G1_TOLERANCE
    assert np.array_equal(a[:, 6:12], b[:, 6:12])


@pytest.mark.parametrize("name", TRAJECTORIES)
def test_G3_trajectory_equals_matlab_rk4_run(name):
    stage, munk_kept = _scenario(name)
    table = _trajectory(name)
    states = _integrate(_builder(stage, munk_kept), table)
    error = np.abs(states - _columns(table, "x", 12)).max(axis=0)
    tolerance = _trajectory_tolerances()[name]
    over = [f"x{i + 1}: {error[i]:.3g} > {tolerance[i]:.3g}"
            for i in range(12) if error[i] > tolerance[i]]
    assert not over, f"{name}: " + "; ".join(over)


PERTURBATIONS = [
    ("body_density", lambda v: v * 1.01),
    ("body_center_of_gravity", lambda v: [v[0], v[1], v[2] + 0.01]),
    ("time_constants", lambda v: [2 * v[0], v[1], v[2]]),
    ("water_density", lambda v: v + 1.0),
    ("oswald_efficiency", lambda v: v * 1.01),
]


@pytest.mark.parametrize("name, change", PERTURBATIONS, ids=[p[0] for p in PERTURBATIONS])
def test_G4_perturbed_parameter_is_caught_by_G1(name, change, tmp_path):
    p = _parameter_set(HULL)
    p[name] = change(p[name])
    table = _read_csv(DERIVATIVE_FILE)
    current = np.column_stack([table["Vc"], table["betaVc"], table["w_c"]])
    xdot, _ = _evaluate(_builder(HULL, values=p, tmp_path=tmp_path), _columns(table, "x", 12), current,
                        tau_ext=_columns(table, "tau", 6))
    assert np.abs(xdot - _columns(table, "xdot", 12)).max() > 10 * G1_TOLERANCE


def test_G4_perturbed_center_of_gravity_is_caught_by_G3(tmp_path):
    p = _parameter_set(HULL)
    cg = p["body_center_of_gravity"]
    p["body_center_of_gravity"] = [cg[0], cg[1], cg[2] + 0.01]
    name = "hull_roll_pitch_decay"
    table = _trajectory(name)
    error = np.abs(_integrate(_builder(HULL, values=p, tmp_path=tmp_path), table)
                   - _columns(table, "x", 12)).max(axis=0)
    assert (error > 10 * _trajectory_tolerances()[name]).any()


def test_G4_perturbed_fin_coefficient_is_caught_by_G1(tmp_path):
    p = _parameter_set(FULL)
    p["rudder_lift_coefficient"] *= 1.01
    table = _read_csv(DERIVATIVE_FILE)
    current = np.column_stack([table["Vc"], table["betaVc"], table["w_c"]])
    xdot, _ = _evaluate(_builder(FULL, values=p, tmp_path=tmp_path), _columns(table, "x", 12), current,
                        u=_columns(table, "ui", 3))
    assert np.abs(xdot - _columns(table, "xdot", 12)).max() > 10 * G1_TOLERANCE


@pytest.mark.xfail(strict=True, raises=AssertionError,
                   reason="MSS REMUS turns near 1.9 deg/s; Prestero 2001 p. 53 prints about 10 deg/s")
def test_G5_prestero_turn_rate_at_4deg_rudder_and_1p54_m_per_s():
    table = _trajectory("rudder_step_4deg_at_1p54ms")
    states = _integrate(_builder(FULL), table)
    last = table["t"] >= table["t"][-1] - 10.0
    rate = np.mean(np.abs(states[last, 5])) / DEG
    assert abs(rate - 10.0) <= 5.0, f"steady yaw rate {rate:.3g} deg/s"


@pytest.mark.parametrize("stage", [HULL, FULL], ids=["hull", "full"])
def test_one_density_reaches_every_block(stage, tmp_path):
    p = _parameter_set(stage)
    table = _read_csv(DERIVATIVE_FILE)
    x = _columns(table, "x", 12)
    current = np.column_stack([table["Vc"], table["betaVc"], table["w_c"]])
    tau = _columns(table, "tau", 6) if stage == HULL else None
    u = _columns(table, "ui", 3) if stage == FULL else None
    one = _evaluate(_builder(stage), x, current, tau_ext=tau, u=u)[0]
    scaled = dict(p)
    for key in ("water_density", "body_density"):
        scaled[key] *= 2
    two = _evaluate(_builder(stage, values=scaled, tmp_path=tmp_path), x, current,
                    tau_ext=(2 * tau if tau is not None else None), u=u)[0]
    scale = np.maximum(np.abs(one).max(axis=1, keepdims=True), 1.0)
    assert (np.abs(two - one) / scale).max() <= 1e-12


@pytest.mark.parametrize("stage", [HULL, FULL], ids=["hull", "full"])
@pytest.mark.parametrize("munk_kept", [False, True], ids=["munk_removed", "kirchhoff_full"])
def test_neutrally_buoyant_body_at_rest_stays_at_rest(stage, munk_kept):
    xdot, _ = _evaluate(_builder(stage, munk_kept), np.zeros((1, 12)), np.zeros((1, 3)))
    assert np.abs(xdot).max() <= 1e-12


@pytest.mark.parametrize("munk_kept", [False, True], ids=["munk_removed", "kirchhoff_full"])
def test_energy_never_increases_without_input_or_current(munk_kept):
    p = _parameter_set(HULL)
    weight = _parameter_file()["one_value"]["mass"] * _parameter_file()["gravity"]["g"]
    r_g = np.asarray(p["body_center_of_gravity"])
    r_b = np.asarray(p["center_of_buoyancy"])
    rng = np.random.default_rng(20261007)
    n = 1000
    x = np.column_stack([rng.uniform(-1, 3, n), rng.uniform(-1, 1, (n, 2)),
                         rng.uniform(-0.5, 0.5, (n, 3)), rng.uniform(-10, 10, (n, 3)),
                         rng.uniform(-np.pi, np.pi, n), rng.uniform(-1.2, 1.2, n),
                         rng.uniform(-np.pi, np.pi, n)])
    xdot, mass = _evaluate(_builder(HULL, munk_kept), x, np.zeros((n, 3)))
    worst = -np.inf
    for k in range(n):
        nu, nu_dot = x[k, :6], xdot[k, :6]
        rotation = MatrixTransforms.Rzyx_explicit(x[k, 9:12])
        power = nu @ mass @ nu_dot + (rotation @ MatrixTransforms.skew(nu[3:6]) @
                                      (weight * r_b - weight * r_g))[2]
        scale = 1.0 + abs(nu @ mass @ nu_dot)
        worst = max(worst, power / scale)
    assert worst <= 1e-9, f"energy grows: dE/dt / scale = {worst:.3g}"


@pytest.mark.parametrize("munk_kept", [False, True], ids=["munk_removed", "kirchhoff_full"])
def test_surge_only_motion_produces_no_sway_roll_or_yaw(munk_kept):
    rng = np.random.default_rng(7)
    x = np.zeros((50, 12))
    x[:, 0] = rng.uniform(-1, 3, 50)
    x[:, 6:9] = rng.uniform(-10, 10, (50, 3))
    x[:, 11] = rng.uniform(-np.pi, np.pi, 50)
    tau = np.zeros((50, 6))
    tau[:, 0] = rng.uniform(-50, 50, 50)
    xdot, _ = _evaluate(_builder(HULL, munk_kept), x, np.zeros((50, 3)), tau_ext=tau)
    assert np.abs(xdot[:, [1, 3, 5]]).max() <= 1e-12


def test_default_keeps_added_mass_coriolis_when_fins_are_modelled():
    """The generator's default Coriolis part is the physics form, and the user vehicle (fins and propeller) uses it:
    the user tree and the gate tree with ``KirchhoffFull`` give the same derivative."""
    default = make_trees.torpedo(full=True)
    assert default["children"]["added_mass_coriolis"]["plugin"] == "KirchhoffFull"
    assert make_trees.user_trees()["remus100"]["children"]["added_mass_coriolis"]["plugin"] == "KirchhoffFull"
    table = _read_csv(DERIVATIVE_FILE)
    x = _columns(table, "x", 12)[:20]
    u = _columns(table, "ui", 3)[:20]
    current = np.column_stack([table["Vc"], table["betaVc"], table["w_c"]])[:20]
    kept = _evaluate(_builder(FULL, True), x, current, u=u)[0]
    assert np.abs(kept - _columns(table, "xdot_munk", 12)[:20]).max() <= G1_TOLERANCE


def test_a_parameter_the_tree_leaves_out_takes_the_plugin_default(tmp_path):
    """Every part declares its parameters with defaults (the frozen REMUS 100 set, cited in the part): a tree whose
    parts carry no ``parameters.py`` values gives the same vehicle."""
    tree = _tree(HULL)
    for item in tree["children"].values():
        for entry in item if isinstance(item, list) else [item]:
            entry["params"] = {k: v for k, v in entry["params"].items() if k == "open_parameters"}
            for child in entry["children"].values():
                child["params"] = {}
    table = _read_csv(HULL_DERIVATIVE_FILE)
    current = np.column_stack([table["Vc"], table["betaVc"], table["w_c"]])[:50]
    x = _columns(table, "x", 12)[:50]
    tau = _columns(table, "tau_ext", 6)[:50]
    bare = _evaluate(vc.build_variant(tmp_path, tree), x, current, tau_ext=tau)[0]
    assert np.abs(bare - _columns(table, "xdot", 12)[:50]).max() <= G1_TOLERANCE


def test_G4_perturbed_weight_through_site_gravity_is_caught_by_G1(tmp_path):
    tree = _tree(HULL)
    p = _parameter_set(HULL)
    site = tree["children"]["site"]
    tree["children"]["site"] = make_trees.node("SiteGiven", {
        "gravity": 1.01 * _parameter_file()["gravity"]["g"], "water_density": p["water_density"],
        "kinematic_viscosity": site["params"]["kinematic_viscosity"]})
    table = _read_csv(DERIVATIVE_FILE)
    current = np.column_stack([table["Vc"], table["betaVc"], table["w_c"]])
    xdot, _ = _evaluate(vc.build_variant(tmp_path, tree), _columns(table, "x", 12), current,
                        tau_ext=_columns(table, "tau", 6))
    assert np.abs(xdot - _columns(table, "xdot", 12)).max() > 10 * G1_TOLERANCE
