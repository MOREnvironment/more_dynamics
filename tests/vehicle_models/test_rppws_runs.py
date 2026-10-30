"""The frozen REMUS 100 runs through the library's ``.rppws``: the vehicle is
built by rpp's own ``ComponentContextBuilder`` from the named part ``REMUS 100
(MSS)`` (configuration ``remus100`` of ``script_descriptions/vehicle_simulation.json``,
through the ``VehicleSimulation`` of ``scripts/vehicles/vehicle_simulation.py``),
not assembled by ``vehicle_contract`` from the frozen parameter file as
``test_spheroid_auv.py`` does. Three gates, one run per test:

* the vehicle of the part, stepped with its own RK4 at the frozen step and
  commanded through its input descriptions, equals the MATLAB trajectory with
  the tolerance of ``test_spheroid_auv.test_G3_trajectory_equals_matlab_rk4_run``
  (the same table, read by the same function);
* the command schedule of the script description is the frozen run's ``ui``
  columns, sample by sample;
* the script's CVODES run of that schedule agrees with the RK4 run within the
  bound derived below.

A frozen run belongs here when the part's own selectors reproduce it: the full
vehicle (four fins, rudder pair and stern-plane pair, and the propulsor), the
added-mass Coriolis without the Munk couplings, still water. The runs of the bare hull and of the Munk-kept form are
other compositions of the plugin, and the run in a current needs a part with a
current form; ``EXCLUDED`` names the one that is left out and why.

Integrator bound. CVODES works with ``abstol`` 1e-8 and ``reltol`` 1e-6
(``more_simulation/simulation.py``, ``_simulate``) and keeps the weighted
root-mean-square norm of the error of each sample interval under 1, the weight of
state ``i`` being ``abstol + reltol * |x_i|``; that norm allows one state of ``m``
states to take ``sqrt(m)`` times its weight. With ``N`` intervals, the error of
``x_i`` of scale ``s_i = max |x_i|`` is therefore bounded by
``sqrt(m) * N * (abstol + reltol * s_i)`` if every interval used its whole
allowance (``m`` = 12 states of the vehicle). The RK4 run (step 0.05 s, the frozen
step) matches MATLAB to 1e-14, far under that bound, so the bound is one on the
CVODES run. (A first bound without the ``sqrt(m)`` factor, written before the
runs were compared, failed on three runs by up to 1.45 times; the factor is
CVODES's own norm, not a fit.)

Author:    Enio Krizman
Date:      2026-10-16
"""

import importlib.util
import inspect
import json
import re
from pathlib import Path

import numpy as np
import pytest

import vehicle_contract as vc
import test_spheroid_auv as g3

LIBRARY = Path(__file__).resolve().parents[2]
SCRIPT = LIBRARY / "scripts" / "vehicles" / "vehicle_simulation.py"
DESCRIPTION = LIBRARY / ".rppws" / "script_descriptions" / "vehicle_simulation.json"
CONFIGURATION = "remus100"
PART = "REMUS 100 (MSS)"
# the frozen file's command columns (MSS order: rudder, stern plane, shaft speed) as the vehicle's named inputs: the
# four fins in the order of the vehicle's actuators (rudder, rudder, stern plane, stern plane), then the propulsor
FROZEN_COMMAND_COLUMNS = {"fin_angle_command": ("ui_01", "ui_01", "ui_02", "ui_02"), "drive_command": ("ui_03",)}
CVODES_ABSTOL, CVODES_RELTOL = 1e-8, 1e-6  # more_simulation/simulation.py, _simulate (checked below)
CURRENT_COLUMNS = ("Vc", "betaVc", "w_c")


def _has_current(name):
    table = g3._trajectory(name)
    return any(np.abs(table[c]).max() > 0.0 for c in CURRENT_COLUMNS)


FULL_RUNS = [n for n in g3.TRAJECTORIES if g3._scenario(n) == (g3.FULL, False)]
EXCLUDED = {n: "needs the three current inputs (current_form 'yaw_rate_terms'); the part's current_form is 'none'"
            for n in FULL_RUNS if _has_current(n)}
RUNS = [n for n in FULL_RUNS if n not in EXCLUDED]

_CACHE = {}


def _script():
    """The script module, loaded from its file (it is not a package)."""
    if "module" not in _CACHE:
        vc.plugins()  # skips with the registry message when rpp's generated types are missing
        try:
            spec = importlib.util.spec_from_file_location("vehicle_simulation_script", SCRIPT)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
        except Exception as exc:  # no rpp, or more_simulation is not installed
            pytest.skip(f"{vc.REGISTER_MESSAGE} ({type(exc).__name__}: {str(exc)[:80]})")
        _CACHE["module"] = module
    return _CACHE["module"]


def _simulation(command=None, configuration=CONFIGURATION, **options):
    """A ``VehicleSimulation`` built by rpp's builder from the named part (``command`` a schedule name)."""
    module = _script()
    try:
        return module.VehicleSimulation(rpp_workspace=str(LIBRARY), rpp_configuration=configuration,
                                        command=command, **options)
    except ValueError:
        raise
    except Exception as exc:  # more_dynamics is not registered here
        pytest.skip(f"{vc.REGISTER_MESSAGE} ({type(exc).__name__}: {str(exc)[:80]})")


def _vehicle(configuration=CONFIGURATION):
    if ("vehicle", configuration) not in _CACHE:
        _CACHE[("vehicle", configuration)] = _simulation(configuration=configuration).vessels[0]
    return _CACHE[("vehicle", configuration)]


def _frozen_inputs(vehicle, table):
    """The frozen file's commands as the vehicle's stacked inputs, assembled in the order of its input
    descriptions by name; several inputs of one name (the four fins) take that name's columns one after the other
    (a name the table has no column for, or a column no input has, is an error)."""
    descriptions = list(vehicle.getInputDescriptions())
    assert sorted(set(d.name for d in descriptions)) == sorted(FROZEN_COMMAND_COLUMNS), [d.name for d in descriptions]
    columns, taken = [], {}
    for d in descriptions:
        start = taken.get(d.name, 0)
        part = FROZEN_COMMAND_COLUMNS[d.name][start:start + d.size]
        assert len(part) == d.size, d.name
        taken[d.name] = start + d.size
        columns += part
    assert taken == {name: len(c) for name, c in FROZEN_COMMAND_COLUMNS.items()}
    return np.column_stack([table[c] for c in columns])


def _rk4_run(vehicle, table):
    """The frozen run stepped by the vehicle's own RK4 at the frozen step, states in MSS order [nu; eta]."""
    inputs = _frozen_inputs(vehicle, table)
    steps = int(round((table["t"][1] - table["t"][0]) / g3.STEP))
    x = g3._mss_to_ours(g3._columns(table, "x", 12)[0])
    states = np.zeros((len(table["t"]), 12))
    for row in range(len(table["t"])):
        states[row] = g3._ours_to_mss(x)
        if row < len(table["t"]) - 1:
            for _ in range(steps):
                x = vehicle.step_vector(x, inputs[row], g3.STEP)
    return states


def test_the_runs_of_the_part_are_all_the_full_vehicle_runs_but_the_one_in_a_current():
    assert len(RUNS) == 7, RUNS
    assert set(EXCLUDED) == {"rudder_step_10deg_in_current"}, EXCLUDED
    assert all(set(g3._trajectory(n)) >= {"ui_01", "ui_02", "ui_03"} for n in RUNS)


def test_the_vehicle_is_the_named_part_built_by_rpp():
    vehicle = _vehicle()
    assert type(vehicle).__name__ == "SpheroidAuv"
    # four fins (rudder, rudder, stern plane, stern plane) and the propulsor: each actuator names its command after
    # its plugin, so the four fins carry one name (a schedule gives that name the four values in actuator order)
    assert [d.name for d in vehicle.getInputDescriptions()] == ["fin_angle_command"] * 4 + ["drive_command"]
    assert [d.size for d in vehicle.getInputDescriptions()] == [1, 1, 1, 1, 1]
    description = json.loads(DESCRIPTION.read_text())
    part_id = description["Configurations"][CONFIGURATION]["Components"]["vessels"][0]["Id"]
    found = [p for p in (LIBRARY / ".rppws" / "parts").rglob("description.json")
             if json.loads(p.read_text()).get("Id") == part_id]
    assert len(found) == 1 and json.loads(found[0].read_text())["Name"] == PART


@pytest.mark.parametrize("name", RUNS)
def test_G3_through_rpp_trajectory_equals_matlab_rk4_run(name):
    table = g3._trajectory(name)
    states = _rk4_run(_vehicle(), table)
    error = np.abs(states - g3._columns(table, "x", 12)).max(axis=0)
    tolerance = g3._trajectory_tolerances()[name]
    over = [f"x{i + 1}: {error[i]:.3g} > {tolerance[i]:.3g}" for i in range(12) if error[i] > tolerance[i]]
    assert not over, f"{name}: " + "; ".join(over)


def _commands():
    return json.loads(DESCRIPTION.read_text())["Configurations"][CONFIGURATION]["Commands"]


def test_the_script_description_has_one_command_per_run_of_the_part():
    assert sorted(_commands()) == sorted(RUNS)


@pytest.mark.parametrize("name", RUNS)
def test_the_command_schedule_is_the_frozen_run_sample_by_sample(name):
    table = g3._trajectory(name)
    schedule = _commands()[name]
    vehicle = _vehicle()
    names = [d.name for d in vehicle.getInputDescriptions()]
    sizes = [d.size for d in vehicle.getInputDescriptions()]
    rows = _script().command_rows(schedule, names, sizes, table["t"])
    assert np.array_equal(rows, _frozen_inputs(vehicle, table))
    assert schedule["SampleTime"] == table["t"][1] - table["t"][0]
    assert schedule["Duration"] == table["t"][-1]
    start = g3._columns(table, "x", 12)[0]
    assert schedule["InitialSurgeSpeed"] == start[0] and not start[1:].any()


@pytest.mark.parametrize("name", RUNS)
def test_the_script_cvodes_run_agrees_with_the_rk4_run_within_the_integrator_bound(name):
    table = g3._trajectory(name)
    reference = g3._columns(table, "x", 12)
    rk4 = _rk4_run(_vehicle(), table)
    result = _simulation(name).run()[0]
    assert np.array_equal(result["time"], table["t"])
    cvodes = g3._ours_to_mss(result["states"][:, -12:])
    n = len(table["t"]) - 1
    bound = np.sqrt(12) * n * (CVODES_ABSTOL + CVODES_RELTOL * np.abs(reference).max(axis=0))
    error = np.abs(cvodes - rk4).max(axis=0)
    assert (error <= bound).all(), f"{name}: " + "; ".join(
        f"x{i + 1}: {error[i]:.3g} > {bound[i]:.3g}" for i in range(12) if error[i] > bound[i])


def test_the_integrator_bound_uses_the_tolerances_of_the_simulation_class():
    source = inspect.getsource(_script().Simulation._simulate)
    assert re.search(r'"abstol":\s*1e-8', source) and re.search(r'"reltol":\s*1e-6', source)
    assert (CVODES_ABSTOL, CVODES_RELTOL) == (1e-8, 1e-6)


def test_a_schedule_is_refused_when_it_does_not_match_the_vehicles_input_descriptions():
    command_rows = _script().command_rows
    names, sizes, times = ["fin_deflection_command", "shaft_speed_command"], [2, 1], np.array([0.0, 0.25])
    good = {"Steps": [{"Time": 0.0, "fin_deflection_command": [0.1, 0.0], "shaft_speed_command": [900.0]}]}
    assert command_rows(good, names, sizes, times).tolist() == [[0.1, 0.0, 900.0]] * 2
    missing = {"Steps": [{"Time": 0.0, "fin_deflection_command": [0.1, 0.0]}]}
    extra = {"Steps": [{**good["Steps"][0], "thruster_command": [1.0]}]}
    wrong_size = {"Steps": [{**good["Steps"][0], "fin_deflection_command": [0.1]}]}
    late = {"Steps": [{**good["Steps"][0], "Time": 1.0}]}
    for bad in (missing, extra, wrong_size, late):
        with pytest.raises(ValueError):
            command_rows(bad, names, sizes, times)
    # the vehicle's own order, not the schedule's: swapped names give the same row
    swapped = {"Steps": [{"shaft_speed_command": [900.0], "fin_deflection_command": [0.1, 0.0], "Time": 0.0}]}
    assert command_rows(swapped, names, sizes, times).tolist() == command_rows(good, names, sizes, times).tolist()
    assert command_rows(good, names[::-1], sizes[::-1], times).tolist() == [[900.0, 0.1, 0.0]] * 2


def test_inputs_of_one_name_take_that_names_values_one_after_the_other_in_actuator_order():
    """Four fins and a propulsor: ``fin_angle_command`` four times, ``drive_command`` once."""
    command_rows = _script().command_rows
    names, sizes, times = ["fin_angle_command"] * 4 + ["drive_command"], [1] * 5, np.array([0.0, 0.25, 0.5])
    schedule = {"Steps": [{"Time": 0.0, "fin_angle_command": [1.0, 2.0, 3.0, 4.0], "drive_command": [9.0]},
                          {"Time": 0.25, "fin_angle_command": [5.0, 6.0, 7.0, 8.0], "drive_command": [9.5]}]}
    assert command_rows(schedule, names, sizes, times).tolist() == [[1, 2, 3, 4, 9], [5, 6, 7, 8, 9.5],
                                                                    [5, 6, 7, 8, 9.5]]
    short = {"Steps": [{"Time": 0.0, "fin_angle_command": [1.0, 2.0], "drive_command": [9.0]}]}
    with pytest.raises(ValueError, match="fin_angle_command"):
        command_rows(short, names, sizes, times)
    # a vehicle with a two-value input and a one-value input of the same name takes three values
    mixed = {"Steps": [{"Time": 0.0, "u": [1.0, 2.0, 3.0]}]}
    assert command_rows(mixed, ["u", "u"], [2, 1], np.array([0.0])).tolist() == [[1.0, 2.0, 3.0]]


def test_a_command_the_configuration_does_not_have_is_refused():
    with pytest.raises(ValueError, match="no command"):
        _simulation("not_a_run")


# ---------------------------------------------------------------------------------------------------------------
# Otter (MSS): the commanded runs of the unmodified otter.m, through the two linked thrusters
# ---------------------------------------------------------------------------------------------------------------
OTTER_CONFIGURATION, OTTER_PART = "otter", "Otter (MSS)"
OTTER_DATA = Path(__file__).resolve().parent / "data" / "otter"
OTTER_TAG = "_mss_ed2ce7b"
OTTER_RUNS = sorted(f.stem[len("otter_trajectory_"):-len(OTTER_TAG)] for f in OTTER_DATA.glob(f"otter_trajectory_*{OTTER_TAG}.csv")
                    if "tolerances" not in f.name)
OTTER_COLUMNS = {"drive_command": ("ui_01", "ui_02")}  # port (n(1)), starboard (n(2)) shaft speeds, rad/s


def _otter_table(name):
    data = np.genfromtxt(OTTER_DATA / f"otter_trajectory_{name}{OTTER_TAG}.csv", delimiter=",", names=True)
    return {key: np.atleast_1d(data[key]) for key in data.dtype.names}


def _otter_tolerances():
    table = {}
    with open(OTTER_DATA / f"otter_trajectory_tolerances{OTTER_TAG}.csv", newline="") as handle:
        import csv
        for row in csv.DictReader(handle):
            table.setdefault(row["scenario"], np.zeros(12))[int(row["state"]) - 1] = float(row["tolerance"])
    return table


def _otter_inputs(vehicle, table, columns=OTTER_COLUMNS):
    """The frozen commands as the vehicle's stacked inputs: the two ``drive_command`` inputs take the name's two columns
    one after the other, in actuator order (port, starboard)."""
    descriptions = list(vehicle.getInputDescriptions())
    assert [d.name for d in descriptions] == ["drive_command", "drive_command"] and [d.size for d in descriptions] == [1, 1]
    return np.column_stack([table[c] for c in columns["drive_command"]])


def _otter_rk4_run(vehicle, table, columns=OTTER_COLUMNS, refine=1):
    """The frozen run stepped by the vehicle's own RK4 at the frozen step (MSS state order), commands held per row;
    ``refine`` > 1 takes that many times as many steps of the same method (a finer run of the same equations)."""
    inputs = _otter_inputs(vehicle, table, columns)
    steps = int(round((table["t"][1] - table["t"][0]) / g3.STEP)) * refine
    x = g3._mss_to_ours(np.array([table[f"x_{i:02d}"][0] for i in range(1, 13)]))
    states = np.zeros((len(table["t"]), 12))
    for row in range(len(table["t"])):
        states[row] = g3._ours_to_mss(x)
        if row < len(table["t"]) - 1:
            for _ in range(steps):
                x = vehicle.step_vector(x, inputs[row], g3.STEP / refine)
    return states


def _otter_reference(table):
    return np.column_stack([table[f"x_{i:02d}"] for i in range(1, 13)])


def test_the_otter_runs_are_the_six_commanded_runs():
    assert OTTER_RUNS == ["astern_from_speed", "commands_beyond_the_limits", "differential_turn", "spin_in_place",
                          "straight_from_rest_thrust_step", "zigzag_differential_steps"]
    assert sorted(_otter_tolerances()) == OTTER_RUNS


def test_the_otter_is_the_named_part_built_by_rpp_with_two_thrusters():
    vehicle = _vehicle(OTTER_CONFIGURATION)
    assert type(vehicle).__name__ == "Catamaran"
    assert [d.name for d in vehicle.getInputDescriptions()] == ["drive_command", "drive_command"]
    description = json.loads(DESCRIPTION.read_text())
    part_id = description["Configurations"][OTTER_CONFIGURATION]["Components"]["vessels"][0]["Id"]
    found = [p for p in (LIBRARY / ".rppws" / "parts").rglob("description.json")
             if json.loads(p.read_text()).get("Id") == part_id]
    assert len(found) == 1 and json.loads(found[0].read_text())["Name"] == OTTER_PART


@pytest.mark.parametrize("name", OTTER_RUNS)
def test_G3_through_rpp_otter_trajectory_equals_matlab_rk4_run(name):
    table = _otter_table(name)
    error = np.abs(_otter_rk4_run(_vehicle(OTTER_CONFIGURATION), table) - _otter_reference(table)).max(axis=0)
    tolerance = _otter_tolerances()[name]
    over = [f"x{i + 1}: {error[i]:.3g} > {tolerance[i]:.3g}" for i in range(12) if error[i] > tolerance[i]]
    assert not over, f"{name}: " + "; ".join(over)


@pytest.mark.parametrize("name", ["differential_turn", "zigzag_differential_steps"])
def test_G4_swapped_thrusters_or_the_corrected_trim_are_detected(name):
    """The gate can fail: the two shaft speeds swapped between port and starboard, and the corrected payload trim in
    place of MSS's (the part's form), each leave the frozen run by far more than its tolerance."""
    table = _otter_table(name)
    tolerance = _otter_tolerances()[name]
    swapped = _otter_rk4_run(_vehicle(OTTER_CONFIGURATION), table, {"drive_command": ("ui_02", "ui_01")})
    assert (np.abs(swapped - _otter_reference(table)).max(axis=0) > 10 * tolerance).any()
    p = vc.plugins()
    thrusters = [(p.Propulsor, vc.pick(p.Propulsor, vc.reference_part(n)))
                 for n in ("Otter port thruster", "Otter starboard thruster")]
    corrected = vc.otter(actuators=thrusters)  # the default form, lever_arm_moment
    run = _otter_rk4_run(corrected, table)
    assert (np.abs(run - _otter_reference(table)).max(axis=0) > 10 * tolerance).any()


@pytest.mark.parametrize("name", OTTER_RUNS)
def test_the_otter_command_schedule_is_the_frozen_run_sample_by_sample(name):
    table = _otter_table(name)
    schedule = json.loads(DESCRIPTION.read_text())["Configurations"][OTTER_CONFIGURATION]["Commands"][name]
    vehicle = _vehicle(OTTER_CONFIGURATION)
    names = [d.name for d in vehicle.getInputDescriptions()]
    sizes = [d.size for d in vehicle.getInputDescriptions()]
    rows = _script().command_rows(schedule, names, sizes, table["t"])
    assert np.array_equal(rows, _otter_inputs(vehicle, table))
    assert schedule["SampleTime"] == table["t"][1] - table["t"][0]
    assert schedule["Duration"] == table["t"][-1]
    assert schedule["InitialSurgeSpeed"] == table["x_01"][0] and not any(table[f"x_{i:02d}"][0] for i in range(2, 13))


def test_the_otter_has_one_command_per_frozen_run():
    commands = json.loads(DESCRIPTION.read_text())["Configurations"][OTTER_CONFIGURATION]["Commands"]
    assert sorted(commands) == OTTER_RUNS


REFINE = 4  # steps of the finer RK4 per frozen step: its own step error is (1/4)^4 of the frozen step's


@pytest.mark.parametrize("name", OTTER_RUNS)
def test_the_otter_script_cvodes_run_agrees_with_the_finer_rk4_run_within_the_integrator_bound(name):
    """The script's CVODES run of the schedule against the same equations stepped four times finer by RK4, within the
    integrator bound of the REMUS runs above (CVODES's own norm, ``sqrt(12) N (abstol + reltol max|x|)``).

    Why the finer run and not the frozen step: a first version of this test, written like the REMUS one, compared the
    CVODES run with the RK4 run at the frozen step (0.05 s) and failed on all six runs by factors 1.85-7.1 of the bound.
    The cause is the RK4 step error of the frozen step, not CVODES: the Otter's heave, roll and pitch modes
    (natural frequencies of several rad/s) make the frozen step's RK4 differ from a finer run of the same method by
    1.2e-4 to 2.0e-4 in the worst state of each run (the runs four and ten times finer differ by at most 1.5e-6, 2 % of
    the bound), while CVODES stays within 0.04-0.24 of its bound of the ten-times-finer run (measured once, as numbers
    for this note, not asserted). The REMUS runs (slower, Munk couplings removed) do not show it. The G3 gate above compares the frozen step with the frozen step, so MATLAB's and the
    library's RK4 carry the same step error and agree to 1e-14; this test is about CVODES alone."""
    table = _otter_table(name)
    reference = _otter_reference(table)
    fine = _otter_rk4_run(_vehicle(OTTER_CONFIGURATION), table, refine=REFINE)
    result = _simulation(name, configuration=OTTER_CONFIGURATION).run()[0]
    assert np.array_equal(result["time"], table["t"])
    cvodes = g3._ours_to_mss(result["states"][:, -12:])
    n = len(table["t"]) - 1
    bound = np.sqrt(12) * n * (CVODES_ABSTOL + CVODES_RELTOL * np.abs(reference).max(axis=0))
    error = np.abs(cvodes - fine).max(axis=0)
    assert (error <= bound).all(), f"{name}: " + "; ".join(
        f"x{i + 1}: {error[i]:.3g} > {bound[i]:.3g}" for i in range(12) if error[i] > bound[i])


@pytest.mark.parametrize("name", OTTER_RUNS)
def test_the_frozen_step_rk4_error_of_the_otter_runs_is_the_measured_size(name):
    """The frozen runs (MATLAB's and the library's alike) carry the RK4 step error of 0.05 s: the difference to the
    run four times finer is below 3e-4 in every state of every run (measured 1.2e-4 to 2.0e-4); a larger one would mean
    the equations changed, not the step."""
    table = _otter_table(name)
    vehicle = _vehicle(OTTER_CONFIGURATION)
    difference = np.abs(_otter_rk4_run(vehicle, table) - _otter_rk4_run(vehicle, table, refine=REFINE)).max()
    assert difference < 3e-4, difference


# ---------------------------------------------------------------------------------------------------------------
# Grethe: the commanded plant runs (the model file grethe.m @ 679a6e8 (retired 2026-10-30; file deleted, kept as tag archive/grethe-vehicle) with its defaults and its caller's actuator equations)
# ---------------------------------------------------------------------------------------------------------------
GRETHE_CONFIGURATION = "grethe"
GRETHE_DATA = Path(__file__).resolve().parent / "data" / "grethe"
# the runs the schedule can express: the state starts at rest or at a surge speed with the actuators at rest (the two
# free decays start with an attitude, which a schedule has no key for)
GRETHE_RUNS = ["coast_down_from_2p4", "steering_square_wave_5deg", "straight_from_rest_thrust_step", "turning_circle_20deg"]
GRETHE_COLUMNS = {"drive_command": ("ui_01",), "steering_command": ("ui_02",)}  # thrust command N, steering angle command rad


def _grethe_table(name):
    data = np.genfromtxt(GRETHE_DATA / f"grethe_m_plant_trajectory_{name}.csv", delimiter=",", names=True)
    return {key: np.atleast_1d(data[key]) for key in data.dtype.names}


def _grethe_inputs(vehicle, table):
    descriptions = list(vehicle.getInputDescriptions())
    assert [d.name for d in descriptions] == ["drive_command", "steering_command"]
    return np.column_stack([table[c] for d in descriptions for c in GRETHE_COLUMNS[d.name]])


def _grethe_rk4_run(vehicle, table, refine=1):
    """The frozen plant run by the vehicle's own RK4 (the commands held per row; ``refine`` times as many steps), states
    in the file's order ``x_01..12`` then the actual thrust and steering angle."""
    inputs = _grethe_inputs(vehicle, table)
    steps = int(round((table["t"][1] - table["t"][0]) / g3.STEP)) * refine
    x = g3._mss_to_ours(np.r_[[table[f"x_{i:02d}"][0] for i in range(1, 13)], table["act_01"][0], table["act_02"][0]],
                        vehicle)
    states = np.zeros((len(table["t"]), 14))
    for row in range(len(table["t"])):
        states[row] = g3._ours_to_mss(x, vehicle)
        if row < len(table["t"]) - 1:
            for _ in range(steps):
                x = vehicle.step_vector(x, inputs[row], g3.STEP / refine)
    return states


def _grethe_reference(table):
    return np.column_stack([table[f"x_{i:02d}"] for i in range(1, 13)] + [table["act_01"], table["act_02"]])


def test_the_grethe_has_one_command_per_run_the_schedule_can_express():
    commands = json.loads(DESCRIPTION.read_text())["Configurations"][GRETHE_CONFIGURATION]["Commands"]
    assert sorted(commands) == GRETHE_RUNS


@pytest.mark.parametrize("name", GRETHE_RUNS)
def test_the_grethe_command_schedule_is_the_frozen_plant_run_sample_by_sample(name):
    table = _grethe_table(name)
    schedule = json.loads(DESCRIPTION.read_text())["Configurations"][GRETHE_CONFIGURATION]["Commands"][name]
    vehicle = _vehicle(GRETHE_CONFIGURATION)
    names = [d.name for d in vehicle.getInputDescriptions()]
    sizes = [d.size for d in vehicle.getInputDescriptions()]
    rows = _script().command_rows(schedule, names, sizes, table["t"])
    assert np.array_equal(rows, _grethe_inputs(vehicle, table))
    assert schedule["SampleTime"] == table["t"][1] - table["t"][0] and schedule["Duration"] == table["t"][-1]
    assert schedule["InitialSurgeSpeed"] == table["x_01"][0]
    assert not any(table[f"x_{i:02d}"][0] for i in range(2, 13)) and not table["act_01"][0] and not table["act_02"][0]


@pytest.mark.parametrize("name", GRETHE_RUNS)
def test_the_grethe_script_cvodes_run_agrees_with_the_finer_rk4_run_within_the_integrator_bound(name):
    """The script's CVODES run of the schedule against the same plant stepped four times finer by RK4, within the
    integrator bound of the runs above (the Otter section explains why the finer run and not the frozen step)."""
    table = _grethe_table(name)
    reference = _grethe_reference(table)
    fine = _grethe_rk4_run(_vehicle(GRETHE_CONFIGURATION), table, refine=REFINE)
    result = _simulation(name, configuration=GRETHE_CONFIGURATION).run()[0]
    assert np.allclose(result["time"], table["t"], rtol=0.0, atol=1e-9)  # the file's times are k * 0.05 * 10 in floating point
    cvodes = g3._ours_to_mss(result["states"], _vehicle(GRETHE_CONFIGURATION))
    n = len(table["t"]) - 1
    bound = np.sqrt(14) * n * (CVODES_ABSTOL + CVODES_RELTOL * np.abs(reference).max(axis=0))
    error = np.abs(cvodes - fine).max(axis=0)
    assert (error <= bound).all(), f"{name}: " + "; ".join(
        f"state {i + 1}: {error[i]:.3g} > {bound[i]:.3g}" for i in range(14) if error[i] > bound[i])
