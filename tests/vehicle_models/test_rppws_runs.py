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
vehicle (fins and propeller), the added-mass Coriolis without the Munk
couplings, still water. The runs of the bare hull and of the Munk-kept form are
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
# the frozen file's command columns (MSS order: rudder, stern plane, shaft speed) as the vehicle's named inputs
FROZEN_COMMAND_COLUMNS = {"fin_deflection_command": ("ui_01", "ui_02"), "shaft_speed_command": ("ui_03",)}
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


def _simulation(command=None, **options):
    """A ``VehicleSimulation`` built by rpp's builder from the named part (``command`` a schedule name)."""
    module = _script()
    try:
        return module.VehicleSimulation(rpp_workspace=str(LIBRARY), rpp_configuration=CONFIGURATION,
                                        command=command, **options)
    except ValueError:
        raise
    except Exception as exc:  # more_dynamics is not registered here
        pytest.skip(f"{vc.REGISTER_MESSAGE} ({type(exc).__name__}: {str(exc)[:80]})")


def _vehicle():
    if "vehicle" not in _CACHE:
        _CACHE["vehicle"] = _simulation().vessels[0]
    return _CACHE["vehicle"]


def _frozen_inputs(vehicle, table):
    """The frozen file's commands as the vehicle's stacked inputs, assembled in the order of its input
    descriptions by name (a name the table has no column for, or a column no input has, is an error)."""
    names = [d.name for d in vehicle.getInputDescriptions()]
    assert sorted(names) == sorted(FROZEN_COMMAND_COLUMNS), names
    sizes = {d.name: d.size for d in vehicle.getInputDescriptions()}
    assert all(len(FROZEN_COMMAND_COLUMNS[n]) == sizes[n] for n in names)
    return np.column_stack([table[c] for n in names for c in FROZEN_COMMAND_COLUMNS[n]])


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
    assert [d.name for d in vehicle.getInputDescriptions()] == ["fin_deflection_command", "shaft_speed_command"]
    assert [d.size for d in vehicle.getInputDescriptions()] == [2, 1]
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


def test_a_command_the_configuration_does_not_have_is_refused():
    with pytest.raises(ValueError, match="no command"):
        _simulation("not_a_run")


@pytest.mark.parametrize("configuration", ["otter", "grethe"])
def test_the_parts_without_actuators_run_through_the_script_with_zero_input(configuration):
    """Otter (MSS) and Grethe have no force producer child (their parts carry none): the vehicle takes no input,
    the script runs it at zero input, and and a command is refused (the configuration has no schedule)."""
    module = _script()
    simulation = module.VehicleSimulation(rpp_workspace=str(LIBRARY), rpp_configuration=configuration,
                                          duration=5.0, delta_t=0.5)
    assert list(simulation.vessels[0].getInputDescriptions()) == []
    result = simulation.run()[0]
    assert result["inputs"].shape == (11, 0) and np.isfinite(result["states"]).all()
