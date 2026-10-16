"""Run a named vehicle of the library's ``.rppws`` forward in time and print its
final state, in the form of Luka's ``scripts/vehicle_models/hull_vessel_simulation.py``:
a ``Simulation`` of ``more_simulation`` whose ``vessels`` slot takes any
``VehicleModel3D``, the vehicles being the named parts of ``.rppws/parts``
assigned to the configurations of
``.rppws/script_descriptions/vehicle_simulation.json``.

    python vehicle_simulation.py <workspace> --rpp-configuration remus100 --surge-speed 1.5 --duration 20
    python vehicle_simulation.py <workspace> --rpp-configuration remus100 --command zigzag_20_20 --out run.csv

Without ``--command`` every input is zero. ``--command <name>`` takes the
schedule ``Configurations.<configuration>.Commands.<name>`` of the script
description: a table of switching times, each row giving the value of every
input of the vehicle **by the name of its input description** (the vehicle's
own ``IODescription`` names, in its own order); a value holds from its time
until the next row (zero-order hold over each sample interval). ``Duration``,
``SampleTime`` and ``InitialSurgeSpeed`` of the schedule are the defaults of
``--duration``, ``--delta-t`` and ``--surge-speed``. A schedule that names an
input the vehicle does not have, leaves one out or gives it the wrong size is
refused.

The integrator (CVODES) needs the derivative of the vehicle's dynamics; the
hull loads of the REMUS 100 are guarded (``smooth_speed``) so it is finite at
rest and a run may start from rest; ``--surge-speed`` gives an initial surge
speed.

It needs rpp's registry with ``more_dynamics`` registered
(``rpp library register ./more_dynamics --link``).

Author:    Enio Krizman
Date:      2026-10-09
"""

import json
from argparse import ArgumentParser
from pathlib import Path

import casadi as ca
import numpy as np

from rpp_py.data_manager import DataManager

_DATA_MANAGER = DataManager()

from more_common.casadi_graph import RppCasadiGraph
from more_simulation.simulation import Simulation

STATE_NAMES = ("north", "east", "down", "roll", "pitch", "yaw", "u", "v", "w", "p", "q", "r")
DEFAULT_DELTA_T = 0.1  # s, Simulation's own default
DEFAULT_DURATION = 35.0  # s, Simulation's own default
SCHEDULE_KEYS = ("Description", "InitialSurgeSpeed", "Duration", "SampleTime", "Steps")


def load_commands(rpp_workspace, library, script, configuration):
    """``(configuration name, {command: schedule})`` of ``configuration`` (the active one when ``None``) in the
    script description of ``script`` of ``library`` in the workspace."""
    path = DataManager(workspace_path=rpp_workspace).get_script_description_path_from_library(library, script)
    description = json.loads(Path(path).read_text(encoding="utf-8"))
    name = configuration or description["ActiveConfiguration"]
    return name, description["Configurations"][name].get("Commands", {})


def command_rows(schedule, names, sizes, times):
    """The stacked input rows of ``schedule`` at ``times``: per row the value of the last step whose ``Time`` is
    not after the sample (to ``1e-9`` s), each input taken by its name, in the order of ``names``."""
    steps = schedule["Steps"]
    given = [{key for key in step if key != "Time"} for step in steps]
    for step, keys in zip(steps, given):
        if keys != set(names):
            raise ValueError(f"command step at t={step['Time']}: inputs {sorted(keys)} but the vehicle's input "
                             f"descriptions are {list(names)}")
        for name, size in zip(names, sizes):
            if len(step[name]) != size:
                raise ValueError(f"command step at t={step['Time']}: input {name!r} has {len(step[name])} values, "
                                 f"the vehicle's description has {size}")
    switch = np.array([step["Time"] for step in steps], dtype=float)
    if len(steps) == 0 or switch[0] > times[0] + 1e-9 or np.any(np.diff(switch) <= 0.0):
        raise ValueError("command steps must start at the first sample and have increasing times")
    rows = np.zeros((len(times), int(sum(sizes))))
    for k, time in enumerate(times):
        step = steps[int(np.searchsorted(switch, time + 1e-9, side="right")) - 1]
        rows[k, :] = np.concatenate([np.asarray(step[name], dtype=float) for name in names])
    return rows


class VehicleSimulation(Simulation):
    COMPONENTS = {
        "vessels": "List[more_dynamics::VehicleModel3D]",
    }
    RPP_SCRIPT_LIBRARY = "more_dynamics"
    RPP_SCRIPT_NAME = "vehicle_simulation"

    def __init__(
        self,
        rpp_workspace: str,
        rpp_configuration: str | None = None,
        delta_t: float | None = None,
        duration: float | None = None,
        surge_speed: float | None = None,
        command: str | None = None,
    ):
        self.command = command
        self.schedule = None
        if command is not None:
            configuration, commands = load_commands(
                rpp_workspace, self.RPP_SCRIPT_LIBRARY, self.RPP_SCRIPT_NAME, rpp_configuration)
            if command not in commands:
                raise ValueError(f"configuration {configuration!r} has no command {command!r} "
                                 f"(it has: {sorted(commands)})")
            self.schedule = commands[command]
            unknown = set(self.schedule) - set(SCHEDULE_KEYS)
            if unknown:
                raise ValueError(f"command {command!r}: unknown keys {sorted(unknown)}")
        schedule = self.schedule or {}
        self.surge_speed = schedule.get("InitialSurgeSpeed", 0.0) if surge_speed is None else surge_speed
        super().__init__(
            rpp_workspace=rpp_workspace,
            rpp_configuration=rpp_configuration,
            delta_t=schedule.get("SampleTime", DEFAULT_DELTA_T) if delta_t is None else delta_t,
            duration=schedule.get("Duration", DEFAULT_DURATION) if duration is None else duration,
        )

    def run(self):
        """Every vessel with the schedule as its inputs (zero inputs without one, as ``Simulation.run``)."""
        if self.schedule is None:
            return super().run()
        num_steps = int(self.duration / self.delta_t)
        times = np.arange(num_steps + 1, dtype=float) * self.delta_t
        self.results = []
        for vessel_index, vessel in enumerate(self.vessels):
            graph = RppCasadiGraph(vessel.graph())
            descriptions = vessel.getInputDescriptions()
            rows = command_rows(self.schedule, [d.name for d in descriptions], [d.size for d in descriptions], times)
            result = self._simulate(graph, ca.DM(rows), num_steps, vessel_index)
            result["input_names"] = [d.name for d in descriptions]
            result["input_sizes"] = [d.size for d in descriptions]
            self.results.append(result)
        return self.results

    def _extract_initial_conditions(self, graph):
        initial = np.array(super()._extract_initial_conditions(graph), dtype=float).ravel()
        initial[-6] = self.surge_speed  # u, the first of the six body velocities that close the state
        return ca.DM(initial)


def write_csv(path, result):
    """``t``, the twelve vehicle states and the vehicle's inputs (named by its input descriptions) per sample."""
    names = list(STATE_NAMES)
    for name, size in zip(result.get("input_names", []), result.get("input_sizes", [])):
        names += [name] if size == 1 else [f"{name}_{i + 1}" for i in range(size)]
    table = np.column_stack([result["time"], result["states"][:, -12:], result["inputs"]])
    np.savetxt(path, table, delimiter=",", header="t," + ",".join(names), comments="")


def main() -> None:
    parser = ArgumentParser()
    parser.add_argument("rpp_workspace")
    parser.add_argument("--rpp-configuration", default=None)
    parser.add_argument("--command", default=None, help="a schedule of Configurations.<configuration>.Commands")
    parser.add_argument("--duration", type=float, default=None, help="simulated time, s (default: the command's)")
    parser.add_argument("--delta-t", type=float, default=None, help="step between samples, s (default: the command's)")
    parser.add_argument("--surge-speed", type=float, default=None, help="initial surge speed, m/s")
    parser.add_argument("--out", default=None, help="write the trajectory to this CSV (one file per vessel)")
    arguments = parser.parse_args()
    simulation = VehicleSimulation(
        rpp_workspace=arguments.rpp_workspace,
        rpp_configuration=arguments.rpp_configuration,
        delta_t=arguments.delta_t,
        duration=arguments.duration,
        surge_speed=arguments.surge_speed,
        command=arguments.command,
    )
    results = simulation.run()
    for index, result in enumerate(results):
        final = result["states"][-1][-12:]
        print("t      " + "  ".join(f"{name:>8}" for name in STATE_NAMES))
        print(f"{result['time'][-1]:5.1f}  " + "  ".join(f"{value:8.3f}" for value in final))
        if arguments.out:
            out = Path(arguments.out)
            write_csv(out if len(results) == 1 else out.with_name(f"{out.stem}_{index + 1}{out.suffix}"), result)


if __name__ == "__main__":
    main()
