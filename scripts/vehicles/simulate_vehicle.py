"""Run a ``MarineCraft6DOF`` vehicle of ``scripts/vehicles/.rppws`` forward in
time with zero commands, from a given surge speed, and print its state.

    python simulate_vehicle.py remus100 --surge-speed 1.5 --duration 20

The vehicle is built by rpp's own ``ComponentContextBuilder`` from the script
description of the same name (``.rppws/script_descriptions/<vehicle>.json``)
and advanced with its own ``step`` (RK4 sub-steps, zero-order hold). It needs
rpp's registry with ``more_dynamics`` registered
(``rpp library register ./more_dynamics --link``).

Author:    Enio Krizman
Date:      2026-10-08
"""

import math
from argparse import ArgumentParser
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
STATE_NAMES = ("north", "east", "down", "roll", "pitch", "yaw", "u", "v", "w", "p", "q", "r")


def build(vehicle):
    from rpp_py.data_manager import DataManager
    from rpp_plugin_registrator.library_manager import LibraryManager
    data_manager = DataManager(library_manager=LibraryManager())
    from rpp_py.context_builder import ComponentContextBuilder
    description = HERE / ".rppws" / "script_descriptions" / f"{vehicle}.json"
    context = ComponentContextBuilder(data_manager=data_manager).build_script_from_description_path(
        str(description), configuration=vehicle)
    context.initialize()
    return context.get_component("vessels")[0]


def main():
    parser = ArgumentParser()
    parser.add_argument("vehicle", choices=("remus100", "otter"))
    parser.add_argument("--surge-speed", type=float, default=0.0, help="initial surge speed, m/s")
    parser.add_argument("--duration", type=float, default=10.0, help="simulated time, s")
    parser.add_argument("--print-every", type=float, default=1.0, help="s between printed states")
    arguments = parser.parse_args()
    vessel = build(arguments.vehicle)
    commands = np.zeros(sum(d.size for d in vessel.getInputDescriptions()))
    x = np.zeros(12)
    x[6] = arguments.surge_speed
    print("t      " + "  ".join(f"{name:>8}" for name in STATE_NAMES))
    steps = int(round(arguments.duration / arguments.print_every))
    for k in range(steps + 1):
        print(f"{k * arguments.print_every:5.1f}  " + "  ".join(f"{v:8.3f}" for v in x))
        if k < steps:
            x = vessel.step_vector(x, commands, arguments.print_every)
    print(f"inputs held at zero: {[d.name for d in vessel.getInputDescriptions()]}")


if __name__ == "__main__":
    main()
