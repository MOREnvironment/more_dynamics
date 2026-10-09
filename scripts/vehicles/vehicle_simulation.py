"""Run a named vehicle of the library's ``.rppws`` forward in time with zero
commands and print its final state, in the form of Luka's
``scripts/vehicle_models/hull_vessel_simulation.py``: a ``Simulation`` of
``more_simulation`` whose ``vessels`` slot takes any ``VehicleModel3D``, the
vehicles being the parts REMUS 100 and Otter of ``.rppws/parts`` assigned to
the configurations of ``.rppws/script_descriptions/vehicle_simulation.json``.

    python vehicle_simulation.py <workspace> --rpp-configuration remus100 --surge-speed 1.5 --duration 20

The integrator (CVODES) needs the derivative of the vehicle's dynamics: the
slender-body lift of the REMUS hull has none at exactly zero relative
velocity, so a REMUS run starts with a surge speed.

It needs rpp's registry with ``more_dynamics`` registered
(``rpp library register ./more_dynamics --link``).

Author:    Enio Krizman
Date:      2026-10-09
"""

from argparse import ArgumentParser

import casadi as ca
import numpy as np

from rpp_py.data_manager import DataManager

_DATA_MANAGER = DataManager()

from more_simulation.simulation import Simulation

STATE_NAMES = ("north", "east", "down", "roll", "pitch", "yaw", "u", "v", "w", "p", "q", "r")


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
        delta_t: float = 0.1,
        duration: float = 35.0,
        surge_speed: float = 0.0,
    ):
        self.surge_speed = surge_speed
        super().__init__(
            rpp_workspace=rpp_workspace,
            rpp_configuration=rpp_configuration,
            delta_t=delta_t,
            duration=duration,
        )

    def _extract_initial_conditions(self, graph):
        initial = np.array(super()._extract_initial_conditions(graph), dtype=float).ravel()
        initial[-6] = self.surge_speed  # u, the first of the six body velocities that close the state
        return ca.DM(initial)


def main() -> None:
    parser = ArgumentParser()
    parser.add_argument("rpp_workspace")
    parser.add_argument("--rpp-configuration", default=None)
    parser.add_argument("--duration", type=float, default=35.0, help="simulated time, s")
    parser.add_argument("--delta-t", type=float, default=0.1, help="step between samples, s")
    parser.add_argument("--surge-speed", type=float, default=0.0, help="initial surge speed, m/s")
    arguments = parser.parse_args()
    simulation = VehicleSimulation(
        rpp_workspace=arguments.rpp_workspace,
        rpp_configuration=arguments.rpp_configuration,
        delta_t=arguments.delta_t,
        duration=arguments.duration,
        surge_speed=arguments.surge_speed,
    )
    for result in simulation.run():
        final = result["states"][-1][-12:]
        print("t      " + "  ".join(f"{name:>8}" for name in STATE_NAMES))
        print(f"{result['time'][-1]:5.1f}  " + "  ".join(f"{value:8.3f}" for value in final))


if __name__ == "__main__":
    main()
