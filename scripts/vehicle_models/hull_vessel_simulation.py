from pathlib import Path

from rpp_py.data_manager import DataManager

_DATA_MANAGER = DataManager()

from more_simulation.plotting import plot_results
from more_simulation.simulation import Simulation


class HullVesselSimulation(Simulation):
    COMPONENTS = {
        "vessels": "List[more_dynamics::VehicleModel3D]",
    }

    def __init__(self, delta_t: float = 0.1, duration: float = 35.0):
        super().__init__(
            delta_t=delta_t,
            duration=duration,
            script_path=Path(__file__),
        )


def main() -> None:
    simulation = HullVesselSimulation()
    results = simulation.run()
    plot_results(results)


if __name__ == "__main__":
    main()
