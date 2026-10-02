from argparse import ArgumentParser

from rpp_py.data_manager import DataManager

_DATA_MANAGER = DataManager()

from more_simulation.plotting import plot_results
from more_simulation.simulation import Simulation


class HullVesselSimulation(Simulation):
    COMPONENTS = {
        "vessels": "List[more_dynamics::VehicleModel3D]",
    }
    RPP_SCRIPT_LIBRARY = "more_dynamics"
    RPP_SCRIPT_NAME = "hull_vessel_simulation"

    def __init__(
        self,
        rpp_workspace: str,
        rpp_configuration: str | None = None,
        delta_t: float = 0.1,
        duration: float = 35.0,
    ):
        super().__init__(
            rpp_workspace=rpp_workspace,
            rpp_configuration=rpp_configuration,
            delta_t=delta_t,
            duration=duration,
        )


def main() -> None:
    parser = ArgumentParser()
    parser.add_argument("rpp_workspace")
    parser.add_argument("--rpp-configuration", default=None)
    arguments = parser.parse_args()
    simulation = HullVesselSimulation(
        rpp_workspace=arguments.rpp_workspace,
        rpp_configuration=arguments.rpp_configuration,
    )
    results = simulation.run()
    plot_results(results)


if __name__ == "__main__":
    main()
