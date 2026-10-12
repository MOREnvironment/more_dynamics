import casadi as ca
import numpy as np
import pytest

from rpp_py.data_manager import DataManager

_DATA_MANAGER = DataManager()

from more_common.casadi_graph import RppCasadiGraph  # noqa: E402
from more_dynamics.plugins.force_producers.jet_nozzle import (  # noqa: E402
    JetNozzle,
)
from more_dynamics.plugins.hydrodynamics import (  # noqa: E402
    linear_surface_hydrodynamics,
)
from more_dynamics.plugins.hydrostatics import (  # noqa: E402
    linear_surface_hydrostatics,
)
from more_dynamics.plugins.vehicle_models.hull_vessel import (  # noqa: E402
    HullVessel,
)
from rpp_plugin_types.more_dynamics import VehicleModel3D  # noqa: E402
from rpp_py.context import ComponentContext  # noqa: E402
from rpp_py.parameter_handler import ParameterHandler  # noqa: E402


def _context(plugin, parameters=None, **kwargs) -> ComponentContext:
    resolved = ParameterHandler.resolve_params(plugin.PARAMETERS, {})
    resolved.params.update(parameters or {})
    return ComponentContext(instance=plugin, params=resolved, **kwargs)


@pytest.fixture(scope="module")
def vessel() -> HullVessel:
    hull = HullVessel()
    _context(
        hull,
        {
            "mass": 200.0,
            "inertia": [55.0, 200.0, 200.0],
            "added_mass": [11.3, 300.0, 200.0, 11.0, 160.0, 340.0],
            "center_of_gravity": [0.0, 0.0, 0.1],
        },
        subcomponents={
            "actuators": [
                _context(
                    JetNozzle(),
                    {
                        "location": [-1.6, 0.0, 0.0],
                        "max_thrust": 6000.0,
                        "min_angle": -0.5235987756,
                        "max_angle": 0.5235987756,
                    },
                )
            ],
            "sensors": [],
            "hydrostatics": [
                _context(
                    linear_surface_hydrostatics.LinearSurfaceHydrostatics(),
                    {
                        "length": 3.5,
                        "beam": 1.2,
                        "draft": 0.16,
                        "block_coefficient": 0.29,
                        "waterplane_coefficient": 0.7,
                        "center_of_gravity": [0.0, 0.0, 0.1],
                        "longitudinal_center_of_flotation": -0.15,
                    },
                )
            ],
            "hydrodynamics": [
                _context(
                    linear_surface_hydrodynamics.LinearSurfaceHydrodynamics(),
                    {
                        "damping_coefficients": [
                            240.0, 1800.0, 7000.0, 900.0, 6500.0, 8000.0,
                        ]
                    },
                )
            ],
        },
        spec=HullVessel.COMPONENTS,
    ).initialize()
    return hull


def _command(*values: float) -> list:
    command = VehicleModel3D.Command()
    command.data.extend(values)
    return [command]


def _state(odometry) -> np.ndarray:
    return HullVessel._state_from_odometry(odometry)


def _drive(vessel, command, steps, delta_t=0.05):
    odometry = HullVessel._odometry_from_state(np.zeros(12))
    for step in range(steps):
        odometry = vessel.step(odometry, command, step * delta_t, delta_t)
    return _state(odometry)


def test_odometry_conversion_round_trips_the_vessel_state():
    state = np.array(
        [1.0, -2.0, 0.5, 0.1, -0.2, 2.5, 3.0, 0.4, -0.1, 0.01, 0.02, -0.3]
    )

    np.testing.assert_allclose(
        _state(HullVessel._odometry_from_state(state)), state, atol=1e-12
    )


def test_step_without_command_keeps_a_resting_vessel_at_rest(vessel):
    np.testing.assert_allclose(
        _drive(vessel, _command(0.0, 0.0), 40), np.zeros(12), atol=1e-9
    )


def test_step_reaches_the_cruise_speed_of_the_thrust(vessel):
    state = _drive(vessel, _command(0.2, 0.0), 400)

    assert state[6] == pytest.approx(0.2 * 6000.0 / 240.0, rel=1e-3)
    assert state[0] > 50.0
    np.testing.assert_allclose(state[[1, 5, 7, 11]], 0.0, atol=1e-6)


def test_step_turns_with_the_nozzle(vessel):
    straight = _drive(vessel, _command(0.2, 0.0), 200)
    turning = _drive(vessel, _command(0.2, 0.5), 200)

    assert abs(turning[5]) > 0.3
    assert abs(turning[1]) > 1.0
    assert abs(straight[5]) < 1e-6


def test_step_matches_the_graph_once_the_actuator_has_settled(vessel):
    """With the nozzle at its steady state the graph is the same model."""
    graph = RppCasadiGraph(vessel.graph())
    command = np.array([0.2, 0.5])
    delta_t = 0.05
    state_symbol = graph.step.sx_in(0)
    input_symbol = graph.step.sx_in(1)
    integrator = ca.integrator(
        "reference",
        "cvodes",
        {
            "x": state_symbol,
            "p": input_symbol,
            "ode": graph.step(state_symbol, input_symbol),
        },
        0,
        delta_t,
        {"abstol": 1e-10, "reltol": 1e-10},
    )
    state = ca.DM.zeros(graph.num_states)
    state[0] = 0.2
    state[1] = 0.5 * 0.5235987756
    for _ in range(200):
        state = integrator(x0=state, p=command)["xf"]

    np.testing.assert_allclose(
        _drive(vessel, _command(*command), 200, delta_t),
        np.asarray(state.full()).reshape(-1)[-12:],
        atol=1e-4,
    )


def test_step_ignores_missing_command_values(vessel):
    np.testing.assert_allclose(
        _drive(vessel, _command(0.2), 100),
        _drive(vessel, _command(0.2, 0.0), 100),
    )


def test_step_requires_initialization():
    with pytest.raises(RuntimeError, match="must be initialized"):
        HullVessel().step(
            HullVessel._odometry_from_state(np.zeros(12)), [], 0.0, 0.1
        )
