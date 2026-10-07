"""Drive a jetski-like HullVessel and measure it with noisy MORE sensors.

The vessel is assembled from the registered plugins the same way a workspace
does it, and every sensor is sampled through its payload noise description.
"""
import os
from pathlib import Path

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
from more_sensors.models import SensorSampler  # noqa: E402
from more_sensors.plugins.dvl import DVL  # noqa: E402
from more_sensors.plugins.gnss import GNSS  # noqa: E402
from more_sensors.plugins.imu import IMU  # noqa: E402
from more_sensors.plugins.magnetometer import (  # noqa: E402
    Magnetometer,
)
from more_sensors.plugins.pose_sensor import PoseSensor  # noqa: E402
from more_sensors.plugins.pressure import Pressure  # noqa: E402
from more_sensors.plugins.sbl import SBL  # noqa: E402
from rpp_py.context import ComponentContext  # noqa: E402
from rpp_py.parameter_handler import ParameterHandler  # noqa: E402


# Run with MORE_STORE_PLOTS=1 to write the sensor plots to the data folder.
STORE_DATA = os.environ.get("MORE_STORE_PLOTS", "") == "1"
DELTA_T = 0.1
DURATION = 40.0

SENSOR_PARAMETERS = {
    PoseSensor: {},
    IMU: {
        "noise_enabled": True,
        "random_seed": 42,
        "rate_hz": 10.0,
        "gyro_bias": [0.002, -0.001, 0.003],
        "gyro_white_noise_std_per_sample": [0.001, 0.001, 0.001],
        "gyro_bias_random_walk_std": [1e-4, 1e-4, 1e-4],
        "orientation_white_noise_std_per_sample": [0.002, 0.002, 0.005],
        "location": [0.3, 0.0, 0.2],
        "accel_bias": [0.05, -0.02, 0.03],
        "accel_white_noise_std_per_sample": [0.02, 0.02, 0.02],
    },
    Magnetometer: {
        "noise_enabled": True,
        "random_seed": 47,
        "white_noise_std_per_sample": [0.3e-6, 0.3e-6, 0.3e-6],
    },
    GNSS: {
        "noise_enabled": True,
        "random_seed": 44,
        "rate_hz": 1.0,
        "location": [0.5, 0.0, 1.2],
        "fix_location": [45.8007257, 15.9721655, 0.0],
        "position_white_noise_std_per_sample": [0.5, 0.5, 1.0],
    },
    DVL: {
        "noise_enabled": True,
        "random_seed": 43,
        "rate_hz": 5.0,
        "location": [1.0, 0.0, -0.3],
        "scale_factor": [1.01, 1.0, 1.0],
        "bias": [0.02, 0.0, 0.0],
        "white_noise_std_per_sample": [0.01, 0.01, 0.02],
    },
    Pressure: {
        "noise_enabled": True,
        "random_seed": 45,
        "location": [0.0, 0.0, -0.5],
        "white_noise_std_per_sample": [50.0],
    },
    SBL: {
        "noise_enabled": True,
        "random_seed": 46,
        "rate_hz": 2.0,
        "dropout_probability": 0.2,
        "array_position": [-20.0, 0.0, 0.0],
        "white_noise_std_per_sample": [0.1, 0.1, 0.2],
    },
}


def _context(plugin, parameters=None, **kwargs) -> ComponentContext:
    resolved = ParameterHandler.resolve_params(plugin.PARAMETERS, {})
    unknown = set(parameters or {}) - set(resolved.params)
    assert not unknown, f"undeclared parameters: {unknown}"
    resolved.params.update(parameters or {})
    return ComponentContext(instance=plugin, params=resolved, **kwargs)


def _jetski() -> HullVessel:
    """Assemble the Jetski of the more_simulation workspace with sensors."""
    jet = _context(
        JetNozzle(),
        {
            "location": [-1.6, 0.0, 0.0],
            "max_thrust": 6000.0,
            "thrust_coefficient": 1.0,
            "thrust_rise_time": 0.7,
            "thrust_fall_time": 0.25,
            "nozzle_velocity_radpersec": 1.5,
            "min_angle": -0.5235987756,
            "max_angle": 0.5235987756,
        },
    )
    hydrodynamics = _context(
        linear_surface_hydrodynamics.LinearSurfaceHydrodynamics(),
        {
            "damping_coefficients": [
                240.0,
                1800.0,
                7000.0,
                900.0,
                6500.0,
                8000.0,
            ]
        },
    )
    hydrostatics = _context(
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
    sensors = [
        _context(plugin_type(), parameters)
        for plugin_type, parameters in SENSOR_PARAMETERS.items()
    ]

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
            "actuators": [jet],
            "sensors": sensors,
            "hydrostatics": [hydrostatics],
            "hydrodynamics": [hydrodynamics],
        },
        spec=HullVessel.COMPONENTS,
    ).initialize()
    return hull


def _inputs(num_steps: int) -> np.ndarray:
    """Throttle up, cruise straight, carve a turn, then throttle down.

    The throttle stays low: this hull model capsizes above half throttle.
    """
    time = np.arange(num_steps) * DELTA_T
    thrust = np.zeros(num_steps)
    angle = np.zeros(num_steps)
    thrust[(time >= 1.0) & (time < 25.0)] = 0.2
    thrust[(time >= 25.0) & (time < 35.0)] = 0.12
    angle[(time >= 10.0) & (time < 20.0)] = 0.4
    return np.column_stack([thrust, angle])


def _simulate(
    graph: RppCasadiGraph,
    inputs: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Integrate the vessel graph into state and acceleration histories."""
    state_symbol = graph.step.sx_in(0)
    input_symbol = graph.step.sx_in(1)
    integrator = ca.integrator(
        "jetski",
        "cvodes",
        {
            "x": state_symbol,
            "p": input_symbol,
            "ode": graph.step(state_symbol, input_symbol),
        },
        0,
        DELTA_T,
        {"abstol": 1e-8, "reltol": 1e-6},
    )
    state = ca.DM.zeros(graph.num_states)
    vessel_states = []
    vessel_accelerations = []
    for control_input in [*inputs, inputs[-1]]:
        vessel_states.append(np.asarray(state.full()).reshape(-1)[-12:])
        vessel_accelerations.append(
            np.asarray(graph.step(state, control_input).full()).reshape(-1)[
                -6:
            ]
        )
        state = integrator(x0=state, p=control_input)["xf"]
    return np.asarray(vessel_states), np.asarray(vessel_accelerations)


@pytest.fixture(scope="module")
def drive():
    hull = _jetski()
    graph = RppCasadiGraph(hull.graph())
    num_steps = int(DURATION / DELTA_T)
    vessel_states, vessel_accelerations = _simulate(
        graph, _inputs(num_steps)
    )
    time = np.arange(num_steps + 1) * DELTA_T

    measurements = {}
    for sensor in hull.sensors:
        sampler = SensorSampler(sensor.graph())
        samples = [
            (
                sampler.sample(state, instant, acceleration),
                sampler.truth(state, acceleration),
                instant,
            )
            for state, acceleration, instant in zip(
                vessel_states, vessel_accelerations, time
            )
        ]
        measurements[type(sensor).__name__] = [
            sample for sample in samples if sample[0] is not None
        ]
    return hull, graph, time, vessel_states, measurements


def _errors(samples) -> np.ndarray:
    return np.array([measured - truth for measured, truth, _ in samples])


def test_hull_accepts_every_sensor(drive):
    hull, graph, _, _, _ = drive

    assert len(hull.sensors) == len(SENSOR_PARAMETERS)
    assert graph.num_outputs == 12 + 6 + 9 + 3 + 3 + 3 + 1 + 3


def test_jetski_accelerates_and_turns(drive):
    _, _, _, vessel_states, _ = drive

    assert np.all(np.isfinite(vessel_states))
    assert 4.0 < vessel_states[:, 6].max() < 6.0
    assert 0.5 < abs(vessel_states[-1, 5]) < 1.5
    assert np.hypot(*vessel_states[-1, 0:2]) > 100.0
    assert np.abs(vessel_states[:, 3:5]).max() < 0.05


def test_sensors_report_at_their_configured_rate(drive):
    _, _, time, _, measurements = drive

    assert len(measurements["PoseSensor"]) == len(time)
    assert len(measurements["Pressure"]) == len(time)
    assert len(measurements["IMU"]) == len(time)
    assert len(measurements["DVL"]) == pytest.approx(DURATION * 5.0, abs=1)
    assert len(measurements["GNSS"]) == pytest.approx(DURATION * 1.0, abs=1)


def test_sbl_drops_out_with_the_configured_probability(drive):
    _, _, _, _, measurements = drive

    due = DURATION * 2.0 + 1
    assert len(measurements["SBL"]) / due == pytest.approx(0.8, abs=0.12)


def test_pose_sensor_without_noise_reports_the_truth(drive):
    _, _, _, vessel_states, measurements = drive

    measured = np.array([sample[0] for sample in measurements["PoseSensor"]])
    np.testing.assert_allclose(measured, vessel_states[:, 0:6])


def test_sensor_noise_matches_its_configuration(drive):
    _, _, _, _, measurements = drive

    gyro_error = _errors(measurements["IMU"])[:, 3:6]
    np.testing.assert_allclose(
        gyro_error.mean(axis=0), [0.002, -0.001, 0.003], atol=6e-4
    )
    np.testing.assert_allclose(gyro_error.std(axis=0), 0.001, rtol=0.25)

    accel_error = _errors(measurements["IMU"])[:, 6:9]
    np.testing.assert_allclose(
        accel_error.mean(axis=0), [0.05, -0.02, 0.03], atol=6e-3
    )

    pressure_error = _errors(measurements["Pressure"])
    np.testing.assert_allclose(pressure_error.std(), 50.0, rtol=0.15)

    truth = np.array([sample[1] for sample in measurements["DVL"]])
    dvl_error = _errors(measurements["DVL"]) - 0.01 * truth * [1, 0, 0]
    np.testing.assert_allclose(
        dvl_error.mean(axis=0), [0.02, 0.0, 0.0], atol=5e-3
    )


def test_gnss_error_is_metres_at_the_fix_location(drive):
    _, _, _, _, measurements = drive

    error = _errors(measurements["GNSS"])
    north_error = error[:, 0] * 111_000.0
    altitude_error = error[:, 2]

    assert 0.2 < north_error.std() < 1.0
    assert 0.4 < altitude_error.std() < 2.0


def test_imu_measures_gravity_and_the_turn(drive):
    _, _, time, _, measurements = drive

    truth = np.array([sample[1] for sample in measurements["IMU"]])
    cruising = (time > 5.0) & (time < 9.0)
    turning = (time > 15.0) & (time < 19.0)

    np.testing.assert_allclose(truth[cruising, 8].mean(), 9.80665, atol=0.2)
    assert abs(truth[turning, 7].mean()) > 0.2
    assert abs(truth[cruising, 7].mean()) < 0.05


def test_magnetometer_follows_the_heading(drive):
    _, _, _, vessel_states, measurements = drive

    truth = np.array([sample[1] for sample in measurements["Magnetometer"]])

    np.testing.assert_allclose(
        np.arctan2(truth[:, 0], truth[:, 1]),
        vessel_states[:, 5],
        atol=0.01,
    )


def _plot(name, title, samples, labels, unit, scale=1.0):
    import matplotlib.pyplot as plt

    instants = [instant for _, _, instant in samples]
    measured = np.array([sample[0] for sample in samples]) * scale
    truth = np.array([sample[1] for sample in samples]) * scale
    figure, axes = plt.subplots(
        len(labels), 1, sharex=True, figsize=(9, 2.4 * len(labels) + 1)
    )
    for axis, column, label in zip(
        np.atleast_1d(axes), labels, labels.values()
    ):
        axis.plot(instants, measured[:, column], ".", ms=3, label="measured")
        axis.plot(instants, truth[:, column], lw=1, label="truth")
        axis.set_ylabel(f"{label} ({unit})")
        axis.grid(True, alpha=0.3)
    np.atleast_1d(axes)[0].legend(loc="upper right")
    np.atleast_1d(axes)[0].set_title(title)
    np.atleast_1d(axes)[-1].set_xlabel("time (s)")
    figure.tight_layout()
    figure.savefig(Path(__file__).parent / "data" / name)
    plt.close(figure)


def test_store_drive_plots(drive):
    if not STORE_DATA:
        pytest.skip("set MORE_STORE_PLOTS=1 to write the plots")
    import matplotlib.pyplot as plt

    _, _, _, vessel_states, measurements = drive
    data_folder = Path(__file__).parent / "data"
    data_folder.mkdir(parents=True, exist_ok=True)

    figure = plt.figure(figsize=(7, 6))
    plt.plot(vessel_states[:, 0], vessel_states[:, 1], label="truth")
    sbl = np.array([sample[0] for sample in measurements["SBL"]])
    plt.plot(sbl[:, 0] - 20.0, sbl[:, 1], ".", ms=4, label="SBL fix")
    plt.xlabel("x / east (m)")
    plt.ylabel("y / north (m)")
    plt.axis("equal")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.title("Jetski trajectory")
    figure.savefig(data_folder / "jetski_trajectory.png")
    plt.close(figure)

    axes = {0: "x", 1: "y", 2: "z"}
    _plot(
        "jetski_imu_gyro.png", "IMU angular velocity",
        measurements["IMU"], {3: "x", 4: "y", 5: "z"}, "rad/s",
    )
    _plot(
        "jetski_imu_accel.png", "IMU specific force",
        measurements["IMU"], {6: "x", 7: "y", 8: "z"}, "m/s^2",
    )
    _plot(
        "jetski_magnetometer.png", "Magnetometer",
        measurements["Magnetometer"], axes, "uT", scale=1e6,
    )
    _plot(
        "jetski_dvl.png", "DVL body velocity",
        measurements["DVL"], {0: "surge", 1: "sway", 2: "heave"}, "m/s",
    )
    _plot(
        "jetski_gnss.png", "GNSS fix",
        measurements["GNSS"], {0: "latitude", 1: "longitude"}, "deg",
    )
    _plot(
        "jetski_pressure.png", "Pressure",
        measurements["Pressure"], {0: "pressure"}, "Pa",
    )
    _plot(
        "jetski_sbl.png", "SBL fix in the array frame",
        measurements["SBL"], axes, "m",
    )
