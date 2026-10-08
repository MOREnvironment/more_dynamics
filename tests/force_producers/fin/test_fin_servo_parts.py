"""Gate tests of the servo part: one device, ``ideal``, ``lag_rate_angle``,
``output_saturated`` and ``no_limits`` as its named settings of the three
switches ``dynamics``, ``rate_limit``, ``angle_limit`` (``fin.servo.servo``).

Written 2026-10-08. The servo is one consolidated part, not one file per
form; there is no rate-only servo (no part, no "1:1" or "DUNE" claim
anywhere). The servo type is
shared by a fin and by a propulsor's steering. Coupling inputs ``command``
(1, rad) and ``servo_state``; outputs ``deflection`` (1, rad) and
``servo_state_dot``.

Settings and their sources (keys in References):

* ``ideal`` (``dynamics="none"``, ``angle_limit="on_command"``):
  ``delta = sat(command, +-max_deflection)``, no state (MSS ``remus100.m``
  109, 113-114). G1 against MATLAB running MSS.
* ``no_limits`` (``dynamics="none"``, ``angle_limit="none"``): the output is
  the command itself, no state, no limit at all.
* ``lag_rate_angle`` (``dynamics="first_order_lag"``, ``rate_limit=True``,
  ``angle_limit="on_command"``): one state, the deflection;
  ``d(delta)/dt = sat((sat(command, +-max_deflection) - delta) /
  time_constant, +-max_rate)`` — amplitude limit, first-order lag, rate
  limit (Murray-Smith 2016, Fig. 1 and text, p. 246; time constant = 1/G_r,
  p. 246). The command is an angle, so the command gain G_a is 1 (p. 246:
  "If the required rudder deflection delta_c is the variable subjected to
  limiting, the gain factor G_a is unity").
* ``output_saturated`` (``dynamics="first_order_lag"``, ``rate_limit=True``,
  ``angle_limit="on_output"``): one state, the integrator's own output,
  not itself clamped; ``d(x)/dt = sat((command - x) / time_constant,
  +-max_rate)``, ``delta = sat(x, +-max_deflection)`` (Sarhadi 2026, Fig. 4,
  p. 4: lag -> rate limit -> integrator -> amplitude limit on the output).
  Equals ``lag_rate_angle`` only while the command never exceeds the limit
  and the state starts inside it; a command held beyond the limit winds the
  state up past the stop, and the output leaves the stop late after a
  reversal.

Invalid switch combinations are refused, naming the switch: a rate limit
needs the lag's own state (``dynamics="none"`` has none), and
``angle_limit="on_output"`` needs it too (there is no state distinct from
the output to wind up without one).

The closed-form step response used below is a derivation of this module
from the ``lag_rate_angle`` equation: from rest, with target ``d =
sat(command)``, the rate limit holds while ``|d - delta| > time_constant *
max_rate``, i.e. up to ``t1 = (|d| - time_constant * max_rate) / max_rate``;
after ``t1`` (or from the start when ``|d| <= time_constant * max_rate``)
the lag is linear: ``d - delta = (d - delta(t1)) exp(-(t - t1) /
time_constant)``.

The wind-up release times used below (0.667 s, 5.00 s) are derived here,
not published values: command held at twice the angle limit for 30 s, then
reversed to minus twice the limit; the output-side state winds up unclamped
and the output leaves the stop ``max_deflection / max_rate`` seconds after
the reversal, deep inside the rate-limited phase: 20 deg / 30 deg/s =
0.667 s with the servo of Sarhadi (2026, Fig. 4, p. 4) and 35 deg / 7 deg/s
= 5.00 s with the Zeefakkel actuator of Murray-Smith (2016, p. 246)
(`tests/force_producers/fin/data/SOURCE.md`).

References
----------
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2
    with the fixes of 2026-10-07. https://github.com/cybergalactic/MSS, MIT
    licence, revision ``cc07579``: ``CRAFT/AUV/models/remus100.m`` 109,
    113-114.
[Murray-Smith 2016] Murray-Smith, D. J. (2016). Inverse simulation methods
    applied to investigations of actuator nonlinearities in ship steering.
    *Simulation Notes Europe* 26(4), 245-256. Fig. 1 and text, p. 246;
    Table 1, p. 247 (R.O.V. Zeefakkel: time constant 3 s, +-35 deg,
    +-7 and +-10 deg/s, p. 246).
[Sarhadi 2026] Sarhadi, P. (2026). Simple yet effective anti-windup
    techniques for amplitude and rate saturation: an AUV case study.
    arXiv:2601.01302v2, Fig. 4, p. 4 (read as a page image).

Author:    Enio Krizman
Date:      2026-10-08
"""

import casadi as ca
import numpy as np
import pytest

from fin_parts_contract import (
    DECLARED,
    G1_TOLERANCE,
    G2_TOLERANCE,
    N_RANDOM,
    call,
    module,
    mss_constant,
    mss_fins,
    part,
    part_declared,
    rng,
)

FORMS = ("ideal", "no_limits", "lag_rate_angle", "output_saturated")
DEG = np.pi / 180.0


def _servo(form, **values):
    return part("servo", form, values)


def _integrate(servo, delta0, command, times):
    """delta(t) of a servo with one state under a constant command (CVODES)."""
    x, c = ca.SX.sym("x"), ca.SX.sym("c")
    rate = servo(command=c, servo_state=x)["servo_state_dot"]
    integrator = ca.integrator("servo", "cvodes", {"x": x, "p": c, "ode": rate}, 0.0, list(times),
                               {"abstol": 1e-13, "reltol": 1e-13, "max_num_steps": 200000})
    return np.array(integrator(x0=delta0, p=command)["xf"], dtype=float).reshape(-1)


def _step_response(command, max_deflection, max_rate, time_constant, t):
    """Closed form of the lag_rate_angle equation from rest (module docstring)."""
    target = np.clip(command, -max_deflection, max_deflection)
    sign, size = np.sign(target), abs(target)
    t1 = max(0.0, (size - time_constant * max_rate) / max_rate)
    t = np.asarray(t, dtype=float)
    ramp = sign * max_rate * t
    lag = target - (target - sign * max_rate * t1) * np.exp(-(t - t1) / time_constant)
    return np.where(t <= t1, ramp, lag)


# --------------------------------------------------------------------------
# Declarations (rule: every parameter a primitive with unit, meaning, range)
# --------------------------------------------------------------------------
@pytest.mark.parametrize("form", FORMS)
def test_servo_declares_its_primitives(form):
    declared = part_declared("servo", form)
    expected = DECLARED[("servo", form)]
    assert [d.name for d in declared] == list(expected)
    for d in declared:
        shape, unit, minimum, exclusive = expected[d.name]
        assert (d.shape, d.unit, d.minimum, d.minimum_exclusive) == (shape, unit, minimum, exclusive), d
        assert d.meaning


@pytest.mark.parametrize("form", [f for f in FORMS if DECLARED[("servo", f)]])
def test_servo_refuses_values_outside_their_ranges(form):
    from more_transformations.more_casadi_transformations import check_values

    declared = part_declared("servo", form)
    good = {d.name: 0.5 for d in declared}
    check_values(declared, good)
    for d in declared:
        for bad in (0.0, -0.1, float("nan")):
            with pytest.raises(ValueError, match=d.name):
                check_values(declared, {**good, d.name: bad})


# --------------------------------------------------------------------------
# Invalid switch combinations (rule: refused, naming the switch)
# --------------------------------------------------------------------------
@pytest.mark.parametrize("dynamics, rate_limit, angle_limit, word", [
    ("bogus", False, "on_command", "dynamics"),
    ("none", False, "bogus", "angle_limit"),
    ("none", True, "on_command", "rate_limit"),
    ("none", False, "on_output", "angle_limit"),
])
def test_invalid_switch_combinations_are_refused(dynamics, rate_limit, angle_limit, word):
    servo_module = module("fin.servo.servo")
    with pytest.raises(ValueError, match=word):
        servo_module.servo_parameters(dynamics=dynamics, rate_limit=rate_limit, angle_limit=angle_limit)
    with pytest.raises(ValueError, match=word):
        servo_module.servo_casadi(dynamics=dynamics, rate_limit=rate_limit, angle_limit=angle_limit)


# --------------------------------------------------------------------------
# ideal / no_limits
# --------------------------------------------------------------------------
def test_G1_ideal_servo_equals_mss_saturation():
    """deflection = MATLAB's delta_r, delta_s (remus100.m 113-114) on every
    row, with max_deflection = MATLAB's delta_max (line 109)."""
    ref = mss_fins()
    servo = _servo("ideal", max_deflection=mss_constant(ref, "delta_max"))
    for command, expected in ((ref["ui1"], ref["delta_r"]), (ref["ui2"], ref["delta_s"])):
        got = np.array([call(servo, command=c, servo_state=np.zeros(0))["deflection"][0] for c in command])
        assert np.max(np.abs(got - expected)) <= G1_TOLERANCE
    assert np.any(np.abs(ref["ui1"]) > mss_constant(ref, "delta_max"))   # the limit is reached


def test_ideal_servo_has_no_state():
    servo = _servo("ideal", max_deflection=0.4)
    out = call(servo, command=0.3, servo_state=np.zeros(0))
    assert servo.size1_in("servo_state") == 0 and out["servo_state_dot"].size == 0


def test_no_limits_servo_passes_the_command_through():
    servo = _servo("no_limits")
    assert servo.size1_in("servo_state") == 0
    g = rng()
    for command in g.uniform(-5.0, 5.0, 50):
        out = call(servo, command=command, servo_state=np.zeros(0))
        assert abs(out["deflection"][0] - command) <= G2_TOLERANCE


# --------------------------------------------------------------------------
# lag_rate_angle
# --------------------------------------------------------------------------
def _lag_rate_angle_rate(command, delta, max_deflection, max_rate, time_constant):
    """Murray-Smith 2016, Fig. 1, p. 246 (G_a = 1)."""
    target = np.clip(command, -max_deflection, max_deflection)
    return np.clip((target - delta) / time_constant, -max_rate, max_rate)


def test_G2_lag_rate_angle_equals_the_block_diagram():
    p = {"max_deflection": 35 * DEG, "max_rate": 7 * DEG, "time_constant": 3.0}
    servo = _servo("lag_rate_angle", **p)
    g = rng()
    for command, delta in zip(g.uniform(-1.2, 1.2, N_RANDOM), g.uniform(-0.7, 0.7, N_RANDOM)):
        out = call(servo, command=command, servo_state=[delta])
        assert abs(out["deflection"][0] - delta) <= G2_TOLERANCE              # the state is the deflection
        expected = _lag_rate_angle_rate(command, delta, **p)
        assert abs(out["servo_state_dot"][0] - expected) <= G2_TOLERANCE, (command, delta)


@pytest.mark.parametrize("rate_deg", [7.0, 10.0])
@pytest.mark.parametrize("command_deg", [5.0, 30.0, -30.0, 50.0])
def test_G5_lag_rate_angle_step_response_on_the_zeefakkel_actuator(rate_deg, command_deg):
    """Time constant 3 s, +-35 deg, +-7 or +-10 deg/s (Murray-Smith 2016,
    p. 246): 5 deg stays linear; 30 deg runs into the rate limit; 50 deg is
    cut at 35 deg first. Integrated vs the closed form (module docstring)."""
    p = {"max_deflection": 35 * DEG, "max_rate": rate_deg * DEG, "time_constant": 3.0}
    servo = _servo("lag_rate_angle", **p)
    times = [0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0]
    got = _integrate(servo, 0.0, command_deg * DEG, times)
    expected = _step_response(command_deg * DEG, **p, t=times)
    assert np.max(np.abs(got - expected)) <= 1e-8, (got, expected)


def test_lag_rate_angle_with_an_unbounded_rate_is_a_first_order_lag():
    """max_rate -> infinity: d(delta)/dt = (sat(command) - delta) / time_constant."""
    p = {"max_deflection": 0.4, "max_rate": 1e12, "time_constant": 0.2}
    servo = _servo("lag_rate_angle", **p)
    g = rng()
    for command, delta in zip(g.uniform(-0.6, 0.6, 100), g.uniform(-0.4, 0.4, 100)):
        rate = call(servo, command=command, servo_state=[delta])["servo_state_dot"][0]
        assert abs(rate - (np.clip(command, -0.4, 0.4) - delta) / 0.2) <= G2_TOLERANCE


def test_lag_rate_angle_with_distant_limits_equals_the_unlimited_lag():
    """max_deflection and max_rate set far beyond the command/state range
    used: the command-side limited, rate-limited lag equals the same lag
    with no limits at all (``dynamics="first_order_lag"``,
    ``rate_limit=False``, ``angle_limit="none"``) to machine precision,
    since neither limit ever binds."""
    from more_transformations.more_casadi_transformations import freeze

    servo_module = module("fin.servo.servo")
    far = {"max_deflection": 1.0e6, "max_rate": 1.0e9, "time_constant": 0.2}
    limited = _servo("lag_rate_angle", **far)
    unlimited_block = servo_module.servo_casadi(dynamics="first_order_lag", rate_limit=False, angle_limit="none")
    unlimited_declared = servo_module.servo_parameters(dynamics="first_order_lag", rate_limit=False,
                                                        angle_limit="none")
    unlimited = freeze(unlimited_block, unlimited_declared, {"time_constant": far["time_constant"]})
    g = rng()
    for command, delta in zip(g.uniform(-1.0, 1.0, 100), g.uniform(-1.0, 1.0, 100)):
        a = call(limited, command=command, servo_state=[delta])
        b = call(unlimited, command=command, servo_state=[delta])
        assert abs(a["servo_state_dot"][0] - b["servo_state_dot"][0]) <= G2_TOLERANCE
        assert abs(a["deflection"][0] - b["deflection"][0]) <= G2_TOLERANCE


def test_lag_rate_angle_with_a_tiny_time_constant_approaches_no_dynamics():
    """time_constant -> 0 with an unbinding rate limit: the lag converges to
    its target orders of magnitude faster than any simulated step, so after
    1000 time constants the state equals the instantaneous saturation
    (``ideal`` with the same max_deflection) to the stated tolerance (the
    residual after 1000 time constants of exponential decay, e^-1000, is
    far below double precision; the bound used here is loose by comparison
    on purpose, since the integrator's own tolerance dominates)."""
    p = {"max_deflection": 0.3, "max_rate": 1.0e9, "time_constant": 1.0e-6}
    servo = _servo("lag_rate_angle", **p)
    ideal = _servo("ideal", max_deflection=p["max_deflection"])
    g = rng()
    for command in g.uniform(-0.6, 0.6, 20):
        state = _integrate(servo, 0.0, command, [1000 * p["time_constant"]])[0]
        expected = call(ideal, command=command, servo_state=np.zeros(0))["deflection"][0]
        assert abs(state - expected) <= 1e-6


# --------------------------------------------------------------------------
# output_saturated (Sarhadi 2026, Fig. 4, p. 4)
# --------------------------------------------------------------------------
def _output_saturated_rate(command, state, max_rate, time_constant, **_):
    return np.clip((command - state) / time_constant, -max_rate, max_rate)


def test_command_beyond_the_angle_limit_saturates_the_output_but_not_the_state():
    """A command held well past max_deflection drives the output to the
    stop while the internal state keeps climbing past it (Fig. 4: the
    amplitude limit is on the output, not fed back into the integrator)."""
    p = {"max_deflection": 20 * DEG, "max_rate": 30 * DEG, "time_constant": 0.1}
    servo = _servo("output_saturated", **p)
    x = _integrate(servo, 0.0, 5 * p["max_deflection"], [5.0])[0]
    out = call(servo, command=5 * p["max_deflection"], servo_state=[x])
    assert abs(out["deflection"][0] - p["max_deflection"]) <= 1e-9
    assert x > p["max_deflection"] + 1e-6                           # the state, unlike the output, is not clamped


def test_initial_state_beyond_the_limit_is_clamped_on_the_output_only():
    p = {"max_deflection": 20 * DEG, "max_rate": 30 * DEG, "time_constant": 0.1}
    servo = _servo("output_saturated", **p)
    state0 = 1.5 * p["max_deflection"]
    out = call(servo, command=0.0, servo_state=[state0])
    assert abs(out["deflection"][0] - p["max_deflection"]) <= 1e-9       # output clamped
    expected_rate = max(-p["max_rate"], min(p["max_rate"], (0.0 - state0) / p["time_constant"]))
    assert abs(out["servo_state_dot"][0] - expected_rate) <= G2_TOLERANCE  # state free to move off the stop


def test_G2_output_saturated_equals_its_equation():
    p = {"max_deflection": 20 * DEG, "max_rate": 30 * DEG, "time_constant": 0.1}
    servo = _servo("output_saturated", **p)
    g = rng()
    commands, states = g.uniform(-1.5, 1.5, N_RANDOM), g.uniform(-1.0, 1.0, N_RANDOM)
    rates = np.empty(N_RANDOM)
    for i, (command, state) in enumerate(zip(commands, states)):
        out = call(servo, command=command, servo_state=[state])
        assert abs(out["deflection"][0] - max(-p["max_deflection"], min(p["max_deflection"], state))) <= G2_TOLERANCE
        expected_rate = _output_saturated_rate(command, state, **p)
        assert abs(out["servo_state_dot"][0] - expected_rate) <= G2_TOLERANCE, (command, state)
        rates[i] = out["servo_state_dot"][0]
    assert np.any(rates > 0) and np.any(rates < 0)                      # both rate signs are exercised above


@pytest.mark.parametrize(
    "name, time_constant, max_deflection_deg, max_rate_deg, release_s",
    [
        ("Sarhadi 2026 Fig. 4 p. 4 (tau 0.1 s, +-20 deg, +-30 deg/s)", 0.1, 20.0, 30.0, 20.0 / 30.0),
        ("Murray-Smith 2016 p. 246 (tau 3 s, +-35 deg, +-7 deg/s)", 3.0, 35.0, 7.0, 35.0 / 7.0),
    ],
)
def test_G5_output_saturated_wind_up_release_after_a_saturating_reversal(name, time_constant, max_deflection_deg,
                                                                          max_rate_deg, release_s):
    """Command held at twice the angle limit for 30 s (the unclamped state
    winds up close to it), then reversed to minus twice the limit: the
    output leaves the +max_deflection stop max_deflection / max_rate
    seconds later (20 / 30 = 0.667 s for Sarhadi 2026, Fig. 4, p. 4;
    35 / 7 = 5.00 s for Murray-Smith 2016, p. 246). The state is still deep in the rate-limited phase there
    (|state + 2 max_deflection| = 3 max_deflection >> max_rate *
    time_constant for both parameter sets), so the crossing time does not
    depend on time_constant."""
    max_deflection, max_rate = max_deflection_deg * DEG, max_rate_deg * DEG
    servo = _servo("output_saturated", max_deflection=max_deflection, max_rate=max_rate, time_constant=time_constant)
    x_wound_up = _integrate(servo, 0.0, 2 * max_deflection, [30.0])[0]
    assert abs(x_wound_up - 2 * max_deflection) <= 1e-3                     # converged, unclamped
    x_before, x_after = _integrate(servo, x_wound_up, -2 * max_deflection, [release_s - 0.05, release_s + 0.05])
    deflection_before = call(servo, command=-2 * max_deflection, servo_state=[x_before])["deflection"][0]
    deflection_after = call(servo, command=-2 * max_deflection, servo_state=[x_after])["deflection"][0]
    assert abs(deflection_before - max_deflection) <= 1e-9                  # still at the stop
    assert deflection_after < max_deflection - 1e-6                        # left the stop


def test_output_saturated_equals_lag_rate_angle_when_the_command_never_exceeds_the_limit():
    """Equal only when: the command-side form's gain is 1 (true here, both
    read the command directly); the state starts inside +-max_deflection;
    and the command never exceeds max_deflection. Under those conditions
    the output-side amplitude limit never engages, so both equations
    reduce to the same unsaturated lag + rate law."""
    p = {"max_deflection": 20 * DEG, "max_rate": 30 * DEG, "time_constant": 0.1}
    output_saturated = _servo("output_saturated", **p)
    lag_rate_angle = _servo("lag_rate_angle", **p)
    times = [0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0]
    g = rng()
    for command, state0 in zip(g.uniform(-p["max_deflection"], p["max_deflection"], 20),
                               g.uniform(-p["max_deflection"], p["max_deflection"], 20)):
        x_output = _integrate(output_saturated, state0, command, times)
        x_lag = _integrate(lag_rate_angle, state0, command, times)
        assert max(abs(x_output - x_lag)) <= 1e-8, (command, state0)


# --------------------------------------------------------------------------
# Open-parameter interface and gradient (every settable settings)
# --------------------------------------------------------------------------
@pytest.mark.parametrize("form", ["lag_rate_angle", "output_saturated"])
def test_open_parameters_and_the_gradient_through_time_constant(form):
    open_servo = part("servo", form)          # every own parameter left open
    names = [d.name for d in part_declared("servo", form)]
    assert set(open_servo.name_in()) == {"command", "servo_state", *names}
    command, state = ca.SX.sym("command"), ca.SX.sym("state")
    # max_deflection and max_rate fixed well clear of the command/state/rate
    # used below, so the function stays in its unsaturated, smooth branch.
    fixed = {n: 10.0 for n in names if n != "time_constant"}
    tau = ca.SX.sym("tau")
    out = open_servo(command=command, servo_state=state, time_constant=tau, **fixed)
    f = ca.Function("f", [command, state, tau], [out["servo_state_dot"], ca.jacobian(out["servo_state_dot"], tau)])

    def rate_at(tau_value):
        return float(np.array(f(0.21, 0.2, tau_value)[0]).reshape(-1)[0])

    gradient = float(np.array(f(0.21, 0.2, 0.2)[1]).reshape(-1)[0])
    h = 1e-6
    central = (rate_at(0.2 + h) - rate_at(0.2 - h)) / (2 * h)
    assert abs(gradient) > 1e-6
    assert abs(gradient - central) <= 1e-6
