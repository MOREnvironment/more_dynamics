"""The jet as the first composition of the propulsor: the shared repository's jet nozzle exactly, and the fixed jet.

Written 2026-10-10, before the parts exist (contract: ``propulsor_contract``).

The shared repository's ``JetNozzle`` is left as it is and ours stands alone:
nothing here imports it. "Equal to his" is read from rows recorded once from
his plugin (``data/jet_nozzle_reference.csv``, ``data/SOURCE.md`` says how);
``test_propulsor_live_reference.py`` runs his live plugin when it is importable.

* His forms: drive ``normalized_lag`` + steering ``jet_nozzle_clamp`` +
  mounting ``jet_nozzle_direction`` + inflow ``none`` + law
  ``normalized_linear`` + housing ``none`` equal his output (the wrench) and his
  dynamics (thrust and nozzle-angle derivatives) on the 600 frozen rows, three
  parameter sets (trim 0, trim 0.2 rad, a negative trim with a yaw bias and
  asymmetric angle limits), tolerance 1e-12. Command order [thrust, angle],
  state order [thrust, nozzle angle]: his.
* His behaviours recorded, not endorsed: a full astern effort pushes
  backwards as hard as ahead (a waterjet cannot do this without a reversing
  bucket, which turns the jet down and forward: Varyani and Krishnankutty 2009,
  J. Mar. Eng. Technol. 13, p. 14), and with trim his thrust vector is longer
  than the commanded thrust.
* The fixed jet (the default): our servo for the nozzle angle and our rotation
  for the direction. At trim 0 its wrench equals his on every frozen row whose
  angle state is inside the limits (the fix changes nothing there); at trim 0.2
  rad the difference is quantified and must be positive, and the thrust vector
  has the commanded length at every angle. Its dynamics equal the his-forms
  composition's under the conditions of ``test_propulsor_steering_parts``
  (symmetric limits, 1 s time constant, rate limit equal to his nozzle
  velocity, state inside the limits).

Author:    Enio Krizman
Date:      2026-10-10
"""

import numpy as np
import pytest

from propulsor_contract import (EXACT_TOLERANCE, JET_CSV, JET_CSV_SHA256, DATA_DIR, jet_forms_fixed, jet_forms_his,
                                jet_parts, jet_reference, jet_values_fixed, jet_values_his, max_diff,
                                propulsor_from_forms, rng, run, sha256)

NU_R = np.zeros(6)
RHO = 1025.0
SETS = sorted(jet_parts())


def _his(name):
    return propulsor_from_forms(jet_forms_his(), jet_values_his(jet_parts()[name]))


def _fixed(name, **changes):
    numbers = {**jet_parts()[name], **changes}
    return propulsor_from_forms(jet_forms_fixed(), jet_values_fixed(numbers))


def test_the_frozen_rows_are_the_recorded_ones():
    assert sha256(DATA_DIR / JET_CSV) == JET_CSV_SHA256


@pytest.mark.parametrize("name", SETS)
def test_his_forms_equal_his_output_and_dynamics_on_the_frozen_rows(name):
    ref = jet_reference()[name]
    f = _his(name)
    assert f.size_in("command") == (2, 1) and f.size_in("state") == (2, 1)
    worst_out = worst_dot = 0.0
    for state, command, output, dot in zip(ref.state, ref.command, ref.output, ref.dot):
        out = run(f, command, NU_R, RHO, state=state)
        worst_out = max(worst_out, max_diff(out["tau"], output))
        worst_dot = max(worst_dot, max_diff(out["state_dot"], dot))
    assert worst_out <= EXACT_TOLERANCE and worst_dot <= EXACT_TOLERANCE, (worst_out, worst_dot)


@pytest.mark.parametrize("name", SETS)
def test_his_rows_cover_the_clamps_and_both_directions(name):
    ref = jet_reference()[name]
    numbers = jet_parts()[name]
    assert np.any(np.abs(ref.state[:, 0]) > 1.0) and np.any(ref.state[:, 1] > numbers["max_angle"])
    assert np.any(ref.state[:, 1] < numbers["min_angle"])
    assert np.any(ref.command[:, 0] > 1.0) and np.any(ref.command[:, 1] < -1.0)
    assert np.any(ref.dot[:, 0] > 0) and np.any(ref.dot[:, 0] < 0)


def test_the_density_fed_by_the_vehicle_does_not_change_his_jet():
    ref = jet_reference()["trimmed"]
    f = _his("trimmed")
    for state, command in list(zip(ref.state, ref.command))[:20]:
        a = run(f, command, NU_R, 900.0, state=state)["tau"]
        b = run(f, command, NU_R, 1100.0, state=state)["tau"]
        assert max_diff(a, b) == 0.0


def test_his_jet_is_indifferent_to_the_boat_speed():
    f = _his("trimmed")
    a = run(f, [1.0, 0.3], NU_R, RHO, state=[1.0, 0.2])["tau"]
    b = run(f, [1.0, 0.3], [8.0, 1.0, 0.0, 0.1, 0.0, 0.2], RHO, state=[1.0, 0.2])["tau"]
    assert max_diff(a, b) == 0.0


def test_his_jet_pushes_astern_as_hard_as_ahead_a_recorded_behaviour_not_ours():
    """Effort -1 gives -T_max c straight astern (a waterjet without a reversing bucket cannot)."""
    numbers = jet_parts()["straight"]
    f = _his("straight")
    ahead = run(f, [0, 0], NU_R, RHO, state=[1.0, 0.0])["tau"][0]
    astern = run(f, [0, 0], NU_R, RHO, state=[-1.0, 0.0])["tau"][0]
    assert abs(ahead - numbers["max_thrust"] * numbers["thrust_coefficient"]) <= 1e-12 and abs(astern + ahead) <= 1e-12


def test_his_thrust_vector_is_longer_than_the_commanded_thrust_when_trimmed():
    numbers = jet_parts()["trimmed"]
    f = _his("trimmed")
    full = numbers["max_thrust"] * numbers["thrust_coefficient"]
    lengths = [np.linalg.norm(run(f, [0, 0], NU_R, RHO, state=[1.0, a])["tau"][:3]) / full
               for a in np.linspace(numbers["min_angle"], numbers["max_angle"], 21)]
    assert max(lengths) > 1.003 and min(lengths) >= 1.0 - 1e-15          # up to ~0.7 %, never shorter


# --------------------------------------------------------------------------
# The fixed jet
# --------------------------------------------------------------------------
def _inside(ref, numbers):
    return (ref.state[:, 1] >= numbers["min_angle"]) & (ref.state[:, 1] <= numbers["max_angle"])


def _angle_command(command, numbers):
    """His normalised angle command as radians: the one thing our servo takes differently."""
    c = np.clip(command, -1.0, 1.0)
    return c * numbers["max_angle"] if c >= 0 else -c * numbers["min_angle"]


def test_at_zero_trim_the_fixed_jet_equals_his_wrench_wherever_the_nozzle_is_inside_its_limits():
    ref = jet_reference()["straight"]
    numbers = jet_parts()["straight"]
    f = _fixed("straight")
    mask = _inside(ref, numbers)
    assert mask.sum() > 100
    worst = 0.0
    for state, command, output in zip(ref.state[mask], ref.command[mask], ref.output[mask]):
        out = run(f, [command[0], _angle_command(command[1], numbers)], NU_R, RHO, state=state)
        worst = max(worst, max_diff(out["tau"], output))
    assert worst <= EXACT_TOLERANCE, worst


def test_the_fixed_jets_thrust_vector_has_the_commanded_length_at_every_trim_and_angle():
    numbers = jet_parts()["trimmed"]
    f = _fixed("trimmed")
    full = numbers["max_thrust"] * numbers["thrust_coefficient"]
    for trim in (0.0, 0.2, -0.3, 0.7):
        g = _fixed("trimmed", trim_angle=trim)
        for a in np.linspace(numbers["min_angle"], numbers["max_angle"], 15):
            tau = run(g, [0, 0], NU_R, RHO, state=[1.0, a])["tau"]
            assert abs(np.linalg.norm(tau[:3]) - full) <= 1e-9, (trim, a)
    assert f.size_in("state") == (2, 1)


def test_the_fix_is_real_at_trim_0_2_rad_the_difference_is_quantified():
    ref = jet_reference()["trimmed"]
    numbers = jet_parts()["trimmed"]
    his, fixed = _his("trimmed"), _fixed("trimmed")
    mask = _inside(ref, numbers)
    worst = worst_relative = 0.0
    full = numbers["max_thrust"] * numbers["thrust_coefficient"]
    for state, command in zip(ref.state[mask], ref.command[mask]):
        c = [command[0], _angle_command(command[1], numbers)]
        a = run(his, command, NU_R, RHO, state=state)["tau"]
        b = run(fixed, c, NU_R, RHO, state=state)["tau"]
        worst = max(worst, float(np.linalg.norm(a[:3] - b[:3])))
        worst_relative = max(worst_relative, float(np.linalg.norm(a[:3] - b[:3]) / full))
    assert worst > 1.0, worst                    # newtons at full thrust: the largest force difference over the rows
    assert worst_relative < 0.02, worst_relative  # a few per mille to ~1 %: a correction, not another jet


def test_the_fixed_jets_servo_dynamics_equal_the_his_forms_under_the_reduction_conditions():
    """Symmetric limits, 1 s time constant, rate limit = his nozzle velocity, state inside the limits, angle command
    his normalised command times the limit: the nozzle-angle derivative equals his clamp's; the thrust state is the
    same drive part."""
    g = rng()
    for numbers in ({**jet_parts()["straight"], "min_angle": -0.6, "max_angle": 0.6, "nozzle_velocity": 0.3},
                    {**jet_parts()["biased"], "min_angle": -0.5, "max_angle": 0.5, "nozzle_velocity": 1.5}):
        his = propulsor_from_forms(jet_forms_his(), jet_values_his(numbers))
        ours = propulsor_from_forms(jet_forms_fixed(), jet_values_fixed(numbers))
        worst = 0.0
        for _ in range(300):
            state = [g.uniform(-1.2, 1.2), g.uniform(-numbers["max_angle"], numbers["max_angle"])]
            command = [g.uniform(-1.2, 1.2), g.uniform(-1.0, 1.0)]
            a = run(his, command, NU_R, RHO, state=state)["state_dot"]
            b = run(ours, [command[0], command[1] * numbers["max_angle"]], NU_R, RHO, state=state)["state_dot"]
            worst = max(worst, max_diff(a, b))
        assert worst <= EXACT_TOLERANCE, (numbers, worst)


def test_the_fixed_jet_is_the_default_jet_and_names_no_vehicle():
    forms = jet_forms_fixed()
    assert forms["steering"] == "servo_lag_rate_angle" and forms["mounting"] == "steered_rotation"
    assert all("jet_nozzle" not in v for k, v in forms.items() if k in ("steering", "mounting"))
