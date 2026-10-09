"""Type-interface and vehicle-feed checks of the fin skeleton and the
force-producer set.

Written 2026-10-08, before the parts exist (contract: ``fin_parts_contract``).

* Every part form has the coupling inputs and outputs of its type, names,
  sizes and order, and then its own parameters; nothing else.
* ``lifting_fin_casadi`` refuses, naming the slot: a missing or unknown
  slot, a part without a coupling name of its type (input or output), a
  coupling of the wrong size, and a part with an input ``water_density``.
  (The rpp context builder instantiates a child by its plugin name and does
  not compare its type with the slot's, so the skeleton is where a wrong
  part is caught.)
* ``check_lifting_fin_values`` refuses a chord or lift axis that is not a
  unit vector and axes that are not perpendicular: the sign and size of
  the lift are the lift axis.
* The vehicle feeds ``nu_r`` and ``water_density`` (owner's answers of
  2026-10-08): the skeleton and the set take both as inputs, no part or
  skeleton declares ``water_density``, the set refuses a producer without
  either input (naming ``producers[i]``), and the set has one density for
  all its producers, so two values cannot enter.
* The set stacks its producers' states and their rates in list order and
  refuses a command map of the wrong shape.

Author:    Enio Krizman
Date:      2026-10-08
"""

import casadi as ca
import numpy as np
import pytest

from fin_parts_contract import (
    COUPLINGS,
    DECLARED,
    E_X,
    E_Y,
    FIN_OUTPUTS,
    NEG_Y,
    NEG_Z,
    PART_MODULES,
    PRODUCER_INPUTS,
    PRODUCER_OUTPUTS,
    SKELETON_DECLARED,
    call,
    fin_from_forms,
    geometry,
    max_diff,
    module,
    part,
    part_block,
    part_declared,
    producer_set,
    rng,
    skeleton,
    tau,
)


def _fake(name, inputs, outputs):
    """A stand-in part: inputs {name: rows}, outputs {name: expression of the inputs' first entries}."""
    syms = {n: ca.SX.sym(n, rows) for n, rows in inputs.items()}
    first = syms[next(iter(syms))][0] if syms else ca.SX(1.0)
    outs = {n: ca.repmat(first, rows, 1) for n, rows in outputs.items()}
    return ca.Function(name, list(syms.values()), list(outs.values()), list(syms), list(outs))


def _good_parts():
    return {"servo": part("servo", "ideal", {"max_deflection": 0.35}),
            "inflow": part("inflow", "translational"),
            "flow_angle": part("flow_angle", "none"),
            "interference": part("interference", "none"),
            "section": part("section", "quadratic_drag", {"lift_slope": 0.5})}


def _deflection_only_fin(position=(-0.8, 0, 0), lift=NEG_Y, servo="ideal", servo_values=None):
    return fin_from_forms({"servo": servo, "inflow": "translational", "flow_angle": "none",
                           "interference": "none", "section": "quadratic_drag"},
                          {"servo": servo_values or {"max_deflection": 0.35}, "section": {"lift_slope": 0.5}},
                          geometry(list(position), lift, 0.0133, E_X))


# --------------------------------------------------------------------------
# Parts: couplings of their type
# --------------------------------------------------------------------------
@pytest.mark.parametrize("slot, form", sorted(PART_MODULES))
def test_part_has_the_couplings_of_its_type_then_its_parameters(slot, form):
    f = part_block(slot, form)
    inputs, outputs = COUPLINGS[slot]
    declared = [d.name for d in part_declared(slot, form)]
    assert declared == list(DECLARED[(slot, form)])
    assert f.name_in() == [*inputs, *declared], f.name_in()
    assert f.name_out() == list(outputs), f.name_out()
    for name, rows in inputs.items():
        if rows is not None:
            assert f.size_in(name) == (rows, 1), (name, f.size_in(name))
    for name, rows in outputs.items():
        if rows is not None:
            assert f.size_out(name) == (rows, 1), (name, f.size_out(name))
    assert "water_density" not in f.name_in()


@pytest.mark.parametrize("slot, form", sorted(PART_MODULES))
def test_part_is_one_casadi_function(slot, form):
    """One CasADi function per part (all-CasADi blocks)."""
    assert isinstance(part_block(slot, form), ca.Function)


# --------------------------------------------------------------------------
# Skeleton
# --------------------------------------------------------------------------
def test_skeleton_declares_its_geometry_and_no_density():
    declared = skeleton().lifting_fin_parameters()
    assert [d.name for d in declared] == list(SKELETON_DECLARED)
    for d in declared:
        shape, unit, minimum, exclusive = SKELETON_DECLARED[d.name]
        assert (d.shape, d.unit, d.minimum, d.minimum_exclusive) == (shape, unit, minimum, exclusive), d
        assert d.meaning


def test_skeleton_signature():
    f = skeleton().lifting_fin_casadi(_good_parts())
    names = [d.name for d in skeleton().lifting_fin_parameters()]
    assert f.name_in() == [*PRODUCER_INPUTS, *names], f.name_in()
    assert f.name_out() == FIN_OUTPUTS
    assert f.size_in("command") == (1, 1) and f.size_in("nu_r") == (6, 1) and f.size_in("water_density") == (1, 1)
    assert f.size_in("state") == (0, 1) and f.size_out("tau") == (6, 1)


def test_skeleton_state_is_the_servo_state():
    parts = {**_good_parts(), "servo": part("servo", "lag_rate_angle",
                                            {"max_deflection": 0.35, "max_rate": 0.5, "time_constant": 0.1})}
    f = skeleton().lifting_fin_casadi(parts)
    assert f.size_in("state") == (1, 1) and f.size_out("state_dot") == (1, 1)


@pytest.mark.parametrize("change, words", [
    (lambda p: p.pop("section"), ["section"]),
    (lambda p: p.update(rudder_trim=p["interference"]), ["rudder_trim"]),
    (lambda p: p.update(inflow=_fake("bad_inflow", {"nu_r": 6}, {"fin_velocity": 3})), ["inflow", "fin_position"]),
    (lambda p: p.update(flow_angle=_fake("bad_angle", {"fin_velocity": 3, "chord_axis": 3, "lift_axis": 3},
                                         {"flow_angle": 1})), ["flow_angle", "speed_squared"]),
    (lambda p: p.update(section=_fake("bad_section", {"angle_of_attack": 2},
                                      {"lift_coefficient": 1, "drag_coefficient": 1})), ["section", "angle_of_attack"]),
    (lambda p: p.update(inflow=_fake("bad_inflow", {"nu_r": 3, "fin_position": 3}, {"fin_velocity": 3})), ["inflow", "nu_r"]),
    (lambda p: p.update(section=_fake("dense_section", {"angle_of_attack": 1, "water_density": 1},
                                      {"lift_coefficient": 1, "drag_coefficient": 1})), ["section", "water_density"]),
])
def test_skeleton_refuses_a_part_that_does_not_fit_its_slot(change, words):
    parts = _good_parts()
    change(parts)
    with pytest.raises(ValueError) as error:
        skeleton().lifting_fin_casadi(parts)
    for word in words:
        assert word in str(error.value), (word, str(error.value))


@pytest.mark.parametrize("change, name", [
    ({"chord_axis": [1.0, 0.1, 0.0]}, "chord_axis"),
    ({"lift_axis": [0.0, -2.0, 0.0]}, "lift_axis"),
    ({"lift_axis": [0.6, -0.8, 0.0]}, "lift_axis"),         # unit, not perpendicular to the chord
    ({"fin_area": 0.0}, "fin_area"),
    ({"fin_position": [-0.8, 0.0]}, "fin_position"),
])
def test_skeleton_refuses_invalid_geometry(change, name):
    values = {**geometry([-0.8, 0, 0], NEG_Y, 0.0133, E_X), **change}
    with pytest.raises(ValueError, match=name):
        skeleton().check_lifting_fin_values(values)


def test_skeleton_accepts_a_rotated_orthonormal_frame():
    c = np.array([np.cos(0.3), np.sin(0.3), 0.0])
    n = np.array([-np.sin(0.3), np.cos(0.3), 0.0])
    skeleton().check_lifting_fin_values(geometry([-0.8, 0, 0], n.tolist(), 0.0133, c.tolist()))


# --------------------------------------------------------------------------
# Set and vehicle feed
# --------------------------------------------------------------------------
def test_set_signature_has_one_density():
    fin_set = producer_set([_deflection_only_fin(), _deflection_only_fin(lift=NEG_Z)], np.eye(2))
    assert fin_set.name_in() == PRODUCER_INPUTS, fin_set.name_in()
    assert fin_set.name_out()[:2] == PRODUCER_OUTPUTS
    assert sum(name.endswith("water_density") for name in fin_set.name_in()) == 1


def test_one_fed_density_scales_every_lifting_fin():
    fin_set = producer_set([_deflection_only_fin(), _deflection_only_fin((-0.7, 0, 0), NEG_Z)], np.eye(2))
    g = rng()
    for command, nu_r in zip(g.uniform(-0.4, 0.4, (50, 2)), g.uniform(-2, 2, (50, 6))):
        a, b = tau(fin_set, command, nu_r, 1000.0), tau(fin_set, command, nu_r, 1030.0)
        assert max_diff(b, 1.03 * a) <= 1e-12 * max(1.0, float(np.max(np.abs(b))))


@pytest.mark.parametrize("missing", ["nu_r", "water_density"])
def test_set_refuses_a_producer_without_a_fed_input(missing):
    inputs = {"command": 1, "state": 0, "nu_r": 6, "water_density": 1}
    inputs.pop(missing)
    syms = {n: ca.SX.sym(n, rows) for n, rows in inputs.items()}
    bad = ca.Function("bad", list(syms.values()), [ca.SX.zeros(6), ca.SX(0, 1)], list(syms), ["tau", "state_dot"])
    with pytest.raises(ValueError) as error:
        module("shared.force_producer_set").force_producer_set_casadi([_deflection_only_fin(), bad], 2)
    assert "producers[1]" in str(error.value) and missing in str(error.value), str(error.value)


def test_set_refuses_a_command_map_of_the_wrong_shape():
    with pytest.raises(ValueError, match="command_map"):
        producer_set([_deflection_only_fin(), _deflection_only_fin(lift=NEG_Z)], np.ones((3, 2)))


def test_set_stacks_states_in_list_order():
    lag = {"max_deflection": 0.35, "max_rate": 5.0, "time_constant": 0.1}
    first = _deflection_only_fin(servo="lag_rate_angle", servo_values=lag)
    second = _deflection_only_fin((-0.7, 0, 0), NEG_Z, servo="ideal")
    third = _deflection_only_fin((-0.7, 0, 0), E_Y, servo="lag_rate_angle", servo_values={**lag, "time_constant": 0.2})
    fin_set = producer_set([first, second, third], np.eye(3))
    assert fin_set.size_in("state") == (2, 1)
    out = call(fin_set, command=[0.2, 0.1, -0.3], state=[0.0, 0.0], nu_r=[1.5, 0, 0, 0, 0, 0], water_density=1026.0)
    assert max_diff(out["state_dot"], [0.2 / 0.1, -0.3 / 0.2]) <= 1e-12
