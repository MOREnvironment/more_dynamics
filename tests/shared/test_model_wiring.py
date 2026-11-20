"""The mechanics shared by every model block: ``function_from``, ``restrict``,
``with_passthrough`` and ``as_producer`` (``models/shared/wiring.py``,
``as_producer.py``). The checks are structural (names, sizes, dependence) and
by value on small closed forms; no frame or sign enters.
"""
import casadi as ca
import numpy as np
import pytest

from more_dynamics.models.shared.as_producer import as_producer
from more_dynamics.models.shared.prescribed_wrench import PRESCRIBED_WRENCH_PARAMETERS, prescribed_wrench
from more_dynamics.models.shared.wiring import function_from, restrict, with_passthrough


def _block():
    """Inputs a (2), b (1), c (3); outputs s = a + b, t = 2 c (depends on c only)."""
    a, b, c = ca.SX.sym("a", 2), ca.SX.sym("b"), ca.SX.sym("c", 3)
    return function_from("block", {"a": a, "b": b, "c": c}, {"s": a + b, "t": 2 * c})


def test_function_from_keeps_names_and_order():
    f = _block()
    assert f.name_in() == ["a", "b", "c"]
    assert f.name_out() == ["s", "t"]
    out = f(a=[1.0, 2.0], b=10.0, c=[1.0, 2.0, 3.0])
    assert np.allclose(out["s"].full().ravel(), [11.0, 12.0])
    assert np.allclose(out["t"].full().ravel(), [2.0, 4.0, 6.0])


def test_function_from_allows_a_primitive_handed_on_under_its_own_name():
    x = ca.SX.sym("hull_mass")
    f = function_from("pass", {"hull_mass": x}, {"hull_mass": x})
    assert f(hull_mass=3.0)["hull_mass"] == 3.0


def test_restrict_drops_the_inputs_the_kept_outputs_do_not_read():
    cut = restrict(_block(), ["t"], "cut")
    assert cut.name_in() == ["c"]
    assert cut.name_out() == ["t"]
    assert np.allclose(cut(c=[1.0, 1.0, 1.0])["t"].full().ravel(), [2.0, 2.0, 2.0])


def test_restrict_keeps_every_input_an_output_reads():
    cut = restrict(_block(), ["s"], "cut")
    assert cut.name_in() == ["a", "b"]


def test_restrict_renames_an_output_to_the_slot_coupling_name():
    cut = restrict(_block(), ["t"], "cut", rename={"t": "coupling"})
    assert cut.name_out() == ["coupling"]


def test_restrict_unknown_output_is_refused():
    with pytest.raises(Exception):
        restrict(_block(), ["missing"], "cut")


def test_with_passthrough_hands_the_listed_inputs_on_unchanged():
    f = with_passthrough(_block(), "wide", ["b"])
    assert f.name_out() == ["s", "t", "b"]
    assert f(a=[0.0, 0.0], b=7.0, c=[0.0, 0.0, 0.0])["b"] == 7.0


def test_as_producer_renames_the_command_input_and_keeps_listed_outputs():
    u, k = ca.SX.sym("u", 2), ca.SX.sym("k")
    block = function_from("blk", {"u": u, "k": k}, {"tau": ca.vertcat(u, 0, 0, 0, 0) * k, "extra": k})
    producer = as_producer(block, "u", "producer")
    assert producer.name_in() == ["command", "k"]
    assert producer.name_out() == ["tau"]
    tau = producer(command=[1.0, 2.0], k=3.0)["tau"].full().ravel()
    assert np.allclose(tau, [3.0, 6.0, 0.0, 0.0, 0.0, 0.0])


def test_prescribed_wrench_is_the_command_unchanged():
    f = prescribed_wrench()
    wrench_command = [1.0, -2.0, 3.0, -4.0, 5.0, -6.0]
    assert f.name_in() == ["command"] and f.name_out() == ["tau"]
    assert np.allclose(f(command=wrench_command)["tau"].full().ravel(), wrench_command)
    assert PRESCRIBED_WRENCH_PARAMETERS == ()
