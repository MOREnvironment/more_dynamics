"""``PayloadBuilder``, payload names, the unset-value refusal and
``frozen_block`` (``plugins/shared/payload_io.py``).

A plugin carries no vehicle numbers: a declared value is ``unset()`` until the
composition gives it, and the plugin's configuration callback refuses with one
message that lists every missing value. The decorator ``refuses_unset`` wraps
``on_configure(context)`` of an ``rpp_py`` plugin, whose result is a
``CallbackResult`` (a stand-in with the same fields is used here, so the tests
do not import ``rpp_py``). The context is a plain stand-in with
``get_parameter``; no component tree is needed.
"""
from dataclasses import dataclass
from types import SimpleNamespace

import casadi as ca
import numpy as np
import pytest

from more_common.casadi_graph import RppCasadiGraph
from more_dynamics.models.shared.wiring import function_from
from more_dynamics.plugins.shared.payload_io import (
    PAYLOAD_NAMES, PayloadBuilder, UnsetValues, frozen_block, payload_name, read_declared, read_value,
    refuses_unset)
from more_transformations.more_casadi_transformations import Parameter, unset

LENGTH = Parameter("length", (1, 1), "m", "waterline length", minimum=0.0, minimum_exclusive=True)
MASS = Parameter("hull_mass", (1, 1), "kg", "dry mass", minimum=0.0, minimum_exclusive=True)
MATRIX = Parameter("stiffness", (2, 2), "N/m", "stiffness matrix, BODY")


class Context:
    def __init__(self, **values):
        self.values = values

    def get_parameter(self, name):
        return self.values[name]


def _payload():
    return SimpleNamespace(inputDescription=[], outputDescription=[], stateDescription=[], dynamics=b"", output=b"")


# ---- names -------------------------------------------------------------------------------------------------------

def test_payload_names_translate_model_symbols_and_keep_others():
    assert payload_name("eta") == "pose" == PAYLOAD_NAMES["eta"]
    assert payload_name("nu_r") == "relative_velocity"
    assert payload_name("something_else") == "something_else"


def test_payload_name_rename_wins_over_the_table():
    assert payload_name("eta", {"eta": "position"}) == "position"


def test_payload_names_are_one_word_per_quantity():
    assert len(set(PAYLOAD_NAMES.values())) == len(PAYLOAD_NAMES)


# ---- reading values ----------------------------------------------------------------------------------------------

def test_read_declared_returns_every_declared_value():
    values = read_declared(Context(length=5.0, hull_mass=10.0), (LENGTH, MASS))
    assert values == {"length": 5.0, "hull_mass": 10.0}


def test_read_declared_refuses_an_unset_value_with_unit_and_meaning():
    with pytest.raises(UnsetValues) as caught:
        read_declared(Context(length=unset(), hull_mass=10.0), (LENGTH, MASS))
    text = str(caught.value)
    assert "'length' [m]" in text and "waterline length" in text and "hull_mass" not in text


def test_zero_and_empty_text_are_values_not_the_marker():
    assert read_declared(Context(length=0.0), (LENGTH,)) == {"length": 0.0}


def test_read_value_refuses_unset_and_passes_a_value():
    assert read_value(Context(count=2), "count", "1", "number of hulls") == 2
    with pytest.raises(UnsetValues, match=r"'count' \[1\]"):
        read_value(Context(count=unset()), "count", "1", "number of hulls")


# ---- the refusal decorator ---------------------------------------------------------------------------------------

@dataclass(frozen=True)
class CallbackResult:
    """Same fields as ``rpp_py.lifecycle.CallbackResult``."""
    success: bool = True
    message: str = ""


class Hull:
    @refuses_unset
    def on_configure(self, context):
        read_declared(context, (LENGTH, MASS, MATRIX))
        self.built = True
        return CallbackResult()


def test_decorator_passes_the_result_of_a_complete_configuration():
    plugin = Hull()
    context = Context(length=5.0, hull_mass=1.0, stiffness=[[1.0, 0.0], [0.0, 1.0]])
    assert plugin.on_configure(context) == CallbackResult(True, "")
    assert plugin.built


def test_decorator_passes_a_failed_result_through():
    class Refusing:
        @refuses_unset
        def on_configure(self, context):
            return CallbackResult(False, "no")
    assert Refusing().on_configure(Context()) == CallbackResult(False, "no")


def test_decorator_lists_every_missing_value_at_once():
    context = Context(length=5.0, hull_mass=unset(), stiffness=unset())
    with pytest.raises(UnsetValues) as caught:
        Hull().on_configure(context)
    text = str(caught.value)
    assert text.startswith("Hull: not set in the composition -- missing parameter(s): ")
    assert "'hull_mass' [kg]" in text and "'stiffness' [N/m]" in text and "'length'" not in text


def test_decorator_runs_the_callback_a_second_time_to_collect_the_names():
    calls = []

    @refuses_unset
    def configure(self, context):
        calls.append(1)
        read_declared(context, (LENGTH,))
    with pytest.raises(UnsetValues):
        configure(object(), Context(length=unset()))
    assert len(calls) == 2  # the first refusal, then the collecting pass (a side effect happens twice)


def test_decorator_keeps_the_callback_name_and_docstring():
    def on_configure(self, context):
        """doc"""
    wrapped = refuses_unset(on_configure)
    assert wrapped.__name__ == "on_configure" and wrapped.__doc__ == "doc"


def test_decorator_does_not_hide_other_errors():
    @refuses_unset
    def configure(self, context):
        raise ValueError("a different problem")
    with pytest.raises(ValueError, match="a different problem"):
        configure(object(), Context())


# ---- freezing ----------------------------------------------------------------------------------------------------

def test_frozen_block_fixes_the_declared_parameters_to_the_plugin_values():
    x, k = ca.SX.sym("x"), ca.SX.sym("hull_mass")
    block = function_from("scale", {"x": x, "hull_mass": k}, {"y": k * x})
    frozen = frozen_block(Context(hull_mass=4.0), block, (MASS,))
    assert frozen.name_in() == ["x"]
    assert frozen(x=3.0)["y"] == 12.0


def test_frozen_block_without_declared_parameters_is_the_block_itself():
    x = ca.SX.sym("x")
    block = function_from("id", {"x": x}, {"y": x})
    assert frozen_block(Context(), block, ()) is block


def test_frozen_block_refuses_an_unset_value():
    x, k = ca.SX.sym("x"), ca.SX.sym("hull_mass")
    block = function_from("scale", {"x": x, "hull_mass": k}, {"y": k * x})
    with pytest.raises(UnsetValues):
        frozen_block(Context(hull_mass=unset()), block, (MASS,))


# ---- the builder -------------------------------------------------------------------------------------------------

def _built_graph():
    builder = PayloadBuilder(_payload())
    eta = builder.input("pose", 6, "position and Euler angles, NED")
    again = builder.input("pose", 6)
    assert again is eta
    gain = builder.input("gain", (2, 2), "gain matrix")
    z = builder.state("z", 2, [1.0, 2.0], "filter state")
    builder.output("scaled", gain @ z, "gain times state")
    builder.output("first", eta[0], "north position")
    builder.state_dot(-z + gain @ z)
    return builder.payload()


def test_builder_declares_names_sizes_and_initial_conditions():
    payload = _built_graph()
    assert [(d.name, d.size) for d in payload.inputDescription] == [("pose", 6), ("gain", 4)]
    assert [(d.name, d.size) for d in payload.outputDescription] == [("scaled", 2), ("first", 1)]
    assert [(d.name, d.size, list(d.ic)) for d in payload.stateDescription] == [("z", 2, [1.0, 2.0])]


def test_builder_output_and_dynamics_take_one_stacked_state_and_input_and_a_matrix_is_column_major():
    graph = RppCasadiGraph(_built_graph())
    gain = np.array([[1.0, 2.0], [3.0, 4.0]])
    stacked = np.concatenate([np.arange(6.0), gain.ravel(order="F")])  # pose, then vec(gain)
    state = np.array([1.0, 1.0])
    out = np.asarray(graph.output(state, stacked)).ravel()
    assert np.allclose(out, [*(gain @ state), 0.0])
    dot = np.asarray(graph.step(state, stacked)).ravel()
    assert np.allclose(dot, -state + gain @ state)


def test_builder_refuses_one_name_asked_with_two_shapes():
    builder = PayloadBuilder(_payload())
    builder.input("pose", 6)
    with pytest.raises(ValueError, match="pose"):
        builder.input("pose", (3, 2))


def test_builder_without_dynamics_has_no_step():
    builder = PayloadBuilder(_payload())
    u = builder.input("u", 1)
    builder.output("y", 2 * u)
    payload = builder.payload()
    assert payload.dynamics == b""
    assert RppCasadiGraph(payload).step is None


def test_builder_call_names_a_blocks_inputs_by_payload_name_and_shares_symbols():
    nu, eta = ca.SX.sym("nu", 6), ca.SX.sym("eta", 6)
    block = function_from("blk", {"nu": nu, "eta": eta}, {"out": nu + eta})
    builder = PayloadBuilder(_payload())
    out = builder.call(block)
    builder.output("out", out["out"])
    names = [d.name for d in builder.payload().inputDescription]
    assert names == ["velocity", "pose"]


def test_builder_call_given_input_does_not_become_a_payload_input():
    nu, eta = ca.SX.sym("nu", 6), ca.SX.sym("eta", 6)
    block = function_from("blk", {"nu": nu, "eta": eta}, {"out": nu + eta})
    builder = PayloadBuilder(_payload())
    out = builder.call(block, given={"eta": ca.SX.zeros(6)})
    builder.output("out", out["out"])
    assert [d.name for d in builder.payload().inputDescription] == ["velocity"]
