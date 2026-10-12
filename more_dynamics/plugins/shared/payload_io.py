"""Building and reading ``CasadyPayload`` (vehicle_model.capnp 51-58) the way
Luka's plugins do: ``output(state, input)`` and ``dynamics(state, input)``
take one stacked ``state`` and one stacked ``input`` vector; the names and
sizes of the entries live only in ``inputDescription`` / ``outputDescription``
/ ``stateDescription`` (``jet_nozzle.py`` 50-156,
``linear_surface_hydrostatics.py`` 54-77).

A matrix crosses the payload as its column-major ``vec`` (an ``IODescription``
carries a size, not a shape); a consumer that declares a 6x6 input of size 36
reshapes it back the same way.

Payload names are words, as in Luka's descriptions (``pose``, ``velocity``,
``restoring_force``); a model block keeps its own symbols (``eta``, ``nu_r``,
``M_RB``) and the plugin names them once, through ``PAYLOAD_NAMES``.

A plugin carries no vehicle's numbers: every value of a vehicle, a device or a
site is declared with ``unset()`` (an empty list, ``more_transformations``) and
comes from the composition's named part. ``read_declared`` and ``read_value``
refuse an unset value, and ``refuses_unset`` (the decorator of a plugin's
``initialize``) turns the refusal into one message that lists every value of
the component that its selected forms read and the composition left unset,
with unit and meaning. A convention of this library, no source equation.

Author:    Enio Krizman
Date:      2026-10-08
"""

import contextvars
import functools

import casadi as ca

from more_common.casadi_graph import graph_to_bytes
from more_transformations.more_casadi_transformations import Parameter, freeze, is_unset


class UnsetValues(ValueError):
    """A declared value the composition has not set (its default is still ``unset()``)."""


_COLLECTING = contextvars.ContextVar("unset_values_collector", default=None)


def _placeholder(parameter):
    """Some value inside the declared range, for the pass that only lists what is missing."""
    low, high = parameter.minimum, parameter.maximum
    candidates = [1.0, 0.5, 0.0, 2.0, 10.0]
    if low is not None:
        candidates += [low + 1.0, low + 0.5]
    if high is not None:
        candidates += [high - 1.0, high - 0.5]
    for value in candidates:
        if ((low is None or value > low or (value == low and not parameter.minimum_exclusive))
                and (high is None or value < high or (value == high and not parameter.maximum_exclusive))):
            return [[value] * parameter.shape[1] for _ in range(parameter.shape[0])]
    return [[1.0] * parameter.shape[1] for _ in range(parameter.shape[0])]


def read_declared(context, declared):
    """``{name: value}`` of the plugin's values of the ``declared`` parameters, refusing an unset one
    (``UnsetValues``, naming it with its unit and meaning)."""
    declared = tuple(declared)
    values = {d.name: context.get_parameter(d.name) for d in declared}
    unset = [d for d in declared if is_unset(values[d.name])]
    if unset:
        collector = _COLLECTING.get()
        if collector is None:
            raise UnsetValues(_missing_text(unset))
        for d in unset:
            collector.setdefault(d.name, d)
            values[d.name] = _placeholder(d)
    return values


def read_value(context, name, unit, meaning, placeholder=1.0):
    """The plugin's value of ``name`` where it is read directly (a count, a given matrix) and not through a
    model block; an unset one is refused like ``read_declared``'s."""
    value = context.get_parameter(name)
    if is_unset(value):
        declaration = Parameter(name, (1, 1), unit, meaning)
        collector = _COLLECTING.get()
        if collector is None:
            raise UnsetValues(_missing_text([declaration]))
        collector.setdefault(name, declaration)
        return placeholder
    return value


def _missing_text(declared):
    return "missing parameter(s): " + ", ".join(f"{d.describe()} ({d.meaning})" for d in declared)


def refuses_unset(initialize):
    """Decorator of a plugin's ``initialize``: when a value its selected forms read is still unset, the error
    names the plugin and lists every such value at once (a second pass lists the rest, with a placeholder in
    each, before the plugin refuses; no graph is made)."""
    @functools.wraps(initialize)
    def wrapper(self, context):
        try:
            return initialize(self, context)
        except UnsetValues as first:
            error = first
        collected = {}
        token = _COLLECTING.set(collected)
        try:
            initialize(self, context)
        except Exception:  # a placeholder may not fit a later check; only the names are wanted
            pass
        finally:
            _COLLECTING.reset(token)
        listed = _missing_text(list(collected.values())) if collected else str(error)
        raise UnsetValues(f"{type(self).__name__}: not set in the composition -- {listed}; set them in the "
                          "part's params/parameters.py (a reference part shows working values)") from None
    return wrapper


def frozen_block(context, block, declared):
    """``block`` with every declared parameter frozen to the plugin's value of it."""
    declared = tuple(declared)
    if not declared:
        return block
    return freeze(block, declared, read_declared(context, declared))


# model-layer symbol -> payload name (one name per quantity across every part)
PAYLOAD_NAMES = {
    "eta": "pose", "nu": "velocity", "nu_r": "relative_velocity",
    "M_RB": "rigid_body_mass_matrix", "C_RB": "rigid_body_coriolis_matrix",
    "M_A": "added_mass_matrix", "C_A": "added_mass_coriolis_matrix",
    "G": "restoring_matrix", "G_CF": "restoring_matrix_at_flotation", "D": "damping_matrix",
}


def payload_name(model_name, rename=None):
    rename = rename or {}
    return rename.get(model_name, PAYLOAD_NAMES.get(model_name, model_name))


class PayloadBuilder:
    """Collects named inputs, states and outputs, then serialises Luka's two
    positional functions."""

    def __init__(self, payload):
        self._payload = payload
        self._inputs, self._states, self._outputs, self._state_dot = [], [], [], None

    def input(self, name, shape, description=""):
        """A named input; asking again for the same name and shape gives the one symbol (several models of one
        plugin read the same quantity)."""
        rows, cols = (shape, 1) if isinstance(shape, int) else shape
        for known, symbol, _ in self._inputs:
            if known == name:
                if symbol.size() != (rows, cols):
                    raise ValueError(f"input {name!r} asked as {rows}x{cols} and as {symbol.size1()}x{symbol.size2()}")
                return symbol
        symbol = ca.SX.sym(name, rows, cols)
        self._inputs.append((name, symbol, description))
        return symbol

    def state(self, name, size, ic, description=""):
        symbol = ca.SX.sym(name, size)
        self._states.append((name, symbol, list(ic), description))
        return symbol

    def output(self, name, expression, description=""):
        self._outputs.append((name, ca.vec(ca.SX(expression)), description))

    def state_dot(self, expression):
        self._state_dot = ca.vec(ca.SX(expression))

    def call(self, function, rename=None, given=None):
        """Call a named model block, every input it still has becoming a
        payload input named by ``payload_name`` (unless ``given``)."""
        given = dict(given or {})
        args = {}
        for i, name in enumerate(function.name_in()):
            args[name] = given[name] if name in given else self.input(payload_name(name, rename),
                                                                       function.size_in(i))
        return function.call(args)

    def payload(self):
        payload = self._payload
        for name, symbol, description in self._inputs:
            payload.inputDescription.append(_io(name, symbol.numel(), description))
        for name, expression, description in self._outputs:
            payload.outputDescription.append(_io(name, expression.numel(), description))
        for name, symbol, ic, description in self._states:
            payload.stateDescription.append(_state(name, symbol.numel(), ic, description))
        n_in = sum(s.numel() for _, s, _ in self._inputs)
        n_state = sum(s.numel() for _, s, _, _ in self._states)
        state, inputs = ca.SX.sym("state", n_state), ca.SX.sym("input", n_in)
        substitute = _splitter([s for _, s, _, _ in self._states], state) + \
            _splitter([s for _, s, _ in self._inputs], inputs)
        originals = [s for _, s, _, _ in self._states] + [s for _, s, _ in self._inputs]
        outputs = ca.vertcat(*[e for _, e, _ in self._outputs])
        output = ca.substitute([outputs], originals, substitute)[0] if originals else outputs
        payload.output = graph_to_bytes(ca.Function("output", [state, inputs], [output], ["state", "input"],
                                                    ["output"]))
        if self._state_dot is not None:
            state_dot = ca.substitute([self._state_dot], originals, substitute)[0]
            payload.dynamics = graph_to_bytes(ca.Function("dynamics", [state, inputs], [state_dot],
                                                          ["state", "input"], ["state_dot"]))
        return payload


def _io(name, size, description):
    from rpp_schema.more_dynamics.IODescription import IODescription
    return IODescription(size, name=name, description=description)


def _state(name, size, ic, description):
    from rpp_schema.more_dynamics.StateDescription import StateDescription
    return StateDescription(size, ic=ic, name=name, description=description)


def _splitter(symbols, stacked):
    parts, k = [], 0
    for s in symbols:
        n = s.numel()
        parts.append(ca.reshape(stacked[k:k + n], s.size1(), s.size2()))
        k += n
    return parts
