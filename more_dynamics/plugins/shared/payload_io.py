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

Author:    Enio Krizman
Date:      2026-10-08
"""

import casadi as ca

from more_common.casadi_graph import RppCasadiGraph, graph_to_bytes
from more_transformations.more_casadi_transformations import freeze
from rpp_py.parameter_description import ParameterDescription

# Every part declares it: the names of its own parameters left as graph inputs
# (identification, a time-varying current in a gate). Default: none, every
# parameter frozen in initialize().
OPEN_PARAMETERS = ParameterDescription("open_parameters", [])


def frozen_block(context, block, declared):
    """``block`` with every declared parameter frozen to its composition
    value except the ones the composition lists in ``open_parameters``."""
    declared = tuple(declared)
    names = {d.name for d in declared}
    opened = list(context.get_parameter("open_parameters") or [])
    unknown = [n for n in opened if n not in names]
    if unknown:
        raise ValueError(f"open_parameters {unknown} are not parameters of block {block.name()!r}: {sorted(names)}")
    kept = tuple(d for d in declared if d.name not in opened)
    if not kept:
        return block
    return freeze(block, kept, {d.name: context.get_parameter(d.name) for d in kept})

# model-layer symbol -> payload name (one name per quantity across every part)
PAYLOAD_NAMES = {
    "eta": "pose", "nu": "velocity", "nu_r": "relative_velocity",
    "M_RB": "rigid_body_mass_matrix", "C_RB": "rigid_body_coriolis_matrix",
    "M_A": "added_mass_matrix", "C_A": "added_mass_coriolis_matrix",
    "G": "restoring_matrix", "G_CF": "restoring_matrix_at_flotation", "D": "damping_matrix",
}


MODEL_NAMES = {payload: model for model, payload in PAYLOAD_NAMES.items()}


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


def named_entries(descriptions, vector):
    """``{name: slice of vector}`` in description order (sizes add up)."""
    out, k = {}, 0
    for d in descriptions:
        out[d.name] = vector[k:k + d.size]
        k += d.size
    if k != vector.numel():
        raise ValueError(f"descriptions cover {k} values, the vector has {vector.numel()}")
    return out


def named_function(payload, name, rename=None, state_name=None):
    """A child's payload as one ``ca.Function`` with named inputs and outputs:
    what a composite (a fin, the cross-flow strip) hands to its model-layer
    assembly, which speaks the model symbols (``rename``: payload name ->
    model name). With ``state_name`` the child's stacked state is the input
    ``state_name`` and its derivative the output ``<state_name>_dot`` (empty
    when the child has no state)."""
    rename = rename or {}
    graph = RppCasadiGraph(payload)
    state = ca.SX.sym("state", graph.num_states)
    ins = {d.name: ca.SX.sym(d.name, d.size) for d in payload.inputDescription}
    stacked = ca.vertcat(*ins.values()) if ins else ca.SX(0, 1)
    outs = named_entries(payload.outputDescription, graph.output(state, stacked))
    inputs = dict(ins)
    if state_name is not None:
        inputs = {state_name: state} | inputs
        outs[f"{state_name}_dot"] = graph.step(state, stacked) if graph.step is not None else ca.SX(0, 1)
    elif graph.num_states:
        raise ValueError(f"{name}: a child with states needs a state name")
    return ca.Function(name, list(inputs.values()), list(outs.values()),
                       [rename.get(n, n) for n in inputs], [rename.get(n, n) for n in outs])
