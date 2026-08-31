from types import SimpleNamespace

import casadi as ca

from more_common.casadi_graph import RppCasadiGraph, graph_to_bytes


def test_graph_without_dynamics_has_no_step_function():
    value = ca.SX.sym("value", 1)
    output = ca.Function("output", [value], [value])
    payload = SimpleNamespace(
        inputDescription=[],
        outputDescription=[],
        stateDescription=[],
        dynamics=b"",
        output=graph_to_bytes(output),
    )

    graph = RppCasadiGraph(payload)

    assert graph.step is None
    assert graph.output(2.0) == 2.0
