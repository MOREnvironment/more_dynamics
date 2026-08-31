import casadi as ca

from rpp_plugin_types.more_dynamics import HydrodynamicsModel
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription
from rpp_schema.more_dynamics.IODescription import IODescription

from more_common.casadi_graph import graph_to_bytes
from more_dynamics.models.hydrodynamics import (
    linear_surface_hydrodynamics_casadi,
    preprocess_linear_surface_hydrodynamics,
)


class LinearSurfaceHydrodynamics(HydrodynamicsModel):
    PARAMETERS = [
        ParameterDescription(
            "damping_coefficients",
            [
                0.0,
                20.02547625,
                7201.948346468731,
                1562.418311261956,
                14535.734088287643,
                3655.04992515,
            ],
        ),
    ]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        properties = preprocess_linear_surface_hydrodynamics(
            context.get_parameter("damping_coefficients")
        )
        self._model = linear_surface_hydrodynamics_casadi(properties)

    def graph(self) -> HydrodynamicsModel.CasadyPayload:
        if self._model is None:
            raise RuntimeError(
                "LinearSurfaceHydrodynamics must be initialized before graph()"
            )

        payload = HydrodynamicsModel.CasadyPayload()
        payload.inputDescription.append(IODescription(6, name="velocity"))
        payload.outputDescription.append(
            IODescription(6, name="hydrodynamic_force")
        )

        state = ca.SX.sym("state", 0)
        velocity = ca.SX.sym("input", 6)
        payload.output = graph_to_bytes(
            ca.Function("output", [state, velocity], [self._model(velocity)])
        )
        return payload
