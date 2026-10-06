import casadi as ca

from rpp_plugin_types.more_dynamics import HydrodynamicsModel
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription
from rpp_schema.more_dynamics.IODescription import IODescription

from more_common.casadi_graph import graph_to_bytes
from more_dynamics.models.hydrodynamics import (
    crossflow_surface_hydrodynamics_casadi,
    preprocess_crossflow_surface_hydrodynamics,
)


class CrossflowSurfaceHydrodynamics(HydrodynamicsModel):
    """Linear surface-vessel damping plus Fossen/Hoerner cross-flow drag."""

    PARAMETERS = [
        ParameterDescription(
            "damping_coefficients",
            [240.0, 1800.0, 7000.0, 900.0, 6500.0, 8000.0],
        ),
        ParameterDescription("length", 3.5),
        ParameterDescription("beam", 1.2),
        ParameterDescription("draft", 0.16),
        ParameterDescription("water_density", 1025.0),
        ParameterDescription("crossflow_scale", 1.0),
        ParameterDescription("strip_count", 20),
    ]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        properties = preprocess_crossflow_surface_hydrodynamics(
            damping_coefficients=context.get_parameter("damping_coefficients"),
            length=context.get_parameter("length"),
            beam=context.get_parameter("beam"),
            draft=context.get_parameter("draft"),
            water_density=context.get_parameter("water_density"),
            crossflow_scale=context.get_parameter("crossflow_scale"),
            strip_count=context.get_parameter("strip_count"),
        )
        self._model = crossflow_surface_hydrodynamics_casadi(properties)

    def graph(self) -> HydrodynamicsModel.CasadyPayload:
        if self._model is None:
            raise RuntimeError(
                "CrossflowSurfaceHydrodynamics must be initialized before graph()"
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
