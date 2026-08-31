import casadi as ca

from rpp_plugin_types.more_dynamics import HydrostaticsModel
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription
from rpp_schema.more_dynamics.IODescription import IODescription

from more_dynamics.models.hydrostatics import (
    linear_surface_hydrostatics_casadi,
    preprocess_linear_surface_hydrostatics,
)
from more_common.casadi_graph import graph_to_bytes


class LinearSurfaceHydrostatics(HydrostaticsModel):
    PARAMETERS = [
        ParameterDescription("length", 5.2),
        ParameterDescription("beam", 2.15),
        ParameterDescription("draft", 0.3),
        ParameterDescription("block_coefficient", 0.233),
        ParameterDescription("waterplane_coefficient", 0.8),
        ParameterDescription("water_density", 1025.0),
        ParameterDescription("gravity", 9.81),
        ParameterDescription("center_of_gravity", [0.0, 0.0, 0.025]),
        ParameterDescription("longitudinal_center_of_flotation", -0.3),
        ParameterDescription("coefficient_scales", [1.0, 0.7, 1.0]),
        ParameterDescription("reference_point", [0.0, 0.0, 0.0]),
    ]

    def __init__(self) -> None:
        self._properties = None
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        self._properties = preprocess_linear_surface_hydrostatics(
            length=context.get_parameter("length"),
            beam=context.get_parameter("beam"),
            draft=context.get_parameter("draft"),
            block_coefficient=context.get_parameter("block_coefficient"),
            waterplane_coefficient=context.get_parameter(
                "waterplane_coefficient"
            ),
            water_density=context.get_parameter("water_density"),
            gravity=context.get_parameter("gravity"),
            center_of_gravity=context.get_parameter("center_of_gravity"),
            longitudinal_center_of_flotation=context.get_parameter(
                "longitudinal_center_of_flotation"
            ),
            coefficient_scales=context.get_parameter("coefficient_scales"),
            reference_point=context.get_parameter("reference_point"),
        )
        self._model = linear_surface_hydrostatics_casadi(self._properties)

    def graph(self) -> HydrostaticsModel.CasadyPayload:
        if self._model is None:
            raise RuntimeError(
                "LinearSurfaceHydrostatics must be initialized before graph()"
            )

        payload = HydrostaticsModel.CasadyPayload()
        payload.inputDescription.append(IODescription(6, name="pose"))
        payload.outputDescription.append(
            IODescription(6, name="restoring_force")
        )
        payload.outputDescription.append(
            IODescription(1, name="displaced_volume")
        )
        payload.outputDescription.append(
            IODescription(1, name="wetted_surface_area")
        )

        state = ca.SX.sym("state", 0)
        pose = ca.SX.sym("input", 6)
        payload.output = graph_to_bytes(
            ca.Function("output", [state, pose], [self._model(pose)])
        )
        return payload
