"""Rigid body of a hull plus a point payload (otter.m 94-98, 122-128).

Author:    Enio Krizman
Date:      2026-10-08
"""
from rpp_plugin_types.more_dynamics import RigidBodyModel
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_dynamics.models.vehicles.hull_parts.rigid_body.rigid_body_parts import hull_with_point_payload, hull_with_point_payload_parameters
from more_dynamics.plugins.payload_io import OPEN_PARAMETERS, PayloadBuilder, frozen_block, payload_name


class HullWithPointPayload(RigidBodyModel):
    PARAMETERS = [
        ParameterDescription("hull_mass", 55.0),  # otter.m 94-98
        ParameterDescription("payload_mass", 25.0),  # generate_catamaran_mss.m 61-62 (payload mass)
        ParameterDescription("hull_center_of_gravity", [0.2, 0.0, -0.2]),  # otter.m 94-98
        ParameterDescription("payload_position", [0.05, 0.0, -0.35]),  # generate_catamaran_mss.m 61-62 (payload position)
        ParameterDescription("radii_of_gyration", [0.4, 0.25, 0.25]),  # otter.m 94-98
        OPEN_PARAMETERS,
    ]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        self._model = frozen_block(context, hull_with_point_payload(), hull_with_point_payload_parameters())

    def graph(self) -> RigidBodyModel.CasadyPayload:
        if self._model is None:
            raise RuntimeError("HullWithPointPayload must be initialized before graph()")
        io = PayloadBuilder(RigidBodyModel.CasadyPayload())
        out = io.call(self._model)
        for name in self._model.name_out():
            io.output(payload_name(name), out[name])
        return io.payload()
