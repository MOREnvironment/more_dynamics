"""Rigid body of a homogeneous prolate spheroid (spheroid.m 35-42).

Author:    Enio Krizman
Date:      2026-10-08
"""
from rpp_plugin_types.more_dynamics import RigidBodyModel
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_dynamics.models.vehicles.hull_parts.rigid_body.rigid_body_parts import homogeneous_spheroid, homogeneous_spheroid_parameters
from more_dynamics.plugins.shared.payload_io import OPEN_PARAMETERS, PayloadBuilder, frozen_block, payload_name


class HomogeneousSpheroid(RigidBodyModel):
    PARAMETERS = [
        ParameterDescription("body_density", 1054.7872613500267),  # mass 31.9 kg (remus100.m:3, given) / (4/3 pi a b^2); MSS: 1025 at spheroid.m:35
        ParameterDescription("body_center_of_gravity", [0, 0, 0.02]),  # remus100.m:137 r_bG
        OPEN_PARAMETERS,
    ]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        self._model = frozen_block(context, homogeneous_spheroid(), homogeneous_spheroid_parameters())

    def graph(self) -> RigidBodyModel.CasadyPayload:
        if self._model is None:
            raise RuntimeError("HomogeneousSpheroid must be initialized before graph()")
        io = PayloadBuilder(RigidBodyModel.CasadyPayload())
        out = io.call(self._model)
        for name in self._model.name_out():
            io.output(payload_name(name), out[name])
        return io.payload()
