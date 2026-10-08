"""Restoring of two surface hulls at their equilibrium draft, nabla = m / rho,
T = nabla / (2 Cb L B_pont) (otter.m 121-122), as Luka's
``HydrostaticsModel``: first output ``restoring_force = -g`` (added by the
vehicle, ``vehicle_model.capnp`` 11-13); then the restoring matrices, the
displaced volume and the draft for the parts that consume them.

Author:    Enio Krizman
Date:      2026-10-08
"""
from rpp_plugin_types.more_dynamics import HydrostaticsModel
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_dynamics.models.vehicles.hull_parts.restoring.restoring_parts import (
    TWIN_HULL_METACENTRIC_EQUILIBRIUM_DRAFT_PARAMETERS, twin_hull_metacentric_equilibrium_draft)
from more_dynamics.plugins.shared.payload_io import OPEN_PARAMETERS, PayloadBuilder, frozen_block, payload_name


class TwinHullMetacentricEquilibriumDraft(HydrostaticsModel):
    PARAMETERS = [
        ParameterDescription("longitudinal_inertia_factor", 0.8),  # otter.m 179 (I_L = 0.8 ...)
        ParameterDescription("longitudinal_center_of_flotation", -0.2),  # otter.m 177, 192
        ParameterDescription("reference_point", [0.0, 0.0, 0.0]),  # the body-frame origin (CO)
        OPEN_PARAMETERS,
    ]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        self._model = frozen_block(context, twin_hull_metacentric_equilibrium_draft(),
                                   TWIN_HULL_METACENTRIC_EQUILIBRIUM_DRAFT_PARAMETERS)

    def graph(self) -> HydrostaticsModel.CasadyPayload:
        if self._model is None:
            raise RuntimeError("TwinHullMetacentricEquilibriumDraft must be initialized before graph()")
        io = PayloadBuilder(HydrostaticsModel.CasadyPayload())
        out = io.call(self._model)
        io.output("restoring_force", -out["g"], "signed generalized force added by the vehicle, BODY, N and N m")
        for name in self._model.name_out():
            if name != "g":
                io.output(payload_name(name), out[name])
        return io.payload()
