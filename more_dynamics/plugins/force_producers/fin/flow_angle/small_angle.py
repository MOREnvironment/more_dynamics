"""Fin flow angle in the small-angle form: beta = u_n / u_c and the section
speed squared u_c^2 (Prestero 2001, eqs. 4.42-4.43, pp. 32-33). Inputs
``fin_velocity``, ``chord_axis``, ``lift_axis`` (from the fin), outputs
``flow_angle`` and ``speed_squared``.

References
----------
[Prestero 2001] Prestero, T. (2001). Verification of a six-degree of freedom
    simulation model for the REMUS autonomous underwater vehicle. MIT/WHOI
    MSc thesis. Eqs. 4.42-4.43, pp. 32-33.

Author:    Enio Krizman
Date:      2026-10-08
"""
from rpp_plugin_types.more_dynamics import FinFlowAngle
from rpp_py.context import ComponentContext

from more_dynamics.models.force_producers.fin.flow_angle.small_angle import small_angle_casadi, small_angle_parameters
from more_dynamics.plugins.payload_io import OPEN_PARAMETERS, PayloadBuilder, frozen_block, payload_name


class FinFlowAngleSmallAngle(FinFlowAngle):
    PARAMETERS = [OPEN_PARAMETERS]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        self._model = frozen_block(context, small_angle_casadi(), small_angle_parameters())

    def graph(self) -> FinFlowAngle.CasadyPayload:
        if self._model is None:
            raise RuntimeError("FinFlowAngleSmallAngle must be initialized before graph()")
        io = PayloadBuilder(FinFlowAngle.CasadyPayload())
        out = io.call(self._model)
        for name in self._model.name_out():
            io.output(payload_name(name), out[name])
        return io.payload()
