"""Fin flow angle absent: beta = 0 and the section speed squared from the
fin-velocity components along the chord and lift axes (MSS ``remus100.m``
234-245 @ cc07579). Inputs ``fin_velocity``, ``chord_axis``, ``lift_axis``,
outputs ``flow_angle`` and ``speed_squared``.

References
----------
[MSS] Fossen, T. I. MSS, MIT, CRAFT/AUV/models/remus100.m 234-245 @ cc07579.

Author:    Enio Krizman
Date:      2026-10-08
"""
from rpp_plugin_types.more_dynamics import FinFlowAngle
from rpp_py.context import ComponentContext

from more_dynamics.models.force_producers.fin.flow_angle.none import none_casadi, none_parameters
from more_dynamics.plugins.payload_io import OPEN_PARAMETERS, PayloadBuilder, frozen_block, payload_name


class FinFlowAngleNone(FinFlowAngle):
    PARAMETERS = [OPEN_PARAMETERS]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        self._model = frozen_block(context, none_casadi(), none_parameters())

    def graph(self) -> FinFlowAngle.CasadyPayload:
        if self._model is None:
            raise RuntimeError("FinFlowAngleNone must be initialized before graph()")
        io = PayloadBuilder(FinFlowAngle.CasadyPayload())
        out = io.call(self._model)
        for name in self._model.name_out():
            io.output(payload_name(name), out[name])
        return io.payload()
