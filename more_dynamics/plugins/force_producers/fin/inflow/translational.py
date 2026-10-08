"""Fin inflow from the vehicle's translation alone, v_fin = v_r[0:3] (MSS
``remus100.m`` 234-235 @ cc07579). Inputs ``relative_velocity`` and
``fin_position`` (unused), output ``fin_velocity``.

References
----------
[MSS] Fossen, T. I. MSS, MIT, CRAFT/AUV/models/remus100.m 234-235 @ cc07579.

Author:    Enio Krizman
Date:      2026-10-08
"""
from rpp_plugin_types.more_dynamics import FinInflow
from rpp_py.context import ComponentContext

from more_dynamics.models.force_producers.fin.inflow.translational import (translational_casadi,
                                                                                translational_parameters)
from more_dynamics.plugins.payload_io import OPEN_PARAMETERS, PayloadBuilder, frozen_block, payload_name


class FinInflowTranslational(FinInflow):
    PARAMETERS = [OPEN_PARAMETERS]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        self._model = frozen_block(context, translational_casadi(), translational_parameters())

    def graph(self) -> FinInflow.CasadyPayload:
        if self._model is None:
            raise RuntimeError("FinInflowTranslational must be initialized before graph()")
        io = PayloadBuilder(FinInflow.CasadyPayload())
        out = io.call(self._model)
        for name in self._model.name_out():
            io.output(payload_name(name), out[name])
        return io.payload()
