"""No fin-body interference: both factors, ``deflection_factor`` and
``flow_angle_factor``, equal one (Prestero 2001, eq. 4.41, p. 32; MSS
``remus100.m`` 238-245 @ cc07579). No inputs.

References
----------
[Prestero 2001] Prestero, T. (2001). Verification of a six-degree of freedom
    simulation model for the REMUS autonomous underwater vehicle. MIT/WHOI
    MSc thesis. Eq. 4.41, p. 32.
[MSS] Fossen, T. I. MSS, MIT, CRAFT/AUV/models/remus100.m 238-245 @ cc07579.

Author:    Enio Krizman
Date:      2026-10-08
"""
from rpp_plugin_types.more_dynamics import FinInterference
from rpp_py.context import ComponentContext

from more_dynamics.models.force_producers.fin.interference.none import none_casadi, none_parameters
from more_dynamics.plugins.shared.payload_io import OPEN_PARAMETERS, PayloadBuilder, frozen_block, payload_name


class FinInterferenceNone(FinInterference):
    PARAMETERS = [OPEN_PARAMETERS]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        self._model = frozen_block(context, none_casadi(), none_parameters())

    def graph(self) -> FinInterference.CasadyPayload:
        if self._model is None:
            raise RuntimeError("FinInterferenceNone must be initialized before graph()")
        io = PayloadBuilder(FinInterference.CasadyPayload())
        out = io.call(self._model)
        for name in self._model.name_out():
            io.output(payload_name(name), out[name])
        return io.payload()
