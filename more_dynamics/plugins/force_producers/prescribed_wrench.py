"""Stage-1 force producer: the commanded wrench acts on the vehicle unchanged
(no actuator model), the external-wrench boundary of ``hydroVessel.m``.

Author:    Enio Krizman
Date:      2026-10-08
"""
from rpp_plugin_types.more_dynamics import ForceProducer
from rpp_py.context import ComponentContext

from more_dynamics.models.force_producers.force_producer_parts import prescribed_wrench
from more_dynamics.plugins.payload_io import OPEN_PARAMETERS, PayloadBuilder


class PrescribedWrench(ForceProducer):
    PARAMETERS = [OPEN_PARAMETERS]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        self._model = prescribed_wrench()

    def graph(self) -> ForceProducer.CasadyPayload:
        if self._model is None:
            raise RuntimeError("PrescribedWrench must be initialized before graph()")
        io = PayloadBuilder(ForceProducer.CasadyPayload())
        out = io.call(self._model, rename={"command": "desired_wrench"})
        io.output("generated_force", out["tau"], "wrench on the vehicle, BODY, N and N m")
        return io.payload()
