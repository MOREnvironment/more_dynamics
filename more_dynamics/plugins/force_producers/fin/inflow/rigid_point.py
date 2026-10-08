"""Fin inflow at a point fixed to the rigid body, retaining the rotational
terms: v_fin = v_r + omega x r_fin (Prestero 2001, eq. 4.40, p. 32). Inputs
``relative_velocity`` (from the vehicle) and ``fin_position`` (from the fin),
output ``fin_velocity``.

References
----------
[Prestero 2001] Prestero, T. (2001). Verification of a six-degree of freedom
    simulation model for the REMUS autonomous underwater vehicle. MIT/WHOI
    MSc thesis. Eq. 4.40, p. 32.

Author:    Enio Krizman
Date:      2026-10-08
"""
from rpp_plugin_types.more_dynamics import FinInflow
from rpp_py.context import ComponentContext

from more_dynamics.models.force_producers.fin.inflow.rigid_point import rigid_point_casadi, rigid_point_parameters
from more_dynamics.plugins.payload_io import OPEN_PARAMETERS, PayloadBuilder, frozen_block, payload_name


class FinInflowRigidPoint(FinInflow):
    PARAMETERS = [OPEN_PARAMETERS]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        self._model = frozen_block(context, rigid_point_casadi(), rigid_point_parameters())

    def graph(self) -> FinInflow.CasadyPayload:
        if self._model is None:
            raise RuntimeError("FinInflowRigidPoint must be initialized before graph()")
        io = PayloadBuilder(FinInflow.CasadyPayload())
        out = io.call(self._model)
        for name in self._model.name_out():
            io.output(payload_name(name), out[name])
        return io.payload()
