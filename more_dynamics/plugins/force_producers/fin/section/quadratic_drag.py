"""Fin section with linear lift and quadratic deflection drag, the MSS
comparison form: C_L = c_L_alpha alpha, C_D = c_L_alpha alpha^2 (MSS
``remus100.m`` 238-245 @ cc07579). The lift slope defaults to the REMUS fin's
(Prestero 2001, Table A.5, p. 103). Input ``angle_of_attack``, outputs
``lift_coefficient`` and ``drag_coefficient``.

References
----------
[Prestero 2001] Prestero, T. (2001). Verification of a six-degree of freedom
    simulation model for the REMUS autonomous underwater vehicle. MIT/WHOI
    MSc thesis. Table A.5, p. 103.
[MSS] Fossen, T. I. MSS, MIT, CRAFT/AUV/models/remus100.m 238-245 @ cc07579.

Author:    Enio Krizman
Date:      2026-10-08
"""
from rpp_plugin_types.more_dynamics import FinSection
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_dynamics.models.force_producers.fin.section.quadratic_drag import (quadratic_drag_casadi,
                                                                                   quadratic_drag_parameters)
from more_dynamics.plugins.payload_io import OPEN_PARAMETERS, PayloadBuilder, frozen_block, payload_name


class FinSectionQuadraticDrag(FinSection):
    PARAMETERS = [
        ParameterDescription("lift_slope", 3.12),  # 1/rad, Prestero 2001, Table A.5, p. 103 (c_L_alpha)
        OPEN_PARAMETERS,
    ]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        self._model = frozen_block(context, quadratic_drag_casadi(), quadratic_drag_parameters())

    def graph(self) -> FinSection.CasadyPayload:
        if self._model is None:
            raise RuntimeError("FinSectionQuadraticDrag must be initialized before graph()")
        io = PayloadBuilder(FinSection.CasadyPayload())
        out = io.call(self._model)
        for name in self._model.name_out():
            io.output(payload_name(name), out[name])
        return io.payload()
