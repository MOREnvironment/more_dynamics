"""Fin section with linear lift and a constant zero-lift drag: C_L = c_L_alpha
alpha, C_D = C_D0 (Prestero 2001, eq. 4.37, p. 31). The lift slope defaults to
the REMUS fin's (Table A.5, p. 103); the zero-lift drag defaults to zero
because the cited equation carries no fin drag, and is replaced by a measured
value when one exists. Input ``angle_of_attack``, outputs ``lift_coefficient``
and ``drag_coefficient``.

References
----------
[Prestero 2001] Prestero, T. (2001). Verification of a six-degree of freedom
    simulation model for the REMUS autonomous underwater vehicle. MIT/WHOI
    MSc thesis. Eq. 4.37, p. 31; Table A.5, p. 103.

Author:    Enio Krizman
Date:      2026-10-08
"""
from rpp_plugin_types.more_dynamics import FinSection
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_dynamics.models.force_producers.fin.section.linear_section import (linear_section_casadi,
                                                                                   linear_section_parameters)
from more_dynamics.plugins.shared.payload_io import OPEN_PARAMETERS, PayloadBuilder, frozen_block, payload_name


class FinSectionLinearSection(FinSection):
    PARAMETERS = [
        ParameterDescription("lift_slope", 3.12),  # 1/rad, Prestero 2001, Table A.5, p. 103 (c_L_alpha)
        ParameterDescription("zero_lift_drag", 0.0),  # Prestero 2001, eq. 4.37, p. 31 carries no fin drag
        OPEN_PARAMETERS,
    ]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        self._model = frozen_block(context, linear_section_casadi(), linear_section_parameters())

    def graph(self) -> FinSection.CasadyPayload:
        if self._model is None:
            raise RuntimeError("FinSectionLinearSection must be initialized before graph()")
        io = PayloadBuilder(FinSection.CasadyPayload())
        out = io.call(self._model)
        for name in self._model.name_out():
            io.output(payload_name(name), out[name])
        return io.payload()
