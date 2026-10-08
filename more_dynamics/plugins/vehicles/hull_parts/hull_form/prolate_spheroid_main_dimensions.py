"""Hull form of a prolate spheroid from length and diameter; produces the
shared dimensions every other part reads (remus100.m 131-132, 220-221).
``section_beam``: the beam of a cross-flow strip section, here the diameter.

Author:    Enio Krizman
Date:      2026-10-08
"""
from rpp_plugin_types.more_dynamics import HullForm
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_dynamics.models.vehicles.hull_parts.hull_form.hull_form import (PROLATE_SPHEROID_MAIN_DIMENSIONS_PARAMETERS,
                                                     prolate_spheroid_main_dimensions)
from more_dynamics.plugins.shared.payload_io import OPEN_PARAMETERS, PayloadBuilder, frozen_block, payload_name


class ProlateSpheroidMainDimensions(HullForm):
    PARAMETERS = [
        ParameterDescription("length", 1.6),  # remus100.m:131 L_auv, passed as L at :221
        ParameterDescription("diameter", 0.19),  # remus100.m:132 D_auv, passed as B at :221
        OPEN_PARAMETERS,
    ]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        self._model = frozen_block(context, prolate_spheroid_main_dimensions(),
                                   PROLATE_SPHEROID_MAIN_DIMENSIONS_PARAMETERS)

    def graph(self) -> HullForm.CasadyPayload:
        if self._model is None:
            raise RuntimeError("ProlateSpheroidMainDimensions must be initialized before graph()")
        io = PayloadBuilder(HullForm.CasadyPayload())
        out = io.call(self._model)
        for name in self._model.name_out():
            io.output(payload_name(name), out[name])
        return io.payload()
