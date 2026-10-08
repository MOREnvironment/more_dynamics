"""Hull form of two pontoons: the primitives, unchanged (otter.m 92-93,
104-107). ``section_beam``: the beam of a cross-flow strip section, here one
pontoon's beam (otter.m 245 at cc07579: crossFlowDrag(L, B_pont, T, nu_r)).

Author:    Enio Krizman
Date:      2026-10-08
"""
from rpp_plugin_types.more_dynamics import HullForm
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_dynamics.models.vehicles.hull_parts.hull_form.hull_form import TWIN_PONTOONS_PARAMETERS, twin_pontoons
from more_dynamics.plugins.shared.payload_io import OPEN_PARAMETERS, PayloadBuilder, frozen_block, payload_name


class TwinPontoons(HullForm):
    PARAMETERS = [
        ParameterDescription("length", 2.0),  # otter.m 92-93, 104-107
        ParameterDescription("beam", 1.08),  # otter.m 92-93, 104-107
        ParameterDescription("pontoon_beam", 0.25),  # otter.m 104-107
        ParameterDescription("pontoon_lateral_offset", 0.395),  # otter.m 104-107
        ParameterDescription("pontoon_block_coefficient", 0.4),  # otter.m 104-107
        ParameterDescription("pontoon_waterplane_coefficient", 0.75),  # otter.m 104-107
        OPEN_PARAMETERS,
    ]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        self._model = frozen_block(context, twin_pontoons(), TWIN_PONTOONS_PARAMETERS)

    def graph(self) -> HullForm.CasadyPayload:
        if self._model is None:
            raise RuntimeError("TwinPontoons must be initialized before graph()")
        io = PayloadBuilder(HullForm.CasadyPayload())
        out = io.call(self._model)
        for name in self._model.name_out():
            io.output(payload_name(name), out[name])
        return io.payload()
