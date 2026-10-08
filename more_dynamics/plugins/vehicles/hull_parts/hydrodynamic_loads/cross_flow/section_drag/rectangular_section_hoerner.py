"""2-D drag of a rectangular section from Hoerner's table, Cd_2D(B / 2T)
(Hoerner.m 47-51), the section law of ``CrossFlowStrip`` on a pontoon.

Author:    Enio Krizman
Date:      2026-10-08
"""
from rpp_plugin_types.more_dynamics import SectionDragModel
from rpp_py.context import ComponentContext

from more_dynamics.models.vehicles.hull_parts.hydrodynamic_loads.cross_flow.section_drag.rectangular_section_hoerner import rectangular_section_hoerner
from more_dynamics.plugins.shared.payload_io import OPEN_PARAMETERS, PayloadBuilder


class RectangularSectionHoerner(SectionDragModel):
    PARAMETERS = [OPEN_PARAMETERS]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        self._model = rectangular_section_hoerner()

    def graph(self) -> SectionDragModel.CasadyPayload:
        if self._model is None:
            raise RuntimeError("RectangularSectionHoerner must be initialized before graph()")
        io = PayloadBuilder(SectionDragModel.CasadyPayload())
        out = io.call(self._model, rename={"beam": "section_beam"})
        io.output("section_drag_coefficient", out["section_drag_coefficient"])
        return io.payload()
