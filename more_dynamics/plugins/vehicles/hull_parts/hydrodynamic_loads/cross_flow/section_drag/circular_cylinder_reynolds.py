"""2-D drag of a circular section, C_D(Re) kappa(L/D) (cylinderDrag.m
78-110), the section law of ``CrossFlowStrip`` on a torpedo hull.

Author:    Enio Krizman
Date:      2026-10-08
"""
from rpp_plugin_types.more_dynamics import SectionDragModel
from rpp_py.context import ComponentContext

from more_dynamics.models.vehicles.hull_parts.hydrodynamic_loads.cross_flow.section_drag.circular_cylinder_reynolds import circular_cylinder_reynolds
from more_dynamics.plugins.shared.payload_io import OPEN_PARAMETERS, PayloadBuilder


class CircularCylinderReynolds(SectionDragModel):
    PARAMETERS = [OPEN_PARAMETERS]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        self._model = circular_cylinder_reynolds()

    def graph(self) -> SectionDragModel.CasadyPayload:
        if self._model is None:
            raise RuntimeError("CircularCylinderReynolds must be initialized before graph()")
        io = PayloadBuilder(SectionDragModel.CasadyPayload())
        out = io.call(self._model, rename={"beam": "section_beam"})
        io.output("section_drag_coefficient", out["section_drag_coefficient"])
        return io.payload()
