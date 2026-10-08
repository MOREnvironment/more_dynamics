"""Cross-flow drag by the strip integral (crossFlowDrag.m 54-69; Fossen 2011,
eqs. 6.91-6.92, p. 127) as Luka's ``HydrodynamicsModel``, its 2-D section
drag law taken from its own child slot ``section`` — a composite, like
``HullVessel`` and ``CascadeController2D``: the child's payload is read
through its descriptions and handed to the model-layer strip integral.

Author:    Enio Krizman
Date:      2026-10-08
"""
from rpp_plugin_types.more_dynamics import HydrodynamicsModel
from rpp_py.context import ComponentContext

from more_dynamics.models.vehicles.hull_parts.hydrodynamic_loads.cross_flow.cross_flow_strip_parts import cross_flow_strip_with_section
from more_dynamics.plugins.payload_io import MODEL_NAMES, PayloadBuilder, named_function
from rpp_py.parameter_description import ParameterDescription


class CrossFlowStrip(HydrodynamicsModel):
    COMPONENTS = {"section": "more_dynamics::SectionDragModel"}
    PARAMETERS = [ParameterDescription("open_parameters", [])]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        section = context.get_component("section")
        if section == []:
            raise ValueError("CrossFlowStrip: slot 'section' is empty")
        self._model = cross_flow_strip_with_section(
            named_function(section.graph(), "section_drag", rename=MODEL_NAMES))

    def graph(self) -> HydrodynamicsModel.CasadyPayload:
        if self._model is None:
            raise RuntimeError("CrossFlowStrip must be initialized before graph()")
        io = PayloadBuilder(HydrodynamicsModel.CasadyPayload())
        out = io.call(self._model)
        io.output("hydrodynamic_force", out["tau"], "signed generalized force added by the vehicle, BODY")
        return io.payload()
