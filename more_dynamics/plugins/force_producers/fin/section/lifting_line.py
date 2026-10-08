"""Lifting-line fin section with induced drag: reserved. The part refuses to
initialise until the source is read in full and its form and parameters are
fixed.

Author:    Enio Krizman
Date:      2026-10-08
"""
from rpp_plugin_types.more_dynamics import FinSection
from rpp_py.context import ComponentContext


class FinSectionLiftingLine(FinSection):
    PARAMETERS = []

    def initialize(self, context: ComponentContext) -> None:
        raise NotImplementedError("FinSectionLiftingLine: waiting on source: lifting-line induced drag")

    def graph(self) -> FinSection.CasadyPayload:
        raise NotImplementedError("FinSectionLiftingLine: waiting on source: lifting-line induced drag")
