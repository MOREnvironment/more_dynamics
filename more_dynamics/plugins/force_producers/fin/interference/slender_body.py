"""Slender-body fin-body interference factors (Pitts et al. 1957): reserved.
The part refuses to initialise until the source is read in full and its
form and parameters are fixed.

Author:    Enio Krizman
Date:      2026-10-08
"""
from rpp_plugin_types.more_dynamics import FinInterference
from rpp_py.context import ComponentContext


class FinInterferenceSlenderBody(FinInterference):
    PARAMETERS = []

    def initialize(self, context: ComponentContext) -> None:
        raise NotImplementedError("FinInterferenceSlenderBody: waiting on source: Pitts 1957 interference factors")

    def graph(self) -> FinInterference.CasadyPayload:
        raise NotImplementedError("FinInterferenceSlenderBody: waiting on source: Pitts 1957 interference factors")
