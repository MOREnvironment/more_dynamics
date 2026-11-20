"""Stand-ins for the generated description types of ``rpp_schema.more_dynamics``.

``PayloadBuilder`` builds ``IODescription`` and ``StateDescription`` entries
from the types rpp generates when the library is registered
(``rpp library register ./more_dynamics --link``). Where they are not on the
path, minimal stand-ins with the same constructor signatures are put in
``sys.modules`` for the test, so the shared-helper tests do not depend on the
registration. Where the generated types exist they are used.
"""
import sys
import types

import pytest


class _IODescription:
    def __init__(self, size, min=None, max=None, name="", description=""):
        self.size, self.min, self.max = size, list(min or []), list(max or [])
        self.name, self.description = name, description


class _StateDescription:
    def __init__(self, size, ic=None, name="", description=""):
        self.size, self.ic, self.name, self.description = size, list(ic or []), name, description


@pytest.fixture(autouse=True)
def generated_description_types(monkeypatch):
    try:
        import rpp_schema.more_dynamics.IODescription  # noqa: F401
        import rpp_schema.more_dynamics.StateDescription  # noqa: F401
        return
    except ImportError:
        pass
    modules = {"rpp_schema": types.ModuleType("rpp_schema"),
               "rpp_schema.more_dynamics": types.ModuleType("rpp_schema.more_dynamics"),
               "rpp_schema.more_dynamics.IODescription": types.ModuleType("IODescription"),
               "rpp_schema.more_dynamics.StateDescription": types.ModuleType("StateDescription")}
    modules["rpp_schema.more_dynamics.IODescription"].IODescription = _IODescription
    modules["rpp_schema.more_dynamics.StateDescription"].StateDescription = _StateDescription
    for name, module in modules.items():
        monkeypatch.setitem(sys.modules, name, module)
