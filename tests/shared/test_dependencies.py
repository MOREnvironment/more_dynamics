"""The library declares the packages its shared helpers import: the parameter
framework (``Parameter``, ``freeze``, ``unset``) and the CasADi transforms
come from ``more_transformations``."""
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def test_package_json_lists_more_transformations():
    package = json.loads((ROOT / "package.json").read_text())
    assert "more_transformations" in package["Dependencies"]


def test_setup_requires_more_transformations_with_a_lower_bound():
    setup = (ROOT / "setup.py").read_text()
    assert re.search(r'"more_transformations>=\d+\.\d+\.\d+"', setup)


def test_the_parameter_framework_is_importable_from_more_transformations():
    from more_transformations.more_casadi_transformations import MatrixTransforms, Parameter, freeze, is_unset, unset
    assert is_unset(unset()) and not is_unset(0.0)
    assert all((MatrixTransforms, Parameter, freeze))
