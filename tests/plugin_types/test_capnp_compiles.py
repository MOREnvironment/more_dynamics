"""Check the registered part types against the Cap'n Proto compiler.

Author:    Enio Krizman
Date:      2026-10-08
"""

import json
import os
from pathlib import Path
import re
import shutil
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[2]
TYPES = ROOT / "more_dynamics" / "plugin_types"
SCHEMA = "vehicle_parts.capnp"
INTERFACES = (
    "SiteModel", "CurrentModel", "HullForm", "RigidBodyModel", "AddedMassModel", "AddedMassCoriolisModel",
    "SectionDragModel", "ActuatorServo", "FinInflow", "FinFlowAngle", "FinInterference", "FinSection",
)


def test_declared_interfaces_and_registry_entries():
    registry = json.loads((ROOT / "plugins.json").read_text())
    paths = {entry["Path"] for entry in registry["PluginTypes"]}
    text = (TYPES / SCHEMA).read_text()
    declared = re.findall(r"^interface\s+(\w+)\s+\$Anot\.plugin\(\"(\w+)\"\)", text, re.M)
    assert [name for name, _ in declared] == [plugin for _, plugin in declared]
    assert set(name for name, _ in declared) == set(INTERFACES)
    for name in INTERFACES:
        assert re.search(rf"interface\s+{name}\b.*\n\s*graph @0 \(\) -> \(graph :Vm\.CasadyPayload\);", text), name
    assert f"more_dynamics/plugin_types/{SCHEMA}" in paths
    assert "more_dynamics/plugin_types/force_producer_parts.capnp" not in paths
    assert "BlockPayload" not in text and "Wire" not in text  # no payload beside Luka's CasadyPayload


def test_capnp_compiles(tmp_path):
    compiler = shutil.which("capnp")
    if compiler is None:
        pytest.skip("capnp compiler is not on PATH")
    common_root = os.environ.get("RPP_COMMON_CAPNP_DIR")
    if not common_root:
        pytest.skip("RPP_COMMON_CAPNP_DIR is not set; it must contain rpp_common/anot.capnp and msgs.capnp")
    common = Path(common_root) / "rpp_common"
    for name in ("anot.capnp", "msgs.capnp"):
        if not (common / name).is_file():
            pytest.skip(f"RPP_COMMON_CAPNP_DIR has no rpp_common/{name}")

    (tmp_path / "rpp_common").mkdir()
    for name in ("anot.capnp", "msgs.capnp"):
        shutil.copy2(common / name, tmp_path / "rpp_common" / name)
    for name in (SCHEMA, "vehicle_model.capnp"):
        shutil.copy2(TYPES / name, tmp_path / name)
    # the schema imports "more_dynamics/vehicle_model.capnp": the library folder is reachable under its own name,
    # as in rpp's registry
    (tmp_path / "more_dynamics").symlink_to(tmp_path, target_is_directory=True)
    environment = dict(os.environ, PWD=str(tmp_path))
    result = subprocess.run(
        [compiler, "compile", "-o-", SCHEMA], cwd=tmp_path,
        env=environment, capture_output=True, timeout=60,
    )
    assert result.returncode == 0, result.stderr.decode(errors="replace")
    assert result.stdout
