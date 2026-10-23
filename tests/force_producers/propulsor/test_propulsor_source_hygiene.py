"""What the propulsor's own source must not contain: Luka's modules, vehicle names, model traces, local transforms.

Written 2026-10-10, before the parts exist (contract: ``propulsor_contract``).

The scans read the propulsor's library folders: ``models/propulsor`` (its
parts and skeleton) and the plugin ``plugins/force_producers/propulsor.py``.

* **Ours stands alone** (owner, 2026-10-08: *"leave his as is but we make it in
  ours as well because my plan is to delete his parts later"*): no import of
  the shared repository's jet-nozzle or thruster plugins, hull-vessel,
  linear-surface or cross-flow-surface modules. The one place a test may
  import his plugin is the live-reference test, only when it is present.
* **Generic**: code outside docstrings and comments (which the
  syntax tree drops) names no vehicle, craft, motor make or author-form
  (``remus``, ``otter``, ``grethe``, ``marie``, ``hugin``, ``lauv``, ``nps``,
  ``torqeedo``, ``prestero``, ``mss``); equation citations live in docstrings
  and comments.
* **No trace of how the code was made**: no tool or author other than the
  owner named, no job, review or decision ids, anywhere in the files including
  comments (the scan's word list is stored reversed, so this file carries none
  of them).
* **Transforms**: no local rotation matrix, skew matrix or Euler
  convention: the parts import them from ``more_transformations``.
* **Density**: no part and no plugin declares or reads a
  ``water_density`` parameter; it reaches the propulsor as an input.
* **No new plugin type**: the skeleton is a plugin of the shared
  repository's ``ForceProducer`` type, its parts are model functions chosen by
  selectors; no ``Propulsor*`` interface is added to ``plugin_types``.

Author:    Enio Krizman
Date:      2026-10-10
"""

import ast
import re
from pathlib import Path

import pytest

LIBRARY = Path(__file__).resolve().parents[3]
MODEL_DIR = LIBRARY / "more_dynamics" / "models" / "propulsor"
PLUGIN_FILE = LIBRARY / "more_dynamics" / "plugins" / "force_producers" / "propulsor.py"

HIS_MODULES = ("plugins.force_producers.jet_nozzle", "plugins.force_producers.thruster", "models.hull_vessel",
               "models.hydrodynamics", "models.hydrostatics", "plugins.vehicle_models.hull_vessel",
               "plugins.hydrodynamics.linear_surface", "plugins.hydrostatics.linear_surface")
HIS_CLASSES = ("JetNozzle", "HullVessel", "LinearSurfaceHydrostatics", "LinearSurfaceHydrodynamics")
VEHICLE_WORDS = ("remus", "otter", "grethe", "marie", "hugin", "lauv", "nps", "torqeedo", "mariner", "cybership",
                 "prestero", "mss", "marble", "jetski", "jet_ski")
TRACE_WORDS = tuple(w[::-1] for w in (
    'edualc',
    'xedoc',
    'tpg',
    'ciporhtna',
    'ianepo',
    'tpgtahc',
    'tolipoc',
    'yb detareneg',
    'derohtua-oc',
    'a ssap',
    'b ssap',
    'dnuor riaper',
    'ffodnah',
    'feirb',
    '21u',
    '84-a',
    '85-a',
    '86-a',
    '29-e',
    '57-e',
    '401-e',
    '111-e',
    'ssm-d',
    'ssm-f',
    ']mgm[',
    'sledom_cireneg_erom'))


def _sources():
    files = sorted(MODEL_DIR.rglob("*.py")) if MODEL_DIR.exists() else []
    if PLUGIN_FILE.exists():
        files.append(PLUGIN_FILE)
    if not files:
        pytest.fail(f"not ported yet: {MODEL_DIR} and {PLUGIN_FILE}")
    return files


def _code(path):
    """The syntax tree back to text: no comments, no docstrings."""
    tree = ast.parse(path.read_text())
    for node in ast.walk(tree):
        body = getattr(node, "body", None)
        if isinstance(body, list) and body and isinstance(body[0], ast.Expr) \
                and isinstance(getattr(body[0], "value", None), ast.Constant) and isinstance(body[0].value.value, str):
            body.pop(0)
    return ast.unparse(tree)


def test_the_propulsor_has_its_model_folder_and_its_plugin():
    assert MODEL_DIR.is_dir(), f"not ported yet: {MODEL_DIR}"
    assert PLUGIN_FILE.is_file(), f"not ported yet: {PLUGIN_FILE}"
    for slot in ("drive", "steering", "mounting", "inflow", "law", "housing"):
        assert (MODEL_DIR / slot).is_dir() or (MODEL_DIR / f"{slot}.py").is_file(), slot


def test_no_import_of_the_shared_repositorys_jet_hull_or_surface_modules():
    offenders = []
    for path in _sources():
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            names = []
            if isinstance(node, ast.Import):
                names = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom):
                names = [(node.module or "") + "." + a.name for a in node.names] + [node.module or ""]
            for name in names:
                if any(module in name for module in HIS_MODULES) or any(name.endswith(c) for c in HIS_CLASSES):
                    offenders.append((str(path.relative_to(LIBRARY)), name))
    assert not offenders, offenders


def test_no_vehicle_or_author_form_in_the_code():
    offenders = []
    for path in _sources():
        words = set(re.findall(r"[a-z]+", _code(path).lower()))
        tokens = set(re.findall(r"[a-z_]+", _code(path).lower()))
        for word in VEHICLE_WORDS:
            if word in words or any(t.startswith(word + "_") or t.endswith("_" + word) for t in tokens):
                offenders.append((str(path.relative_to(LIBRARY)), word))
    assert not offenders, offenders


def test_no_ai_trace_and_no_job_ids_anywhere_in_the_files():
    offenders = []
    for path in _sources():
        text = path.read_text().lower()
        for word in TRACE_WORDS:
            if re.search(r"(?<![a-z0-9])" + re.escape(word) + r"(?![a-z0-9])", text):
                offenders.append((str(path.relative_to(LIBRARY)), word))
    assert not offenders, offenders


def test_every_file_ends_its_docstring_with_the_owner_as_author():
    for path in _sources():
        if path.name == "__init__.py" and not path.read_text().strip():
            continue
        doc = ast.get_docstring(ast.parse(path.read_text())) or ""
        assert "Author:    Enio Krizman" in doc, str(path.relative_to(LIBRARY))
        assert re.search(r"Date:\s+20\d\d-\d\d-\d\d", doc), str(path.relative_to(LIBRARY))


def test_no_local_rotation_skew_or_euler_convention():
    offenders = []
    for path in _sources():
        code = _code(path)
        if re.search(r"def\s+(skew|rot[xyz]|rzyx|rotation|euler|cross_matrix|smtrx)\w*\(", code, re.IGNORECASE) \
                or re.search(r"\[\s*\[?\s*0(\.0)?\s*,\s*-\s*\w+\[2\]\s*,", code):
            offenders.append(str(path.relative_to(LIBRARY)))
    assert not offenders, offenders


def test_no_part_and_no_plugin_declares_or_reads_a_density_parameter():
    offenders = []
    for path in _sources():
        code = _code(path)
        if re.search(r"Parameter\(\s*['\"]water_density['\"]", code) or \
                re.search(r"ParameterDescription\(\s*['\"]water_density['\"]", code) or \
                re.search(r"get_parameter\(\s*['\"]water_density['\"]", code):
            offenders.append(str(path.relative_to(LIBRARY)))
    assert not offenders, offenders


def test_no_propulsor_plugin_type_is_added():
    types = (LIBRARY / "more_dynamics" / "plugin_types").glob("*.capnp")
    text = "\n".join(p.read_text() for p in types)
    assert not re.search(r"interface\s+Propulsor\w*", text)
    assert not re.search(r"plugin\(\"Propulsor\w+\"\)", text)
