"""Checks of the frozen reference data of the fin parts
(``tests/data/force_producers/fin_parts/``).

Written 2026-10-08. These tests need no part and run now:

* every frozen file has the sha256 its ``SOURCE.md`` row names;
* the MATLAB file is what its generator says: ``remus100.m``'s propeller
  terms are zero, its ``tau`` is the fins' share (lines 249-254 with the
  propeller at rest), its constants are constant, and its rows reach the
  cases the gates need (zero, astern and sideways speed; angles beyond the
  limit);
* every number typed from Prestero 2001 names its table and page and equals
  its printed text;
* the generator finds MSS only through ``MSS_DIR`` and writes only under
  ``OUT_DIR`` (nothing relative to one machine);
* the cited lines of ``remus100.m`` are unchanged at the pinned revision
  (skipped when ``MSS_DIR`` is not set).

References
----------
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2
    with the fixes of 2026-10-07. https://github.com/cybergalactic/MSS, MIT
    licence, revision ``cc07579``: ``CRAFT/AUV/models/remus100.m`` 98-254.
[Prestero 2001] Prestero, T. (2001). *Verification of a six-degree of
    freedom simulation model for the REMUS autonomous underwater vehicle*.
    MSc thesis, MIT/WHOI. App. A, Tables A.1, A.5, pp. 102-103; App. C,
    Table C.10, p. 111.

Author:    Enio Krizman
Date:      2026-10-08
"""

import re

import numpy as np
import pytest

from fin_parts_contract import (
    DATA_DIR,
    G1_TOLERANCE,
    MSS_FINS_CSV,
    PRESTERO_JSON,
    external,
    mss_columns,
    mss_constant,
    mss_fins,
    prestero,
    sha256,
)

GENERATOR = "generate_remus100_fins_mss.m"
FROZEN = (MSS_FINS_CSV, GENERATOR, PRESTERO_JSON)
REMUS = "CRAFT/AUV/models/remus100.m"
CITED_LINES = {
    98: "rho = 1026;",
    109: "delta_max = deg2rad(20);",
    113: "delta_r = sat(ui(1), delta_max);",
    114: "delta_s = sat(ui(2), delta_max);",
    179: "S_fin = 0.00665;",
    182: "CL_delta_r = 0.5;",
    183: "A_r = 2 * S_fin;",
    184: "x_r = -a;",
    187: "CL_delta_s = 0.7;",
    188: "A_s = 2 * S_fin;",
    189: "x_s = -a;",
    234: "U_rh = sqrt( nu_r(1)^2 + nu_r(2)^2 );",
    235: "U_rv = sqrt( nu_r(1)^2 + nu_r(3)^2 );",
    238: "X_r = -0.5 * rho * U_rh^2 * A_r * CL_delta_r * delta_r^2;",
    239: "X_s = -0.5 * rho * U_rv^2 * A_s * CL_delta_s * delta_s^2;",
    242: "Y_r = -0.5 * rho * U_rh^2 * A_r * CL_delta_r * delta_r;",
    245: "Z_s = -0.5 * rho * U_rv^2 * A_s * CL_delta_s * delta_s;",
    249: "tau(1) = (1-t_prop) * X_prop + X_r + X_s;",
    250: "tau(2) = Y_r;",
    251: "tau(3) = Z_s;",
    252: "tau(4) = K_prop / 10;",
    253: "tau(5) = -x_s * Z_s;",
    254: "tau(6) = x_r * Y_r;",
}


def _source_hashes():
    text = (DATA_DIR / "SOURCE.md").read_text()
    hashes = {}
    for line in text.splitlines():
        found = re.search(r"\b([0-9a-f]{64})\b", line)
        name = re.search(r"`([^`]+\.(?:csv|m|json))`", line)
        if found and name:
            hashes[name.group(1)] = found.group(1)
    return hashes


@pytest.mark.parametrize("name", FROZEN)
def test_frozen_file_has_the_hash_its_source_row_names(name):
    assert _source_hashes().get(name) == sha256(DATA_DIR / name), name


def test_matlab_file_is_the_fins_share_of_remus100_tau():
    ref = mss_fins()
    assert len(ref["case_id"]) == 1024
    assert np.all(ref["X_prop"] == 0.0) and np.all(ref["K_prop"] == 0.0) and np.all(ref["ui3"] == 0.0)
    x_r, x_s = mss_constant(ref, "x_r"), mss_constant(ref, "x_s")
    fins = np.column_stack([ref["X_r"] + ref["X_s"], ref["Y_r"], ref["Z_s"], np.zeros_like(ref["Y_r"]),
                            -x_s * ref["Z_s"], x_r * ref["Y_r"]])
    assert np.max(np.abs(fins - mss_columns(ref, "tau", 6))) <= G1_TOLERANCE


def test_matlab_file_constants_are_remus100s():
    ref = mss_fins()
    for name in ("rho", "delta_max", "A_r", "A_s", "CL_delta_r", "CL_delta_s", "x_r", "x_s"):
        mss_constant(ref, name)
    assert mss_constant(ref, "rho") == 1026.0 and mss_constant(ref, "A_r") == mss_constant(ref, "A_s")


def test_matlab_file_reaches_the_cases_the_gates_need():
    ref = mss_fins()
    nu_r = mss_columns(ref, "nu_r", 6)
    limit = mss_constant(ref, "delta_max")
    assert np.any((ref["U_rh"] == 0.0) & (ref["U_rv"] == 0.0))            # at rest
    assert np.any(nu_r[:, 0] < 0.0)                                        # astern
    assert np.any((nu_r[:, 0] == 0.0) & (nu_r[:, 1] != 0.0))               # sideways only
    assert np.any(np.abs(ref["ui1"]) > limit) and np.any(np.abs(ref["ui2"]) > limit)


def test_prestero_numbers_name_their_page_and_equal_their_print():
    data = prestero()
    for name, entry in data["parameters"].items():
        assert re.search(r"Table [AC]\.\d+, p\. \d+", entry["where"]), name
        assert float(entry["printed"]) == entry["value"], name
    assert re.search(r"Table C\.10, p\. 111", data["control_fin_coefficients"]["where"])


def test_generator_finds_mss_only_through_its_variables():
    text = (DATA_DIR / GENERATOR).read_text()
    assert "getenv('MSS_DIR')" in text and "getenv('OUT_DIR')" in text
    assert not re.search(r"/home/|/Users/|[A-Z]:\\\\", text)


def test_cited_remus100_lines_are_unchanged():
    path = external("MSS_DIR", REMUS)
    lines = path.read_text().splitlines()
    for number, text in CITED_LINES.items():
        assert lines[number - 1].strip().startswith(text), (number, lines[number - 1])
