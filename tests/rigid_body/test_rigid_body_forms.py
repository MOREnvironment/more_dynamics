"""Gate tests for the rigid-body + added-mass block, **all forms**.

The rules they encode: the block is generic, with no
vehicle name or number in its code; the newest MSS is the reference, not the
truth, and every deviation from it is stated; there is no
``legacy_otter_inertia`` flag (its test is in ``test_rigid_body_block.py``);
skew and H come only from ``more_transformations``. The catamaran path
(``hull_with_payload``) keeps its gates in ``test_rigid_body_block.py``; this
file adds the other forms, the cross-form and physical-sign tests, the
transforms-import test and the generic test.

Gates (the test names carry them): G1 the block against MATLAB running MSS
(frozen CSVs), 1e-9 absolute; G2 the block against an independent transcription
or its own frozen form, 1e-10 absolute; G4 a perturbed model must differ
from the reference by more than 10x the G1 tolerance; G5 a printed number
reproduced at its printed precision.

Contract of the block (all CasADi, no numpy in the block)
---------------------------------------------------------
``more_dynamics.models.rigid_body`` exports:

* ``rigid_body_parameters(mass_properties="hull_with_payload")``: the declared
  parameters of one form (name, shape, SI unit, meaning, admissible range).
* ``rigid_body_casadi(*, mass_properties="hull_with_payload", coriolis="co",
  stabilize_added_mass_coriolis=False) -> ca.Function``. A form is chosen by
  these keywords, never by a vehicle name; unknown ``mass_properties`` or
  ``coriolis`` -> ``ValueError``. Inputs ``nu``, ``nu_r`` (6x1) and one input
  per declared parameter, by name; outputs by name: ``M``, ``C_RB``, ``C_A``,
  ``M_RB``, ``M_A`` (6x6), ``mass``, ``center_of_gravity`` (3x1, CO -> CG),
  ``inertia`` (3x3, about the CG), plus ``displaced_volume`` and
  ``wetted_surface`` (``L B + 2 T B``, used by the surge damping block) for
  ``"displacement_hull"`` and ``lamb_k_factors`` (3x1) for ``"spheroid"``.
* ``check_rigid_body_values(values, *, mass_properties)``: the numbers of one
  form checked where they enter (missing or foreign name, shape, finite,
  range, prolate spheroid, ``M`` symmetric positive definite) ->
  ``ValueError`` naming the parameter; returns ``ca.DM`` values.

  ``mass_properties`` (how m, I about the CG, r_g and M_A are obtained):

  - ``"hull_with_payload"`` (default; the catamaran path): ``length, beam,
    water_density, hull_mass, payload_mass, hull_center_of_gravity,
    payload_position, added_mass_coefficients, radii_of_gyration``.
  - ``"displacement_hull"``: ``water_density,
    length, beam, draft, block_coefficient, radii_of_gyration,
    hull_center_of_gravity, added_mass_coefficients``. ``m = rho Cb L B T``,
    ``I = m diag((R_s * [B, L, L])^2)`` about the CG, M_A by scaled derivatives
    ``c * [A11(m, L, rho), m, m, I11, I22, I33]`` (the ``otter.m`` 152-157
    pattern: the rotational terms scale the inertia **about the CG**).
  - ``"spheroid"`` (MSS ``spheroid.m`` + ``imlay61.m``):
    ``semi_major_axis, semi_minor_axis, body_density, water_density,
    roll_added_inertia_ratio, body_center_of_gravity``. Body: ``m =
    body_density 4/3 pi a b^2``, ``I = diag(2/5 m b^2, 1/5 m (a^2+b^2),
    1/5 m (a^2+b^2))``. Added mass (Lamb k-factors) from the displaced fluid,
    ``m_f`` and ``I_f`` with ``water_density``: ``diag(k1 m_f, k2 m_f, k2 m_f,
    r44 I_f,x, k' I_f,y, k' I_f,y)``. Two densities because MSS uses 1025 in
    ``spheroid.m`` (line 35) and 1026 in ``imlay61.m`` (line 31).
  - ``"explicit"``: ``body_mass,
    principal_inertia`` (3x1, the principal moments about the CG: a
    non-diagonal inertia cannot be given), ``body_center_of_gravity``, ``origin`` (``r_g =
    body_center_of_gravity - origin``),
    ``added_mass_derivatives`` (``[X_du, Y_dv, Z_dw, K_dp, M_dq, N_dr]``,
    ``M_A = -diag``).

  ``coriolis`` (the C_RB parametrisation; all give the same ``C_RB(nu) nu``):

  - ``"co"`` (default; MSS ``rbody.m`` / ``spheroid.m``):
    ``H^T diag(m S(w), -S(I w)) H``.
  - ``"book"``: ``[[m S(w), -m S(w) S(r_g)],
    [m S(r_g) S(w), -S(I_O w)]]``, ``I_O = I - m S(r_g)^2``. Differs from
    ``"co"`` only in the moment block, by ``-m |r_g|^2 S(w)``.
  - ``"lagrangian"`` (MSS ``m2c.m``):
    ``m2c(M_RB, nu)``, depends on the linear velocity too.

  ``C_A(nu_r)`` is always ``m2c(M_A, nu_r)`` (MSS ``m2c.m`` with its
  symmetrisation). ``stabilize_added_mass_coriolis=True``
  zeroes the eight entries of MSS ``remus100.m`` lines 207-210.
* 3-DOF forms (a **reduction of the 6-DOF block**; ``planar.py``):
  ``PLANAR_DOFS == (0, 1, 5)``; ``planar_added_mass_matrix(X_du, Y_dv, Y_dr,
  N_dr)`` (3x3, ``-[[X_du,0,0],[0,Y_dv,Y_dr],[0,Y_dr,N_dr]]``);
  ``planar_coriolis_casadi() -> ca.Function`` named ``"planar_coriolis"``,
  inputs ``["nu", "mass_matrix"]`` (3x1 ``[u, v, r]``, 3x3), output ``["C"]``
  (3x3), MSS ``m2c.m`` lines 50-54, with the declaration
  ``planar_coriolis_parameters()``; ``planar_casadi(*, mass_properties,
  coriolis, stabilize_added_mass_coriolis) -> ca.Function``, inputs ``nu``,
  ``nu_r`` (3x1) and the declared parameters, outputs ``["M", "C_RB",
  "C_A"]`` (3x3), the surge, sway, yaw reduction of ``rigid_body_casadi`` at
  ``[u, v, 0, 0, 0, r]``.
* Transforms from ``more_transformations.more_casadi_transformations`` only
  (``skew``/``H_matrix``); no numpy in the block.

Every test calls the block with numbers, ``rigid_body_casadi(...)(nu=...,
nu_r=..., **values)``, and reads the named outputs; ``_build`` returns them
(the state-independent ones as a namespace with the attribute names the
assertions use) together with the block and its checked values.

MSS (the checkout named by ``MSS_DIR``; read at ``72656d1``, release 2.0.2;
the code of every file below is unchanged since ``99bf0b3``): ``LIBRARY/modeling/rbody.m`` 31-45, ``m2c.m`` 33-56,
``spheroid.m`` 32-52, ``imlay61.m`` 30-62, ``addedMassSurge.m`` 32-34,
``LIBRARY/kinematics/Hmtrx.m`` 16-18, ``CRAFT/AUV/models/remus100.m`` 131-137,
199-210, ``CRAFT/USV/models/otter.m`` 152-159.

Frozen references (``tests/rigid_body/data/SOURCE.md``, section of 2026-10-06):
``rigid_body_rbody_mss_current.csv``, ``rigid_body_spheroid_mss_current.csv``,
``rigid_body_hull_mss_current.csv``, ``rigid_body_m2c_3dof_mss_current.csv``
(MATLAB R2026a running MSS), ``rigid_body_m2c_3dof_coupled_mss_cc07579.csv``
(``m2c.m`` 3-DOF on fully coupled mass matrices, MSS ``cc07579``), ``spheroid_matlab_reference_mss_current.csv``
(columns of the 2026-10-05 ``remus100.m`` CSV).

The physical tests below are written from Newton-Euler mechanics in the test
itself (CG velocity, force and moment about the CG moved to the CO: Fossen,
T. I. (2011), *Handbook of Marine Craft Hydrodynamics and Motion Control*,
Wiley, eqs. 3.14-3.18, pp. 47-48, and 3.33-3.40, pp. 50-51), not from model
code.

Outside this repo (nothing relative to one machine): ``MSS_DIR`` (MSS checkout) for the
pinned-line test; unset -> it skips and names the variable. Every other test
needs only ``tests/rigid_body/data/``.
"""

import ast
import contextlib
import importlib
import inspect
import io
import os
import subprocess
import sys
import warnings
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

DATA_DIR = Path(__file__).resolve().parent / "data"
# Outside this repo (nothing relative to one machine): found only through these variables,
# default not set; a test that needs one skips with a message naming it. The
# frozen CSVs in this folder's data/ are the only hard dependency.
MSS_DIR_VARIABLE = "MSS_DIR"                                  # MSS checkout root
CONTRACT_MODULE = "more_dynamics.models.rigid_body"

G1_TOLERANCE = 1e-9   # gate G1 (module docstring)
G2_TOLERANCE = 1e-10  # gate G2
G4_FACTOR = 10.0      # gate G4
SEED = 20261006
N_RANDOM_STATES = 1000
N_FROZEN = 50         # cases 1-50 of every reference CSV = nu of inputs.csv

# MSS lines whose constants are typed below (pinned by test_mss_cited_lines_are_unchanged).
MSS_LINES = {
    "CRAFT/AUV/models/remus100.m": {
        131: "L_auv = 1.6;",
        132: "D_auv = 0.19;",
        134: "a = 1.0096 * L_auv/2;    % Scaled spheroid semi-axes a and b to obtain m = 31.9 kg",
        135: "b = 1.0096 * D_auv/2;",
        136: "r44 = 0.3;",
        137: "r_bG = [ 0 0 0.02 ]';",
        199: "[MRB,CRB] = spheroid(a,b,nu(4:6),r_bG);",
        200: "[MA,CA] = imlay61(a, b, nu_r, r44);",
        207: "CA(5,3) = 0; CA(3,5) = 0;",
        208: "CA(5,1) = 0; CA(1,5) = 0;",
        209: "CA(6,1) = 0; CA(1,6) = 0;",
        210: "CA(6,2) = 0; CA(2,6) = 0;",
    },
    "LIBRARY/modeling/spheroid.m": {35: "rho = 1025;", 36: "m = 4/3 * pi * rho * a * b^2;"},
    "LIBRARY/modeling/imlay61.m": {31: "rho = 1026;", 59: "MA = diag([m*k1 m*k2 m*k2 MA_44 k_prime*Iy k_prime*Iy]);"},
    "LIBRARY/modeling/rbody.m": {35: "I_G   = m * diag([R44^2, R55^2, R66^2]);"},
    "LIBRARY/modeling/m2c.m": {35: "M = 0.5 * (M + M');", 54: "p = M * nu;"},
}
REMUS = {  # remus100.m lines 131-137
    "L_auv": 1.6, "D_auv": 0.19, "scale": 1.0096, "r44": 0.3, "r_bG": [0.0, 0.0, 0.02],
}
SPHEROID_M_DENSITY = 1025.0  # spheroid.m line 35
IMLAY61_DENSITY = 1026.0     # imlay61.m line 31
REMUS_PRINTED_MASS = 31.9    # remus100.m line 134 comment ("to obtain m = 31.9 kg")

# remus100.m 207-210, 0-based (row, col) pairs.
STABILIZED_PAIRS = ((4, 2), (4, 0), (5, 0), (5, 1))

VEHICLE_NAMES = ("otter", "remus", "catamaran", "grethe", "marie", "lauv", "hugin",
                 "mariner", "torqeedo", "cybership")


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def _contract():
    try:
        return importlib.import_module(CONTRACT_MODULE)
    except ModuleNotFoundError as exc:
        if exc.name == "casadi":
            pytest.fail(f"casadi is not installed: {exc}")
        pytest.fail(f"block not ported yet: {exc}")


def _build(mass_properties, coriolis="co", stabilize=False, **keywords):
    """The block of the chosen form with the checked numbers, and its
    state-independent named outputs under the attribute names the assertions
    below use (``rigid_body_mass_matrix`` = output ``M_RB``, ...;
    ``added_mass_derivatives`` = ``-diag(M_A)``)."""
    block = _contract()
    function = block.rigid_body_casadi(
        mass_properties=mass_properties,
        coriolis=coriolis,
        stabilize_added_mass_coriolis=stabilize,
    )
    values = block.check_rigid_body_values(keywords, mass_properties=mass_properties)
    out = function(nu=np.zeros(6), nu_r=np.zeros(6), **values)

    def optional(name, scalar=False):
        if name not in out:
            return None
        return float(out[name]) if scalar else np.array(out[name], dtype=float).ravel()

    constants = SimpleNamespace(
        mass_properties=mass_properties,
        coriolis=coriolis,
        stabilize_added_mass_coriolis=stabilize,
        mass=float(out["mass"]),
        center_of_gravity=np.array(out["center_of_gravity"], dtype=float).ravel(),
        inertia=np.array(out["inertia"], dtype=float),
        rigid_body_mass_matrix=np.array(out["M_RB"], dtype=float),
        added_mass_matrix=np.array(out["M_A"], dtype=float),
        total_mass_matrix=np.array(out["M"], dtype=float),
        added_mass_derivatives=-np.diag(np.array(out["M_A"], dtype=float)),
        lamb_k_factors=optional("lamb_k_factors"),
        displaced_volume=optional("displaced_volume", scalar=True),
        wetted_surface=optional("wetted_surface", scalar=True),
    )
    return constants, (function, values)


def _evaluate(function, nu, nu_r):
    block, values = function
    out = block(nu=np.asarray(nu, float), nu_r=np.asarray(nu_r, float), **values)
    return {name: np.array(out[name], dtype=float) for name in ("M", "C_RB", "C_A")}


def _planar_coriolis(mass_matrix):
    """``planar_coriolis_casadi()`` with ``mass_matrix`` frozen in: a
    function of ``nu = [u, v, r]`` only."""
    from more_transformations.more_casadi_transformations import freeze

    block = _contract()
    return freeze(block.planar_coriolis_casadi(), block.planar_coriolis_parameters(),
                  {"mass_matrix": np.asarray(mass_matrix, float)})


def _max_diff(a, b):
    return float(np.max(np.abs(np.asarray(a, float) - np.asarray(b, float))))


def _skew(v):
    x, y, z = v
    return np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])


class _Csv:
    """One reference CSV: echoed inputs, states and row-major matrices."""

    def __init__(self, name):
        path = DATA_DIR / name
        self.header = path.read_text().splitlines()[0].split(",")
        self.values = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)

    def column(self, name):
        return self.values[:, self.header.index(name)]

    def vector(self, prefix, n):
        return self.values[:, [self.header.index(f"{prefix}_{i}") for i in range(1, n + 1)]]

    def matrix(self, prefix, size=6):
        cols = [self.header.index(f"{prefix}_{i:02d}") for i in range(1, size * size + 1)]
        return self.values[:, cols].reshape(-1, size, size)  # row-major

    def rows(self, set_id):
        return np.flatnonzero(self.column("set") == set_id)


def _rbody():
    return _Csv("rigid_body_rbody_mss_current.csv")


def _spheroid():
    return _Csv("rigid_body_spheroid_mss_current.csv")


def _hull():
    return _Csv("rigid_body_hull_mss_current.csv")


def _dof3():
    return _Csv("rigid_body_m2c_3dof_mss_current.csv")


def _dof3_coupled():
    return _Csv("rigid_body_m2c_3dof_coupled_mss_cc07579.csv")


def _explicit_keywords(csv, row, added_mass_derivatives=(0.0,) * 6):
    """rbody.m line 35: I_G = m diag(R44^2, R55^2, R66^2)."""
    m = csv.column("m")[row]
    radii = np.array([csv.column(f"R{i}")[row] for i in (44, 55, 66)])
    return {
        "body_mass": m,
        "principal_inertia": m * radii**2,
        "body_center_of_gravity": csv.vector("r_bG", 3)[row],
        "origin": np.zeros(3),
        "added_mass_derivatives": np.asarray(added_mass_derivatives, float),
    }


def _spheroid_keywords(csv, row, body_density=SPHEROID_M_DENSITY, water_density=IMLAY61_DENSITY):
    return {
        "semi_major_axis": csv.column("a")[row],
        "semi_minor_axis": csv.column("b")[row],
        "body_density": body_density,
        "water_density": water_density,
        "roll_added_inertia_ratio": csv.column("r44")[row],
        "body_center_of_gravity": csv.vector("r_bG", 3)[row],
    }


def _hull_keywords(csv, row):
    return {
        "water_density": csv.column("rho")[row],
        "length": csv.column("L")[row],
        "beam": csv.column("B")[row],
        "draft": csv.column("T")[row],
        "block_coefficient": csv.column("Cb")[row],
        "radii_of_gyration": csv.vector("Rs", 3)[row],
        "hull_center_of_gravity": csv.vector("r_bG", 3)[row],
        "added_mass_coefficients": csv.vector("c", 6)[row],
    }


def _remus_keywords():
    """remus100.m lines 131-137, densities spheroid.m 35 / imlay61.m 31."""
    return {
        "semi_major_axis": REMUS["scale"] * REMUS["L_auv"] / 2,
        "semi_minor_axis": REMUS["scale"] * REMUS["D_auv"] / 2,
        "body_density": SPHEROID_M_DENSITY,
        "water_density": IMLAY61_DENSITY,
        "roll_added_inertia_ratio": REMUS["r44"],
        "body_center_of_gravity": REMUS["r_bG"],
    }


def _states():
    """The 50 frozen nu (inputs.csv) and 1000 seeded states (SEED)."""
    inputs = np.loadtxt(DATA_DIR / "inputs.csv", delimiter=",", skiprows=1, ndmin=2)
    rng = np.random.default_rng(SEED)
    seeded = rng.uniform(-3.0, 3.0, size=(N_RANDOM_STATES, 6))
    return np.vstack([inputs[:, 0:6], seeded])


def _newton_euler_wrench(m, inertia_cg, r_g, nu):
    """Coriolis-centripetal wrench about the CO, written from Newton-Euler.

    CG velocity v_g = v + w x r_g; force m w x v_g; moment about the CG
    w x (I_g w); moved to the CO by adding r_g x force. Needs no model code.
    """
    v, w = np.asarray(nu[:3], float), np.asarray(nu[3:], float)
    v_g = v + np.cross(w, r_g)
    force = m * np.cross(w, v_g)
    moment = np.cross(w, inertia_cg @ w) + np.cross(r_g, force)
    return np.concatenate([force, moment])


def _mss_dir():
    value = os.environ.get(MSS_DIR_VARIABLE)
    if not value:
        pytest.skip(f"{MSS_DIR_VARIABLE} is not set (path to the MSS checkout); MSS line checks skipped")
    return Path(value)


def _package_files():
    module = _contract()
    # the added-mass block sits in its own sibling folder; it is part of this unit's modules
    added_mass = importlib.import_module("more_dynamics.models.added_mass.added_mass")
    if hasattr(module, "__path__"):
        return sorted([Path(p) for base in module.__path__ for p in Path(base).glob("*.py")] + [Path(added_mass.__file__)])
    return [Path(module.__file__), Path(added_mass.__file__)]


# --------------------------------------------------------------------------
# Reference sanity (no port needed)
# --------------------------------------------------------------------------
def test_reference_csv_shapes_and_states():
    frozen_nu = _states()[:N_FROZEN]
    for csv, sets in ((_rbody(), 3), (_spheroid(), 2), (_hull(), 2), (_dof3(), 2)):
        assert len(csv.values) == 200 * sets
        for s in range(1, sets + 1):
            rows = csv.rows(s)
            assert np.array_equal(csv.column("case")[rows], np.arange(1, 201))
            nu = csv.vector("nu", 6)[rows] if "nu_1" in csv.header else None
            if nu is not None:
                np.testing.assert_array_equal(nu[:N_FROZEN], frozen_nu)
    dof3 = _dof3()
    np.testing.assert_array_equal(dof3.vector("nu3", 3)[:N_FROZEN], frozen_nu[:, [0, 1, 5]])
    remus = _Csv("spheroid_matlab_reference_mss_current.csv")
    assert remus.values.shape == (50, 6 + 4 * 36)


def test_reference_internal_consistency():
    """MATLAB columns agree with each other: m2c(MRB) and rbody's CRB give the
    same C nu; CA_stab is CA with exactly the remus100.m 207-210 entries zeroed."""
    rb = _rbody()
    nu = rb.vector("nu", 6)
    lhs = np.einsum("kij,kj->ki", rb.matrix("C_RB_m2c"), nu)
    rhs = np.einsum("kij,kj->ki", rb.matrix("C_RB"), nu)
    assert _max_diff(lhs, rhs) <= G1_TOLERANCE
    for csv in (_spheroid(), _hull()):
        ca_full, ca_stab = csv.matrix("C_A"), csv.matrix("C_A_stab").copy()
        for i, j in STABILIZED_PAIRS:
            assert np.all(ca_stab[:, i, j] == 0.0) and np.all(ca_stab[:, j, i] == 0.0)
            ca_stab[:, i, j], ca_stab[:, j, i] = ca_full[:, i, j], ca_full[:, j, i]
        np.testing.assert_array_equal(ca_stab, ca_full)


def test_mss_cited_lines_are_unchanged():
    mss = _mss_dir()
    for rel, lines in MSS_LINES.items():
        path = mss / rel
        if not path.exists():
            pytest.fail(f"{MSS_DIR_VARIABLE}={mss}: {rel} not found")
        text = path.read_text().splitlines()
        for number, expected in lines.items():
            assert text[number - 1].strip().startswith(expected.strip()), (rel, number, text[number - 1])


# --------------------------------------------------------------------------
# Contract
# --------------------------------------------------------------------------
def test_contract_selection_keywords():
    block = _contract()
    parameters = inspect.signature(block.rigid_body_casadi).parameters
    for name in ("mass_properties", "coriolis", "stabilize_added_mass_coriolis"):
        assert name in parameters, name
    assert parameters["mass_properties"].default == "hull_with_payload"
    assert parameters["coriolis"].default == "co"
    assert parameters["stabilize_added_mass_coriolis"].default is False
    assert "legacy_otter_inertia" not in parameters


def test_contract_rejects_unknown_forms_and_keywords():
    rb = _rbody()
    keywords = _explicit_keywords(rb, rb.rows(1)[0])
    block = _contract()
    with pytest.raises(ValueError):
        block.rigid_body_casadi(mass_properties="otter")
    with pytest.raises(ValueError):
        block.check_rigid_body_values(keywords, mass_properties="otter")
    with pytest.raises(ValueError):
        block.rigid_body_casadi(mass_properties="explicit", coriolis="book_corrected_typo")
    missing = dict(keywords)
    missing.pop("body_mass")
    with pytest.raises((TypeError, ValueError), match=r"missing.*'body_mass' \[kg\]"):
        block.check_rigid_body_values(missing, mass_properties="explicit")
    with pytest.raises((TypeError, ValueError), match=r"unknown.*'draft'"):
        block.check_rigid_body_values({**keywords, "draft": 0.3}, mass_properties="explicit")
    off_diagonal = dict(keywords)
    off_diagonal["principal_inertia"] = np.diag(keywords["principal_inertia"]) + np.array([[0, 1.0, 0], [1.0, 0, 0], [0, 0, 0]])
    # a 3x3 inertia cannot be given: the declared shape is the three principal
    # moments
    with pytest.raises(ValueError, match=r"'principal_inertia' \[kg\*m\^2\] must have shape \(3, 1\)"):
        block.check_rigid_body_values(off_diagonal, mass_properties="explicit")
    sph = _spheroid()
    swapped = _spheroid_keywords(sph, 0)
    swapped["semi_major_axis"], swapped["semi_minor_axis"] = swapped["semi_minor_axis"], swapped["semi_major_axis"]
    with pytest.raises(ValueError, match="'semi_major_axis'"):  # imlay61.m line 42
        block.check_rigid_body_values(swapped, mass_properties="spheroid")
    for name, bad in (("body_mass", -1.0), ("body_mass", float("nan")), ("principal_inertia", [1.0, 0.0, 1.0])):
        with pytest.raises(ValueError, match=f"'{name}'"):
            block.check_rigid_body_values({**keywords, name: bad}, mass_properties="explicit")


@pytest.mark.parametrize("coriolis", ["co", "book", "lagrangian"])
def test_contract_function_and_constants(coriolis):
    rb = _rbody()
    from more_transformations.more_casadi_transformations import freeze

    keywords = _explicit_keywords(rb, rb.rows(1)[0])
    constants, (block, _) = _build("explicit", coriolis, **keywords)
    declared = _contract().rigid_body_parameters("explicit")
    assert block.name_in() == ["nu", "nu_r", *[d.name for d in declared]]
    assert block.name_out() == list(_contract().rigid_body_outputs("explicit"))
    assert block.name_out() == ["M", "C_RB", "C_A", "M_RB", "M_A", "mass", "center_of_gravity", "inertia"]
    function = freeze(block, declared, keywords)  # the numbers frozen in: a function of the state only
    assert function.name_in() == ["nu", "nu_r"]
    assert function.name_out()[:3] == ["M", "C_RB", "C_A"]
    assert [function.size_in(i) for i in range(2)] == [(6, 1)] * 2
    assert [function.size_out(i) for i in range(3)] == [(6, 6)] * 3
    assert constants.mass_properties == "explicit"
    assert constants.coriolis == coriolis
    assert constants.stabilize_added_mass_coriolis is False
    for name in ("rigid_body_mass_matrix", "added_mass_matrix", "total_mass_matrix"):
        assert np.asarray(getattr(constants, name)).shape == (6, 6), name
    assert np.asarray(constants.inertia).shape == (3, 3)
    assert np.asarray(constants.center_of_gravity).shape == (3,)
    np.testing.assert_allclose(
        constants.total_mass_matrix, constants.rigid_body_mass_matrix + constants.added_mass_matrix,
        atol=G2_TOLERANCE, rtol=0.0,
    )
    np.testing.assert_allclose(
        constants.added_mass_derivatives, -np.diag(constants.added_mass_matrix), atol=G2_TOLERANCE, rtol=0.0
    )


# --------------------------------------------------------------------------
# G1 — block vs MATLAB running current MSS
# --------------------------------------------------------------------------
@pytest.mark.parametrize("set_id", [1, 2, 3])
@pytest.mark.parametrize("coriolis", ["co", "lagrangian"])
def test_G1_explicit_matches_rbody_and_m2c(set_id, coriolis):
    """rbody.m (MRB, CRB(nu2)) and m2c(MRB, nu), 200 states per set."""
    rb = _rbody()
    rows = rb.rows(set_id)
    constants, function = _build("explicit", coriolis, **_explicit_keywords(rb, rows[0]))
    assert _max_diff(constants.rigid_body_mass_matrix, rb.matrix("M_RB")[rows[0]]) <= G1_TOLERANCE
    expected = rb.matrix("C_RB" if coriolis == "co" else "C_RB_m2c")
    for k in rows:
        nu = rb.vector("nu", 6)[k]
        out = _evaluate(function, nu, nu)
        assert _max_diff(out["M"], rb.matrix("M_RB")[k]) <= G1_TOLERANCE, k
        assert _max_diff(out["C_RB"], expected[k]) <= G1_TOLERANCE, k


@pytest.mark.parametrize("set_id", [1, 2, 3])
def test_G1_book_form_equals_rbody_plus_documented_term(set_id):
    """No MSS function holds the book form: it equals rbody's CRB minus
    m |r_g|^2 S(w) in the moment block (Fossen 2011, eq. 3.30, p. 50)."""
    rb = _rbody()
    rows = rb.rows(set_id)
    keywords = _explicit_keywords(rb, rows[0])
    _, function = _build("explicit", "book", **keywords)
    m, r_g = keywords["body_mass"], keywords["body_center_of_gravity"]
    for k in rows:
        nu = rb.vector("nu", 6)[k]
        expected = rb.matrix("C_RB")[k].copy()
        expected[3:, 3:] -= m * float(r_g @ r_g) * _skew(nu[3:])
        assert _max_diff(_evaluate(function, nu, nu)["C_RB"], expected) <= G1_TOLERANCE, k


@pytest.mark.parametrize("set_id", [1, 2])
@pytest.mark.parametrize("stabilize", [False, True])
def test_G1_spheroid_matches_spheroid_and_imlay61(set_id, stabilize):
    """spheroid.m (rho 1025) + imlay61.m (rho 1026); set 2 calls imlay61 with
    nargin 3 (r44 = 0); stabilised C_A = remus100.m 207-210."""
    sph = _spheroid()
    rows = sph.rows(set_id)
    constants, function = _build("spheroid", "co", stabilize, **_spheroid_keywords(sph, rows[0]))
    assert _max_diff(constants.rigid_body_mass_matrix, sph.matrix("M_RB")[rows[0]]) <= G1_TOLERANCE
    assert _max_diff(constants.added_mass_matrix, sph.matrix("M_A")[rows[0]]) <= G1_TOLERANCE
    c_a = sph.matrix("C_A_stab" if stabilize else "C_A")
    for k in rows:
        nu = sph.vector("nu", 6)[k]
        out = _evaluate(function, nu, nu)
        assert _max_diff(out["M"], sph.matrix("M_RB")[k] + sph.matrix("M_A")[k]) <= G1_TOLERANCE, k
        assert _max_diff(out["C_RB"], sph.matrix("C_RB")[k]) <= G1_TOLERANCE, k
        assert _max_diff(out["C_A"], c_a[k]) <= G1_TOLERANCE, k


def test_G1_spheroid_matches_remus100_path():
    """The 2026-10-05 CSV (remus100.m as written, 50 cases): MRB, MA, CRB at nu, CA at nu_r
    stabilised (remus100.m 199-210)."""
    ref = _Csv("spheroid_matlab_reference_mss_current.csv")
    nu = _states()[:N_FROZEN]
    nu_r = ref.values[:, [ref.header.index(f"nu_r_{i:02d}") for i in range(1, 7)]]
    constants, function = _build("spheroid", "co", True, **_remus_keywords())
    for k in range(N_FROZEN):
        out = _evaluate(function, nu[k], nu_r[k])
        assert _max_diff(out["M"], ref.matrix("M_RB")[k] + ref.matrix("M_A")[k]) <= G1_TOLERANCE, k
        assert _max_diff(out["C_RB"], ref.matrix("C_RB")[k]) <= G1_TOLERANCE, k
        assert _max_diff(out["C_A"], ref.matrix("C_A")[k]) <= G1_TOLERANCE, k
    assert _max_diff(constants.rigid_body_mass_matrix, ref.matrix("M_RB")[0]) <= G1_TOLERANCE


@pytest.mark.parametrize("set_id", [1, 2])
@pytest.mark.parametrize("stabilize", [False, True])
def test_G1_displacement_hull_matches_rbody_and_added_mass(set_id, stabilize):
    """m = rho Cb L B T; rbody.m; addedMassSurge.m; otter.m 152-159 pattern; m2c.m."""
    hull = _hull()
    rows = hull.rows(set_id)
    constants, function = _build("displacement_hull", "co", stabilize, **_hull_keywords(hull, rows[0]))
    assert abs(constants.mass - hull.column("m")[rows[0]]) <= G1_TOLERANCE
    assert _max_diff(constants.added_mass_matrix, hull.matrix("M_A")[rows[0]]) <= G1_TOLERANCE
    assert abs(constants.added_mass_matrix[0, 0] - (-hull.column("c_1")[rows[0]]) * hull.column("A11")[rows[0]]) <= G1_TOLERANCE
    c_a = hull.matrix("C_A_stab" if stabilize else "C_A")
    for k in rows:
        nu = hull.vector("nu", 6)[k]
        out = _evaluate(function, nu, nu)
        assert _max_diff(out["M"], hull.matrix("M_RB")[k] + hull.matrix("M_A")[k]) <= G1_TOLERANCE, k
        assert _max_diff(out["C_RB"], hull.matrix("C_RB")[k]) <= G1_TOLERANCE, k
        assert _max_diff(out["C_A"], c_a[k]) <= G1_TOLERANCE, k


@pytest.mark.parametrize("set_id", [1, 2])
def test_G1_explicit_added_mass_matches_m2c(set_id):
    """Derivative set given explicitly (M_A_6dof / C_A_6dof): the hull CSV's
    MA and m2c(MA, nu)."""
    hull, rb = _hull(), _rbody()
    rows = hull.rows(set_id)
    derivatives = -np.diag(hull.matrix("M_A")[rows[0]])
    _, function = _build("explicit", "co", **_explicit_keywords(rb, rb.rows(1)[0], derivatives))
    for k in rows:
        nu = hull.vector("nu", 6)[k]
        assert _max_diff(_evaluate(function, nu, nu)["C_A"], hull.matrix("C_A")[k]) <= G1_TOLERANCE, k


@pytest.mark.parametrize("set_id", [1, 2])
def test_G1_planar_forms_match_m2c_3dof(set_id):
    """3-DOF (reduction of the 6-DOF block): m2c.m lines 50-54 on M = -[[Xu,0,0],[0,Yv,Yr],[0,Yr,Nr]]."""
    block = _contract()
    dof3 = _dof3()
    rows = dof3.rows(set_id)
    args = [dof3.column(name)[rows[0]] for name in ("X_du", "Y_dv", "Y_dr", "N_dr")]
    mass_matrix = np.array(block.planar_added_mass_matrix(*args), dtype=float)
    assert _max_diff(mass_matrix, dof3.matrix("M3", 3)[rows[0]]) <= G1_TOLERANCE
    unfrozen = block.planar_coriolis_casadi()
    assert unfrozen.name() == "planar_coriolis"
    assert unfrozen.name_in() == ["nu", "mass_matrix"] and unfrozen.size_in("mass_matrix") == (3, 3)
    function = _planar_coriolis(mass_matrix)
    assert function.name_in() == ["nu"] and function.name_out() == ["C"]
    assert function.size_in(0) == (3, 1) and function.size_out(0) == (3, 3)
    for k in rows:
        c = np.array(function(nu=dof3.vector("nu3", 3)[k])["C"], dtype=float)
        assert _max_diff(c, dof3.matrix("C3", 3)[k]) <= G1_TOLERANCE, k


@pytest.mark.parametrize("set_id", [1, 2])
def test_G1_planar_coriolis_matches_m2c_3dof_on_coupled_mass_matrices(set_id):
    """MSS ``cc07579``, ``m2c.m`` lines 52-57 (p = M nu) on two fully coupled
    symmetric M3, 201 states each: ours equals MATLAB, entries relative to
    max(1, |C|). Case 0 of set 1 is the counterexample [1, 2, 3]: C nu3 =
    [-153, 51, 17]. The branch before ``9aef3ba`` (``72656d1`` lines 52-54,
    transcribed) differs on these rows, so this file sees the surge coupling
    that the uncoupled file above cannot."""
    csv = _dof3_coupled()
    rows = csv.rows(set_id)
    m3 = csv.matrix("M3", 3)[rows[0]]
    assert np.all(np.abs(m3[np.triu_indices(3, 1)]) > 0.0)
    function = _planar_coriolis(m3)
    worst_old = 0.0
    for k in rows:
        nu3 = csv.vector("nu3", 3)[k]
        expected = csv.matrix("C3", 3)[k]
        got = np.array(function(nu=nu3)["C"], dtype=float)
        assert _max_diff(got, expected) <= G1_TOLERANCE * max(1.0, np.abs(expected).max()), k
        u, v, r = nu3
        old = np.array([[0.0, 0.0, -m3[1, 1] * v - m3[1, 2] * r],
                        [0.0, 0.0, m3[0, 0] * u],
                        [m3[1, 1] * v + m3[1, 2] * r, -m3[0, 0] * u, 0.0]])
        worst_old = max(worst_old, _max_diff(old, expected))
    assert worst_old > G4_FACTOR * G1_TOLERANCE, worst_old
    if set_id == 1:
        first = rows[csv.column("case")[rows] == 0][0]
        assert np.array_equal(csv.vector("nu3", 3)[first], [1.0, 2.0, 3.0])
        np.testing.assert_array_equal(csv.matrix("C3", 3)[first] @ [1.0, 2.0, 3.0], [-153.0, 51.0, 17.0])


# --------------------------------------------------------------------------
# G4 — a broken model must fail (one perturbation per form, vs the MATLAB CSV)
# --------------------------------------------------------------------------
def _worst(csv, rows, function, c_a_prefix="C_A", c_rb_prefix="C_RB", with_added_mass=True):
    worst = 0.0
    for k in rows:
        nu = csv.vector("nu", 6)[k]
        out = _evaluate(function, nu, nu)
        expected_m = csv.matrix("M_RB")[k] + (csv.matrix("M_A")[k] if with_added_mass else 0.0)
        worst = max(worst, _max_diff(out["M"], expected_m), _max_diff(out["C_RB"], csv.matrix(c_rb_prefix)[k]))
        if with_added_mass:
            worst = max(worst, _max_diff(out["C_A"], csv.matrix(c_a_prefix)[k]))
    return worst


def _g4_cases():
    return {
        "explicit_mass_plus_1_percent": ("explicit", lambda kw: {**kw, "body_mass": kw["body_mass"] * 1.01}),
        "explicit_cg_x_plus_1_cm": ("explicit", lambda kw: {**kw, "body_center_of_gravity": kw["body_center_of_gravity"] + [0.01, 0, 0]}),
        "spheroid_body_density_plus_1_percent": ("spheroid", lambda kw: {**kw, "body_density": kw["body_density"] * 1.01}),
        "spheroid_r44_times_zero": ("spheroid", lambda kw: {**kw, "roll_added_inertia_ratio": 0.0}),
        "spheroid_water_density_equals_body": ("spheroid", lambda kw: {**kw, "water_density": kw["body_density"]}),
        "hull_block_coefficient_plus_1_percent": ("displacement_hull", lambda kw: {**kw, "block_coefficient": kw["block_coefficient"] * 1.01}),
        "hull_cg_z_plus_1_cm": ("displacement_hull", lambda kw: {**kw, "hull_center_of_gravity": kw["hull_center_of_gravity"] + [0, 0, 0.01]}),
        "hull_added_mass_sway_times_zero": ("displacement_hull", lambda kw: {**kw, "added_mass_coefficients": kw["added_mass_coefficients"] * [1, 0, 1, 1, 1, 1]}),
    }


@pytest.mark.parametrize("name", sorted(_g4_cases()))
def test_G4_perturbed_form_is_detected(name):
    form, perturb = _g4_cases()[name]
    if form == "explicit":
        csv = _rbody()
        rows = csv.rows(1)
        keywords = _explicit_keywords(csv, rows[0])
        added = False
    elif form == "spheroid":
        csv = _spheroid()
        rows = csv.rows(1)
        keywords = _spheroid_keywords(csv, rows[0])
        added = True
    else:
        csv = _hull()
        rows = csv.rows(2)
        keywords = _hull_keywords(csv, rows[0])
        added = True
    _, control = _build(form, **keywords)
    assert _worst(csv, rows, control, with_added_mass=added) <= G1_TOLERANCE, "control"
    _, function = _build(form, **perturb(keywords))
    worst = _worst(csv, rows, function, with_added_mass=added)
    assert worst > G4_FACTOR * G1_TOLERANCE, (name, worst)


def test_G4_planar_perturbation_is_detected():
    block = _contract()
    dof3 = _dof3()
    rows = dof3.rows(2)
    args = [dof3.column(name)[rows[0]] for name in ("X_du", "Y_dv", "Y_dr", "N_dr")]
    for perturbed, expect_detect in ((args, False), ([args[0], args[1], 0.0, args[3]], True)):
        function = _planar_coriolis(block.planar_added_mass_matrix(*perturbed))
        worst = max(_max_diff(np.array(function(nu=dof3.vector("nu3", 3)[k])["C"], float), dof3.matrix("C3", 3)[k])
                    for k in rows)
        assert (worst > G4_FACTOR * G1_TOLERANCE) is expect_detect, (perturbed, worst)


# --------------------------------------------------------------------------
# 3-DOF as a reduction of the 6-DOF block: planar_casadi
# --------------------------------------------------------------------------
PLANAR_SEED = 20261007
N_PLANAR_STATES = 300


def _planar_forms():
    """Two mass-property forms whose CG couples the planar rows: ``explicit``
    on rbody set 1 (r_g = [0.12, -0.04, 0.21], so M_RB(1,6) = -m y_g != 0)
    and ``displacement_hull`` on hull set 2 (r_g = [0.15, 0, -0.10], sway-yaw
    coupling m x_g); each with the ``co`` and ``lagrangian`` C_RB, stabilised
    C_A and not."""
    rb, hull = _rbody(), _hull()
    forms = []
    for coriolis in ("co", "lagrangian"):
        for stabilize in (False, True):
            forms.append((f"explicit-1-{coriolis}-{'stab' if stabilize else 'full'}", "explicit", coriolis,
                          stabilize, _explicit_keywords(rb, rb.rows(1)[0], -np.diag(hull.matrix("M_A")[0]))))
            forms.append((f"hull-2-{coriolis}-{'stab' if stabilize else 'full'}", "displacement_hull", coriolis,
                          stabilize, _hull_keywords(hull, hull.rows(2)[0])))
    return forms


def _planar_states():
    """Seeded ``[u, v, r]`` and ``[u_r, v_r, r_r]``, uniform in [-3, 3]."""
    rng = np.random.default_rng(PLANAR_SEED)
    return rng.uniform(-3.0, 3.0, size=(N_PLANAR_STATES, 2, 3))


def _embed3(nu3):
    """``[u, v, r] -> [u, v, 0, 0, 0, r]``."""
    return np.array([nu3[0], nu3[1], 0.0, 0.0, 0.0, nu3[2]])


def _planar(form, coriolis, stabilize, full):
    """``planar_casadi`` built with the same keywords as the 6-DOF block
    ``full`` (the pair returned by ``_build``), with the same checked values."""
    function = _contract().planar_casadi(
        mass_properties=form, coriolis=coriolis, stabilize_added_mass_coriolis=stabilize
    )
    return function, full[1]


def _planar_worst(planar, full, index):
    """Largest difference of M, C_RB nu3 and C_A nu_r3 of ``planar`` from the
    surge, sway, yaw part of the 6-DOF ``full`` at the embedded states."""
    worst = 0.0
    for nu3, nu_r3 in _planar_states():
        nu6, nu_r6 = _embed3(nu3), _embed3(nu_r3)
        six = _evaluate(full, nu6, nu_r6)
        out = planar[0](nu=nu3, nu_r=nu_r3, **planar[1])
        three = {name: np.array(out[name], dtype=float) for name in ("M", "C_RB", "C_A")}
        worst = max(
            worst,
            _max_diff(three["M"], six["M"][np.ix_(index, index)]),
            _max_diff(three["C_RB"] @ nu3, (six["C_RB"] @ nu6)[index]),
            _max_diff(three["C_A"] @ nu_r3, (six["C_A"] @ nu_r6)[index]),
        )
    return worst


@pytest.mark.parametrize("label,form,coriolis,stabilize,keywords", _planar_forms(),
                         ids=[f[0] for f in _planar_forms()])
def test_G2_planar_casadi_is_the_reduction_of_the_6dof_block(label, form, coriolis, stabilize, keywords):
    """``planar_casadi(constants)`` against ``rigid_body_casadi(constants)`` at
    ``[u, v, 0, 0, 0, r]`` / ``[u_r, v_r, 0, 0, 0, r_r]``: M is rows and columns
    ``PLANAR_DOFS`` of the 6-DOF M; C_RB nu3 and C_A nu_r3 are the surge, sway,
    yaw rows of the 6-DOF forces (G2, absolute)."""
    block = _contract()
    assert tuple(block.PLANAR_DOFS) == (0, 1, 5)
    index = list(block.PLANAR_DOFS)
    constants, full = _build(form, coriolis, stabilize, **keywords)
    planar = _planar(form, coriolis, stabilize, full)
    declared = [d.name for d in block.rigid_body_parameters(form)]
    assert planar[0].name_in() == ["nu", "nu_r", *declared] and planar[0].name_out() == ["M", "C_RB", "C_A"]
    assert [planar[0].size_in(i) for i in range(2)] == [(3, 1)] * 2
    assert [planar[0].size_out(i) for i in range(3)] == [(3, 3)] * 3
    reduced = np.array(planar[0](nu=np.zeros(3), nu_r=np.zeros(3), **planar[1])["M"], dtype=float)
    np.testing.assert_array_equal(reduced,
                                  np.asarray(constants.total_mass_matrix)[np.ix_(index, index)])
    assert _planar_worst(planar, full, index) <= G2_TOLERANCE, label


@pytest.mark.parametrize("form", ["explicit", "displacement_hull"])
def test_G4_planar_casadi_perturbed_sway_coefficient_is_detected(form):
    """One planar coefficient changed (sway added-mass derivative x1.01 for
    ``explicit``, sway added-mass coefficient x1.01 for ``displacement_hull``):
    the reduction of the perturbed constants is compared with the unperturbed
    6-DOF block and must differ by more than G4_FACTOR x G1."""
    block = _contract()
    label, _, coriolis, stabilize, keywords = next(f for f in _planar_forms() if f[1] == form)
    index = list(block.PLANAR_DOFS)
    constants, full = _build(form, coriolis, stabilize, **keywords)
    assert _planar_worst(_planar(form, coriolis, stabilize, full), full, index) <= G2_TOLERANCE, "control"
    key = "added_mass_derivatives" if form == "explicit" else "added_mass_coefficients"
    perturbed = {**keywords, key: np.asarray(keywords[key], float) * [1, 1.01, 1, 1, 1, 1]}
    _, perturbed_full = _build(form, coriolis, stabilize, **perturbed)
    worst = _planar_worst(_planar(form, coriolis, stabilize, perturbed_full), full, index)
    assert worst > G4_FACTOR * G1_TOLERANCE, (label, worst)


def test_physics_planar_coriolis_keeps_surge_coupling_of_a_general_mass_matrix():
    """A generally coupled symmetric positive-definite M3 (seeded, every
    off-diagonal entry non-zero): ``planar_coriolis_casadi(M3)`` equals the
    Kirchhoff form written here from the planar momentum p = M3 [u, v, r]
    (no model code), C3 = [[0, 0, -p_y], [0, 0, p_x], [p_y, -p_x, 0]].

    MSS ``m2c.m``, transcribed on the same states:
    - from ``cc07579`` (lines 54-57, p = M nu; corrected in ``9aef3ba`` /
      ``a3406cf``): equal to ours, 1e-10 relative — MSS now confirms the
      reduction;
    - history, before ``9aef3ba`` (lines 52-54 at ``72656d1``): kept only
      M(1,1) u in p_x and M(2,2) v + M(2,3) r in p_y, and must differ, so the
      test cannot pass on that branch."""
    block = _contract()
    rng = np.random.default_rng(PLANAR_SEED)
    for k in range(50):
        factor = rng.uniform(-1.0, 1.0, size=(3, 3))
        m3 = factor @ factor.T + 3.0 * np.eye(3)
        assert np.all(np.abs(m3[np.triu_indices(3, 1)]) > 0.0)
        function = _planar_coriolis(m3)
        for nu3 in rng.uniform(-3.0, 3.0, size=(20, 3)):
            p_x, p_y, _ = m3 @ nu3
            kirchhoff = np.array([[0.0, 0.0, -p_y], [0.0, 0.0, p_x], [p_y, -p_x, 0.0]])
            got = np.array(function(nu=nu3)["C"], dtype=float)
            assert _max_diff(got, kirchhoff) <= G2_TOLERANCE * max(1.0, np.abs(kirchhoff).max()), k
            p = m3 @ nu3  # m2c.m 54 at cc07579
            corrected = np.array([[0.0, 0.0, -p[1]], [0.0, 0.0, p[0]], [p[1], -p[0], 0.0]])  # m2c.m 55-57
            assert _max_diff(got, corrected) <= G2_TOLERANCE * max(1.0, np.abs(corrected).max()), k
            u, v, r = nu3
            branch = np.array([[0.0, 0.0, -m3[1, 1] * v - m3[1, 2] * r],
                               [0.0, 0.0, m3[0, 0] * u],
                               [m3[1, 1] * v + m3[1, 2] * r, -m3[0, 0] * u, 0.0]])
            assert _max_diff(branch, kirchhoff) > G4_FACTOR * G1_TOLERANCE, k


# --------------------------------------------------------------------------
# G5 — printed number
# --------------------------------------------------------------------------
def test_G5_remus_spheroid_mass_printed_in_remus100():
    """remus100.m line 134: semi-axes scaled 'to obtain m = 31.9 kg' (rho 1025,
    spheroid.m line 35-36); printed precision 0.1 kg."""
    constants, _ = _build("spheroid", **_remus_keywords())
    assert round(float(constants.mass), 1) == REMUS_PRINTED_MASS


# --------------------------------------------------------------------------
# Cross-form consistency (block only)
# --------------------------------------------------------------------------
@pytest.mark.parametrize("set_id", [1, 2, 3])
def test_cross_form_coriolis_forms_give_the_same_force(set_id):
    """co, book and lagrangian: same M_RB and same C_RB(nu) nu; co - book =
    m |r_g|^2 S(w) in the moment block only; lagrangian is a different matrix
    whenever the linear velocity is non-zero."""
    rb = _rbody()
    keywords = _explicit_keywords(rb, rb.rows(set_id)[0])
    built = {c: _build("explicit", c, **keywords) for c in ("co", "book", "lagrangian")}
    m, r_g = keywords["body_mass"], keywords["body_center_of_gravity"]
    for c in ("book", "lagrangian"):
        assert _max_diff(built[c][0].rigid_body_mass_matrix, built["co"][0].rigid_body_mass_matrix) <= G2_TOLERANCE
    differs = False
    for k, nu in enumerate(_states()):
        c = {name: _evaluate(f, nu, nu)["C_RB"] for name, (_, f) in built.items()}
        assert _max_diff(c["book"] @ nu, c["co"] @ nu) <= G2_TOLERANCE, k
        assert _max_diff(c["lagrangian"] @ nu, c["co"] @ nu) <= G2_TOLERANCE, k
        term = np.zeros((6, 6))
        term[3:, 3:] = m * float(r_g @ r_g) * _skew(nu[3:])
        assert _max_diff(c["co"] - c["book"], term) <= G2_TOLERANCE, k
        differs |= _max_diff(c["lagrangian"], c["co"]) > G4_FACTOR * G1_TOLERANCE
    assert differs


def test_cross_form_spheroid_equals_explicit_with_its_mass_properties():
    sph = _spheroid()
    spheroid_constants, spheroid_function = _build("spheroid", **_spheroid_keywords(sph, sph.rows(2)[0]))
    _, explicit_function = _build(
        "explicit",
        body_mass=spheroid_constants.mass,
        principal_inertia=np.diag(spheroid_constants.inertia),
        body_center_of_gravity=spheroid_constants.center_of_gravity,
        origin=np.zeros(3),
        added_mass_derivatives=spheroid_constants.added_mass_derivatives,
    )
    for k, nu in enumerate(_states()[:200]):
        a, b = _evaluate(spheroid_function, nu, nu), _evaluate(explicit_function, nu, nu)
        for name in ("M", "C_RB", "C_A"):
            assert _max_diff(a[name], b[name]) <= G2_TOLERANCE, (name, k)


def test_cross_form_displacement_hull_equals_explicit_with_its_mass_properties():
    hull = _hull()
    hull_constants, hull_function = _build("displacement_hull", **_hull_keywords(hull, hull.rows(2)[0]))
    _, explicit_function = _build(
        "explicit",
        body_mass=hull_constants.mass,
        principal_inertia=np.diag(hull_constants.inertia),
        body_center_of_gravity=hull_constants.center_of_gravity,
        origin=np.zeros(3),
        added_mass_derivatives=hull_constants.added_mass_derivatives,
    )
    for k, nu in enumerate(_states()[:200]):
        a, b = _evaluate(hull_function, nu, nu), _evaluate(explicit_function, nu, nu)
        for name in ("M", "C_RB", "C_A"):
            assert _max_diff(a[name], b[name]) <= G2_TOLERANCE, (name, k)


def test_cross_form_stabilization_zeroes_exactly_eight_entries():
    sph = _spheroid()
    keywords = _spheroid_keywords(sph, sph.rows(1)[0])
    _, full = _build("spheroid", "co", False, **keywords)
    _, stabilized = _build("spheroid", "co", True, **keywords)
    mask = np.ones((6, 6), dtype=bool)
    for i, j in STABILIZED_PAIRS:
        mask[i, j] = mask[j, i] = False
    assert (~mask).sum() == 8
    for k, nu in enumerate(_states()[:200]):
        a, b = _evaluate(full, nu, nu), _evaluate(stabilized, nu, nu)
        assert np.all(b["C_A"][~mask] == 0.0), k
        assert _max_diff(a["C_A"][mask], b["C_A"][mask]) <= G2_TOLERANCE, k
        assert _max_diff(a["C_RB"], b["C_RB"]) == 0.0 and _max_diff(a["M"], b["M"]) == 0.0


# --------------------------------------------------------------------------
# Physical-sign tests (no model code; MSS is checked, not trusted)
# --------------------------------------------------------------------------
def _all_forms():
    rb, sph, hull = _rbody(), _spheroid(), _hull()
    forms = []
    for set_id in (1, 2, 3):
        for coriolis in ("co", "book", "lagrangian"):
            forms.append((f"explicit-{set_id}-{coriolis}", "explicit", coriolis,
                          _explicit_keywords(rb, rb.rows(set_id)[0], -np.diag(hull.matrix("M_A")[0]))))
    for set_id in (1, 2):
        forms.append((f"spheroid-{set_id}", "spheroid", "co", _spheroid_keywords(sph, sph.rows(set_id)[0])))
        forms.append((f"hull-{set_id}", "displacement_hull", "co", _hull_keywords(hull, hull.rows(set_id)[0])))
    forms.append(("remus100", "spheroid", "co", _remus_keywords()))
    return forms


@pytest.mark.parametrize("label,form,coriolis,keywords", _all_forms(), ids=[f[0] for f in _all_forms()])
def test_physics_coriolis_force_equals_newton_euler(label, form, coriolis, keywords):
    """C_RB(nu) nu equals the Coriolis-centripetal wrench written from
    Newton-Euler about the CG and moved to the CO (no model code); M_RB's
    rotational block is the parallel-axis inertia I + m(|r|^2 1 - r r^T)."""
    constants, function = _build(form, coriolis, **keywords)
    m, inertia, r_g = constants.mass, np.asarray(constants.inertia), np.asarray(constants.center_of_gravity)
    parallel_axis = inertia + m * (float(r_g @ r_g) * np.eye(3) - np.outer(r_g, r_g))
    scale = max(1.0, np.abs(parallel_axis).max())
    assert _max_diff(constants.rigid_body_mass_matrix[3:, 3:], parallel_axis) <= G2_TOLERANCE * scale
    assert np.all(np.diag(constants.rigid_body_mass_matrix)[3:] >= np.diag(inertia) - G2_TOLERANCE * scale)
    for k, nu in enumerate(_states()[:300]):
        exact = _newton_euler_wrench(m, inertia, r_g, nu)
        got = _evaluate(function, nu, nu)["C_RB"] @ nu
        assert _max_diff(got, exact) <= G2_TOLERANCE * max(1.0, np.abs(exact).max()), (label, k)


@pytest.mark.parametrize("label,form,coriolis,keywords", _all_forms(), ids=[f[0] for f in _all_forms()])
def test_physics_no_work_and_positive_kinetic_energy(label, form, coriolis, keywords):
    """Coriolis forces do no work (nu^T C nu = 0, with and without the
    stabilised entries); kinetic energy 1/2 nu^T M nu > 0 (M symmetric, PD)."""
    for stabilize in (False, True):
        constants, function = _build(form, coriolis, stabilize, **keywords)
        total = np.asarray(constants.total_mass_matrix)
        assert _max_diff(total, total.T) <= G2_TOLERANCE * max(1.0, np.abs(total).max())
        assert np.linalg.eigvalsh(0.5 * (total + total.T)).min() > 0.0
        assert np.linalg.eigvalsh(np.asarray(constants.rigid_body_mass_matrix)).min() > 0.0
        for k, nu in enumerate(_states()[:300]):
            out = _evaluate(function, nu, nu)
            for name in ("C_RB", "C_A"):
                work = float(nu @ out[name] @ nu)
                assert abs(work) <= 1e-9 * max(1.0, np.abs(out[name]).max()), (label, name, k)


@pytest.mark.parametrize("coriolis", ["co", "book", "lagrangian"])
def test_physics_centripetal_force_points_to_the_rotation_axis(coriolis):
    """Pure yaw rotation with the CG ahead of the CO: the centripetal force
    m w x (w x r_g) needed to keep the CG on its circle points back towards the
    axis (negative x), magnitude m r^2 x_g."""
    m, x_g, yaw_rate = 50.0, 0.3, 0.8
    _, function = _build("explicit", coriolis, body_mass=m, principal_inertia=[2.0, 3.0, 4.0],
                         body_center_of_gravity=[x_g, 0.0, 0.0], origin=np.zeros(3),
                         added_mass_derivatives=np.zeros(6))
    nu = np.array([0.0, 0.0, 0.0, 0.0, 0.0, yaw_rate])
    force = (_evaluate(function, nu, nu)["C_RB"] @ nu)[:3]
    assert force[0] < 0.0
    np.testing.assert_allclose(force, [-m * yaw_rate**2 * x_g, 0.0, 0.0], atol=1e-12, rtol=1e-12)


def test_physics_munk_moment_is_destabilizing_for_a_prolate_spheroid():
    """Lamb: a slender body has more added mass sideways than along (k2 > k1).
    Moving forward and sideways (u, v > 0) the added-mass yaw moment on the body,
    -(C_A nu)_6 = -(A22 - A11) u v, turns the nose away from the velocity
    (negative yaw), i.e. it increases the drift angle; likewise in pitch,
    -(C_A nu)_5 = (A33 - A11) u w. The stabilised form removes both."""
    constants, function = _build("spheroid", **_remus_keywords())
    _, stabilized = _build("spheroid", "co", True, **_remus_keywords())
    a11, a22, a33 = np.diag(constants.added_mass_matrix)[:3]
    assert a22 > a11 > 0.0 and a33 > a11
    u, v, w = 1.5, 0.4, 0.3
    yaw = -(_evaluate(function, np.zeros(6), [u, v, 0, 0, 0, 0])["C_A"] @ [u, v, 0, 0, 0, 0])[5]
    pitch = -(_evaluate(function, np.zeros(6), [u, 0, w, 0, 0, 0])["C_A"] @ [u, 0, w, 0, 0, 0])[4]
    assert yaw < 0.0 and pitch > 0.0
    np.testing.assert_allclose(yaw, -(a22 - a11) * u * v, rtol=1e-12)
    np.testing.assert_allclose(pitch, (a33 - a11) * u * w, rtol=1e-12)
    yaw_s = (_evaluate(stabilized, np.zeros(6), [u, v, 0, 0, 0, 0])["C_A"] @ [u, v, 0, 0, 0, 0])[5]
    pitch_s = (_evaluate(stabilized, np.zeros(6), [u, 0, w, 0, 0, 0])["C_A"] @ [u, 0, w, 0, 0, 0])[4]
    assert yaw_s == 0.0 and pitch_s == 0.0


def test_physics_added_mass_is_positive_for_negative_derivatives():
    """Hydrodynamic derivatives X_du ... N_dr < 0 mean the fluid adds inertia:
    M_A = -diag(derivatives) has a positive diagonal (SNAME sign convention)."""
    hull = _hull()
    for set_id in (1, 2):
        constants, _ = _build("displacement_hull", **_hull_keywords(hull, hull.rows(set_id)[0]))
        assert np.all(constants.added_mass_derivatives < 0.0)
        assert np.all(np.diag(constants.added_mass_matrix) > 0.0)


# --------------------------------------------------------------------------
# Generic: two parameter sets per form (above), no vehicle in the code
# --------------------------------------------------------------------------
def _code_identifiers_and_strings(path):
    tree = ast.parse(path.read_text())
    docstrings = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            if node.body and isinstance(node.body[0], ast.Expr) and isinstance(getattr(node.body[0], "value", None), ast.Constant):
                docstrings.add(id(node.body[0].value))
    found = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Name):
            found.append(node.id)
        elif isinstance(node, ast.Attribute):
            found.append(node.attr)
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            found.append(node.name)
        elif isinstance(node, ast.arg):
            found.append(node.arg)
        elif isinstance(node, ast.keyword) and node.arg:
            found.append(node.arg)
        elif isinstance(node, ast.Constant) and isinstance(node.value, str) and id(node) not in docstrings:
            found.append(node.value)
    return found


def test_generic_block_names_no_vehicle():
    for path in _package_files():
        for token in _code_identifiers_and_strings(path):
            lowered = token.lower()
            for name in VEHICLE_NAMES:
                assert name not in lowered, (path.name, token)


# --------------------------------------------------------------------------
# Transforms only from more_transformations
# --------------------------------------------------------------------------
LOCAL_TRANSFORM_NAMES = {"skew", "_skew", "_skew_casadi", "H_matrix", "transform_matrix", "Smtrx", "Hmtrx"}


def test_transforms_no_local_definitions_and_imports_from_more_transformations():
    numpy_import = casadi_import = False
    for path in _package_files():
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                assert node.name not in LOCAL_TRANSFORM_NAMES, (path.name, node.name)
            if isinstance(node, ast.ImportFrom) and node.module:
                if node.module.startswith("more_transformations.more_casadi_transformations"):
                    casadi_import = True
                elif node.module.startswith("more_transformations"):
                    numpy_import = True
            if isinstance(node, ast.Import):
                for alias in node.names:
                    if alias.name.startswith("more_transformations.more_casadi_transformations"):
                        casadi_import = True
                    elif alias.name.startswith("more_transformations"):
                        numpy_import = True
    assert not numpy_import, "import from the numpy more_transformations (the block is all CasADi)"
    assert casadi_import, "no import from more_transformations.more_casadi_transformations (graph)"


# --------------------------------------------------------------------------
# The all-CasADi block: committed block, frozen path, gradient, no numpy
# --------------------------------------------------------------------------
COMMITTED_REVISION = "5d67caa"  # last revision of the block with numpy pre-processing
COMMITTED_FILES = ("__init__.py", "added_mass.py", "constants.py", "kinetics.py", "planar.py")
COMMITTED_SEED = 20261007
N_COMMITTED_STATES = 300
GRADIENT_TOLERANCE = 1e-8
GRADIENT_STEP = 1e-4  # kg, central difference

# config/dataclass/plant/asv_catamaran_params.py 8-31 (as PARAMETERS of
# test_rigid_body_block.py, line per key there).
HULL_WITH_PAYLOAD_VALUES = {
    "length": 2.0,
    "beam": 1.08,
    "water_density": 1025.0,
    "hull_mass": 55.0,
    "payload_mass": 25.0,
    "hull_center_of_gravity": [0.2, 0.0, -0.2],
    "payload_position": [0.05, 0.0, -0.35],
    "added_mass_coefficients": [-1.0, -1.5, -1.0, -0.2, -0.8, -1.7],
    "radii_of_gyration": [0.4, 0.25, 0.25],
}


def _form_values():
    """One parameter set per form (the frozen CSV sets and the catamaran)."""
    rb, sph, hull = _rbody(), _spheroid(), _hull()
    return {
        "hull_with_payload": HULL_WITH_PAYLOAD_VALUES,
        "displacement_hull": _hull_keywords(hull, hull.rows(2)[0]),
        "spheroid": _spheroid_keywords(sph, sph.rows(1)[0]),
        "explicit": _explicit_keywords(rb, rb.rows(1)[0], -np.diag(hull.matrix("M_A")[0])),
    }


def _committed_block(tmp_path):
    """The block as committed at ``COMMITTED_REVISION``, read with ``git
    show`` into a temporary package (method: the old modules are imported,
    no values are frozen in a file). Skips when git or the revision is not
    available (a source archive without history)."""
    root = Path(__file__).resolve().parents[2]
    package = tmp_path / "rigid_body_committed"
    package.mkdir()
    for name in COMMITTED_FILES:
        try:
            shown = subprocess.run(
                ["git", "-C", str(root), "show",
                 f"{COMMITTED_REVISION}:more_dynamics/models/rigid_body/{name}"],
                capture_output=True, text=True, check=True,
            )
        except (OSError, subprocess.CalledProcessError) as exc:
            pytest.skip(f"revision {COMMITTED_REVISION} not readable with git here: {exc}")
        (package / name).write_text(shown.stdout)
    sys.path.insert(0, str(tmp_path))
    try:
        return importlib.import_module("rigid_body_committed")
    finally:
        sys.path.remove(str(tmp_path))


def _committed_keywords(form, values):
    """The keyword names of the committed block for the same numbers."""
    old = dict(values)
    if form == "displacement_hull":
        old["center_of_gravity"] = old.pop("hull_center_of_gravity")
    elif form == "spheroid":
        old["center_of_gravity"] = old.pop("body_center_of_gravity")
    elif form == "explicit":
        old["mass"] = old.pop("body_mass")
        old["inertia"] = np.diag(old.pop("principal_inertia"))
        old["center_of_gravity"] = old.pop("body_center_of_gravity")
    return old


def test_G2_block_equals_the_committed_block_on_all_forms(tmp_path):
    """Every form x Coriolis form x stabilisation: M, C_RB, C_A and the
    state-independent quantities of the block equal the block committed at
    ``COMMITTED_REVISION`` (numpy pre-processing) on 300 seeded states."""
    committed = _committed_block(tmp_path)
    rng = np.random.default_rng(COMMITTED_SEED)
    nus = rng.uniform(-3.0, 3.0, size=(N_COMMITTED_STATES, 6))
    nu_rs = rng.uniform(-3.0, 3.0, size=(N_COMMITTED_STATES, 6))
    for form, values in _form_values().items():
        for coriolis in ("co", "book", "lagrangian"):
            for stabilize in (False, True):
                label = (form, coriolis, stabilize)
                constants, function = _build(form, coriolis, stabilize, **values)
                old_constants = committed.preprocess_rigid_body(
                    mass_properties=form, coriolis=coriolis,
                    stabilize_added_mass_coriolis=stabilize, **_committed_keywords(form, values))
                old_function = committed.rigid_body_casadi(old_constants)
                for name in ("mass", "center_of_gravity", "inertia", "rigid_body_mass_matrix",
                             "added_mass_matrix", "total_mass_matrix", "added_mass_derivatives",
                             "lamb_k_factors", "displaced_volume", "wetted_surface"):
                    old_value = getattr(old_constants, name)
                    if old_value is None:
                        assert getattr(constants, name) is None, (label, name)
                    else:
                        assert _max_diff(getattr(constants, name), old_value) <= G2_TOLERANCE, (label, name)
                for k in range(N_COMMITTED_STATES):
                    new = _evaluate(function, nus[k], nu_rs[k])
                    old = old_function(nu=nus[k], nu_r=nu_rs[k])
                    for name in ("M", "C_RB", "C_A"):
                        assert _max_diff(new[name], np.array(old[name], dtype=float)) <= G2_TOLERANCE, (label, name, k)


@pytest.mark.parametrize("form", ["hull_with_payload", "displacement_hull", "spheroid", "explicit"])
def test_frozen_block_equals_the_block_called_with_numbers(form):
    """``freeze`` (what a plugin carries: a function of the state only)
    against the block called with the numbers, 300 seeded states; entry by
    entry within G2 relative to max(1, |value|)."""
    from more_transformations.more_casadi_transformations import freeze

    block = _contract()
    values = _form_values()[form]
    _, (function, checked) = _build(form, "co", True, **values)
    frozen = freeze(function, block.rigid_body_parameters(form), values)
    assert frozen.name_in() == ["nu", "nu_r"]
    assert frozen.name_out() == list(block.rigid_body_outputs(form))
    rng = np.random.default_rng(COMMITTED_SEED)
    for k in range(N_COMMITTED_STATES):
        nu, nu_r = rng.uniform(-3.0, 3.0, size=(2, 6))
        a, b = frozen(nu=nu, nu_r=nu_r), function(nu=nu, nu_r=nu_r, **checked)
        for name in frozen.name_out():
            unfrozen = np.array(b[name], dtype=float)
            error = np.abs(np.array(a[name], dtype=float) - unfrozen) / np.maximum(1.0, np.abs(unfrozen))
            assert error.max() <= G2_TOLERANCE, (form, name, k)


def test_gradient_of_the_wrench_with_respect_to_hull_mass():
    """Identification path: the block called with ``hull_mass`` left as a
    symbol; d(M nu_dot + C_RB(nu) nu + C_A(nu_r) nu_r)/d hull_mass from CasADi
    equals a central difference of the block called with numbers."""
    import casadi as ca

    function = _contract().rigid_body_casadi()
    rng = np.random.default_rng(COMMITTED_SEED)
    m0 = HULL_WITH_PAYLOAD_VALUES["hull_mass"]

    def wrench(hull_mass, nu, nu_r, nu_dot):
        out = function(nu=nu, nu_r=nu_r, **{**HULL_WITH_PAYLOAD_VALUES, "hull_mass": hull_mass})
        return out["M"] @ nu_dot + out["C_RB"] @ nu + out["C_A"] @ nu_r

    symbol = ca.SX.sym("hull_mass")
    for k in range(20):
        nu, nu_r = rng.uniform(-3.0, 3.0, size=(2, 6))
        nu_dot = rng.uniform(-1.0, 1.0, size=6)
        gradient = ca.Function("gradient", [symbol], [ca.jacobian(wrench(symbol, nu, nu_r, nu_dot), symbol)])
        exact = np.array(gradient(m0), dtype=float).ravel()
        upper = np.array(wrench(m0 + GRADIENT_STEP, nu, nu_r, nu_dot), dtype=float).ravel()
        lower = np.array(wrench(m0 - GRADIENT_STEP, nu, nu_r, nu_dot), dtype=float).ravel()
        difference = (upper - lower) / (2.0 * GRADIENT_STEP)
        assert _max_diff(exact, difference) <= GRADIENT_TOLERANCE, k
        assert np.abs(exact).max() > 1e-2, k  # the gradient is not trivially zero


def test_block_imports_no_numpy():
    """No module of the block imports numpy or scipy, or the numpy
    ``more_transformations`` modules (AST scan of every file)."""
    for path in _package_files():
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            modules = []
            if isinstance(node, ast.Import):
                modules = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module:
                modules = [node.module]
            for module in modules:
                assert module.split(".")[0] not in ("numpy", "scipy"), (path.name, module)
                if module.startswith("more_transformations"):
                    assert module.startswith("more_transformations.more_casadi_transformations"), (path.name, module)
