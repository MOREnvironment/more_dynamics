"""Gate tests for the rigid-body + added-mass block, **all forms**.

Written 2026-10-06, before the forms were ported. The rules they encode
(the owner's decisions of 2026-10-05/06): the block is generic, with no
vehicle name or number in its code; the newest MSS is the reference, not the
truth, and every deviation from it is stated; there is no
``legacy_otter_inertia`` flag (its test is in ``test_rigid_body_block.py``);
skew and H come only from ``more_transformations``. The catamaran path
(``hull_with_payload``) keeps its gates in ``test_rigid_body_block.py``; this
file adds the other forms, the cross-form and physical-sign tests, the
transforms-import test and the generic test.

Gates (the test names carry them): G1 the block against MATLAB running MSS
(frozen CSVs), 1e-9 absolute; G2 the block against the numpy source or an
independent transcription, 1e-10 absolute; G4 a perturbed model must differ
from the reference by more than 10x the G1 tolerance; G5 a printed number
reproduced at its printed precision.

Contract of the block (nothing renamed since the first version)
---------------------------------------------------------------
``more_dynamics.models.rigid_body`` exports:

* ``preprocess_rigid_body(*, mass_properties="hull_with_payload", coriolis="co",
  stabilize_added_mass_coriolis=False, **form_keywords) -> RigidBodyConstants``
  (numpy). A form is chosen by these explicit arguments, never by a vehicle
  name. Unknown ``mass_properties`` or ``coriolis`` -> ``ValueError``; a missing
  or foreign keyword for the chosen form -> ``TypeError`` or ``ValueError``.

  ``mass_properties`` (how m, I about the CG, r_g and M_A are obtained):

  - ``"hull_with_payload"`` (default; the checked catamaran path, unchanged):
    ``length, beam, water_density, hull_mass, payload_mass,
    hull_center_of_gravity, payload_position, added_mass_coefficients,
    radii_of_gyration``.
  - ``"displacement_hull"`` (source ``get_params_hull``): ``water_density,
    length, beam, draft, block_coefficient, radii_of_gyration,
    center_of_gravity, added_mass_coefficients``. ``m = rho Cb L B T``,
    ``I = m diag((R_s * [B, L, L])^2)`` about the CG, M_A by scaled derivatives
    ``c * [A11(m, L, rho), m, m, I11, I22, I33]`` (the ``otter.m`` 152-157
    pattern: the rotational terms scale the inertia **about the CG**).
  - ``"spheroid"`` (source ``get_I_gb`` + ``get_added_mass_derivates`` +
    ``M_A_lamb_6dof``; MSS ``spheroid.m`` + ``imlay61.m``):
    ``semi_major_axis, semi_minor_axis, body_density, water_density,
    roll_added_inertia_ratio, center_of_gravity``. Body: ``m = body_density
    4/3 pi a b^2``, ``I = diag(2/5 m b^2, 1/5 m (a^2+b^2), 1/5 m (a^2+b^2))``.
    Added mass (Lamb k-factors) from the displaced fluid, ``m_f`` and ``I_f``
    with ``water_density``: ``diag(k1 m_f, k2 m_f, k2 m_f, r44 I_f,x,
    k' I_f,y, k' I_f,y)``. Two densities because MSS uses 1025 in
    ``spheroid.m`` (line 35) and 1026 in ``imlay61.m`` (line 31); the numpy
    source uses one ``rho`` for both.
  - ``"explicit"`` (source ``RigidBody6DOF`` + ``M_A_6dof``): ``mass, inertia``
    (3x3 about the CG, **diagonal**, as the source line 43-44; otherwise
    ``ValueError``), ``center_of_gravity``, ``added_mass_derivatives``
    (``[X_du, Y_dv, Z_dw, K_dp, M_dq, N_dr]``, ``M_A = -diag``), optional
    ``origin`` (default zeros; ``r_g = center_of_gravity - origin``, source
    ``get_r_bg``).

  ``coriolis`` (the C_RB parametrisation; all give the same ``C_RB(nu) nu``):

  - ``"co"`` (default; MSS ``rbody.m`` / ``spheroid.m``, source ``get_C_RB_co``,
    ``C_RB``, ``C_RB_``, ``C_RB_explicit``, ``C_RB_runtime``):
    ``H^T diag(m S(w), -S(I w)) H``. Source ``get_C_RB_book_corrected`` equals
    it algebraically (dropped as a separate form by the owner's decision of
    2026-10-06; the source function maps here).
  - ``"book"`` (source ``get_C_RB_book``): ``[[m S(w), -m S(w) S(r_g)],
    [m S(r_g) S(w), -S(I_O w)]]``, ``I_O = I - m S(r_g)^2``. Differs from
    ``"co"`` only in the moment block, by ``-m |r_g|^2 S(w)``.
  - ``"lagrangian"`` (MSS ``m2c.m``, source ``get_C_RB_lagrangian``):
    ``m2c(M_RB, nu)``, depends on the linear velocity too.

  ``C_A(nu_r)`` is always ``m2c(M_A, nu_r)`` (MSS ``m2c.m`` with its
  symmetrisation; source ``C_A_6dof``, ``C_A_6dof_lagrangian``,
  ``AddedMass.get_C_RB_lagrangian``). ``stabilize_added_mass_coriolis=True``
  zeroes the eight entries of MSS ``remus100.m`` lines 207-210 (source
  ``stabilize_C_A``).
* ``RigidBodyConstants`` (frozen dataclass), at least: ``mass_properties``,
  ``coriolis``, ``stabilize_added_mass_coriolis``, ``mass``,
  ``center_of_gravity`` (3,), ``inertia`` (3x3, about the CG),
  ``rigid_body_mass_matrix``, ``added_mass_matrix``, ``total_mass_matrix``
  (6x6), ``added_mass_derivatives`` (6,) ``= -diag(M_A)``; ``lamb_k_factors``
  (3,) for ``"spheroid"``; ``displaced_volume`` and ``wetted_surface``
  (``L B + 2 T B``, used by the surge damping block) for
  ``"displacement_hull"``.
* ``rigid_body_casadi(constants) -> ca.Function`` unchanged: inputs
  ``["nu", "nu_r"]`` (6x1), outputs ``["M", "C_RB", "C_A"]`` (6x6).
* 3-DOF forms (kept as a **reduction of the 6-DOF block**, owner's decision
  of 2026-10-06; ``planar.py``): ``PLANAR_DOFS == (0, 1, 5)``;
  ``planar_reduction(matrix) -> np.ndarray`` (rows and columns ``PLANAR_DOFS``
  of a 6x6 matrix); ``planar_added_mass_matrix(X_du, Y_dv, Y_dr, N_dr) ->
  np.ndarray`` (3x3, ``-[[X_du,0,0],[0,Y_dv,Y_dr],[0,Y_dr,N_dr]]``);
  ``planar_coriolis_casadi(mass_matrix) -> ca.Function`` named
  ``"planar_coriolis"``, input ``["nu"]`` (3x1 ``[u, v, r]``), output ``["C"]``
  (3x3), MSS ``m2c.m`` lines 50-54; ``planar_casadi(constants) ->
  ca.Function``, inputs ``["nu", "nu_r"]`` (3x1), outputs ``["M", "C_RB",
  "C_A"]`` (3x3), the surge, sway, yaw reduction of
  ``rigid_body_casadi(constants)`` at ``[u, v, 0, 0, 0, r]``.
* Transforms from ``more_transformations`` only: ``skew``/``H_matrix`` from
  ``more_transformations.matrix_transforms`` (numpy) and
  ``more_transformations.more_casadi_transformations`` (graph).

Numpy source read in full (``more_generic_models/more_generic_models/``,
``524e336``), ``dynamics/plant/matrices/``:

* ``rigid_body_kinetics.py``: ``RigidBody6DOF.__init__`` 23-57 (diagonal ``I_g``
  required 40-44; ``M_RB = M_RB_book``), ``get_r_bg`` 72-77 (warns when the CG is
  not below the CO), ``normalize_zero`` 80-81, ``_construct_Igb`` 87-123,
  ``_compute_Ibb`` 127-138 (``I_O = I_g - m S^2``), ``get_M_RB_book`` 143-154,
  ``get_M_RB_co`` 156-176, ``get_M_RB`` 179-206, ``get_C_RB_book`` 210-233,
  ``get_C_RB_book_corrected`` 236-267, ``get_C_RB_co`` 270-300,
  ``get_C_RB_lagrangian`` 303-352 (6-DOF and 3-DOF ``m2c``), ``get_C_RB``
  357-402.
* ``rigid_body_kinetics_spheroid.py``: ``get_I_gb`` 31-59, ``M_RB`` 62-88 (uses
  the diagonal of ``I`` only), ``C_RB`` 95-127, ``C_RB_`` 136-159,
  ``C_RB_explicit`` 162-197, ``get_RB_matrices`` 203-207.
* ``rigid_body_kinetics_surface_vessel.py``: ``get_params_hull`` 54-94 (returns
  ``Ig_co = Ig_CG - m S(r_bg)^2`` and stores it as ``self.I_gb``, 83 and 91),
  ``M_RB`` 97-121, ``C_RB`` 124-138, ``C_RB_runtime`` 147-170 (uses
  ``self.I_gb``).
* ``added_mass.py``: ``normalize_zero`` 17-18, ``M_A_6dof`` 24-41,
  ``M_A_lamb_6dof`` 45-78, ``M_A_3dof`` 81-88, ``C_A_6dof`` 97-113,
  ``C_A_6dof_lagrangian`` 116-122, ``C_A_3dof`` 126-152, ``stabilize_C_A``
  156-171, ``get_added_mass_derivates`` 174-267, ``get_added_mass_derivates_hull``
  291-310 (same formula as ``_catamaran`` 269-289), ``get_C_RB_lagrangian``
  363-388 (no symmetrisation).
* Callers: ``plant/auv_spheroid/auv_spheroid.py`` 158-228 (spheroid + Lamb +
  stabilised C_A), ``plant/asv_hull/asv_hull.py`` 84-177 (``M_RB`` with
  ``Ig_cg`` but ``C_RB_runtime`` and the rotational added mass with
  ``Ig_co``), ``plant/asv_catamaran/asv_catamaran.py`` (hull_with_payload).

MSS (the checkout named by ``MSS_DIR``; read at ``72656d1``, release 2.0.2;
the code of every file below is unchanged since ``99bf0b3``): ``LIBRARY/modeling/rbody.m`` 31-45, ``m2c.m`` 33-56,
``spheroid.m`` 32-52, ``imlay61.m`` 30-62, ``addedMassSurge.m`` 32-34,
``LIBRARY/kinematics/Hmtrx.m`` 16-18, ``CRAFT/AUV/models/remus100.m`` 131-137,
199-210, ``CRAFT/USV/models/otter.m`` 152-159.

Frozen references (``tests/data/rigid_body/SOURCE.md``, section of 2026-10-06):
``rigid_body_rbody_mss_current.csv``, ``rigid_body_spheroid_mss_current.csv``,
``rigid_body_hull_mss_current.csv``, ``rigid_body_m2c_3dof_mss_current.csv``
(MATLAB R2026a running MSS), ``spheroid_matlab_reference_mss_current.csv``
(columns of the 2026-10-05 ``remus100.m`` CSV).

The physical tests below are written from Newton-Euler mechanics in the test
itself (CG velocity, force and moment about the CG moved to the CO: Fossen,
T. I. (2011), *Handbook of Marine Craft Hydrodynamics and Motion Control*,
Wiley, eqs. 3.14-3.18, pp. 47-48, and 3.33-3.40, pp. 50-51), not from model
code.

Outside this repo (nothing relative to one machine): ``MSS_DIR`` (MSS checkout) for the
pinned-line test, ``MORE_GENERIC_MODELS_DIR`` (repository root of
``more_generic_models``) for the G2 tests; unset -> those tests skip and name
the variable. Every other test needs only ``tests/data/rigid_body/``.
"""

import ast
import contextlib
import importlib
import inspect
import io
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "rigid_body"
# Outside this repo (nothing relative to one machine): found only through these variables,
# default not set; a test that needs one skips with a message naming it. The
# frozen CSVs in tests/data/ are the only hard dependency.
MSS_DIR_VARIABLE = "MSS_DIR"                                  # MSS checkout root
MORE_GENERIC_MODELS_DIR_VARIABLE = "MORE_GENERIC_MODELS_DIR"  # repo root of more_generic_models
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
    "LIBRARY/modeling/m2c.m": {33: "M = 0.5 * (M + M');"},
}
REMUS = {  # remus100.m lines 131-137
    "L_auv": 1.6, "D_auv": 0.19, "scale": 1.0096, "r44": 0.3, "r_bG": [0.0, 0.0, 0.02],
}
SPHEROID_M_DENSITY = 1025.0  # spheroid.m line 35
IMLAY61_DENSITY = 1026.0     # imlay61.m line 31
REMUS_PRINTED_MASS = 31.9    # remus100.m line 134 comment ("to obtain m = 31.9 kg")

# remus100.m 207-210 (source stabilize_C_A 160-169), 0-based (row, col) pairs.
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
    block = _contract()
    constants = block.preprocess_rigid_body(
        mass_properties=mass_properties,
        coriolis=coriolis,
        stabilize_added_mass_coriolis=stabilize,
        **keywords,
    )
    return constants, block.rigid_body_casadi(constants)


def _evaluate(function, nu, nu_r):
    out = function(nu=np.asarray(nu, float), nu_r=np.asarray(nu_r, float))
    return {name: np.array(out[name], dtype=float) for name in ("M", "C_RB", "C_A")}


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


def _explicit_keywords(csv, row, added_mass_derivatives=(0.0,) * 6):
    """rbody.m line 35: I_G = m diag(R44^2, R55^2, R66^2)."""
    m = csv.column("m")[row]
    radii = np.array([csv.column(f"R{i}")[row] for i in (44, 55, 66)])
    return {
        "mass": m,
        "inertia": m * np.diag(radii**2),
        "center_of_gravity": csv.vector("r_bG", 3)[row],
        "added_mass_derivatives": np.asarray(added_mass_derivatives, float),
    }


def _spheroid_keywords(csv, row, body_density=SPHEROID_M_DENSITY, water_density=IMLAY61_DENSITY):
    return {
        "semi_major_axis": csv.column("a")[row],
        "semi_minor_axis": csv.column("b")[row],
        "body_density": body_density,
        "water_density": water_density,
        "roll_added_inertia_ratio": csv.column("r44")[row],
        "center_of_gravity": csv.vector("r_bG", 3)[row],
    }


def _hull_keywords(csv, row):
    return {
        "water_density": csv.column("rho")[row],
        "length": csv.column("L")[row],
        "beam": csv.column("B")[row],
        "draft": csv.column("T")[row],
        "block_coefficient": csv.column("Cb")[row],
        "radii_of_gyration": csv.vector("Rs", 3)[row],
        "center_of_gravity": csv.vector("r_bG", 3)[row],
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
        "center_of_gravity": REMUS["r_bG"],
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


def _source(module, name):
    """The numpy source, from the repository named by MORE_GENERIC_MODELS_DIR."""
    root = os.environ.get(MORE_GENERIC_MODELS_DIR_VARIABLE)
    if not root:
        pytest.skip(
            f"{MORE_GENERIC_MODELS_DIR_VARIABLE} is not set (repository root of "
            "more_generic_models); G2 against the numpy source skipped"
        )
    if root not in sys.path:
        sys.path.insert(0, root)
    try:
        mod = importlib.import_module(f"more_generic_models.dynamics.plant.{module}")
    except ImportError as exc:
        pytest.skip(f"numpy source not importable from {MORE_GENERIC_MODELS_DIR_VARIABLE}={root}: {exc}")
    return getattr(mod, name)


def _quiet(call, *args, **kwargs):
    """Source dispatchers print and RigidBody6DOF warns (get_r_bg)."""
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        return call(*args, **kwargs)


def _package_files():
    module = _contract()
    if hasattr(module, "__path__"):
        return sorted(Path(p) for base in module.__path__ for p in Path(base).glob("*.py"))
    return [Path(module.__file__)]


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
    parameters = inspect.signature(block.preprocess_rigid_body).parameters
    for name in ("mass_properties", "coriolis", "stabilize_added_mass_coriolis"):
        assert name in parameters, name
    assert parameters["mass_properties"].default == "hull_with_payload"
    assert parameters["coriolis"].default == "co"
    assert parameters["stabilize_added_mass_coriolis"].default is False
    assert "legacy_otter_inertia" not in parameters  # flag dropped (owner, 2026-10-05)


def test_contract_rejects_unknown_forms_and_keywords():
    rb = _rbody()
    keywords = _explicit_keywords(rb, rb.rows(1)[0])
    block = _contract()
    with pytest.raises(ValueError):
        block.preprocess_rigid_body(mass_properties="otter", **keywords)
    with pytest.raises(ValueError):
        block.preprocess_rigid_body(mass_properties="explicit", coriolis="book_corrected_typo", **keywords)
    missing = dict(keywords)
    missing.pop("mass")
    with pytest.raises((TypeError, ValueError)):
        block.preprocess_rigid_body(mass_properties="explicit", **missing)
    with pytest.raises((TypeError, ValueError)):
        block.preprocess_rigid_body(mass_properties="explicit", draft=0.3, **keywords)
    off_diagonal = dict(keywords)
    off_diagonal["inertia"] = keywords["inertia"] + np.array([[0, 1.0, 0], [1.0, 0, 0], [0, 0, 0]])
    with pytest.raises(ValueError):  # source rigid_body_kinetics.py 43-44
        block.preprocess_rigid_body(mass_properties="explicit", **off_diagonal)
    sph = _spheroid()
    swapped = _spheroid_keywords(sph, 0)
    swapped["semi_major_axis"], swapped["semi_minor_axis"] = swapped["semi_minor_axis"], swapped["semi_major_axis"]
    with pytest.raises(ValueError):  # imlay61.m line 42; source added_mass.py 181-182
        block.preprocess_rigid_body(mass_properties="spheroid", **swapped)


@pytest.mark.parametrize("coriolis", ["co", "book", "lagrangian"])
def test_contract_function_and_constants(coriolis):
    rb = _rbody()
    constants, function = _build("explicit", coriolis, **_explicit_keywords(rb, rb.rows(1)[0]))
    assert function.name_in() == ["nu", "nu_r"]
    assert function.name_out() == ["M", "C_RB", "C_A"]
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
    m |r_g|^2 S(w) in the moment block (identity of the source forms; Fossen 2011, eq. 3.30, p. 50)."""
    rb = _rbody()
    rows = rb.rows(set_id)
    keywords = _explicit_keywords(rb, rows[0])
    _, function = _build("explicit", "book", **keywords)
    m, r_g = keywords["mass"], keywords["center_of_gravity"]
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
    mass_matrix = block.planar_added_mass_matrix(*args)
    assert _max_diff(mass_matrix, dof3.matrix("M3", 3)[rows[0]]) <= G1_TOLERANCE
    function = block.planar_coriolis_casadi(mass_matrix)
    assert function.name_in() == ["nu"] and function.name_out() == ["C"]
    assert function.size_in(0) == (3, 1) and function.size_out(0) == (3, 3)
    for k in rows:
        c = np.array(function(nu=dof3.vector("nu3", 3)[k])["C"], dtype=float)
        assert _max_diff(c, dof3.matrix("C3", 3)[k]) <= G1_TOLERANCE, k


# --------------------------------------------------------------------------
# G2 — block vs numpy source (50 frozen + 1000 seeded states)
# --------------------------------------------------------------------------
@pytest.mark.parametrize("set_id", [1, 2, 3])
def test_G2_explicit_matches_RigidBody6DOF_all_modes(set_id):
    RigidBody6DOF = _source("matrices.rigid_body_kinetics", "RigidBody6DOF")
    AddedMass = _source("matrices.added_mass", "AddedMass")
    rb, hull = _rbody(), _hull()
    derivatives = -np.diag(hull.matrix("M_A")[hull.rows(1 + (set_id - 1) % 2)[0]])
    keywords = _explicit_keywords(rb, rb.rows(set_id)[0], derivatives)
    source = _quiet(RigidBody6DOF, keywords["mass"], keywords["inertia"], r_cg=keywords["center_of_gravity"])
    modes = {"co": "co", "book": "book", "lagrangian": "lagrangian"}
    functions = {}
    for coriolis in modes:
        constants, functions[coriolis] = _build("explicit", coriolis, **keywords)
        assert _max_diff(constants.rigid_body_mass_matrix, source.M_RB_book) <= G2_TOLERANCE
        assert _max_diff(constants.rigid_body_mass_matrix, source.M_RB_co) <= G2_TOLERANCE
        assert _max_diff(constants.rigid_body_mass_matrix, _quiet(source.get_M_RB, "co")) <= G2_TOLERANCE
        assert _max_diff(constants.added_mass_matrix, AddedMass.M_A_6dof(derivatives)) <= G2_TOLERANCE
    _, stabilized = _build("explicit", "co", True, **keywords)
    for k, nu in enumerate(_states()):
        for coriolis, mode in modes.items():
            out = _evaluate(functions[coriolis], nu, nu)
            assert _max_diff(out["C_RB"], _quiet(source.get_C_RB, nu, mode)) <= G2_TOLERANCE, (coriolis, k)
            assert _max_diff(out["C_A"], AddedMass.C_A_6dof(nu, *derivatives)) <= G2_TOLERANCE, (coriolis, k)
            assert _max_diff(out["C_A"], AddedMass.C_A_6dof_lagrangian(nu, AddedMass.M_A_6dof(derivatives))) <= G2_TOLERANCE
        # book_corrected (dropped as a form, owner 2026-10-06) is the "co" matrix
        assert _max_diff(_evaluate(functions["co"], nu, nu)["C_RB"],
                         _quiet(source.get_C_RB, nu, "book_corrected")) <= G2_TOLERANCE, k
        expected = AddedMass.stabilize_C_A(AddedMass.C_A_6dof_lagrangian(nu, AddedMass.M_A_6dof(derivatives)))
        assert _max_diff(_evaluate(stabilized, nu, nu)["C_A"], expected) <= G2_TOLERANCE, k


def test_G2_explicit_origin_is_get_r_bg():
    """source get_r_bg (72-77): r_bg = r_cg - r_co."""
    RigidBody6DOF = _source("matrices.rigid_body_kinetics", "RigidBody6DOF")
    rb = _rbody()
    keywords = _explicit_keywords(rb, rb.rows(1)[0])
    origin = np.array([0.05, -0.02, 0.01])
    shifted = {**keywords, "center_of_gravity": keywords["center_of_gravity"] + origin, "origin": origin}
    constants, _ = _build("explicit", **shifted)
    np.testing.assert_allclose(constants.center_of_gravity, RigidBody6DOF.get_r_bg(
        shifted["center_of_gravity"], origin), atol=G2_TOLERANCE, rtol=0.0)
    reference, _ = _build("explicit", **keywords)
    assert _max_diff(constants.rigid_body_mass_matrix, reference.rigid_body_mass_matrix) <= G2_TOLERANCE


@pytest.mark.parametrize("set_id", [1, 2])
def test_G2_spheroid_matches_source(set_id):
    """Source uses one rho for body and added mass: both densities = rho here."""
    Spheroid = _source("matrices.rigid_body_kinetics_spheroid", "RigidBody6DOFSpheroid")
    AddedMass = _source("matrices.added_mass", "AddedMass")
    sph = _spheroid()
    row = sph.rows(set_id)[0]
    rho = (1025.0, 1026.0)[set_id - 1]
    keywords = _spheroid_keywords(sph, row, body_density=rho, water_density=rho)
    a, b, r44, r_g = (keywords[k] for k in ("semi_major_axis", "semi_minor_axis",
                                            "roll_added_inertia_ratio", "center_of_gravity"))
    source = Spheroid(np.asarray(r_g, float))
    inertia, mass, _ = source.get_I_gb(rho, a, b)
    derivatives, k_factors = AddedMass.get_added_mass_derivates(mass, a, b)
    m_a = AddedMass.M_A_lamb_6dof(mass, inertia, k_factors, r44)
    constants, function = _build("spheroid", "co", False, **keywords)
    _, stabilized = _build("spheroid", "co", True, **keywords)
    assert abs(constants.mass - mass) <= G2_TOLERANCE * mass
    assert _max_diff(constants.inertia, inertia) <= G2_TOLERANCE
    assert _max_diff(constants.rigid_body_mass_matrix, Spheroid.M_RB(mass, inertia, r_g)) <= G2_TOLERANCE
    assert _max_diff(constants.added_mass_matrix, m_a) <= G2_TOLERANCE
    assert _max_diff(constants.lamb_k_factors, k_factors) <= G2_TOLERANCE
    # Fossen derivative set (get_added_mass_derivates) = Lamb with K_dp = 0
    expected_derivatives = derivatives.copy()
    expected_derivatives[3] = -r44 * inertia[0, 0]
    assert _max_diff(constants.added_mass_derivatives, expected_derivatives) <= G2_TOLERANCE
    for k, nu in enumerate(_states()):
        out = _evaluate(function, nu, nu)
        assert _max_diff(out["C_RB"], Spheroid.C_RB(mass, inertia, nu, r_g)) <= G2_TOLERANCE, k
        assert _max_diff(out["C_RB"], source.C_RB_(nu)) <= G2_TOLERANCE, k
        assert _max_diff(out["C_RB"], Spheroid.C_RB_explicit(mass, inertia, nu, r_g)) <= G2_TOLERANCE, k
        m_rb, c_rb = Spheroid.get_RB_matrices(mass, inertia, nu, r_g)
        assert _max_diff(out["M"], m_rb + m_a) <= G2_TOLERANCE and _max_diff(out["C_RB"], c_rb) <= G2_TOLERANCE
        assert _max_diff(out["C_A"], AddedMass.C_A_6dof_lagrangian(nu, m_a)) <= G2_TOLERANCE, k
        expected = AddedMass.stabilize_C_A(AddedMass.C_A_6dof_lagrangian(nu, m_a))
        assert _max_diff(_evaluate(stabilized, nu, nu)["C_A"], expected) <= G2_TOLERANCE, k


def test_G2_spheroid_fossen_derivative_set_equals_lamb_without_roll():
    """Source get_added_mass_derivates (Fossen eqs. 8.74-8.77, unread) equals
    the Lamb form with r44 = 0: block derivatives with r44 = 0 = source set."""
    Spheroid = _source("matrices.rigid_body_kinetics_spheroid", "RigidBody6DOFSpheroid")
    AddedMass = _source("matrices.added_mass", "AddedMass")
    keywords = {**_remus_keywords(), "body_density": 1025.0, "water_density": 1025.0,
                "roll_added_inertia_ratio": 0.0}
    _, mass, _ = Spheroid(np.zeros(3)).get_I_gb(1025.0, keywords["semi_major_axis"], keywords["semi_minor_axis"])
    derivatives, _ = AddedMass.get_added_mass_derivates(mass, keywords["semi_major_axis"], keywords["semi_minor_axis"])
    constants, _ = _build("spheroid", **keywords)
    assert _max_diff(constants.added_mass_derivatives, derivatives) <= G2_TOLERANCE


@pytest.mark.parametrize("set_id", [1, 2])
def test_G2_displacement_hull_matches_source_functions(set_id):
    """get_params_hull mass, volume, wetted surface, I about the CG; the static
    M_RB / C_RB with I about the CG; the scaled derivatives with I about the CG."""
    SurfaceVessel = _source("matrices.rigid_body_kinetics_surface_vessel", "RigidBody6DOFSurfaceVessel")
    AddedMass = _source("matrices.added_mass", "AddedMass")
    hull = _hull()
    keywords = _hull_keywords(hull, hull.rows(set_id)[0])
    r_g = keywords["center_of_gravity"]
    _, ig_cg, mass, volume, wetted = SurfaceVessel().get_params_hull(
        keywords["water_density"], keywords["length"], keywords["beam"], keywords["draft"],
        keywords["block_coefficient"], keywords["radii_of_gyration"], r_g)
    derivatives = AddedMass().get_added_mass_derivates_hull(
        mass, keywords["length"], keywords["water_density"], ig_cg, keywords["added_mass_coefficients"])
    m_a = AddedMass.M_A_6dof(derivatives)
    constants, function = _build("displacement_hull", "co", False, **keywords)
    _, stabilized = _build("displacement_hull", "co", True, **keywords)
    assert abs(constants.mass - mass) <= G2_TOLERANCE * mass
    assert abs(constants.displaced_volume - volume) <= G2_TOLERANCE
    assert abs(constants.wetted_surface - wetted) <= G2_TOLERANCE
    assert _max_diff(constants.inertia, ig_cg) <= G2_TOLERANCE
    assert _max_diff(constants.rigid_body_mass_matrix, SurfaceVessel.M_RB(mass, ig_cg, r_g)) <= G2_TOLERANCE
    assert _max_diff(constants.added_mass_matrix, m_a) <= G2_TOLERANCE
    for k, nu in enumerate(_states()):
        out = _evaluate(function, nu, nu)
        assert _max_diff(out["C_RB"], SurfaceVessel.C_RB(mass, ig_cg, nu, r_g)) <= G2_TOLERANCE, k
        assert _max_diff(out["C_A"], AddedMass.C_A_6dof_lagrangian(nu, m_a)) <= G2_TOLERANCE, k
        expected = AddedMass.stabilize_C_A(AddedMass.C_A_6dof_lagrangian(nu, m_a))
        assert _max_diff(_evaluate(stabilized, nu, nu)["C_A"], expected) <= G2_TOLERANCE, k


def test_G2_displacement_hull_source_wiring_double_shift_is_asserted():
    """Source get_params_hull stores Ig_co = Ig_CG - m S(r_bg)^2 as I_gb (lines
    83, 91); C_RB_runtime then moves it with H a second time, while the hull
    plant builds M_RB from Ig_CG (asv_hull.py 98). The block follows rbody.m.
    Asserted here: the source's C_RB_runtime differs from the block for an
    offset CG and equals it for r_bg = 0; Newton-Euler (no model) agrees with
    the block, not with the source."""
    SurfaceVessel = _source("matrices.rigid_body_kinetics_surface_vessel", "RigidBody6DOFSurfaceVessel")
    hull = _hull()
    keywords = _hull_keywords(hull, hull.rows(2)[0])  # r_bG = [0.15, 0, -0.10]
    vessel = SurfaceVessel()
    ig_co, ig_cg, mass, _, _ = vessel.get_params_hull(
        keywords["water_density"], keywords["length"], keywords["beam"], keywords["draft"],
        keywords["block_coefficient"], keywords["radii_of_gyration"], keywords["center_of_gravity"])
    _, function = _build("displacement_hull", **keywords)
    nu = _states()[0]
    block_c = _evaluate(function, nu, nu)["C_RB"]
    source_c = vessel.C_RB_runtime(nu)
    assert _max_diff(block_c, source_c) > G4_FACTOR * G1_TOLERANCE
    exact = _newton_euler_wrench(mass, ig_cg, keywords["center_of_gravity"], nu)
    assert _max_diff(block_c @ nu, exact) <= 1e-9 * max(1.0, np.abs(exact).max())
    assert _max_diff(source_c @ nu, exact) > 1e-6 * np.abs(exact).max()
    zero_offset = {**keywords, "center_of_gravity": np.zeros(3)}
    vessel0 = SurfaceVessel()
    vessel0.get_params_hull(keywords["water_density"], keywords["length"], keywords["beam"], keywords["draft"],
                            keywords["block_coefficient"], keywords["radii_of_gyration"], np.zeros(3))
    _, function0 = _build("displacement_hull", **zero_offset)
    assert _max_diff(_evaluate(function0, nu, nu)["C_RB"], vessel0.C_RB_runtime(nu)) <= G2_TOLERANCE
    assert _max_diff(ig_co, ig_cg) > 0.0


def test_G2_planar_forms_match_source():
    """3-DOF (reduction of the 6-DOF block): M_A_3dof, C_A_3dof, RigidBody6DOF m2c 3-DOF branch."""
    RigidBody6DOF = _source("matrices.rigid_body_kinetics", "RigidBody6DOF")
    AddedMass = _source("matrices.added_mass", "AddedMass")
    block = _contract()
    dof3 = _dof3()
    states = _states()[:, [0, 1, 5]]
    for set_id in (1, 2):
        args = [dof3.column(name)[dof3.rows(set_id)[0]] for name in ("X_du", "Y_dv", "Y_dr", "N_dr")]
        mass_matrix = block.planar_added_mass_matrix(*args)
        assert _max_diff(mass_matrix, AddedMass.M_A_3dof(*args)) <= G2_TOLERANCE
        function = block.planar_coriolis_casadi(mass_matrix)
        for k, nu3 in enumerate(states):
            c = np.array(function(nu=nu3)["C"], dtype=float)
            assert _max_diff(c, AddedMass.C_A_3dof(args[0], args[1], args[2], *nu3)) <= G2_TOLERANCE, k
            assert _max_diff(c, RigidBody6DOF.get_C_RB_lagrangian(mass_matrix, nu3)) <= G2_TOLERANCE, k


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
        "explicit_mass_plus_1_percent": ("explicit", lambda kw: {**kw, "mass": kw["mass"] * 1.01}),
        "explicit_cg_x_plus_1_cm": ("explicit", lambda kw: {**kw, "center_of_gravity": kw["center_of_gravity"] + [0.01, 0, 0]}),
        "spheroid_body_density_plus_1_percent": ("spheroid", lambda kw: {**kw, "body_density": kw["body_density"] * 1.01}),
        "spheroid_r44_times_zero": ("spheroid", lambda kw: {**kw, "roll_added_inertia_ratio": 0.0}),
        "spheroid_water_density_equals_body": ("spheroid", lambda kw: {**kw, "water_density": kw["body_density"]}),
        "hull_block_coefficient_plus_1_percent": ("displacement_hull", lambda kw: {**kw, "block_coefficient": kw["block_coefficient"] * 1.01}),
        "hull_cg_z_plus_1_cm": ("displacement_hull", lambda kw: {**kw, "center_of_gravity": kw["center_of_gravity"] + [0, 0, 0.01]}),
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
        function = block.planar_coriolis_casadi(block.planar_added_mass_matrix(*perturbed))
        worst = max(_max_diff(np.array(function(nu=dof3.vector("nu3", 3)[k])["C"], float), dof3.matrix("C3", 3)[k])
                    for k in rows)
        assert (worst > G4_FACTOR * G1_TOLERANCE) is expect_detect, (perturbed, worst)


# --------------------------------------------------------------------------
# 3-DOF as a reduction of the 6-DOF block (owner, 2026-10-06): planar_casadi
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


def _planar_worst(planar, full, index):
    """Largest difference of M, C_RB nu3 and C_A nu_r3 of ``planar`` from the
    surge, sway, yaw part of the 6-DOF ``full`` at the embedded states."""
    worst = 0.0
    for nu3, nu_r3 in _planar_states():
        nu6, nu_r6 = _embed3(nu3), _embed3(nu_r3)
        six = _evaluate(full, nu6, nu_r6)
        out = planar(nu=nu3, nu_r=nu_r3)
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
    planar = block.planar_casadi(constants)
    assert planar.name_in() == ["nu", "nu_r"] and planar.name_out() == ["M", "C_RB", "C_A"]
    assert [planar.size_in(i) for i in range(2)] == [(3, 1)] * 2
    assert [planar.size_out(i) for i in range(3)] == [(3, 3)] * 3
    np.testing.assert_array_equal(block.planar_reduction(constants.total_mass_matrix),
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
    assert _planar_worst(block.planar_casadi(constants), full, index) <= G2_TOLERANCE, "control"
    key = "added_mass_derivatives" if form == "explicit" else "added_mass_coefficients"
    perturbed = {**keywords, key: np.asarray(keywords[key], float) * [1, 1.01, 1, 1, 1, 1]}
    perturbed_constants, _ = _build(form, coriolis, stabilize, **perturbed)
    worst = _planar_worst(block.planar_casadi(perturbed_constants), full, index)
    assert worst > G4_FACTOR * G1_TOLERANCE, (label, worst)


def test_physics_planar_coriolis_keeps_surge_coupling_of_a_general_mass_matrix():
    """A generally coupled symmetric positive-definite M3 (seeded, every
    off-diagonal entry non-zero): ``planar_coriolis_casadi(M3)`` equals the
    Kirchhoff form written here from the planar momentum p = M3 [u, v, r]
    (no model code), C3 = [[0, 0, -p_y], [0, 0, p_x], [p_y, -p_x, 0]]. The
    3-DOF branch of MSS ``m2c.m`` (lines 52-54 at ``72656d1``) keeps only
    M(1,1) u in p_x and M(2,2) v + M(2,3) r in p_y; it is transcribed below
    and must differ, so the test cannot pass on that branch."""
    block = _contract()
    rng = np.random.default_rng(PLANAR_SEED)
    for k in range(50):
        factor = rng.uniform(-1.0, 1.0, size=(3, 3))
        m3 = factor @ factor.T + 3.0 * np.eye(3)
        assert np.all(np.abs(m3[np.triu_indices(3, 1)]) > 0.0)
        function = block.planar_coriolis_casadi(m3)
        for nu3 in rng.uniform(-3.0, 3.0, size=(20, 3)):
            p_x, p_y, _ = m3 @ nu3
            kirchhoff = np.array([[0.0, 0.0, -p_y], [0.0, 0.0, p_x], [p_y, -p_x, 0.0]])
            got = np.array(function(nu=nu3)["C"], dtype=float)
            assert _max_diff(got, kirchhoff) <= G2_TOLERANCE * max(1.0, np.abs(kirchhoff).max()), k
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
    m, r_g = keywords["mass"], keywords["center_of_gravity"]
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
        mass=spheroid_constants.mass,
        inertia=spheroid_constants.inertia,
        center_of_gravity=spheroid_constants.center_of_gravity,
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
        mass=hull_constants.mass,
        inertia=hull_constants.inertia,
        center_of_gravity=hull_constants.center_of_gravity,
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
    _, function = _build("explicit", coriolis, mass=m, inertia=np.diag([2.0, 3.0, 4.0]),
                         center_of_gravity=[x_g, 0.0, 0.0], added_mass_derivatives=np.zeros(6))
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
    assert numpy_import, "no import from more_transformations (numpy, preprocessing)"
    assert casadi_import, "no import from more_transformations.more_casadi_transformations (graph)"
