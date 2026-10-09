"""Gate tests for the surface hydrostatics block (linear restoring ``G eta``).

Written 2026-10-06, before the block existed. This is **our** surface module;
the other developer's ``models/hydrostatics/linear_surface.py`` is a
different block, is not edited and is not imported here.

Gates (the test names carry them): G1 the block against MATLAB running MSS
(frozen CSVs), 1e-9 absolute; G2 against an independent transcription or the
block's own frozen form, 1e-10; G4 a perturbed model is detected; G5 a printed number of Fossen (2011) reproduced.

Contract of the block
---------------------
``more_dynamics.models.restoring.surface`` exports (all-CasADi block,
named parameters; geometry in, no vehicle numbers inside):

* ``surface_restoring_matrices(water_density, gravity, displacement_volume,
  waterplane_area, transverse_metacentric_height,
  longitudinal_metacentric_height, longitudinal_center_of_flotation,
  reference_point) -> (G, G_CO, G_CF)`` (CasADi, 6x6 each; numbers in, ``DM``
  out); MSS ``Gmtrx.m`` 29-37:
  ``G_CF = diag(0, 0, rho g A_wp, rho g nabla GMT, rho g nabla GML, 0)``,
  ``G_CO = H(r_bF)^T G_CF H(r_bF)`` with ``r_bF = [LCF 0 0]``,
  ``G = H(r_bP)^T G_CO H(r_bP)``; ``H`` from ``more_transformations`` ``H_matrix`` (``Hmtrx.m``).
* ``metacentric_heights(displacement_volume, center_of_buoyancy_above_keel,
  transverse_waterplane_inertia, longitudinal_waterplane_inertia,
  center_of_gravity, draft) -> 6x1 [BM_T, BM_L, KM_T, KM_L, GM_T,
  GM_L]`` (CasADi): ``BM = I/nabla``, ``KG = draft - z_g``, ``GM = KB + BM - KG``.
* ``surface_hydrostatics_parameters(*, hull_count, gravity_source="latitude",
  center_of_buoyancy_area="waterplane")`` declares ``length``, ``beam``,
  ``draft``, ``displacement_volume``, ``waterplane_area_coefficient``,
  ``center_of_gravity``, ``hull_lateral_offset`` (two hulls only),
  ``longitudinal_inertia_factor``, ``longitudinal_center_of_flotation``,
  ``reference_point``, ``water_density`` and ``latitude`` (degrees, gravity
  from ``more_transformations``) or ``gravity`` (``gravity_source="value"``,
  9.81 for the MSS checks, ``Gmtrx.m`` line 26).
* ``surface_hydrostatics_casadi(*, <the same selectors>) -> ca.Function``:
  inputs ``eta`` (6x1, NED, rad) and the declared parameters by name;
  outputs ``g = G eta``, ``G``, ``G_CO``, ``G_CF``, ``waterplane_area``,
  ``transverse_waterplane_inertia``, ``longitudinal_waterplane_inertia``,
  ``center_of_buoyancy_above_keel``, ``transverse_metacentric_height``,
  ``longitudinal_metacentric_height``, ``acceleration_of_gravity``. Per hull:
  ``A_hull = Cw L B``, Munro-Smith
  ``I_T,hull = (1/12) L B^3 6 Cw^3 / ((1 + Cw)(1 + 2 Cw))``. Whole craft:
  ``A_wp = n A_hull``, ``I_T = n I_T,hull + n A_hull y^2``,
  ``I_L = n c_L (1/12) B L^3``; Morrish
  ``KB = (1/3)(5 T/2 - (nabla/n) / A_KB)`` with ``A_KB = A_hull`` (default
  ``"waterplane"``) or ``A_KB = L B`` (``"length_times_beam"``, the
  ``otter.m`` behaviour up to MSS ``72656d1`` (line 177), which drops ``Cw``;
  MSS ``cf349d4`` changed ``otter.m`` to the waterplane area (line 178 at
  ``cc07579``), the default here, because a wall-sided hull then has KB = T/2
  exactly). Only the selectors ``gravity_source`` and
  ``center_of_buoyancy_area`` have defaults.

``g`` is the left-hand-side restoring vector, as in the submerged block: the
equation of motion is ``M nu_dot + ... + G eta = tau`` (``otter.m`` line 262 at
``cc07579``, ``... - G * eta``), so the applied force is ``-g``.

Conventions (hidden assumptions)
--------------------------------
* NED, z down; the CO on the waterline (``KG = T - z_g``, MSS
  ``exShipHydrostatics.m`` line 41 ``r_bB = [.. T-KB]``); small roll, pitch
  and heave (Fossen 2011, eqs. 4.22-4.24, pp. 64-65).
* ``LCF`` and ``reference_point`` are measured from the CO in BODY.
* Twin hulls are identical and placed at ``+-y``; the waterplane of each is
  symmetric about its own centreline.

Frozen reference: ``tests/restoring/data/`` (``SOURCE.md``, section of 2026-10-06).
"""

import ast
import importlib
import inspect
import os
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

DATA_DIR = Path(__file__).resolve().parent / "data"
HYDRODYNAMICS_DATA = Path(__file__).resolve().parents[1] / "hydrodynamic_loads" / "data"   # the catamaran reference of the hydrodynamic loads
# Outside this repo only through environment variables (nothing relative to one machine):
# MSS_DIR = an MSS checkout. Unset -> the tests that need it skip; the gates on the frozen
# CSVs in this folder's data/ run everywhere.
CONTRACT_MODULE = "more_dynamics.models.restoring.surface"
GMTRX_CSV = "surface_gmtrx_mss_current.csv"
CHAIN_CSV = "surface_chain_mss_current.csv"

G1_TOLERANCE = 1e-9   # G1: block vs MATLAB running MSS, absolute
G4_FACTOR = 10.0      # G4: a perturbed model must differ by more than this x G1
SEED = 20261006
N_RANDOM_STATES = 1000
N_RANDOM_SETS = 50
# Ship-scale G entries reach 1e9-1e10 N/m, N m/rad: an absolute 1e-9 is below
# double precision there. Those G are checked relative (not a gate), their
# length-scale coefficients (KB, BM, GM) at G1.
SHIP_G_RELATIVE = 1e-12

CHAIN_CSV_CC07579 = "surface_chain_mss_cc07579.csv"  # otter.m KB with the waterplane area (cf349d4)

DEFAULT_KB_AREA = "waterplane"          # Morrish with the hull's waterplane area (otter.m 178 from cf349d4)
MSS_OTTER_KB_AREA = "length_times_beam"  # otter.m line 177 up to MSS 72656d1 (behind the flag)

VEHICLE_NAMES = ("remus", "otter", "grethe", "marie", "hugin", "lauv",
                 "mariner", "torqeedo", "cybership", "prestero", "osv")

OTTER = "CRAFT/USV/models/otter.m"
EXSHIP = "mssExamples/exShipHydrostatics.m"
OSV = "CRAFT/SHIP/models/osv.m"
GMTRX = "LIBRARY/modeling/Gmtrx.m"
CITED_LINES = {
    (GMTRX, 1): "function G = Gmtrx(nabla,A_wp,GMT,GML,x_F,r_bP)",
    (GMTRX, 25): "rho = 1025;",
    (GMTRX, 26): "g   = 9.81;",
    (GMTRX, 29): "r_bF = [x_F, 0, 0]';",
    (GMTRX, 32): "G33_CF  = rho * g * A_wp;",
    (GMTRX, 33): "G44_CF  = rho * g * nabla * GMT;",
    (GMTRX, 34): "G55_CF  = rho * g * nabla * GML;",
    (GMTRX, 36): "G_CO = Hmtrx(r_bF)' * G_CF * Hmtrx(r_bF);",
    (GMTRX, 37): "G = Hmtrx(r_bP)' * G_CO * Hmtrx(r_bP);",
    ("LIBRARY/kinematics/Hmtrx.m", 17): "H = [eye(3)     S'",
    (OTTER, 90): "g   = 9.81;",
    (OTTER, 104): "B_pont  = 0.25;",
    (OTTER, 105): "y_pont  = 0.395;",
    (OTTER, 106): "Cw_pont = 0.75;",
    (OTTER, 174): "Aw_pont = Cw_pont * L * B_pont;",
    (OTTER, 175): "I_T = 2 * (1/12)*L*B_pont^3 * (6*Cw_pont^3/((1+Cw_pont)*(1+2*Cw_pont)))...",
    (OTTER, 176): "+ 2 * Aw_pont * y_pont^2;",
    (OTTER, 177): "I_L = 0.8 * 2 * (1/12) * B_pont * L^3;",
    (OTTER, 178): "KB = (1/3) * (5*T/2 - 0.5*nabla/Aw_pont);",
    (OTTER, 183): "KG = T - rg(3);",
    (OTTER, 187): "G33 = rho * g * (2 * Aw_pont);",
    (OTTER, 192): "LCF = -0.2;",
    (OTTER, 262): "M \\ ( tau + tau_damp + tau_crossflow - C * nu_r - G * eta )",
    (EXSHIP, 32): "Awp = Cw * B * L;",
    (EXSHIP, 35): "KB = (1/3) * (5*T/2 - nabla/Awp);",
    (EXSHIP, 38): "k_munro_smith =  (6 * Cw^3) / ((1+Cw) * (1+2*Cw));",
    (EXSHIP, 41): "r_bB = [-0.5 0 T-KB]';",
    (EXSHIP, 48): "I_L = 0.7 * (L^3 * B) / 12;",
    (EXSHIP, 55): "GM_T = BM_T - BG;",
    (OSV, 84): "vessel.KB = (1/3) * (5*vessel.T/2 - vessel.nabla/vessel.Awp);",
    (OSV, 92): "vessel.I_L = 0.7 * (vessel.L^3 * vessel.B) / 12;",
}


# --------------------------------------------------------------------------
# Helpers: MSS text, contract, data
# --------------------------------------------------------------------------
def _env_dir(variable):
    value = os.environ.get(variable)
    if not value:
        pytest.skip(f"{variable} is not set (nothing relative to one machine); set it to run this check")
    path = Path(value).expanduser()
    if not path.is_dir():
        pytest.skip(f"{variable}={value} is not a directory")
    return path


def _mss_line(rel, number):
    path = _env_dir("MSS_DIR") / rel
    if not path.exists():
        pytest.skip(f"{rel} not found under MSS_DIR")
    return path.read_text().splitlines()[number - 1].strip()


def _contract():
    try:
        return importlib.import_module(CONTRACT_MODULE)
    except ModuleNotFoundError as exc:
        if exc.name == "casadi":
            pytest.fail(f"casadi is not installed: {exc}")
        pytest.fail(f"block not ported yet: {exc}")


def _split(params):
    """Selectors and numbers of a parameter set (test-side): ``hull_count`` and
    ``center_of_buoyancy_area`` are selectors; ``gravity`` given means
    ``gravity_source="value"`` (``latitude`` then unused), otherwise
    ``"latitude"``; a single hull has no ``hull_lateral_offset`` parameter
    (the sets carry 0 there)."""
    values = dict(params)
    selectors = {
        "hull_count": int(values.pop("hull_count")),
        "center_of_buoyancy_area": values.pop("center_of_buoyancy_area", DEFAULT_KB_AREA),
        "gravity_source": "value" if "gravity" in values else "latitude",
    }
    if selectors["gravity_source"] == "value":
        values.pop("latitude", None)
    if selectors["hull_count"] == 1:
        assert values.pop("hull_lateral_offset", 0.0) == 0.0, "a single hull sits on the centreline"
    return selectors, values


def _build(params):
    """The block with the checked numbers, and its state-independent named
    outputs under the attribute names the assertions below use
    (``stiffness_matrix`` = ``G``, ``stiffness_at_origin`` = ``G_CO``,
    ``stiffness_at_flotation`` = ``G_CF``, ``gravity`` =
    ``acceleration_of_gravity``)."""
    from more_transformations.more_casadi_transformations import check_values

    block = _contract()
    selectors, numbers = _split(params)
    function = block.surface_hydrostatics_casadi(**selectors)
    values = check_values(block.surface_hydrostatics_parameters(**selectors), numbers)
    out = function(eta=np.zeros(6), **values)
    scalar = lambda name: float(out[name])
    constants = SimpleNamespace(
        gravity=scalar("acceleration_of_gravity"),
        waterplane_area=scalar("waterplane_area"),
        transverse_waterplane_inertia=scalar("transverse_waterplane_inertia"),
        longitudinal_waterplane_inertia=scalar("longitudinal_waterplane_inertia"),
        center_of_buoyancy_above_keel=scalar("center_of_buoyancy_above_keel"),
        transverse_metacentric_height=scalar("transverse_metacentric_height"),
        longitudinal_metacentric_height=scalar("longitudinal_metacentric_height"),
        stiffness_at_flotation=np.array(out["G_CF"], dtype=float),
        stiffness_at_origin=np.array(out["G_CO"], dtype=float),
        stiffness_matrix=np.array(out["G"], dtype=float),
    )
    return constants, (function, values)


def _evaluate(function, eta):
    block, values = function
    out = block(eta=np.asarray(eta, float), **values)
    return (np.array(out["g"], dtype=float).reshape(-1),
            np.array(out["G"], dtype=float))


def _load_csv(name):
    path = DATA_DIR / name
    header = path.read_text().splitlines()[0].split(",")
    data = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
    return [dict(zip(header, row)) for row in data]


def _g_matrix(row):
    return np.array([row[f"G_{i:02d}"] for i in range(1, 37)]).reshape(6, 6)  # row-major


def _max_diff(a, b):
    return float(np.max(np.abs(np.asarray(a, float) - np.asarray(b, float))))


def _h(r):
    """Hmtrx.m 16-18, test-side (the block must take it from more_transformations)."""
    x, y, z = r
    s = np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])
    return np.block([[np.eye(3), s.T], [np.zeros((3, 3)), np.eye(3)]])


def _gravity_library(latitude):
    from more_transformations.ecef_ned_transforms import ECEFNEDtransform
    return ECEFNEDtransform.gravity(latitude)


# --------------------------------------------------------------------------
# Parameter sets (geometry, no vehicle number in the block)
# --------------------------------------------------------------------------
def _chain_parameters(row, center_of_buoyancy_area):
    """Block parameters for one MSS chain row (generator columns)."""
    return {
        "length": row["L"],
        "beam": row["B_hull"],
        "draft": row["T"],
        "displacement_volume": row["nabla"],
        "waterplane_area_coefficient": row["Cw"],
        "center_of_gravity": np.array([0.0, 0.0, row["r_bg_3"]]),
        "hull_count": int(row["hull_count"]),
        "hull_lateral_offset": row["y_hull"],
        "longitudinal_inertia_factor": row["I_L_factor"],
        "longitudinal_center_of_flotation": row["x_F"],
        "reference_point": np.zeros(3),
        "water_density": row["rho"],
        "latitude": 63.446827,
        "gravity": row["g"],
        "center_of_buoyancy_area": center_of_buoyancy_area,
    }


def _twin_hull_parameters():
    """otter.m payload case 1 (mp = 25 kg), MSS behaviour of line 177 up to
    revision 72656d1 (the flag)."""
    row = _load_csv(CHAIN_CSV)[0]
    assert row["kind"] == 1 and row["mp"] == 25
    return _chain_parameters(row, MSS_OTTER_KB_AREA)


def _monohull_parameters():
    """A small monohull (test construction, not reference data): 5 m, 2 m
    beam, 0.4 m draft, Cb 0.5, Cw 0.8; CF aft, reference point off the CO."""
    length, beam, draft = 5.0, 2.0, 0.4
    return {
        "length": length,
        "beam": beam,
        "draft": draft,
        "displacement_volume": 0.5 * length * beam * draft,
        "waterplane_area_coefficient": 0.8,
        "center_of_gravity": np.array([-0.2, 0.0, -0.3]),
        "hull_count": 1,
        "hull_lateral_offset": 0.0,
        "longitudinal_inertia_factor": 0.7,
        "longitudinal_center_of_flotation": -0.3,
        "reference_point": np.array([0.4, 0.0, -0.2]),
        "water_density": 1025.0,
        "latitude": 63.446827,
    }


def _random_eta():
    rng = np.random.default_rng(SEED)
    eta = rng.uniform(-0.3, 0.3, size=(N_RANDOM_STATES, 6))
    eta[:, :3] = rng.uniform(-1.0, 1.0, size=(N_RANDOM_STATES, 3))
    eta[:, 5] = rng.uniform(-np.pi, np.pi, N_RANDOM_STATES)
    return eta


def _random_parameter_sets():
    """50 seeded small-craft sets, monohull and twin hull, both KB areas.
    Sizes keep |G| below ~1e5, where the frozen G2 tolerance (1e-10 absolute)
    is above double-precision rounding."""
    rng = np.random.default_rng(SEED + 1)
    sets = []
    for k in range(N_RANDOM_SETS):
        n = 1 if k % 2 == 0 else 2
        length = rng.uniform(1.5, 5.0)
        beam = rng.uniform(0.3, 1.5) if n == 1 else rng.uniform(0.15, 0.4)
        draft = rng.uniform(0.1, 0.4)
        cw = rng.uniform(0.6, 0.95)
        cb = rng.uniform(0.35, 0.95) * cw
        sets.append({
            "length": length,
            "beam": beam,
            "draft": draft,
            "displacement_volume": n * cb * length * beam * draft,
            "waterplane_area_coefficient": cw,
            "center_of_gravity": np.array([rng.uniform(-0.5, 0.5), 0.0, rng.uniform(-0.6, 0.1)]),
            "hull_count": n,
            "hull_lateral_offset": 0.0 if n == 1 else rng.uniform(0.3, 1.0),
            "longitudinal_inertia_factor": rng.uniform(0.6, 1.0),
            "longitudinal_center_of_flotation": rng.uniform(-0.5, 0.2),
            "reference_point": rng.uniform(-0.5, 0.5, 3),
            "water_density": rng.uniform(1000.0, 1030.0),
            "latitude": rng.uniform(-80.0, 80.0),
            "center_of_buoyancy_area": DEFAULT_KB_AREA if n == 1 or k % 4 == 1 else MSS_OTTER_KB_AREA,
        })
    return sets


def _morrish_kb(p, area):
    """Morrish KB (exShipHydrostatics.m 35 per hull; otter.m 177 for the flag)."""
    n = p["hull_count"]
    a_kb = (p["waterplane_area_coefficient"] * p["length"] * p["beam"]
            if area == DEFAULT_KB_AREA else p["length"] * p["beam"])
    return (2.5 * p["draft"] - p["displacement_volume"] / n / a_kb) / 3.0


# --------------------------------------------------------------------------
# Reference sanity and pinned lines (no port needed)
# --------------------------------------------------------------------------
def test_cited_mss_lines_are_unchanged():
    for (rel, number), text in CITED_LINES.items():
        line = _mss_line(rel, number)
        assert line.startswith(text), (rel, number, line)


def test_reference_shapes_and_kinds():
    gm = _load_csv(GMTRX_CSV)
    chain = _load_csv(CHAIN_CSV)
    assert len(gm) == 50 and len(chain) == 22
    assert [r["kind"] for r in chain] == [1.0] * 20 + [2.0, 3.0]
    assert all(r["rho"] == 1025.0 and r["g"] == 9.81 for r in chain)
    # off-origin reference points are present (the CO -> P transform is exercised)
    assert min(np.linalg.norm([r["r_bP_1"], r["r_bP_2"], r["r_bP_3"]]) for r in gm) > 0.0


def test_chain_case_1_equals_the_A26_catamaran_reference():
    """Same otter.m G as ``hydrodynamics/catamaran_matlab_reference_mss_current.csv``."""
    path = HYDRODYNAMICS_DATA / "catamaran_matlab_reference_mss_current.csv"
    header = path.read_text().splitlines()[0].split(",")
    a26 = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
    g26 = a26[0, [header.index(f"G_{i:02d}") for i in range(1, 37)]].reshape(6, 6)
    assert _max_diff(_g_matrix(_load_csv(CHAIN_CSV)[0]), g26) <= G1_TOLERANCE


def test_chain_case_1_equals_the_catamaran_reference_at_cc07579():
    """Both files made by MSS ``cc07579`` (``otter.m`` 178, the waterplane KB):
    chain row 1 G = ``catamaran_matlab_reference_mss_cc07579.csv`` G."""
    path = HYDRODYNAMICS_DATA / "catamaran_matlab_reference_mss_cc07579.csv"
    header = path.read_text().splitlines()[0].split(",")
    cat = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
    g_cat = cat[0, [header.index(f"G_{i:02d}") for i in range(1, 37)]].reshape(6, 6)
    assert _max_diff(_g_matrix(_load_csv(CHAIN_CSV_CC07579)[0]), g_cat) <= G1_TOLERANCE


def test_transcription_of_gmtrx_reproduces_matlab():
    """Gmtrx.m 29-37 transcribed (test side) equals the MATLAB CSV; checks the
    transcription before it is used as the MSS reference below."""
    for k, r in enumerate(_load_csv(GMTRX_CSV)):
        g_cf = 1025.0 * 9.81 * np.diag([0, 0, r["A_wp"], r["nabla"] * r["GMT"], r["nabla"] * r["GML"], 0])
        h_f, h_p = _h([r["x_F"], 0, 0]), _h([r["r_bP_1"], r["r_bP_2"], r["r_bP_3"]])
        assert _max_diff(h_p.T @ h_f.T @ g_cf @ h_f @ h_p, _g_matrix(r)) <= G1_TOLERANCE, k


@pytest.mark.parametrize("csv", [CHAIN_CSV, CHAIN_CSV_CC07579])
def test_transcription_of_chain_reproduces_matlab(csv):
    """The otter.m / exShipHydrostatics.m / osv.m chain transcribed equals
    MATLAB: ``otter.m`` 177 (``L B``) up to ``72656d1`` (``CHAIN_CSV``), the
    waterplane area from ``cc07579`` (``CHAIN_CSV_CC07579``, ``otter.m`` 178)."""
    for k, r in enumerate(_load_csv(csv)):
        area = MSS_OTTER_KB_AREA if r["kind"] == 1 and csv == CHAIN_CSV else DEFAULT_KB_AREA
        p = _chain_parameters(r, area)
        n, cw, length, beam = r["hull_count"], r["Cw"], r["L"], r["B_hull"]
        a_hull = cw * length * beam
        i_t = n * length * beam ** 3 / 12 * 6 * cw ** 3 / ((1 + cw) * (1 + 2 * cw)) + n * a_hull * r["y_hull"] ** 2
        i_l = n * r["I_L_factor"] * beam * length ** 3 / 12
        kb = _morrish_kb(p, area)
        gm_t = kb + i_t / r["nabla"] - (r["T"] - r["r_bg_3"])
        gm_l = kb + i_l / r["nabla"] - (r["T"] - r["r_bg_3"])
        for name, value in {"A_wp": n * a_hull, "I_T": i_t, "I_L": i_l, "KB": kb,
                            "GM_T": gm_t, "GM_L": gm_l}.items():
            assert abs(value - r[name]) <= G1_TOLERANCE, (k, name, value, r[name])


# --------------------------------------------------------------------------
# G1 — block vs current-MSS MATLAB
# --------------------------------------------------------------------------
def test_G1_restoring_matrices_match_gmtrx():
    """``surface_restoring_matrices`` vs MATLAB ``Gmtrx.m``, 50 cases with the
    reference point off the CO; ``G_CO``/``G_CF`` are checked by moving back."""
    block = _contract()
    for k, r in enumerate(_load_csv(GMTRX_CSV)):
        r_bp = np.array([r["r_bP_1"], r["r_bP_2"], r["r_bP_3"]])
        g, g_co, g_cf = (np.array(m, dtype=float) for m in block.surface_restoring_matrices(
            water_density=1025.0, gravity=9.81, displacement_volume=r["nabla"],
            waterplane_area=r["A_wp"], transverse_metacentric_height=r["GMT"],
            longitudinal_metacentric_height=r["GML"],
            longitudinal_center_of_flotation=r["x_F"], reference_point=r_bp,
        ))
        assert _max_diff(g, _g_matrix(r)) <= G1_TOLERANCE, k
        assert _max_diff(_h(r_bp).T @ g_co @ _h(r_bp), _g_matrix(r)) <= G1_TOLERANCE, k
        h_f = _h([r["x_F"], 0, 0])
        assert _max_diff(h_f.T @ g_cf @ h_f, g_co) <= G1_TOLERANCE, k
        assert _max_diff(g_cf, np.diag(np.diag(g_cf))) == 0.0, k


@pytest.mark.parametrize("kind, csv", [(1, CHAIN_CSV), (2, CHAIN_CSV), (3, CHAIN_CSV), (1, CHAIN_CSV_CC07579)])
def test_G1_coefficient_chain_matches_mss(kind, csv):
    """The block's coefficients vs MATLAB: otter.m (twin hull, 20 payloads;
    up to ``72656d1`` with the line-177 KB flag, ``CHAIN_CSV``; from
    ``cc07579`` with the default, ``CHAIN_CSV_CC07579``: the default is the G1
    of current MSS), exShipHydrostatics.m and osv.m (monohull, default).
    Length-scale terms at G1; ship-scale G relative (module note)."""
    for k, r in enumerate(_load_csv(csv)):
        if r["kind"] != kind:
            continue
        area = MSS_OTTER_KB_AREA if kind == 1 and csv == CHAIN_CSV else DEFAULT_KB_AREA
        constants, function = _build(_chain_parameters(r, area))
        assert constants.gravity == 9.81
        for field, column in {
            "waterplane_area": "A_wp",
            "transverse_waterplane_inertia": "I_T",
            "longitudinal_waterplane_inertia": "I_L",
            "center_of_buoyancy_above_keel": "KB",
            "transverse_metacentric_height": "GM_T",
            "longitudinal_metacentric_height": "GM_L",
        }.items():
            assert abs(getattr(constants, field) - r[column]) <= G1_TOLERANCE, (kind, k, field)
        g_ref = _g_matrix(r)
        _, g_out = _evaluate(function, np.zeros(6))
        if kind == 1:
            assert _max_diff(constants.stiffness_matrix, g_ref) <= G1_TOLERANCE, (kind, k)
            assert _max_diff(g_out, g_ref) <= G1_TOLERANCE, (kind, k)
        else:
            scale = np.max(np.abs(g_ref))
            assert _max_diff(constants.stiffness_matrix, g_ref) <= SHIP_G_RELATIVE * scale, (kind, k)
            assert _max_diff(g_out, g_ref) <= SHIP_G_RELATIVE * scale, (kind, k)


def test_G1_metacentric_heights_match_mss_chain():
    """``metacentric_heights`` on the MATLAB chain's own KB, I_T, I_L."""
    block = _contract()
    for k, r in enumerate(_load_csv(CHAIN_CSV)):
        out = block.metacentric_heights(
            displacement_volume=r["nabla"], center_of_buoyancy_above_keel=r["KB"],
            transverse_waterplane_inertia=r["I_T"], longitudinal_waterplane_inertia=r["I_L"],
            center_of_gravity=np.array([0.0, 0.0, r["r_bg_3"]]), draft=r["T"],
        )
        out = np.asarray(out, float)
        assert out.shape == (6, 1)
        out = out.reshape(-1)
        expected = [r["BM_T"], r["BM_L"], r["KB"] + r["BM_T"], r["KB"] + r["BM_L"], r["GM_T"], r["GM_L"]]
        assert _max_diff(out, expected) <= G1_TOLERANCE, k


# --------------------------------------------------------------------------
# G2 — block signature and frozen form
# --------------------------------------------------------------------------
def test_G2_function_signature():
    from more_transformations.more_casadi_transformations import freeze

    block = _contract()
    constants, (function, values) = _build(_monohull_parameters())
    declared = block.surface_hydrostatics_parameters(hull_count=1)
    assert function.name_in() == ["eta", *[d.name for d in declared]]
    assert function.name_out() == list(block.SURFACE_HYDROSTATICS_OUTPUTS)
    assert function.name_out()[:2] == ["g", "G"]
    frozen = freeze(function, declared, values)
    assert frozen.name_in() == ["eta"]
    assert frozen.size_in(0) == (6, 1)
    assert frozen.size_out("g") == (6, 1)
    assert frozen.size_out("G") == (6, 6)
    twin = block.surface_hydrostatics_parameters(hull_count=2, gravity_source="value")
    assert "hull_lateral_offset" in [d.name for d in twin] and "gravity" in [d.name for d in twin]
    assert "hull_lateral_offset" not in [d.name for d in declared] and "latitude" in [d.name for d in declared]
    for field in ("gravity", "waterplane_area", "transverse_waterplane_inertia",
                  "longitudinal_waterplane_inertia", "center_of_buoyancy_above_keel",
                  "transverse_metacentric_height", "longitudinal_metacentric_height",
                  "stiffness_at_flotation", "stiffness_at_origin", "stiffness_matrix"):
        assert hasattr(constants, field), field


# --------------------------------------------------------------------------
# G4 — a broken model must fail
# --------------------------------------------------------------------------
def _perturb_cw(p):
    return {**p, "waterplane_area_coefficient": 1.01 * p["waterplane_area_coefficient"]}


@pytest.mark.parametrize("kind, csv", [(1, CHAIN_CSV), (2, CHAIN_CSV), (1, CHAIN_CSV_CC07579)])
def test_G4_waterplane_coefficient_plus_1_percent_is_detected(kind, csv):
    row = next(r for r in _load_csv(csv) if r["kind"] == kind)
    area = MSS_OTTER_KB_AREA if kind == 1 and csv == CHAIN_CSV else DEFAULT_KB_AREA
    p = _chain_parameters(row, area)
    control, _ = _build(p)
    assert abs(control.transverse_metacentric_height - row["GM_T"]) <= G1_TOLERANCE, "control"
    broken, _ = _build(_perturb_cw(p))
    worst = max(abs(broken.waterplane_area - row["A_wp"]),
                abs(broken.transverse_metacentric_height - row["GM_T"]))
    assert worst > G4_FACTOR * G1_TOLERANCE, worst


def test_G4_perturbed_lcf_is_detected_by_gmtrx_reference():
    block = _contract()
    r = _load_csv(GMTRX_CSV)[0]
    g, _, _ = block.surface_restoring_matrices(
        water_density=1025.0, gravity=9.81, displacement_volume=r["nabla"],
        waterplane_area=r["A_wp"], transverse_metacentric_height=r["GMT"],
        longitudinal_metacentric_height=r["GML"],
        longitudinal_center_of_flotation=r["x_F"] + 0.01,
        reference_point=np.array([r["r_bP_1"], r["r_bP_2"], r["r_bP_3"]]))
    assert _max_diff(g, _g_matrix(r)) > G4_FACTOR * G1_TOLERANCE


# --------------------------------------------------------------------------
# G5 — published numbers (Fossen 2011, Example 4.2, p. 67; read in full)
# --------------------------------------------------------------------------
def _barge_parameters():
    """Fossen (2011), Example 4.2: barge 100 m x 8 m, draft 5 m, KG = 3.0 m,
    so with the CO on the waterline z_g = T - KG = 2.0 m; box hull: Cw = 1,
    I_L without reduction (eq. 4.38), nabla = 100 x 8 x 5 (eq. 4.39)."""
    return {
        "length": 100.0, "beam": 8.0, "draft": 5.0, "displacement_volume": 4000.0,
        "waterplane_area_coefficient": 1.0, "center_of_gravity": np.array([0.0, 0.0, 2.0]),
        "hull_count": 1, "hull_lateral_offset": 0.0, "longitudinal_inertia_factor": 1.0,
        "longitudinal_center_of_flotation": 0.0, "reference_point": np.zeros(3),
        "water_density": 1025.0, "latitude": 63.446827,
    }


def test_G5_fossen_2011_example_4_2_barge():
    """Printed: KB = 2.5 m, I_L = 666 666.7 m^4 (4.38), BM_L = 166.7 m (4.41),
    GM_L = 166.2 m (4.43), I_T = 4 266.7 m^4 (4.37). The Munro-Smith factor
    is 1 at Cw = 1, so I_T is the printed box value."""
    constants, _ = _build(_barge_parameters())
    assert round(constants.center_of_buoyancy_above_keel, 1) == 2.5
    assert round(constants.longitudinal_waterplane_inertia, 1) == 666666.7
    assert round(constants.transverse_waterplane_inertia, 1) == 4266.7
    assert round(constants.longitudinal_waterplane_inertia / 4000.0, 1) == 166.7
    assert round(constants.longitudinal_metacentric_height, 1) == 166.2


def test_G5_fossen_2011_example_4_2_printed_BM_T_is_a_misprint():
    """The book prints BM_T = 2.08 m (4.40) and GM_T = 1.58 m (4.42), but its
    own I_T = 4 266.7 m^4 (4.37) and nabla = 4 000 m^3 (4.39) give
    BM_T = 1.07 m and GM_T = 0.57 m (Fossen 2011, p. 67). Pinned so nobody
    "repairs" towards the printed value."""
    assert round(4266.7 / 4000.0, 2) == 1.07  # the book's own numbers
    constants, _ = _build(_barge_parameters())
    assert round(constants.transverse_metacentric_height, 2) == 0.57
    assert round(constants.transverse_metacentric_height, 2) != 1.58


# --------------------------------------------------------------------------
# Physical signs (no model needed; MSS is checked, not trusted)
# --------------------------------------------------------------------------
def _level_params(**changes):
    p = _monohull_parameters()
    p.update({"longitudinal_center_of_flotation": 0.0, "reference_point": np.zeros(3)})
    p.update(changes)
    return p


def test_physical_heave_down_gives_upward_restoring_force():
    """z > 0 (deeper, NED): extra displaced volume A_wp z, so the applied
    force is -rho g A_wp z (Fossen 2011 eq. 4.13-4.15), i.e. upward. The block
    returns ``g = +rho g A_wp z`` and the force on the craft is ``-g``."""
    p = _level_params()
    constants, function = _build(p)
    a_wp = p["waterplane_area_coefficient"] * p["length"] * p["beam"]
    for z in (0.01, 0.05, 0.1):
        g, _ = _evaluate(function, [0, 0, z, 0, 0, 0])
        expected_applied_z = -p["water_density"] * constants.gravity * a_wp * z
        assert abs(-g[2] - expected_applied_z) <= G1_TOLERANCE * abs(expected_applied_z), z
        assert -g[2] < 0.0, "upward in NED is negative z"
        assert _max_diff(np.delete(g, 2), 0.0) <= G1_TOLERANCE


def test_physical_roll_and_pitch_give_righting_moments():
    """GM > 0: a roll (pitch) angle produces an applied moment of the opposite
    sign (Fossen 2011 eq. 4.30-4.31, Definition 4.2)."""
    p = _level_params()
    constants, function = _build(p)
    assert constants.transverse_metacentric_height > 0.0
    assert constants.longitudinal_metacentric_height > 0.0
    for angle in (0.02, -0.05):
        g, _ = _evaluate(function, [0, 0, 0, angle, 0, 0])
        assert -g[3] * angle < 0.0, ("roll", angle)
        g, _ = _evaluate(function, [0, 0, 0, 0, angle, 0])
        assert -g[4] * angle < 0.0, ("pitch", angle)


def test_physical_heave_pitch_coupling_sign_with_flotation_centre_aft():
    """Derived from first principles, not from any model: a waterplane strip at
    x immerses by ``z - x theta`` (bow-up theta lifts the bow, z down), so the
    applied heave force is ``-rho g sum (z - x theta) dA`` and the applied
    pitch moment ``-sum x dF``. For a rectangular waterplane centred at
    ``x_F`` this gives ``G_CO[2,4] = G_CO[4,2] = -rho g A_wp x_F`` and
    ``G_CO[4,4] = rho g (nabla GM_L + A_wp x_F^2)``. (Fossen 2011 eq. 4.28 prints
    ``-Z_theta = +rho g int x dA``; the sign here follows the derivation and
    ``Gmtrx.m``.)"""
    x_f = -0.3
    p = _level_params(longitudinal_center_of_flotation=x_f)
    constants, _ = _build(p)
    rho_g = p["water_density"] * constants.gravity
    # strip integration of a rectangular waterplane of the same area
    a_wp = constants.waterplane_area
    xs = np.linspace(-0.5, 0.5, 2001)[:-1] + 0.5 / 2000
    xs = x_f + xs * p["length"]
    da = a_wp / xs.size
    force_per_theta = rho_g * np.sum(xs) * da           # dF_applied/dtheta
    moment_per_z = rho_g * np.sum(xs) * da               # dM_applied/dz = -sum x (-rho g dA)
    g_co = constants.stiffness_at_origin
    assert abs(g_co[2, 4] - (-force_per_theta)) <= 1e-6 * abs(force_per_theta)
    assert abs(g_co[4, 2] - (-moment_per_z)) <= 1e-6 * abs(moment_per_z)
    assert g_co[2, 4] > 0.0  # CF aft (x_F < 0): bow-up pitch immerses the stern more
    expected_55 = rho_g * (p["displacement_volume"] * constants.longitudinal_metacentric_height + a_wp * x_f ** 2)
    assert abs(g_co[4, 4] - expected_55) <= 1e-9 * expected_55


def test_physical_wall_sided_hull_has_its_centre_of_buoyancy_at_half_draft():
    """A wall-sided hull (vertical sides, Cb = Cw) is a prism of height T: its
    centroid is at T/2 above the keel, whatever the waterplane shape. Morrish
    with the hull's waterplane area gives exactly T/2 (the default, and
    ``otter.m`` 178 from MSS ``cf349d4``); the ``otter.m`` line-177 form up to
    ``72656d1`` (L B instead of A_wp) does not unless Cw = 1. The default must
    pass."""
    for n in (1, 2):
        for cw in (0.6, 0.75, 0.9, 1.0):
            p = {**_monohull_parameters(), "hull_count": n,
                 "hull_lateral_offset": 0.0 if n == 1 else 0.8,
                 "waterplane_area_coefficient": cw}
            p["displacement_volume"] = n * cw * p["length"] * p["beam"] * p["draft"]
            default, _ = _build(p)
            assert abs(default.center_of_buoyancy_above_keel - p["draft"] / 2) <= G1_TOLERANCE, (n, cw)
            mss, _ = _build({**p, "center_of_buoyancy_area": MSS_OTTER_KB_AREA})
            if cw < 1.0:
                assert mss.center_of_buoyancy_above_keel > p["draft"] / 2 + 1e-3, (n, cw)


def test_restoring_matrix_is_symmetric_at_every_reference_point():
    """Fossen 2011 eq. 4.26: G = G^T."""
    for p in _random_parameter_sets()[:10]:
        constants, _ = _build(p)
        for mat in (constants.stiffness_matrix, constants.stiffness_at_origin):
            assert _max_diff(mat, mat.T) <= G1_TOLERANCE


# --------------------------------------------------------------------------
# Morrish KB area: waterplane by default, MSS otter.m behaviour behind the flag
# --------------------------------------------------------------------------
def test_kb_default_differs_from_otter_line_177_by_the_documented_kb_term():
    """The default differs from ``otter.m`` before MSS ``cf349d4`` (line 177,
    up to ``72656d1``) by the KB term alone, and equals ``otter.m`` at
    ``cc07579`` (line 178)."""
    p = _twin_hull_parameters()
    flag, _ = _build(p)
    default, _ = _build({**p, "center_of_buoyancy_area": DEFAULT_KB_AREA})
    delta = _morrish_kb(p, DEFAULT_KB_AREA) - _morrish_kb(p, MSS_OTTER_KB_AREA)
    assert delta < 0.0  # Morrish with A_wp puts the CB lower for Cw < 1
    assert abs(default.transverse_metacentric_height
               - (flag.transverse_metacentric_height + delta)) <= G1_TOLERANCE
    assert abs(default.waterplane_area - flag.waterplane_area) == 0.0
    current = _load_csv(CHAIN_CSV_CC07579)[0]
    assert current["kind"] == 1 and current["mp"] == 25
    assert abs(default.center_of_buoyancy_above_keel - current["KB"]) <= G1_TOLERANCE
    assert abs(default.transverse_metacentric_height - current["GM_T"]) <= G1_TOLERANCE


def test_center_of_buoyancy_area_rejects_unknown_values():
    block = _contract()
    with pytest.raises(ValueError):
        block.surface_hydrostatics_casadi(hull_count=1, center_of_buoyancy_area="box")


def test_single_hull_must_be_on_the_centreline_and_hull_count_is_one_or_two():
    """A single hull has no lateral offset to give (the parameter is not
    declared, so a value is refused by name); 3 hulls is refused."""
    from more_transformations.more_casadi_transformations import check_values

    block = _contract()
    _, numbers = _split(_monohull_parameters())
    with pytest.raises(ValueError, match="hull_lateral_offset"):
        check_values(block.surface_hydrostatics_parameters(hull_count=1),
                     {**numbers, "hull_lateral_offset": 0.5})
    with pytest.raises(ValueError):
        block.surface_hydrostatics_casadi(hull_count=3)


# --------------------------------------------------------------------------
# Gravity and transforms from more_transformations
# --------------------------------------------------------------------------
def test_gravity_from_latitude_by_default_and_explicit_override():
    p = _monohull_parameters()
    constants, _ = _build(p)
    assert constants.gravity == _gravity_library(p["latitude"])
    overridden, _ = _build({**p, "gravity": 9.81})
    assert overridden.gravity == 9.81
    ratio = overridden.stiffness_at_flotation[2, 2] / constants.stiffness_at_flotation[2, 2]
    assert abs(ratio - 9.81 / constants.gravity) <= 1e-14


def _l0_classes():
    from more_transformations.ecef_ned_transforms import ECEFNEDtransform as GravityNumpy
    from more_transformations.matrix_transforms import MatrixTransforms as MatrixNumpy
    from more_transformations.more_casadi_transformations.ecef_ned_transforms import (
        ECEFNEDtransform as GravityCasadi)
    from more_transformations.more_casadi_transformations.matrix_transforms import (
        MatrixTransforms as MatrixCasadi)
    return (MatrixNumpy, MatrixCasadi), (GravityNumpy, GravityCasadi)


def test_transforms_and_gravity_are_taken_from_more_transformations(monkeypatch):
    """H (Hmtrx.m) and gravity (gravity.m) come from more_transformations, numpy
    or CasADi form; counted by wrapping those functions while the block is built."""
    matrices, gravities = _l0_classes()
    calls = {"H_matrix": 0, "gravity": 0}

    def wrap(cls, name):
        original = getattr(cls, name)

        def spy(*args, **kwargs):
            calls[name] += 1
            return original(*args, **kwargs)
        monkeypatch.setattr(cls, name, staticmethod(spy))

    for cls in matrices:
        wrap(cls, "H_matrix")
    for cls in gravities:
        wrap(cls, "gravity")
    _build(_monohull_parameters())
    assert calls["H_matrix"] >= 2, calls   # CF -> CO and CO -> P
    assert calls["gravity"] >= 1, calls


def _code_without_docstrings(module):
    tree = ast.parse(Path(module.__file__).read_text())
    for node in ast.walk(tree):
        body = getattr(node, "body", None)
        if isinstance(body, list) and body and isinstance(body[0], ast.Expr) \
                and isinstance(getattr(body[0], "value", None), ast.Constant) \
                and isinstance(body[0].value.value, str):
            body.pop(0)
    return tree


def test_transforms_no_local_copies_or_gravity_constants():
    block = _contract()
    tree = _code_without_docstrings(block)
    local_defs = {n.name.lower() for n in ast.walk(tree)
                  if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))}
    forbidden = {"_skew", "skew", "smtrx", "_h", "h_matrix", "hmtrx", "_transform", "gravity", "_gravity"}
    assert not (local_defs & forbidden), local_defs & forbidden
    literals = {n.value for n in ast.walk(tree) if isinstance(n, ast.Constant) and isinstance(n.value, float)}
    assert not literals & {9.81, 9.7803253359, 1025.0}, literals


# --------------------------------------------------------------------------
# Generic by construction (no vehicle name or number)
# --------------------------------------------------------------------------
def test_generic_block_reads_no_vehicle_name_and_has_no_vehicle_defaults():
    block = _contract()
    code = ast.unparse(_code_without_docstrings(block)).lower()
    assert not [n for n in VEHICLE_NAMES if n in code]
    signature = inspect.signature(block.surface_hydrostatics_casadi)
    defaults = {n: q.default for n, q in signature.parameters.items()
                if q.default is not inspect.Parameter.empty}
    assert defaults == {"gravity_source": "latitude", "center_of_buoyancy_area": DEFAULT_KB_AREA}, defaults
    for name in ("surface_restoring_matrices", "metacentric_heights"):
        sig = inspect.signature(getattr(block, name))
        assert all(q.default is inspect.Parameter.empty for q in sig.parameters.values()), name
    # two parameter sets give two different blocks from the same code
    _, a = _build(_monohull_parameters())
    _, b = _build(_twin_hull_parameters())
    eta = _random_eta()[0]
    assert _max_diff(_evaluate(a, eta)[0], _evaluate(b, eta)[0]) > G4_FACTOR * G1_TOLERANCE


# --------------------------------------------------------------------------
# The all-CasADi block: committed block, frozen path, gradient, no numpy
# --------------------------------------------------------------------------
COMMITTED_REVISION = "5d67caa"  # last revision of the block with numpy pre-processing
COMMITTED_FILES = ("__init__.py", "linear_surface.py", "submerged.py", "surface.py")
COMMITTED_SEED = 20261007
N_COMMITTED_STATES = 300
G2_TOLERANCE = 1e-10
GRADIENT_TOLERANCE = 1e-8
GRADIENT_STEP = 0.1  # m^3, central difference (g is quadratic in nabla: the difference is exact)
COEFFICIENTS = ("waterplane_area", "transverse_waterplane_inertia", "longitudinal_waterplane_inertia",
                "center_of_buoyancy_above_keel", "transverse_metacentric_height",
                "longitudinal_metacentric_height", "gravity", "stiffness_at_flotation",
                "stiffness_at_origin", "stiffness_matrix")


def _committed_package(tmp_path):
    """The hydrostatics package as committed at ``COMMITTED_REVISION``, read
    with ``git show`` into a temporary package (the old modules are imported,
    no values are frozen in a file). Skips when git or the revision is not
    available (a source archive without history)."""
    import subprocess

    root = Path(__file__).resolve().parents[2]
    package = tmp_path / "hydrostatics_committed_surface"
    package.mkdir()
    for name in COMMITTED_FILES:
        try:
            shown = subprocess.run(
                ["git", "-C", str(root), "show",
                 f"{COMMITTED_REVISION}:more_dynamics/models/hydrostatics/{name}"],
                capture_output=True, text=True, check=True,
            )
        except (OSError, subprocess.CalledProcessError) as exc:
            pytest.skip(f"revision {COMMITTED_REVISION} not readable with git here: {exc}")
        (package / name).write_text(shown.stdout)
    sys.path.insert(0, str(tmp_path))
    try:
        return importlib.import_module("hydrostatics_committed_surface.surface")
    finally:
        sys.path.remove(str(tmp_path))


def _committed_sets():
    """Every selector choice: monohull and twin hull, latitude and given
    gravity, both KB areas (the seeded small-craft sets, the two named sets)."""
    sets = {"monohull": _monohull_parameters(), "twin_hull_mss_flag": _twin_hull_parameters(),
            "monohull_flag": {**_monohull_parameters(), "center_of_buoyancy_area": MSS_OTTER_KB_AREA},
            "monohull_flag_value": {**_monohull_parameters(), "center_of_buoyancy_area": MSS_OTTER_KB_AREA,
                                    "gravity": 9.81}}
    for k, p in enumerate(_random_parameter_sets()[:12]):
        sets[f"random_{k}"] = p if k % 3 else {**p, "gravity": 9.81}
    return sets


def test_G2_block_equals_the_committed_block(tmp_path):
    """Every selector choice: g, G and the coefficients equal the block
    committed at ``COMMITTED_REVISION`` (numpy pre-processing) on 300 seeded
    states, 1e-10."""
    committed = _committed_package(tmp_path)
    rng = np.random.default_rng(COMMITTED_SEED)
    etas = rng.uniform(-0.3, 0.3, size=(N_COMMITTED_STATES, 6))
    seen = set()
    for name, p in _committed_sets().items():
        selectors, _ = _split(p)
        seen.add(tuple(sorted(selectors.items())))
        constants, function = _build(p)
        # Renamed to beam only from this round (A-50, shared names): the
        # committed revision still calls it hull_beam.
        old_p = {**{k: v for k, v in p.items() if k != "beam"}, "hull_beam": p["beam"]}
        old_constants = committed.preprocess_surface_hydrostatics(**old_p)
        old_function = committed.surface_hydrostatics_casadi(old_constants)
        for field in COEFFICIENTS:
            assert _max_diff(getattr(constants, field), getattr(old_constants, field)) <= G2_TOLERANCE, (name, field)
        for k, eta in enumerate(etas):
            g, stiffness = _evaluate(function, eta)
            old = old_function(eta=eta)
            assert _max_diff(g, np.array(old["g"], dtype=float).reshape(-1)) <= G2_TOLERANCE, (name, k)
            assert _max_diff(stiffness, np.array(old["G"], dtype=float)) <= G2_TOLERANCE, (name, k)
    assert len(seen) == 8, seen  # 2 hull counts x 2 gravity sources x 2 KB areas


def test_frozen_block_equals_the_block_called_with_numbers():
    """``freeze`` (what a plugin carries: a function of ``eta`` only) against
    the block called with the numbers, 300 seeded states; entry by entry
    within 1e-14 relative to max(1, |value|)."""
    from more_transformations.more_casadi_transformations import freeze

    block = _contract()
    for name, p in {"monohull": _monohull_parameters(), "twin_hull": _twin_hull_parameters()}.items():
        selectors, numbers = _split(p)
        _, (function, checked) = _build(p)
        frozen = freeze(function, block.surface_hydrostatics_parameters(**selectors), numbers)
        assert frozen.name_in() == ["eta"]
        rng = np.random.default_rng(COMMITTED_SEED)
        for k in range(N_COMMITTED_STATES):
            eta = rng.uniform(-0.3, 0.3, 6)
            a, b = frozen(eta=eta), function(eta=eta, **checked)
            for output in frozen.name_out():
                unfrozen = np.array(b[output], dtype=float)
                error = np.abs(np.array(a[output], dtype=float) - unfrozen) / np.maximum(1.0, np.abs(unfrozen))
                assert error.max() <= 1e-14, (name, output, k)


def test_gradient_of_g_with_respect_to_the_displacement_volume():
    """Identification path: the block called with ``displacement_volume``
    left as a symbol; dg/d nabla from CasADi equals a central difference of
    the block called with numbers, 20 seeded states (monohull, latitude
    gravity)."""
    import casadi as ca

    p = _monohull_parameters()
    selectors, numbers = _split(p)
    function = _contract().surface_hydrostatics_casadi(**selectors)
    symbol = ca.SX.sym("displacement_volume")
    nabla = numbers["displacement_volume"]
    rng = np.random.default_rng(COMMITTED_SEED)

    def g_of(volume, eta):
        return function(eta=eta, **{**numbers, "displacement_volume": volume})["g"]

    for k in range(20):
        eta = rng.uniform(-0.3, 0.3, 6)
        gradient = ca.Function("gradient", [symbol], [ca.jacobian(g_of(symbol, eta), symbol)])
        exact = np.array(gradient(nabla), dtype=float).ravel()
        upper = np.array(g_of(nabla + GRADIENT_STEP, eta), dtype=float).ravel()
        lower = np.array(g_of(nabla - GRADIENT_STEP, eta), dtype=float).ravel()
        difference = (upper - lower) / (2.0 * GRADIENT_STEP)
        assert _max_diff(exact, difference) <= GRADIENT_TOLERANCE, k
        assert np.abs(exact).max() > 1e-2, k  # the gradient is not trivially zero


def test_block_imports_no_numpy():
    """``surface.py`` imports neither numpy nor scipy, nor the numpy
    ``more_transformations`` modules (AST scan)."""
    for node in ast.walk(ast.parse(Path(_contract().__file__).read_text())):
        modules = []
        if isinstance(node, ast.Import):
            modules = [alias.name for alias in node.names]
        elif isinstance(node, ast.ImportFrom) and node.module:
            modules = [node.module]
        for module in modules:
            assert module.split(".")[0] not in ("numpy", "scipy"), module
            if module.startswith("more_transformations"):
                assert module.startswith("more_transformations.more_casadi_transformations"), module
