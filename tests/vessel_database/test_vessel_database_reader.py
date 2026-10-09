"""Gate tests for the vessel-database (BEM) reader — **not yet built** (U10b).

Written 2026-10-09 (U10a), before the reader exists, so every "Reader contract"
test below fails with a named `pytest.fail` until U10b lands the module this
file imports. The "Fixture sanity" tests above them exercise only the frozen
data this job freezes (the copied `LAUV_marie.mat` and the MATLAB dump of it)
and the already-landed `more_transformations` helper; they pass today and stay
green through U10b, so a regression in the fixtures themselves is caught apart
from the reader's own correctness.

Why this fixture, why this scope (read first; see `SOURCE.md` beside this
file for every number's provenance)
----------------------------------------------------------------------------
* **File:** `data/LAUV_marie.mat`, copied byte-identical from MSS `cc07579`
  (`HYDRO/vessels_capytaine/LAUV_marie/LAUV_marie.mat`, MIT licence, sha256
  pinned in `SOURCE.md` and below), T. I. Fossen's idealised submerged body
  for NTNU AUR-Lab's LAUV *Marie* (E-37 Q2; A-35 Sec.1-2).
* **First intake only** (ADR 0004 Sec.2, accepted E-34 a): constant matrices
  — `A` at one frequency limit -> `M_A`, and the hydrostatic quantities a
  submerged body needs (`mass, water_density, gravity, center_of_gravity,
  center_of_buoyancy, displaced_volume`). `B(omega)` and the RAOs are **out of
  scope**; they enter later through a fluid-memory block of the hydrodynamic
  kind. No test below asks the reader to touch `B(omega)` beyond leaving it
  alone (`test_never_clips_or_rebuilds_damping`).
* **Placement** (E-53, the owner's verbatim answer *"then 2 how luka
  imagined it"*): the reader is a **script**, `scripts/vessel_database/
  vessel_database_reader.py` (repo-root `scripts/`, the same place as
  `scripts/thruster/`, `scripts/vehicle_models/`, `scripts/vehicles/`), not a
  `more_dynamics.models` block — ADR 0005 (all-CasADi blocks, no file I/O in
  `models/`) ruled the model-kind placement E-37 first proposed out (E-53
  Sec."Raised"). It reads the `.mat` **once** and writes plain numbers into a
  named vehicle part's `parameters.py`, with a provenance line in that part's
  `SOURCE.md` (E-53 Sec."Read as" item 1); **the model never opens a file**
  (rule 19, ADR 0005 — no numpy/file I/O inside a CasADi block). Import, for
  a repo-root script with no `__init__.py` (Python 3 implicit namespace
  package, the same as the three `scripts/<kind>/` folders already there):
  ``import scripts.vessel_database.vessel_database_reader``. This test file
  inserts the repo root onto `sys.path` itself (below) so it imports the same
  way whether pytest is invoked from the repo root (`python -m pytest tests
  -q`, the documented command) or from inside `tests/`.
* **Names below are this job's proposal (rule 19)**, not a landed contract:
  U10b may rename a function during porting, but must keep the refusals,
  the CO-shift math, the "C not handed on for a submerged body" rule and the
  "never touch B(omega)" rule enforced below, because each is an owner
  decision or a Fossen clarification, not a style choice.

Contract proposed for `scripts/vessel_database/vessel_database_reader.py`
(every symbol a "Reader contract" test below imports and exercises)
----------------------------------------------------------------------------
Exceptions (all subclass ``VesselDatabaseError``):
    ``VesselDatabaseError``, ``UnsupportedFileFormatError`` (v7.3/HDF5 .mat —
    A-35 Sec.4 Inputs, "a v7.3/HDF5 file is refused with a message"),
    ``MissingFieldError`` (A-35 Sec.4 Validation item 1),
    ``FrequencySentinelError`` (unknown omega or a missing 0/10 sentinel —
    Validation item 1, Inputs bullet 2), ``VelocityNotZeroError`` (no
    zero-speed column and none named — Inputs bullet 3), ``AxesMismatchError``
    (``hydrodynamic_axes`` != ``"MSS FSD"`` without a caller override).

``load_vessel_mat(path) -> VesselFile``
    ``scipy.io.loadmat(path, struct_as_record=False, squeeze_me=True)``,
    raising ``UnsupportedFileFormatError`` for a file scipy cannot load as a
    MATLAB v5 structure. ``VesselFile`` carries the raw fields untouched:
    ``freqs, velocities`` (1-D), ``A, B, C`` (6x6xnfreqxnvel), ``MRB`` (6x6),
    ``hydrodynamic_axes, hydrodynamic_reference, hydrodynamic_source`` (str),
    and ``main`` (mapping: ``m, rho, g, nabla, k44, k55, k66, GM_T, GM_L,
    Lpp, Lwl, T, B, C_B, submerged, submergenceDepth, CG, CB, name``).

``validate_vessel_file(vessel: VesselFile) -> None``
    A-35 Sec.4 "Validation". **Refuses** (raises, naming the field/value) on
    items 1-5: missing fields; wrong matrix shapes; `MRB`/`A` not symmetric
    (1e-8 relative); `MRB + A` not positive definite; `MRB` inconsistent with
    a CG reference when ``hydrodynamic_reference == "CG"``. **Warns**
    (``warnings.warn``, never raises) on item 6: a negative `B(omega)`
    diagonal entry; a `GM_T`/`GM_L` vs `C44/(rho g nabla)` mismatch beyond
    0.5%. Never silently repairs either.

``select_added_mass(vessel, *, frequency="infinite", velocity_index=None) -> np.ndarray``
    6x6, **at the file's own reference point** (no CO shift here — kept a
    separate step, A-35 Sec.4 Outputs). ``frequency``: ``"infinite"``
    (default, the last `freqs` entry when it equals 10, the MSS sentinel —
    Capytaine's `omega = infinity` solution, not a literal 10 rad/s point,
    `capytaine_vessel.py:183-186`), ``"zero"`` (the first entry when it
    equals 0), or an exact finite ``float`` matched against the `freqs` grid
    **with no interpolation** (A-35 Sec.4 Inputs bullet 2) — refuses by
    naming the requested value when absent. ``velocity_index``: ``None``
    auto-selects the one column where ``velocities == 0``; refuses by name
    when no zero-speed column exists and none is given explicitly (first
    intake is zero-speed only, Inputs bullet 3).

``shift_added_mass_to_co(matrix_at_reference, *, co_to_reference_point) -> np.ndarray``
    ``H(r)^T @ matrix_at_reference @ H(r)``, ``r = co_to_reference_point``
    (the vector from our body origin, CO, to the file's own reference point —
    A-35 Sec.4 Outputs: *"H(r)^T A_ref H(r) with r = CO -> reference point"*;
    `H` **only** from ``more_transformations.matrix_transforms.
    MatrixTransforms.H_matrix`` (Fossen 2011 eq. 3.24, p. 49) — never a local
    skew/H definition (Transforms rule, E-25/E-62). **The loader never
    guesses the origin** (A-35 Sec.4 Inputs bullet 4): the caller always
    states ``co_to_reference_point``; this function takes no default.

``hydrostatic_quantities(vessel) -> HydrostaticQuantities``
    ``mass, center_of_gravity, center_of_buoyancy`` (both CO-relative, so the
    caller's ``co_to_reference_point`` enters here too for the CG/CB fields),
    ``displaced_volume, water_density, gravity`` — exactly what
    ``preprocess_submerged_hydrostatics`` needs for ``W = m g``,
    ``B = rho g nabla`` (A-35 Sec.4 Outputs bullet 2). **Never** derived from
    ``GM_T``/``GM_L`` (X-5 reply 2: those are `C(4,4)/(m g)`, not
    `C(4,4)/(rho g nabla)` — a 0.16% difference that is not a fault, F-MSS-12
    closed as *explained*).

``select_restoring(vessel, *, body_kind) -> np.ndarray | None``
    ``body_kind="submerged"`` returns ``None`` — a submerged body's `C` is
    derived information that drops the net-buoyancy residual the
    hydrostatics block needs, so it is **not handed on** (ADR 0004 Sec.2;
    A-35 Sec.4 Outputs: *"For a submerged body the exported C is not handed
    on"*). ``body_kind="surface"`` returns ``vessel.C`` at the file's
    reference point, untouched — **never reconstructed from `GM_T`/`GM_L`**.

``read_bem_vessel(path, *, added_mass_frequency, body_kind, co_to_reference_point, generator_revision, velocity_index=None) -> BemVesselReading``
    Orchestrates the above into one frozen dataclass: ``added_mass_matrix``
    (at CO), ``restoring_matrix`` (``None`` for a submerged body),
    hydrostatics (the `HydrostaticQuantities` fields), ``rigid_body_mass_
    matrix_export`` (validation only — **not** shifted to CO, A-35 Sec.4
    Outputs: *"our M_RB comes from our mass properties"*), ``provenance``
    (a `BemProvenance`: ``kind="bem"``, ``generator`` = `vessel.
    hydrodynamic_source`, ``file_name``, ``sha256``, ``generator_revision``
    (git revision of the file's generating repository — the `.mat` itself
    carries none, so it is caller-supplied), ``frequency_limit``,
    ``reference_point``, ``axes``, ``submergence_depth``,
    ``co_to_reference_point``). No `damping_matrix` field anywhere on this
    dataclass — `B(omega)` is out of scope by construction, not merely
    unused.

``write_vessel_part_parameters(reading, *, parameters_path, source_md_path, part_name) -> None``
    E-53 placement a: writes the reading's numbers as plain Python literals
    into ``parameters_path`` (a `.rppws` part's ``parameters.py``, rpp's own
    form — no YAML, ADR 0004 addendum) and appends one `SOURCE.md` line per
    written value naming its `kind: bem`, the file, its sha256, the
    frequency limit and the CO shift used (A-35 Sec.4 "What the parameter
    record must say"). The function itself never reopens the `.mat` — the
    composition is self-contained afterwards (E-53's whole point).

Frozen reference: `data/vessel_database_lauv_marie_mss_reference.csv` (made
by `data/generate_vessel_database_mss.m`, run against MSS `cc07579`;
`data/SOURCE.md` names every column, origin, revision and cross-check).
"""

from __future__ import annotations

import csv
import hashlib
import importlib
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

# scripts/ has no __init__.py anywhere in this repo (an implicit namespace
# package, like scripts/thruster/, scripts/vehicle_models/, scripts/vehicles/
# already there); inserting the repo root here makes the import work the same
# way regardless of pytest's invocation directory (nothing relative to a
# machine: this is relative to this test file, inside the same repo, rule 9).
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

DATA_DIR = Path(__file__).resolve().parent / "data"
MAT_FILE = DATA_DIR / "LAUV_marie.mat"
REFERENCE_CSV = DATA_DIR / "vessel_database_lauv_marie_mss_reference.csv"
CONTRACT_MODULE = "scripts.vessel_database.vessel_database_reader"

# SOURCE.md pins the same two digests; both must change together if either
# frozen file is ever regenerated (rule 4: no number typed from memory for
# the *values*, but a sha256 is exactly a check that nothing was retyped).
MAT_FILE_SHA256 = "7eed3a0c4643ad446b6e45db9f74659052cc56866b3bc3a9c540bd5befcadb61"
REFERENCE_CSV_SHA256 = "7e38579daebcd940dcb754a7fcd4789799b40fe5d24e2cafa15b522721401426"

# A-66 case 2 table (30_checks/2026-10-09_A66_spheroid_auv_soundness.md),
# independently re-derived in 70_reviews/2026-10-09_A66_passB_codex.md: BEM
# A(infinity) at our CO, published to 6 significant figures.
A66_CASE2_BEM_AT_CO_DIAG = np.array(
    [0.650876, 33.9008, 33.8897, 2.13e-5, 9.76159, 9.76464]
)
A66_CASE2_TOLERANCE = 2e-4  # matches the ledger's own printed precision


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _reference_row() -> dict:
    with open(REFERENCE_CSV, newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 1, "one frozen vehicle, one row (SOURCE.md)"
    return rows[0]


def _matrix(row: dict, prefix: str) -> np.ndarray:
    """A 6x6 matrix from `<prefix>_01..<prefix>_36`, row-major (SOURCE.md)."""
    values = [float(row[f"{prefix}_{k:02d}"]) for k in range(1, 37)]
    return np.array(values, dtype=float).reshape(6, 6)


def _co_to_cg(row: dict) -> np.ndarray:
    """CO -> CG, `main.CG` as the file states it (never hand-typed)."""
    return np.array(
        [float(row["main_CG_1"]), float(row["main_CG_2"]), float(row["main_CG_3"])]
    )


def _h_matrix(r_bg: np.ndarray) -> np.ndarray:
    """`more_transformations`'s own H, not a local definition (Transforms rule)."""
    from more_transformations.matrix_transforms import MatrixTransforms

    return MatrixTransforms.H_matrix(r_bg)


def _contract():
    try:
        return importlib.import_module(CONTRACT_MODULE)
    except ModuleNotFoundError as exc:
        pytest.fail(
            "vessel-database reader not ported yet (U10b): "
            f"{exc}. Expected `{CONTRACT_MODULE}` "
            f"({CONTRACT_MODULE.replace('.', '/')}.py), see this file's module "
            "docstring for the proposed contract."
        )


# --------------------------------------------------------------------------
# Fixture sanity — exercise only this job's frozen data; green today.
# --------------------------------------------------------------------------
class TestFixtureSanity:
    def test_mat_file_is_the_credited_mss_copy(self):
        assert MAT_FILE.is_file(), "LAUV_marie.mat missing from tests/vessel_database/data/"
        assert _sha256(MAT_FILE) == MAT_FILE_SHA256, (
            "LAUV_marie.mat changed: re-freeze with a new revision note in SOURCE.md "
            "(rule 5 — supersede, never silently overwrite)"
        )

    def test_reference_csv_is_frozen(self):
        assert REFERENCE_CSV.is_file(), "run data/generate_vessel_database_mss.m (SOURCE.md)"
        assert _sha256(REFERENCE_CSV) == REFERENCE_CSV_SHA256

    def test_reference_csv_states_the_mss_sentinel_convention(self):
        row = _reference_row()
        assert row["hydrodynamic_axes"] == "MSS FSD"
        assert row["hydrodynamic_reference"] == "CG"
        assert int(row["n_freq"]) == 36
        assert int(row["freq_zero_idx"]) == 1
        assert int(row["freq_inf_idx"]) == 36
        assert float(row["freq_inf_value"]) == 10.0

    def test_restoring_matrix_is_constant_over_frequency(self):
        """A-35 Sec.3 item 5: checked once here so the reader's "C not handed
        on for a submerged body" shortcut (ADR 0004) is not hiding a
        frequency-dependence the generator would otherwise have caught."""
        row = _reference_row()
        assert int(row["c_varies_over_freq"]) == 0

    def test_co_shift_of_the_frozen_added_mass_matches_a66_case2(self):
        """Reproduce the published BEM-at-CO diagonal directly from this
        job's frozen CSV and `more_transformations`'s H, with no reader
        involved yet — the required cross-check against the independent
        shape-vs-formula finding (see this folder's `SOURCE.md`)."""
        row = _reference_row()
        a_inf_cg = _matrix(row, "A_inf")
        r = _co_to_cg(row)
        H = _h_matrix(r)
        a_inf_co = H.T @ a_inf_cg @ H
        np.testing.assert_allclose(
            np.diag(a_inf_co), A66_CASE2_BEM_AT_CO_DIAG, atol=A66_CASE2_TOLERANCE
        )

    def test_negative_damping_entries_are_present_and_unclipped_in_the_fixture(self):
        """X-5 reply 2 / item 12: three frequencies carry a small negative
        diagonal B(omega) entry, kept raw on purpose. Pinned here so a test
        of the reader (below) can assert it never "fixes" them."""
        row = _reference_row()
        freqs = [float(row[f"freq_{k:02d}"]) for k in range(1, 37)]
        for target in (3.2, 3.6, 5.5):
            idx = freqs.index(target) + 1  # 1-based, matches B_diag_f%02d
            diag = [float(row[f"B_diag_f{idx:02d}_dof{d}"]) for d in range(1, 7)]
            assert min(diag) < 0.0, f"expected a negative entry at {target} rad/s"


# --------------------------------------------------------------------------
# Reader contract — red until U10b lands `scripts/vessel_database/
# vessel_database_reader.py` (see the module docstring for the full contract).
# --------------------------------------------------------------------------
class TestReaderContract:
    def test_load_vessel_mat_exposes_frozen_fields_unmodified(self):
        reader = _contract()
        vessel = reader.load_vessel_mat(MAT_FILE)
        row = _reference_row()
        assert vessel.hydrodynamic_axes == row["hydrodynamic_axes"]
        assert vessel.hydrodynamic_reference == row["hydrodynamic_reference"]
        np.testing.assert_allclose(vessel.MRB, _matrix(row, "MRB"))
        assert vessel.freqs[0] == 0.0
        assert vessel.freqs[-1] == 10.0
        assert float(vessel.velocities) == 0.0 or 0.0 in np.atleast_1d(vessel.velocities)

    def test_load_vessel_mat_refuses_an_unsupported_file(self, tmp_path):
        reader = _contract()
        bogus = tmp_path / "not_a_matlab_file.mat"
        bogus.write_bytes(b"this is not a MATLAB v5 .mat file at all")
        with pytest.raises(reader.UnsupportedFileFormatError):
            reader.load_vessel_mat(bogus)

    def test_select_added_mass_infinite_matches_frozen(self):
        reader = _contract()
        vessel = reader.load_vessel_mat(MAT_FILE)
        row = _reference_row()
        a_inf = reader.select_added_mass(vessel, frequency="infinite")
        np.testing.assert_allclose(a_inf, _matrix(row, "A_inf"))

    def test_select_added_mass_zero_matches_frozen(self):
        reader = _contract()
        vessel = reader.load_vessel_mat(MAT_FILE)
        row = _reference_row()
        a_zero = reader.select_added_mass(vessel, frequency="zero")
        np.testing.assert_allclose(a_zero, _matrix(row, "A_zero"))

    def test_select_added_mass_refuses_an_unknown_frequency(self):
        reader = _contract()
        vessel = reader.load_vessel_mat(MAT_FILE)
        with pytest.raises(reader.FrequencySentinelError, match=r"3\.3"):
            reader.select_added_mass(vessel, frequency=3.3)  # not on the 34-point grid

    def test_select_added_mass_refuses_without_no_interpolation(self):
        """A-35 Sec.4 Inputs bullet 2: "no interpolation" — a value between
        two grid points is refused, not silently interpolated."""
        reader = _contract()
        vessel = reader.load_vessel_mat(MAT_FILE)
        with pytest.raises(reader.FrequencySentinelError):
            reader.select_added_mass(vessel, frequency=0.15)  # between 0.1 and 0.2

    def test_shift_to_co_matches_a66_case2(self):
        reader = _contract()
        vessel = reader.load_vessel_mat(MAT_FILE)
        row = _reference_row()
        a_inf_cg = reader.select_added_mass(vessel, frequency="infinite")
        a_inf_co = reader.shift_added_mass_to_co(
            a_inf_cg, co_to_reference_point=_co_to_cg(row)
        )
        np.testing.assert_allclose(
            np.diag(a_inf_co), A66_CASE2_BEM_AT_CO_DIAG, atol=A66_CASE2_TOLERANCE
        )

    def test_shift_to_co_requires_the_vector_explicitly(self):
        """"The loader never guesses the origin" (A-35 Sec.4 Inputs bullet 4)."""
        reader = _contract()
        import inspect

        params = inspect.signature(reader.shift_added_mass_to_co).parameters
        assert "co_to_reference_point" in params
        assert params["co_to_reference_point"].default is inspect.Parameter.empty

    def test_hydrostatic_quantities_match_frozen_main_fields(self):
        reader = _contract()
        vessel = reader.load_vessel_mat(MAT_FILE)
        row = _reference_row()
        hydro = reader.hydrostatic_quantities(vessel)
        assert hydro.mass == pytest.approx(float(row["main_m"]))
        assert hydro.water_density == pytest.approx(float(row["main_rho"]))
        assert hydro.gravity == pytest.approx(float(row["main_g"]))
        assert hydro.displaced_volume == pytest.approx(float(row["main_nabla"]))

    def test_hydrostatic_quantities_never_use_gm(self):
        """X-5 reply 2: GM_T/GM_L are C(4,4)/(m g), not the geometric
        CG-CB separation; never read when deriving mass/CG/CB/volume."""
        reader = _contract()
        import inspect

        source = inspect.getsource(reader.hydrostatic_quantities)
        assert "GM_T" not in source and "GM_L" not in source

    def test_restoring_matrix_not_exposed_for_a_submerged_body(self):
        reader = _contract()
        vessel = reader.load_vessel_mat(MAT_FILE)
        assert reader.select_restoring(vessel, body_kind="submerged") is None

    def test_rigid_body_mass_matrix_export_is_not_shifted(self):
        """A-35 Sec.4 Outputs: "validation only; our M_RB comes from our mass
        properties" — the export stays at the file's reference point."""
        reader = _contract()
        row = _reference_row()
        reading = reader.read_bem_vessel(
            MAT_FILE,
            added_mass_frequency="infinite",
            body_kind="submerged",
            co_to_reference_point=_co_to_cg(row),
            generator_revision="cc07579",
        )
        np.testing.assert_allclose(
            reading.rigid_body_mass_matrix_export, _matrix(row, "MRB")
        )

    def test_read_bem_vessel_has_no_damping_field(self):
        """B(omega) is out of scope by construction (ADR 0004 Sec.2), not
        merely unused — nothing on the reading can be mistaken for a damping
        matrix and later clipped by a caller."""
        reader = _contract()
        row = _reference_row()
        reading = reader.read_bem_vessel(
            MAT_FILE,
            added_mass_frequency="infinite",
            body_kind="submerged",
            co_to_reference_point=_co_to_cg(row),
            generator_revision="cc07579",
        )
        field_names = {f.lower() for f in vars(reading)}
        assert not any("damping" in f or f == "b" for f in field_names)

    def test_read_bem_vessel_provenance_names_everything_a35_requires(self):
        reader = _contract()
        row = _reference_row()
        reading = reader.read_bem_vessel(
            MAT_FILE,
            added_mass_frequency="infinite",
            body_kind="submerged",
            co_to_reference_point=_co_to_cg(row),
            generator_revision="cc07579",
        )
        provenance = reading.provenance
        assert provenance.kind == "bem"
        assert provenance.sha256 == MAT_FILE_SHA256
        assert provenance.generator_revision == "cc07579"
        assert provenance.frequency_limit in ("infinite", 10.0)
        assert provenance.reference_point == "CG"
        assert provenance.axes == "MSS FSD"
        assert provenance.submergence_depth == pytest.approx(float(row["main_submergenceDepth"]))

    def test_validate_warns_not_raises_on_the_known_negative_damping(self):
        reader = _contract()
        vessel = reader.load_vessel_mat(MAT_FILE)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            reader.validate_vessel_file(vessel)
        assert any("negative" in str(w.message).lower() for w in caught)

    def test_write_vessel_part_parameters_round_trips(self, tmp_path):
        reader = _contract()
        row = _reference_row()
        reading = reader.read_bem_vessel(
            MAT_FILE,
            added_mass_frequency="infinite",
            body_kind="submerged",
            co_to_reference_point=_co_to_cg(row),
            generator_revision="cc07579",
        )
        parameters_path = tmp_path / "parameters.py"
        source_md_path = tmp_path / "SOURCE.md"
        reader.write_vessel_part_parameters(
            reading,
            parameters_path=parameters_path,
            source_md_path=source_md_path,
            part_name="LAUV Marie (BEM)",
        )
        parameters_text = parameters_path.read_text()
        assert "class ComponentParameters" in parameters_text  # Luka's own form
        assert "bem" in source_md_path.read_text().lower()
        assert MAT_FILE_SHA256[:12] in source_md_path.read_text()

    def test_write_vessel_part_parameters_never_reopens_the_mat_file(self):
        """E-53's whole point: "the model never opens a file"."""
        reader = _contract()
        import inspect

        source = inspect.getsource(reader.write_vessel_part_parameters)
        assert "loadmat" not in source and ".mat" not in source
