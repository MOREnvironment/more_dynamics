"""Vessel-database (BEM) reader: constant hydrodynamic quantities (added mass
at a frequency limit, the restoring matrix, hydrostatic terms) read once from
an MSS ``vessel`` structure (``HYDRO/wamit2vessel.m`` 50-94 field list) and
shifted to our body origin (CO) by the standard 6x6 transformation.

Scope (first intake): the constant matrices and scalars a BEM-derived
vehicle composition needs once, at load time -- ``A`` at one frequency limit
(``"infinite"`` or ``"zero"``, or an exact grid point; never interpolated
between grid points), ``M_RB``, and the hydrostatic quantities ``mass,
water_density, gravity, displaced_volume, center_of_gravity,
center_of_buoyancy``. The frequency-dependent radiation damping
``B(omega)`` and the response-amplitude operators are **out of scope**: a
fluid-memory block reads them later, directly from the same file; this
module never clips, zeroes or otherwise repairs a negative ``B(omega))``
diagonal entry (T. I. Fossen, personal communication, 2026-10-07: mesh
noise, negligible once spectrum-weighted, not to be removed by hand) and
carries no ``damping``/``B`` field on its result for exactly that reason --
nothing downstream can mistake an absent value for a zeroed one.

The file's own hydrodynamic reference point need not be the body's CO:
``select_added_mass`` returns the matrix exactly as stored, and
``shift_added_mass_to_co`` moves it with the standard transformation
(``H(r)^T A_ref H(r)``, Fossen 2011 eq. 3.24 p. 49) given the vector from CO
to that reference point explicitly -- never a default, since a wrong origin
would silently corrupt every downstream dynamics term. A submerged body's
restoring matrix ``C`` is hydrostatic-equilibrium information a submerged-
body preprocessor does not read (it derives its own net-buoyancy residual
from ``mass``, ``water_density``, ``gravity`` and ``displaced_volume``
instead): ``select_restoring`` returns ``None`` for
``body_kind="submerged"`` rather than handing on a matrix nothing should
read. The hydrostatic quantities are read from the file's mass and volume
fields directly, never reconstructed from its metacentric heights ``GM_T``/
``GM_L`` (T. I. Fossen, personal communication, 2026-10-07: those equal
``C(4,4)/(m g)``, a quantity distinct from ``C(4,4)/(rho g nabla)`` by
construction, not a data fault).

``write_vessel_part_parameters`` writes the result's numbers as plain
Python literals into a composition part's parameters file (the form the
registry's own editor writes and reads: one class attribute per value, no
nested structure) and appends one provenance line per value to a companion
file beside it; once written, the composition is self-contained and this
module is not read again to build or run it.

Equations
---------
* Reference-point shift: ``M_CO = H(r)^T M_ref H(r)``, ``r`` the vector from
  CO to the file's own reference point, ``H`` built from the 3x3 skew of
  ``r`` (Fossen 2011, eq. 3.24, p. 49; MSS ``Hmtrx.m`` 16-18).
* Hydrostatic restoring sanity (validation only, never used to derive a
  value): ``C(4,4) == rho g nabla (z_G - z_B)`` at the file's own reference
  point (Fossen 2011, eq. 4.9 and Table 4.1, Ch. 4).

References
----------
[Fossen 2011] Fossen, T. I. (2011). *Handbook of Marine Craft Hydrodynamics
    and Motion Control*, 1st ed. John Wiley & Sons, Chichester. Ch. 3,
    eq. 3.24, p. 49; Ch. 4, eq. 4.9 and Table 4.1, pp. 59-80.
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2.
    https://github.com/cybergalactic/MSS, MIT licence, revision ``cc07579``:
    ``HYDRO/wamit2vessel.m`` 50-94 (the ``vessel`` structure's field list),
    ``HYDRO/computeManeuveringModel.m`` (frequency/velocity grid use).
[MSS-Capytaine] Fossen, T. I. *MSS-Capytaine*, release 1.0.
    https://github.com/cybergalactic/MSS-Capytaine, revision ``bfc9dc9``:
    ``src/capytaine_vessel.py`` 183-186, 239-241 (the ``[0, ..., 10]``
    frequency-sentinel convention: the ``10`` rad/s grid point holds the
    zero-speed infinite-frequency radiation solution, not a literal result
    at 10 rad/s).

Author:    Enio Krizman
Date:      2026-10-09
"""

from __future__ import annotations

import hashlib
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import scipy.io

from more_transformations.matrix_transforms import MatrixTransforms

_MSS_FREQUENCY_SENTINEL = 10.0  # capytaine_vessel.py:183-186,239-241


# ---------------------------------------------------------------------------
# Exceptions
# ---------------------------------------------------------------------------
class VesselDatabaseError(Exception):
    """Base class for every refusal this module raises."""


class UnsupportedFileFormatError(VesselDatabaseError):
    """The file is not readable as a MATLAB v5 ``vessel`` structure."""


class MissingFieldError(VesselDatabaseError):
    """A required field is absent, or has the wrong shape."""


class FrequencySentinelError(VesselDatabaseError):
    """The requested frequency is not on the file's own grid (no interpolation)."""


class VelocityNotZeroError(VesselDatabaseError):
    """No zero-speed column exists, and none was given explicitly."""


class AxesMismatchError(VesselDatabaseError):
    """``hydrodynamic_axes`` is not ``"MSS FSD"`` and no override was given."""


# ---------------------------------------------------------------------------
# Data carried from the file, untouched
# ---------------------------------------------------------------------------
@dataclass
class MainFields:
    m: Optional[float] = None
    rho: Optional[float] = None
    g: Optional[float] = None
    nabla: Optional[float] = None
    k44: Optional[float] = None
    k55: Optional[float] = None
    k66: Optional[float] = None
    GM_T: Optional[float] = None
    GM_L: Optional[float] = None
    Lpp: Optional[float] = None
    Lwl: Optional[float] = None
    T: Optional[float] = None
    B: Optional[float] = None
    C_B: Optional[float] = None
    submerged: Optional[float] = None
    submergenceDepth: Optional[float] = None
    CG: Optional[np.ndarray] = None
    CB: Optional[np.ndarray] = None
    name: Optional[str] = None


@dataclass
class VesselFile:
    """The MSS ``vessel`` structure (``wamit2vessel.m`` 50-94), read as is.

    A field absent from the source structure is carried through as `None`
    rather than guessed, for `validate_vessel_file` to refuse by name.
    """

    freqs: Optional[np.ndarray]
    velocities: Optional[np.ndarray]
    A: Optional[np.ndarray]
    B: Optional[np.ndarray]
    C: Optional[np.ndarray]
    MRB: Optional[np.ndarray]
    hydrodynamic_axes: Optional[str]
    hydrodynamic_reference: Optional[str]
    hydrodynamic_source: Optional[str]
    main: Optional[MainFields]


@dataclass
class HydrostaticQuantities:
    mass: float
    center_of_gravity: np.ndarray
    center_of_buoyancy: np.ndarray
    displaced_volume: float
    water_density: float
    gravity: float


@dataclass
class BemProvenance:
    kind: str
    generator: str
    file_name: str
    sha256: str
    generator_revision: str
    frequency_limit: object
    reference_point: str
    axes: str
    submergence_depth: float
    co_to_reference_point: np.ndarray


@dataclass
class BemVesselReading:
    """One frozen reading. No `damping`/`B` field anywhere: `B(omega)` is
    out of scope by construction, not merely unused."""

    added_mass_matrix: np.ndarray
    restoring_matrix: Optional[np.ndarray]
    hydrostatics: HydrostaticQuantities
    rigid_body_mass_matrix_export: np.ndarray
    provenance: BemProvenance


# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
def _as_float_array(value):
    return None if value is None else np.asarray(value, dtype=float)


def _read_main_fields(main_struct) -> MainFields:
    def scalar(name):
        value = getattr(main_struct, name, None)
        return None if value is None else float(value)

    def vector(name):
        value = getattr(main_struct, name, None)
        return None if value is None else np.asarray(value, dtype=float)

    return MainFields(
        m=scalar("m"),
        rho=scalar("rho"),
        g=scalar("g"),
        nabla=scalar("nabla"),
        k44=scalar("k44"),
        k55=scalar("k55"),
        k66=scalar("k66"),
        GM_T=scalar("GM_T"),
        GM_L=scalar("GM_L"),
        Lpp=scalar("Lpp"),
        Lwl=scalar("Lwl"),
        T=scalar("T"),
        B=scalar("B"),
        C_B=scalar("C_B"),
        submerged=scalar("submerged"),
        submergenceDepth=scalar("submergenceDepth"),
        CG=vector("CG"),
        CB=vector("CB"),
        name=getattr(main_struct, "name", None),
    )


def load_vessel_mat(path) -> VesselFile:
    """Read an MSS ``vessel`` structure (``wamit2vessel.m`` 50-94). Refuses
    a v7.3/HDF5 file, or anything else `scipy` cannot read as a MATLAB v5
    structure, by name."""
    path = Path(path)
    try:
        raw = scipy.io.loadmat(str(path), struct_as_record=False, squeeze_me=True)
    except Exception as exc:  # scipy raises several distinct error types
        raise UnsupportedFileFormatError(
            f"{path}: not readable as a MATLAB v5 vessel-database file ({exc})"
        ) from exc

    vessel_struct = raw.get("vessel")
    if vessel_struct is None:
        raise MissingFieldError(
            f"{path}: no top-level 'vessel' structure (wamit2vessel.m field list)"
        )

    main_struct = getattr(vessel_struct, "main", None)
    main = _read_main_fields(main_struct) if main_struct is not None else None

    return VesselFile(
        freqs=_as_float_array(getattr(vessel_struct, "freqs", None)),
        velocities=_as_float_array(getattr(vessel_struct, "velocities", None)),
        A=_as_float_array(getattr(vessel_struct, "A", None)),
        B=_as_float_array(getattr(vessel_struct, "B", None)),
        C=_as_float_array(getattr(vessel_struct, "C", None)),
        MRB=_as_float_array(getattr(vessel_struct, "MRB", None)),
        hydrodynamic_axes=getattr(vessel_struct, "hydrodynamic_axes", None),
        hydrodynamic_reference=getattr(vessel_struct, "hydrodynamic_reference", None),
        hydrodynamic_source=getattr(vessel_struct, "hydrodynamic_source", None),
        main=main,
    )


# ---------------------------------------------------------------------------
# Validate
# ---------------------------------------------------------------------------
def _require(condition: bool, message: str) -> None:
    if not condition:
        raise MissingFieldError(message)


def _warn_on_damping_and_metacentric_height(vessel: VesselFile) -> None:
    B = vessel.B
    if B is not None:
        diagonal = (
            np.stack([B[i, i, :] for i in range(6)])
            if B.ndim == 3
            else np.stack([B[i, i, :, 0] for i in range(6)])
        )
        if np.any(diagonal < 0.0):
            warnings.warn(
                "vessel.B has a negative diagonal entry (mesh noise, kept "
                "raw on purpose; see this module's docstring)",
                stacklevel=2,
            )

    main = vessel.main
    C = vessel.C
    if main is not None and C is not None and main.GM_T and main.nabla and main.rho and main.g:
        c44 = C[3, 3, 0] if C.ndim == 3 else C[3, 3, 0, 0]
        roll_restoring = main.rho * main.g * main.nabla * main.GM_T  # Fossen 2011, Table 4.1
        relative_difference = abs(c44 - roll_restoring) / abs(roll_restoring)
        if relative_difference > 0.005:
            warnings.warn(
                "vessel.main.GM_T does not match C(4,4)/(rho g nabla) to "
                "within 0.5% (distinct quantities by construction; see "
                "this module's docstring)",
                stacklevel=2,
            )


def validate_vessel_file(vessel: VesselFile) -> None:
    """Refuse a structurally unsound file; warn, never repair, a known
    physically-explained peculiarity.

    Raises (naming the field) on: a missing required field; a matrix of the
    wrong shape; `M_RB` or `A` not symmetric to 1e-8 relative tolerance;
    `M_RB + A` not positive definite at the zero- and infinite-frequency
    limits; `M_RB` inconsistent with a CG reference (a non-zero
    mass-inertia coupling block) when ``hydrodynamic_reference == "CG"``.
    Warns, never raises or repairs, on: a negative `B(omega)` diagonal
    entry; a `GM_T`/`GM_L` vs `C(4,4)/(rho g nabla)` mismatch beyond 0.5%.
    """
    _require(vessel.MRB is not None, "vessel.MRB is missing")
    _require(vessel.A is not None, "vessel.A is missing")
    _require(vessel.B is not None, "vessel.B is missing")
    _require(vessel.C is not None, "vessel.C is missing")
    _require(vessel.freqs is not None, "vessel.freqs is missing")
    _require(vessel.velocities is not None, "vessel.velocities is missing")
    _require(vessel.hydrodynamic_axes is not None, "vessel.hydrodynamic_axes is missing")
    _require(
        vessel.hydrodynamic_reference is not None, "vessel.hydrodynamic_reference is missing"
    )
    _require(vessel.main is not None, "vessel.main is missing")

    main = vessel.main
    for name in ("m", "rho", "g", "nabla", "CG", "CB", "submergenceDepth"):
        _require(getattr(main, name) is not None, f"vessel.main.{name} is missing")

    MRB = vessel.MRB
    A = vessel.A
    _require(MRB.shape == (6, 6), f"vessel.MRB has shape {MRB.shape}, expected (6, 6)")
    _require(
        A.ndim in (3, 4) and A.shape[:2] == (6, 6),
        f"vessel.A has shape {A.shape}, expected (6, 6, n_freq[, n_vel])",
    )

    if not np.allclose(MRB, MRB.T, rtol=1e-8, atol=1e-9):
        raise VesselDatabaseError("vessel.MRB is not symmetric to 1e-8 relative tolerance")

    if vessel.hydrodynamic_reference == "CG":
        coupling = MRB[0:3, 3:6]  # zero at the CG reference (Fossen 2011, eq. 3.26, p. 50)
        scale = max(float(np.max(np.abs(np.diag(MRB)[0:3]))), 1.0)
        if np.max(np.abs(coupling)) > 1e-6 * scale:
            raise VesselDatabaseError(
                "vessel.MRB has a non-zero mass-inertia coupling block, "
                "inconsistent with hydrodynamic_reference == 'CG'"
            )

    for freq_index in sorted({0, A.shape[2] - 1}):
        a_slice = A[:, :, freq_index] if A.ndim == 3 else A[:, :, freq_index, 0]
        if not np.allclose(a_slice, a_slice.T, rtol=1e-8, atol=1e-9):
            raise VesselDatabaseError(
                f"vessel.A is not symmetric to 1e-8 relative tolerance at "
                f"frequency index {freq_index}"
            )
        eigenvalues = np.linalg.eigvalsh(MRB + a_slice)
        if np.min(eigenvalues) <= 0.0:
            raise VesselDatabaseError(
                f"vessel.MRB + vessel.A is not positive definite at frequency index {freq_index}"
            )

    _warn_on_damping_and_metacentric_height(vessel)


# ---------------------------------------------------------------------------
# Select
# ---------------------------------------------------------------------------
def _resolve_velocity_index(vessel: VesselFile, velocity_index: Optional[int]) -> int:
    if velocity_index is not None:
        return int(velocity_index)
    velocities = np.atleast_1d(vessel.velocities)
    zero_indices = np.flatnonzero(velocities == 0.0)
    if zero_indices.size == 0:
        raise VelocityNotZeroError(
            "vessel.velocities has no zero-speed column; pass velocity_index explicitly"
        )
    return int(zero_indices[0])


def _resolve_frequency_index(vessel: VesselFile, frequency) -> int:
    freqs = vessel.freqs
    if frequency == "infinite":
        if freqs[-1] != _MSS_FREQUENCY_SENTINEL:
            raise FrequencySentinelError(
                f"vessel.freqs does not end with the {_MSS_FREQUENCY_SENTINEL} "
                "rad/s infinite-frequency sentinel"
            )
        return len(freqs) - 1
    if frequency == "zero":
        if freqs[0] != 0.0:
            raise FrequencySentinelError("vessel.freqs does not start with the 0 rad/s sentinel")
        return 0
    value = float(frequency)
    matches = np.flatnonzero(np.isclose(freqs, value, rtol=1e-9, atol=1e-12))
    if matches.size == 0:
        raise FrequencySentinelError(
            f"requested frequency {value} rad/s is not on the vessel's frequency grid "
            "(no interpolation)"
        )
    return int(matches[0])


def _slice_matrix(array: np.ndarray, freq_index: int, vel_index: int) -> np.ndarray:
    if array.ndim == 3:
        return np.array(array[:, :, freq_index], dtype=float)
    if array.ndim == 4:
        return np.array(array[:, :, freq_index, vel_index], dtype=float)
    raise VesselDatabaseError(f"unexpected matrix rank {array.ndim}, expected 3 or 4")


def select_added_mass(
    vessel: VesselFile, *, frequency="infinite", velocity_index=None
) -> np.ndarray:
    """6x6 added-mass matrix at the file's own reference point, at one
    frequency limit or an exact grid point (never interpolated)."""
    vel_index = _resolve_velocity_index(vessel, velocity_index)
    freq_index = _resolve_frequency_index(vessel, frequency)
    return _slice_matrix(vessel.A, freq_index, vel_index)


def shift_added_mass_to_co(matrix_at_reference: np.ndarray, *, co_to_reference_point) -> np.ndarray:
    """``H(r)^T M_ref H(r)``, ``r`` the vector from CO to the file's own
    reference point (Fossen 2011, eq. 3.24, p. 49). The caller always
    states ``r``; a wrong origin would silently corrupt every downstream
    term, so this function takes no default."""
    H = MatrixTransforms.H_matrix(co_to_reference_point)  # Fossen 2011, eq. 3.24, p. 49
    return H.T @ np.asarray(matrix_at_reference, dtype=float) @ H


def hydrostatic_quantities(vessel: VesselFile) -> HydrostaticQuantities:
    """``mass, center_of_gravity, center_of_buoyancy, displaced_volume,
    water_density, gravity`` read directly from the file's mass and volume
    fields. Never derived from the file's metacentric heights: those equal
    ``C(4,4)/(m g)``, a quantity distinct by construction from
    ``C(4,4)/(rho g nabla)`` (T. I. Fossen, personal communication,
    2026-10-07)."""
    main = vessel.main
    return HydrostaticQuantities(
        mass=float(main.m),
        center_of_gravity=np.asarray(main.CG, dtype=float),
        center_of_buoyancy=np.asarray(main.CB, dtype=float),
        displaced_volume=float(main.nabla),
        water_density=float(main.rho),
        gravity=float(main.g),
    )


def select_restoring(vessel: VesselFile, *, body_kind: str) -> Optional[np.ndarray]:
    """The restoring matrix `C` at the file's own reference point, for a
    surface body; `None` for a submerged body, whose equilibrium `C` is
    derived information a submerged-body preprocessor does not read (that
    preprocessor derives its own net-buoyancy residual from
    `hydrostatic_quantities` instead -- Fossen 2011, eq. 4.9, Ch. 4)."""
    if body_kind == "submerged":
        return None
    if body_kind == "surface":
        vel_index = _resolve_velocity_index(vessel, None)
        return _slice_matrix(vessel.C, 0, vel_index)
    raise ValueError(f"unknown body_kind {body_kind!r}, expected 'submerged' or 'surface'")


# ---------------------------------------------------------------------------
# Orchestrate
# ---------------------------------------------------------------------------
def _sha256_of(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_bem_vessel(
    path,
    *,
    added_mass_frequency,
    body_kind: str,
    co_to_reference_point,
    generator_revision: str,
    velocity_index=None,
) -> BemVesselReading:
    """Orchestrate a load, a validation, an added-mass CO shift and the
    hydrostatic quantities into one frozen reading. `rigid_body_mass_
    matrix_export` stays at the file's own reference point (validation
    only; the vehicle's own mass properties are the inertia term it
    actually uses)."""
    path = Path(path)
    vessel = load_vessel_mat(path)
    validate_vessel_file(vessel)

    added_mass_at_reference = select_added_mass(
        vessel, frequency=added_mass_frequency, velocity_index=velocity_index
    )
    added_mass_matrix = shift_added_mass_to_co(
        added_mass_at_reference, co_to_reference_point=co_to_reference_point
    )
    restoring_matrix = select_restoring(vessel, body_kind=body_kind)
    hydrostatics = hydrostatic_quantities(vessel)

    provenance = BemProvenance(
        kind="bem",
        generator=vessel.hydrodynamic_source,
        file_name=path.name,
        sha256=_sha256_of(path),
        generator_revision=generator_revision,
        frequency_limit=added_mass_frequency,
        reference_point=vessel.hydrodynamic_reference,
        axes=vessel.hydrodynamic_axes,
        submergence_depth=float(vessel.main.submergenceDepth),
        co_to_reference_point=np.asarray(co_to_reference_point, dtype=float),
    )

    return BemVesselReading(
        added_mass_matrix=added_mass_matrix,
        restoring_matrix=restoring_matrix,
        hydrostatics=hydrostatics,
        rigid_body_mass_matrix_export=np.asarray(vessel.MRB, dtype=float),
        provenance=provenance,
    )


# ---------------------------------------------------------------------------
# Write
# ---------------------------------------------------------------------------
def _format_number(value: float) -> str:
    return repr(float(value))


def _format_vector(values) -> str:
    return "[" + ", ".join(_format_number(v) for v in np.asarray(values, dtype=float)) + "]"


def _format_matrix(values) -> str:
    rows = [
        "    [" + ", ".join(_format_number(v) for v in row) + "]"
        for row in np.asarray(values, dtype=float)
    ]
    return "[\n" + ",\n".join(rows) + ",\n]"


def write_vessel_part_parameters(
    reading: BemVesselReading,
    *,
    parameters_path,
    source_md_path,
    part_name: str,
) -> None:
    """Write `reading`'s numbers as plain Python literals into a
    composition part's parameters file (one class attribute per value, the
    form the registry's own editor writes and reads) and append one
    provenance line per value to a companion file beside it. Never reopens
    the file the reading came from: once written, the composition is
    self-contained."""
    provenance = reading.provenance
    hydrostatics = reading.hydrostatics
    short_hash = provenance.sha256[:12]

    parameters_lines = [
        f"# {part_name}: vessel-database (BEM) hydrodynamic quantities,",
        f"# read once from {provenance.file_name} ({provenance.generator},",
        f"# revision {provenance.generator_revision}). Source: SOURCE.md beside this file.",
        "",
        "",
        "class ComponentParameters:",
        f"    mass = {_format_number(hydrostatics.mass)}  # bem {short_hash}",
        f"    water_density = {_format_number(hydrostatics.water_density)}  # bem {short_hash}",
        f"    gravity = {_format_number(hydrostatics.gravity)}  # bem {short_hash}",
        f"    displaced_volume = {_format_number(hydrostatics.displaced_volume)}  # bem {short_hash}",
        f"    center_of_gravity = {_format_vector(hydrostatics.center_of_gravity)}  # bem {short_hash}",
        f"    center_of_buoyancy = {_format_vector(hydrostatics.center_of_buoyancy)}  # bem {short_hash}",
        f"    added_mass_matrix = {_format_matrix(reading.added_mass_matrix)}  # bem {short_hash}",
    ]
    parameters_path = Path(parameters_path)
    parameters_path.parent.mkdir(parents=True, exist_ok=True)
    parameters_path.write_text("\n".join(parameters_lines) + "\n")

    source_lines = [
        f"# Vessel-database (BEM) parameters -- {part_name}",
        "",
        f"mass -- bem -- {provenance.file_name} (sha256 {short_hash}...), "
        f"{provenance.generator} {provenance.generator_revision}",
        f"water_density -- bem -- {provenance.file_name} (sha256 {short_hash}...)",
        f"gravity -- bem -- {provenance.file_name} (sha256 {short_hash}...)",
        f"displaced_volume -- bem -- {provenance.file_name} (sha256 {short_hash}...)",
        f"center_of_gravity -- bem -- {provenance.file_name} (sha256 {short_hash}...), "
        f"reference {provenance.reference_point}",
        f"center_of_buoyancy -- bem -- {provenance.file_name} (sha256 {short_hash}...), "
        f"reference {provenance.reference_point}",
        f"added_mass_matrix -- bem -- {provenance.file_name} (sha256 {short_hash}...), "
        f"A({provenance.frequency_limit}) ({provenance.axes}, reference "
        f"{provenance.reference_point}), shifted to CO by H(r), "
        f"r = co_to_reference_point = {provenance.co_to_reference_point.tolist()}, "
        f"zero speed, submergence {provenance.submergence_depth} m",
    ]
    source_md_path = Path(source_md_path)
    source_md_path.parent.mkdir(parents=True, exist_ok=True)
    source_md_path.write_text("\n".join(source_lines) + "\n")
