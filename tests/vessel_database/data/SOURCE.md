# Frozen reference — vessel-database (BEM) reader, `LAUV_marie`

**Frozen** files: a changed byte in a CSV needs a decision record; a new MSS revision that
changes an output adds a file `*_mss_<revision>.csv` beside the old one, never writes over it.
The generator below lives in this repository, finds MSS only through `MSS_DIR` (and MATLAB
through `MATLAB_BIN`), and writes under `OUT_DIR` (default `<tempdir>/more_mss_references/vessel_database/`),
never into this folder; a re-run is checked with `cmp` against the file of the same name here.
Nothing here is relative to one machine.

## Files

| File | What it holds | Made by | Made at | sha256 |
|---|---|---|---|---|
| `LAUV_marie.mat` | MSS `vessel` structure (`wamit2vessel.m` 50-94 field list), T. I. Fossen's idealised submerged body for NTNU AUR-Lab's LAUV *Marie* (2.15 m, 0.15 m diameter, 34 kg, 5 m centre depth), produced by the BEM solver [MSS-Capytaine](https://github.com/cybergalactic/MSS-Capytaine) `bfc9dc9` (release 1.0) and shipped **inside MSS itself** | T. I. Fossen, MSS, MIT licence, `HYDRO/vessels_capytaine/LAUV_marie/LAUV_marie.mat` | MSS `cc07579` (2.0.2 + Fossen's post-release fixes, pulled 2026-10-07, A-41) | `7eed3a0c4643ad446b6e45db9f74659052cc56866b3bc3a9c540bd5befcadb61` |
| `vessel_database_lauv_marie_mss_reference.csv` | one frozen row: `A` at the zero- and infinite-frequency sentinels (both at the file's native reference point, the CG) and `MRB`, `C` (6x6 each, row-major `*_01..*_36`); the `main` scalar/vector fields (`m, rho, g, nabla, k44, k55, k66, GM_T, GM_L, Lpp, Lwl, T, B, C_B, submerged, submergenceDepth, CG, CB`); `hydrodynamic_axes`, `hydrodynamic_reference`, `hydrodynamic_source`, `name`; the raw diagonal of `B(omega)` at every one of the 36 frequencies, zero speed, unmodified (Fossen X-5 reply 2: never clip) | `generate_vessel_database_mss.m` (this folder) | MSS `cc07579`, run 2026-10-09 | `7e38579daebcd940dcb754a7fcd4789799b40fe5d24e2cafa15b522721401426` |
| `generate_vessel_database_mss.m` | the generator above; loads the `.mat` with MATLAB's own `load` (no MSS function call — the file is read exactly as committed) | — | — | — |

**Licence (rule 3, E-37 Q2):** MIT (MSS `LICENSE`, Copyright T. I. Fossen). The file is Fossen's
idealised hydrodynamic example, not measured Marie data (his README, quoted below); copied with
credit under the MIT terms.

**Provenance the loader (U10b) must cite beside every value it takes from this file** (A-35 Sec.4
"What the parameter record must say"), e.g.:

> `added_mass_matrix — bem — MSS-Capytaine, shipped in MSS cc07579, HYDRO/vessels_capytaine/LAUV_marie/LAUV_marie.mat (sha256 7eed3a0c…), A(infinity) (freqs index 36, 10 rad/s sentinel), CG -> CO by H(r), r = [0, 0, 0.01] m, zero speed, 5 m submergence, idealised body (Fossen's README), loaded <date>`

## What the fixture says about itself (`LAUV_marie/README.md`, MSS `cc07579`)

> The exact hull offsets, mass distribution, centers, appendages, and propeller geometry are not
> public. The model therefore uses a documented generic tapered-cylinder hull without fins,
> antennas, sonar heads, or propeller. It is a hydrodynamic example rather than a validated
> digital twin of Marie.

Two independent comparisons against a dimension-matched Lamb spheroid (same length and diameter)
confirm this directly: the added-mass diagonal is 25-55% low against this file's BEM numbers, by
amounts a strip-theory check on the hull's actual offsets already explains — this hull's block
coefficient (0.69) is a tapered cylinder's, not a spheroid's (pi/6 = 0.52 for any spheroid), which
accounts for the volume-driven (sway/heave) and volume-squared-times-length-squared-driven
(pitch/yaw) split in the gap. The mismatch is **geometry, not a defect in either formula**.

## Conventions carried into the loader's test contract (rule 19 — do not reinvent)

- **Axes and reference point.** `hydrodynamic_axes = "MSS FSD"` (forward-starboard-down, z down);
  `hydrodynamic_reference = "CG"` — `vessel.A`, `vessel.MRB` and `vessel.C` are **all given at the
  body's centre of gravity**, not at a hull origin (`CO`). `main.CG = [-0, 0, 0.01]` m is the CG's
  position **in CO-relative body-fixed coordinates**: the CG sits 0.01 m below the hull origin (z
  down). This matches `MSS-Capytaine/vessels_capytaine/LAUV_marie/config.json`
  `"center_of_mass_m": [0.0, 0.0, -0.01]` (Capytaine's own, z-up, pre-export convention) and A-66's
  passB finding (`70_reviews/2026-10-09_A66_passB_codex.md`: *"With `r_CG/CO=(0,0,+0.01) m` in
  forward-starboard-down axes"*).
- **Frequency sentinel** (`capytaine_vessel.py:183-186,239-241`; matched by MSS's own `plotBv.m:20-22`):
  `vessel.freqs = [0, <34 finite omegas>, 10]`, 36 entries; the entry at 10 rad/s holds Capytaine's
  `omega = infinity` radiation solution, not a literal result at 10 rad/s. `A(0)` and `A(infinity)`
  read as `vessel.A(:,:,1)` and `vessel.A(:,:,36)` respectively (zero speed, the file's only speed
  column, `vessel.velocities = 0`).
- **CO shift** (A-35 Sec.4 "Outputs"; Fossen 2011 eq. 3.24 p. 49, `more_transformations.matrix_transforms.H_matrix` / its CasADi mirror): `added_mass_matrix_at_CO = H(r)^T @ A_at_CG @ H(r)`, `r` = the
  CO -> CG vector = `main.CG` read from the file (**never a hand-typed 0.01**). Verified against
  A-66 case 2 and its passB independent re-derivation: `H([0,0,0.01])^T @ A_inf_CG @ H([0,0,0.01])`
  reproduces the ledger's CO diagonal (surge 0.650876, sway 33.900835, heave 33.889732, roll
  2.13e-5, pitch 9.76159, yaw 9.764639 kg or kg*m^2) to the digits both ledgers report.
- **Restoring matrix `C`.** Constant over frequency for this submerged body (checked by the
  generator, `c_varies_over_freq = 0` in the CSV: max deviation < 1e-12 across all 36 slices). Per
  ADR 0004 Sec.2 and A-35 Sec.4, a **submerged** body's `C` is **not** handed on to the
  hydrostatics preprocessor (it is derived information that would drop the net-buoyancy residual
  `preprocess_submerged_hydrostatics` needs as `W = m g`, `B = rho g nabla`); only `mass, rho, g,
  nabla, center_of_gravity, center_of_buoyancy` are forwarded. `C44` is still used as a **validation**
  check: `C44 == rho * g * nabla * (z_G - z_B)` at the CG (A-35 Sec.4 item 5) — here
  `rho * g * nabla * (CG_3 - CB_3) = 1025 * 9.81 * 0.033224317683393435 * (0.01 - 2.05406869e-06)
  = 3.3401019830...`, matching the CSV's `C_22` (row 4, col 4, row-major) `3.3401019830143164` to the
  ratio-1.0000000 precision A-35 Sec.4 item 5 already reports; the file's
  own `main.GM_T`, `main.GM_L` are **not** used to reconstruct `C` (X-5 reply 2: they are `C(4,4)/(m g)`,
  not `C(4,4)/(rho g nabla)` — a 0.16% difference Fossen explained directly, not a fault, F-MSS-12
  closed as *explained*).
- **`B(omega)` is out of scope for this first intake** (ADR 0004 Sec.2: "`B(omega)` and the RAOs
  enter later through a fluid-memory block") and is **never clipped**: three frequencies carry a
  small negative diagonal entry (3.2, 3.6, 5.5 rad/s — X-5 item 12, Fossen's reply 2: mesh noise,
  negligible once spectrum-weighted, kept raw on purpose). The frozen CSV's `B_diag_f*_dof*`
  columns exist only so a test can assert the loader never touches, zeroes or clips them.
- **Mass-matrix sanity already checked** (A-35 Sec.4 items 3-4, reproduced by this file's own
  numbers): `MRB` diagonal is `[m, m, m, m k44^2, m k55^2, m k66^2]` with zero off-diagonal
  coupling (consistent with being evaluated at the CG); `A(infinity)` is symmetric to the
  generator's own double precision and its translational block is positive definite.

## Readers

`more_dynamics/tests/vessel_database/test_vessel_database_reader.py` (red until U10b lands the
reader named in that file's module docstring).
