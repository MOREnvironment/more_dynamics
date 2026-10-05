# Frozen reference — rigid-body + added-mass block (Otter-based catamaran)

**Frozen:** 2026-10-05 (job A-4a, verifier). **Append-only:** a changed byte in either CSV needs an ADR.

## Files

| File | Origin (relative to `more_generic_models/more_generic_models/`) | How it was made | sha256 |
|---|---|---|---|
| `inputs.csv` | `test/plant/asv_catamaran/test_mss_matlab/catamaran_dynamics_inputs.csv` | byte copy | `b653df012854130ebf7af4e916a78fda74abb0ffb0a2b3f9fe82fd593eb1cd48` |
| `matlab_reference.csv` | `test/plant/asv_catamaran/compare_results/compare_dynamics_consistency/catamaran_dynamics_full_debug_matlab.csv` (sha256 `addfe98dc853d8c7e9553eb4e60dc673ebf0d8fd5edcb3d5a41f05d0d44baafc`) | the 180 columns matching `^(M_RB\|M_A\|C_RB\|C_A\|C_total)_\d\d$`, in source order, each field copied as text (no float re-formatting) | `126106d05e78fb04dfcd2ec7844271f566fa016455e3322060298adc319a8ccb` |

Re-extraction (must be byte-equal to `matlab_reference.csv`):

```sh
F=catamaran_dynamics_full_debug_matlab.csv
cols=$(head -1 $F | tr ',' '\n' | grep -n -E '^(M_RB|M_A|C_RB|C_A|C_total)_[0-9]{2}$' | cut -d: -f1 | paste -sd,)
cut -d, -f$cols $F | cmp - matlab_reference.csv
```

## Revisions

- `more_generic_models`: HEAD `e6814a0f80bc140ad9bacea975ae275ddff49bbb`. Both CSVs and the MATLAB script are clean at commit `4d4bffc` (2026-02-18). The three source modules the block ports (`asv_catamaran.py`, `rigid_body_kinetics_surface_vessel.py`, `added_mass.py`) carry one uncommitted line each: the import moved from `more_generic_models.transforms.matrix_transforms` to `more_transformations.matrix_transforms`. The physics is as committed.
- MATLAB generator: `test/plant/asv_catamaran/test_mss_matlab/test_dynamics_consistency.m` (sha256 `66f65393601890cbdef29736cb8fa952d951d5223593f6e66735ea4c3c370b0c`). It carries its own inline copy of the Otter model (`otter_debug`, lines 132–226), **not** a call to MSS `otter.m`.
- MSS anchor for G5: `source-sim/MSS/CRAFT/USV/models/otter.m`, MSS HEAD `99bf0b3`, file last changed `e1dff2a` (2026-06-27), sha256 `66c144fea629cee4bea0ad55b62a1bc5e55b46e68ac8e94860ddc7f0aff12248`.

## Layout and conventions

- 50 rows, one per case. `inputs.csv` columns `x1..x12` = `[nu (6); eta (6)]`, `tau1..tau6` (unused by this block).
- A matrix `X` is stored as `X_01..X_36` **row-major** (`flatten_matrix` transposes before `(:)`, lines 122–129): `X_01..X_06` is the first row.
- `C_RB` is evaluated at `nu`; `C_A` at `nu_r = nu - nu_c` with the current of the generator script: `Vc = 0.3` m/s, `beta_c = deg2rad(30)` (lines 14–15), `nu_c = [Vc cos(beta_c - psi), Vc sin(beta_c - psi), 0, 0, 0, 0]` (lines 147–149), `psi = x12`. `nu_r` is not frozen here; the test recomputes it from these lines.
- `C_total = C_RB + C_A`. `M_total` is not frozen (outside the brief's column list); the test compares `M` with `M_RB + M_A`.
- Payload `mp = 25` kg at `rp = [0.05 0 -0.35]'` m (lines 12–13).

## Known divergence from current MSS

The generator computes the inertia as `Ig = Ig_CG - m*S(rg_total)^2 - mp*S(rp)^2` (line 156), then shifts it again with `H(rg_total)` (lines 158–161). MSS `otter.m` has computed it about the combined CG since its revision of 2026-04-20 (*"Correct payload lever arm moment of inertia"*, otter.m line 78, commit `880b2ef`; formula at lines 123–129). The reference therefore reproduces the **pre-2026-04-20** Otter inertia. The rotational block of `M_RB` differs from current MSS by up to 6.99 kg·m² (roll 21.55 vs 15.53 kg·m²), and `K_pdot`, `M_qdot`, `N_rdot` follow. See ledger `agents-more/30_checks/2026-10-05_A4a_rigid_body_gates.md`.

## Current-MSS regeneration (job A-26, 2026-10-05) — appended, nothing above changed

**Generated:** 2026-10-05 with MATLAB `26.1.0.3312084 (R2026a) Update 4`, MSS `source-sim/MSS` HEAD `99bf0b30e9ae0dca3515d02e0bbd06c54cc0c2f7` (clean). Ledger: `agents-more/30_checks/2026-10-05_A26_matlab_references_current_mss.md`. The `*_mss_current.csv` files are the **current-MSS (default) path**. The legacy files above stay as they are (template path).

How they were made: new generators named `*_mss_current.m` sit beside the old ones in `more_generic_models` and are not committed. They run MSS's own `remus100.m` / `otter.m` on the same `inputs.csv`. A copy of the MSS file is made in `tempdir` at run time; it adds only one extra output, the function's workspace, so the intermediate terms can be read. On every case that copy equals the unmodified MSS function (`xdot`, `M`): max |diff| = 0. External `tau` is added in the model's own state equation, because the MSS functions take actuator commands, which are set to zero.

| File | Origin | How it was made | sha256 |
|---|---|---|---|
| `matlab_reference_mss_current.csv` | `test/plant/asv_catamaran/test_mss_matlab/catamaran_dynamics_full_debug_mss_current.csv` (sha256 `7db813baed8f5a72a3958524fbeb15e4efe46c43c6a82859d28ef9dd95bdcb79`), from `test_dynamics_consistency_mss_current.m` (sha256 `1c344c17bf5a9c5e4b4709d59b02f560f9fef68556859daa7ddc7b2182bd5a1c`) | same 180 columns and re-extraction command as above | `ce33d1401f586aaac6e3c09bae1b4b25157ec357311645ef8804140372229d06` |

This is the **corrected Otter inertia** (E-11 c default path). `otter.m:123–129` gives `M_RB` rotational diagonal 15.52682 / 21.275 / 16.0125 kg·m². The legacy file has 21.55470 / 28.26826 / 16.97788. `M_A` and `C_RB` equal the A-4 test transcription of current `otter.m` (`_otter_m_current_matrices`, `_otter_m_crb`): differences 4.4e-15 and 8.1e-13. The translational blocks are unchanged.

## Readers (job A-30, 2026-10-05) — appended

- `matlab_reference_mss_current.csv` is read by `tests/rigid_body/test_rigid_body_block.py` (G1 block and numpy source, G2 cases, G4); since owner decision E-24 dropped `legacy_otter_inertia`, **`matlab_reference.csv` stays on disk for history and no test reads it**.
