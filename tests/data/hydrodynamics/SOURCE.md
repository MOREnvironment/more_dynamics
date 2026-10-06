# Frozen reference — hydrodynamics blocks (linear damping, cross-flow drag, lift/drag)

**Frozen:** 2026-10-05 (job A-22, verifier). **Append-only:** a changed byte in any CSV needs an ADR.

## Files

| File | Origin (relative to `more_generic_models/more_generic_models/`) | How it was made | sha256 |
|---|---|---|---|
| `inputs.csv` | `test/plant/auv_spheroid/test_mss_matlab/spheroid_auv_dynamics_inputs.csv` | byte copy; byte-identical to `test/plant/asv_catamaran/test_mss_matlab/catamaran_dynamics_inputs.csv`, so it serves both vehicles | `b653df012854130ebf7af4e916a78fda74abb0ffb0a2b3f9fe82fd593eb1cd48` |
| `spheroid_matlab_reference.csv` | `test/plant/auv_spheroid/compare_results/compare_dynamics_consistency/spheroid_auv_dynamics_full_debug_matlab.csv` (sha256 `326cfb58…dab423b1`) | the 128 columns matching `^(nu_r_0[1-6]\|alpha\|U_r\|tau_lift_drag_0[1-6]\|tau_crossflow_0[1-6]\|M_RB_[0-9]{2}\|M_A_[0-9]{2}\|D_[0-9]{2})$`, source order, fields copied as text | `2ba042b128b366fc059916d964d929db35e43d9c93303c0c853ab3fc869a390a` |
| `catamaran_matlab_reference.csv` | `test/plant/asv_catamaran/compare_results/compare_dynamics_consistency/catamaran_dynamics_full_debug_matlab.csv` (sha256 `addfe98d…d0d44baafc`) | the 126 columns matching `^(nu_r_0[1-6]\|tau_damp_0[1-6]\|tau_crossflow_0[1-6]\|M_total_[0-9]{2}\|G_[0-9]{2}\|D_[0-9]{2})$`, source order, fields copied as text | `7c6e805825cced99449fdb2b8f6612b2395f26fb90e875503caa4bbc93b8b8fb` |

Re-extraction (each must be byte-equal):

```sh
cut_cols(){ cols=$(head -1 $1 | tr ',' '\n' | grep -n -E "$2" | cut -d: -f1 | paste -sd,); cut -d, -f$cols $1; }
cut_cols spheroid_auv_dynamics_full_debug_matlab.csv '^(nu_r_0[1-6]|alpha|U_r|tau_lift_drag_0[1-6]|tau_crossflow_0[1-6]|M_RB_[0-9]{2}|M_A_[0-9]{2}|D_[0-9]{2})$' | cmp - spheroid_matlab_reference.csv
cut_cols catamaran_dynamics_full_debug_matlab.csv '^(nu_r_0[1-6]|tau_damp_0[1-6]|tau_crossflow_0[1-6]|M_total_[0-9]{2}|G_[0-9]{2}|D_[0-9]{2})$' | cmp - catamaran_matlab_reference.csv
```

## Revisions

- `more_generic_models`: HEAD `e6814a0f80bc140ad9bacea975ae275ddff49bbb`; origin CSVs and both generators clean at `4d4bffc` (2026-02-18).
- Spheroid generator: `test/plant/auv_spheroid/test_mss_matlab/test_dynamics_consitency.m` (sha256 `821ca50d…ebbcd30`), local `remus100_core` (lines 160–374).
- Catamaran generator: `test/plant/asv_catamaran/test_mss_matlab/test_dynamics_consistency.m` (sha256 `66f65393…3c370b0c`), local `otter_debug` (lines 132–226).
- MSS (`source-sim/MSS`, HEAD `99bf0b3`): `remus100.m` `1264bbd` (2026-08-26), `otter.m` `e1dff2a`, `Dmtrx.m` `ac792a9`, `forceLiftDrag.m` `ac792a9`, `coeffLiftDrag.m` `cb09f7c`, `crossFlowDrag.m` `6a2a064` (2026-08-26, 20 strip midpoints), `HYDRO/cylinderDrag.m` `8acb535`, `HYDRO/Hoerner.m` `2be7156`, `INS/functions/gravity.m` `6d4b553`.

## Layout and conventions

- 50 rows; a matrix `X` is `X_01..X_36` **row-major** (generators' `flatten_matrix`: transpose then `(:)`).
- `nu_r = nu - [Vc cos(beta_c - psi), Vc sin(beta_c - psi), 0, 0, 0, 0]` with `Vc = 0.3`, `beta_c = 30 deg` (spheroid generator lines 21–22, 219–225; catamaran lines 14–15, 147–149). Frozen, not recomputed.
- **The stored references are the template (legacy) path, not current MSS.** The spheroid generator differs from current `remus100.m` in exactly the template's departures: rho = 1025 in every function (its own `forceLiftDragg`, line 467; `imlay611`, line 407), the sway damping fades with speed (line 268, `D(2,2)*exp(-3*U_r)`), and cross-flow drag uses the pre-2026-08-26 `crossFlowDrag.m` (21 end points; inferred in ledger `agents-more/30_checks/2026-10-05_A2_crossflow_residual.md`). Reynolds number on the length, as current MSS.
- The catamaran generator uses `g = gravity(mu)` (line 137), where `otter.m` uses `g = 9.81` (line 90): it enters `Xu = -24.4 g / Umax` (line 187) and the restoring matrix `G`.
- `G` is at the CO; the test moves it to the centre of flotation with `H(-r_f)`, `r_f = [LCF 0 0]`, `LCF = -0.2` (generator line 182; `otter.m` line 192). `D` and `tau_damp` follow `otter.m`'s sign: diagonal entries are the negative derivatives `Xu..Nr`, and `tau_damp` is **added** to the right-hand side.
- Spheroid `D` follows `Dmtrx.m`'s sign (positive diagonal) and is **subtracted** (`- D*nu_r`).

## Current-MSS regeneration (job A-26, 2026-10-05) — appended, nothing above changed

**Generated:** 2026-10-05 with MATLAB `26.1.0.3312084 (R2026a) Update 4`, MSS `source-sim/MSS` HEAD `99bf0b30e9ae0dca3515d02e0bbd06c54cc0c2f7` (clean). Ledger: `agents-more/30_checks/2026-10-05_A26_matlab_references_current_mss.md`. The `*_mss_current.csv` files are the **current-MSS (default) path**. The legacy files above stay as they are (template path).

How they were made: new generators named `*_mss_current.m` sit beside the old ones in `more_generic_models` and are not committed. They run MSS's own `remus100.m` / `otter.m` on the same `inputs.csv`. A copy of the MSS file is made in `tempdir` at run time; it adds only one extra output, the function's workspace, so the intermediate terms can be read. On every case that copy equals the unmodified MSS function (`xdot`, `M`): max |diff| = 0. External `tau` is added in the model's own state equation, because the MSS functions take actuator commands, which are set to zero.

| File | Origin | How it was made | sha256 |
|---|---|---|---|
| `spheroid_matlab_reference_mss_current.csv` | `test/plant/auv_spheroid/test_mss_matlab/spheroid_auv_dynamics_full_debug_mss_current.csv` (sha256 `f490c79bc7b23fbd3aa13d3de02480d13d3cff762375e91026064e8a6c54b5c9`), from `test_dynamics_consistency_mss_current.m` (sha256 `51ca41176581c5d390ac5fc12c66eb82750ea5406193273d3371cb6c17a33c3e`) | same 128 columns and `cut_cols` regex as above | `63abfa7e47411827945a595bebdf3bedd9b684b81c5c5a15861501ebc0d0e27d` |
| `catamaran_matlab_reference_mss_current.csv` | `test/plant/asv_catamaran/test_mss_matlab/catamaran_dynamics_full_debug_mss_current.csv` (sha256 `7db813baed8f5a72a3958524fbeb15e4efe46c43c6a82859d28ef9dd95bdcb79`), from `test_dynamics_consistency_mss_current.m` (sha256 `1c344c17bf5a9c5e4b4709d59b02f560f9fef68556859daa7ddc7b2182bd5a1c`) | same 126 columns and `cut_cols` regex as above | `6bb2ccab765c11cff535d0a7be0579cac94a97dde4fc6da15ff82d9faf3d9c1b` |

Inputs: the same `inputs.csv` (the generators read their own byte-identical copy). What differs from the legacy files, all by current MSS:

- **Spheroid** (`remus100.m` `1264bbd`):
  - `M_A` uses `imlay61.m`'s rho = 1026, so it is ×1026/1025.
  - `tau_lift_drag` uses `forceLiftDrag.m`'s rho = 1026, so it is ×1026/1025.
  - `D` fades **surge only** (`remus100.m:216`) and is built on the new `M`.
  - `tau_crossflow` uses `crossFlowDrag.m` `6a2a064`: 20 strip midpoints, Reynolds number on the length, full `cylinderDrag.m` table.
- **Catamaran** (`otter.m` `e1dff2a`):
  - `g = 9.81` (`otter.m:90`), so `G` is ×0.99880 and `Xu = -24.4·9.81/Umax`.
  - The inertia is about the combined CG (`otter.m:123–129`). This changes `M_total` and the rotational added mass, and through them `D` (Kp, Mq, Nr).
  - `tau_crossflow` uses 20 strip midpoints (Hoerner).
  - `D = diag([Xu Yv Zw Kp Mq Nr])`, the same sign convention as above (negative derivatives).

## Readers (job A-30, 2026-10-05) — appended

- `tests/hydrodynamics/test_hydrodynamics_block.py`: the `*_mss_current.csv` files for the default paths (section "G1-MSS"; spheroid cross-flow on the MSS flag `cross_flow_reynolds_length="length"`, register row D-MSS-1), the legacy `*_matlab_reference.csv` files for the template flags, the numpy source's G1, G2 and G4.

## Job U3a (verifier, 2026-10-06) — appended, nothing above changed

**MSS** `source-sim/MSS` HEAD `ac77394b74a2184317e92c6b73d8deca97b48359` (clean), which carries issue #81's fix: `HYDRO/cylinderDrag.m` 79-80, `Re = U_crossflow * B / nu_water`. Between `99bf0b3` and `ac77394`, the only other file these tests use that changed is `LIBRARY/modeling/Dmtrx.m`, in a comment only (line 4, `98506da`). New: `LIBRARY/modeling/XuuITTC.m` (`97fae93`). **MATLAB** `26.1.0.3312084 (R2026a) Update 4`. Ledger: `agents-more/30_checks/2026-10-06_U3_hydrodynamics_gates.md`.

| File | Origin | How it was made | sha256 |
|---|---|---|---|
| `spheroid_matlab_reference_mss_ac77394.csv` | the A-26 generator `test/plant/auv_spheroid/test_mss_matlab/test_dynamics_consistency_mss_current.m` (sha256 `51ca4117…c33e`, unchanged) run from a scratch copy with its input CSV, `MSS_DIR` = the checkout above; full output sha256 `1fba1e7f…896a` | same 128 columns and `cut_cols` regex as above | `5303cdc0a205d2d8fb06bd861718a8f0243a33def8c99e78951436516de63041` |
| `surge_damping_mss_ac77394.csv` | `generate_u3_surge_floating_mss.m` (this folder) | MSS `forceSurgeDamping.m` (both branches) and `addedMassSurge.m` called directly; set 1 = `osv.m` 62-80, 129; sets 2-5 seeded (`rng(20261006,'twister')`); 40 speeds per set and branch + 4 probe rows at Rn = 100 (`probe = 1`, physical tests only) | `752fbc4d90b38e32db0943b7603ffcdfa2f81ccfee26801aabd86311543ebe80` |
| `xuu_ittc_mss_ac77394.csv` | same generator | MSS `XuuITTC.m` called directly, same five sets, 40 speeds | `3ccbe4b4cbd5b0349c6abe05adc5ff04dbd3f90f6294ac30bd2aecba20539721` |
| `floating_damping_mss_ac77394.csv` | same generator | MSS `Dmtrx.m` surface branch (G a matrix) on 20 seeded craft; matrices row-major | `6cf64656315ebf354996980445670b0e453474ac85804358796ca437ccee8bc4` |
| `generate_u3_surge_floating_mss.m` | written by U3a | needs `MSS_DIR` (errors without it); run with `MSS_DIR=<MSS> "$MATLAB_BIN" -batch "run('generate_u3_surge_floating_mss.m')"` from this folder; writes `%.17g` (exact doubles) | `1268d96827d666cbab9b07ba3ba2055c80f987f4feb4a5b01388725572f60e3a` |
| `cited_lines_snapshot.json` | `snapshot_cited_lines.py` (this folder) with `MSS_DIR` = the checkout above and `MORE_GENERIC_MODELS_DIR` = `more_generic_models` at `524e336` | the stripped text of all 134 cited lines and the 4 MSS tables the tests read; the gates read this file, so they run without either checkout (agents-more rule 9) | `785a29b6116644ff6d116edb181e4576df933fa2fcce73af66bea4b9230d3b21` |
| `snapshot_cited_lines.py` | written by U3a | needs both variables; refuses a cited line that no longer starts with its pinned text | `f580674b414db24d4c30063edaf99699f95f27e518a52110355967e259f3aee6` |

Checks on the new spheroid file: the header equals that of `spheroid_matlab_reference_mss_current.csv`, and only `tau_crossflow_02/03/05/06` differ, by up to 175.5 N (test `test_ac77394_reference_differs_from_99bf0b3_only_in_cross_flow`). Byte-unchanged, re-checked: `inputs.csv`, both legacy files, both A-26 files (sha256 as in the tables above).

**Readers (U3a):** `tests/hydrodynamics/test_hydrodynamics_block.py`. The ac77394 file is the default-path G1 of the cross-flow block. `spheroid_matlab_reference_mss_current.csv` (99bf0b3) remains the G1 file for damping and lift/drag, and the history file for cross-flow (pre-fix Re on the length, transcription only). The block test on the dropped `cross_flow_reynolds_length` flag was removed (owner, 2026-10-06).

**Paths (agents-more rule 9):** nothing here or in the test is relative to a workspace. MSS and `more_generic_models` are reached only through `MSS_DIR` and `MORE_GENERIC_MODELS_DIR`, and MATLAB through `MATLAB_BIN`. The A-26 generators in `more_generic_models` still default to a sibling-folder path when `MSS_DIR` is unset; U3a set it explicitly (their owner's edit, listed in the ledger).
