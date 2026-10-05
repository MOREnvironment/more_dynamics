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
