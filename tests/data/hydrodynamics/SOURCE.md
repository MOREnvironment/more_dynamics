# Frozen reference — hydrodynamics blocks (linear damping, cross-flow drag, lift/drag)

**Frozen:** 2026-10-05. **Append-only:** a changed byte in any CSV needs an ADR.

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
- **The stored references are the template (legacy) path, not current MSS.** The spheroid generator differs from current `remus100.m` in exactly the template's departures: rho = 1025 in every function (its own `forceLiftDragg`, line 467; `imlay611`, line 407), the sway damping fades with speed (line 268, `D(2,2)*exp(-3*U_r)`), and cross-flow drag uses the pre-2026-08-26 `crossFlowDrag.m` (21 end points; inferred from the size of the residual against current MSS, 2026-10-05). Reynolds number on the length, as current MSS.
- The catamaran generator uses `g = gravity(mu)` (line 137), where `otter.m` uses `g = 9.81` (line 90): it enters `Xu = -24.4 g / Umax` (line 187) and the restoring matrix `G`.
- `G` is at the CO; the test moves it to the centre of flotation with `H(-r_f)`, `r_f = [LCF 0 0]`, `LCF = -0.2` (generator line 182; `otter.m` line 192). `D` and `tau_damp` follow `otter.m`'s sign: diagonal entries are the negative derivatives `Xu..Nr`, and `tau_damp` is **added** to the right-hand side.
- Spheroid `D` follows `Dmtrx.m`'s sign (positive diagonal) and is **subtracted** (`- D*nu_r`).

## Current-MSS regeneration (2026-10-05) — appended, nothing above changed

**Generated:** 2026-10-05 with MATLAB `26.1.0.3312084 (R2026a) Update 4`, MSS `source-sim/MSS` HEAD `99bf0b30e9ae0dca3515d02e0bbd06c54cc0c2f7` (clean). The `*_mss_current.csv` files are the **current-MSS (default) path**. The legacy files above stay as they are (template path).

How they were made: new generators named `*_mss_current.m` sit beside the old ones in `more_generic_models` (committed there in `524e336`, 2026-10-06). They run MSS's own `remus100.m` / `otter.m` on the same `inputs.csv`. A copy of the MSS file is made in `tempdir` at run time; it adds only one extra output, the function's workspace, so the intermediate terms can be read. On every case that copy equals the unmodified MSS function (`xdot`, `M`): max |diff| = 0. External `tau` is added in the model's own state equation, because the MSS functions take actuator commands, which are set to zero.

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

## Readers (2026-10-05) — appended

- `tests/hydrodynamics/test_hydrodynamics_block.py`: the `*_mss_current.csv` files for the default paths (section "G1-MSS"; spheroid cross-flow on the then MSS flag `cross_flow_reynolds_length="length"`, the Reynolds number on the length, removed 2026-10-06), the legacy `*_matlab_reference.csv` files for the template flags, the numpy source's G1, G2 and G4.

## Surge, ITTC and floating-damping references; cross-flow at `ac77394` (2026-10-06) — appended, nothing above changed

**MSS** `source-sim/MSS` HEAD `ac77394b74a2184317e92c6b73d8deca97b48359` (clean), which carries issue #81's fix: `HYDRO/cylinderDrag.m` 79-80, `Re = U_crossflow * B / nu_water`. Between `99bf0b3` and `ac77394`, the only other file these tests use that changed is `LIBRARY/modeling/Dmtrx.m`, in a comment only (line 4, `98506da`). New: `LIBRARY/modeling/XuuITTC.m` (`97fae93`). **MATLAB** `26.1.0.3312084 (R2026a) Update 4`.

| File | Origin | How it was made | sha256 |
|---|---|---|---|
| `spheroid_matlab_reference_mss_ac77394.csv` | the 2026-10-05 generator `test/plant/auv_spheroid/test_mss_matlab/test_dynamics_consistency_mss_current.m` (sha256 `51ca4117…c33e`, unchanged) run from a scratch copy with its input CSV, `MSS_DIR` = the checkout above; full output sha256 `1fba1e7f…896a` | same 128 columns and `cut_cols` regex as above | `5303cdc0a205d2d8fb06bd861718a8f0243a33def8c99e78951436516de63041` |
| `surge_damping_mss_ac77394.csv` | `generate_u3_surge_floating_mss.m` (this folder) | MSS `forceSurgeDamping.m` (both branches) and `addedMassSurge.m` called directly; set 1 = `osv.m` 62-80, 129; sets 2-5 seeded (`rng(20261006,'twister')`); 40 speeds per set and branch + 4 probe rows at Rn = 100 (`probe = 1`, physical tests only) | `752fbc4d90b38e32db0943b7603ffcdfa2f81ccfee26801aabd86311543ebe80` |
| `xuu_ittc_mss_ac77394.csv` | same generator | MSS `XuuITTC.m` called directly, same five sets, 40 speeds | `3ccbe4b4cbd5b0349c6abe05adc5ff04dbd3f90f6294ac30bd2aecba20539721` |
| `floating_damping_mss_ac77394.csv` | same generator | MSS `Dmtrx.m` surface branch (G a matrix) on 20 seeded craft; matrices row-major | `6cf64656315ebf354996980445670b0e453474ac85804358796ca437ccee8bc4` |
| `generate_u3_surge_floating_mss.m` | written 2026-10-06 | needs `MSS_DIR` (errors without it); run with `MSS_DIR=<MSS> "$MATLAB_BIN" -batch "run('generate_u3_surge_floating_mss.m')"` from this folder; writes `%.17g` (exact doubles) | `1268d96827d666cbab9b07ba3ba2055c80f987f4feb4a5b01388725572f60e3a` |
| `cited_lines_snapshot.json` | `snapshot_cited_lines.py` (this folder) with `MSS_DIR` = the checkout above and `MORE_GENERIC_MODELS_DIR` = `more_generic_models` at `524e336` | the stripped text of all 134 cited lines and the 4 MSS tables the tests read; the gates read this file, so they run without either checkout (nothing relative to a machine) | `785a29b6116644ff6d116edb181e4576df933fa2fcce73af66bea4b9230d3b21` |
| `snapshot_cited_lines.py` | written 2026-10-06 | needs both variables; refuses a cited line that no longer starts with its pinned text | `f580674b414db24d4c30063edaf99699f95f27e518a52110355967e259f3aee6` |

Checks on the new spheroid file: the header equals that of `spheroid_matlab_reference_mss_current.csv`, and only `tau_crossflow_02/03/05/06` differ, by up to 175.5 N (test `test_ac77394_reference_differs_from_99bf0b3_only_in_cross_flow`). Byte-unchanged, re-checked: `inputs.csv`, both legacy files, both 2026-10-05 current-MSS files (sha256 as in the tables above).

**Readers (2026-10-06):** `tests/hydrodynamics/test_hydrodynamics_block.py`. The ac77394 file is the default-path G1 of the cross-flow block. `spheroid_matlab_reference_mss_current.csv` (99bf0b3) remains the G1 file for damping and lift/drag, and the history file for cross-flow (pre-fix Re on the length, transcription only). The block test on the dropped `cross_flow_reynolds_length` flag was removed (owner, 2026-10-06).

**Paths:** nothing here or in the test is relative to a workspace. MSS and `more_generic_models` are reached only through `MSS_DIR` and `MORE_GENERIC_MODELS_DIR`, and MATLAB through `MATLAB_BIN`. The 2026-10-05 generators in `more_generic_models` then fell back to a sibling-folder path when `MSS_DIR` was unset (it was set explicitly for these runs); since 2026-10-07 they require `MSS_DIR`.

## MSS release 2.0.2 (`72656d1`), re-checked 2026-10-07 — appended, nothing above changed

MSS moved to `72656d1b3b962bf533003ded2049347371153062` (release 2.0.2, T. I. Fossen; commits `108ceda` 2026-10-06 to `72656d1` 2026-10-07). Every MATLAB generator named above was re-run on 2026-10-07 with MATLAB `26.1.0.3312084 (R2026a) Update 4` and `MSS_DIR` at that revision, from a scratch copy with its input files, and its output byte-compared with the frozen file (extractions with the commands above). This table lists every file now in the folder; it is the one the folder is checked against.

| File | Made at | Status at MSS `72656d1` (2026-10-07) | sha256 now |
|---|---|---|---|
| `inputs.csv` | `more_generic_models` input file (byte copy) | unchanged input | `b653df012854130ebf7af4e916a78fda74abb0ffb0a2b3f9fe82fd593eb1cd48` |
| `spheroid_matlab_reference.csv` | template generator (inline `remus100.m` copy), not MSS | not an MSS product; not re-anchored | `2ba042b128b366fc059916d964d929db35e43d9c93303c0c853ab3fc869a390a` |
| `catamaran_matlab_reference.csv` | template generator (inline `otter.m` copy), not MSS | not an MSS product; not re-anchored | `7c6e805825cced99449fdb2b8f6612b2395f26fb90e875503caa4bbc93b8b8fb` |
| `spheroid_matlab_reference_mss_current.csv` | MSS `99bf0b3` | **differs** — the `72656d1` run changes only `tau_crossflow_02/03/05/06`, by up to 175.5 / 167.8 / 58.41 / 61.98 (N, N, N·m, N·m); every other column byte-equal. Cause: the cross-flow Reynolds number on the diameter (`HYDRO/cylinderDrag.m:80`, fixed in `ac77394`, 2026-10-05), not 2.0.2. Stays anchored at `99bf0b3` as the pre-fix history file; no `*_mss_2_0_2.csv` was added because the `72656d1` output is byte-identical to `spheroid_matlab_reference_mss_ac77394.csv` | `63abfa7e47411827945a595bebdf3bedd9b684b81c5c5a15861501ebc0d0e27d` |
| `spheroid_matlab_reference_mss_ac77394.csv` | MSS `ac77394` | **byte-identical** → also MSS `72656d1` | `5303cdc0a205d2d8fb06bd861718a8f0243a33def8c99e78951436516de63041` |
| `catamaran_matlab_reference_mss_current.csv` | MSS `99bf0b3` | **byte-identical** → also MSS `72656d1` | `6bb2ccab765c11cff535d0a7be0579cac94a97dde4fc6da15ff82d9faf3d9c1b` |
| `surge_damping_mss_ac77394.csv` | MSS `ac77394` | **stays anchored at `ac77394`**: `LIBRARY/modeling/forceSurgeDamping.m` was deleted in MSS commit `108ceda` (2026-10-06, "MSS Toolbox 2.0 release"); not regenerated. (Control run with that file taken from `ac77394` beside MSS `72656d1`: byte-identical, `addedMassSurge.m` changed only a comment line.) | `752fbc4d90b38e32db0943b7603ffcdfa2f81ccfee26801aabd86311543ebe80` |
| `xuu_ittc_mss_ac77394.csv` | MSS `ac77394` | **byte-identical** (`XuuITTC.m` unchanged in 2.0.2; run with `forceSurgeDamping.m` from `ac77394` beside MSS `72656d1` so the seeded draws reach this section unchanged) → also MSS `72656d1` | `3ccbe4b4cbd5b0349c6abe05adc5ff04dbd3f90f6294ac30bd2aecba20539721` |
| `floating_damping_mss_ac77394.csv` | MSS `ac77394` | **byte-identical** (`Dmtrx.m` unchanged in 2.0.2; same run) → also MSS `72656d1` | `6cf64656315ebf354996980445670b0e453474ac85804358796ca437ccee8bc4` |
| `osv_surge_damping_mss_2_0_2.csv` | MSS `72656d1`, 2026-10-07 (new, below) | new | `77aad05a4c4498ef0b0dcc62e1a16beffa16337c2132dfebd65458d24b66dd0f` |
| `generate_osv_surge_mss_2_0_2.m` | generator (new, below) | needs MSS 2.0.2 or later (stops if `forceSurgeDamping` is on the path); re-run in place byte-identical | `8a81bcab934319b7f948009c4d9f01c42cc6a71bcc149b12675b9c02ee95a490` |
| `generate_u3_surge_floating_mss.m` | generator | header comments edited 2026-10-07 (revision note, no code change); runs only where `forceSurgeDamping.m` exists (MSS `ac77394` or older) | `e9e0610298e2559740ca73b0b9b39a822d019b23c23b201847e5c74afe0c95b5` |
| `snapshot_cited_lines.py` | script | comments, exit message and the `made_by` text edited 2026-10-07 (no logic change) | `a3b8761e01f68b3bfc9cb1b994d3bc36c7444410b0a99498fcaba41552b5a260` |
| `cited_lines_snapshot.json` | re-pinned at MSS `72656d1` on 2026-10-07 by the cited-line re-pin (its `revisions` field) | not regenerated or checked here; sha256 as read at the end of this check | `3e41f30877d5e33bdc15d09218c36c0eaadc8dacb0d380ec89693418864bbf6c` |

Generators in `more_generic_models`, headers edited on 2026-10-07 (revision note; `MSS_DIR` now required instead of a sibling-folder fallback; outputs re-run byte-identical): spheroid `test/plant/auv_spheroid/test_mss_matlab/test_dynamics_consistency_mss_current.m` sha256 `0df62b81cde0b10f79ef67936a45513cef79595edb8633eed1a58684015f2f0a`; catamaran `test/plant/asv_catamaran/test_mss_matlab/test_dynamics_consistency_mss_current.m` sha256 `b3074f217a7d317e6fc9e41764bec2ccb318954d3be906a14c4823236498fb72`.

What 2.0.2 changed in the MSS files these data come from: `remus100.m`, `otter.m`, `crossFlowDrag.m`, `addedMassSurge.m`, `imlay61.m`, `m2c.m` — comments, layout and the book reference only (Fossen 2021 → *Fossen (2027), Handbook of Marine Craft Hydrodynamics and Motion Control, 3rd ed.*, not on disk); `Dmtrx.m`, `XuuITTC.m`, `forceLiftDrag.m`, `coeffLiftDrag.m`, `HYDRO/cylinderDrag.m`, `HYDRO/Hoerner.m`, `gravity.m` unchanged. **Line numbers moved:** `otter.m` every line −1; `remus100.m` −1 before line 203, +2 from line 205; `crossFlowDrag.m` code +12 to +14 (new References block in its header: `rho`/`nStrips` 24–25 → 36–37, `Hoerner` 31 → 44, `cylinderDrag` 34 → 48); `osv.m` −1 before line 176. The line numbers written above stay those of the revision named beside them. **`forceSurgeDamping.m` is gone** from MSS; `CRAFT/SHIP/models/osv.m`, its only MSS caller, was rewritten (next section).

### `osv_surge_damping_mss_2_0_2.csv` (220 × 20) — surge damping of `osv.m` at `72656d1`

MSS 2.0.2 damps the OSV in surge with (`CRAFT/SHIP/models/osv.m` 178–185 @ `72656d1`)

`D_nl_11 = exp(−k_u·|u_r|)·D11 − Xuu·|u_r|`, `k_u = 3`, surge force `X = −D_nl_11·u_r` (surge row of `−(CRB + CA + D_nonlinear)·nu_r`, `osv.m` 203–207; `D` from `Dmtrx.m` is diagonal),

with `D11 = M(1,1)/T1` (`LIBRARY/modeling/Dmtrx.m` 62) and `Xuu` from `LIBRARY/modeling/XuuITTC.m` (ITTC-1957 line, Reynolds number bounded below at 1e5, Mumford wetted area). The generator `generate_osv_surge_mss_2_0_2.m` (this folder) evaluates it with the MSS functions on the speed grid of `surge_damping_mss_ac77394.csv` (sets 1–5, branch 0: 44 speeds per set including the 4 probe rows), with `B, T, C_B` of each set from `xuu_ittc_mss_ac77394.csv`. Set 1 is `osv.m` itself: `D11 = 59768.875` (`Dmtrx` on `osv.m`'s own `MRB, MA, G`, the call of lines 133–134), `M11 = 5976887.5`, and the formula equals `D_nonlinear(1,1)` read from `osv.m` (instrumented copy, equal to the unmodified function to 0) to 0. Sets 2–5 (test construction): `D11 = Dmtrx(...)(1,1)` with `M(1,1) = m + A11`, `A11` the `addedMassSurge.m` column of the grid file.

Columns: `set, probe, m, L, B, T, C_B, rho, T1, M11, D11, k_u, u_r, Xuu, D_nl_11, X`, then for set 1 only (NaN otherwise) `D_nl_11_osv_first, D_nl_11_osv_second, nu_dot_1_osv_first, nu_dot_1_osv_second`: the unmodified `osv.m` called twice on the same state `x = [u_r, 0, …, 0]` after `clear osv`. `osv.m` keeps its parameters in a `persistent` struct and sets `vessel.D(1,1) = 0` (line 188) after line 185 has read it, so **the linear term exists only on the first call**: from the second call on, `D_nl_11 = −Xuu·|u_r|`. At `u_r = 0.1` m/s the coefficient is 44578.2 on the first call and 300.3 on the second (kg/s), at 1 m/s 5049.6 and 2073.9, at 8 m/s 12489.536 on both (they differ by 2.3e-6). The first-call value is the formula the file's comment states (line 182); the `*_second` columns record what a simulation sees. Units: m, kg, kg/m³, s, m/s, kg/m (`Xuu`), kg/s (`D11`, `D_nl_11`), N (`X`), m/s² (`nu_dot_1`). `%.17g`.

Run: `MSS_DIR=<MSS checkout at 72656d1 or later> "$MATLAB_BIN" -batch "run('generate_osv_surge_mss_2_0_2.m')"` from this folder (no machine path). Readers: none yet — a candidate third surge form beside the two of `surge_damping.py`, for the owner to decide.
