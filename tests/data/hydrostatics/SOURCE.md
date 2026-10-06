# Frozen reference — submerged hydrostatics block (restoring forces g(eta))

**Frozen:** 2026-10-05 (job A-22, verifier). **Append-only:** a changed byte in either CSV needs an ADR.

## Files

| File | Origin (relative to `more_generic_models/more_generic_models/`) | How it was made | sha256 |
|---|---|---|---|
| `inputs.csv` | `test/plant/auv_spheroid/test_mss_matlab/spheroid_auv_dynamics_inputs.csv` | byte copy (byte-identical to the catamaran inputs used in `../rigid_body/`) | `b653df012854130ebf7af4e916a78fda74abb0ffb0a2b3f9fe82fd593eb1cd48` |
| `spheroid_matlab_reference.csv` | `test/plant/auv_spheroid/compare_results/compare_dynamics_consistency/spheroid_auv_dynamics_full_debug_matlab.csv` (sha256 `326cfb58f97dd87577c738ba74fa1375be1d580911eeb40ff99ae0d2dab423b1`) | the 7 columns matching `^(M_RB_01\|g_0[1-6])$`, in source order, fields copied as text | `22a64de9587f77ddecc4fb74cc296d70f7f6f3c376fdfb11365757c19de21a7c` |

Re-extraction (must be byte-equal):

```sh
F=spheroid_auv_dynamics_full_debug_matlab.csv
cols=$(head -1 $F | tr ',' '\n' | grep -n -E '^(M_RB_01|g_0[1-6])$' | cut -d: -f1 | paste -sd,)
cut -d, -f$cols $F | cmp - spheroid_matlab_reference.csv
```

## Revisions

- `more_generic_models`: HEAD `e6814a0f80bc140ad9bacea975ae275ddff49bbb`. Both origin CSVs and the generator are clean at `4d4bffc` (2026-02-18).
- Generator: `test/plant/auv_spheroid/test_mss_matlab/test_dynamics_consitency.m` (sha256 `821ca50d389a48ab1b64b44fd9dc6437bd98b59cddd6e4ec649925ee6ebbcd30`), local function `remus100_core` (lines 160–374), an inline copy of MSS `remus100.m` without actuators. `g` is `gRvect(W,B,R,r_bG,r_bB)` (line 282) with `W = B = MRB(1,1) * gravity(mu)` (lines 204–205, 259–261), `R` from `eulerang(phi,theta,psi)` (line 280), `r_bG = [0 0 0.02]'`, `r_bB = [0 0 0]'` (lines 238–239).
- MSS anchor (`source-sim/MSS`, HEAD `99bf0b3`): `LIBRARY/modeling/gRvect.m` (`d7c1e96`, 2026-08-21), `LIBRARY/modeling/gvect.m` (`d7c1e96`), `LIBRARY/kinematics/Rzyx.m` (`ac792a9`), `INS/functions/gravity.m` (`6d4b553`), `CRAFT/AUV/models/remus100.m` (`1264bbd`, 2026-08-26).

## Layout and conventions

- 50 rows. `inputs.csv` columns `x1..x12` = `[nu (6); eta (6)]`, `tau1..tau6` (unused here). `eta = [x y z phi theta psi]`, NED, z down.
- `M_RB_01` is the rigid-body mass `m` of `spheroid(a,b,...)` with rho = 1025; the test forms `W = m * gravity(mu)` from it (`gravity.m` lines 11–12, `mu` = 63.446827 deg, generator line 204). It is frozen here only so the weight can be rebuilt exactly as MATLAB built it.
- The vehicle is neutrally buoyant (`B = W`), so the reference exercises only the moment rows; the force rows are zero. The block's force rows are tested against `gRvect.m` with `W != B` inside the test.

## Current-MSS regeneration (job A-26, 2026-10-05) — appended, nothing above changed

**Generated:** 2026-10-05 with MATLAB `26.1.0.3312084 (R2026a) Update 4`, MSS `source-sim/MSS` HEAD `99bf0b30e9ae0dca3515d02e0bbd06c54cc0c2f7` (clean). Ledger: `agents-more/30_checks/2026-10-05_A26_matlab_references_current_mss.md`. The `*_mss_current.csv` files are the **current-MSS (default) path**. The legacy files above stay as they are (template path).

How they were made: new generators named `*_mss_current.m` sit beside the old ones in `more_generic_models` and are not committed. They run MSS's own `remus100.m` / `otter.m` on the same `inputs.csv`. A copy of the MSS file is made in `tempdir` at run time; it adds only one extra output, the function's workspace, so the intermediate terms can be read. On every case that copy equals the unmodified MSS function (`xdot`, `M`): max |diff| = 0. External `tau` is added in the model's own state equation, because the MSS functions take actuator commands, which are set to zero.

| File | Origin | How it was made | sha256 |
|---|---|---|---|
| `spheroid_matlab_reference_mss_current.csv` | `test/plant/auv_spheroid/test_mss_matlab/spheroid_auv_dynamics_full_debug_mss_current.csv` (sha256 `f490c79bc7b23fbd3aa13d3de02480d13d3cff762375e91026064e8a6c54b5c9`), from `test_dynamics_consistency_mss_current.m` (sha256 `51ca41176581c5d390ac5fc12c66eb82750ea5406193273d3371cb6c17a33c3e`) | same 7 columns and re-extraction command as above, applied to the new CSV | `22a64de9587f77ddecc4fb74cc296d70f7f6f3c376fdfb11365757c19de21a7c` |

**Byte-identical to `spheroid_matlab_reference.csv`.** Current MSS changed nothing on this path: `spheroid.m` still uses rho = 1025 for `m`, `gravity(mu)` is the same, and `gRvect.m` only had its variables renamed (`d7c1e96`). So the legacy hydrostatics reference already reproduces current MSS.

## Readers (job A-30, 2026-10-05) — appended

- `tests/hydrostatics/test_submerged_hydrostatics_block.py`: `spheroid_matlab_reference_mss_current.csv` for the block's G1 and G4, `spheroid_matlab_reference.csv` for the numpy source's G1; a test pins that the two are byte-identical.

## Surface hydrostatics references (job U2a, 2026-10-06) — appended, nothing above changed

**Generated:** 2026-10-06 with MATLAB `26.1.0.3312084 (R2026a) Update 4`, MSS HEAD `ac77394` (clean; the hydrostatics files below are unchanged since `99bf0b3` except `osv.m`, whose lines 79–102 were read). Ledger: `agents-more/30_checks/2026-10-06_U2_hydrostatics_gates.md`. Generator: `test_surface_hydrostatics_mss_current.m` **in this folder** (not in `more_generic_models`, so the repo needs nothing outside itself). Numbers written with `%.17g` (every double round-trips), not `writetable`'s 15 digits.

Re-generation (agents-more rule 9: no machine path; `MSS_DIR` is required, `MATLAB_BIN` is the matlab executable): `MSS_DIR=<MSS checkout> "$MATLAB_BIN" -batch "run('<abs path>/test_surface_hydrostatics_mss_current.m')"`. The random draws are seeded (`rng(20261006,'twister')`); a re-run on 2026-10-06 was byte-identical.

| File | What | sha256 |
|---|---|---|
| `surface_gmtrx_mss_current.csv` | `LIBRARY/modeling/Gmtrx.m` (`14be3e4`) on 50 random small-craft argument sets: columns `nabla, A_wp, GMT, GML, x_F, r_bP_1..3`, then `G_01..G_36` (row-major). `r_bP` off the CO in every case. `Gmtrx.m` fixes `rho = 1025`, `g = 9.81` (lines 25–26). Ranges keep \|G\| ≲ 1e6, where the G1 tolerance (1e-9 absolute) is above double rounding. | `ed0209d192ed211f9f9404cd99e6d1b2c6ffcef004760904155d0d5d9ce9a53f` |
| `surface_chain_mss_current.csv` | the MSS coefficient chain, 22 rows: `kind` 1 = `CRAFT/USV/models/otter.m` (`e1dff2a`) lines 121–194 on 20 payloads (row 1 = A-26 payload `mp = 25`, `rp = [0.05 0 -0.35]`; others `mp` 0–45 kg, random `rp`), read from the function workspace through a run-time instrumented copy (checked equal to the unmodified function, max \|diff\| = 0); `kind` 2 = `mssExamples/exShipHydrostatics.m` (`d7c1e96`, script); `kind` 3 = `CRAFT/SHIP/models/osv.m` (`6b1a4b4`). Columns: inputs `hull_count, L, B_hull, T, nabla, Cw, r_bg_3, y_hull, I_L_factor, x_F, rho, g, mp, rp_1..3`; outputs `A_wp, I_T, I_L, KB, BM_T, BM_L, GM_T, GM_L, G_01..36` (G at the CO). `I_L_factor`: `otter.m` 177 (0.8, per pontoon ×2), `exShipHydrostatics.m` 49 and `osv.m` 93 (0.7). `rho`, `g` for kinds 2–3 are `Gmtrx.m`'s. | `304c0266eacda2ce237345872f94c6612743731c88edc54bc5e4e79983fa5849` |
| `test_surface_hydrostatics_mss_current.m` | the generator (edited for rule 9 after the first run; its re-run output is byte-identical) | `69c0fbe6cd295c13e019a4716b22f05de3e986f3959cf062bd36b1c13430f202` |

Checks: chain row 1's `G` equals the A-26 file `../hydrodynamics/catamaran_matlab_reference_mss_current.csv` (`G_01..36`, constant over its 50 rows) to 0.0; test-side transcriptions of `Gmtrx.m` and of the three chains reproduce both CSVs to ≤ 1e-9 (`test_transcription_of_*`).

Conventions: CO on the waterline, z down (`KG = T - z_g`; `exShipHydrostatics.m` 42 `r_bB = [.. T-KB]`); `x_F`, `r_bP` from the CO in BODY. **`otter.m` line 178 computes KB with `L*B_pont`, not the pontoon waterplane area** (`exShipHydrostatics.m` 36 and `osv.m` 85 use `nabla/Awp`): kind-1 rows carry that form; the block reproduces them on its `center_of_buoyancy_area="length_times_beam"` flag (proposed D-MSS-3, ledger).

Readers: `tests/hydrostatics/test_surface_hydrostatics_block.py` (both CSVs and the A-26 catamaran file for the cross-check).
