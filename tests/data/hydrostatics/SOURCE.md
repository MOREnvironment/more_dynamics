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
