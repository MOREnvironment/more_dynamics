# Frozen reference — force-producer blocks

**Frozen:** 2026-10-05 (job A-22, verifier). **Append-only:** a changed byte in either CSV needs an ADR.

## Files

| File | Origin (relative to `more_generic_models/more_generic_models/`) | How it was made | sha256 |
|---|---|---|---|
| `differential_thruster_matlab.csv` | `test/thruster/thruster_differential/results_thruster_test_matlab.csv` | byte copy | `bf434587340cad6f46d68b7febe0f747ca82f7afe3a68ab3a3c9e3a5be3185bc` |
| `propeller_matlab.csv` | `test/thruster/thruster_wagenigen/compare/remus_propulsion_full.csv` | byte copy | `ba30701bfa877a316f6cc548e007f94e430c8153721715198dce677ce48bb012` |

## Revisions and generators

- `more_generic_models`: HEAD `e6814a0f80bc140ad9bacea975ae275ddff49bbb`.
- `differential_thruster_matlab.csv`: clean at `f953733` (2026-02-19). Written by `test/thruster/thruster_differential/test_diff_dynamical.m` (sha256 `632676cf…40156780`, `writetable` call after the loop). That script does **not** define `g`, `k_pos`, `k_neg`, `n_max`, `n_min`, `l1`, `l2` or `Binv`; it relies on the MATLAB workspace left by `test_diff_trhuster.m` (sha256 `43c2c381…01f9024`), which sets `g = gravity(deg2rad(63.446827))`, `k_pos = 0.02216/2`, `k_neg = 0.01289/2`, `n_max = sqrt((0.5*24.4*g)/k_pos)`, `n_min = -sqrt((0.5*13.6*g)/k_neg)`, `l1 = -0.395`, `l2 = 0.395`. Confirmed numerically: the largest `n_cmd` in the CSV equals that `n_max` (103.993272201449 rad/s); with `otter.m`'s `g = 9.81` it would not.
  Columns: `time, tau_X_cmd, tau_N_cmd` (the allocation request), `n_cmd_left, n_cmd_right` (shaft speed, rad/s, after saturation), `thrust_left, thrust_right` (N), `tau_X` (N), `tau_N` (N·m). No shaft lag: thrust is evaluated on the commanded speed.
- `propeller_matlab.csv`: clean at `d31ffb1` (2026-02-20). Written by `test/thruster/thruster_wagenigen/test_mss_reference.m` (sha256 `73231efd…2c73`), local `remus_propulsion_reference` (lines 75–136): the `remus100.m` linearised propeller with **rho = 1025** (line 83; current `remus100.m` line 99 has 1026) and the printed `KT_0, KQ_0, KT_max, KQ_max` (lines 98–102). Sweep: rpm −1500…1500 step 100, speed `U` 0…2.5 step 0.25 (lines 14–15), 341 rows. Columns `RPM, U, Ja_all, Xprop` (already × (1 − t_prop)), `Kprop` (already / 10), `tau1..tau6`. `Ja_all` uses `|n|` (line 106) and is not an input of the force.
- MSS (`source-sim/MSS`, HEAD `99bf0b3`): `otter.m` `e1dff2a`; `remus100.m` `1264bbd`; `LIBRARY/modeling/wageningen.m` and `utiles/WageningData.mat` / `.txt` `ac792a9`; `INS/functions/gravity.m` `6d4b553`.

## What has no frozen file, and why

- **Fins (rudder and stern planes):** `test/fins/fins_auv_physical/test_mss_matlab/test_fins_remus1000.m` is a copy of `remus100.m` that writes no CSV. The tests transcribe `remus100.m` lines 110, 114–115, 180–190, 232–252 instead, each line cited and pinned.
- **VSIM fins, outboard motor:** no MATLAB file exists; the reference is the numpy source (G2). For the outboard, the reverse-thrust branch is tested **after** the sign fix of owner decision E-20 Q5 a: the fixed force is the negative of the current source's force for `n < 0` (derivation in the test file).
- **Wageningen K_T/K_Q:** the reference is MSS's own copy of the regression table, `WageningData.mat` (read with `scipy.io.loadmat`) and `WageningData.txt`. The paper both cite (Barnitsas, Ray and Kinley 1981, as MSS spells it) is not on disk and was not read.

## Current-MSS regeneration (job A-26, 2026-10-05) — appended, nothing above changed

**Generated:** 2026-10-05 with MATLAB `26.1.0.3312084 (R2026a) Update 4`, MSS `source-sim/MSS` HEAD `99bf0b30e9ae0dca3515d02e0bbd06c54cc0c2f7` (clean). Ledger: `agents-more/30_checks/2026-10-05_A26_matlab_references_current_mss.md`. The `*_mss_current.csv` files are the **current-MSS (default) path**. The legacy files above stay as they are (template path).

How they were made: new generators named `*_mss_current.m` sit beside the old ones in `more_generic_models` and are not committed. They run MSS's own `remus100.m` / `otter.m` on the same `inputs.csv`. A copy of the MSS file is made in `tempdir` at run time; it adds only one extra output, the function's workspace, so the intermediate terms can be read. On every case that copy equals the unmodified MSS function (`xdot`, `M`): max |diff| = 0. External `tau` is added in the model's own state equation, because the MSS functions take actuator commands, which are set to zero.

| File | Origin | How it was made | sha256 |
|---|---|---|---|
| `differential_thruster_mss_current.csv` | `test/thruster/thruster_differential/results_thruster_test_mss_current.csv`, from `test_diff_thruster_mss_current.m` (sha256 `7fc582b0f69d10b9ee86d6b38265ccceb27d2cfe602fa1b54b467ca323403f21`) | byte copy | `ca6c7640e289667ed0cff1372f84ad34bcd1ceca636ee49cdb1bdff85e66d2da` |

The generator is self-contained, which fixes A-22 finding 3:

- `B_prop`, `n_min` and `n_max` come from `[~,~,~,B_prop,n_min,n_max] = otter()`, so `g = 9.81`, `n_max = 103.930864273980` rad/s and `n_min = -101.736665503817` rad/s.
- Allocation follows `SIMotter.m:95, 183–184` (`invQR`, `sign(u).*sqrt(abs(u))`).
- Saturation, thrust and `tau` are read from `otter.m` itself (lines 220–232).
- The sweep is `test_diff_dynamical.m` lines 4–30, unchanged.
- Columns are the same as `differential_thruster_matlab.csv`.

It differs from the legacy file **only on the 283 saturated rows**: the limits use g = 9.81 instead of the latitude gravity. The other 918 rows agree to 1.1e-13.

## Readers (job A-30, 2026-10-05) — appended

- `tests/force_producers/test_force_producers_block.py`: `differential_thruster_mss_current.csv` for the default (`otter.m` g = 9.81, `test_G1_MSS_differential_default_equals_matlab`); `differential_thruster_matlab.csv` for the numpy source's G1 and the block's G1/G4 on the latitude gravity; `propeller_matlab.csv` for the propeller's G1/G4 (no current-MSS file exists).

## Current-MSS fins, propeller and Wageningen (job U4a, 2026-10-06) — appended, nothing above changed

**Generated:** 2026-10-06 with MATLAB `26.1.0.3312084 (R2026a) Update 4`, MSS HEAD `ac77394b74a2184317e92c6b73d8deca97b48359` (clean; no actuator line differs from `99bf0b3`: the diff touches `cylinderDrag.m`, `Dmtrx.m`, `XuuITTC.m` and ship/hydro files only; `remus100.m` last changed in `1264bbd`, `wageningen.m` and `WageningData.mat` in `ac792a9`, `sat.m` read at `ac77394`). Ledger: `agents-more/30_checks/2026-10-06_U4_force_producers_gates.md`. Closes A-30 finding 3 (fins and propeller had no current-MSS file).

| File | Origin (in the `more_generic_models` repository) | How it was made | sha256 |
|---|---|---|---|
| `remus100_actuators_mss_current.csv` | `more_generic_models/test/plant/auv_spheroid/test_mss_matlab/remus100_actuators_mss_current.csv`, from `test_remus100_actuators_mss_current.m` (sha256 `28bd87342dae56d1c92c9af3326003124c2087e563c4669fa27887909e3b59d3`) | byte copy | `789be46937db1412855ff005113b40c8c7c2e31575a6c7f19cef6c1a87ce63de` |
| `wageningen_kt_kq_mss_current.csv` | `more_generic_models/test/thruster/thruster_wagenigen/wageningen_kt_kq_mss_current.csv`, from `test_wageningen_mss_current.m` (sha256 `3ba1dd06d6d346a86d8a064a3b08aae8d6d69c4797bcb21f272e428ea36b526e`) | byte copy | `4d2abfe3b5b1483bdc20ed19285e54793c4507dc973d402a1c8c25d670869ee7` |

- **`remus100_actuators_mss_current.csv`** (1036 × 38). `remus100.m` run as written; an instrumented copy (tempdir, A-26 method) returns its workspace; instrumented vs unmodified `xdot`, `M`: max |diff| = 0 on every case. The actuator lines it reads: propeller 147–178 (`D_prop`, `t_prop`, `Va = 0.944 U_r`, `Ja_max`, printed `KT_0 … KQ_max` 158–162, the linearised law 168–178), saturation 110–116 (`delta_max = deg2rad(20)`, `n_max = 1525`), fins 180–190 (`S_fin`, `CL_delta_r/s`, `A = 2 S_fin`, `x_r = x_s = -a`, `a = 1.0096 L/2` line 135), forces 231–252 (`U_rh`, `U_rv`, `X_r`, `X_s`, `Y_r`, `Z_s`, `tau`), `rho = 1026` line 99. Inputs (test constructions, seed `20261006`, MATLAB `twister`): 36 structured cases at one state (n ∈ {−2000, −1525, −800, −1e−3, 0, 1e−3, 800, 1525, 2000} rpm × fin pairs {(0,0), (±20°), (0.6, −0.6), (−0.2, 0.1)} rad) + 1000 random (u ∈ [−1, 3], v, w ∈ [−1, 1] m/s, p, q, r ∈ [−0.5, 0.5] rad/s, φ, θ ∈ [−0.5, 0.5], ψ ∈ [−π, π], fins ∈ [−0.6, 0.6] rad, n ∈ [−2000, 2000] rpm); current `Vc = 0.3` m/s, `betaVc = 30°`, `w_c = 0.1` m/s. 42 % of the fin commands and 24 % of the speeds lie beyond the limits. Columns: `case_id`, `x1..x12`, `ui1..ui3` (commands), `nu_r_01..06`, `U_r`, `n_p` (rps, saturated), `delta_r`, `delta_s` (saturated), `X_prop`, `K_prop`, `X_r`, `X_s`, `Y_r`, `Z_s`, `tau_01..06` (`remus100.m`'s actuator `tau`, before its state equation).
- **`wageningen_kt_kq_mss_current.csv`** (3756 × 6). MSS `LIBRARY/modeling/wageningen.m` called as written on J = `linspace(-0.5, 2.0, 251)`, J = 0.6632 and 1000 seeded J in [−0.5, 2.0] (seed `20261006`), for three geometries (P/D, A_E/A_O, z): (1, 0.718, 3) `remus100.m` line 157; (0.83, 0.55, 4) test construction; (0.83, 0.718, 3) the outboard parameter set of `electrical_outboard_motor_params.py`. Columns `PD, AEAO, z, J, KT, KQ`.
- **How to regenerate (agents-more rule 9: nothing relative to one machine).** Both scripts find MSS only through `MSS_DIR` and stop with an error naming it when it is unset: `MSS_DIR=<MSS root> "$MATLAB_BIN" -batch "run('<script>')"`. Regenerated that way on 2026-10-06: both CSVs byte-identical to the files above.
- **Readers:** `tests/force_producers/test_force_producers_block.py` section 7. The test file finds MSS and `more_generic_models` only through `MSS_DIR` and `MORE_GENERIC_MODELS_DIR`; unset, the pinned-line check and the two `WageningData` readers skip naming the variable, and every CSV gate runs on the pinned line text frozen in `CITED_LINES`.
