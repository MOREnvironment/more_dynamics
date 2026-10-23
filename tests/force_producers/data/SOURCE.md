# Frozen reference — force producers (differential thruster, propeller, Wageningen K_T/K_Q, fins, outboard)

**Frozen** files: a changed byte in a CSV needs a decision record; a new MSS revision that changes an output adds a file `*_mss_<revision>.csv` beside the old one. Every generator lives in this repository, finds MSS only through `MSS_DIR` (MATLAB through `MATLAB_BIN`) and writes under `OUT_DIR` (default `<tempdir>/more_mss_references/<block>/`), never into this folder; a re-run is checked with `cmp` against the file of the same name here.

## Files (rewritten 2026-10-07; MSS checked at `cc07579b63630342802a198c34c74ac7e096d5fb`)

| File | What it holds | Made by | Made at MSS | At MSS `cc07579` | sha256 |
|---|---|---|---|---|---|
| `differential_thruster_mss_current.csv` | 1201 × 9: `time, tau_X_cmd, tau_N_cmd, n_cmd_left, n_cmd_right, thrust_left, thrust_right, tau_X, tau_N` — a requested-wrench sweep (120 s, 0.1 s, ramps, steps and reversals) allocated as `SIMotter.m` does (`Binv = invQR(B_prop)`, `n_c = sign(u).*sqrt(abs(u))`), saturated and passed through `otter.m`'s thrust law; `B_prop`, `n_min`, `n_max` from `otter()` (g = 9.81: `n_max = 103.930864273980`, `n_min = −101.736665503817` rad/s) | `generate_differential_thruster_mss.m` | `99bf0b3` | identical | `ca6c7640e289667ed0cff1372f84ad34bcd1ceca636ee49cdb1bdff85e66d2da` |
| `remus100_actuators_mss_current.csv` | 1036 × 38: MSS `remus100.m`'s own `X_prop`, `K_prop`, `X_r`, `X_s`, `Y_r`, `Z_s` and `tau` (read from its workspace, instrumented copy equal to the unmodified function to 0) on 36 structured + 1000 seeded cases (seed 20261006; n ∈ [−2000, 2000] rpm, fins ∈ [−0.6, 0.6] rad, current 0.3 m/s at 30 deg, w_c = 0.1 m/s) | `generate_remus100_actuators_mss.m` | `ac77394` | identical | `d392eee824970e63edb18e1a1863453d1546e6c9270ab9566351c1b06b40a474` |
| `wageningen_kt_kq_mss_current.csv` | 3756 × 6: `PD, AEAO, z, J, KT, KQ` — MSS `LIBRARY/modeling/wageningen.m` on J = linspace(−0.5, 2.0, 251), J = 0.6632 and 1000 seeded J, for three geometries: (1, 0.718, 3) of `remus100.m`, (0.83, 0.55, 4) test construction, (0.83, 0.718, 3) the outboard set | `generate_wageningen_kt_kq_mss.m` | `ac77394` | identical | `4d2abfe3b5b1483bdc20ed19285e54793c4507dc973d402a1c8c25d670869ee7` |
| `differential_thruster_matlab.csv` | the same 9 columns from the **template** thruster (constants written out as `otter.m` / `SIMotter.m` stood in 2026-02, latitude gravity `gravity(63.446827 deg)` in the speed limits): differs from the MSS file only on the 283 saturated rows | `generate_differential_thruster_template.m` | MSS of 2026-02 | identical | `bf434587340cad6f46d68b7febe0f747ca82f7afe3a68ab3a3c9e3a5be3185bc` |
| `propeller_matlab.csv` | 341 × 11: `RPM, U, Ja_all, Xprop, Kprop, tau1..6` — the **template** REMUS 100 linearised propeller (written out from `remus100.m` as of 2026-02) with rho = 1025 (`remus100.m` uses 1026) | `generate_propeller_template.m` | MSS of 2026-02 | identical | `ba30701bfa877a316f6cc548e007f94e430c8153721715198dce677ce48bb012` |
| `parameter_sets.json` | two parameter sets the block tests use as inputs: `outboard_rpm` and `outboard_throttle` (an electric outboard of the Grethe class: T_max 520 N, P_max 2000 W, η 0.6, 1300 rpm, D 0.305 m, P/D 0.83, A_E/A_O 0.718, z 3, at `[−2.3 0 0.55]` m; values unverified, not identified on the boat). The third set, `vsim_fins` (DUNE's simulator fin forces), left with its block on 2026-10-10 (E-111) for `agents-more/30_checks/reference_blocks/dune_vsim_fins/data/` | written 2026-10-07 from the parameter values the tests used until then; `vsim_fins` removed 2026-10-10 | — | — | `93e17e514cc7293714bca1554a59bbcdf9ccf83d1121efae6323da8371cc3f00` |
| `generate_differential_thruster_mss.m`, `generate_remus100_actuators_mss.m`, `generate_wageningen_kt_kq_mss.m` | the three MSS generators | — | — | re-run 2026-10-07 at `72656d1` and `cc07579`: byte-identical outputs | `7da2ea37bdfee91edd4d04c4ae0feb485590841121153652a7dd563505e18c3f`, `950815365197fc8706dc7cacbc591e5bb364ae19fd4ec6ce196e0880f288519c`, `ca27d0ff5813a475fb72d1bfcbbf9d1921840d3ed6876edcbb15b93e363422b1` |
| `generate_differential_thruster_template.m`, `generate_propeller_template.m` | the two template generators (MATLAB, first written 2026-02-19 / 2026-02-20) | — | — | re-run 2026-10-07 at `cc07579`: byte-identical outputs | `6b383b860d1654974cabe585ef9000061b4d749f9e0cd6de9da39fd3694a28c7`, `67bc3c353f6d70ed0ea902d0b0d4bd4897e95bbcc3067edd39161be8247e99c0` |

## What has no independent frozen reference (stated, not hidden)

- **VSIM fins**: the block and its tests left the library on 2026-10-10 (E-111; DUNE's simulated parts are comparison blocks only); signature and physical tests only, in `agents-more/30_checks/reference_blocks/dune_vsim_fins/`.
- **Outboard throttle mode**: no frozen reference of its force law yet; the rpm mode's K_T/K_Q is G1 against `wageningen_kt_kq_mss_current.csv` (outboard geometry rows) and its force law is written in the test from the stated equations.

## Layout and conventions

- Thruster index 0 = left (port, y < 0), 1 = right; shaft speeds rad/s (Otter) or rpm (REMUS, outboard), as named in each file's columns; forces N, moments N·m.
- `remus100.m`'s fin and propeller constants are MSS's own (rho 1026, fin limit 20 deg, x = −a); the template propeller file uses rho 1025.
- Numbers: `writetable` (15 significant digits) in the MSS and template sweep files; JSON numbers as written by Python's `json` (exact doubles).

## MSS line numbers

Written at the revision named beside them. From `72656d1` to `cc07579`: `otter.m` +1 from line 79 (thrust law 223–227 → 224–228, wrench 231 → 232); `SIMotter.m`, `remus100.m`, `wageningen.m`, `WageningData.mat`, `gravity.m` unchanged.

## Readers

`tests/force_producers/test_force_producers_block.py`. `differential_thruster_mss_current.csv` is also byte-copied into `more_control`'s `tests/data/allocators/` (that repository is cloned alone).

## Re-pin to MSS `99fdf4d` (2026-10-23, MATLAB R2026a Update 4; `99fdf4dbdfe0b3f124e26a689c6aabe0beb5654f`, 41 commits after `cc07579`)

`remus100_actuators_mss_current.csv` rewritten in place: `remus100.m` reads `rho` 1025 (was 1026, `23c307e`), so the propeller force and moment move by 9.7e-4 relative (max |Delta X_prop| 0.113 N). `differential_thruster_mss_current.csv`, `wageningen_kt_kq_mss_current.csv` and the two template files are byte-identical. Note: `SIMotter.m` now allocates reverse thrust with `k_neg` (2026-10-07 revision, "Corrected reverse-thrust control allocation"); the allocation of `generate_differential_thruster_mss.m` still transcribes the earlier `sign(u) .* sqrt(abs(u))` and is not changed by this re-pin (reported in the ledger).
