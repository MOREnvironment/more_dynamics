# Frozen reference — fin parts and fin compositions (skeleton `LiftingFin`, set `ForceProducerSet`)

**Frozen** files: a changed byte needs a decision record; a new MSS revision that changes an output adds a file `*_mss_<revision>.csv` beside the old one. The generator lives here, finds MSS only through `MSS_DIR` (MATLAB through `MATLAB_BIN`) and writes under `OUT_DIR` (default `<tempdir>/more_mss_references/`), never into this folder; a re-run is checked with `cmp` against the file of the same name here. `tests/force_producers/fin/test_fin_parts_reference_data.py` checks every hash below.

## Files (made 2026-10-08; MSS at `cc07579b63630342802a198c34c74ac7e096d5fb`)

| File | What it holds | Made by | sha256 |
|---|---|---|---|
| `remus100_fins_mss_cc07579.csv` | 1024 × 49: MSS `remus100.m`'s own workspace on 24 structured cases (no current; speed zero, ahead 1.5, astern −1.0, sideways 0.4, vertical 0.3 m/s and one general state; fin commands 0, ±20°, ±0.6 rad, (−0.2, 0.1)) and 1000 seeded cases (seed 20261008; current 0.3 m/s at 30°, w_c = 0.1 m/s; fin commands in [−0.6, 0.6] rad). Propeller command 0 in every case, so `X_prop = K_prop = 0` and `tau` is the fins' share. Columns: `case_id`, `x1..x12`, `ui1..ui3`, `Vc`, `betaVc`, `w_c`, `nu_r_01..06`, `delta_r`, `delta_s`, `U_rh`, `U_rv`, `X_r`, `X_s`, `Y_r`, `Z_s`, `X_prop`, `K_prop`, `tau_01..06`, `rho`, `delta_max`, `A_r`, `A_s`, `CL_delta_r`, `CL_delta_s`, `x_r`, `x_s`. Instrumented copy equal to the unmodified function (xdot, M) to 0 on every case; re-run byte-identical (2026-10-08). | `generate_remus100_fins_mss.m` (MATLAB R2026a Update 4) | `f5dc8cfd01487025d66cfbe1cf3a91f952aa71d43908237d247c0303a012dc0f` |
| `generate_remus100_fins_mss.m` | the generator (runs `remus100.m` as written; a copy instrumented at run time in tempdir only returns its workspace; no MSS line typed) | — | `cbb83e1010c81429114093cbe5db94fe2ef2d674bf8a7549229e67ba79c5e87e` |
| `prestero_2001_remus_fins.json` | numbers typed from the printed pages of Prestero (2001), each with its table and page: ρ (Table A.1, p. 102); S_fin, b_fin, x_finpost, AR_e, ã, c_Lα (Table A.5, p. 103); the twelve control-fin coefficients (Table C.10, p. 111); half a unit of the last printed digit | typed 2026-10-08 from the saved PDF (`source-sim/0_literature/1_vehicle_models/Prestero_2001_REMUS_6DOF_Model_Verification_MIT_WHOI.pdf`, PDF pages 103, 104, 112) | `fa1f9e0df946d877601e27b951df88f84aba013b9bef8952d0958c8eae967694` |

Also read by these tests, frozen in the parent folder: `../parameter_sets.json` (`vsim_fins`: DUNE's LAUV fin forces and positions, `etc/common/vsim-models.ini` 54–62) for the set of one-fin leaves.

## Typed numbers and their pages (rule: no number from memory)

| Value | Printed | Where |
|---|---|---|
| water density | `+1.03e+003` kg/m³ | Prestero 2001, Table A.1, p. 102 |
| fin planform area S_fin | `+6.65e-003` m² | Table A.5, p. 103 |
| fin span b_fin | `+8.57e-002` m | Table A.5, p. 103 |
| fin post x_finpost | `-6.38e-001` m | Table A.5, p. 103 |
| effective aspect ratio AR_e | `+2.21e+000` | Table A.5, p. 103 |
| lift-slope parameter ã | `+9.00e-001` | Table A.5, p. 103 |
| fin lift slope c_Lα | `+3.12e+000` /rad | Table A.5, p. 103 |
| Y_uudr, Z_uuds, M_uuds, N_uudr, Y_uvf, Z_uwf, Y_urf, Z_uqf, M_uwf, N_uvf, M_uqf, N_urf | +9.64, −9.64, −6.15, −6.15, −9.64, −9.64, +6.15, −6.15, −6.15, +6.15, −3.93, −3.93 | Table C.10, p. 111 |
| Zeefakkel actuator: time constant 3 s, ±35°, ±7 and ±10 °/s | text | Murray-Smith 2016, p. 246 (typed in `test_fin_servo_parts.py`) |
| DUNE servo simulator defaults: 90°, 333.3 °/s | `defaultValue("90.0")`, `defaultValue("333.3")` | DUNE `src/Simulators/Servos/Task.cpp` 96–108 @ `555ef0b` (typed in `test_fin_servo_parts.py`) |
| Sarhadi 2026 output-saturated actuator block (lag → rate limit → integrator → amplitude limit on the output, state not clamped) | Fig. 4, p. 4 | Sarhadi, P. (2026). Simple yet effective anti-windup techniques for amplitude and rate saturation: an AUV case study. arXiv:2601.01302v2 (typed in `test_fin_servo_output_saturated.py`) |
| Sarhadi 2026 / Murray-Smith 2016 wind-up release times (0.667 s, 5.00 s) | scratch | `agents-more/60_working/A49_scripts/a49_envelope_output.txt`, `SERVO2` lines (own construction, not published; typed in `test_fin_servo_output_saturated.py`) |

**Known inconsistency of the source (not a typing slip):** Table C.10's scale `Y_uudr = 9.64` does not follow from eq. 4.44 (`ρ c_Lα S_fin`) with Tables A.1 and A.5 (1030 × 3.12 × 0.00665 = 21.4, or 10.7 with the ½ of eq. 4.43); the ratios of the twelve entries do follow from x_finpost. The tests use the table's ratios only and pin the scale difference.

## Layout and conventions

- BODY axes FRD; forces N, moments N·m about the CO; angles rad; `nu_r` relative to the water.
- MSS sign of the fins: a positive rudder angle gives −Y (lift axis −y), a positive stern-plane angle gives −Z (lift axis −z). Prestero's rudder lifts along +y, his stern planes along −z (eq. 4.43).
- Numbers: `writetable` (15 significant digits); JSON numbers as typed.
