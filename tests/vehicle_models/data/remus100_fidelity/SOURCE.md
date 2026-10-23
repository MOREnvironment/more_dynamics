# Reference — higher-fidelity forms of the torpedo-AUV class on the REMUS 100 hull

**Frozen:** `remus100_fidelity_turn_and_stability.csv` (a changed byte needs a decision record). Its generator, `generate_remus100_fidelity_references.py`, lives here, writes under `OUT_DIR` (default `<tempdir>/more_references/vehicles/remus100_fidelity/`) and never into this folder; a re-run is checked against the file here with `cmp`.

Made 2026-10-13 with numpy (workspace virtual environment, Python 3.12) and `more_transformations` on its branch `casadi-transforms`; a second run the same day wrote the file byte-identical.

- `remus100_fidelity_turn_and_stability.csv` `ea5e3eaec206f6679bc8c9a00405d1749fb4ebd533f17af5fbaff6a921e95086`

Run: `python generate_remus100_fidelity_references.py` from this folder (the workspace environment, `more_transformations` importable).

## Files

| File | What it holds |
|---|---|
| `reference_torpedo_auv.py` | numpy reference functions, one per form (citations beside each line), and `TorpedoAuv`, a whole vehicle from them; the straight-running equilibrium (Newton), the motion-state Jacobian (central differences, step 1e-7) and the sway–yaw roots (Fossen 2011, eqs. 12.17–12.21, p. 351) |
| `fidelity_parts.py` | the parameter sets: the named part "REMUS 100 (Munk + body lift)" and the comparison cases, all built from the one-value MSS set `../remus100/remus100_parameters_consistent.json` plus the part file |
| `remus100_munk_body_lift_parameters.json` | the named part "REMUS 100 (Munk + body lift)": only what differs from the one-value MSS set, every value with its kind and source; the cruise shaft speed by reference to a frozen MSS cell |
| `generate_remus100_fidelity_references.py` | the generator of the CSV |
| `remus100_fidelity_turn_and_stability.csv` | per case: the straight-running equilibrium at the cruise shaft speed and zero fins (`cruise_u`, `equilibrium_residual`), the two sway–yaw roots (real, imaginary; largest real part first), the largest real part over all motion states, and the 4° rudder turn (positive command held from t = 0, 60 s from 1.54 m/s, RK4 h = 0.05 s, `SIMremus100.m` 60): the mean yaw rate over the last 10 s and its spread, rad/s |

## Cases (one mechanism changed at a time; the pressure centre at 0.60 / 0.65 / 0.70 of the length)

| Case | Forms | Sway–yaw roots (1/s) | 4° turn (°/s) |
|---|---|---|---|
| `mss_forms` | every MSS form, MSS's fin pairs (deflection only, calibrated 0.5 / 0.7) | −0.052, −1.000 | +1.8574 = the frozen MSS run `rudder_step_4deg_at_1p54ms_consistent.csv` (last 10 s, ≤ 1e-9 rad/s) |
| `munk_kept_alone` | MSS forms, every added-mass Coriolis term kept | **+1.672**, −2.724 | +26.137 (a spin) |
| `geometric_fin_slope_effective_angle` | MSS hull forms, four effective-angle fins, lift slope from the fin geometry | −1.016, −3.362 | +1.0329 |
| `paired_mss_fin_slope` | Munk kept, point-force body lift, effective-angle fins with MSS's calibrated coefficients | **+0.659 / +0.570 / +0.476** | +15.85 / +14.13 / +12.38 (spins) |
| `paired_geometric_fin_slope` | the named part: Munk kept, point-force body lift, effective-angle fins, lift slope from the fin geometry | −1.117, −3.963 / −1.171, −3.909 / −1.227, −3.852 | **+6.3051 / +5.9095 / +5.5260** |

The named part declares the fraction 0.65 (Prestero 2001, eq. 4.35, p. 31; the owner's choice of the one value, 2026-10-13); 0.60 and 0.70 bound Hoerner's range as Prestero gives it (p. 31) and are run as a sensitivity, not as a tolerance.

**What these numbers are and are not.** A second implementation of the stated equations on one synthetic parameter set (the MSS REMUS geometry with Prestero's fin geometry): a form check and one straight-line stability check at one trim speed, open loop, still water. They are not a validation of any vehicle. Prestero's field record (Fig. 7-7, p. 57: about 10 °/s after a 4° rudder, five seconds, depth under closed-loop control, fin misalignment up to 5°, p. 61) is not reproduced and is not a gate.

## Values (all published or derived; none fitted)

| Quantity | Value | Source | Kind |
|---|---|---|---|
| hull, mass, site, hull loads, propeller | the one-value MSS set | `../remus100/SOURCE.md` | published / derived there |
| body-lift slope per degree c_yβ | 0.003 /deg | Prestero 2001, eq. 4.31, p. 30 (after Hoerner 1965, p. 13-3, not read; valid 6.7 ≤ l/d ≤ 10; here l/d = 8.42) | published |
| body-lift slope c_ydβ | 1.4475 /rad | (l/d) c_yβ 180/π, Prestero 2001, eq. 4.32 | derived |
| pressure-centre fraction | 0.65 (0.60, 0.70 as sensitivity) | Prestero 2001, eq. 4.35, p. 31; Hoerner's range through Prestero p. 31 | published |
| pressure centre x_cp | −0.240 m at 0.65 | L/2 − f L, the hull centred on the body origin (`remus100.m` 131, 134) | derived |
| fin area S_fin (one fin) | 0.00665 m² | Prestero 2001, Table 2.2, p. 17 (= `remus100.m` 179) | published |
| fin span b_fin | 0.0857 m | Prestero 2001, Table 2.2, p. 17 | published |
| empirical factor ā | 0.9 | Prestero 2001, Table 2.2, p. 17, eq. 4.38, p. 31 (after Hoerner 1965 p. 3-2, not read) | published |
| fin lift slope c_Lα | 3.1158 /rad per fin | eqs. 4.38–4.39 (AR_e = 2 b²/S = 2.209); Table 2.2 prints 3.12 and 2.21 | derived |
| fin positions | (−L/2, 0, 0) | `remus100.m` 184, 189 | derived |
| cruise shaft speed | 920.088 rpm | `../remus100/remus100_trajectory_rudder_step_4deg_at_1p54ms_consistent.csv`, row 1, `ui_03` (holds 1.54 m/s straight, Prestero 2001, Table 8.1, p. 60) | derived |

## Published anchors used by the tests (typed with their page)

Prestero 2001: Table 2.1 p. 15 (l = 1.33 m, d = 0.191 m), Table 2.2 p. 17 (c_Lα 3.12, AR_e 2.21), Table A.1 p. 102 (c_ydβ 1.20, x_cp −0.321 m), Table D.1 p. 112 (M_wl / Z_wl = −11.1 / −34.5); the copy read is the filed thesis (the signed archival scan; a 2002 LaTeX file prints three other values elsewhere in the tables, none used here).

## Not in this reference (not derived or not settled; stated, never guessed)

- The hull-profile forms (Myring, tabulated offsets), strip-theory added mass and the profile cross-flow damping: their reproductions of Prestero's tables are reported, not checked, and what the profile damping form replaces of the time-constant damping is not defined.
- A part on Prestero's own geometry and coefficients (also: his propulsion as a constant thrust has no form here, his weight is printed in N with no g, his buoyancy is printed as 306 N and 308 N).
- The behaviour of `slender_body_hoerner` for a hull outside 6.7 ≤ l/d ≤ 10 (open; owner, 2026-10-13).
