# Reference — higher-fidelity forms of the torpedo-AUV class on the REMUS 100 hull

**Frozen:** `remus100_fidelity_turn_and_stability.csv` (a changed byte needs a decision record). Its generator, `generate_remus100_fidelity_references.py`, lives here, writes under `OUT_DIR` (default `<tempdir>/more_references/vehicles/remus100_fidelity/`) and never into this folder; a re-run is checked against the file here with `cmp`.

Made with numpy (workspace virtual environment, Python 3.12) and `more_transformations`; a second run wrote the file byte-identical. The cases `published_forms`, `munk_kept_alone` and `geometric_fin_slope_effective_angle` use the base propeller roll scale 0.1; the cases of the part "REMUS 100 (estimated)" use its roll scale 0.35, which moves only the 4-degree-turn columns.

- `remus100_fidelity_turn_and_stability.csv` `3b95ce7f7c3bda95bbbd8e8b47e246ea64da128ac83ba9ba0329ef98e2d785fc`

Run: `python generate_remus100_fidelity_references.py` from this folder (the workspace environment, `more_transformations` importable).

## Files

| File | What it holds |
|---|---|
| `reference_torpedo_auv.py` | numpy reference functions, one per form (equations cited beside each), and `TorpedoAuv`, a whole vehicle from them; the straight-running equilibrium (Newton), the shaft speed that holds a given speed (`shaft_speed_holding`), the motion-state Jacobian (central differences, step 1e-7) and the sway–yaw roots (Fossen 2011, eqs. 12.17–12.21, p. 351) |
| `fidelity_parts.py` | the parameter sets: the named part "REMUS 100 (estimated)" and the comparison cases, all built from the one-value base set `../remus100/remus100_parameters_consistent.json` plus the part file |
| `remus100_estimated_parameters.json` | the named part "REMUS 100 (estimated)": only what differs from the base set, every value with its kind and source; the cruise shaft speed as a derived value |
| `generate_remus100_fidelity_references.py` | the generator of the CSV |
| `remus100_fidelity_turn_and_stability.csv` | per case: the straight-running equilibrium at the cruise shaft speed and zero fins (`cruise_u`, `equilibrium_residual`), the two sway–yaw roots (real, imaginary; largest real part first), the largest real part over all motion states, and the 4° rudder turn (positive command held from t = 0, 60 s from 1.54 m/s, RK4 h = 0.05 s): the mean yaw rate over the last 10 s and its spread, rad/s |

## Cases (one mechanism changed at a time; the pressure centre at 0.60 / 0.65 / 0.70 of the length)

| Case | Forms | Sway–yaw roots (1/s) | 4° turn (°/s) |
|---|---|---|---|
| `published_forms` | every form of the published reference model (shortcut Coriolis, wing lift and drag, linear roll damping), its fin pairs (deflection only, lift coefficients 0.5 / 0.7) | −0.052, −1.000 | +1.8567 |
| `munk_kept_alone` | the same with every added-mass Coriolis term kept | **+1.672**, −2.724 | +26.116 (a spin) |
| `geometric_fin_slope_effective_angle` | published hull forms, four effective-angle fins, lift slope from the fin geometry | −1.017, −3.359 | +1.0323 |
| `paired_published_fin_slope` | Munk kept, point-force body lift, effective-angle fins with the published lift coefficients | **+0.660 / +0.570 / +0.476** | +15.86 / +14.14 / +12.38 (spins) |
| `paired_geometric_fin_slope` | the named part: Munk kept, point-force body lift, effective-angle fins, lift slope from the fin geometry, roll scale 0.35 | −1.116, −3.962 / −1.170, −3.907 / −1.226, −3.851 | **+6.305 / +5.909 / +5.525** |

The named part declares the fraction 0.65 (Prestero 2001, eq. 4.35, p. 31); 0.60 and 0.70 bound Hoerner's range as Prestero gives it (p. 31) and are run as a sensitivity, not as a tolerance.

The cruise shaft speed (920.0884 rpm) is the one at which the reference model with the published forms runs straight at 1.54 m/s (Prestero 2001, Table 8.1, p. 60); `shaft_speed_holding` finds it by bisection and a test checks that it still does. A comparison with a run of another simulator is a local check, not part of this folder's tests.

**What these numbers are and are not.** A second implementation of the stated equations on one synthetic parameter set (the REMUS 100 hull of the base set with Prestero's fin geometry): a form check and one straight-line stability check at one trim speed, open loop, still water. They are not a validation of any vehicle. Prestero's field record (Fig. 7-7, p. 57: about 10 °/s after a 4° rudder, five seconds, depth under closed-loop control, fin misalignment up to 5°, p. 61) is not reproduced and is not a gate.

## Values (published, derived or estimated; none fitted here)

| Quantity | Value | Source | Kind |
|---|---|---|---|
| hull, mass, site, hull loads, propeller | the one-value base set | `../remus100/SOURCE.md` | published / derived there |
| axial drag coefficient c_d | 0.42 | Allen et al. (2000), cited second-hand, not read | published |
| body-lift slope per degree c_yβ | 0.003 /deg | Prestero 2001, eq. 4.31, p. 30 (after Hoerner 1965, p. 13-3, not read; valid 6.7 ≤ l/d ≤ 10; here l/d = 8.42) | published |
| body-lift slope c_ydβ | 1.4475 /rad | (l/d) c_yβ 180/π, Prestero 2001, eq. 4.32 | derived |
| pressure-centre fraction | 0.65 (0.60, 0.70 as sensitivity) | Prestero 2001, eq. 4.35, p. 31; Hoerner's range through Prestero p. 31 | published |
| pressure centre x_cp | −0.240 m at 0.65 | L/2 − f L, the hull centred on the body origin | derived |
| fin area S_fin (one fin) | 0.00665 m² | Prestero 2001, Table 2.2, p. 17 | published |
| fin span b_fin | 0.0857 m | Prestero 2001, Table 2.2, p. 17 | published |
| aspect-ratio factor k_AR | 2 | Prestero 2001, eq. 4.39, p. 31 | published |
| empirical factor ā | 0.9 | Prestero 2001, Table 2.2, p. 17, eq. 4.38, p. 31 (after Hoerner 1965 p. 3-2, not read) | published |
| fin lift slope c_Lα | 3.1158 /rad per fin | eqs. 4.38–4.39 (AR_e = 2 b²/S = 2.209); Table 2.2 prints 3.12 and 2.21 | derived |
| fin positions | (−L/2, 0, 0) | the fin post at the tail of the hull centred on the origin | derived |
| propeller roll scale | 0.35 | Prestero 2001, Sec. 4.6.2, p. 33; calibrated in the generic-models REMUS actuator set to a steady roll of −5.3° at 1500 rpm and 1.51 m/s; not re-calibrated here | estimate |
| cruise shaft speed | 920.0884 rpm | `shaft_speed_holding` on the published-forms case at 1.54 m/s | derived |

## Published anchors used by the tests (typed with their page)

Prestero 2001: Table 2.1 p. 15 (l = 1.33 m, d = 0.191 m), Table 2.2 p. 17 (c_Lα 3.12, AR_e 2.21), Table A.1 p. 102 (c_ydβ 1.20, x_cp −0.321 m), Table D.1 p. 112 (M_wl / Z_wl = −11.1 / −34.5); the copy read is the filed thesis (the signed archival scan; a 2002 LaTeX file prints three other values elsewhere in the tables, none used here). The text of eqs. 4.5–4.6, 4.10–4.11 and 4.27–4.39 was read from the text layer of the filed thesis; Tables A.1 and D.1 have no text layer and are typed from the scan.

## Not in this reference (not derived or not settled; stated, never guessed)

- The hull-profile forms (Myring, tabulated offsets), strip-theory added mass and the profile cross-flow damping.
- A part on Prestero's own geometry and coefficients (his propulsion as a constant thrust has no form here, his weight is printed in N with no g, his buoyancy is printed as 306 N and 308 N).
- The sign of the propeller roll reaction against the physical shaft-torque direction: the reference applies `-scale K` (Prestero 2001, Sec. 4.6.2, p. 33) as the plugin does; it does not enter the straight-running yaw-sway block.
