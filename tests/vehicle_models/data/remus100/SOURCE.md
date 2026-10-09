# Frozen reference — torpedo-AUV class on the REMUS 100 parameter set (MSS `remus100.m`, one value per quantity)

**Frozen** files: a changed byte in a CSV or in the parameter file needs a decision record; a new MSS revision that changes an output adds a file `*_mss_<revision>.csv` beside the old one. The generator lives in this folder, finds MSS only through `MSS_DIR` (MATLAB through `MATLAB_BIN`) and writes under `OUT_DIR` (default `<tempdir>/more_mss_references/vehicles/remus100/`), never into this folder; a re-run is checked with `cmp` against the file of the same name here.

**Two reference sets.** The **default** (the class is gated against it) is the one-value set, `*_consistent.*`, made by `generate_remus100_consistent.m` from the marked copy `remus100_consistent.m`: one water density, one geometry, the mass given (sections below). The **history** set is the output of unmodified MSS `remus100.m`, made by `generate_remus100_mss.m` (`*_mss.csv`, the unsuffixed trajectory files, `remus100_parameters.json`); it is kept, never overwritten, and one test shows by how much it differs. Why: within one model MSS uses two water densities and two geometries (registered as MSS faults; the owner, 2026-10-07: *"we dont match mss if doesnt make sense"*); our model holds one value per physical quantity.

History set: made 2026-10-07 by MATLAB R2026a Update 4 running MSS at `cc07579b63630342802a198c34c74ac7e096d5fb` (`CRAFT/AUV/models/remus100.m`, `CRAFT/AUV/SIMremus100.m`, `LIBRARY/numericalMethods/rk4.m`). Re-run the same day from this folder: all 21 outputs byte-identical.

Default set: made 2026-10-07 by the same MATLAB and MSS revision; a second run the same day wrote all 24 outputs byte-identical.

Run: `MSS_DIR=<MSS checkout> "$MATLAB_BIN" -batch "run('<this folder>/generate_remus100_mss.m')"`.

## The one-value reference (the default)

**The model.** `remus100_consistent.m` is a marked copy of MSS `remus100.m` at `cc07579` (T. I. Fossen, MIT licence, carried in the copy). Every changed line ends with `%<added>` and says `one value:`; a line starting with `%<removed> ` is the original commented out; the generator rebuilds the original from the copy and stops unless it equals MSS byte for byte (so does `test_marked_copies_rebuild_their_mss_files` when `MSS_DIR` is set). Four marked copies of MSS library functions take the quantity they had fixed inside as an input: `spheroid_consistent.m` (the mass; MSS: 1025 kg/m³ × volume, `spheroid.m` 35–36), `imlay61_consistent.m` (ρ; MSS 1026, `imlay61.m` 31), `forceLiftDrag_consistent.m` (ρ; MSS 1026, `forceLiftDrag.m` 26), `crossFlowDrag_consistent.m` (ρ; MSS 1025, `crossFlowDrag.m` 36).

| Quantity | One value | MSS `remus100.m` uses |
|---|---|---|
| water density | **1026 kg/m³** (`remus100.m` 98) in the added mass, hull lift and drag, cross-flow, propeller and fins | 1025 (body mass, cross-flow) and 1026 (added mass, lift and drag, propeller, fins) |
| length, diameter | **L = 1.6 m, D = 0.19 m** (131–132) everywhere: a = L/2, b = D/2 for the mass, the added mass, `CD_0` (143–144) and the fin positions x_r = x_s = −a (184, 189), as for S (133), lift and drag (220) and cross-flow (221) | 1.0096·L/2, 1.0096·D/2 for mass, added mass, `CD_0`, fins (134–135); L, D for S, lift and drag, cross-flow |
| mass | **31.9 kg, given** (line 3: "the mass of the vehicle is 31.9 kg"; the generator reads line 3 and stops unless the copy holds the same number); inertia of a homogeneous spheroid of that mass (`spheroid.m` 39–42) | 1025 kg/m³ × the scaled spheroid volume = 31.9005 kg |
| weight, buoyancy | **W = m g_mu, B = W** as `remus100.m` 214 states it (neutrally buoyant, the REMUS assumption); not B = ρ g ∇, which would be 31.03 kg of water against 31.9 kg of vehicle (0.87 kg heavy) | the same line |

The propeller parameters are named `propeller_diameter` and `max_shaft_speed` after the part declarations; their values and MSS places are unchanged. The parameter set gives the mass through `body_density` = mass / (4/3 π a b²) = 1054.79 kg/m³, the input the rigid-body block's `"spheroid"` form declares (it has no mass input).

**Options of the copy** (6th argument, a struct; the defaults are the one-value model): `opt.munk` keeps the added-mass Coriolis couplings of lines 207–210 (the `xdot_munk` columns); `opt.tau_ext` replaces the actuator wrench by an external one (the bare hull, stage 1; `ui = 0`); `opt.rho`, `opt.mass` change the one density and the mass (the invariance check); `opt.geometry = 'mss'` gives MSS's scaled semi-axes and the mass of the displaced water at the one ρ — MSS with one density only, used for history (`xdot_one_density`, `*_one_density.csv`). The fourth output is the actuator wrench.

**Density invariance (measured, MATLAB).** With one density reaching every call and the mass scaling with it, every force, mass and damping term is proportional to ρ (Dmtrx's roll and pitch frequencies √(W Δz / M) do not depend on it; the cylinder drag coefficient depends on viscosity, not ρ), so ν̇ cannot depend on ρ. On the 640 derivative cases, `max |Δν̇| / max |ν̇|` with (ρ, m) scaled together: **0 exactly** for the factors 0.5 and 2 (powers of two scale every number exactly), **1.6e-14** for 1000/1026 (rounding of ρ and m); with MSS geometry and the mass from ρ, ρ = 1000 against 1026: 1.1e-14; the bare hull with (ρ, m, `tau_ext`) scaled by 1000/1026 on the 600 hull cases: 8.6e-15. Unmodified MSS is not invariant (A-47: 1.6e-3).

**The split is exact:** the hull copy (`opt.tau_ext`) fed with the full copy's own actuator wrench gives the full copy's `xdot` with difference **0** on the 640 cases.

**How far the history set is from the default (measured; `test_history_unmodified_mss_differs_by_one_density_and_one_geometry`).**

| Change | Derivative, 640 cases: max over cases of `max|Δν̇|/max|ν̇|` (and `‖Δν̇‖/‖ν̇‖`) | Turn, 15° rudder, 925 rpm, from 1.54 m/s, 60 s (`rudder_step_15deg_at_925rpm`) |
|---|---|---|
| density alone: MSS → one density (MSS geometry) | 1.88e-3 (1.28e-3) | heading −0.179°, turn rate 0.087240 → 0.087188 rad/s, position 6.2 cm |
| geometry and given mass: one density → one value | 5.85e-2 (4.74e-2) | heading **+2.905°**, turn rate → 0.088027 rad/s, position 0.94 m, speed 1.5078 → 1.5153 m/s |
| total: MSS → one value | 6.0e-2 (4.84e-2) | heading +2.727°, position 0.87 m |

The density row reproduces the A-47 scan (−0.179°, 0.1 %). The geometry row is our choice (the hull's own L and D everywhere, mass given) and is +2.9°, not the −4.06° of A-47's other variant (the *scaled* geometry everywhere).

**Files (default set).**

| File | What it holds |
|---|---|
| `remus100_parameters_consistent.json` | the parameter set by the class's parameter names, one `water_density` (no `cross_flow_water_density`): `value`, `unit`, `block`, `block_parameter`, `declared` (`new`: every block is in the declared-parameter form), `place`, `stage`, `kind` (`published` = a number MSS prints; `derived` = computed from published ones, the formula in `place`); `one_value` (density, geometry, the given mass 31.9 kg and its line, B = W), `gravity`, `selectors` (with `"form": "submerged"` for `linear_damping_parameters` and `"open_water_coefficients": "given"` for `propeller_parameters`: the set gives the printed K_T, K_Q). Every value read from the copy's workspace or MSS text by the generator. |
| `remus100_derivative_consistent.csv` | 640 × 59: the cases of `remus100_derivative_mss.csv` (same seed; the generator stops unless states, inputs and currents are equal): `case_id, x, ui, Vc, betaVc, w_c, tau_01..06` (the copy's actuator wrench), `xdot`, `xdot_munk`, `xdot_one_density` |
| `remus100_hull_derivative_consistent.csv` | 600 × 46: the cases of `remus100_hull_derivative_mss.csv` (checked equal): `xdot`, `xdot_munk` of the hull copy |
| `remus100_mass_matrix_consistent.csv` | `M` of the empty call `[~,~,M] = remus100_consistent()` |
| `remus100_trajectory_<scenario>_consistent.csv` | the 16 scenarios of the history set (same inputs and lengths) plus `rudder_step_15deg_at_925rpm`; the Prestero scenario at **920.088018654 rpm** (the speed holding 1.54 m/s straight in the copy, `fzero`, surge acceleration −6.4e-14 m/s²): steady yaw rate, last 10 s, **1.85743°/s** at 1.535 m/s |
| `remus100_trajectory_tolerances_consistent.csv` | as the history tolerance file, for the 17 default scenarios |
| `remus100_trajectory_rudder_step_15deg_at_925rpm_{mss,one_density}.csv` | the same turn by unmodified `remus100.m` and by the one-density model: history, no gate |
| `generate_remus100_consistent.m`, `remus100_consistent.m`, `spheroid_consistent.m`, `imlay61_consistent.m`, `forceLiftDrag_consistent.m`, `crossFlowDrag_consistent.m` | the generator and the five marked copies |

Run: `MSS_DIR=<MSS checkout> "$MATLAB_BIN" -batch "run('<this folder>/generate_remus100_consistent.m')"`.

## Parameter set — one value per quantity (the default)

| Name | Value | Unit | Block (all in the declared-parameter form) | Place (MSS lines at `cc07579`) | Stage | Kind |
|---|---|---|---|---|---|---|
| `semi_major_axis` | `0.8` | m | rigid_body | remus100.m:131 L_auv / 2 (one geometry; MSS scales by 1.0096 at :134) | hull | published |
| `semi_minor_axis` | `0.095` | m | rigid_body | remus100.m:132 D_auv / 2 (one geometry; MSS scales by 1.0096 at :135) | hull | published |
| `body_density` | `1054.7872613500267` | kg/m^3 | rigid_body | mass 31.9 kg (remus100.m:3, given) / (4/3 pi a b^2); MSS: 1025 at spheroid.m:35 | hull | derived |
| `water_density` | `1026` | kg/m^3 | rigid_body, lift_drag, cross_flow, propeller, fins | remus100.m:98 rho = 1026, the one value (MSS: 1026 at imlay61.m:31 and forceLiftDrag.m:26, 1025 at crossFlowDrag.m:36) | hull | published |
| `roll_added_inertia_ratio` | `0.3` | 1 | rigid_body | remus100.m:136 r44 | hull | published |
| `body_center_of_gravity` | `[0, 0, 0.02]` | m | rigid_body | remus100.m:137 r_bG | hull | published |
| `center_of_buoyancy` | `[0, 0, 0]` | m | submerged_hydrostatics, submerged_linear_damping | remus100.m:138 r_bB | hull | published |
| `weight` | `313.31493718007266` | N | submerged_hydrostatics, submerged_linear_damping | remus100.m:214 W = m * g_mu, m = 31.9 kg given (remus100.m:3), g_mu = gravity(mu) (remus100.m:96-97) | hull | derived |
| `buoyancy` | `313.31493718007266` | N | submerged_hydrostatics | remus100.m:214 B = W (neutral, as REMUS states) | hull | derived |
| `time_constants` | `[20, 20, 1]` | s | submerged_linear_damping | remus100.m:192, 193, 196 [T1 T2 T6] (Dmtrx.m call at :217) | hull | published |
| `damping_ratios` | `[0.3, 0.8]` | 1 | submerged_linear_damping | remus100.m:194, 195 [zeta4 zeta5] | hull | published |
| `span` | `0.19` | m | lift_drag | remus100.m:132 D_auv, passed as b at :220 | hull | published |
| `planform_area` | `0.2128` | m^2 | lift_drag | remus100.m:133 S = 0.7 * L_auv * D_auv | hull | derived |
| `parasitic_drag_coefficient` | `0.05595961914206819` | 1 | lift_drag | remus100.m:143-144 CD_0 = Cd * pi * b^2 / S, b = D_auv / 2 (one geometry) | hull | derived |
| `oswald_efficiency` | `0.3` | 1 | lift_drag | coeffLiftDrag.m:54 e = 0.3 | hull | published |
| `length` | `1.6` | m | cross_flow | remus100.m:131 L_auv, passed as L at :221 | hull | published |
| `beam` | `0.19` | m | cross_flow | remus100.m:132 D_auv, passed as B at :221 | hull | published |
| `draft` | `0.19` | m | cross_flow | remus100.m:132 D_auv, passed as T at :221 | hull | published |
| `propeller_diameter` | `0.14` | m | propeller | remus100.m:148 D_prop | actuators | published |
| `max_shaft_speed` | `1525` | rpm | propeller | remus100.m:110 n_max | actuators | published |
| `thrust_deduction` | `0.1` | 1 | propeller | remus100.m:149 t_prop | actuators | published |
| `wake_fraction` | `0.05600000000000005` | 1 | propeller | remus100.m:150 Va = 0.944 * U_r, w = 1 - 0.944 | actuators | published |
| `pitch_diameter_ratio` | `1` | 1 | propeller | remus100.m:156 wageningen(0,1,0.718,3), second argument | actuators | published |
| `blade_area_ratio` | `0.718` | 1 | propeller | remus100.m:155-156 blade-area ratio 0.718 | actuators | published |
| `blade_count` | `3` | 1 | propeller | remus100.m:155-156 3 blades | actuators | published |
| `max_advance_number` | `0.6632` | 1 | propeller | remus100.m:153 Ja_max | actuators | published |
| `roll_moment_scale` | `0.1` | 1 | propeller | remus100.m:252 tau(4) = K_prop / 10 | actuators | published |
| `position` | `[0, 0, 0]` | m | propeller | remus100.m:249-252 thrust on x_b through the CO (no moment arm) | actuators | published |
| `orientation` | `[0, 0, 0]` | rad | propeller | remus100.m:249, 252 shaft along x_b (thrust in tau(1), torque in tau(4)) | actuators | published |
| `thrust_torque_coefficients` | `[0.4566, 0.07, 0.1798, 0.0312]` | 1 | propeller | remus100.m:157, 158, 160, 161 [KT_0 KQ_0 KT_max KQ_max] | actuators | published |
| `rudder_area` | `0.0133` | m^2 | fins | remus100.m:183 A_r = 2 * S_fin (S_fin :179) | actuators | published |
| `stern_plane_area` | `0.0133` | m^2 | fins | remus100.m:188 A_s = 2 * S_fin (S_fin :179) | actuators | published |
| `rudder_lift_coefficient` | `0.5` | 1/rad | fins | remus100.m:182 CL_delta_r | actuators | published |
| `stern_plane_lift_coefficient` | `0.7` | 1/rad | fins | remus100.m:187 CL_delta_s | actuators | published |
| `rudder_position` | `-0.8` | m | fins | remus100.m:184 x_r = -a, a = L_auv / 2 (one geometry) | actuators | derived |
| `stern_plane_position` | `-0.8` | m | fins | remus100.m:189 x_s = -a, a = L_auv / 2 (one geometry) | actuators | derived |
| `max_deflection` | `0.3490658503988659` | rad | fins | remus100.m:109 delta_max = deg2rad(20) | actuators | published |

Selectors (keywords of the blocks, not parameters): as the history set, plus `form = "submerged"` (linear damping) and `open_water_coefficients = "given"` (propeller). Gravity as the history set: g = 9.8217848645790813 m/s² at `mu = deg2rad(63.446827)`; `weight = buoyancy = 31.9 kg · g`.

## History: unmodified MSS `remus100.m` and its three marked copies

MSS `remus100.m` (T. I. Fossen, MIT licence) is run as written. Three copies beside the generator change only the lines they mark, and carry MSS's licence text:

| Copy | Change against `remus100.m` at `cc07579` |
|---|---|
| `remus100_hull.m` | the control inputs and the actuator forces (lines 108–115, 146–189, 233–254) replaced by an external generalized force `tau_ext` (6×1; N and N·m; BODY axes at the CO), the boundary of MSS `CRAFT/hydroVessel.m`; lines 1, 84, 90 adapted to the new argument |
| `remus100_munk.m` | the added-mass Coriolis couplings zeroed at lines 207–210 kept (`C_A = m2c(M_A, nu_r)` in full) |
| `remus100_hull_munk.m` | both changes |

A line starting with `%<removed> ` is the original line commented out; a line ending with `%<added>` is new. The generator first rebuilds the original from each copy (drop the added lines, un-comment the removed ones) and stops unless it equals MSS's `remus100.m` byte for byte; `test_marked_copies_rebuild_mss_remus100` repeats the check when `MSS_DIR` is set.

**The hull/actuator split is exact:** on the 640 cases of `remus100_derivative_mss.csv`, `remus100_hull(x, tau, …)` with the wrench `tau` that `remus100.m` itself computes (read from its workspace through a copy written to `tempdir` at run time, only its function line changed and one capture line added; equal to the unmodified function to 0) gives `remus100.m`'s `xdot` with a maximum difference of **0**.

## Files (history set)

| File | What it holds |
|---|---|
| `remus100_parameters.json` | the REMUS 100 parameter set by the class's parameter names: `value` (SI; propeller speed rpm), `unit`, `block` (the block or blocks that declare it), `block_parameter` (the name inside the block), `declared` (`new` = declared-parameter form, `old` = argument of the block's `preprocess_*`), `place` (MSS file and line at `cc07579`), `stage` (`hull` or `actuators`), `kind` (`published`); plus `selectors`, `gravity` and `source`. Every value read by the generator from the MSS workspace (`remus100.m`, `spheroid.m`, `imlay61.m`, `forceLiftDrag.m`, `coeffLiftDrag.m`, `crossFlowDrag.m`) or, for numbers that `remus100.m` writes inside an expression or a comment, parsed from its text (lines 150, 156, 252); no number typed. Numbers: `%.17g` (exact doubles). |
| `remus100_derivative_mss.csv` | 640 × 47: `case_id, x_01..x_12, ui_01..ui_03, Vc, betaVc, w_c, tau_01..tau_06, xdot_01..xdot_12, xdot_munk_01..xdot_munk_12`. (a) `xdot` of `remus100.m`; its own actuator wrench `tau`; (c) `xdot_munk` of `remus100_munk.m`. Cases: 36 structured at one state (rpm −2000, −1525, −800, −1e-3, 0, 1e-3, 800, 1525, 2000 × fins (0, 0), (20°, −20°), (0.6, −0.6) rad, (−0.2, 0.1) rad; current 0.3 m/s at 30°, `w_c` 0.1 m/s); 4 at rest (no input; fins ±0.2 rad; 1000 rpm in current 0.5 m/s at 30°, 0.1 m/s; fins and 1000 rpm); 600 seeded (`rng(20261007, 'twister')`: u ∈ [−1, 3], v, w ∈ [−1, 1] m/s, p, q, r ∈ [−0.5, 0.5] rad/s, x, y ∈ [−10, 10], z ∈ [0, 50] m, φ, ψ ∈ [−π, π], θ ∈ [−1.2, 1.2] rad; fins ∈ [−0.6, 0.6] rad (the limit is 0.349), rpm ∈ [−2000, 2000]; current `Vc` ∈ [0, 1] m/s, `betaVc` ∈ [−π, π], `w_c` ∈ [−0.2, 0.2] m/s on every second case, off on the others). |
| `remus100_hull_derivative_mss.csv` | 600 × 46: `case_id, x_01..x_12, tau_ext_01..tau_ext_06, Vc, betaVc, w_c, xdot_01..xdot_12, xdot_munk_01..xdot_munk_12`. Stage 1: `remus100_hull.m` and `remus100_hull_munk.m` on seeded states (`rng(20261008)`, the same ranges) and wrenches (forces ∈ [−50, 50] N, moments ∈ [−10, 10] N·m; current on every second case); case 1 at rest, no force, no current. |
| `remus100_mass_matrix_mss.csv` | 6 × 6 `M` of the empty call `[~,~,M] = remus100()` (`remus100.m` 83–85, 212). |
| `remus100_trajectory_<scenario>.csv` | `t`, the inputs held from that row to the next (`ui_01..ui_03` = rudder rad, stern plane rad, propeller rpm; or `tau_ext_01..tau_ext_06`), `Vc, betaVc, w_c`, `x_01..x_12` at `t`. Integration: MSS `rk4.m` at `h = 0.05` s (`SIMremus100.m` 60, 348: the Euler-angle path); inputs decided and held over 5 steps (0.25 s), one row per 5 steps, except the `hull_zoh_*` files (one row per step). Scenarios below. |
| `remus100_trajectory_tolerances.csv` | `scenario, state, offset_plus, offset_minus, reordered, tolerance`: per scenario and state, the largest difference over the rows between the reference run and (1) the same run with every derivative entry shifted by +1e-9, (2) by −1e-9 (the G1 tolerance: the largest derivative error G1 admits, applied with one sign throughout), (3) the same Runge–Kutta 4 with its weighted sum reordered (`x + h (k1/6 + k2/3 + k3/3 + k4/6)`). `tolerance` = the largest of the three; (3) is at most 7e-15, (1)–(2) between 6e-10 and 3e-7. |
| `generate_remus100_mss.m`, `remus100_hull.m`, `remus100_munk.m`, `remus100_hull_munk.m` | the generator and the three marked copies |

## Scenarios (open loop; inputs recorded in each file; both sets)

In the default set the model of each row is the marked copy with the same option (`remus100.m` → defaults, `remus100_munk.m` → `opt.munk`, `remus100_hull.m` → `opt.tau_ext`), and the Prestero rpm is 920.088018654. All start at the stated state (every entry not named is 0) with no current unless stated. `deg` = π/180.

| Scenario | Model | Initial state | Inputs | Length |
|---|---|---|---|---|
| `straight_from_rest` | `remus100.m` | at rest | 1300 rpm (the commanded speed of `SIMremus100.m` 98), fins 0 | 30 s |
| `rudder_step_20deg` | `remus100.m` | u = 1.5 m/s | 1300 rpm; rudder 20 deg from t = 5 s | 40 s |
| `rudder_step_4deg_at_1p54ms` | `remus100.m` | u = 1.54 m/s | rudder 4 deg from t = 10 s; propeller 924.958670961 rpm = the speed that holds 1.54 m/s straight in `remus100.m` (surge acceleration −4.3e-14 m/s², found with `fzero`). Prestero (2001) p. 53: rudder roughly 4 deg; p. 60 Table 8.1: u = 1.54 m/s. Steady yaw rate in MSS, last 10 s: **1.84604 deg/s** at 1.535 m/s | 60 s |
| `stern_plane_step_10deg` | `remus100.m` | u = 1.5 m/s | 1300 rpm; stern plane 10 deg (nose down) from t = 5 s | 30 s |
| `zigzag_20_20` | `remus100.m` | u = 1.5 m/s | 1300 rpm; rudder +20 deg from t = 5 s, reversed each time the heading has moved 20 deg past the initial heading (decided every 0.25 s) | 60 s |
| `dive_and_level_off` | `remus100.m` | u = 1.5 m/s | 1300 rpm; stern plane +15 deg from t = 5 s until z ≥ 3 m, then −15 deg until θ ≥ 0, then 0 (decided every 0.25 s) | 40 s |
| `rudder_step_10deg_in_current` | `remus100.m` | u = 1.5 m/s | current `Vc` 0.5 m/s, `betaVc` 30 deg, `w_c` 0.1 m/s (`SIMremus100.m` 90–92); 1300 rpm; rudder 10 deg from t = 5 s | 30 s |
| `rudder_step_15deg_at_925rpm` (default set; history runs `_mss`, `_one_density`) | the copy (`remus100.m`, the one-density model) | u = 1.54 m/s | rudder 15 deg and 925 rpm from t = 0 (the turn of the A-47 scan) | 60 s |
| `munk_straight_perturbed` | `remus100_munk.m` | u = 1.5, v = 0.01 m/s, r = 0.01 rad/s | 1300 rpm, fins 0 | 30 s |
| `munk_rudder_step_20deg` | `remus100_munk.m` | u = 1.5 m/s | 1300 rpm; rudder 20 deg from t = 5 s | 40 s |
| `hull_glide` | `remus100_hull.m` | u = 2 m/s | `tau_ext` = 0 | 30 s |
| `hull_sway_yaw_release` | `remus100_hull.m` | u = 1, v = 0.3 m/s, r = 0.2 rad/s | `tau_ext` = 0 | 20 s |
| `hull_roll_pitch_decay` | `remus100_hull.m` | u = 1 m/s, φ = 20 deg, θ = 15 deg | `tau_ext` = 0 | 20 s |
| `hull_munk_glide` | `remus100_hull_munk.m` | u = 2, v = 0.01 m/s, r = 0.01 rad/s | `tau_ext` = 0 | 20 s |
| `hull_zoh_rudder_step_20deg`, `hull_zoh_stern_plane_step_10deg`, `hull_zoh_zigzag_20_20` | `remus100_hull.m` | as the scenario of the same name | `tau_ext` = the actuator wrench of `remus100.m` at the start of each step (the inputs of that scenario), held over the step — so these differ from the full runs by the hold (end heading of the zig-zag 0.0317 vs 0.0315 rad) | as named |

## Parameter set of the unmodified MSS run (history, `remus100_parameters.json`; MSS lines at `cc07579`)

| Name | Value | Unit | Block (name inside) | Form | MSS place | Stage |
|---|---|---|---|---|---|---|
| `semi_major_axis` | `0.8076800000000001` | m | rigid_body | new | remus100.m:134 a = 1.0096 * L_auv/2 | hull |
| `semi_minor_axis` | `0.09591200000000001` | m | rigid_body | new | remus100.m:135 b = 1.0096 * D_auv/2 | hull |
| `body_density` | `1025` | kg/m^3 | rigid_body | new | spheroid.m:35 rho = 1025 (called at remus100.m:199) | hull |
| `water_density` | `1026` | kg/m^3 | rigid_body, lift_drag, propeller, fins | new, old, old, old | imlay61.m:31 rho = 1026 (remus100.m:200); equal to forceLiftDrag.m:26 and remus100.m:98 | hull |
| `roll_added_inertia_ratio` | `0.3` | 1 | rigid_body | new | remus100.m:136 r44 | hull |
| `body_center_of_gravity` | `[0, 0, 0.02]` | m | rigid_body | new | remus100.m:137 r_bG | hull |
| `center_of_buoyancy` | `[0, 0, 0]` | m | submerged_hydrostatics, submerged_linear_damping | old, old | remus100.m:138 r_bB | hull |
| `weight` | `313.31999801537864` | N | submerged_hydrostatics, submerged_linear_damping | old, old | remus100.m:214 W = m * g_mu (m = MRB(1,1); g_mu = gravity(mu), remus100.m:96-97) | hull |
| `buoyancy` | `313.31999801537864` | N | submerged_hydrostatics | old | remus100.m:214 B = W | hull |
| `time_constants` | `[20, 20, 1]` | s | submerged_linear_damping | old | remus100.m:192, 193, 196 [T1 T2 T6] (Dmtrx.m call at :217) | hull |
| `damping_ratios` | `[0.3, 0.8]` | 1 | submerged_linear_damping | old | remus100.m:194, 195 [zeta4 zeta5] | hull |
| `span` | `0.19` | m | lift_drag | old | remus100.m:132 D_auv, passed as b at :220 | hull |
| `planform_area` | `0.2128` | m^2 | lift_drag | old | remus100.m:133 S = 0.7 * L_auv * D_auv | hull |
| `parasitic_drag_coefficient` | `0.05703920106809604` | 1 | lift_drag | old | remus100.m:143-144 CD_0 = Cd * pi * b^2 / S | hull |
| `oswald_efficiency` | `0.3` | 1 | lift_drag | old | coeffLiftDrag.m:54 e = 0.3 | hull |
| `length` | `1.6` | m | cross_flow | old | remus100.m:131 L_auv, passed as L at :221 | hull |
| `beam` | `0.19` | m | cross_flow | old | remus100.m:132 D_auv, passed as B at :221 | hull |
| `draft` | `0.19` | m | cross_flow | old | remus100.m:132 D_auv, passed as T at :221 | hull |
| `cross_flow_water_density` | `1025` | kg/m^3 | cross_flow (`water_density`) | old | crossFlowDrag.m:36 rho = 1025 (remus100.m:221) | hull |
| `propeller_diameter` | `0.14` | m | propeller | old | remus100.m:148 D_prop | actuators |
| `max_shaft_speed` | `1525` | rpm | propeller | old | remus100.m:110 n_max | actuators |
| `thrust_deduction` | `0.1` | 1 | propeller | old | remus100.m:149 t_prop | actuators |
| `wake_fraction` | `0.05600000000000005` | 1 | propeller | old | remus100.m:150 Va = 0.944 * U_r, w = 1 - 0.944 | actuators |
| `pitch_diameter_ratio` | `1` | 1 | propeller | old | remus100.m:156 wageningen(0,1,0.718,3), second argument | actuators |
| `blade_area_ratio` | `0.718` | 1 | propeller | old | remus100.m:155-156 blade-area ratio 0.718 | actuators |
| `blade_count` | `3` | 1 | propeller | old | remus100.m:155-156 3 blades | actuators |
| `max_advance_number` | `0.6632` | 1 | propeller | old | remus100.m:153 Ja_max | actuators |
| `roll_moment_scale` | `0.1` | 1 | propeller | old | remus100.m:252 tau(4) = K_prop / 10 | actuators |
| `position` | `[0, 0, 0]` | m | propeller | old | remus100.m:249-252 thrust on x_b through the CO (no moment arm) | actuators |
| `orientation` | `[0, 0, 0]` | rad | propeller | old | remus100.m:249, 252 shaft along x_b (thrust in tau(1), torque in tau(4)) | actuators |
| `thrust_torque_coefficients` | `[0.4566, 0.07, 0.1798, 0.0312]` | 1 | propeller | old | remus100.m:157, 158, 160, 161 [KT_0 KQ_0 KT_max KQ_max] | actuators |
| `rudder_area` | `0.0133` | m^2 | fins | old | remus100.m:183 A_r = 2 * S_fin (S_fin :179) | actuators |
| `stern_plane_area` | `0.0133` | m^2 | fins | old | remus100.m:188 A_s = 2 * S_fin (S_fin :179) | actuators |
| `rudder_lift_coefficient` | `0.5` | 1/rad | fins | old | remus100.m:182 CL_delta_r | actuators |
| `stern_plane_lift_coefficient` | `0.7` | 1/rad | fins | old | remus100.m:187 CL_delta_s | actuators |
| `rudder_position` | `-0.8076800000000001` | m | fins | old | remus100.m:184 x_r = -a | actuators |
| `stern_plane_position` | `-0.8076800000000001` | m | fins | old | remus100.m:189 x_s = -a | actuators |
| `max_deflection` | `0.3490658503988659` | rad | fins | old | remus100.m:109 delta_max = deg2rad(20) | actuators |

Selectors (keywords of the blocks, not parameters): `mass_properties = "spheroid"`, `coriolis = "co"` (`spheroid.m` 46–52), cross-flow `drag_model = "cylinder"` (`remus100.m` 221) with `strip_grid = "midpoint"` (`crossFlowDrag.m` 56), `sway_damping_fade = false` and `smooth_speed = false` (MSS exactly, `remus100.m` 127, 218; with `true` the regularisation `smooth_speed_epsilon` would be a parameter: regularisation, not physics, not part of this set), propeller `open_water_model = "linearized"` (`remus100.m` 163–177) with `clip_advance_ratio = false`, fins `convention = "starboard_down_positive"` (`remus100.m` 242, 245). `stabilize_added_mass_coriolis` is chosen per test (`true` = MSS, `remus100.m` 207–210).

Gravity: `weight = MRB(1,1) · gravity(mu)`, `mu = deg2rad(63.446827)` (`remus100.m` 96–97, 214; `INS/functions/gravity.m` 11–12): g = 9.8217848645790813 m/s², m = weight / g = 31.900515266356955 kg (`spheroid.m` 35–36 with the scaled semi-axes). The class takes `weight` and `buoyancy` as values; nothing in the class computes them from the mass.

Density: MSS uses 1025 kg/m³ for the body (`spheroid.m` 35) and the cross-flow (`crossFlowDrag.m` 36), 1026 kg/m³ for the added mass (`imlay61.m` 31), the hull lift and drag (`forceLiftDrag.m` 26), the propeller and the fins (`remus100.m` 98). The history set keeps the cross-flow value apart as `cross_flow_water_density`; the default set has one `water_density`.

## Numbers and conventions

- `%.17g` in every CSV and in the JSON: each number reads back to the double MATLAB held.
- Units SI; angles rad; propeller speed rpm; NED positions, Euler ZYX angles, BODY velocities (`remus100.m` 7–24); current as `remus100.m` 40–43: `v_c = [Vc cos(betaVc − psi), Vc sin(betaVc − psi), w_c]` in BODY.

## sha256

History set:


- `generate_remus100_mss.m` `3b06c0952a1e780ce0f2534493a46fcc8109c687c2a6e594031a7193ca654b56`
- `remus100_derivative_mss.csv` `2eedf8ea294abedaa2f376f4d88d109b9cfc57c0d97b86814660e943057b4379`
- `remus100_hull_derivative_mss.csv` `446982c6e7e8c837d94bd345a548fddc86fcf60f53f6b29873c81ae8bc80d641`
- `remus100_hull.m` `33a411b89c4dc6f9bca8d0a6edcbf7d13f21cc4b0c82a083c0276427dd5eafe0`
- `remus100_hull_munk.m` `726ee6fee77ee8958862306b4229d59a44b1b5bc2485ee881a69e717deaa184e`
- `remus100_mass_matrix_mss.csv` `1789d81f5751c3493e3a6f08f30c66f1928b7df33c4dccf9997dbc7542c85ba0`
- `remus100_munk.m` `e6cb656d34c5022481a4e5cd73bc20e0fc1a3fb5786475cfaeb9d71acbbe596f`
- `remus100_parameters.json` `ba10fa9821e36b88fd2c8be819ecfdc129373c9030cc8ec8f6e4f771ffdab2a5`
- `remus100_trajectory_dive_and_level_off.csv` `3af0e7b7e93c8d5c2f08882c138da61f15584386ab5fdc198ea78799ba0e809b`
- `remus100_trajectory_hull_glide.csv` `58d3abfcfbbdba9f7a257bc9c4855d6470ae1e89a48e534cef7b221dc0af9c67`
- `remus100_trajectory_hull_munk_glide.csv` `549b3db15cae1178469e1af8a5b2a7c8c1398fc32598f276e7b77c7df8488fb9`
- `remus100_trajectory_hull_roll_pitch_decay.csv` `1ef01efa85d844aa6df000a66ea171743c92c371c0461662fba74ba2b964e8db`
- `remus100_trajectory_hull_sway_yaw_release.csv` `b747d7be6279410c2a9b942d07d155119fb7c05a513d9f44d4e05cc651c29131`
- `remus100_trajectory_hull_zoh_rudder_step_20deg.csv` `b3e48fb25864dcdc09e9308e7acf6aa76eec3a9ccd59314213000332fb387cb5`
- `remus100_trajectory_hull_zoh_stern_plane_step_10deg.csv` `f88225927be7637fede0dc4e632a7cb1a56c05b96b621e9b6ce873015166960d`
- `remus100_trajectory_hull_zoh_zigzag_20_20.csv` `fa1595d5fc944f726426b06c532ba2a91e3d6efd268b9b3fba39287cdfb4f9bf`
- `remus100_trajectory_munk_rudder_step_20deg.csv` `d970de52138ef038ed122a328d0bffda3b244f3d0c8100d8230ec1517efedd9b`
- `remus100_trajectory_munk_straight_perturbed.csv` `a4672924e4ea4e574265dc7edd536f48635428923497cb8c323f9223d70984f2`
- `remus100_trajectory_rudder_step_10deg_in_current.csv` `52715d057ec6c31254e7429c1c23448d8b38e455672ed0998058b2d6f9febbdc`
- `remus100_trajectory_rudder_step_20deg.csv` `0e6ed69548e3475bdd96a4f321fa6184a688b1fc61680822916feaa311a328af`
- `remus100_trajectory_rudder_step_4deg_at_1p54ms.csv` `00c6cab65f15f0f45c1341d35499d1380a5b593d0040ad98d9277d75d52d1385`
- `remus100_trajectory_stern_plane_step_10deg.csv` `ad48b7f80e91bb98af0e069e630c2e83830396e98c93f2f951de2575a96a3389`
- `remus100_trajectory_straight_from_rest.csv` `4d4325ae0e27128c5bdcd66a25bf8b141e3ed6cc548bb322d8c1ac3f954d8651`
- `remus100_trajectory_tolerances.csv` `89efdcc934395073f795fa65c13aadc881085ffee79dee5236f128377b6334c8`
- `remus100_trajectory_zigzag_20_20.csv` `b4a85a3374a9c3b8ce47ad5ac6a80bc6b50d79a5079ac9753e4fe1da862a2fc0`

Default set:

- `crossFlowDrag_consistent.m` `6432d2e2609a8f19b5755e67a489166e52e67480fe9660e6d9fb6ef917a9a5bd`
- `forceLiftDrag_consistent.m` `2743899a6c8c3a396a344ffee3e4751011af4c089cfd53b5293e334aba3c31ad`
- `generate_remus100_consistent.m` `c8b0b86a4a2ec0ef495da9ef6ee70d52e904e9fe7cb5da0c98c1733fa4c00dc4`
- `imlay61_consistent.m` `b57b3cf10b4fc8ea2efc58bfb2845398854f3149ec26948062d7a8bfaef7162b`
- `remus100_consistent.m` `82023c064f1b07f77a3cb7afdb3f869dd3a8d5d069e9c4d41415d8224b3b50b7`
- `remus100_derivative_consistent.csv` `7c1524b5acfcc441561d67d5f2feb0fbf01ec4171522ec3136932ac736c370a0`
- `remus100_hull_derivative_consistent.csv` `a16e7e633983dbf23d81758d835f0b9ac47acb2c66240f09675cdb62b81112e2`
- `remus100_mass_matrix_consistent.csv` `d0a6402a62fb399b22152fdb4c38246a6f5a5f702bb9b03c2785486d292efff6`
- `remus100_parameters_consistent.json` `6e85a98ca9efff777c233ec47aff197b12e4750983d3a17af4babf0e041b0cff`
- `remus100_trajectory_dive_and_level_off_consistent.csv` `3e3e8dac0bb525e0f21c85fd8cdc6994b4cf53e6c5ddffbf674867a16d8ca368`
- `remus100_trajectory_hull_glide_consistent.csv` `5c52be0765ead46d31bc77cca98c3bf6a385e06514cf24fe715f10d5872dbba8`
- `remus100_trajectory_hull_munk_glide_consistent.csv` `9aa07e2e16e34ae21df30b243367c630aab6756c4e3e24f7f119479a7a294356`
- `remus100_trajectory_hull_roll_pitch_decay_consistent.csv` `46fc54c9395e44b6986f9b0fd0865b36d93555f73e26f2a2789a6dbded91d6be`
- `remus100_trajectory_hull_sway_yaw_release_consistent.csv` `c83c969d790f9d068f41ca6493f25c01dc80094d0083a04c6ad68e10f65f690d`
- `remus100_trajectory_hull_zoh_rudder_step_20deg_consistent.csv` `2c39d07750ff0c5362a682bc814ef286c5c956a0b8935ad11ae13aef9a9de195`
- `remus100_trajectory_hull_zoh_stern_plane_step_10deg_consistent.csv` `e14ff763f033c83237d129f05feb30f7df5bfd4b027f7a3754d9c6032c2c2fbf`
- `remus100_trajectory_hull_zoh_zigzag_20_20_consistent.csv` `b733fc0e3b45653452f7a0e3d594e3c2ac2ed12db68018599846eb60c697033d`
- `remus100_trajectory_munk_rudder_step_20deg_consistent.csv` `b04f12a04bf8a95d5eac3d3fd8004f6fc16a17d59608843e37c7b0744c5622f8`
- `remus100_trajectory_munk_straight_perturbed_consistent.csv` `b0dfd2a9ba165ce352b8e05947e70d9c422772e261678956d0545cbf20b70a56`
- `remus100_trajectory_rudder_step_10deg_in_current_consistent.csv` `003e93cc2783692c716ac537162194b0baeaea8c7f852a1107638dc503976127`
- `remus100_trajectory_rudder_step_15deg_at_925rpm_consistent.csv` `c6278ef930685ace8cd6c0d9bec0a0c0d0d9bde65845c2d7d3e466cfbddcd2e8`
- `remus100_trajectory_rudder_step_15deg_at_925rpm_mss.csv` `7aa709ce34a62e226ddcd97dbafd3f17357418b6cc9737e75e1abf8f1115d111`
- `remus100_trajectory_rudder_step_15deg_at_925rpm_one_density.csv` `fec64606279ef986828bfdc5a9ea735324fea7d2f406f4f3734471f5996c3172`
- `remus100_trajectory_rudder_step_20deg_consistent.csv` `9a4ce160cd8ae25459e9f26ed36f0b62b89bae9697376396e5b41ff4fa32d360`
- `remus100_trajectory_rudder_step_4deg_at_1p54ms_consistent.csv` `5f8b0581abb903eefb25dbdd788539e6f43989aa67d7f00aa2f57ca950e8c213`
- `remus100_trajectory_stern_plane_step_10deg_consistent.csv` `ad6589d2e65df9c5cb34ec99e831c495e91be31961284b1a2f97c998616201ab`
- `remus100_trajectory_straight_from_rest_consistent.csv` `bbd8c92c90bdd89b3b5d532b8ba070c41e2e7916f54905bb8e1cb743d624e5f9`
- `remus100_trajectory_tolerances_consistent.csv` `30bf56f25e57fd4672e94bb9202deb4ab86e00b22f01ce93a531c3cc84492d16`
- `remus100_trajectory_zigzag_20_20_consistent.csv` `292fbb1c1e73f0c06a3cecc4b31ec396d19a1c66448f56b49313c2e0d8d25f6b`
- `spheroid_consistent.m` `853516efd11df409425d7dda9f1ad834009b808b54a8d6b5f3c191308680dea1`

## Readers

`tests/vehicles/test_torpedo_class.py`.

Author: Enio Krizman · Date: 2026-10-07


