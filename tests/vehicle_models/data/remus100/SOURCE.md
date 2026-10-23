# Frozen reference — torpedo-AUV class on the REMUS 100 parameter set (MSS `remus100.m`, one value per quantity)

**Frozen** files: a changed byte in a CSV or in the parameter file needs a decision record; a new MSS revision that changes an output adds a file `*_mss_<revision>.csv` beside the old one. The generator lives in this folder, finds MSS only through `MSS_DIR` (MATLAB through `MATLAB_BIN`) and writes under `OUT_DIR` (default `<tempdir>/more_mss_references/vehicles/remus100/`), never into this folder; a re-run is checked with `cmp` against the file of the same name here.

**Two reference sets.** The **default** (the class is gated against it) is the one-value set, `*_consistent.*`, made by `generate_remus100_consistent.m` from the marked copy `remus100_consistent.m`: one water density, one geometry, the mass given (sections below). The **history** set is the output of unmodified MSS `remus100.m`, made by `generate_remus100_mss.m` (`*_mss.csv`, the unsuffixed trajectory files, `remus100_parameters.json`); it is kept, never overwritten, and one test shows by how much it differs. Why: within one model MSS uses two water densities and two geometries (registered as MSS faults; the owner, 2026-10-07: *"we dont match mss if doesnt make sense"*); our model holds one value per physical quantity.

History set: made 2026-10-07 by MATLAB R2026a Update 4 running MSS at `99fdf4dbdfe0b3f124e26a689c6aabe0beb5654f` (`CRAFT/AUV/models/remus100.m`, `CRAFT/AUV/SIMremus100.m`, `LIBRARY/numericalMethods/rk4.m`). Re-run the same day from this folder: all 21 outputs byte-identical.

Default set: made 2026-10-07 by the same MATLAB and MSS revision; a second run the same day wrote all 24 outputs byte-identical.

Run: `MSS_DIR=<MSS checkout> "$MATLAB_BIN" -batch "run('<this folder>/generate_remus100_mss.m')"`.

## The one-value reference (the default)

**The model.** `remus100_consistent.m` is a marked copy of MSS `remus100.m` at `99fdf4d` (T. I. Fossen, MIT licence, carried in the copy). Every changed line ends with `%<added>` and says `one value:`; a line starting with `%<removed> ` is the original commented out; the generator rebuilds the original from the copy and stops unless it equals MSS byte for byte (so does `test_marked_copies_rebuild_their_mss_files` when `MSS_DIR` is set). Four marked copies of MSS library functions take the quantity they had fixed inside as an input: `spheroid_consistent.m` (the mass; MSS: 1025 kg/m³ × volume, `spheroid.m` 35–36), `imlay61_consistent.m` (ρ; MSS 1026, `imlay61.m` 31), `forceLiftDrag_consistent.m` (ρ; MSS 1026, `forceLiftDrag.m` 26), `crossFlowDrag_consistent.m` (ρ; MSS 1025, `crossFlowDrag.m` 36).

| Quantity | One value | MSS `remus100.m` uses |
|---|---|---|
| water density | **1025 kg/m³** (`remus100.m` 98–99, `mssConstants.m` 12) in the added mass, hull lift and drag, cross-flow, propeller and fins | the same: 1025 everywhere since `23c307e` (F-MSS-19 fixed upstream) |
| length, diameter | **L = 1.6 m, D = 0.19 m** (132–133) everywhere: a = L/2, b = D/2 for the mass, the added mass, `CD_0` (144–145) and the fin positions x_r = x_s = −a (185, 190), as for S (134), lift and drag (221) and cross-flow (222) | 1.0096·L/2, 1.0096·D/2 for mass, added mass, `CD_0`, fins (135–136); L, D for S, lift and drag, cross-flow |
| mass | **31.9 kg, given** (line 3: "the mass of the vehicle is 31.9 kg"; the generator reads line 3 and stops unless the copy holds the same number); inertia of a homogeneous spheroid of that mass (`spheroid.m` 39–42) | 1025 kg/m³ × the scaled spheroid volume = 31.9005 kg |
| weight, buoyancy | **W = m g_mu, B = W** as `remus100.m` 214 states it (neutrally buoyant, the REMUS assumption); not B = ρ g ∇, which would be 31.03 kg of water against 31.9 kg of vehicle (0.87 kg heavy) | the same line |

The propeller parameters are named `propeller_diameter` and `max_shaft_speed` after the part declarations; their values and MSS places are unchanged. The parameter set gives the mass through `body_density` = mass / (4/3 π a b²) = 1054.79 kg/m³, the input the rigid-body block's `"spheroid"` form declares (it has no mass input).

**Options of the copy** (6th argument, a struct; the defaults are the one-value model): `opt.munk` keeps the added-mass Coriolis couplings of lines 207–210 (the `xdot_munk` columns); `opt.tau_ext` replaces the actuator wrench by an external one (the bare hull, stage 1; `ui = 0`); `opt.rho`, `opt.mass` change the one density and the mass (the invariance check); `opt.geometry = 'mss'` gives MSS's scaled semi-axes and the mass of the displaced water at the one ρ — MSS with one density only, used for history (`xdot_one_density`, `*_one_density.csv`). The fourth output is the actuator wrench.

**Density invariance (measured, MATLAB).** With one density reaching every call and the mass scaling with it, every force, mass and damping term is proportional to ρ (Dmtrx's roll and pitch frequencies √(W Δz / M) do not depend on it; the cylinder drag coefficient depends on viscosity, not ρ), so ν̇ cannot depend on ρ. On the 640 derivative cases, `max |Δν̇| / max |ν̇|` with (ρ, m) scaled together: **0 exactly** for the factors 0.5 and 2 (powers of two scale every number exactly), **1.2e-14** for 1000/1025 (rounding of ρ and m); with MSS geometry and the mass from ρ, ρ = 1000 against 1025: 1.5e-14; the bare hull with (ρ, m, `tau_ext`) scaled by 1000/1025 on the 600 hull cases: 1.1e-14. Unmodified MSS is not invariant (A-47: 1.6e-3).

**The split is exact:** the hull copy (`opt.tau_ext`) fed with the full copy's own actuator wrench gives the full copy's `xdot` with difference **0** on the 640 cases.

**How far the history set is from the default (measured; `test_history_unmodified_mss_differs_by_one_density_and_one_geometry`).**

| Change | Derivative, 640 cases: max over cases of `max|Δν̇|/max|ν̇|` (and `‖Δν̇‖/‖ν̇‖`) | Turn, 15° rudder, 925 rpm, from 1.54 m/s, 60 s (`rudder_step_15deg_at_925rpm`) |
|---|---|---|
| density alone: MSS → one density (MSS geometry) | 0 (0): one density in MSS since `23c307e` | heading +0.000°, turn rate 0.087188 → 0.087188 rad/s, position 0.0 m |
| geometry and given mass: one density → one value | 5.60e-2 (4.68e-2) | heading **+2.797°**, turn rate → 0.087996 rad/s, position 0.90 m, speed 1.5078 → 1.5153 m/s |
| total: MSS → one value | 5.60e-2 (4.68e-2) | heading +2.797°, position 0.90 m |

The density row was −0.179° / 0.1 % at `cc07579` (A-47 scan, two water densities in one MSS model, F-MSS-19); at `99fdf4d` MSS holds one density and the row is zero. The geometry row is our choice (the hull's own L and D everywhere, mass given) and is +2.8°, not the −4.06° of A-47's other variant (the *scaled* geometry everywhere).

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

| Name | Value | Unit | Block (all in the declared-parameter form) | Place (MSS lines at `99fdf4d`) | Stage | Kind |
|---|---|---|---|---|---|---|
| `semi_major_axis` | `0.8` | m | rigid_body | remus100.m:132 L_auv / 2 (one geometry; MSS scales by 1.0096 at :135) | hull | published |
| `semi_minor_axis` | `0.095` | m | rigid_body | remus100.m:133 D_auv / 2 (one geometry; MSS scales by 1.0096 at :136) | hull | published |
| `body_density` | `1054.7872613500267` | kg/m^3 | rigid_body | mass 31.9 kg (remus100.m:3, given) / (4/3 pi a b^2); MSS: mssConstants.m:12 = 1025 at spheroid.m:35-36 | hull | derived |
| `water_density` | `1025` | kg/m^3 | rigid_body, lift_drag, cross_flow, propeller, fins | remus100.m:98-99 rho = constant.rho_water, the one value (mssConstants.m:12 = 1025; read the same way at imlay61.m:31-32, forceLiftDrag.m:26-27 and crossFlowDrag.m:36-37) | hull | published |
| `roll_added_inertia_ratio` | `0.3` | 1 | rigid_body | remus100.m:137 r44 | hull | published |
| `body_center_of_gravity` | `[0, 0, 0.02]` | m | rigid_body | remus100.m:138 r_bG | hull | published |
| `center_of_buoyancy` | `[0, 0, 0]` | m | submerged_hydrostatics, submerged_linear_damping | remus100.m:139 r_bB | hull | published |
| `weight` | `313.31493718007266` | N | submerged_hydrostatics, submerged_linear_damping | remus100.m:215 W = m * g_mu, m = 31.9 kg given (remus100.m:3), g_mu = gravity(mu) (remus100.m:96-97) | hull | derived |
| `buoyancy` | `313.31493718007266` | N | submerged_hydrostatics | remus100.m:215 B = W (neutral, as REMUS states) | hull | derived |
| `time_constants` | `[20, 20, 1]` | s | submerged_linear_damping | remus100.m:193, 194, 197 [T1 T2 T6] (Dmtrx.m call at :218) | hull | published |
| `damping_ratios` | `[0.3, 0.8]` | 1 | submerged_linear_damping | remus100.m:195, 196 [zeta4 zeta5] | hull | published |
| `span` | `0.19` | m | lift_drag | remus100.m:133 D_auv, passed as b at :221 | hull | published |
| `planform_area` | `0.2128` | m^2 | lift_drag | remus100.m:134 S = 0.7 * L_auv * D_auv | hull | derived |
| `parasitic_drag_coefficient` | `0.05595961914206819` | 1 | lift_drag | remus100.m:144-145 CD_0 = Cd * pi * b^2 / S, b = D_auv / 2 (one geometry) | hull | derived |
| `oswald_efficiency` | `0.3` | 1 | lift_drag | coeffLiftDrag.m:54 e = 0.3 | hull | published |
| `length` | `1.6` | m | cross_flow | remus100.m:132 L_auv, passed as L at :222 | hull | published |
| `beam` | `0.19` | m | cross_flow | remus100.m:133 D_auv, passed as B at :222 | hull | published |
| `draft` | `0.19` | m | cross_flow | remus100.m:133 D_auv, passed as T at :222 | hull | published |
| `propeller_diameter` | `0.14` | m | propeller | remus100.m:149 D_prop | actuators | published |
| `max_shaft_speed` | `1525` | rpm | propeller | remus100.m:111 n_max | actuators | published |
| `thrust_deduction` | `0.1` | 1 | propeller | remus100.m:150 t_prop | actuators | published |
| `wake_fraction` | `0.05600000000000005` | 1 | propeller | remus100.m:151 Va = 0.944 * U_r, w = 1 - 0.944 | actuators | published |
| `pitch_diameter_ratio` | `1` | 1 | propeller | remus100.m:157 wageningen(0,1,0.718,3), second argument | actuators | published |
| `blade_area_ratio` | `0.718` | 1 | propeller | remus100.m:156-157 blade-area ratio 0.718 | actuators | published |
| `blade_count` | `3` | 1 | propeller | remus100.m:156-157 3 blades | actuators | published |
| `max_advance_number` | `0.6632` | 1 | propeller | remus100.m:154 Ja_max | actuators | published |
| `roll_moment_scale` | `0.1` | 1 | propeller | remus100.m:253 tau(4) = K_prop / 10 | actuators | published |
| `position` | `[0, 0, 0]` | m | propeller | remus100.m:250-253 thrust on x_b through the CO (no moment arm) | actuators | published |
| `orientation` | `[0, 0, 0]` | rad | propeller | remus100.m:250, 253 shaft along x_b (thrust in tau(1), torque in tau(4)) | actuators | published |
| `thrust_torque_coefficients` | `[0.4566, 0.07, 0.1798, 0.0312]` | 1 | propeller | remus100.m:158, 159, 161, 162 [KT_0 KQ_0 KT_max KQ_max] | actuators | published |
| `rudder_area` | `0.0133` | m^2 | fins | remus100.m:184 A_r = 2 * S_fin (S_fin :180) | actuators | published |
| `stern_plane_area` | `0.0133` | m^2 | fins | remus100.m:189 A_s = 2 * S_fin (S_fin :180) | actuators | published |
| `rudder_lift_coefficient` | `0.5` | 1/rad | fins | remus100.m:183 CL_delta_r | actuators | published |
| `stern_plane_lift_coefficient` | `0.7` | 1/rad | fins | remus100.m:188 CL_delta_s | actuators | published |
| `rudder_position` | `-0.8` | m | fins | remus100.m:185 x_r = -a, a = L_auv / 2 (one geometry) | actuators | derived |
| `stern_plane_position` | `-0.8` | m | fins | remus100.m:190 x_s = -a, a = L_auv / 2 (one geometry) | actuators | derived |
| `max_deflection` | `0.3490658503988659` | rad | fins | remus100.m:110 delta_max = deg2rad(20) | actuators | published |

Selectors (keywords of the blocks, not parameters): as the history set, plus `form = "submerged"` (linear damping) and `open_water_coefficients = "given"` (propeller). Gravity as the history set: g = 9.8217848645790813 m/s² at `mu = deg2rad(63.446827)`; `weight = buoyancy = 31.9 kg · g`.

## History: unmodified MSS `remus100.m` and its three marked copies

MSS `remus100.m` (T. I. Fossen, MIT licence) is run as written. Three copies beside the generator change only the lines they mark, and carry MSS's licence text:

| Copy | Change against `remus100.m` at `99fdf4d` |
|---|---|
| `remus100_hull.m` | the control inputs and the actuator forces (lines 108–115, 146–189, 233–254) replaced by an external generalized force `tau_ext` (6×1; N and N·m; BODY axes at the CO), the boundary of MSS `CRAFT/hydroVessel.m`; lines 1, 84, 90 adapted to the new argument |
| `remus100_munk.m` | the added-mass Coriolis couplings zeroed at lines 207–210 kept (`C_A = m2c(M_A, nu_r)` in full) |
| `remus100_hull_munk.m` | both changes |

A line starting with `%<removed> ` is the original line commented out; a line ending with `%<added>` is new. The generator first rebuilds the original from each copy (drop the added lines, un-comment the removed ones) and stops unless it equals MSS's `remus100.m` byte for byte; `test_marked_copies_rebuild_mss_remus100` repeats the check when `MSS_DIR` is set.

**The hull/actuator split is exact:** on the 640 cases of `remus100_derivative_mss.csv`, `remus100_hull(x, tau, …)` with the wrench `tau` that `remus100.m` itself computes (read from its workspace through a copy written to `tempdir` at run time, only its function line changed and one capture line added; equal to the unmodified function to 0) gives `remus100.m`'s `xdot` with a maximum difference of **0**.

## Files (history set)

| File | What it holds |
|---|---|
| `remus100_parameters.json` | the REMUS 100 parameter set by the class's parameter names: `value` (SI; propeller speed rpm), `unit`, `block` (the block or blocks that declare it), `block_parameter` (the name inside the block), `declared` (`new` = declared-parameter form, `old` = argument of the block's `preprocess_*`), `place` (MSS file and line at `99fdf4d`), `stage` (`hull` or `actuators`), `kind` (`published`); plus `selectors`, `gravity` and `source`. Every value read by the generator from the MSS workspace (`remus100.m`, `spheroid.m`, `imlay61.m`, `forceLiftDrag.m`, `coeffLiftDrag.m`, `crossFlowDrag.m`) or, for numbers that `remus100.m` writes inside an expression or a comment, parsed from its text (lines 150, 156, 252); no number typed. Numbers: `%.17g` (exact doubles). |
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

## Parameter set of the unmodified MSS run (history, `remus100_parameters.json`; MSS lines at `99fdf4d`)

This file is kept unmodified and hash-pinned below. Its `block_parameter` names for the propeller (`diameter`, `max_speed`) are the keyword names of the block's earlier `preprocess_*` form (`declared` = `old`), as `generate_remus100_mss.m` writes them; the block now declares `propeller_diameter` and `max_shaft_speed` (rpm), and the live set `remus100_parameters_consistent.json` uses those names. No test reads this file against the block's declarations.

| Name | Value | Unit | Block (name inside) | Form | MSS place | Stage |
|---|---|---|---|---|---|---|
| `semi_major_axis` | `0.8076800000000001` | m | rigid_body | new | remus100.m:135 a = 1.0096 * L_auv/2 | hull |
| `semi_minor_axis` | `0.09591200000000001` | m | rigid_body | new | remus100.m:136 b = 1.0096 * D_auv/2 | hull |
| `body_density` | `1025` | kg/m^3 | rigid_body | new | spheroid.m:35-36 rho = constant.rho_water, mssConstants.m:12 = 1025 (called at remus100.m:200) | hull |
| `water_density` | `1025` | kg/m^3 | rigid_body, lift_drag, propeller, fins | new, old, old, old | imlay61.m:31-32 rho = constant.rho_water, mssConstants.m:12 = 1025 (remus100.m:201); equal to forceLiftDrag.m:26-27 and remus100.m:98-99 | hull |
| `roll_added_inertia_ratio` | `0.3` | 1 | rigid_body | new | remus100.m:137 r44 | hull |
| `body_center_of_gravity` | `[0, 0, 0.02]` | m | rigid_body | new | remus100.m:138 r_bG | hull |
| `center_of_buoyancy` | `[0, 0, 0]` | m | submerged_hydrostatics, submerged_linear_damping | old, old | remus100.m:139 r_bB | hull |
| `weight` | `313.31999801537864` | N | submerged_hydrostatics, submerged_linear_damping | old, old | remus100.m:215 W = m * g_mu (m = MRB(1,1); g_mu = gravity(mu), remus100.m:96-97) | hull |
| `buoyancy` | `313.31999801537864` | N | submerged_hydrostatics | old | remus100.m:215 B = W | hull |
| `time_constants` | `[20, 20, 1]` | s | submerged_linear_damping | old | remus100.m:193, 194, 197 [T1 T2 T6] (Dmtrx.m call at :218) | hull |
| `damping_ratios` | `[0.3, 0.8]` | 1 | submerged_linear_damping | old | remus100.m:195, 196 [zeta4 zeta5] | hull |
| `span` | `0.19` | m | lift_drag | old | remus100.m:133 D_auv, passed as b at :221 | hull |
| `planform_area` | `0.2128` | m^2 | lift_drag | old | remus100.m:134 S = 0.7 * L_auv * D_auv | hull |
| `parasitic_drag_coefficient` | `0.05703920106809604` | 1 | lift_drag | old | remus100.m:144-145 CD_0 = Cd * pi * b^2 / S | hull |
| `oswald_efficiency` | `0.3` | 1 | lift_drag | old | coeffLiftDrag.m:54 e = 0.3 | hull |
| `length` | `1.6` | m | cross_flow | old | remus100.m:132 L_auv, passed as L at :222 | hull |
| `beam` | `0.19` | m | cross_flow | old | remus100.m:133 D_auv, passed as B at :222 | hull |
| `draft` | `0.19` | m | cross_flow | old | remus100.m:133 D_auv, passed as T at :222 | hull |
| `cross_flow_water_density` | `1025` | kg/m^3 | cross_flow (`water_density`) | old | crossFlowDrag.m:36-37 rho = constant.rho_water, mssConstants.m:12 = 1025 (remus100.m:222) | hull |
| `propeller_diameter` | `0.14` | m | propeller | old | remus100.m:149 D_prop | actuators |
| `max_shaft_speed` | `1525` | rpm | propeller | old | remus100.m:111 n_max | actuators |
| `thrust_deduction` | `0.1` | 1 | propeller | old | remus100.m:150 t_prop | actuators |
| `wake_fraction` | `0.05600000000000005` | 1 | propeller | old | remus100.m:151 Va = 0.944 * U_r, w = 1 - 0.944 | actuators |
| `pitch_diameter_ratio` | `1` | 1 | propeller | old | remus100.m:157 wageningen(0,1,0.718,3), second argument | actuators |
| `blade_area_ratio` | `0.718` | 1 | propeller | old | remus100.m:156-157 blade-area ratio 0.718 | actuators |
| `blade_count` | `3` | 1 | propeller | old | remus100.m:156-157 3 blades | actuators |
| `max_advance_number` | `0.6632` | 1 | propeller | old | remus100.m:154 Ja_max | actuators |
| `roll_moment_scale` | `0.1` | 1 | propeller | old | remus100.m:253 tau(4) = K_prop / 10 | actuators |
| `position` | `[0, 0, 0]` | m | propeller | old | remus100.m:250-253 thrust on x_b through the CO (no moment arm) | actuators |
| `orientation` | `[0, 0, 0]` | rad | propeller | old | remus100.m:250, 253 shaft along x_b (thrust in tau(1), torque in tau(4)) | actuators |
| `thrust_torque_coefficients` | `[0.4566, 0.07, 0.1798, 0.0312]` | 1 | propeller | old | remus100.m:158, 159, 161, 162 [KT_0 KQ_0 KT_max KQ_max] | actuators |
| `rudder_area` | `0.0133` | m^2 | fins | old | remus100.m:184 A_r = 2 * S_fin (S_fin :180) | actuators |
| `stern_plane_area` | `0.0133` | m^2 | fins | old | remus100.m:189 A_s = 2 * S_fin (S_fin :180) | actuators |
| `rudder_lift_coefficient` | `0.5` | 1/rad | fins | old | remus100.m:183 CL_delta_r | actuators |
| `stern_plane_lift_coefficient` | `0.7` | 1/rad | fins | old | remus100.m:188 CL_delta_s | actuators |
| `rudder_position` | `-0.8076800000000001` | m | fins | old | remus100.m:185 x_r = -a | actuators |
| `stern_plane_position` | `-0.8076800000000001` | m | fins | old | remus100.m:190 x_s = -a | actuators |
| `max_deflection` | `0.3490658503988659` | rad | fins | old | remus100.m:110 delta_max = deg2rad(20) | actuators |

Selectors (keywords of the blocks, not parameters): `mass_properties = "spheroid"`, `coriolis = "co"` (`spheroid.m` 46–52), cross-flow `drag_model = "cylinder"` (`remus100.m` 221) with `strip_grid = "midpoint"` (`crossFlowDrag.m` 56), `sway_damping_fade = false` and `smooth_speed = false` (MSS exactly, `remus100.m` 127, 218; with `true` the regularisation `smooth_speed_epsilon` would be a parameter: regularisation, not physics, not part of this set), propeller `open_water_model = "linearized"` (`remus100.m` 163–177) with `clip_advance_ratio = false`, fins `convention = "starboard_down_positive"` (`remus100.m` 242, 245). `stabilize_added_mass_coriolis` is chosen per test (`true` = MSS, `remus100.m` 207–210).

Gravity: `weight = MRB(1,1) · gravity(mu)`, `mu = deg2rad(63.446827)` (`remus100.m` 96–97, 214; `INS/functions/gravity.m` 11–12): g = 9.8217848645790813 m/s², m = weight / g = 31.900515266356955 kg (`spheroid.m` 35–36 with the scaled semi-axes). The class takes `weight` and `buoyancy` as values; nothing in the class computes them from the mass.

Density: since MSS `23c307e` (`mssConstants.m` 12) every call in `remus100.m` reads one water density, 1025 kg/m³ (body `spheroid.m` 35–36, cross-flow `crossFlowDrag.m` 36–37, added mass `imlay61.m` 31–32, lift and drag `forceLiftDrag.m` 26–27, propeller and fins `remus100.m` 98–99). The history set therefore equals the default set in density; only geometry and mass remain different. `cross_flow_water_density` stays in the history file as the same number.

## Numbers and conventions

- `%.17g` in every CSV and in the JSON: each number reads back to the double MATLAB held.
- Units SI; angles rad; propeller speed rpm; NED positions, Euler ZYX angles, BODY velocities (`remus100.m` 7–24); current as `remus100.m` 40–43: `v_c = [Vc cos(betaVc − psi), Vc sin(betaVc − psi), w_c]` in BODY.

## sha256

History set:


- `generate_remus100_mss.m` `3b51468ebafba53e37b19a46f865384b55007cb0c2522b191bd840fd733d7dfa`
- `remus100_derivative_mss.csv` `f41d86378cf79ec338cb638d2ff37c772c8af77f536c1d6f5169543071df8a21`
- `remus100_hull_derivative_mss.csv` `7fcf15ce0e8375912ae99d9d1f29eeba2d970a3bfcd6f57867deac08855a61f9`
- `remus100_hull.m` `b8934fb6fff3ee032d567c77f75e126142e0c617babf2fbb30af95fd8644253f`
- `remus100_hull_munk.m` `285db255b3bd5a5b67ed1721e7970a4dbf26df923e6ed9c4073a26bd8c2a4bb0`
- `remus100_mass_matrix_mss.csv` `28039234999b9413cb75cf343a240d0fcb6e1d6e80609036b1f07c302183fcef`
- `remus100_munk.m` `bb4eac8889b195a842d809902b7d883f834dd780c151cdf46f93550139885d54`
- `remus100_parameters.json` `f2782457d1ea684fc5f0c8bca5bcfa6c4d3ea4ea06c09ee3d1f96e0ff29e1ddb`
- `remus100_trajectory_dive_and_level_off.csv` `6d8ae2dd1f05c814a899e033b8a78fe9f30a3b65adaa53fe4a229fc8933f4ee1`
- `remus100_trajectory_hull_glide.csv` `774a3f7d3c9d321a8a379fa23693ce0efe9c277cbb795ada98685123c390eb62`
- `remus100_trajectory_hull_munk_glide.csv` `a539247c594b2d6633c1c86d33de1e7a58538a91f46ddeeaea294212d6769978`
- `remus100_trajectory_hull_roll_pitch_decay.csv` `180bd2ee23bc807f3bc60f43c5283fef57d3cb1024df05c0a991c8608b33cf99`
- `remus100_trajectory_hull_sway_yaw_release.csv` `9bc46e5d98e56c30c5af1168dbe701e08e3dc1c0f44ce78d99740faa31105142`
- `remus100_trajectory_hull_zoh_rudder_step_20deg.csv` `62a79321eece81747dea7414ffe50b0eb2a682a5dfa17e2c3d0add9d9cdd6499`
- `remus100_trajectory_hull_zoh_stern_plane_step_10deg.csv` `9e6ab41268348d9bbe500d792335163e68be69268426e14354ab4c9e27ec4691`
- `remus100_trajectory_hull_zoh_zigzag_20_20.csv` `a482ee38d4a3af812fd9834feb2fee2b1c2551f6a228fe73f23cb490a3280b62`
- `remus100_trajectory_munk_rudder_step_20deg.csv` `f09ebb869039f304544b892b566b2ba26b941307c3b935398f461648931ad899`
- `remus100_trajectory_munk_straight_perturbed.csv` `8d3240fe9b9628a8e048b8728b71ae255683744548b3614b2d9e18d56ae09981`
- `remus100_trajectory_rudder_step_10deg_in_current.csv` `2ef207ec8ecbf61a08b36b25d3779508405ba63dc75d1ccd950129ba1c16c850`
- `remus100_trajectory_rudder_step_20deg.csv` `ca67686dc75a66693ba205e554a468f4161c2e2b0db1e80589d15861aa5e811b`
- `remus100_trajectory_rudder_step_4deg_at_1p54ms.csv` `063d37a271d7e87ba113a3587d178e46a644c034389e40819c084ce4e0dd3ea8`
- `remus100_trajectory_stern_plane_step_10deg.csv` `ae34f953a7554ba5460ecb9968d8e02759a1501016d3d83019ac67953c0434c3`
- `remus100_trajectory_straight_from_rest.csv` `1972df1a5b12bc40fc6caf0119623423b37200b2341aca00a9789e43da02d379`
- `remus100_trajectory_tolerances.csv` `fa02baa7b1acea09d3f924a8fd7ae35c8e3f0c98595bd49c5b4fa47fbab311f0`
- `remus100_trajectory_zigzag_20_20.csv` `dbc2c4b16f4f195763488bcf5fa085e3a6ebce45d4d8e072832f5e37f85ebc5d`

Default set:

- `crossFlowDrag_consistent.m` `13b83c72a4f9299ff785bdd4a45396af45041fba0c05b68ac0396a782c3aa804`
- `forceLiftDrag_consistent.m` `3b803b7c88f92f7293c12110fe03eee36751ee03c46fb984ef1d887dff729075`
- `generate_remus100_consistent.m` `b1887e88a0cd441a4137b772161e52417ce3a0d6b07606a38dd46f4e11b8b338`
- `imlay61_consistent.m` `d0438c28d856bc977a412941512d7f7badc9456cf7dd5fd3a57df102958fdcbb`
- `remus100_consistent.m` `02289e99a01a063707ed7cc0280de6de32c148c0d8878803793bf76af185222b`
- `remus100_derivative_consistent.csv` `7c6158fff64a66b1e7fb9848e9c4c13d1e4a2e0080e8a6f08bd8188b2ff4d7b1`
- `remus100_hull_derivative_consistent.csv` `76e2a3eb5f0b476fce7515eff2701b61b5ca0a34e5dab6b196dac9bf4d85ac8c`
- `remus100_mass_matrix_consistent.csv` `2deb95df017037a8140f4bcfa07c31ec1e00dfd2d6b8c8a8c1c09c88f2dbbfac`
- `remus100_parameters_consistent.json` `94cbaa62cd83e59b726c949e78b70b6cccb0c6f515629c9a00e822c045634f4c`
- `remus100_trajectory_dive_and_level_off_consistent.csv` `5dd39313313beca1d372b0911864eca0d09373230dcab966abe036a34fd79297`
- `remus100_trajectory_hull_glide_consistent.csv` `b31e119cb68fbc1b64f96b6574de575f6b8b29a2d0d41ddaac7e29708a9c1bc4`
- `remus100_trajectory_hull_munk_glide_consistent.csv` `0e5e659711b5a93c2df9074aba6afe9101ec2d793020262cea9c1519f6ff2816`
- `remus100_trajectory_hull_roll_pitch_decay_consistent.csv` `bf31f709902d0f4edc54afe78f665313c207fb455c864c4626615c87200b7906`
- `remus100_trajectory_hull_sway_yaw_release_consistent.csv` `a525882d45f1eeb1762d93b0806d94b72f8a6a9c6a32f31113f6c70f7de81579`
- `remus100_trajectory_hull_zoh_rudder_step_20deg_consistent.csv` `da27a5a15c37fa414ee87233b92aafae3484fd176beb17e0f0d2b9b6c4108e2a`
- `remus100_trajectory_hull_zoh_stern_plane_step_10deg_consistent.csv` `f552a31608d6f0f239f629d3f1474c58c93cbb64573001d1d5a55b20df8f7a7e`
- `remus100_trajectory_hull_zoh_zigzag_20_20_consistent.csv` `c3bfc4b5e0e4a0cbfee4f7c1185dd264bb833d52cec45e839002b9da95538c41`
- `remus100_trajectory_munk_rudder_step_20deg_consistent.csv` `5eb11fa5a9451644bcdbc5b2a476d4835962051ef64f9e8d5779be8b13e04371`
- `remus100_trajectory_munk_straight_perturbed_consistent.csv` `1699b48578a5e1fd854b059c0f695b943b76a4785adb8d8c2ea88dd33dd2c04b`
- `remus100_trajectory_rudder_step_10deg_in_current_consistent.csv` `12768ae53c4f1bac31d01caecea39ffae87ae68e424e29914601bae7a70ae59d`
- `remus100_trajectory_rudder_step_15deg_at_925rpm_consistent.csv` `29527944e5e8e487905a15db947990dd9a7a2e253643b7ccdfdb4f224628edbd`
- `remus100_trajectory_rudder_step_15deg_at_925rpm_mss.csv` `03fd043f6241b407ba5971c18668d9885443a6d8784f818f5a32405226bd4d39`
- `remus100_trajectory_rudder_step_15deg_at_925rpm_one_density.csv` `03fd043f6241b407ba5971c18668d9885443a6d8784f818f5a32405226bd4d39`
- `remus100_trajectory_rudder_step_20deg_consistent.csv` `fdabce6078c008e593b06d2920056b27bb865c1edd3328e322dfeb8cf9d81236`
- `remus100_trajectory_rudder_step_4deg_at_1p54ms_consistent.csv` `5b7419cbf2d6dffb8ead92507c61f8ed49e77b96ce8fb384a78ca4defc594c26`
- `remus100_trajectory_stern_plane_step_10deg_consistent.csv` `3ea45fac9a9731413cf25748594637bb85a6a68b335038bc187833f89094e170`
- `remus100_trajectory_straight_from_rest_consistent.csv` `bc74899d7298bec604e01720f677109468b0b28279aadc2a029a96fbc78e7b3a`
- `remus100_trajectory_tolerances_consistent.csv` `ad303c6e69435cb44683a3ed7d41b477503c240c71e64d19f122a500c74ed229`
- `remus100_trajectory_zigzag_20_20_consistent.csv` `ffcfa2613d418c6e1aecad03cb24d480c4bfbfdd222f70c072ed3869b118e9c3`
- `spheroid_consistent.m` `3adca8c5a3ab576716e9e2f38fed1a54763d9f9b576c0bff858d55148525fa7d`

## Readers

`tests/vehicles/test_torpedo_class.py`.

Author: Enio Krizman · Date: 2026-10-07


