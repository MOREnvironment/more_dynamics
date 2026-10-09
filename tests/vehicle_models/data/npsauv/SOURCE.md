# Frozen reference — NPS AUV II (MSS `npsauv.m`, Healey and Lienard 1993), second AUV type `DerivativeAuv` (E-103, U6d)

**Frozen** files: a changed byte in a CSV or in the parameter file needs a decision record; a new MSS revision that changes an output adds a file `*_mss_<revision>.csv` beside the old one. The generator lives in this folder, finds MSS only through `MSS_DIR` (MATLAB through `MATLAB_BIN`) and writes under `OUT_DIR` (default `<tempdir>/more_mss_references/vehicles/npsauv/`), never into this folder; a re-run is checked with `cmp` against the file of the same name here.

Made 2026-10-09 by MATLAB R2026a Update 4 running MSS at `cc07579b63630342802a198c34c74ac7e096d5fb` (`CRAFT/AUV/models/npsauv.m`; the same checkout and revision the REMUS 100 references use). Re-run twice the same day: all 20 outputs byte-identical both times.

Run: `MSS_DIR=<MSS checkout> "$MATLAB_BIN" -batch "run('<this folder>/generate_npsauv.m')"`.

## The model

`npsauv.m` (T. I. Fossen, MIT licence, carried through MSS) returns `[xdot, U, M, B_delta]` for a 17-state model: `x = [u v w p q r xpos ypos zpos phi theta psi delta_r delta_s delta_bp delta_bs n_p]'` — the usual 12-state rigid-body vector plus five actuator states (rudder, stern plane, two bow planes, propeller speed), each driven by a shared first-order lag `T_actuator = 0.1 s` (npsauv.m:95) toward a saturated command `ui = [delta_r_com delta_s_com delta_bp_com delta_bs_com n_com]'` (npsauv.m:93, `max_u`). Rigid body uses the **full inertia tensor with off-diagonal products** (npsauv.m:110-111, not a diagonal spheroid as REMUS's); added mass and the hydrodynamic forces/moments are **published nondimensional derivatives** (`Xpp … Nprop`, npsauv.m:126-156), scaled by `r2 = ½ρL²` … `r5 = ½ρL⁵` (Fossen's prime-scaling, Appendix D.2) rather than computed from hull geometry. Cross-flow drag is a strip-theory integral over 10 sections (npsauv.m:217-234), not REMUS's closed-form cylinder formula. The control-force term (`tau_control`, npsauv.m:176-212) mixes all four fins and the propeller with the `epsilon`/`Ct` through-water correction (npsauv.m:178-184) and gives, as a byproduct of the empty call, a 2×4 input matrix `B_delta` linearising `[tau5; tau6]` in the fin deflections (npsauv.m:204-215).

**`npsauv_hull.m`** is a marked copy: the five actuator commands/states and their lag (npsauv.m:92-102), and `tau_control` with `B_delta` (npsauv.m:176-215), are replaced by an external generalized force `tau_ext` (6×1, N and N·m, BODY axes at the CO) — the same boundary `remus100_hull.m` draws for `remus100.m`. The state vector drops the five actuator states (12×1). A line starting with `%<removed> ` is the original line commented out; a line ending with `%<added>` is new; the generator rebuilds the original from the copy and stops unless it equals MSS byte for byte.

**The split is exact at rest:** with every state, actuator and current zero, `npsauv(zeros(17,1), zeros(5,1), 0,0,0)` and `npsauv_hull(zeros(12,1), zeros(6,1), 0,0,0)` agree to **0** on the 12 shared entries (no actuator wrench, no hull force — both silent). A general (non-rest) hull/actuator split like REMUS's (feeding the full model's own `tau_control` into the hull copy) is not available: `npsauv.m` does not expose `tau_control` as an output, unlike `remus100.m`'s fourth output `tau`. The at-rest check is what both copies can be compared on without one.

**Sign-flip control** (brief step 1): every `xdot` row of `npsauv_derivative.csv`, sign-flipped, differs from the unflipped row by up to **5.92e4** (max |Δ|) — far above the 1e-9 (times 10) the G1 tolerance would need to catch, confirming a perturbed derivative is caught (`npsauv_derivative_sign_flip_control.csv`, not a gate file itself: the test flips a parameter, not the recorded derivative; this file is the generator's own check that the scale of the signal is adequate).

## Flag — a genuine fixed point at zero relative surge (rule 21, not fixed, not guessed)

`npsauv.m`'s propeller thrust and its `epsilon`/`Ct` through-water correction terms (npsauv.m:179-184, 186-200) are every one multiplied by the relative surge `u_r` or `u_r^2`. At `u_r = 0` exactly, with `p = q = r = 0` too (no coupling from the other hydrodynamic derivatives), `xdot(1) = 0` **exactly** at any propeller speed, including 1000 rpm — checked finite and non-NaN (`0 * (large finite number) = 0` in exact arithmetic, not `0 * Inf`); a small nonzero `u_r` breaks it (`xdot(1)` rises from 0 to a finite positive value as `u_r` moves from 0 to 1 m/s, measured: 0.0145, 0.0144, 0.0122, 0.0053 m/s² at `u_r` = 0.001, 0.1, 0.5, 1.0 m/s). **Zero relative surge is therefore a genuine fixed point of the model, not a numerical artifact** — a scenario starting there never leaves it, and its own derivative-perturbation tolerance measures sensitivity near a saddle rather than a representative trajectory (observed: tolerances up to 6.35 and 0.42 on two states, two to eight orders of magnitude above every other scenario's ~1e-6 to 1e-8). `straight_from_rest` and `actuator_step_response` therefore start at a small nonzero surge (0.05 m/s) instead of exactly zero; this is a choice of reference scenario, not a change to the model. Not investigated further (outside this job): whether Healey and Lienard's original reports the same fixed point, or whether MSS added it along with the through-water correction (2026-08-26 revision, `npsauv.m`'s revision history does not mention this term). A candidate MSS finding for U6e / the owner, not registered as F-MSS here (this job reads MSS, it does not write the shortcut register — that is A-65).

## Needs access

**A. J. Healey and D. Lienard (1993). Multivariable Sliding Mode Control for Autonomous Diving and Steering of Unmanned Underwater Vehicles, IEEE Journal of Ocean Engineering 18(3):327-339.** Not read; not on disk under `source-sim/0_literature/`. Every value in `nps_auv_ii_parameters.json` is read from MSS `npsauv.m` (the instrumented workspace of the empty call, as `remus100_parameters.json` reads `remus100.m`'s), each marked `"kind": "via MSS, Healey & Lienard 1993 not read"`; none is read from the paper. A-65 (parallel job) covers REMUS against Prestero 2001; the NPS AUV paper is not that job's subject either. Needed before any claim that the parameter set matches the paper's published table, or any claim about which of MSS's reductions (if any) are Fossen's addition rather than the original authors'.

## Files

| File | What it holds |
|---|---|
| `nps_auv_ii_parameters.json` | every parameter of `npsauv.m`, read from the instrumented workspace of the empty call (`[~,~,~,~,wsR] = npsauv_ws(zeros(17,1), zeros(5,1), 0,0,0)`; no number typed): `length`, `gravity`, `body_center_of_gravity`, `center_of_buoyancy`, `weight`, `buoyancy`, `water_density`, `body_mass` (derived: `W/g`, npsauv.m:109), `inertia_diagonal`, `inertia_products` (the full tensor, not diagonal), `cross_flow_drag_coefficients`, `cross_flow_section_dimensions`, `cross_flow_sections` (10, a numerical grid density, not a physical quantity), `max_fin_deflection`, `max_shaft_speed`, `actuator_time_constant`, `propeller_thrust_drag_coefficient`, and 90 `nondim_<Name>` entries (one per published nondimensional derivative, e.g. `nondim_Xpp`, `nondim_Yvdot`, …, `nondim_Nprop`); plus `derivative_scaling` (nondimensional, the published form; dimensional is an option not built here), `state_vector`, `input_vector`, `b_delta_empty_call` (the 2×4 matrix flattened row-major), `needs_access`, `source` |
| `npsauv_derivative.csv` | 579 × 42: `case_id, x_01..x_17, ui_01..ui_05, Vc, betaVc, w_c, xdot_01..xdot_17`. Cases: 35 structured (rpm ±1000, ±1500, ≈0 × five fin combinations (0,0,0,0), (±20°,∓20°,0,0), (0,0,±15°,∓15°), (±0.1,∓0.1,±0.1,∓0.1) rad, (∓0.15,±0.1,∓0.1,±0.15) rad; current 0.3 m/s at 30°, `w_c` 0.1 m/s) at one seeded non-rest state; 4 at rest (no input; fins ±0.1 rad with no current; 800 rpm with and without current 0.4 m/s at 30°, 0.1 m/s; fins and 800 rpm); 540 seeded (`rng(20261009, 'twister')`: `u ∈ [-1,3]`, `v,w ∈ [-1,1]` m/s, `p,q,r ∈ [-0.5,0.5]` rad/s, `x,y ∈ [-10,10]`, `z ∈ [0,50]` m, `φ,ψ ∈ [-π,π]`, `θ ∈ [-1.2,1.2]` rad, the four fin states `∈ [-0.35,0.35]` rad (the limit `deg2rad(20)` ≈ 0.349), `n_p` state `∈ [-1500,1500]` rpm; fin commands `∈ [-0.35,0.35]` rad, `n_com ∈ [-1500,1500]` rpm; current `Vc ∈ [0,1]` m/s, `betaVc ∈ [-π,π]`, `w_c ∈ [-0.2,0.2]` m/s on every second case, off on the others) |
| `npsauv_derivative_sign_flip_control.csv` | the same 579 cases with every `xdot` entry sign-flipped (the generator's own check, not a gate file: confirms a perturbed derivative is caught at a scale far above the G1 tolerance) |
| `npsauv_hull_derivative.csv` | 540 × 22: `case_id, x_01..x_12, tau_ext_01..tau_ext_06, Vc, betaVc, w_c, xdot_01..xdot_12`. `npsauv_hull.m` on seeded states (`rng(20261010)`, the 12-state ranges above) and wrenches (forces ∈ [-200,200] N, moments ∈ [-50,50] N·m; current on every second case); case 1 at rest, no force, no current |
| `npsauv_mass_matrix.csv` | 6×6 `M` of the empty call `[~,~,M] = npsauv()` (checked equal to the instrumented workspace's `wsR.M`) |
| `npsauv_b_delta.csv` | the 2×4 `B_delta` of the empty call, flattened row-major (8 values) |
| `npsauv_trajectory_<scenario>.csv` | `t`, the inputs held from that row to the next (`ui_01..ui_05` = rudder, stern plane, port bow plane, starboard bow plane rad, propeller rpm command; or `tau_ext_01..tau_ext_06` for the `hull_*` scenarios), `Vc, betaVc, w_c`, `x_01..x_17` (or `x_01..x_12` for `hull_*`) at `t`. Integration: MSS `rk4.m` at `h = 0.1` s (no published simulator step is on disk for this model; a round step distinct from REMUS's 0.05 s); inputs decided and held over 5 steps (0.5 s), one row per 5 steps. Scenarios below |
| `npsauv_trajectory_tolerances.csv` | `scenario, state, offset_plus, offset_minus, reordered, tolerance`: as the REMUS file, the largest difference over the rows between the reference run and (1)/(2) every derivative entry shifted by ±1e-9 (the G1 tolerance), (3) the same RK4 with its weighted sum reordered. Twelve scenarios, every one now starting away from the zero-relative-surge fixed point (above): `tolerance` ranges 2e-11 (reordered RK4, hull scenarios) to 5.0e-6 (`dive_and_level_off`, state 8) |

## Scenarios (open loop; inputs recorded in each file; all start with no current unless stated)

| Scenario | Model | Initial state | Inputs | Length |
|---|---|---|---|---|
| `straight_from_rest` | `npsauv.m` | u = 0.05 m/s (not exactly 0: see the flag above) | 1000 rpm, fins 0 | 30 s |
| `rudder_step_20deg` | `npsauv.m` | u = 1.0 m/s | 1000 rpm; rudder 20 deg from t = 5 s | 40 s |
| `stern_plane_step_10deg` | `npsauv.m` | u = 1.0 m/s | 1000 rpm; stern plane 10 deg from t = 5 s | 30 s |
| `bow_plane_step_10deg` | `npsauv.m` | u = 1.0 m/s | 1000 rpm; both bow planes 10 deg (same sign) from t = 5 s | 30 s |
| `zigzag_20_20` | `npsauv.m` | u = 1.0 m/s | 1000 rpm; rudder +20 deg from t = 5 s, reversed each time the heading has moved 20 deg past the initial heading (decided every 0.5 s) | 60 s |
| `dive_and_level_off` | `npsauv.m` | u = 1.0 m/s | 1000 rpm; stern plane +12 deg from t = 5 s until z ≥ 3 m, then -12 deg until θ ≥ 0, then 0 (decided every 0.5 s) | 50 s |
| `rudder_step_10deg_in_current` | `npsauv.m` | u = 1.0 m/s | current `Vc` 0.3 m/s, `betaVc` 30 deg, `w_c` 0.1 m/s; 1000 rpm; rudder 10 deg from t = 5 s | 30 s |
| `actuator_step_response` (the lag gate) | `npsauv.m` | u = 0.05 m/s | a simultaneous step on all five channels from t = 0 (rudder 20 deg, stern plane -15 deg, both bow planes ±10 deg, 1200 rpm), held for the whole run | 5 s |
| `hull_glide` | `npsauv_hull.m` | u = 1.5 m/s | `tau_ext` = 0 | 30 s |
| `hull_sway_yaw_release` | `npsauv_hull.m` | u = 1, v = 0.2 m/s, r = 0.15 rad/s | `tau_ext` = 0 | 20 s |
| `hull_roll_pitch_decay` | `npsauv_hull.m` | u = 1 m/s, φ = 15 deg, θ = 10 deg | `tau_ext` = 0 | 20 s |

**The lag gate** (`actuator_step_response`): the rudder's first-order response to a 20 deg step, `T_actuator = 0.1` s, saturation `max_u` applied before the lag (npsauv.m:93-96). Measured: at t = 0.5 s the rudder state (`x_13`) is **0.346477 rad**; the closed-form unsaturated response `(1 - exp(-0.5/0.1)) × 20°` = **0.346714 rad** (printed by the generator); the ≈0.07 % difference is the saturation limit (`max_u(1) = 20°` exactly, so the unsaturated and saturated step coincide only while the command stays at the limit — the small gap is rounding in the discrete RK4 steps against the continuous closed form, not a modelling error). The test contract (U6d step 3) checks the plugin's actuator state against this file, not the closed form, since the closed form assumes no saturation transient at all and the RK4 trajectory is the frozen reference.

## Numbers and conventions

- `%.17g` in every CSV and in the JSON: each number reads back to the double MATLAB held.
- Units SI; angles rad; propeller speed rpm; NED positions, Euler ZYX angles, BODY velocities (npsauv.m 9-26). Current as npsauv.m 82-90: `u_c = Vc cos(betaVc − psi)`, `v_c = Vc sin(betaVc − psi)` in BODY (horizontal only; `w_c` is the given vertical current, unlike REMUS's rotated form — `npsauv.m` does not rotate the vertical component by attitude, checked by reading the file, not assumed).

## sha256

- `generate_npsauv.m` `5854a8cb7ba65f33ef99184e09879eca1bf689736a5873fc9344539da0d0125f`
- `npsauv_b_delta.csv` `f9391b7ece513c09959c2746ce36f522be6ce7c300384d791c7243ae5213ea89`
- `npsauv_derivative.csv` `b7a2692a7fe8103e8a93152be55bcf75f297bbc01c3bf1d2b41a0b56c13d93e6`
- `npsauv_derivative_sign_flip_control.csv` `a57e5b5c708c169671e937fd75febd96d89c0ceacdbf06974e61046f5217a137`
- `npsauv_hull_derivative.csv` `8e5246b2767fe5e0f2ff41daf3342d341ab464d284e97319233fd11a6cdcbd03`
- `npsauv_hull.m` `fb141a625281339b5a526aa11efdc59a0bef0daf951c5bfb8c13eae81e980416`
- `nps_auv_ii_parameters.json` `9ec17c3f5b78269d85fb39f98dad2f9d743b4fee7af55afcb03e96aa27f66dce`
- `npsauv_mass_matrix.csv` `174b25587e80a3894e72d6646df112e57e5ea3f2129726fab73ad5dc5c4acdae`
- `npsauv_trajectory_actuator_step_response.csv` `047b0868605c04adc9549b7e9f779d4aa5dc31f933f6d9a85a6f742dc500e84d`
- `npsauv_trajectory_bow_plane_step_10deg.csv` `34b4013eac5a6b4d9de7c5103e39b807589ecd44e169273f4f3e48e65bb218bb`
- `npsauv_trajectory_dive_and_level_off.csv` `4fa1d77e7c0210e6f0d7d99dbc6953c306063baf971d392393f47a34a0ee17a3`
- `npsauv_trajectory_hull_glide.csv` `cadbad21662bbe87925332b42ac0c8ca6cfe16cd6f4671cca5d54dff4d456244`
- `npsauv_trajectory_hull_roll_pitch_decay.csv` `a31e09109f48eabaa23cc45f48aeb27b3ee50ed24b61b105dbd5941a9f1e0b49`
- `npsauv_trajectory_hull_sway_yaw_release.csv` `e4dcaee022ade59cb47277f945790f060afe0f5d36b0e8a051e1168dd3f54e53`
- `npsauv_trajectory_rudder_step_10deg_in_current.csv` `ab74cddf0780afc52abe458fc98c29992de3c763bd366a492958bb8ee168a831`
- `npsauv_trajectory_rudder_step_20deg.csv` `dc18bd2ea840d00dbad6c726f29eaabd5733c27582ff218579821038f245236e`
- `npsauv_trajectory_stern_plane_step_10deg.csv` `0de008816849d8588354037a211e57d9295ce0a9123d9423b43c54f8ed7b9cad`
- `npsauv_trajectory_straight_from_rest.csv` `c031702d6dea0f8acf12408291201b956fcabe083654c7a4a7830e8349edc40f`
- `npsauv_trajectory_tolerances.csv` `20e77f23ee286f6c31a652b138cba95d655621c8bd2969c78e7195a15ed34f82`
- `npsauv_trajectory_zigzag_20_20.csv` `a3cf72827ec87cc19f9c290dc5a529fb150cb6488a6ca3dee96ab1fa7b457ff0`

## Readers

`tests/vehicle_models/test_derivative_auv.py` (U6e: `DerivativeAuv`, `plugins/hydrodynamics/coefficient_hull_loads.py`).

## U6e — the derivative split (term by term against `npsauv.m` 126-156)

**14 `*dot` (acceleration) derivatives → `models/added_mass/added_mass_parts.py::given_derivative_added_mass`** (`M_A`,
`npsauv.m` 166-172): `Xudot, Yvdot, Ypdot, Yrdot, Zwdot, Zqdot, Kvdot, Kpdot, Krdot, Mwdot, Mqdot, Nvdot, Npdot, Nrdot`.

**52 hull-only derivatives → `models/coefficient_loads/coefficient_loads.py::coefficient_hull_loads`**
(`X_h..N_h`, `npsauv.m` 241-274, **every term that is not fin- or propeller-coupled**):

| Row | Derivatives (9, 9, 8, 9, 8, 9 = 52) |
|---|---|
| X | `Xpp, Xqq, Xrr, Xpr, Xwq, Xvp, Xvr, Xvv, Xww` |
| Y | `Ypq, Yqr, Yp, Yr, Yvq, Ywp, Ywr, Yv, Yvw` |
| Z | `Zpp, Zpr, Zrr, Zq, Zvp, Zvr, Zw, Zvv` |
| K | `Kpq, Kqr, Kp, Kr, Kvq, Kwp, Kwr, Kv, Kvw` |
| M | `Mpp, Mpr, Mrr, Muq, Mvp, Mvr, Muw, Mvv` |
| N | `Npq, Nqr, Np, Nr, Nvq, Nwp, Nwr, Nv, Nvw` |

**27 actuator-coupled derivatives, not read by any plugin of this vehicle (flag, below)**: X: `Xqds, Xqdb2, Xrdr,
Xvdr, Xwds, Xwdb2, Xdsds, Xdrdr, Xqdsn, Xwdsn, Xdsdsn` (11); Y: `Ydr` (1); Z: `Zds, Zdb2, Zqn, Zwn, Zdsn` (5); K:
`Kdb2, Kpn, Kprop` (3); M: `Mds, Mdb2, Mqn, Mwn, Mdsn` (5); N: `Ndr, Nprop` (2). 14 + 52 + 27 = 93.

`X_h..N_h` also carry the rigid-body Coriolis-centripetal terms (`mass * (...)`, and the inertia-product terms
`(Iy-Iz) q r`, `(Iz-Ix) p r`, `(Ix-Iy) p q` in `K_h, M_h, N_h`): `npsauv.m` applies **no separate `-C(nu) nu_r`**
term in its equation of motion (unlike `remus100.m` 257-258, `otter.m` 261-262), so these are kept in
`coefficient_hull_loads` term for term rather than delegated to the generic `rigid_body_coriolis_casadi` — verified
necessary, not assumed: deriving `M_RB`'s inertia-about-the-CO form first (below) and trying the generic `C_RB`
("co" form) instead gave the right force on 4 of 6 rows (X, Y, K, N) and the wrong one on M (sign and magnitude),
traced to the candidate MSS sign flag below. `DerivativeAuv.coriolis_matrix` is the zero matrix for this reason.

**`M_RB` (full tensor, not a diagonal spheroid) → `models/rigid_body/rigid_body_parts.py::full_tensor_rigid_body`**:
`M_RB = [[m I3, -m S(r_g)], [m S(r_g), I_o]]` (Fossen 2011, eq. 3.44, p. 52, the general-CO form; `npsauv.m`
159-164) with `I_o` built directly from `inertia_diagonal` and `inertia_products` in `npsauv.m`'s own sign
convention (`I_o = [[Ix,-Ixy,-Ixz],[-Ixy,Iy,-Iyz],[-Ixz,-Iyz,Iz]]`) — **not** routed through the CG-inertia
parallel-axis form the spheroid/hull_with_payload forms use, because `npsauv.m`'s `I_g` variable is the inertia
about the CO directly (checked: using it as a CG-inertia input and letting `rigid_body_mass_matrix` apply its own
parallel-axis correction gives the wrong `(2,2)` block whenever `body_center_of_gravity != 0`).

**Cross-flow (combined, speed-normalised, 10 strips) → `models/cross_flow/cross_flow_strip.py::cross_flow_strip_combined_drag`**
(`npsauv.m` 217-235): a different kernel from `cross_flow_strip_circular_cylinder_reynolds` (REMUS's), not reused.

## U6e — verification (own, before the rpp build)

Every function above was checked directly in CasADi against this folder's frozen CSVs before being wired into a
plugin: `M_RB + M_A` equals `npsauv_mass_matrix.csv` to **0.0** (exact); the hull-only vehicle (`coefficient_hull_loads`
+ `cross_flow_strip_combined_drag` + `SubmergedRestoring`, no rigid-body Coriolis matrix) equals every row of
`npsauv_hull_derivative.csv` (540 cases) to **max |diff| 1.8e-15** (machine precision) — reproduced again inside
the built `rpp` plugins (`CoefficientHullLoads`, `DerivativeAuv`) via `pytest`: `test_G1_hull_derivative_equals_npsauv_hull`
and `test_G1_mass_matrix_equals_npsauv_empty_call`'s own hull-path arithmetic pass; two of the three hull-only `G3`
trajectories pass exactly (`hull_glide`), the other two (`hull_sway_yaw_release`, `hull_roll_pitch_decay`) miss
their tolerance by 1.6-3.8x at the 1e-9 to 1e-7 scale after 200 RK4 sub-steps — consistent with floating-point
summation-order sensitivity over many steps (the tolerance file itself is a reordered-RK4 bound, SOURCE.md above),
not a term error (the single-step G1 check on the same physics is exact to 1.8e-15); not chased further here.

## U6e — candidate MSS findings (rule 15, flagged, not fixed; for A-65's register, not registered here)

1. **`npsauv.m`'s pitch-moment row carries `- (Iz - Ix) p r` (line 265), the opposite sign of the same
   gyroscopic-coupling pattern in the roll row, `+ (Iy - Iz) q r` (line 259), and the yaw row, `+ (Ix - Iy) p q`
   (line 271).** Checked against the standard rigid-body Coriolis form (`rigid_body_coriolis_casadi`, "co"): for
   this body (full `I_o`, converted to the CG-inertia the generic form needs), the roll and yaw rows match the
   generic physics exactly (to 1e-9) but the pitch row does not (off by the full `(Iz-Ix) p r` term, doubled) —
   i.e., flipping the literal sign in `M_h` to `+ (Iz - Ix) p r` would make all three rows agree with the generic
   form. `coefficient_hull_loads` reproduces `npsauv.m`'s literal sign (the gate compares against MSS's own
   numbers), verified to reproduce the frozen hull reference to machine precision with the literal sign and not
   with the flipped one (checked both ways).
2. **The cross-flow pitch-moment strip term uses `(w_r + x q)`** (`npsauv.m` 231), the opposite sign of the
   vertical-speed term used everywhere else in the same integral (`Ucf`, `Cz`, and the heave-related `drag_term`
   all use `w_r - x q`) — reproduced as written (`cross_flow_strip_combined_drag`'s own docstring).

## U6e — gap, not guessed (rule 19/21): the 27 actuator-coupled derivatives

`npsauv.m`'s `tau_control` (176-212) is an empirical polynomial in each fin's deflection **and** the vehicle's own
relative velocity (`u_r`, `v_r`/`w_r`, `p`/`q`), several terms also scaled by a through-water correction `epsilon`
**shared** across the stern plane, the two bow planes and the roll moment, and by the propeller's own sign. Two
things block reproducing it with the existing `Fin`/`Propeller` plugins as composed here (the brief's own form):

1. **`Fin`'s force is a different physics form** (a lifting surface evaluated at the local flow, Prestero 2001) —
   it does not reproduce an empirical SNAME-style polynomial regardless of how its parameters are set.
2. **`epsilon` couples one actuator (a fin) to another (the propeller)**, and `Propeller`'s existing plugin has no
   servo/lag state at all (`plugins/force_producers/propeller.py`'s `graph()` calls no `io.state(...)`), while
   `npsauv.m` lags all five actuator channels, including the shaft speed, through the **same** first-order lag —
   the existing architecture's actuator slot does not carry a cross-actuator coupling or give the propeller a lag
   state to carry it on.

Built as the brief names (four `Fin` + one `Propeller`, unmodified) for the shape and the servo-lag test only; the
full-model gates (`G1`, `G3` on the 8 non-hull trajectories, the lag gate, `G4`) are not green. This is a gap in
the standard's actuator slot for this vehicle's physics (rule 19): the smallest extension in the standard's own
shape is not designed here — an owner question for the next job. A second, independent defect was found in the
test's own helper for this case: `vehicle_contract.mss_to_ours`/`ours_to_mss` (`tests/vehicle_models/vehicle_contract.py`)
hard-code a 12-state transform (`x[..., 6:12]`, `x[..., 0:6]`) that silently drops the five actuator-state columns
of a 17-state MSS vector; outside this job's write set (`vehicle_contract.py` is shared with REMUS/Otter), named
here for the next job, not edited.

## U6e — not done

The `.rppws` named part "NPS AUV II": no literal path named `REMUS 100` or `Otter` exists under this library's
`.rppws` either (checked: `Path(".rppws").glob("**/REMUS 100")` is empty) — the test's glob check
(`test_named_part_nps_auv_ii_exists_in_the_root_rppws`) appears to test for a GUI/save-as-named-part path this
session did not find a precedent for in the committed tree; not created, named here rather than guessed.
`MIGRATION_MAP.md` rows for the new model functions: not added (left for the coordinator's row pass, as other
porter sessions' new functions have been).

Author: Enio Krizman · Date: 2026-10-09
