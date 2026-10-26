# Grethe

A named part of `more_dynamics::Monohull` (single-hull surface craft): NTNU Mariner 5 USV "Grethe", Maritime Robotics
build 40402, Pioner 17 ft hull, one steerable electric outboard. **Every value and form is the model file's**:
`CRAFT/USV/models/grethe.m` of the MSS fork branch `grethe-vehicle`, revision `679a6e8` (the file may still be revised;
the part changes when it does). The model file is a best effort of a non-hydrodynamicist: most numbers are engineering
estimates, and the kind of each is written beside it in the file and below. It is not validated against the boat (the
field level, G7a, is a later test, and `X-9` lists what to measure).

The library part and the model file are one model: the part, built by rpp's builder with its outboard and the actuator
equations of the file's caller (`SIMgrethe.m`), equals the file's plant to 1e-9 in the derivative of 300 seeded states and
on six open-loop runs (`tests/vehicle_models/test_monohull_grethe.py`), and every option of the file (`opt.leg`,
`opt.ramp`, `opt.x_ref`, the payload trim) equals a selector of the part's forms.

## Values: where each comes from and of what kind

Kinds are the model file's: **document** (a vendor drawing, data sheet or manual), **published** (a printed number or
formula), **logged** (fitted to her control-station logs of 2026: a fit to one data set, not a validation),
**fitted** (chosen to match a logged quantity), **estimate** (a starting value with no document behind it), **approved
guess** (an estimate the owner accepted as a placeholder, no source), **derived** (computed from other numbers).

| Value (part file) | Value | Kind | `grethe.m` line |
|---|---|---|---|
| `length` (root) | 5.2 m | document (A-43 M1, F003 p2, M001 p7) | 238 |
| `beam` (root) | 2.15 m | document | 239 |
| `hull_mass` (root) | 806 kg | document (FAT, +-20 kg) | 240 |
| `hull_block_coefficient` (root) | 0.234 | derived: m / (rho L B T), T = 0.30 m (document), rounded (0.2345 unrounded) | 241 |
| `hull_waterplane_coefficient` (root) | 0.41 | **fitted** to the logged roll period of 1.85 s (a hard-chine hull would be 0.7-0.8; it also lowers the heave stiffness, heave period 1.25 s, not logged) | 242 |
| `hull_center_of_gravity` (root) | [0, 0, 0.025] m | estimate | 244 |
| `radii_of_gyration` (root) | [0.35, 0.25, 0.25] (fractions of B, L, L) | estimate | 245-247 |
| `added_mass_coefficients` (root) | [-1.0, -1.5, -1.0, -0.2, -0.8, -1.7] | estimate (factors on the displaced mass; the file's positive factors, here with the library's sign) | 248 |
| `payload_mass`, `payload_position` (root) | 0 kg, [0, 0, 0] | gap: the 2026 field payload (USBL mount, battery packs) is not weighed; the file's default is no payload | 119, 213 |
| `longitudinal_center_of_flotation` (hydrostatics) | -0.3 m | estimate | 249 |
| `longitudinal_inertia_factor` (hydrostatics) | 1.0 | estimate | 250 |
| `reference_point` (hydrostatics) | [0, 0, 0] | convention: the restoring matrix at the body origin (the waterline midpoint, E-101 Q4 a); the file's `opt.x_ref` default | 349 |
| `wetted_surface_method` (hydrostatics) | regression_table | published: the hydrostatic term of the Radojcic et al. (2014) regression | 353, 433-441 |
| `time_constants` (hull loads) | [8, 4] s ([T_sway, T_yaw]) | estimate (T_yaw chosen with N_v for a yaw-stable open loop) | 251-252 |
| `sway_yaw_coupling` (hull loads) | [0, 2500] ([Y_r, N_v]) | estimate: Y_r = 0 not identified; N_v = 2500 N m s/m **chosen** (about 7 deg/s at 20 deg and 450 N), no source | 363-364 |
| `damping_ratios` (hull loads) | [0.3, 0.05, 0.5] ([heave, roll, pitch]) | heave and pitch: approved guess, no source (E-135); roll 0.05: logged, band 0.03-0.08 | 253-255 |
| `max_forward_thrust`, `max_speed` (hull loads) | 0 N, 2.40 m/s | the surge resistance is the ITTC and residual forms alone (a thrust here would count a second resistance, E-135); 2.40 m/s is the logged top speed, band 2.32-2.47 m/s, kept here as a record, not read by the damping | 361 (`Xu = 0`) |
| `time_constant` (hull loads) | 15 s | estimate (the linear surge time constant T_1) | 256 |
| `form_factor` | 0.1 | published (Fossen 2011, p. 125; `XuuITTC.m` 33) | 257 |
| `crossover_speed` | 2.0 m/s | published (MSS value) | 258 |
| `ittc_reynolds_floor` | 1e5 | published (`XuuITTC.m`) | 259 |
| `surge_added_mass_factor` | 2.7 | published (`addedMassSurge.m` 32-33, MSS) | 315 |
| `surge_resistance` (hull loads) | ittc_residual_ramped | ITTC 1957 friction line plus the RA14 residual coefficient, scaled to the hull by Froude similarity (a **construction** of the file, not an equation of the paper); below the regression's range (Fn_V < 0.6, 1.805 m/s) the residual is ramped from zero by a smoothstep (a construction, no source; without it the resistance steps from 57 N to 115 N at 1.805 m/s) | 386-388, 443-495 |
| `gravity`, `water_density`, `kinematic_viscosity` (site) | 9.81, 1025, 1e-6 | published / document | 235-237 |
| outboard (`Grethe outboard`, linked) | see that part's `SOURCE.md` | document / estimate | 262-271, 391-408 |

## Selectors (this part's choices)

| Selector | Value | Where | Why |
|---|---|---|---|
| `added_mass_mass_basis` | `displaced` | root | the added mass follows the displaced mass, hull plus payload (E-101 Q5 a; `grethe.m` 290-291, 315-322) |
| `current_form` | `none` | root | still water: the reference data holds current fixed |
| `site_form` | `given` | the part's own `Site` | no site derivation |
| `wetted_surface_method` | `regression_table` | hydrostatics | A-64 Section 3.3 R11 |
| `payload_trim_form` | `lever_arm_moment` (the plugin's default) | hydrostatics | the payload weight is in the displacement; only the moment of its lever arm shifts the equilibrium (`grethe.m` 414-419; E-109) |
| `manoeuvring_damping` | `time_constants_coupled` | hull loads | `Y_v = -M22 / T_sway` and `N_r = -M66 / T_yaw` computed, `Y_r`, `N_v` given (`grethe.m` 361-369): the two diagonal terms are not typed again |
| `surge_resistance` | `ittc_residual_ramped` | hull loads | the file's default (`opt.ramp` true) |
| `cross_flow` | `hoerner_strips` | hull loads | the strip integral of `crossFlowDrag.m` with the Hoerner section (`grethe.m` 382) |

## What changed with this part (2026-10-26) and why

The part before followed an earlier design (A-43, A-64, E-101), and the model file was written beside it. The owner chose that
the part follows the file (E-145 Q7: *"mss model on grethe is reference as its easier to test"*). Changed: `hull_waterplane_coefficient`
0.8 -> 0.41; yaw damping T_yaw 16 s and N_v = Y_r = 0 -> T_yaw 4 s and N_v 2500; `reference_point` [-2.24, 0, 0] -> [0, 0, 0];
the surge resistance gets the ramp; the part gets the outboard (it had none) and the payload trim. The earlier yaw damping,
with which the open loop is unstable in yaw, stays as the part "Grethe (yaw-unstable stress test)" for controllers.
`tests/vehicle_models/data/grethe/grethe_parameters.json` is the parameter set of the independent block gates written
before (U8a, numpy references), kept as it was; it is not the part's source.

## Known limits (from the model file's header; none is hidden here)

The residual regression is not valid below Fn_V = 0.6 (ramped, a construction); the roll stiffness is fitted, not physical
(cause not identified); validity Fn_V <= 1.0 (no planing, no hump, no trim with speed); the yaw damping is chosen, not
identified (the real boat is logged turning at about 4 deg/s in missions); with the power cap the model reaches 2.61 m/s at
full thrust where the logs say 2.40 m/s; the steering limit of 35 deg and the lags cannot make the 3-6 s delay before the
yaw responds nor the turns from rest (Y16). To raise a guess: tape and scale at the quay, free decays, full-throttle runs,
coast-downs, steady turning circles, bollard pull, a steering sweep (X-9).

## Values of the hull loads that are not from a document (owner, 2026-10-12)

The numbers the model file chooses and the earlier data file did not carry were, before the placeholders, silently the
Otter's. They are Grethe's own, each labelled in the table above: `max_speed` 2.40 m/s through water is **fitted** from her
2026 control-station logs (A-74 section 1; replace by the full-throttle runs of the sea trials, X-9 G14); `max_forward_thrust`
0 N is a **decision** (A-74 Q1, owner E-135); roll damping 0.05 is **fitted** from her logs (A-74 section 3); heave 0.3 and
pitch 0.5 are an **approved guess** of 2026-10-12 (owner E-135), no source, until the quay free decays (X-9 G12). No August
2026 field data is used (safeguard rule 20).
