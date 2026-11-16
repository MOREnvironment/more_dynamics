# Marie (MGM LAUV_Marie)

A named part of `more_dynamics::SubmergedBemAuv` (fully submerged hull with tabulated hydrodynamics), the research AUV Marie
(NTNU LAUV). **Every number is MGM's**, not MSS's (owner rule 2026-11-13: the check chain is MSS -> MGM -> RPP; this part
matches MGM): `more_generic_models/config/vehicles/LAUV_Marie/{plant,actuators}.yaml` and `config/globals.yaml` at MGM
`develop` e86985cd (read only). The table below gives each parameter, its value and the MGM file and line; a script read the
saved part with rpp's own parameter loader and asserted each row equal to the MGM value (89 rows, all equal).

Made with the orchestrator's `rpp ws component ...` commands (`subcomponent remove`, `subcomponent create`, `subcomponent assign`) in a registered copy of the library; each fin's `params/parameters.py` was then written by hand (see below).
Composition: slots `hydrostatics` (`SubmergedRestoring`), `hydrodynamics` (`BemHullLoads`), `site` (`Site`) own, and `actuators`
five children in order: four **own** `Fin` components `Top fin`, `Keel fin`, `Starboard fin`, `Port fin` (MGM fin index 0, 1, 2, 3),
then the link `Propulsor` to the saved `Marie propulsor`. Inputs: the current (3), then one command per actuator in list order
(`fin_angle_command` x4, `drive_command`): the four fin deflections (rad) in the order top, keel, starboard, port, and the shaft
speed (rpm). The vehicle gives each child its command by list position only (`vehicle_graph.py`, "the vehicle's inputs, in list
order"), so the order of the fins in the slot is what sends the rudder command to top and keel (inputs 4, 5) and the stern-plane
command to starboard and port (inputs 6, 7), as MGM's `fin_command: [0, 0, 1, 1]` (`actuators.yaml:42`) does. The caller gives
a pair the same value to get MGM's two ganged commands.

**No link twice** (RPP decision, Luka): a link is an immutable pointer that carries the target's parameters and plugin type, so two
links to one component in one slot are not allowed. Each physical fin is its own component with its own `params/parameters.py`
(its `fin_position` and `lift_axis`). The values the four fins share are in one file, `data/fins_common.yml` of this part
(next to this file), which every fin's `parameters.py` reads with
`Path(__file__).resolve().parents[3] / "data" / "fins_common.yml"` (file -> `params` -> fin folder -> `subcomponents` -> this part)
and `yaml.safe_load`. `parameters.py` is a Python script on purpose. **Saving from the editor or `parameter set` rewrites a
`parameters.py` as plain values** (`component_storage.py`, `_build_component_parameters_source`): the yml read and the comments are
then gone and the values stay as they were; each fin's file says so in its first comment. Duplicating or re-saving a fin outside
this part breaks the relative path. Both loaders keep `__file__`: rpp_py `ParameterHandler` (`spec_from_file_location`) and the
orchestrator `ComponentParameterStore.load` (develop: `spec_from_file_location`; PR 18 branch: `exec` with `__file__` set), checked
on both.

The saved fin parts `Marie rudder fin` and `Marie stern plane fin` that this part used to link twice are removed (nothing else linked them).
This part builds through `ComponentRuntimeBuilder` (checks below); no hand assembly is needed.

## Relation to `tests/vehicle_models/data/marie/marie_parameters.json` ("Marie (estimated)")

The test data is the same vehicle written from an earlier MGM state. Compared value by value with MGM today: the hull
(`vehicle`), hydrostatics, hull loads and propulsor groups are **equal** to MGM. Three groups differ, and the part follows MGM:

| Test data | MGM today (this part) | Why |
|---|---|---|
| fins: `fin_position` x = -1.075 (tail-end station), `fin_area` 0.00665, `fin_span` 0.0857 | x = -0.83, area 0.0065, span 0.0875 | the test data is Prestero's literature fin set (MGM `actuators_literature_fins.yaml`); MGM's `actuators.yaml` now has the fins measured with a ruler (`actuators.yaml` header) |
| `site_form` "given", `gravity` 9.81 | `latitude` 63.446827 deg (gravity 9.8217848) | MGM derives gravity from the scenario latitude (`scenario.py:213-215`) |
| `smooth_speed`, `smooth_speed_epsilon` in the test data | not in MGM | RPP numerical guard, left at the plugin default |

The test data therefore does not equal MGM's current Marie in the fins and the site; this part is the one to use for MGM parity.

## Parameter table (value, MGM source, file:line)

Files: `plant.yaml`, `actuators.yaml` in `more_generic_models/config/vehicles/LAUV_Marie/`, `globals.yaml` in `more_generic_models/config/`.
Long lists are summarised (count, first, last); the saved part holds them in full and the check compared every entry.

| Part | Parameter | RPP value | MGM source | MGM file:line |
|---|---|---|---|---|
| hull | `frequencies` | 38 values, first 0.0, last 10.0 | `frequencies` | `plant.yaml:62` |
| hull | `added_mass` | 1368 values, first 0.6582080619327014, last 9.76285293551042 | `added_mass` | `plant.yaml:69` |
| hull | `potential_damping` | 1368 values, first 0.0, last 0.0 | `potential_damping` | `plant.yaml:297` |
| hull | `hydrostatic_stiffness_345` | 9 values, first 0.0, last 3.3406509596954157 | `hydrostatic_stiffness_345` | `plant.yaml:525` |
| hull | `database_reference` | [-0.0, 0.0, 0.01] | `database_reference` | `plant.yaml:527` |
| hull | `length` | 2.15 | `length` | `plant.yaml:533` |
| hull | `beam` | 0.14878654394686797 | `beam` | `plant.yaml:534` |
| hull | `height` | 0.15 | `height` | `plant.yaml:535` |
| hull | `omega_p` | 2.7 | `omega_p` | `plant.yaml:28` |
| hull | `spectrum_gamma` | 3.3 | `spectrum_gamma` | `plant.yaml:29` |
| hull | `coriolis_form` | "full" | `coriolis_form` | `plant.yaml:43` |
| hull | `body_mass` | 34.0 | `mass` | `plant.yaml:528` |
| hull | `body_center_of_gravity` | [-0.0, 0.0, 0.01] | `r_cg` | `plant.yaml:530` |
| hull | `inertia_diagonal` | [0.098906, 13.148998160000001, 13.145598160000002] | `mass * radii_of_gyration^2` about the CG, moved to the CO by the parallel axis with `r_cg` (what `M_RB_co` does, `submerged_bem.py` build) | `plant.yaml:529` |
| hull | `inertia_products` | [0.0, 0.0, 0.0] | products of inertia of the diagonal CG inertia moved by `r_cg = (0, 0, 0.01)`: zero | `plant.yaml:530` |
| hull | `kinematics_form` | "euler" | MGM states are NED Euler angles (`rigid_body_6dof/base.py`) | `base.py` (code, no yaml key) |
| hull | `current_form` | "full_attitude" | MGM `CURRENT_FORM = "submerged"`: current rotated by the full attitude (eq. 8.157) | `submerged_bem.py` docstring (code) |
| hull | `integration_max_step` | 0.05 | RPP numerics only (the plugin default); MGM integrates at the runner period | none |
| hydrostatics | `buoyancy_method` | "immersed_body_of_revolution" | `restoring_law` = immersed_body_of_revolution and `neutral_buoyancy` false | `plant.yaml:44`, `:42` |
| hydrostatics | `center_of_buoyancy` | [-1.831323617232447e-16, -4.432096739586757e-19, 4.1081372303963803e-07] | `r_cb` | `plant.yaml:531` |
| hydrostatics | `displaced_volume` | 0.033359804878048786 | `displaced_volume` | `plant.yaml:532` |
| hydrostatics | `hull_profile_x` | 121 values, first -1.075, last 1.075 | `hull_profile_x` | `plant.yaml:542` |
| hydrostatics | `hull_profile_radius` | 121 values, first 0.001, last 0.001 | `hull_profile_radius` | `plant.yaml:563` |
| hydrostatics | `hull_axis_origin` | [0.0, 0.0, 0.0] | `hull_axis_origin` | `plant.yaml:45` |
| hydrostatics | `hull_profile_slices` | 4 | `hull_profile_slices` | `plant.yaml:46` |
| hydrostatics | `hull_profile_cb_tolerance` | 0.001 | `hull_profile_cb_tolerance` | `plant.yaml:47` |
| hull loads | `T_1236` | [50.0, 5.0, 5.0, 5.0] | `T_1236` | `plant.yaml:30` |
| hull loads | `delta_zeta_45` | [0.2, 0.2] | `delta_zeta_45` | `plant.yaml:31` |
| hull loads | `kappa_456` | [5.0, 5.0, 0.0] | `kappa_456` | `plant.yaml:32` |
| hull loads | `surge_damping_blend` | 3.0 | `surge_damping_blend` | `plant.yaml:33` |
| hull loads | `p_factor` | 0.7 | `p_factor` | `plant.yaml:34` |
| hull loads | `c_d` | 0.42 | `c_d` | `plant.yaml:35` |
| hull loads | `e_oswald` | 0.3 | `e_oswald` | `plant.yaml:36` |
| hull loads | `cross_flow_model` | "cylinder" | `cross_flow_model` | `plant.yaml:41` |
| site | `latitude` | 1.1073560310932362 | `mu` 63.446827 deg (to rad); gravity is the WGS-84 normal gravity at it in both (`scenario.py:213-215`) | `globals.yaml:8` |
| site | `water_density` | 1025.0 | `rho` | `globals.yaml:6` |
| site | `kinematic_viscosity` | 1e-06 | literal `1e-6` of the cross-flow Reynolds number (no yaml key) | `dynamics/plant/matrices/hydrodynamics.py:221` |
| site | `site_form` | "latitude" | gravity from the latitude, as MGM `resolve_globals` | `scenario.py:213-215` |
| Top fin | `fin_position` | [-0.83, 0.0, 0.0] | `(fin_x, fin_radius cos(roll), fin_radius sin(roll))` of fin 0; `fin_radius` is 0.0 for every fin, a placeholder (hull diameter at the fins not measured, register item M-0679), kept 0 so the comparison with MGM stays one-to-one | `actuators.yaml:46`, `:45`, `:43`; `lifting_fins.py:92` |
| Top fin | `lift_axis` | [0.0, -1.0, 0.0] | `deflection_sign * R_x(fin_roll_angle) z` = `(-1) (0, -sin 270 deg, cos 270 deg)` of fin 0: MGM has the sign in the linkage, RPP in the axis | `actuators.yaml:43`, `:44`; `lifting_fins.py:78`, `:179` |
| Top fin | command input | input 4 of the vehicle (rudder command) | `fin_command[0]` = 0; list position in the slot | `actuators.yaml:42`; `lifting_fins.py:60`, `:179` |
| Keel fin | `fin_position` | [-0.83, 0.0, 0.0] | `(fin_x, fin_radius cos(roll), fin_radius sin(roll))` of fin 1; `fin_radius` is 0.0 for every fin, a placeholder (hull diameter at the fins not measured, register item M-0679), kept 0 so the comparison with MGM stays one-to-one | `actuators.yaml:46`, `:45`, `:43`; `lifting_fins.py:92` |
| Keel fin | `lift_axis` | [0.0, -1.0, 0.0] | `deflection_sign * R_x(fin_roll_angle) z` = `(+1) (0, -sin 90 deg, cos 90 deg)` of fin 1: MGM has the sign in the linkage, RPP in the axis | `actuators.yaml:43`, `:44`; `lifting_fins.py:78`, `:179` |
| Keel fin | command input | input 5 of the vehicle (rudder command) | `fin_command[1]` = 0; list position in the slot | `actuators.yaml:42`; `lifting_fins.py:60`, `:179` |
| Starboard fin | `fin_position` | [-0.83, 0.0, 0.0] | `(fin_x, fin_radius cos(roll), fin_radius sin(roll))` of fin 2; `fin_radius` is 0.0 for every fin, a placeholder (hull diameter at the fins not measured, register item M-0679), kept 0 so the comparison with MGM stays one-to-one | `actuators.yaml:46`, `:45`, `:43`; `lifting_fins.py:92` |
| Starboard fin | `lift_axis` | [0.0, 0.0, -1.0] | `deflection_sign * R_x(fin_roll_angle) z` = `(-1) (0, -sin 0 deg, cos 0 deg)` of fin 2: MGM has the sign in the linkage, RPP in the axis | `actuators.yaml:43`, `:44`; `lifting_fins.py:78`, `:179` |
| Starboard fin | command input | input 6 of the vehicle (stern plane command) | `fin_command[2]` = 1; list position in the slot | `actuators.yaml:42`; `lifting_fins.py:60`, `:179` |
| Port fin | `fin_position` | [-0.83, 0.0, 0.0] | `(fin_x, fin_radius cos(roll), fin_radius sin(roll))` of fin 3; `fin_radius` is 0.0 for every fin, a placeholder (hull diameter at the fins not measured, register item M-0679), kept 0 so the comparison with MGM stays one-to-one | `actuators.yaml:46`, `:45`, `:43`; `lifting_fins.py:92` |
| Port fin | `lift_axis` | [0.0, 0.0, -1.0] | `deflection_sign * R_x(fin_roll_angle) z` = `(+1) (0, -sin 180 deg, cos 180 deg)` of fin 3: MGM has the sign in the linkage, RPP in the axis | `actuators.yaml:43`, `:44`; `lifting_fins.py:78`, `:179` |
| Port fin | command input | input 7 of the vehicle (stern plane command) | `fin_command[3]` = 1; list position in the slot | `actuators.yaml:42`; `lifting_fins.py:60`, `:179` |
| all four fins (`data/fins_common.yml`) | `fin_area` | 0.0065 | `fin_area` (every entry) | `actuators.yaml:47` |
| all four fins (`data/fins_common.yml`) | `fin_span` | 0.0875 | `fin_span` | `actuators.yaml:48` |
| all four fins (`data/fins_common.yml`) | `aspect_ratio_factor` | 2.0 | `aspect_ratio_factor` | `actuators.yaml:49` |
| all four fins (`data/fins_common.yml`) | `section_efficiency` | 0.9 | `section_efficiency` | `actuators.yaml:50` |
| all four fins (`data/fins_common.yml`) | `max_deflection` | 0.4363323129985824 | `delta_max` 25 deg (to rad) | `actuators.yaml:58` |
| all four fins (`data/fins_common.yml`) | `inflow` | "rigid_point" | `inflow` | `actuators.yaml:55` |
| all four fins (`data/fins_common.yml`) | `flow_angle` | "small_angle" | `flow_angle` | `actuators.yaml:56` |
| all four fins (`data/fins_common.yml`) | `interference` | "none" | `interference` | `actuators.yaml:57` |
| all four fins (`data/fins_common.yml`) | `lift_law` | "finite_wing" | `lift_law` | `actuators.yaml:51` |
| all four fins (`data/fins_common.yml`) | `lift_slope_form` | "low_aspect_ratio" | `lift_slope_form` | `actuators.yaml:52` |
| all four fins (`data/fins_common.yml`) | `section` | "linear_section" | `drag_law none`, `stall_law none`: lift only | `actuators.yaml:53-54` |
| all four fins (`data/fins_common.yml`) | `zero_lift_drag` | 0.0 | `drag_law none` | `actuators.yaml:53` |
| all four fins (`data/fins_common.yml`) | `servo_dynamics` | "none" | header: "fins follow their command (no servo state)" | `actuators.yaml:20` |
| all four fins (`data/fins_common.yml`) | `servo_rate_limit` | false | header: the rate limit sits on the command in `control.yaml` | `actuators.yaml:20` |
| all four fins (`data/fins_common.yml`) | `servo_angle_limit` | "on_command" | range `delta_max` | `actuators.yaml:58` |
| all four fins (`data/fins_common.yml`) | `chord_axis` | [1.0, 0.0, 0.0] | chord = BODY x | `lifting_fins.py:76` |
| propulsor | `power_law_coefficient` | 0.00842703 | `k_thrust` | `actuators.yaml:65` |
| propulsor | `power_law_exponent` | 2.2633 | `p_thrust` | `actuators.yaml:66` |
| propulsor | `astern_thrust_ratio` | 1.0 | `astern_ratio` | `actuators.yaml:67` |
| propulsor | `reaction_torque_per_thrust` | 0.002852 | `reaction_torque_per_thrust` | `actuators.yaml:68` |
| propulsor | `max_command` | 2340.0 | `n_max` (rpm) | `actuators.yaml:64` |
| propulsor | `min_command` | -2340.0 | `-n_max` (the command is ranged to +-n_max) | `actuators.yaml:64` |
| propulsor | `position` | [0.0, 0.0, 0.0] | default (neutral: at CO); `actuators.yaml` does not set it | `power_law/params.py` default `position` |
| propulsor | `direction` | [1.0, 0.0, 0.0] | shaft along +x | `power_law/params.py` default `orientation_rpy` 0 |
| propulsor | `shaft_torque_scale` | 1.0 | full reaction torque | `power_law.py` wrench: `reaction = -shaft_sense c X` |
| propulsor | `drive` | "speed_command_rpm" | selector: MGM `PowerLawPropeller` is a fixed shaft, rpm command, `X = k|n|^p`, no inflow term | `power_law.py` docstring |
| propulsor | `steering` | "fixed" | selector: MGM `PowerLawPropeller` is a fixed shaft, rpm command, `X = k|n|^p`, no inflow term | `power_law.py` docstring |
| propulsor | `mounting` | "fixed_direction" | selector: MGM `PowerLawPropeller` is a fixed shaft, rpm command, `X = k|n|^p`, no inflow term | `power_law.py` docstring |
| propulsor | `inflow` | "none" | selector: MGM `PowerLawPropeller` is a fixed shaft, rpm command, `X = k|n|^p`, no inflow term | `power_law.py` docstring |
| propulsor | `law` | "power_law" | selector: MGM `PowerLawPropeller` is a fixed shaft, rpm command, `X = k|n|^p`, no inflow term | `power_law.py` docstring |
| propulsor | `housing` | "none" | selector: MGM `PowerLawPropeller` is a fixed shaft, rpm command, `X = k|n|^p`, no inflow term | `power_law.py` docstring |

Kinds (from the MGM headers): length 2.15 m and mass 34 kg are manufacturer values; the frequency tables, section dimensions,
centres and radii of gyration are an idealised bare hull (estimate); `displaced_volume` is INTERIM (weight with 0.57 % reserve
buoyancy, +1.9 N); `power_law_coefficient` and `power_law_exponent` are fitted to field data together with MGM's hull
resistance and are not a property of the screw alone (`actuators.yaml` header); the fins are ruler measurements (+-10-20 %)
with three placeholders (`aspect_ratio_factor`, `section_efficiency`, `fin_radius`).

## Check 2: still-water steady surge speed, RPP plugins against MGM

Hull level at 20 m depth, fins at zero, still water; the speed at which the surge acceleration of the full 6-DOF dynamics is
zero (brentq), at the shaft speeds below. RPP: `SubmergedBemAuv` + `SubmergedRestoring` + `BemHullLoads` + `Site` + 4 `Fin` +
`Propulsor`, now built from this part's folder by `ComponentRuntimeBuilder` (four own fins, `fins_common.yml` read). MGM: `World` of the scenario `marie_line_following`, `SubmergedBemPlant` +
`LiftingFins` + `PowerLawPropeller` (`vehicles/LAUV_Marie`), the same method as MGM's `test_power_law_propeller.py`.

| Hull | rpm | RPP u (m/s) | MGM u (m/s) | difference |
|---|---|---|---|---|
| neutral buoyancy (`from_volume` m/rho; MGM `neutral_buoyancy` true, `fully_submerged`) | 1118 | 1.287172975584607 | 1.287172975584607 | 0 |
| neutral buoyancy | 1500 | 1.7973314059424663 | 1.7973314059424663 | 0 |
| as built (floating law, +1.9 N reserve) | 1118 | 1.287179528395304 | 1.287179528395304 | 0 |
| as built | 1500 | 1.7973360834487038 | 1.7973360834487038 | 0 |

(1118 rpm is the shaft speed of the logs, 1.286-1.299 m/s logged.) Also equal to the last digit: the wrench of the rudder pair
(+0.1 rad, 1.5 m/s: Y -4.805332457512563 N, N 3.988425939735427 N m), the stern-plane pair (Z -4.805332457512563 N, M
-3.988425939735427 N m) and the propulsor at 2000 rpm (X 23.572519528271435 N, K -0.06722882569463014 N m).

Re-run after the split into four own fins (part built by the runtime builder, 2026-11-16): the four fin wrenches summed equal MGM's
`LiftingFins` wrench in 18 cases (three velocities including sway, heave and rates, six command pairs of both signs), worst difference
8.9e-16 (largest entry 14.4); each fin alone: a rudder command gives a force only on top and keel, a stern-plane command only on starboard
and port, each half of MGM's pair value. Steady surge speed: neutral buoyancy 1118 rpm 1.287172975584607, 1500 rpm 1.7973314059424694
(MGM 1.287172975584607 and 1.7973314059424663, difference 3e-15); as built 1.287179528395304 and 1.7973360834487067 (MGM 1.287179528395304 and 1.7973360834487038, 3e-15). The MGM speeds are the values recorded above, not
recomputed in this re-run; the MGM fin wrenches were recomputed.

## Gaps: what MGM has for Marie that the plugin types here cannot represent (nothing worked around)

* **Hull lift stall** (`plant.yaml:37-39`, `lift_stall_law` smoothstep, angle and width 20 deg): `BemHullLoads` has no stall
  law, so beyond about 20 deg angle of attack the hull lift and induced drag differ from MGM (MGM blends to the axial drag).
  Not seen in the steady surge check (angle of attack zero).
* **Seabed contact** (`plant.yaml:40` `contact_station_count` 25, bed constants in the environment files): no counterpart in
  `SubmergedBemAuv` or the force producers.
* **Ganged fin commands**: MGM takes two commands (rudder, stern plane) for four fins, with the linkage sign in
  `deflection_sign`; the RPP part takes four deflections (one per own fin component), the sign is folded into each fin's `lift_axis`, and the caller must
  give a pair the same value.
* **Fin command-rate limit** 300 deg/s: in MGM it acts on the command in `control.yaml` (`delta_r_dot_max`, `delta_s_dot_max`);
  the RPP fins have none (`servo_rate_limit` false), as in MGM's actuator.
* **Shaft-rate channel** (`actuators.yaml:63`, `propeller_rate`, read by the `PropellerRateSpeed` sensor): the RPP `Propulsor`
  publishes no such channel.
* **Everything else in the vehicle folder is not dynamics and has no plugin type in `more_dynamics`**: sensors (`sensors.yaml`:
  RateGyro, Accelerometer, MagneticCompass, GaussMarkovGnss with the antenna lever arm, BottomTrackLog DVL, HydrostaticDepth,
  AttitudeHeadingReference, PropellerRateSpeed), navigation (EKF, aided inertial), aiding, control, allocation, guidance,
  supervision, communication, payload and the vehicle API. `VehicleModel3D` holds the hull, its restoring and loads, the
  actuators and the site only.
* **Initial state**: MGM starts a floating hull at its sea equilibrium (`initial.start`); this part holds no initial state.

## Findings about the tooling (not worked around here)

* `rpp ws component rename` refuses a linked subcomponent ("A linked subcomponent cannot be renamed"); the link name "Propulsor" was
  written into the link `description.json`. The slot was emptied (`subcomponent remove` of all five links) and refilled in order, so
  the fins come before the propulsor in the list; `subcomponent create` into the emptied list slot stored the children in the
  parent (the earlier problem of a link assigned into an empty list slot did not occur).
* A fin created by the orchestrator carries all plugin parameters, with empty lists for the unset ones; each fin's `parameters.py`
  was written by hand with the same set of names as the earlier fin parts.
