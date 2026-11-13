# Marie (MGM LAUV_Marie)

A named part of `more_dynamics::SubmergedBemAuv` (fully submerged hull with tabulated hydrodynamics), the research AUV Marie
(NTNU LAUV). **Every number is MGM's**, not MSS's (owner rule 2026-11-13: the check chain is MSS -> MGM -> RPP; this part
matches MGM): `more_generic_models/config/vehicles/LAUV_Marie/{plant,actuators}.yaml` and `config/globals.yaml` at MGM
`develop` e86985cd (read only). The table below gives each parameter, its value and the MGM file and line; a script read the
saved part with rpp's own parameter loader and asserted each row equal to the MGM value (89 rows, all equal).

Made with the orchestrator's `rpp ws component ...` commands (create, `subcomponent create`, `subcomponent assign` for the links,
`parameter set`) in a registered copy of the library; the links were then named as the REMUS 100 part names them ("Rudder fin 1" ...).
Composition, as the REMUS 100 part: slots `hydrostatics` (`SubmergedRestoring`), `hydrodynamics` (`BemHullLoads`), `site`
(`Site`) own, and `actuators` five links in order: `Marie rudder fin` twice (top and keel fin), `Marie stern plane fin` twice
(starboard and port), `Marie propulsor`. Inputs: the current (3), four fin deflections (rad: rudder, rudder, stern plane,
stern plane) and the shaft speed (rpm). The two fins of a pair take the same value to give MGM's two ganged commands.

**The links are the form rpp_py #5 is about**: two links to one fin part in one slot cannot be built by the runtime builder
until #5 is fixed. This part is not built through the builder yet; the checks below build the same plugins from this part's
saved values by hand (`ComponentContext`, as `tests/vehicle_models/vehicle_contract.py` does).

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
| rudder fin | `fin_position` | [-0.83, 0.0, 0.0] | `(fin_x, fin_radius cos(roll), fin_radius sin(roll))` of fin 0 (`lifting_fins.py` `_position`) | `actuators.yaml:46`, `actuators.yaml:45` |
| rudder fin | `lift_axis` | [0.0, -1.0, 0.0] | `deflection_sign * R_x(fin_roll_angle) z` of fin 0: MGM has the sign in the linkage, RPP in the axis | `actuators.yaml:43`, `actuators.yaml:44` |
| rudder fin | `fin_area` | 0.0065 | `fin_area[0]` | `actuators.yaml:47` |
| rudder fin | `fin_span` | 0.0875 | `fin_span[0]` | `actuators.yaml:48` |
| rudder fin | `aspect_ratio_factor` | 2.0 | `aspect_ratio_factor[0]` | `actuators.yaml:49` |
| rudder fin | `section_efficiency` | 0.9 | `section_efficiency[0]` | `actuators.yaml:50` |
| rudder fin | `max_deflection` | 0.4363323129985824 | `delta_max` 25 deg (to rad) | `actuators.yaml:58` |
| rudder fin | `inflow` | "rigid_point" | `inflow` | `actuators.yaml:55` |
| rudder fin | `flow_angle` | "small_angle" | `flow_angle` | `actuators.yaml:56` |
| rudder fin | `interference` | "none" | `interference` | `actuators.yaml:57` |
| rudder fin | `lift_law` | "finite_wing" | `lift_law` | `actuators.yaml:51` |
| rudder fin | `lift_slope_form` | "low_aspect_ratio" | `lift_slope_form` | `actuators.yaml:52` |
| rudder fin | `section` | "linear_section" | `drag_law none`, `stall_law none`: lift only | code / header |
| rudder fin | `zero_lift_drag` | 0.0 | `drag_law none` | code / header |
| rudder fin | `servo_dynamics` | "none" | header `actuators.yaml:20`: "fins follow their command (no servo state)" | code / header |
| rudder fin | `servo_rate_limit` | false | header `actuators.yaml:20`: the rate limit sits on the command in `control.yaml` | code / header |
| rudder fin | `servo_angle_limit` | "on_command" | range `delta_max` | code / header |
| rudder fin | `chord_axis` | [1.0, 0.0, 0.0] | `lifting_fins.py`: chord = BODY x | code / header |
| stern plane fin | `fin_position` | [-0.83, 0.0, 0.0] | `(fin_x, fin_radius cos(roll), fin_radius sin(roll))` of fin 2 (`lifting_fins.py` `_position`) | `actuators.yaml:46`, `actuators.yaml:45` |
| stern plane fin | `lift_axis` | [0.0, 0.0, -1.0] | `deflection_sign * R_x(fin_roll_angle) z` of fin 2: MGM has the sign in the linkage, RPP in the axis | `actuators.yaml:43`, `actuators.yaml:44` |
| stern plane fin | `fin_area` | 0.0065 | `fin_area[2]` | `actuators.yaml:47` |
| stern plane fin | `fin_span` | 0.0875 | `fin_span[2]` | `actuators.yaml:48` |
| stern plane fin | `aspect_ratio_factor` | 2.0 | `aspect_ratio_factor[2]` | `actuators.yaml:49` |
| stern plane fin | `section_efficiency` | 0.9 | `section_efficiency[2]` | `actuators.yaml:50` |
| stern plane fin | `max_deflection` | 0.4363323129985824 | `delta_max` 25 deg (to rad) | `actuators.yaml:58` |
| stern plane fin | `inflow` | "rigid_point" | `inflow` | `actuators.yaml:55` |
| stern plane fin | `flow_angle` | "small_angle" | `flow_angle` | `actuators.yaml:56` |
| stern plane fin | `interference` | "none" | `interference` | `actuators.yaml:57` |
| stern plane fin | `lift_law` | "finite_wing" | `lift_law` | `actuators.yaml:51` |
| stern plane fin | `lift_slope_form` | "low_aspect_ratio" | `lift_slope_form` | `actuators.yaml:52` |
| stern plane fin | `section` | "linear_section" | `drag_law none`, `stall_law none`: lift only | code / header |
| stern plane fin | `zero_lift_drag` | 0.0 | `drag_law none` | code / header |
| stern plane fin | `servo_dynamics` | "none" | header `actuators.yaml:20`: "fins follow their command (no servo state)" | code / header |
| stern plane fin | `servo_rate_limit` | false | header `actuators.yaml:20`: the rate limit sits on the command in `control.yaml` | code / header |
| stern plane fin | `servo_angle_limit` | "on_command" | range `delta_max` | code / header |
| stern plane fin | `chord_axis` | [1.0, 0.0, 0.0] | `lifting_fins.py`: chord = BODY x | code / header |
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
`Propulsor` built by hand from this part's saved values. MGM: `World` of the scenario `marie_line_following`, `SubmergedBemPlant` +
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

## Gaps: what MGM has for Marie that the plugin types here cannot represent (nothing worked around)

* **Hull lift stall** (`plant.yaml:37-39`, `lift_stall_law` smoothstep, angle and width 20 deg): `BemHullLoads` has no stall
  law, so beyond about 20 deg angle of attack the hull lift and induced drag differ from MGM (MGM blends to the axial drag).
  Not seen in the steady surge check (angle of attack zero).
* **Seabed contact** (`plant.yaml:40` `contact_station_count` 25, bed constants in the environment files): no counterpart in
  `SubmergedBemAuv` or the force producers.
* **Ganged fin commands**: MGM takes two commands (rudder, stern plane) for four fins, with the linkage sign in
  `deflection_sign`; the RPP part takes four deflections, the sign is folded into `lift_axis`, and the caller must give a pair the
  same value. Same as the REMUS 100 part.
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

* The 0.2 orchestrator drops a link assigned into an empty list slot: `create_linked_subcomponent_folder`
  (`rpp_orchestrator/component_storage.py`) takes `parent_record.subcomponents.get(slot_name, [])` and appends to that
  fresh list, which is never stored back, so the parent's `description.json` has no `actuators` although the link folders are
  written. To get the well-formed layout of the REMUS part a placeholder fin was created in the slot first, the five links
  assigned, and the placeholder removed.
* `rpp ws component rename` refuses a linked subcomponent ("A linked subcomponent cannot be renamed"); the link names
  "Rudder fin 1/2", "Stern plane fin 1/2", "Propulsor" were written into the link `description.json` files, as in the REMUS part.
