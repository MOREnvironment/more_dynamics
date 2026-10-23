# Defaults

**Plugins carry no vehicle numbers** (owner, 2026-10-12): every value of a vehicle, a device or a site is declared with
`unset()` (an empty list) in the plugin and comes from a named part of the workspace (`.rppws/parts/`) when the
composition is built; a value the selected forms read and the part leaves out makes the build fail, naming it with its unit
and meaning. The only defaults a plugin keeps are text **selectors** (a form or method, with its physics default), on/off
**switches**, and three **numerical settings** that are no vehicle's data: `integration_max_step` (0.05 s, the RK4
sub-step), `smooth_speed_epsilon` (1e-09 m/s, regularisation of the speed guard) and `ittc_reynolds_floor` (1e5, the lower
bound of the ITTC friction line). Each type's out-of-the-box example is a reference named part.

This file records the **provenance of each reference part's values**: one table per part (parameter, value, unit, source, kind).
Kinds: `document` (read off a cited page or file), `estimate` (a starting point with no document behind it, to be
identified), `derived` (computed from another primitive by a stated equation), `selector` (chooses a form), `convention`
(a modelling choice), `fitted` (identified from data, with its residual), `owner decision`, `approved guess` (an owner-approved
placeholder, no source, until a measurement replaces it) and `gap` (a zero or value stated for want of data). A file citation
on a code line (`remus100.m:208`) is not a vehicle name under the everything-generic rule: it stays where the equation is.

The value tables below are generated from the parts' own `params/parameters.py` (the values) and from the earlier per-plugin
tables and the comment beside each value (source and kind); a part's `SOURCE.md` says more.

## REMUS 100 (MSS) (`more_dynamics::SpheroidAuv`)

| Parameter | Value | Unit | Source | Kind |
|---|---|---|---|---|
| `integration_max_step` | 0.05 |  | SIMremus100.m 60 (h = 0.05 s), the step of the MSS vehicle simulations | setting |
| `length` | 1.6 | m | remus100.m:132 L_auv, passed as L at :222 | document |
| `diameter` | 0.19 | m | remus100.m:133 D_auv, passed as B at :222 | document |
| `body_density` | 1054.7872613500267 | kg/m^3 | mass 31.9 kg (remus100.m:3, given) / (4/3 pi a b^2); MSS: 1025 at spheroid.m:35-36 | derived |
| `body_center_of_gravity` | [0, 0, 0.02] | m | remus100.m:138 r_bG | document |
| `roll_added_inertia_ratio` | 0.3 | 1 | remus100.m:137 r44 | document |
| `added_mass_form` | "lamb_spheroid" | - | imlay61.m 31-60 | selector |
| `coriolis_form` | "munk_couplings_removed" | - | remus100.m 208-211 (MSS omits the Munk-moment couplings of the added-mass Coriolis: its reduced model keeps no measured manoeuvring derivatives to balance them; with every term kept this part spins open-loop, SOURCE.md) | selector |
| `current_form` | "none" | - | still water | selector |
| `body_mass` | 31.9 | kg | remus100.m:3 | document |

## REMUS 100 hydrostatics (`more_dynamics::SubmergedRestoring`)

| Parameter | Value | Unit | Source | Kind |
|---|---|---|---|---|
| `buoyancy_method` | "neutral" | - | remus100.m:215 B = W | selector |
| `center_of_buoyancy` | [0, 0, 0] | m | remus100.m:139 r_bB | document |
| `displaced_volume` | 0.030243065278557742 | m^3 | 4/3 pi (L/2)(D/2)^2, remus100.m:132-133 | derived |
| `buoyancy` | 313.31493718007266 | N | m g with m = 31.9 kg (remus100.m:3), g at remus100.m:96-97 | derived |

## REMUS 100 fin pairs (`more_dynamics::FinPairsDeflectionOnly`)

| Parameter | Value | Unit | Source | Kind |
|---|---|---|---|---|
| `rudder_area` | 0.0133 | m^2 | remus100.m:184 A_r = 2 * S_fin (S_fin :180) | document |
| `stern_plane_area` | 0.0133 | m^2 | remus100.m:189 A_s = 2 * S_fin (S_fin :180) | document |
| `rudder_lift_coefficient` | 0.5 | 1/rad | remus100.m:183 CL_delta_r | document |
| `stern_plane_lift_coefficient` | 0.7 | 1/rad | remus100.m:188 CL_delta_s | document |
| `rudder_position` | -0.8 | m | remus100.m:185 x_r = -a, a = L_auv / 2 (one geometry) | derived |
| `stern_plane_position` | -0.8 | m | remus100.m:190 x_s = -a, a = L_auv / 2 (one geometry) | derived |
| `max_deflection` | 0.3490658503988659 | rad | remus100.m:110 delta_max = deg2rad(20) | document |

## REMUS 100 hull loads (`more_dynamics::AuvHullLoads`)

| Parameter | Value | Unit | Source | Kind |
|---|---|---|---|---|
| `time_constants` | [20, 20, 1] | s | remus100.m:193, 194, 197 [T1 T2 T6] (Dmtrx.m call at :218) | document |
| `damping_ratios` | [0.3, 0.8] | 1 | remus100.m:195, 196 [zeta4 zeta5] | document |
| `planform_area` | 0.2128 | m^2 | remus100.m:134 S = 0.7 * L_auv * D_auv | derived |
| `parasitic_drag_coefficient` | 0.05595961914206819 | 1 | remus100.m:144-145 CD_0 = Cd * pi * b^2 / S, b = D_auv / 2 (one geometry) | derived |
| `oswald_efficiency` | 0.3 | 1 | coeffLiftDrag.m:54 e = 0.3 | document |
| `planform_fraction` | 0.7 | 1 | remus100.m:134 (S = 70% of the rectangle L_auv * D_auv) | document |
| `hull_frontal_drag_coefficient` | 0.42 | 1 | remus100.m:143 (Cd), Allen et al. (2000) as MSS cites it, not read | document |

## REMUS 100 propeller (`more_dynamics::Propeller`)

| Parameter | Value | Unit | Source | Kind |
|---|---|---|---|---|
| `propeller_diameter` | 0.14 | m | remus100.m:149 D_prop | document |
| `max_shaft_speed` | 1525 | rpm | remus100.m:111 n_max | document |
| `thrust_deduction` | 0.1 | 1 | remus100.m:150 t_prop | document |
| `wake_fraction` | 0.05600000000000005 | 1 | remus100.m:151 Va = 0.944 * U_r, w = 1 - 0.944 | derived |
| `pitch_diameter_ratio` | 1 | 1 | remus100.m:157 wageningen(0,1,0.718,3), second argument | document |
| `blade_area_ratio` | 0.718 | 1 | remus100.m:156-157 blade-area ratio 0.718 | document |
| `blade_count` | 3 | 1 | remus100.m:156-157 3 blades | document |
| `max_advance_number` | 0.6632 | 1 | remus100.m:154 Ja_max | document |
| `roll_moment_scale` | 0.1 | 1 | remus100.m:253 tau(4) = K_prop / 10 | document |
| `position` | [0, 0, 0] | m | remus100.m:250-253 thrust on x_b through the CO (no moment arm) | document |
| `orientation` | [0, 0, 0] | rad | remus100.m:250, 253 shaft along x_b (thrust in tau(1), torque in tau(4)) | document |
| `thrust_torque_coefficients` | [0.4566, 0.07, 0.1798, 0.0312] | 1 | remus100.m:158, 159, 161, 162 [KT_0 KQ_0 KT_max KQ_max] | document |

## REMUS 100 (MSS) site (`more_dynamics::Site`)

| Parameter | Value | Unit | Source | Kind |
|---|---|---|---|---|
| `site_form` | 'latitude' | - | remus100.m 96-99 (the site of the MSS simulation: gravity from the latitude) | selector |
| `latitude` | 1.1073560310932362 | rad | remus100.m 96-97 (mu = deg2rad(63.446827)) | document |
| `gravity` | 9.81 | m/s^2 | not read with site_form latitude; mssConstants.m:14 @ 49e03e3 (MSS 2.0.5) | document |
| `water_density` | 1025 | kg/m^3 | remus100.m 98-99 (rho = constant.rho_water), mssConstants.m 12 | document |
| `kinematic_viscosity` | 1e-06 | m^2/s | cylinderDrag.m 78-80; Fossen 2011, p. 125, below eq. 6.85 (20 degC) | document |

## REMUS 100 fin (`more_dynamics::Fin`)

| Parameter | Value | Unit | Source | Kind |
|---|---|---|---|---|
| `fin_position` | [-0.638, 0.0, 0.0] | m | Prestero 2001, Table A.5, p. 103 (x_finpost) | document |
| `chord_axis` | [1.0, 0.0, 0.0] | 1 | unit chord direction, BODY FRD (convention) | convention |
| `lift_axis` | [0.0, 1.0, 0.0] | 1 | unit positive-lift direction, BODY FRD (convention) | convention |
| `fin_area` | 0.00665 | m^2 | Prestero 2001, Table A.5, p. 103 (S_fin) | document |
| `time_constant` | 0.1 | s | Sarhadi 2026, Fig. 4, p. 4 (a demanding choice of that study, not a measurement) | estimate |
| `max_deflection` | 0.3490658503988659 | rad | rad = 20 deg, Sarhadi 2026, Fig. 4, p. 4 | estimate |
| `max_rate` | 0.5235987755982988 | rad/s | rad/s = 30 deg/s, Sarhadi 2026, Fig. 4, p. 4 | estimate |
| `lift_slope` | 3.12 | 1/rad | Prestero 2001, Table A.5, p. 103 (c_L_alpha) | document |
| `zero_lift_drag` | 0.0 | 1 | zero with its reason: Prestero 2001, eq. 4.37, p. 31 carries no fin drag (read by linear_section only) | document |

## Otter (MSS) (`more_dynamics::Catamaran`)

| Parameter | Value | Unit | Source | Kind |
|---|---|---|---|---|
| `integration_max_step` | 0.05 |  | SIMremus100.m 60 (h = 0.05 s), the step of the MSS vehicle simulations | setting |
| `length` | 2.0 | m | otter.m 98-99, 111-114 | document |
| `beam` | 1.08 | m | otter.m 98-99, 111-114 | document |
| `pontoon_beam` | 0.25 | m | otter.m 111-114 | document |
| `pontoon_lateral_offset` | 0.395 | m | otter.m 111-114 | document |
| `pontoon_block_coefficient` | 0.4 | 1 | otter.m 111-114 | document |
| `pontoon_waterplane_coefficient` | 0.75 | 1 | otter.m 111-114 | document |
| `hull_mass` | 55.0 | kg | otter.m 100 | document |
| `payload_mass` | 25.0 | kg | generate_catamaran_mss.m 61-62 (payload mass) | document |
| `hull_center_of_gravity` | [0.2, 0.0, -0.2] | m | otter.m 101 | document |
| `payload_position` | [0.05, 0.0, -0.35] | m | generate_catamaran_mss.m 61-62 (payload position) | document |
| `radii_of_gyration` | [0.4, 0.25, 0.25] | 1 | otter.m 102-104 | document |
| `added_mass_coefficients` | [-1.0, -1.5, -1.0, -0.2, -0.8, -1.7] | 1 | otter.m 161-166 | document |
| `added_mass_form` | "scaled_derivatives_loaded_surge" | - | otter.m 160-168, surge term from m + mP (line 162) | selector |
| `coriolis_form` | "kirchhoff_full" | - | m2c.m 33-48 (every term kept) | selector |
| `current_form` | "none" | - | still water | selector |

## Otter hull loads (`more_dynamics::SurfaceHullLoads`)

| Parameter | Value | Unit | Source | Kind |
|---|---|---|---|---|
| `max_forward_thrust` | 239.364 | N | otter.m, 24.4 kgf = 24.4 * 9.81 N | document |
| `max_speed` | 3.0864 | m/s | otter.m, 6 knots = 6 * 0.5144 m/s | document |
| `time_constants` | [0.5, 1.0] | s | otter.m 106-107, 213, 217 | document |
| `surge_time_constant` | 5.0 | s | otter.m 105, 211 | document |
| `surge_blend_decay` | 3.0 | s/m | otter.m 245 exp(-3 abs(u_r)) | document |
| `manoeuvring_damping` | "time_constants_surge_blend" | - | otter.m 206-217, 244-251 | selector |
| `damping_ratios` | [0.3, 0.2, 0.4] | 1 | otter.m 214-216 | document |
| `yaw_damping_nonlinearity` | 10.0 | s/rad | otter.m 251 | document |
| `surge_resistance` | "none" | - | plugin default (off unless a composition turns it on) | selector |
| `time_constant` | 100.0 | s | osv.m 128, vessel.T1 = 100 s (read with "ittc") | document |
| `form_factor` | 0.1 | 1 | XuuITTC.m 33, k = 0.1; Fossen 2011, p. 125 (read with "ittc") | document |
| `crossover_speed` | 2.0 | m/s | forceSurgeDamping.m 57 @ ac77394, u_cross = 2 m/s (read with "ittc") | document |
| `surge_added_mass_factor` | 2.7 | 1 | addedMassSurge.m 34, A11 = 2.7 rho nabla^(5/3) / L^2 (read with "ittc") | document |
| `ittc_reynolds_floor` | 100000.0 | 1 | XuuITTC.m 35, Re_min = 1e5 (read with "ittc") | document |

## Otter (MSS) site (`more_dynamics::Site`)

| Parameter | Value | Unit | Source | Kind |
|---|---|---|---|---|
| `site_form` | 'given' | - | mssConstants.m 12, 14; otter.m 95-97 (g = constant.g, rho = constant.rho_water) | selector |
| `latitude` | 1.1073560310932362 | rad | not read with site_form given (REMUS 100 site value) | document |
| `gravity` | 9.81 | m/s^2 | mssConstants.m:14; otter.m 96 | document |
| `water_density` | 1025.0 | kg/m^3 | mssConstants.m:12; otter.m 97 | document |
| `kinematic_viscosity` | 1e-06 | m^2/s | cylinderDrag.m 78-80; Fossen 2011, p. 125, below eq. 6.85 (20 degC) | document |

## Otter hydrostatics (`more_dynamics::SurfaceRestoring`)

| Parameter | Value | Unit | Source | Kind |
|---|---|---|---|---|
| `hull_count` | 2 | - | otter.m 181-202: two pontoons | document |
| `longitudinal_inertia_factor` | 0.8 | 1 | otter.m 188 (I_L = 0.8 ...) | document |
| `longitudinal_center_of_flotation` | -0.2 | m | otter.m 186, 201 | document |
| `reference_point` | [0.0, 0.0, 0.0] | m | the body-frame origin (CO) | convention |

## Grethe (`more_dynamics::Monohull`)

| Parameter | Value | Unit | Source | Kind |
|---|---|---|---|---|
| `length` | 5.2 | m | A-43 M1, F003 p2, M001 p7 (overall length) | document |
| `beam` | 2.15 | m | A-43 M2, F003 p2, M001 p7, D002 p1 | document |
| `hull_block_coefficient` | 0.234 | 1 | A-43 M4, 806 / (1025 * 5.2 * 2.15 * 0.3), a mass knob not hull fullness | derived |
| `hull_waterplane_coefficient` | 0.8 | 1 | A-43 M14 | estimate |
| `hull_mass` | 806.0 | kg | A-43 M7, F003 p2 (FAT, +/- 20 kg) | document |
| `payload_mass` | 0.0 | kg | a bare hull; the 2026 field payload (USBL mount, battery packs) is not weighed (A-43 M8, gap) | gap |
| `hull_center_of_gravity` | [0.0, 0.0, 0.025] | m | A-43 M9 | estimate |
| `payload_position` | [0.0, 0.0, 0.0] | m | unused while payload_mass is 0 (A-43 M10) | gap |
| `radii_of_gyration` | [0.35, 0.25, 0.25] | 1 | A-43 M12 | estimate |
| `added_mass_coefficients` | [-1.0, -1.5, -1.0, -0.2, -0.8, -1.7] | 1 | A-43 M11 (open: library yaw is -1.2) | estimate |
| `added_mass_form` | "scaled_derivatives" | - | A-43 M11 | selector |
| `added_mass_mass_basis` | "displaced" | - | E-101 Q5 a | owner decision |
| `coriolis_form` | "kirchhoff_full" | - | every added-mass Coriolis term kept | selector |
| `current_form` | "none" | - | this part is gated in still water; the G1 references hold current fixed | selector |

## Grethe hydrostatics (`more_dynamics::SurfaceRestoring`)

| Parameter | Value | Unit | Source | Kind |
|---|---|---|---|---|
| `hull_count` | 1 | - | A-43 M15, monohull | convention |
| `longitudinal_inertia_factor` | 1.0 | 1 | A-43 M17 | estimate |
| `longitudinal_center_of_flotation` | -0.3 | m | A-43 M18 | estimate |
| `reference_point` | [-2.24, 0.0, 0.0] | m | E-101 Q4 a, CO at the midpoint of the waterline length | owner decision |
| `wetted_surface_method` | "regression_table" | - | A-64 Section 3.3 R11 | owner decision |

## Grethe site (`more_dynamics::Site`)

| Parameter | Value | Unit | Source | Kind |
|---|---|---|---|---|
| `site_form` | 'given' | - | no site derivation for this vehicle | selector |
| `latitude` | 1.1073560310932362 | rad | not read with site_form given | document |
| `gravity` | 9.81 | m/s^2 | estimate: the standard value (A-43 M5) | estimate |
| `water_density` | 1025.0 | kg/m^3 | standard sea water (A-43 M5) | estimate |
| `kinematic_viscosity` | 1e-06 | m^2/s | cylinderDrag.m 78-80; Fossen 2011, p. 125, below eq. 6.85 (20 degC) | document |

## Grethe hull loads (`more_dynamics::SurfaceHullLoads`)

| Parameter | Value | Unit | Source | Kind |
|---|---|---|---|---|
| `manoeuvring_damping` | "linear_coupled" | - | A-64 Section 3.4 G0 set | owner decision |
| `sway_yaw_derivatives` | [-251.875, 0.0, 0.0, -229.86112500000002] | N s/m, N m s | [Y_v, Y_r, N_v, N_r]: Y_v = -M22/T_sway, N_r = -M66/T_yaw (time_constants [8, 16] s); Y_r = N_v = 0, unidentified (E-101 Q2 a) | gap |
| `cross_flow` | "hoerner_strips" | - | A-64 Section 3.4 | selector |
| `surge_resistance` | "ittc_residual" | - | E-101 Q1 a | owner decision |
| `time_constant` | 15.0 | s | A-43 M26 (= M20 T1) | estimate |
| `form_factor` | 0.1 | 1 | Fossen 2011 p. 125; A-43 M27 | published |
| `crossover_speed` | 2.0 | m/s | A-43 M28 | published |
| `surge_added_mass_factor` | 2.7 | 1 | addedMassSurge.m 34; A-43 | published |
| `ittc_reynolds_floor` | 100000.0 | 1 | A-43 M29 | published |
| `max_forward_thrust` | 0 | N | zero on purpose: the surge resistance of this part is the ITTC and residual forms alone, a thrust here would count a second resistance (A-74 Q1, owner E-135); a coast-down-fitted resistance option follows | decision |
| `max_speed` | 2.4 | m/s | fitted from her 2026 control-station logs, through water, logged band 2.32-2.47 m/s (A-74 section 1); replace by the full-throttle runs of the sea trials (X-9 G14) | fitted |
| `damping_ratios` | [0.3, 0.05, 0.5] | 1 | [heave, roll, pitch]: roll 0.05 fitted from her logs, band 0.03-0.08 (A-74 section 3, owner E-135); heave 0.3 and pitch 0.5 approved guess 2026-10-12 (owner E-135), no source, until the quay free decays (X-9 G12) | roll fitted; heave, pitch approved guess |

## Plugins without a named part yet

`DerivativeAuv` and `CoefficientHullLoads` (the NPS AUV II vehicle): the named part is not landed (its own job); until it is, their
reference values are the frozen data file `tests/vehicle_models/data/npsauv/nps_auv_ii_parameters.json` (every entry cited to its
`npsauv.m` line, kind `document`, paper not read: needs-access, its `SOURCE.md`). The vehicle declares the 14 acceleration derivatives
and the inertia; `CoefficientHullLoads` declares the 52 hull derivatives and the cross-flow geometry (a derivative is declared once, on
the part that reads it); the 27 actuator-coupled derivatives are declared by no plugin until an actuator design reads them, and stay in
the data file and its `SOURCE.md`. `cross_flow_sections` (the number of strips, 10 for this vehicle) is vehicle data.

## Named parts

| Named part | `.rppws` location | `SOURCE.md` |
|---|---|---|
| REMUS 100 (MSS) | `.rppws/parts/more_dynamics__spheroid_auv/83497a04-1209-469f-bab2-5d7d41d0dee2/` | yes |
| REMUS 100 fin | `.rppws/parts/more_dynamics__fin/` | yes |
| Otter (MSS) | `.rppws/parts/more_dynamics__catamaran/0705552d-7974-4939-b2d0-e9c666424aae/` | yes |
| Grethe | `.rppws/parts/more_dynamics__monohull/83646192-e837-45e4-b1c3-c36ae72940ab/` | `tests/vehicle_models/data/grethe/SOURCE.md` plus the part's own |
| NPS AUV II (MSS) | not yet landed | n/a |
| Scenario site, Scenario boat A / B | `.rppws/parts/more_dynamics__site/`, `...__catamaran/09d55b0c-...`, `...__catamaran/287972df-...` (two boats of the Otter's numbers linking one `Site`; configuration `shared_site`) | no (a scenario example, not a vehicle) |


HullVessel (Luka's `plugins/vehicle_models/hull_vessel.py`, `.rppws/parts/more_dynamics__hull_vessel/`) and JetNozzle
(`.rppws/parts/more_dynamics__jet_nozzle/`) are Luka's own files and are exempt from this file.
