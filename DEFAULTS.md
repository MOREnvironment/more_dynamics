# Defaults

Rule 20 (AGENTS.md, E-107): plugin and model code never names a vehicle or an author-form.
A plugin's `PARAMETERS` still carry working default *values* (rpp requires one), but the
*provenance* of each value — where it comes from and which named part it was first read
for — lives here, one table per plugin, instead of in the plugin's own docstring or comments.

A file citation on the code line itself (`remus100.m:207`, `otter.m 94-98`) is not a vehicle
name under rule 20 (rule 13's carve-out): those stay in the code. What moves here is the
*prose* that used to say "defaults are the REMUS 100" or "the Otter" beside them.

Columns: **Default** is the plugin's own `PARAMETERS` value (arrays shown whole). **Source**
is the cited file:line or document this number was read from. **Named part** is the
`.rppws` part (E-89, E-106) the value was first read for — the vehicle whose numbers these
are, kept here rather than in the library. **Kind**: `document` (read off a cited page or
file), `estimate` (a starting point with no document behind it, to be identified),
`derived` (computed from another primitive by a stated equation), `selector` (chooses a
form, not a physical quantity), `convention` (a modelling choice, not measured), `owner
decision` or `gap` (named in the owner's own words where the part's `SOURCE.md` uses them,
e.g. Grethe's).

## `more_dynamics::SpheroidAuv` (`plugins/vehicle_models/spheroid_auv.py`)

| Parameter | Default | Unit | Source | Named part | Kind |
|---|---|---|---|---|---|
| `site_form` | `latitude` | - | `remus100.m` 96-97 | REMUS 100 (MSS) | selector |
| `latitude` | 1.1073560310932362 | rad | `remus100.m` 96-97 (mu = deg2rad(63.446827)) | REMUS 100 (MSS) | document |
| `gravity` | 9.81 | m/s^2 | `mssConstants.m`:14 @ `49e03e3` (MSS 2.0.5) | REMUS 100 (MSS) | document |
| `water_density` | 1026 | kg/m^3 | `remus100.m`:98 | REMUS 100 (MSS) | document |
| `kinematic_viscosity` | 1e-06 | m^2/s | `cylinderDrag.m` 78-80; Fossen 2011, p. 125, below eq. 6.85 (20 degC) | - (water property) | document |
| `length` | 1.6 | m | `remus100.m`:131 (L_auv) | REMUS 100 (MSS) | document |
| `diameter` | 0.19 | m | `remus100.m`:132 (D_auv) | REMUS 100 (MSS) | document |
| `body_density_method` | `computed` | - | `spheroid.m`:36 solved for rho_b | REMUS 100 (MSS) | selector |
| `body_mass` | 31.9 | kg | `remus100.m`:3 | REMUS 100 (MSS) | document |
| `body_density` | 1054.7872613500267 | kg/m^3 | given value, read with `"given"`: 31.9 kg / (4/3 pi a b^2) | REMUS 100 (MSS) | derived |
| `body_center_of_gravity` | [0, 0, 0.02] | m | `remus100.m`:137 (r_bG) | REMUS 100 (MSS) | document |
| `roll_added_inertia_ratio` | 0.3 | 1 | `remus100.m`:136 (r44) | REMUS 100 (MSS) | document |
| `added_mass_form` | `lamb_spheroid` | - | `imlay61.m` 31-59 | REMUS 100 (MSS) | selector |
| `coriolis_form` | `kirchhoff_full` | - | `m2c.m` 33-48 | - (physics default) | selector |
| `current_form` | `full_attitude` | - | Fossen 2011, eqs. 8.138-8.141, 8.157, pp. 221-225 | - (physics default) | selector |

## `more_dynamics::Catamaran` (`plugins/vehicle_models/catamaran.py`)

| Parameter | Default | Unit | Source | Named part | Kind |
|---|---|---|---|---|---|
| `site_form` | `given` | - | `otter.m` 90-91 | Otter (MSS) | selector |
| `latitude` | 1.1073560310932362 | rad | `remus100.m` 96-97 (the shared MSS site; unread with `site_form="given"`) | REMUS 100 (MSS) | document |
| `gravity` | 9.81 | m/s^2 | `otter.m` 90-91 | Otter (MSS) | document |
| `water_density` | 1025.0 | kg/m^3 | `otter.m` 90-91 | Otter (MSS) | document |
| `kinematic_viscosity` | 1e-06 | m^2/s | `cylinderDrag.m` 78-80; Fossen 2011, p. 125, below eq. 6.85 | - (water property) | document |
| `length` | 2.0 | m | `otter.m` 92-93, 104-107 | Otter (MSS) | document |
| `beam` | 1.08 | m | `otter.m` 92-93, 104-107 | Otter (MSS) | document |
| `pontoon_beam` | 0.25 | m | `otter.m` 104-107 | Otter (MSS) | document |
| `pontoon_lateral_offset` | 0.395 | m | `otter.m` 104-107 | Otter (MSS) | document |
| `pontoon_block_coefficient` | 0.4 | 1 | `otter.m` 104-107 | Otter (MSS) | document |
| `pontoon_waterplane_coefficient` | 0.75 | 1 | `otter.m` 104-107 | Otter (MSS) | document |
| `hull_mass` | 55.0 | kg | `otter.m` 94-98 | Otter (MSS) | document |
| `payload_mass` | 25.0 | kg | `generate_catamaran_mss.m` 61-62 | Otter (MSS) | document |
| `hull_center_of_gravity` | [0.2, 0.0, -0.2] | m | `otter.m` 94-98 | Otter (MSS) | document |
| `payload_position` | [0.05, 0.0, -0.35] | m | `generate_catamaran_mss.m` 61-62 | Otter (MSS) | document |
| `radii_of_gyration` | [0.4, 0.25, 0.25] | 1 | `otter.m` 94-98 | Otter (MSS) | document |
| `added_mass_coefficients` | [-1.0, -1.5, -1.0, -0.2, -0.8, -1.7] | 1 | `otter.m` 152-157 | Otter (MSS) | document |
| `added_mass_form` | `scaled_derivatives` | - | `otter.m` 152-159 | Otter (MSS) | selector |
| `coriolis_form` | `kirchhoff_full` | - | `m2c.m` 33-48 | - (physics default) | selector |
| `current_form` | `full_attitude` | - | Fossen 2011, eqs. 8.138-8.141, 8.157, pp. 221-225 | - (physics default) | selector |

## `more_dynamics::Monohull` (`plugins/vehicle_models/monohull.py`)

Every value below is Grethe's (NTNU Mariner 5 USV, Maritime Robotics build 40402, Pioner 17
ft hull): the frozen, hand-edited source is `tests/vehicle_models/data/grethe/grethe_parameters.json`
(A-43 Part 2, A-64, E-101), which carries its own per-value `kind` and `place` fields —
reproduced here for the plugin's own declared parameters only (the full set, including
options not yet wired as plugin selectors, is in that file and the `Grethe` named part's
own `SOURCE.md`).

| Parameter | Default | Unit | Source | Named part | Kind |
|---|---|---|---|---|---|
| `site_form` | `given` | - | `otter.m` 90-91 (plugin default, not Grethe-specific) | - (physics default) | selector |
| `latitude` | 1.1073560310932362 | rad | `remus100.m` 96-97 (the shared MSS site; unread with `site_form="given"`) | REMUS 100 (MSS) | document |
| `gravity` | 9.81 | m/s^2 | `grethe_mariner5.yaml hydrostatics.gravity` ("standard") | Grethe | estimate |
| `water_density` | 1025.0 | kg/m^3 | standard sea water | Grethe | document |
| `kinematic_viscosity` | 1e-06 | m^2/s | Fossen 2011, p. 125, below eq. 6.85 (20 degC) | - (water property) | document |
| `length` | 5.2 | m | F003 p2, M001 p7 (overall length) | Grethe | document |
| `beam` | 2.15 | m | F003 p2, M001 p7, D002 p1 | Grethe | document |
| `waterline_length` | 0.0 (not given; 4.63 m is the drawing estimate) | m | D002 p1 drawing measurement at T = 0.30 m; not yet wired | Grethe | estimate |
| `waterline_beam` | 0.0 (not given; 1.97 m is the drawing estimate) | m | D002 p1 drawing measurement at T = 0.30 m; not yet wired | Grethe | estimate |
| `hull_block_coefficient` | 0.233 | 1 | 806 / (1025 * 5.2 * 2.15 * 0.3), a mass knob not hull fullness | Grethe | derived |
| `hull_waterplane_coefficient` | 0.8 | 1 | A-43 M14 | Grethe | estimate |
| `hull_mass` | 806.0 | kg | F003 p2 (FAT, +/- 20 kg) | Grethe | document |
| `payload_mass` | 0.0 | kg | the 2026 field payload (USBL mount, battery packs) is not weighed | Grethe | gap |
| `hull_center_of_gravity` | [0.0, 0.0, 0.025] | m | A-43 M9 | Grethe | estimate |
| `payload_position` | [0.0, 0.0, 0.0] | m | unused while `payload_mass` is 0 | Grethe | gap |
| `radii_of_gyration` | [0.35, 0.25, 0.25] | 1 | A-43 M12 | Grethe | estimate |
| `added_mass_coefficients` | [-1.0, -1.5, -1.0, -0.2, -0.8, -1.7] | 1 | A-43 M11 (the workspace `.rppws` value; open: the library value is -1.2 in yaw) | Grethe | estimate |
| `added_mass_form` | `scaled_derivatives` | - | `otter.m` 152-159 | - (physics default) | selector |
| `added_mass_mass_basis` | `hull` | - | keeps the committed gates green; a payload-carrying part may select `displaced` | - (physics default) | selector |
| `coriolis_form` | `kirchhoff_full` | - | `m2c.m` 33-48 | - (physics default) | selector |
| `current_form` | `full_attitude` | - | Fossen 2011, eqs. 8.138-8.141, 8.157, pp. 221-225 | - (physics default) | selector |

## `more_dynamics::DerivativeAuv` (`plugins/vehicle_models/derivative_auv.py`)

Every value is Healey, A. J. & Lienard, D. B. (1993), via MSS `CRAFT/AUV/models/npsauv.m`
104-172 @ `cc07579`; the paper itself is not read (needs-access, `SOURCE.md`), so every
number here is read from MSS's instrumented workspace, cited to its line.

| Parameter | Default | Unit | Source | Named part | Kind |
|---|---|---|---|---|---|
| `gravity` | 9.81 | m/s^2 | `npsauv.m`:105 (given directly, not from a site) | NPS AUV II (MSS) | document |
| `water_density` | 1025 | kg/m^3 | `npsauv.m`:105 | NPS AUV II (MSS) | document |
| `length` | 5.3 | m | `npsauv.m`:105 | NPS AUV II (MSS) | document |
| `body_mass` | 5443.425076452599 | kg | W / g (`npsauv.m`:109), W at `npsauv.m`:108 | NPS AUV II (MSS) | derived |
| `inertia_diagonal` | [2038, 13587, 13587] | kg*m^2 | `npsauv.m`:110 | NPS AUV II (MSS) | document |
| `inertia_products` | [-13.58, -13.58, -13.58] | kg*m^2 | `npsauv.m`:111 | NPS AUV II (MSS) | document |
| `body_center_of_gravity` | [0, 0, 0.061] | m | `npsauv.m`:106 (r_bg) | NPS AUV II (MSS) | document |
| `weight` | 53400 | N | `npsauv.m`:108 (W) | NPS AUV II (MSS) | document |
| `buoyancy` | 53400 | N | `npsauv.m`:108 (B, neutral) | NPS AUV II (MSS) | document |
| `center_of_buoyancy` | [0, 0, 0] | m | `npsauv.m`:107 (r_bb) | NPS AUV II (MSS) | document |
| `derivative_scaling` | `nondimensional` | - | `npsauv.m` 116-124 (Fossen prime-scaling) | NPS AUV II (MSS) | selector |
| `current_form` | `yaw_rate_terms` | - | `npsauv.m` 83-90 (the gate's own MSS form) | NPS AUV II (MSS) | selector |
| `site_form` | `given` | - | `npsauv.m`:105 (no latitude in this reference) | NPS AUV II (MSS) | selector |
| `latitude` | 0.0 | rad | unread with `site_form="given"` (declared for the shared form) | - | gap |
| `kinematic_viscosity` | 1e-06 | m^2/s | unread by this vehicle's own graph (site's own input) | - | gap |
| `nondim_*` (93 published nondimensional derivatives) | per-derivative, see `npsauv.m` 126-156 | nondimensional (prime system) | `npsauv.m` 126-156 | NPS AUV II (MSS) | document |

## `more_dynamics::Fin` (`plugins/force_producers/fin.py`)

| Parameter | Default | Unit | Source | Named part | Kind |
|---|---|---|---|---|---|
| `fin_position` | [-0.638, 0.0, 0.0] | m | Prestero 2001, Table A.5, p. 103 (x_finpost) | REMUS 100 (MSS) | document |
| `chord_axis` | [1.0, 0.0, 0.0] | 1 | convention (unit chord direction, BODY FRD) | - | convention |
| `lift_axis` | [0.0, 1.0, 0.0] | 1 | convention (unit positive-lift direction, BODY FRD) | - | convention |
| `fin_area` | 0.00665 | m^2 | Prestero 2001, Table A.5, p. 103 (S_fin) | REMUS 100 (MSS) | document |
| `servo_dynamics` | `first_order_lag` | - | Murray-Smith 2016, Fig. 1, p. 246 | - | selector |
| `servo_rate_limit` | `True` | - | Murray-Smith 2016, Fig. 1, p. 246 | - | selector |
| `servo_angle_limit` | `on_command` | - | Murray-Smith 2016, Fig. 1, p. 246 | - | selector |
| `time_constant` | 0.1 | s | Sarhadi 2026, Fig. 4, p. 4 | - (a demanding choice of that study, not a measurement) | estimate |
| `max_deflection` | 0.3490658503988659 (20 deg) | rad | Sarhadi 2026, Fig. 4, p. 4 | - | estimate |
| `max_rate` | 0.5235987755982988 (30 deg/s) | rad/s | Sarhadi 2026, Fig. 4, p. 4 | - | estimate |
| `inflow` | `rigid_point` | - | Prestero 2001, eq. 4.40, p. 32 | - | selector |
| `flow_angle` | `small_angle` | - | Prestero 2001, eqs. 4.42-4.43, pp. 32-33 | - | selector |
| `interference` | `none` | - | Prestero 2001, eq. 4.41, p. 32 (both factors one) | - | selector |
| `section` | `quadratic_drag` | - | `remus100.m` 238-245 | REMUS 100 (MSS) | selector |
| `lift_slope` | 3.12 | 1/rad | Prestero 2001, Table A.5, p. 103 (c_L_alpha) | REMUS 100 (MSS) | document |
| `zero_lift_drag` | 0.0 | 1 | Prestero 2001, eq. 4.37, p. 31 (carries no fin drag; read by `linear_section`) | REMUS 100 (MSS) | document |

## `more_dynamics::FinPairsDeflectionOnly` (`plugins/force_producers/fin_pairs_deflection_only.py`)

| Parameter | Default | Unit | Source | Named part | Kind |
|---|---|---|---|---|---|
| `rudder_area` | 0.0133 | m^2 | `remus100.m`:183 (A_r = 2 * S_fin) | REMUS 100 (MSS) | document |
| `stern_plane_area` | 0.0133 | m^2 | `remus100.m`:188 (A_s = 2 * S_fin) | REMUS 100 (MSS) | document |
| `rudder_lift_coefficient` | 0.5 | 1/rad | `remus100.m`:182 (CL_delta_r) | REMUS 100 (MSS) | document |
| `stern_plane_lift_coefficient` | 0.7 | 1/rad | `remus100.m`:187 (CL_delta_s) | REMUS 100 (MSS) | document |
| `fin_positions_method` | `computed` | - | x_r = x_s = -a (`remus100.m`:184, 189) | REMUS 100 (MSS) | selector |
| `rudder_position` | -0.8 | m | given value, read with `"given"`: `remus100.m`:184 (x_r = -a, a = L_auv / 2) | REMUS 100 (MSS) | derived |
| `stern_plane_position` | -0.8 | m | given value, read with `"given"`: `remus100.m`:189 (x_s = -a, a = L_auv / 2) | REMUS 100 (MSS) | derived |
| `max_deflection` | 0.3490658503988659 (20 deg) | rad | `remus100.m`:109 (delta_max = deg2rad(20)) | REMUS 100 (MSS) | document |

## `more_dynamics::Propeller` (`plugins/force_producers/propeller.py`)

| Parameter | Default | Unit | Source | Named part | Kind |
|---|---|---|---|---|---|
| `propeller_diameter` | 0.14 | m | `remus100.m`:148 (D_prop) | REMUS 100 (MSS) | document |
| `max_shaft_speed` | 1525 | rpm | `remus100.m`:110 (n_max) | REMUS 100 (MSS) | document |
| `thrust_deduction` | 0.1 | 1 | `remus100.m`:149 (t_prop) | REMUS 100 (MSS) | document |
| `wake_fraction` | 0.05600000000000005 | 1 | `remus100.m`:150 (Va = 0.944 * U_r, w = 1 - 0.944) | REMUS 100 (MSS) | derived |
| `pitch_diameter_ratio` | 1 | 1 | `remus100.m`:156 (`wageningen(0,1,0.718,3)`, second argument) | REMUS 100 (MSS) | document |
| `blade_area_ratio` | 0.718 | 1 | `remus100.m` 155-156 | REMUS 100 (MSS) | document |
| `blade_count` | 3 | 1 | `remus100.m` 155-156 | REMUS 100 (MSS) | document |
| `max_advance_number` | 0.6632 | 1 | `remus100.m`:153 (Ja_max) | REMUS 100 (MSS) | document |
| `roll_moment_scale` | 0.1 | 1 | `remus100.m`:252 (tau(4) = K_prop / 10) | REMUS 100 (MSS) | document |
| `position` | [0, 0, 0] | m | `remus100.m` 249-252 (thrust on x_b through the CO, no moment arm) | REMUS 100 (MSS) | document |
| `orientation` | [0, 0, 0] | rad | `remus100.m` 249, 252 (shaft along x_b) | REMUS 100 (MSS) | document |
| `thrust_torque_coefficients` | [0.4566, 0.07, 0.1798, 0.0312] | 1 | `remus100.m` 157, 158, 160, 161 ([KT_0 KQ_0 KT_max KQ_max]) | REMUS 100 (MSS) | document |

## `more_dynamics::AuvHullLoads` (`plugins/hydrodynamics/auv_hull_loads.py`)

| Parameter | Default | Unit | Source | Named part | Kind |
|---|---|---|---|---|---|
| `time_constants` | [20, 20, 1] | s | `remus100.m` 192, 193, 196 ([T1 T2 T6]) | REMUS 100 (MSS) | document |
| `damping_ratios` | [0.3, 0.8] | 1 | `remus100.m` 194, 195 ([zeta4 zeta5]) | REMUS 100 (MSS) | document |
| `oswald_efficiency` | 0.3 | 1 | `coeffLiftDrag.m`:54 (e = 0.3) | REMUS 100 (MSS) | document |
| `planform_area_method` | `computed` | - | S = f L D (`remus100.m`:133) | REMUS 100 (MSS) | selector |
| `planform_fraction` | 0.7 | 1 | `remus100.m`:133 (S = 70% of the rectangle L_auv * D_auv) | REMUS 100 (MSS) | document |
| `planform_area` | 0.2128 | m^2 | given value, read with `"given"`: `remus100.m`:133 | REMUS 100 (MSS) | derived |
| `parasitic_drag_coefficient_method` | `computed` | - | CD_0 = Cd pi (D/2)^2 / S (`remus100.m` 142-144) | REMUS 100 (MSS) | selector |
| `hull_frontal_drag_coefficient` | 0.42 | 1 | `remus100.m`:142 (Cd), Allen et al. (2000) as MSS cites it, not read | REMUS 100 (MSS) | document |
| `parasitic_drag_coefficient` | 0.05595961914206819 | 1 | given value, read with `"given"`: `remus100.m` 143-144 | REMUS 100 (MSS) | derived |
| `smooth_speed` | `True` | - | numerical guard, finite Jacobian at rest | - | selector |
| `smooth_speed_epsilon` | 1e-09 | m/s | numerical regularisation of the guard, not a vehicle quantity | - | convention |

## `more_dynamics::CoefficientHullLoads` (`plugins/hydrodynamics/coefficient_hull_loads.py`)

The 52 hull-only nondimensional derivatives have no generic physics default (an empirical
coefficient table): each defaults to 0.0 and is given entirely by a composition. Today's
only supplying composition is the frozen test fixture
`tests/vehicle_models/data/npsauv/nps_auv_ii_parameters.json` (U6e, the `NPS AUV II`
named part itself is not yet landed in `.rppws`, `test_named_part_nps_auv_ii_exists_in_the_root_rppws`).

| Parameter | Default | Unit | Source | Named part | Kind |
|---|---|---|---|---|---|
| every hull-only derivative (52, `npsauv.m` 126-156) | 0.0 | nondimensional | none (empirical table; no physics default) | NPS AUV II (test fixture, not yet a `.rppws` part) | gap |
| `cross_flow_drag_coefficients` | [0.5, 0.6] | 1 | `npsauv.m`:220 ([Cdy, Cdz]) | NPS AUV II (test fixture) | document |
| `cross_flow_section_dimensions` | [0.53, 0.53] | m | `npsauv.m` 220-221 ([Hx, Bx]) | NPS AUV II (test fixture) | document |
| `cross_flow_sections` | 10 | 1 | `npsauv.m`:218 (number of strips) | NPS AUV II (test fixture) | document |

## `more_dynamics::SurfaceHullLoads` (`plugins/hydrodynamics/surface_hull_loads.py`)

| Parameter | Default | Unit | Source | Named part | Kind |
|---|---|---|---|---|---|
| `max_forward_thrust` | 239.364 | N | `otter.m` (24.4 kgf = 24.4 * 9.81 N) | Otter (MSS) | document |
| `max_speed` | 3.0864 | m/s | `otter.m` (6 knots = 6 * 0.5144 m/s) | Otter (MSS) | document |
| `time_constants` | [1.0, 1.0] | s | `otter.m` 99-101, 202-240 | Otter (MSS) | document |
| `damping_ratios` | [0.3, 0.2, 0.4] | 1 | `otter.m` 99-101, 202-240 | Otter (MSS) | document |
| `yaw_damping_nonlinearity` | 10.0 | s/rad | `otter.m` 99-101, 202-240 | Otter (MSS) | document |
| `manoeuvring_damping` | `time_constants` | - | plugin default | - | selector |
| `sway_yaw_derivatives` | [0.0, 0.0, 0.0, 0.0] | N s/m, N m s | read with `linear_coupled`/`linear_coupled_modulus` | - | gap |
| `modulus_coefficients` | [0.0, 0.0] | - | read with `linear_coupled_modulus`; zero (module docstring) | - | gap |
| `cross_flow` | `hoerner_strips` | - | plugin default | - | selector |
| `surge_resistance` | `none` | - | plugin default (off unless a composition turns it on) | - | selector |
| `time_constant` | 100.0 | s | `osv.m`:128 (vessel.T1 = 100 s) | - (ship-scale MSS calibration) | document |
| `form_factor` | 0.1 | 1 | `XuuITTC.m`:33 (k = 0.1); Fossen 2011, p. 125 | - | document |
| `crossover_speed` | 2.0 | m/s | `forceSurgeDamping.m`:57 @ `ac77394` (u_cross = 2 m/s) | - | document |
| `surge_added_mass_factor` | 2.7 | 1 | `addedMassSurge.m`:33 (A11 = 2.7 rho nabla^(5/3) / L^2) | - | document |
| `ittc_reynolds_floor` | 100000.0 | 1 | `XuuITTC.m`:35 (Re_min = 1e5) | - | document |

## `more_dynamics::SubmergedRestoring` (`plugins/hydrostatics/submerged_restoring.py`)

| Parameter | Default | Unit | Source | Named part | Kind |
|---|---|---|---|---|---|
| `buoyancy_method` | `neutral` | - | `remus100.m`:214 (B = W) | REMUS 100 (MSS) | selector |
| `center_of_buoyancy` | [0, 0, 0] | m | `remus100.m`:138 (r_bB) | REMUS 100 (MSS) | document |
| `displaced_volume` | 0.030243065278557742 | m^3 | `remus100.m` 131-132 (4/3 pi (L/2)(D/2)^2) | REMUS 100 (MSS) | derived |
| `buoyancy` | 313.31493718007266 | N | m g with m = 31.9 kg (`remus100.m`:3), g at `remus100.m` 96-97 | REMUS 100 (MSS) | derived |

## `more_dynamics::SurfaceRestoring` (`plugins/hydrostatics/surface_restoring.py`)

| Parameter | Default | Unit | Source | Named part | Kind |
|---|---|---|---|---|---|
| `hull_count` | 2 | - | `otter.m` 172-193 (two pontoons) | Otter (MSS) | document |
| `longitudinal_inertia_factor` | 0.8 | 1 | `otter.m`:179 (I_L = 0.8 ...) | Otter (MSS) | document |
| `longitudinal_center_of_flotation` | -0.2 | m | `otter.m` 177, 192 | Otter (MSS) | document |
| `reference_point` | [0.0, 0.0, 0.0] | m | convention (the body-frame origin, CO) | - | convention |
| `wetted_surface_method` | `computed` | - | S = n 1.025 L (C_b B + 1.7 T), `XuuITTC.m`:38 (Mumford) | - | selector |
| `wetted_surface` | 1.77 | m^2 | given value, read with `"given"`: the Mumford value at the Otter set (`otter.m` 92-107, mass 80 kg, `XuuITTC.m`:38) | Otter (MSS) | derived |
| `block_coefficient_method` | `given` | - | E-104 (keeps every existing gate; `hull_block_coefficient` a vehicle coupling, draft derived) | - | selector |
| `draft` | 0.3 | m | read with `"computed"` only: document (F003 p2; M001 p7, approx. 0.3 m) | Grethe | document |

## Named parts

| Named part | `.rppws` location | `SOURCE.md` |
|---|---|---|
| REMUS 100 (MSS) | `.rppws/parts/more_dynamics__spheroid_auv/83497a04-1209-469f-bab2-5d7d41d0dee2/` | yes (this job) |
| Otter (MSS) | `.rppws/parts/more_dynamics__catamaran/0705552d-7974-4939-b2d0-e9c666424aae/` | yes (this job) |
| Grethe | `.rppws/parts/more_dynamics__monohull/` (U8b leftover, this job) | `tests/vehicle_models/data/grethe/SOURCE.md` plus the part's own |
| NPS AUV II (MSS) | not yet landed (U6e) | n/a |

HullVessel (Luka's `plugins/vehicle_models/hull_vessel.py`, `.rppws/parts/more_dynamics__hull_vessel/`)
and JetNozzle (`.rppws/parts/more_dynamics__jet_nozzle/`) are Luka's own files (E-28, E-91) and
are exempt from rule 20 and from this file.
