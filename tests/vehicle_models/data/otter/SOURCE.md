# Otter composition inputs

MSS `CRAFT/USV/models/otter.m` at `cc07579` is the source of the Otter values in `otter_mss_cc07579.json` (history; see the last section for `99fdf4d`): gravity and density, lines 90–91; hull dimensions and coefficients, lines 92–107; payload mass and position, `tests/data/hydrodynamics/generate_catamaran_mss.m` lines 61–62; added-mass coefficients, lines 152–157; restoring inputs, lines 177 and 192; damping constants, lines 99–101 and 202–240. `max_forward_thrust` is 24.4 × 9.81 N; `max_speed` is 6 × 0.5144 m/s, both from the MSS file. No new reference outputs are frozen here. The gates read `tests/data/rigid_body/{inputs.csv,matlab_reference_mss_current.csv}` and `tests/data/hydrodynamics/catamaran_matlab_reference_mss_cc07579.csv`.

The composition uses `site_given`, `no_current` or `horizontal_current_full_rotation_rate`, `twin_pontoons`, `hull_with_point_payload`, `scaled_derivatives`, `kirchhoff_full`, `twin_hull_metacentric_equilibrium_draft`, `time_constant_damping_surface`, and `cross_flow_strip+rectangular_section_hoerner`. Its damping wiring maps `restoring_matrix` to `G_CF`; the cross-flow wiring maps `beam` to `pontoon_beam`. The Otter and REMUS compositions use the same `cross_flow_strip` part module with different section-drag parts.

## Re-pin to MSS `99fdf4d` (2026-10-23, MATLAB R2026a Update 4; `99fdf4dbdfe0b3f124e26a689c6aabe0beb5654f`, 41 commits after `cc07579`)

`otter_mss_99fdf4d.json` replaces `otter_mss_cc07579.json` as the composition input of the gates (the old file stays): the lines are those of `otter.m` at `99fdf4d` (gravity and density through `mssConstants.m`, lines 95-97; hull 98-104; `T_surge` 5 s at 105, `T_sway` 0.5 s at 106, `T_yaw` 107; payload added to the surge added mass at 162; added-mass coefficients 161-166; surge law 211-212, 245-246; damping ratios 214-216, yaw 251). New values over the old file: `time_constants` [0.5, 1.0], `surge_time_constant` 5.0, `surge_blend_decay` 3.0 (`exp(-3 |u_r|)`). The composition now selects `added_mass_form = scaled_derivatives_loaded_surge` and `manoeuvring_damping = time_constants_surge_blend`. The gates read `tests/rigid_body/data/matlab_reference_mss_99fdf4d.csv` and `tests/hydrodynamic_loads/data/catamaran_matlab_reference_mss_99fdf4d.csv`.

The part is the MSS form bit for bit: the payload trim `eta_0` (`otter.m` 264-266) with the heave double count (F-MSS-31) and the shifted kinematics (F-MSS-24) are still in MSS at `99fdf4d` and still in this reference; the corrected form is J1b.

## Commanded runs at MSS `ed2ce7b` (2026-10-26, MATLAB R2026a Update 4; Y3)

`generate_otter_commanded_mss.m` runs the unmodified `otter.m` (revision `ed2ce7b`) with commanded propeller shaft speeds
`n = [n_left n_right]` (rad/s): `rk4.m` at the step of `SIMotter.m` (0.05 s), the command decided every 5 steps and held,
payload 25 kg at [0.05, 0, -0.35] m, no current. Files, all with `_mss_ed2ce7b` before the extension:

| File | What it holds |
|---|---|
| `otter_trajectory_<run>.csv` | `t`, `ui_01..02` (shaft speeds, rad/s, port then starboard), `Vc`, `beta_c`, `x_01..12` (MSS order `[nu; eta]`), a row every 0.25 s; runs `straight_from_rest_thrust_step` (16 s), `differential_turn` (60 s), `astern_from_speed` (16 s), `spin_in_place` (40 s), `commands_beyond_the_limits` (30 s, shaft speeds beyond `n_max`, `n_min`), `zigzag_differential_steps` (60 s) |
| `otter_trajectory_tolerances.csv` | per run and state, the measured difference of two integrations of the same run (the derivative with +-1e-9, and the RK4 sum reordered), the G3 tolerance; the same rule as the REMUS 100 runs |
| `otter_commanded_meta.json` | revision, step, hold, payload, shaft-speed limits |

The runs are the unmodified `otter.m`, so they carry its payload trim as MSS has it, the heave double count (F-MSS-31, a
net vertical force of `m_P g` at the shifted equilibrium) and the kinematics at the shifted attitude (F-MSS-24): the
heave starts from the unshifted state and settles at the shift (z = 0.0377 m for 25 kg). The part `Otter (MSS)` selects
the payload trim form that reproduces this. The two straight runs are 16 s: the open-loop yaw is unstable (Munk moment),
so at 40 s a derivative difference of 1e-9 grows to 0.5 rad and the measured tolerance would be 0.4 rad / 2.8 m; at 16 s it
is below 1e-5 in every state. Re-run: `MSS_DIR=<MSS> OUT_DIR=<scratch> "$MATLAB_BIN" -batch "run('<generator>')"`, then
`cmp` each file with its frozen copy; two runs on this machine gave byte-identical files. Readers:
`tests/vehicle_models/test_rppws_runs.py` (G3 through rpp's builder from the named part, the schedule of
`vehicle_simulation.json` sample by sample, CVODES against a finer RK4).
