# REMUS 100 propulsor

A named part of `more_dynamics::Propulsor`: the numbers of the `REMUS 100 propeller` part (same lines of `remus100.m` @ ed2ce7b), in the Propulsor's own names. The `REMUS 100 (MSS)` vehicle links it as its propulsor (the shaft speed command `drive_command`, rpm); the `Propeller` plugin with the `REMUS 100 propeller` part gives the same wrench to 1e-10 on random states (`tests/force_producers/propulsor/test_propulsor_device_compositions.py`) and the same vehicle derivative (`tests/vehicle_models/test_named_vehicles.py`).

| Parameter | Value | Unit | Source | Kind |
|---|---|---|---|---|
| `drive` | speed_command_rpm | - | shaft-speed command, rpm | selector |
| `steering` | fixed | - | a fixed shaft | selector |
| `mounting` | fixed_euler | - | shaft axis from Euler angles | selector |
| `inflow` | speed_magnitude | - | remus100.m:151, 250 | selector |
| `law` | open_water_linearized_given | - | remus100.m:168-178, coefficients as printed | selector |
| `housing` | none | - | no housing load | selector |
| `min_command` | -1525 | rpm | remus100.m:111, 116 (-n_max) | document |
| `max_command` | 1525 | rpm | remus100.m:111 n_max | document |
| `thrust_deduction` | 0.1 | 1 | remus100.m:150 t_prop | document |
| `wake_fraction` | 0.05600000000000005 | 1 | remus100.m:151 Va = 0.944 U_r, w = 1 - 0.944 | derived |
| `position` | [0, 0, 0] | m | remus100.m:250-253 (thrust on x_b through the CO) | document |
| `orientation` | [0, 0, 0] | rad | remus100.m:250, 253 (shaft along x_b) | document |
| `diameter` | 0.14 | m | remus100.m:149 D_prop | document |
| `pitch_diameter_ratio` | 1 | 1 | remus100.m:157 | document |
| `blade_area_ratio` | 0.718 | 1 | remus100.m:156-157 | document |
| `blade_count` | 3 | 1 | remus100.m:156-157 | document |
| `max_advance_number` | 0.6632 | 1 | remus100.m:154 Ja_max | document |
| `thrust_torque_coefficients` | [0.4566, 0.07, 0.1798, 0.0312] | 1 | remus100.m:158, 159, 161, 162 | document |
| `shaft_torque_scale` | 0.1 | 1 | remus100.m:253 K_prop / 10 | document (a comparison scale, not the physical 1) |
