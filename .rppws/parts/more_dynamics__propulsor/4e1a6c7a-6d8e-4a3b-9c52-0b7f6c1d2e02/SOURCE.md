# Otter port thruster

A named part of `more_dynamics::Propulsor`: the port thruster of the Otter, numbers of `otter_legacy.m` @ 6e78e98. The two thrusters are two entries of a vehicle's actuator list (the Otter (MSS) vehicle links both, port first: n(1), n(2) of otter_legacy.m); their summed wrench equals `differential_thruster_casadi` to 1e-10 (`test_propulsor_device_compositions.py`). The shaft-speed limits are derived from the thrust limits as in `otter_legacy.m` 151-152, with the gravity constant of the MSS release (9.81); the frozen MATLAB sweep the gate reads was made with the latitude gravity of an earlier template, so the two agree except on saturated rows. Not wired into the Otter part (the vehicle stays as it is).

| Parameter | Value | Unit | Source | Kind |
|---|---|---|---|---|
| `drive` | speed_command_rad_s | - | the shaft-speed command in rad/s, saturated (otter_legacy.m:236) | selector |
| `steering` | fixed | - | a fixed thruster | selector |
| `mounting` | fixed_direction | - | thrust along x_b at a lever arm (otter_legacy.m:145-146, 228-230) | selector |
| `inflow` | none | - | otter_legacy.m takes nothing from the water for the thrust | selector |
| `law` | quadratic_asymmetric | - | T = k n |n| with k_pos ahead and k_neg astern (otter_legacy.m:238-244) | selector |
| `housing` | none | - | no housing load | selector |
| `min_command` | -101.73666550381692 | rad/s | otter_legacy.m:152 with g = 9.81 (mssConstants.m:14): -sqrt(0.5 * 13.6 g / k_neg) | derived |
| `max_command` | 103.93086427398005 | rad/s | otter_legacy.m:151 with g = 9.81: sqrt(0.5 * 24.4 g / k_pos) | derived |
| `position` | [0.0, -0.395, 0.0] | m | otter_legacy.m:118, 145-146 (y_pont = 0.395) | document |
| `direction` | [1.0, 0.0, 0.0] | 1 | otter_legacy.m:248 | document |
| `positive_thrust_coefficient` | 0.01108 | N s^2/rad^2 | otter_legacy.m:149 | document |
| `negative_thrust_coefficient` | 0.006445 | N s^2/rad^2 | otter_legacy.m:150 | document |
| `shaft_torque_scale` | 0.0 | 1 | otter_legacy.m has no thruster torque | convention |
