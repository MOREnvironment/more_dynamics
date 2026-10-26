# Grethe outboard

A named part of `more_dynamics::Propulsor`: the one steerable electric outboard of Grethe (NTNU Mariner 5 USV, Maritime
Robotics build 40402, a Torqeedo outboard), every value and form as in the model file `grethe.m` (MSS fork branch
`grethe-vehicle`, revision `679a6e8`, `CRAFT/USV/models/grethe.m`; the file may still be revised, and `SIMgrethe.m` 151-154,
257-259 for the actuator equations the model file leaves to its caller). The values are the file's, not new ones: a value
changes here when the file changes, or when a fitted or measured value replaces its estimate. The Grethe vehicle part
links it as its only actuator, and so does the stress-test part.

Inputs, in this order: `drive_command` (the thrust command, N) and `steering_command` (the steering angle command, rad,
positive = thrust force to starboard, bow to port for a unit aft of the CO). States, in this order: `drive_state` (the
actual thrust, N) and `servo_state` (the actual steering angle, rad).

Forms: the thrust is the command through a first-order lag (`thrust_command_lag`) and then limited to the smaller of the
maximum thrust and the power delivered to the water over the horizontal speed (`thrust_power_limited`, with the speed
floored, `horizontal_speed_floored`); the steering is the command through a lag with rate and angle limits
(`servo_lag_rate_angle`); the thrust acts at the motor (`azimuth`), plus the normal force of the leg
(`leg_normal_force`). The model file takes the *actual* thrust and angle as its inputs and leaves the lags to its caller
(`SIMgrethe.m`); here they are states of the part, so the vehicle's `graph()` is the plant of the model file with the
caller's actuator equations (in continuous form, integrated with the vehicle's own RK4; `SIMgrethe.m` takes an Euler step).

| Parameter | Value | Unit | Source | Kind |
|---|---|---|---|---|
| `drive`, `steering`, `mounting`, `inflow`, `law`, `housing` | see `params/parameters.py` | - | the forms above | selector |
| `steering_convention` | port_positive | - | `grethe.m` 38-48, 276 names the same formula `starboard_positive` (the name of the direction of the force); the library names the boat's turn | selector |
| `shaft_torque_scale` | 0.0 | 1 | the model file has no shaft torque on the body | convention |
| `thrust_time_constant` | 0.3 | s | `grethe.m` 268 `T_thrust` | estimate |
| `max_deflection` | 0.6108652381980153 (35 deg) | rad | `grethe.m` 265 `delta_max` | estimate (the logs imply 53-80 deg in turns from rest, Y16 O2, conditional on the inertia estimate; the steering sweep of X-9 decides) |
| `max_rate` | 1.0471975511965976 (60 deg/s) | rad/s | `grethe.m` 266 `ddelta_max` | estimate |
| `time_constant` | 0.5 | s | `grethe.m` 267 `T_delta` | estimate |
| `position` | [-2.24, 0, 0.55] | m | `grethe.m` 264 `r_P`: x from the drawing, z the hub depth | estimate |
| `speed_floor` | 0.1 | m/s | `grethe.m` 391 | numerical floor |
| `max_thrust` | 511.5 | N | `grethe.m` 262 `T_max`: the Torqeedo data sheet, maximum thrust | document |
| `delivered_power` | 1120.0 | W | `grethe.m` 263 `P_prop`: the data sheet, power delivered to the water | document |
| `leg_area` | 0.04 | m^2 | `grethe.m` 269 `A_leg` | estimate (a construction, the leg force) |
| `leg_aspect_ratio` | 1.5 | 1 | `grethe.m` 270 `AR_leg` | estimate (a construction) |

Not in the file and not here: a thrust or steering delay of 3-6 s that the logs show before the yaw responds (Y16 O1; the
file's lags cannot make it), a thrust floor, reverse thrust in transit (the logs say forward only, E-154).
