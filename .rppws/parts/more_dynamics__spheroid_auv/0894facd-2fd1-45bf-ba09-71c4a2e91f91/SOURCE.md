# REMUS 100 (MSS, ocean current)

A named part of `more_dynamics::SpheroidAuv`: **a copy of the part "REMUS 100 (MSS)" with one change, `current_form`
`yaw_rate_terms` instead of `none`.** The vehicle takes the ocean current as three inputs (`current_speed`,
`current_direction`, `current_vertical_speed`, before the actuator inputs) and forms the relative velocity and its
derivative the way `CRAFT/AUV/models/remus100.m` 118-123 does (MSS, MIT, @ `ed2ce7b`): the horizontal current rotated by
the yaw only, the yaw rate in its derivative, the vertical current added to the heave. This is exact for a level craft
and a shortcut otherwise (`more_dynamics/models/current`); the default of the vehicle type stays `full_attitude`.

Every other value and form is the part REMUS 100 (MSS)'s (`../83497a04-.../SOURCE.md`) and it links the same five
actuators. A test asserts that every parameter file of the two parts is equal except `current_form`, so the two cannot
drift apart (`tests/vehicle_models/test_current_parts.py`).
