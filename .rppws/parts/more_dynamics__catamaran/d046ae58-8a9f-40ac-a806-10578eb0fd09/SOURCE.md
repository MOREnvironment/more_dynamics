# Otter (MSS, ocean current)

A named part of `more_dynamics::Catamaran`: **a copy of the part "Otter (MSS)" with one change, `current_form`
`full_rotation_rate` instead of `none`.** The vehicle takes the ocean current as three inputs (`current_speed`,
`current_direction`, `current_vertical_speed`, before the two `drive_command` inputs) and forms the relative velocity
and its derivative the way `CRAFT/USV/models/otter.m` 120-125 does (MSS, MIT, @ `ed2ce7b`): the horizontal current
rotated by the yaw, the full rotation rate in its derivative. Exact for a craft level in roll and pitch; the default of
the vehicle type stays `full_attitude`.

Every other value and form is the part Otter (MSS)'s (`../0705552d-.../SOURCE.md`) and it links the same actuators. A
test asserts that every parameter file of the two parts is equal except `current_form`
(`tests/vehicle_models/test_current_parts.py`).
