# REMUS 100 (MSS, quaternion)

A named part of `more_dynamics::SpheroidAuv`: **a copy of the part "REMUS 100 (MSS, ocean current)" with one change,
`kinematics_form` `unit_quaternion` instead of the default `euler`.** The vehicle carries its attitude as the unit
quaternion `[eta, eps1, eps2, eps3]`: 13 states in the order position, quaternion, velocity (the Euler form has 12:
position, Euler angles, velocity), as the 13-state form of `CRAFT/AUV/models/remus100.m` 225-227, 258-259 (MSS, MIT,
@ `ed2ce7b`), whose state is ordered velocity, position, quaternion. The vehicle's output is still position, Euler
angles and velocity (12 values), so guidance, control and sensors read it as before. The quaternion is renormalised
after every integration step of the vehicle's own step.

MSS's own run in this form (the dialog's default, `SIMremus100.m` 646-660) integrates with forward Euler on the
velocities, backward Euler on the position and the exact quaternion exponential, `SIMremus100.m` 349-356: that scheme
is the test integrator of the gate (`more_mission/tests/mission/mss_chains/test_mss_chain_vehicles.py`), not part of
this part. A test asserts that every parameter file of this part is equal to "REMUS 100 (MSS, ocean current)"'s except
`kinematics_form` (`tests/vehicle_models/test_current_parts.py`).
