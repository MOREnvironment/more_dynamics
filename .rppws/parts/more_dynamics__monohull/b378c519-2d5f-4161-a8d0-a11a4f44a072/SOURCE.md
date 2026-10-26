# Grethe (yaw-unstable stress test)

A named part of `more_dynamics::Monohull`: **a copy of the part "Grethe" with one change, for a controller stress test, not
a vehicle estimate.** Its hull loads take the earlier yaw damping of the Grethe part, `time_constants` [8, 16] s (T_sway,
T_yaw) and `sway_yaw_coupling` [0, 0] (Y_r = N_v = 0), instead of the model file's [8, 4] s and [0, 2500]. With the Munk moment
kept and no weather-cock term the open loop is unstable in yaw: a real pole of about +0.33 rad/s at 1.5 m/s with no steering
(a test, `test_the_part_is_stable_in_yaw_and_the_stress_test_part_is_not`), growing with speed, the plant on which the
first closed-loop attempt did not arrive (owner E-151 G1: kept as a second part, labelled so). It is meant for the exercise
of a heading controller that must stabilise an unstable plant; no value of it is identified, and the real boat is logged
turning at about 4 deg/s in missions, which this plant cannot do open loop.

Every other value and form is the part Grethe's (`../83646192-.../SOURCE.md`, the model file `grethe.m` @ `679a6e8`), and it links
the same outboard (`Grethe outboard`). A test asserts that every parameter file of the two parts is equal except those two
numbers of the hull loads, so the two cannot drift apart.
