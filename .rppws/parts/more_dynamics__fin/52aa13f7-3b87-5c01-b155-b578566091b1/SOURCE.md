# REMUS 100 (estimated) stern plane fin

A named part of `more_dynamics::Fin`: one of the two stern plane fins of "REMUS 100 (estimated)". The vehicle links it twice (a tail pair), so the vehicle has four fins, each with its own command; the command of a pair is given to both fins.

| Parameter | Value | Unit | Source | Kind |
|---|---|---|---|---|
| `fin_position` | [-0.8, 0, 0] | m | the fin post at the tail, x = -L/2 of the 1.6 m hull centred on the origin | derived (typed: `Fin` has no option that takes it from the hull) |
| `chord_axis` | [1, 0, 0] | 1 | BODY FRD | convention |
| `lift_axis` | [0.0, 0.0, -1.0] | 1 | a positive deflection lifts along this axis (negative on z: a positive stern-plane deflection pitches the bow down) | convention |
| `fin_area` | 0.00665 | m^2 | Prestero 2001, Table 2.2, p. 17 | published |
| `servo_dynamics`, `servo_rate_limit`, `servo_angle_limit` | 'none', False, 'on_command' | - | the command is the deflection, saturated | selector |
| `max_deflection` | 0.3490658503988659 | rad | 20 deg | published |
| `inflow` | 'rigid_point' | - | Prestero 2001, eq. 4.40, p. 32 | selector |
| `flow_angle` | 'small_angle' | - | Prestero 2001, eqs. 4.41-4.43, pp. 32-33 | selector |
| `interference` | 'none' | - | Prestero 2001, eq. 4.41, p. 32 | selector |
| `section` | 'lift_times_angle' | - | C_L = a alpha, C_D = C_L alpha | selector |
| `lift_law`, `lift_slope_form` | 'finite_wing', 'low_aspect_ratio' | - | Prestero 2001, eqs. 4.38-4.39, p. 31 | selector |
| `fin_span` | 0.0857 | m | Prestero 2001, Table 2.2, p. 17 | published |
| `aspect_ratio_factor` | 2.0 | 1 | Prestero 2001, eq. 4.39, p. 31 | published |
| `section_efficiency` | 0.9 | 1 | Prestero 2001, Table 2.2, p. 17 and eq. 4.38, p. 31 | published |
| `time_constant`, `max_rate`, `lift_slope`, `zero_lift_drag` | unset | | not read by these selectors | not applicable |

The lift slope per fin is 3.1158 per radian, `[1/(2 pi 0.9) + 1/(pi 2 0.0857^2 / 0.00665)]^-1`; Prestero 2001, Table 2.2, p. 17 prints 3.12 (`tests/vehicle_models/test_spheroid_auv_fidelity_options.py`).
