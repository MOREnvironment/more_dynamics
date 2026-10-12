# REMUS 100 fin

A named part of `more_dynamics::Fin` (one lifting fin: rigid-point inflow, small-angle flow angle,
no interference, quadratic-drag section, first-order servo with rate and angle limits). It is the
reference values the fin tests and a vehicle's fin children are built from; it is not wired into a
vehicle part yet (the REMUS 100 (MSS) vehicle uses the deflection-only fin pair).

| Parameter | Value | Unit | Source | Kind |
|---|---|---|---|---|
| `fin_position` | [-0.638, 0.0, 0.0] | m | Prestero 2001, Table A.5, p. 103 (x_finpost) | document |
| `chord_axis` | [1.0, 0.0, 0.0] | 1 | unit chord direction, BODY FRD | convention |
| `lift_axis` | [0.0, 1.0, 0.0] | 1 | unit positive-lift direction, BODY FRD | convention |
| `fin_area` | 0.00665 | m^2 | Prestero 2001, Table A.5, p. 103 (S_fin) | document |
| `time_constant` | 0.1 | s | Sarhadi 2026, Fig. 4, p. 4: a demanding choice of that study, not a measurement | estimate |
| `max_deflection` | 0.3490658503988659 (20 deg) | rad | Sarhadi 2026, Fig. 4, p. 4 | estimate |
| `max_rate` | 0.5235987755982988 (30 deg/s) | rad/s | Sarhadi 2026, Fig. 4, p. 4 | estimate |
| `lift_slope` | 3.12 | 1/rad | Prestero 2001, Table A.5, p. 103 (c_L_alpha) | document |
| `zero_lift_drag` | 0.0 | 1 | Prestero 2001, eq. 4.37, p. 31 (carries no fin drag; read only by the `linear_section` form) | zero with its reason |

The values were the plugin's former defaults, moved into this part when the plugins became neutral
placeholders; none is new. Prestero 2001 and Sarhadi 2026 are cited as in `DEFAULTS.md`; the values are
as the earlier default carried them (not re-read here). The selectors (`servo_dynamics`,
`servo_rate_limit`, `servo_angle_limit`, `inflow`, `flow_angle`, `interference`, `section`) keep the
plugin's defaults.
