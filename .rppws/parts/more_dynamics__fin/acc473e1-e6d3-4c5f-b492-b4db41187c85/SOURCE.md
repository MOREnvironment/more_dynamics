# REMUS 100 (MSS) stern plane fin

A named part of `more_dynamics::Fin` for the REMUS 100 (MSS) vehicle: one of the two fins of the stern plane (remus100.m 188-190: CL_delta_s 0.7, A_s = 2 S_fin: two fins of S_fin). The vehicle links it
twice (a tail pair), so the vehicle has four fins, each with its own command; MSS's single command per pair
is the same value given to both. Every number is `CRAFT/AUV/models/remus100.m` (MIT, revision `ed2ce7b`),
lines cited beside it in `params/parameters.py`; the forms are the deflection-only form of `remus100.m`
188-190, 240, 246 (inflow `translational`, flow angle `none`, section `quadratic_drag`, an ideal servo whose only
effect is the saturation of the command).

| Parameter | Value | Unit | Source | Kind |
|---|---|---|---|---|
| `fin_position` | [-0.8, 0, 0] | m | `remus100.m`:185, 190, x = -a | derived (a = L_auv / 2, one geometry; typed here because `Fin` has no option that takes it from the hull, so a change of the hull length does not move it: follow-up) |
| `chord_axis` | [1, 0, 0] | 1 | BODY FRD | convention |
| `lift_axis` | [0.0, 0.0, -1.0] | 1 | sign of the force of a positive deflection, `remus100.m`:246 | convention |
| `fin_area` | 0.00665 | m^2 | `remus100.m`:180 S_fin | document |
| `lift_slope` | 0.7 | 1/rad | `remus100.m`:188 | document |
| `max_deflection` | 0.3490658503988659 | rad | `remus100.m`:110 delta_max = deg2rad(20) | document |
| `time_constant`, `max_rate`, `zero_lift_drag` | unset | | not read by these selectors | not applicable |

The MSS pair equations (`remus100.m` 239-246) are matched to rounding by two fins of area S_fin each with
the same deflection (`tests/vehicle_models/test_rppws_runs.py`). Note: the MSS reference used here is the
consistent one (`a = L_auv / 2`, one geometry) of the frozen runs, as in the REMUS 100 (MSS) vehicle.
