# REMUS 100 (estimated)

A named part of `more_dynamics::SpheroidAuv`, a sibling of "REMUS 100 (published reference)" with the forms the published model leaves out. It is an **estimate**, not an identification of the vehicle: the hull, mass, site and hull-load time constants are the one-value REMUS 100 set of the published reference part; the forms and the fin numbers are Prestero (2001).

## Forms (this part's choices)

| Selector | Value | Why |
|---|---|---|
| `coriolis_form` | `kirchhoff_full` | every term of the added-mass Coriolis matrix, the Munk moment of a slender body included (Fossen 2011, eqs. 6.43 and 6.52, pp. 120-121); kept only with the body lift and the geometric fin slope below, and only with the straight-line stability gate green (`tests/vehicle_models/test_spheroid_auv_fidelity_options.py`) |
| `body_lift_form` (hull loads) | `point_force_both_planes` | Prestero 2001, eqs. 4.27-4.36, pp. 30-31; the moment has the sign of the cross product |
| `body_lift_slope_method` (hull loads) | `slender_body_hoerner` | Prestero 2001, eqs. 4.30-4.32, p. 30; l/d = 8.42 is inside the stated 6.7 to 10 |
| `pressure_centre_fraction` (hull loads) | 0.65 | Prestero 2001, eq. 4.35, p. 31 (Hoerner's range 0.6 to 0.7) |
| `hull_drag_form` (hull loads) | `wing_drag` | the drag half of the wing block, so only the lift changes against the published reference part |
| `lift_law` (fins) | `finite_wing`, `low_aspect_ratio` | the lift slope from the fin geometry, Prestero 2001, eqs. 4.38-4.39, p. 31: 3.12 per radian per fin |
| `shaft_torque_scale` (propulsor) | 0.35 | **estimate**: the roll reaction of the propeller, calibrated to a steady roll of -5.3 deg in the generic-models set, not re-calibrated here |

The actuators are linked named parts: the tail rudder pair (`REMUS 100 (estimated) rudder fin`, linked twice), the stern-plane pair (`REMUS 100 (estimated) stern plane fin`, linked twice) and the propeller (`REMUS 100 (estimated) propulsor`): five inputs in this order, four fin deflections (rad) and the shaft speed (rpm, `drive_command`).

`hull_frontal_drag_coefficient` 0.42 of the hull loads is the axial drag coefficient of Allen et al. (2000), cited second-hand, not read; it is read only by `hull_drag_form = quadratic_axial`, which this part does not select.

The pressure-centre fraction of the reference gates at 0.60 and 0.70 bounds Hoerner's range as a sensitivity.
