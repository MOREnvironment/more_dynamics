# Otter (MSS)

A named part of `more_dynamics::Catamaran` (twin-hull surface craft).

Every number in `params/parameters.py` is Fossen's MSS `CRAFT/USV/models/otter_legacy.m` (MIT,
revision `6e78e98`) as cited in the comment beside it, plus `generate_catamaran_mss.m`
(payload mass and position) and `LIBRARY/modeling/cylinderDrag.m` (kinematic viscosity).

## Fidelity selectors (this part's choices)

| Selector | Value | Why |
|---|---|---|
| `site_form` | `given` | the part's own `Site` subcomponent (`site` slot): MSS gives `gravity`/`water_density` directly, no site derivation |
| `added_mass_form` | `scaled_derivatives_loaded_surge` | MSS's scaled-derivative added mass, surge term from the loaded mass `m + mP` (`otter_legacy.m` 168, `Xudot = -addedMassSurge(m+mP,L,rho)`; F-MSS-22 fixed upstream) |
| `coriolis_form` | `kirchhoff_full` | every added-mass Coriolis term kept |
| `current_form` | `none` | this part is gated in still water |

The part's hull-loads subcomponent selects `manoeuvring_damping = time_constants_surge_blend`
(`otter_legacy.m` 212-223, 250-257: surge force `Xu exp(-3|u_r|) u_r + Xuu |u_r| u_r`, `Xu = -M11/T_surge`
with `T_surge` 5 s, `Xuu = -24.4 g / Umax^2`; `T_sway` 0.5 s). The part's restoring selects
`payload_trim_form = lever_arm_and_heave_shift`, so the part is bit-for-bit the MSS form, including
the payload `eta_0` shift of `otter_legacy.m` 265-275 that the registry lists as F-MSS-24 and F-MSS-31 (still present
at `6e78e98`: the payload weight in heave a second time, and the kinematics at the shifted attitude); this part
does not correct them. The library's corrected form is the plugin's default, `lever_arm_moment`.

The actuators are the two thrusters of the Otter, linked from the named parts `Otter port thruster` and
`Otter starboard thruster` (`otter_legacy.m` 32-35, 145-152, 236-248): two inputs named `drive_command`, the shaft speed
of the port (left) and of the starboard (right) propeller in rad/s, in this order (a command schedule gives the
name its two values in this order; `n(1)` and `n(2)` of `otter_legacy.m`).

Every other value (hull geometry, mass, radii of gyration, added-mass coefficients, and the
hydrostatics, hull-loads and site values of the subcomponents) is set in the part itself with
its source in the comment beside it (the plugins carry no vehicle numbers: a value the part
leaves out makes the build fail naming it) -- see `more_dynamics/DEFAULTS.md` for value, unit,
source line and kind.
