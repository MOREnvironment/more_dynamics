# Otter (MSS)

A named part of `more_dynamics::Catamaran` (twin-hull surface craft).

Every number in `params/parameters.py` is Fossen's MSS `CRAFT/USV/models/otter.m` (MIT,
revision `99fdf4d`) as cited in the comment beside it, plus `generate_catamaran_mss.m`
(payload mass and position) and `LIBRARY/modeling/cylinderDrag.m` (kinematic viscosity).

## Fidelity selectors (this part's choices)

| Selector | Value | Why |
|---|---|---|
| `site_form` | `given` | the part's own `Site` subcomponent (`site` slot): MSS gives `gravity`/`water_density` directly, no site derivation |
| `added_mass_form` | `scaled_derivatives_loaded_surge` | MSS's scaled-derivative added mass, surge term from the loaded mass `m + mP` (`otter.m` 162, `Xudot = -addedMassSurge(m+mP,L,rho)`; F-MSS-22 fixed upstream) |
| `coriolis_form` | `kirchhoff_full` | every added-mass Coriolis term kept |
| `current_form` | `none` | this part is gated in still water |

The part's hull-loads subcomponent selects `manoeuvring_damping = time_constants_surge_blend`
(`otter.m` 206-217, 244-251: surge force `Xu exp(-3|u_r|) u_r + Xuu |u_r| u_r`, `Xu = -M11/T_surge`
with `T_surge` 5 s, `Xuu = -24.4 g / Umax^2`; `T_sway` 0.5 s). It is bit-for-bit the MSS form, including
the payload `eta_0` shift of `otter.m` 265-266, 269 that the registry lists as F-MSS-24 and F-MSS-31 (still present
at `99fdf4d`); this part does not correct them.

Every other value (hull geometry, mass, radii of gyration, added-mass coefficients, and the
hydrostatics, hull-loads and site values of the subcomponents) is set in the part itself with
its source in the comment beside it (the plugins carry no vehicle numbers: a value the part
leaves out makes the build fail naming it) -- see `more_dynamics/DEFAULTS.md` for value, unit,
source line and kind.
