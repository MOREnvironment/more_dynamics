# Otter (MSS)

A named part of `more_dynamics::Catamaran` (twin-hull surface craft).

Every number in `params/parameters.py` is Fossen's MSS `CRAFT/USV/models/otter.m` (MIT,
revision `cc07579`) as cited in the comment beside it, plus `generate_catamaran_mss.m`
(payload mass and position) and `LIBRARY/modeling/cylinderDrag.m` (kinematic viscosity).

## Fidelity selectors (this part's choices)

| Selector | Value | Why |
|---|---|---|
| `site_form` | `given` | the part's own `Site` subcomponent (`site` slot): MSS gives `gravity`/`water_density` directly, no site derivation |
| `added_mass_form` | `scaled_derivatives` | MSS's scaled-derivative added mass |
| `coriolis_form` | `kirchhoff_full` | every added-mass Coriolis term kept |
| `current_form` | `none` | this part is gated in still water |

Every other value (hull geometry, mass, radii of gyration, added-mass coefficients, and the
hydrostatics, hull-loads and site values of the subcomponents) is set in the part itself with
its source in the comment beside it (the plugins carry no vehicle numbers: a value the part
leaves out makes the build fail naming it) -- see `more_dynamics/DEFAULTS.md` for value, unit,
source line and kind.
