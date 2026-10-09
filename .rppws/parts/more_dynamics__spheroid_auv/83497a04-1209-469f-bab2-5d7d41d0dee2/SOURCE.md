# REMUS 100 (MSS)

A named part of `more_dynamics::SpheroidAuv` (torpedo-shaped AUV, prolate-spheroid hull).

Every number in `params/parameters.py` is Fossen's MSS `CRAFT/AUV/models/remus100.m` (MIT,
revision `cc07579`) as cited in the comment beside it, plus `LIBRARY/modeling/spheroid.m`
(body density), `imlay61.m` (added mass) and `INS/functions/gravity.m` (latitude site). The
hull geometry and coefficients trace back to Prestero, T. (2001), *Verification of a
six-degree of freedom simulation model for the REMUS autonomous underwater vehicle*,
MIT/WHOI MSc thesis -- not re-read for this part; every value here is read from MSS's
instrumented workspace, not retyped from the thesis (rule 4).

## Fidelity selectors (this part's choices)

| Selector | Value | Why |
|---|---|---|
| `site_form` | `latitude` | plugin default; WGS-84 normal gravity at the MSS site |
| `added_mass_form` | `lamb_spheroid` | plugin default; ideal-fluid prolate spheroid |
| `coriolis_form` | `kirchhoff_full` | plugin default; every added-mass Coriolis term kept |
| `current_form` | `none` | this part is gated in still water |
| `body_density_method` | `computed` (plugin default) | density derived from `body_mass` on the spheroid's own volume |

Every other declared parameter (geometry, mass, roll added-inertia ratio) is the plugin's
own cited default -- see `more_dynamics/DEFAULTS.md` for value, unit, source line and kind.

## Higher-fidelity variants (not yet landed, E-106)

The owner's decision (E-106, 2026-10-09) is to keep this part as the MSS baseline and add
sibling parts under other names as each fidelity option lands and is gated: e.g.
"REMUS 100 (Prestero)" (Prestero's own published coefficients and forms, once U6f/U6g land
the eight fidelity options), later "REMUS 100 (identified)". This part is never renamed to
absorb those -- a new part is added beside it.
