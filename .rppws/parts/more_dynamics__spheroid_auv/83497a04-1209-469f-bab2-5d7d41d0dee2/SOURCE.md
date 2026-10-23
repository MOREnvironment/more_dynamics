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
| `site_form` | `latitude` | the part's own `Site` subcomponent (`site` slot): WGS-84 normal gravity at the MSS site; its density 1025 (mssConstants.m 12) and viscosity 1e-06 are held there too |
| `added_mass_form` | `lamb_spheroid` | ideal-fluid prolate spheroid |
| `coriolis_form` | `munk_couplings_removed` | `remus100.m` 207-210: MSS omits the Munk-moment couplings of the added-mass Coriolis (its reduced model keeps no measured manoeuvring derivatives to balance them); with every term kept (`kirchhoff_full`) this part spins open-loop at zero rudder (design note on the REMUS added-mass options, 2026-10-12) |
| `current_form` | `none` | this part is gated in still water |
| `body_density_method` | `computed` | density derived from `body_mass` on the spheroid's own volume |

Every other value (geometry, mass, roll added-inertia ratio, and the hull-loads, hydrostatics,
fin-pair, propeller and site values of the subcomponents) is set in the part itself, with its
source line in the comment beside it: `more_dynamics/DEFAULTS.md` gives value, unit, source
line and kind for this part. The plugins carry no vehicle numbers (unset placeholders): a
value this part leaves out makes the build fail naming it. Three values were still read from
the former plugin defaults before the placeholders and are now written into the part:
`body_mass` 31.9 kg (`remus100.m`:3, document), and in `REMUS 100 hull loads`
`planform_fraction` 0.7 (`remus100.m`:133, document) and `hull_frontal_drag_coefficient`
0.42 (`remus100.m`:142, Cd of Allen et al. 2000 as MSS cites it, not read).

## Higher-fidelity variants (not yet landed, E-106)

The owner's decision (E-106, 2026-10-09) is to keep this part as the MSS baseline and add
sibling parts under other names as each fidelity option lands and is gated: e.g.
"REMUS 100 (Prestero)" (Prestero's own published coefficients and forms, once U6f/U6g land
the eight fidelity options), later "REMUS 100 (identified)". This part is never renamed to
absorb those -- a new part is added beside it.
