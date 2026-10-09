# Grethe

A named part of `more_dynamics::Monohull` (single-hull surface craft). NTNU Mariner 5 USV,
Maritime Robotics build 40402, Pioner 17 ft hull.

Every number in `params/parameters.py` and in the `hydrostatics`/`hydrodynamics`
subcomponents' own `params/parameters.py` is read from
`tests/vehicle_models/data/grethe/grethe_parameters.json` -- the one hand-edited source
(R10, E-104), which carries its own `kind` and `place` field per value (document / estimate
/ derived / selector / owner decision / gap). Sources: `agents-more/60_working/2026-10-07_A43_vehicle_parameter_dossiers.md`
Part 2 (M1-M50), `agents-more/20_sources/monohull_structure_for_grethe.md` (A-64),
`agents-more/80_owner/inbox/E-101_a64_grethe_monohull_choices_remus_density_at_repin.md`.
No August 2026 field data is used (safeguard rule 20).

## Fidelity selectors (this part's choices, E-101)

| Selector | Value | Where | Why |
|---|---|---|---|
| `added_mass_mass_basis` | `displaced` | root (Monohull) | E-101 Q5 a: a vehicle with a payload displaces more water; the added mass follows |
| `current_form` | `none` | root (Monohull) | still water: the reference data holds current fixed |
| `site_form` | `given` | root (Monohull) | plugin default; no site derivation for Grethe |
| `wetted_surface_method` | `regression_table` | hydrostatics (SurfaceRestoring) | A-64 Section 3.3 R11: the published RA14 hydrostatic wetted-surface regression, one hull |
| `manoeuvring_damping` | `linear_coupled` | hydrodynamics (SurfaceHullLoads) | A-64 Section 3.4 G0 set: `Y_v`, `N_r` from the time constants, `Y_r = N_v = 0` unidentified (E-101 Q2 a) |
| `surge_resistance` | `ittc_residual` | hydrodynamics (SurfaceHullLoads) | E-101 Q1 a: ITTC friction plus the RA14 residual term |
| `cross_flow` | `hoerner_strips` | hydrodynamics (SurfaceHullLoads) | plugin default |

Every other declared parameter is Grethe's own cited default -- see `more_dynamics/DEFAULTS.md`
for value, unit, source and kind, and `tests/vehicle_models/data/grethe/grethe_parameters.json`
for the full set, including values not yet wired as a library option (`waterline_length`,
`waterline_beam`; `SOURCE.md` of that directory).

## Open

- `block_coefficient_method` stays `given` (the plugin default): `hull_block_coefficient` is
  a vehicle coupling, not yet `computed` from a measured draft for this part.
- `yaw_damping_nonlinearity` (hydrodynamics) is left at the plugin's own default (Otter's
  value, 10.0 s/rad, `otter.m` 99-101, 202-240): the reference data marks it
  "published (other vehicle)", unidentified for Grethe.
- `waterline_length`/`waterline_beam` are drawing estimates (4.63 m, 1.97 m) not yet wired
  into any plugin option -- a later job.
