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
| `site_form` | `given` | the part's own `Site` subcomponent (`site` slot) | no site derivation for Grethe |
| `wetted_surface_method` | `regression_table` | hydrostatics (SurfaceRestoring) | A-64 Section 3.3 R11: the published RA14 hydrostatic wetted-surface regression, one hull |
| `manoeuvring_damping` | `linear_coupled` | hydrodynamics (SurfaceHullLoads) | A-64 Section 3.4 G0 set: `Y_v`, `N_r` from the time constants, `Y_r = N_v = 0` unidentified (E-101 Q2 a) |
| `surge_resistance` | `ittc_residual` | hydrodynamics (SurfaceHullLoads) | E-101 Q1 a: ITTC friction plus the RA14 residual term |
| `cross_flow` | `hoerner_strips` | hydrodynamics (SurfaceHullLoads) | Hoerner cross-flow strips, A-64 Section 3.4 |

Every other value is set in the part itself (the plugins carry no vehicle numbers: a value the
part leaves out makes the build fail naming it) -- see `more_dynamics/DEFAULTS.md` for value,
unit, source and kind, and `tests/vehicle_models/data/grethe/grethe_parameters.json`
for the full set, including values not yet wired as a library option (`waterline_length`,
`waterline_beam`; `SOURCE.md` of that directory).

## Open

- `block_coefficient_method` stays `given` (the selector's default): `hull_block_coefficient` is
  a vehicle coupling, not yet `computed` from a measured draft for this part.
- `yaw_damping_nonlinearity` (hydrodynamics) is not set: the selected `linear_coupled`
  damping does not read it (the time-constants form does; the reference data marks the
  Otter's 10.0 s/rad "published (other vehicle)", unidentified for Grethe).
- `waterline_length`/`waterline_beam` are not set: no option reads them yet. The drawing
  estimates (4.63 m, 1.97 m, D002 p1) are in the reference data file; the old "0 = not given"
  convention is retired (a zero is a value, never "not given").

## Values of the hull loads that are not from a document (owner, 2026-10-12)

The three numbers the selected `linear_coupled` damping reads and Grethe's data file does not
carry were, before the placeholders, silently the Otter's (thrust 239.364 N, top speed
3.0864 m/s, damping ratios [0.3, 0.2, 0.4]). They are now Grethe's own, each labelled:

| Value | Set to | Kind | Source |
|---|---|---|---|
| `max_speed` | 2.40 m/s through water | **fitted** from her 2026 control-station logs (steady windows with the speed order above the speed, current removed per day); logged band 2.32-2.47 m/s | A-74 section 1, `a74_top_speed.py`; replace by the straight runs at full throttle of the sea trials (X-9 G14) |
| `max_forward_thrust` | 0 N | **decision**: the surge resistance of this part is the ITTC and residual forms alone; a thrust here would count a second resistance. A coast-down-fitted resistance option follows (held-out error 0.18 m/s) | A-74 Q1, owner E-135 |
| `damping_ratios` roll | 0.05 | **fitted** from her logs, band 0.03-0.08 | A-74 section 3, owner E-135 |
| `damping_ratios` heave 0.3, pitch 0.5 | the values of the earlier test configuration | **approved guess 2026-10-12 (owner E-135), no source**, until the quay free decays (X-9 G12) | owner E-135 |
