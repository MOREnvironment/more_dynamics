# Parameter sets and frozen references — the research AUV Marie, one named part "Marie (estimated)"

**One named part.** Marie is data only: one parameter set of the generic table hull `SubmergedBemAuv` with `SubmergedRestoring` (`buoyancy_method = immersed_body_of_revolution`), `BemHullLoads` and four `Fin`, named by where its data come from — "Marie (estimated)" today; later "Marie (measured)" after a weighing and "Marie (identified)" from the logs. The floating restoring is exact at the surface and fully submerged, so the dive and the surface phase are states of the same part.

## Files

| File | What it holds |
|---|---|
| `marie_parameters.json` | the named part: `part` (its name) and `parameters` grouped by the component that declares them (`vehicle`, `hydrostatics`, `hydrodynamics`, `site`, `fin_common`, `fins`); each entry is a `value` and its `kind` |
| `scaled_torpedo_parameters.json` | **not a vehicle**: the second geometry of the generality check, the part's body scaled by 1.6 with Froude similitude (lengths s, masses s^3, inertias s^5, times s^0.5, frequencies s^-0.5, the 6 x 6 tables as `s^(3 + a_i + a_j)` for the added mass and `s^(2.5 + a_i + a_j)` for the damping, `a` = 0 for a translation, 1 for a rotation) |
| `marie_floating_reference.csv`, `scaled_torpedo_floating_reference.csv` | floating equilibria (depth of the origin below the calm surface, bow-up pitch, immersed volume and centre of buoyancy) of the declared part, the part with the centre of gravity moved aft by 10 and 20 mm (times the scale), and 4 reserves x 4 longitudinal offsets (the reserve set through the mass) |
| `generate_marie_references.py` | writes the second geometry and the two floating references from `marie_parameters.json` (numpy and scipy only, a second implementation: adaptive quadrature of the exact circle segments, eq. 4.5 of Fossen 2011, p. 60) |

## Kinds

* `published`: a value the vehicle's operator publishes (length, mass in air, hull diameter = the section height).
* `fitted to field data`: identified from the vehicle's logs (none in this set yet).
* `estimate`: everything else with a number. The frequency tables, the section breadth, the station table of the hull, the centres and the inertia belong to an **idealised hull** of the vehicle's length and diameter (a bare hull, uniform-cylinder radii of gyration, no fins, antenna or propeller): an example body, not a measurement of the vehicle. The displaced volume is an interim reserve buoyancy of 0.57 % of the weight. The damping time constants and ratios, the speed factors and the lift and drag coefficients are tuned values of a manoeuvring model. The fin area, span and section efficiency belong to another vehicle of the same size class.
* `design choice`: a selector or a numerical setting (forms, the spectrum that weights the tables, the number of slices of the station sums).

The reserve buoyancy controls the floating pitch strongly (at 10 mm aft of the centre of buoyancy: 23.6 deg bow-up at 0.57 %, 7.5 deg at 1 %, 0.8 deg at 2 %); weighing the vehicle in water is the route to a measured value.

## Regenerate

`python3 generate_marie_references.py` in this folder. `marie_parameters.json` itself is written from the vehicle's parameter files outside this library; a changed value needs the references regenerated.
