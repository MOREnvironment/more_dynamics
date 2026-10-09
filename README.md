# more_dynamics
Repository containing more dynamical models

## Plugin types

`more_dynamics/plugin_types/vehicle_model.capnp` (Luka's) defines the types: `VehicleModel3D`, `ForceProducer`, `HydrostaticsModel` and `HydrodynamicsModel`. Every component is a plugin class whose `graph()` returns the `CasadyPayload`; each input, state and output is named in the payload's descriptions. A matrix crosses the payload as its column-major `vec` (a 6x6 matrix is 36 values).

## Layout: models are functions, plugins combine them

`models/` holds CasADi functions with named inputs and outputs, no numbers fixed, no rpp: one folder per **physics**, never per vehicle. `plugins/` holds the plugin classes grouped by component type, as Luka's. A plugin takes one or more models, declares their numbers as `PARAMETERS` with cited defaults (so rpp-orchestrator shows and edits them) and freezes them in `initialize`.

```
more_dynamics/
  models/
    rigid_body/  added_mass/  coriolis/  restoring/  damping/  cross_flow/  lift_drag/  surge_resistance/
    current/  site/  hull_form/  equation_of_motion/
    fin/  servo/  propeller/  thruster/  outboard/  shared/
    hull_vessel/  hydrostatics/  hydrodynamics/        (Luka's, left as they are)
  plugins/
    vehicle_models/     spheroid_auv.py  monohull.py  catamaran.py        (+ Luka's hull_vessel.py)
    hydrostatics/       submerged_restoring.py  surface_restoring.py      (+ Luka's linear_surface_hydrostatics.py)
    hydrodynamics/      auv_hull_loads.py  surface_hull_loads.py          (+ Luka's linear_surface / crossflow_surface)
    force_producers/    fin.py  fin_pairs_deflection_only.py  propeller.py  prescribed_wrench.py  (+ Luka's jet_nozzle.py, thruster.py)
    shared/             payload_io.py  vehicle_graph.py
  plugin_types/         vehicle_model.capnp
tests/                  mirrors the model folders; vehicle tests in tests/vehicle_models/
scripts/vehicles/vehicle_simulation.py
.rppws/                 Luka's parts and the vehicles REMUS 100 and Otter
```

Each test folder keeps its frozen reference data beside it in `data/` (the MATLAB generators, their outputs and a `SOURCE.md` naming the source of every file): `tests/{rigid_body,restoring,hydrodynamic_loads}/data/`, `tests/force_producers/data/`, `tests/force_producers/fin/data/`, `tests/vehicle_models/data/{remus100,otter}/`. References that need MSS skip without `MSS_DIR`; the frozen files run everywhere.

## Vehicle types

A vehicle type is one plugin in the form of `HullVessel`: its own `PARAMETERS` (geometry, mass, centres, forms), the rigid body, the added mass, the Coriolis terms and the current computed inside it from the model functions, and the slots `hydrostatics`, `hydrodynamics` and `actuators` (a list of force producers).

| Type | Hull | Defaults |
|---|---|---|
| `SpheroidAuv` | prolate spheroid, mass from the body density, Lamb added mass | REMUS 100 (Prestero 2001; MSS `remus100.m`) |
| `Catamaran` | two pontoons, a hull mass and a point payload, scaled added-mass derivatives | Otter (MSS `otter.m`) |
| `Monohull` | one hull, the same rigid body and added-mass forms as the catamaran | Grethe: document values and **marked estimates** (the comment of each default says which); none is a measurement |

The graph is the equation of motion of a craft in a current over the children's payloads (`plugins/shared/vehicle_graph.py`): `nu_r = nu - nu_c`, `nu_dot = nu_c_dot + M^-1 (tau - C nu_r)`, `eta_dot = J(eta) nu`, with `tau` the sum of the signed generalized forces of the hydrostatics, the hydrodynamics and every actuator. State `[eta; nu]` after the actuators' states (NED position, ZYX Euler angles, BODY FRD velocity, SI units). A child's input is connected **by the name and size in its description** to a quantity the vehicle computed (`pose`, `velocity`, `relative_velocity`, `mass_matrix`, `water_density`, `gravity`, `length`, ...) or to an output of the hydrostatics (`weight`, `restoring_matrix`, `draft`, ...); a missing or mis-sized one is refused, naming the slot and the child. An actuator input that nothing feeds is a command; the vehicle's inputs are the current values first (when `current_form` is not `none`), then the commands in list order (`getInputDescriptions()`).

### Choosing options

Every option carries its fidelity, source, assumptions and the measurement that raises it; each is documented in its plugin's docstring. Pick each with document values first, after a quay day, and after identification from logs.

| Option | Values | Fidelity | Raised by |
|---|---|---|---|
| `coriolis_form` | `kirchhoff_full` (physics, default); `munk_couplings_removed` (MSS shortcut, `remus100.m` 207-210, a comparison form) | physics form above the shortcut | a frequency-dependent added-mass matrix |
| `current_form` | `none`; `yaw_rate_terms` (MSS `remus100.m` 118-122); `full_rotation_rate` (MSS `otter.m` 113-118) | none < MSS forms | a measured current profile; the form with the full attitude waits on its port |
| `added_mass_form` | `lamb_spheroid` (`SpheroidAuv`); `scaled_derivatives` (`Catamaran`, `Monohull`) | ideal fluid / printed coefficients | a Capytaine or identified matrix |
| `buoyancy_method` (`SubmergedRestoring`) | `neutral` (`B = W`); `from_volume` (`B = rho g nabla`); `given` (`B` as a value) | neutral < volume < measured reserve buoyancy | weigh the vehicle in water at the quay, then a floating trim test |
| `hull_count` (`SurfaceRestoring`) | `1`, `2` | metacentric restoring at the equilibrium draft | a measured draft and the waterplane from the hull lines |
| `surge_resistance` (`SurfaceHullLoads`) | `none`; `ittc` (needs a `wetted_surface` the vehicle gives) | physics form with a printed form factor | a measured resistance curve |

A positively buoyant vehicle (Marie) uses `buoyancy_method = "given"` or `"from_volume"` with a centre of buoyancy apart from the centre of gravity: the net force is up and the trim moment follows; her diving controller and fins carry it.

### Building one

In rpp-orchestrator: the vehicles are parts of this library's `.rppws` (`REMUS 100`, `Otter`), assigned to the configurations `remus100` and `otter` of the script `vehicle_simulation`; make a new vehicle by adding a component of a vehicle type and one child per slot, and edit its parameter file. Run one:

    python scripts/vehicles/vehicle_simulation.py . --rpp-configuration remus100 --surge-speed 1.5 --duration 20

From Python, as Luka's and Matko's tests do, with a `ComponentContext` (the registry must have `more_dynamics`: `rpp library register ./more_dynamics --link`, `workspace/` as the working folder):

```python
from rpp_py.data_manager import DataManager
DataManager()                                   # before the generated plugin types are imported
from rpp_py.context import ComponentContext
from rpp_py.parameter_handler import ParameterHandler
from more_dynamics.plugins.vehicle_models.spheroid_auv import SpheroidAuv
from more_dynamics.plugins.hydrostatics.submerged_restoring import SubmergedRestoring
from more_dynamics.plugins.hydrodynamics.auv_hull_loads import AuvHullLoads
from more_dynamics.plugins.force_producers.fin_pairs_deflection_only import FinPairsDeflectionOnly
from more_dynamics.plugins.force_producers.propeller import Propeller

def context(cls, parameters=None, children=None):
    resolved = ParameterHandler.resolve_params(cls.PARAMETERS, {})
    resolved.params.update(parameters or {})
    return ComponentContext(instance=cls(), params=resolved, subcomponents=children or {}, spec=getattr(cls, "COMPONENTS", {}))

remus = context(SpheroidAuv, {"current_form": "yaw_rate_terms"}, {
    "hydrostatics": [context(SubmergedRestoring)], "hydrodynamics": [context(AuvHullLoads)],
    "actuators": [context(FinPairsDeflectionOnly), context(Propeller)]})
remus.initialize()
vehicle = remus.get_instance()            # vehicle.graph() is one CasADi graph; vehicle.signals(x, u) reads any named quantity
```

### Swapping a force producer: same vehicle, two sets

Only the `actuators` children differ between the REMUS fin pairs (`FinPairsDeflectionOnly`, the rudder and stern-plane pair in the deflection-only form) and two `Fin` plugins (an ideal servo, translational inflow, no flow angle, no interference, a quadratic-drag section): on the frozen MSS rows they give the same wrench (`tests/force_producers/fin/test_fin_plugins.py`). They part the moment one `Fin` setting changes (a lag servo instead of the ideal one, say): change that child's parameters and build again, no code changes.

## Vehicles in `.rppws`

| Name | Type | Configuration | Sources |
|---|---|---|---|
| REMUS 100 | `SpheroidAuv` with `SubmergedRestoring`, `AuvHullLoads`, `FinPairsDeflectionOnly`, `Propeller` | `remus100` | Prestero (2001) and MSS `remus100.m` @ `cc07579`, the line of every value in the comment beside it in `params/parameters.py` |
| Otter | `Catamaran` with `SurfaceRestoring`, `SurfaceHullLoads` | `otter` | MSS `otter.m` @ `cc07579`, the same |

The gates run the same plugins through `tests/vehicle_models/` (REMUS derivative and 17 trajectories against MATLAB running MSS, the Otter block by block). rpp rewrites a parameter file without comments when it is saved from the editor's own writer; "Open Parameter File" opens the text for editing.

## Choosing fin and servo options

The `Fin` plugin takes a servo, an inflow, a flow angle, an interference and a section as selectors, their numbers as parameters. Pick each one with document values first; the measurement or identification route raises it to the next fidelity level once that measurement exists.

### Servo (`servo_dynamics`, `servo_rate_limit`, `servo_angle_limit`; model `models/servo/servo.py`)

One device, three switches, not one file per form: `servo_dynamics` (`"none"` / `"first_order_lag"`), `servo_rate_limit` (`False` / `True`, needs the lag's own state), `servo_angle_limit` (`"none"` / `"on_command"` / `"on_output"`). An invalid combination (a rate limit with no lag; an output-side limit with no lag) is refused, naming the switch.

| Setting (`dynamics` / `rate_limit` / `angle_limit`) | Fidelity | Source | Assumptions | Raised to the next level by |
|---|---|---|---|---|
| `none` / `False` / `on_command` | lowest — instantaneous angle limit, no dynamics | published vehicle model | the servo is fast enough to ignore | a step response showing a measurable rate or lag |
| `none` / `False` / `none` | no limits at all | — | a placeholder for an unlimited actuator, not a vehicle's servo | any of the limited settings, once a limit is known |
| `first_order_lag` / `True` / `on_command` | default — amplitude limit on the command, first-order lag, rate limit | published steering-machine model | the lag's time constant is either given directly or derived from a proportional band; the command enters unscaled | a quay-day step test giving the time constant, rate limit and angle limit directly |
| `first_order_lag` / `True` / `on_output` | comparison setting beside the default — amplitude limit on the integrator's own state instead of on the command | published anti-windup case study | equals the default only while the command stays inside the angle limit and the state starts unwound; diverges after a saturating reversal | not a higher-fidelity route on its own — a comparison setting kept beside the default |

Marie's fins and Grethe's steering: `first_order_lag` / `True` / `on_command`, numbers from a quay-day step test and her logs.

### Fin forms (`models/fin/`)

| Selector | Form | Fidelity | Source | Assumptions | Raised to the next level by |
|---|---|---|---|---|---|
| `inflow` | `translational` | lowest | published vehicle model | the fin sees the vehicle's translation, uncorrected for its own rotation | the `rigid_point` form once a rotation rate matters |
| `inflow` | `rigid_point` | default | published vehicle model | the velocity at a point fixed to the body | the flow measured at the fin |
| `flow_angle` | `none` | MSS form | published simulator model | the flow along the chord | `small_angle` |
| `flow_angle` | `small_angle` | default | published vehicle model | linear ratio, valid away from hover | an `exact` large-angle form once hover-speed manoeuvres are simulated |
| `interference` | `none` | default (no interference modelled) | — | one fin's flow does not affect another's | the slender-body factor once its source is read (`slender_body` refuses until then) |
| `section` | `quadratic_drag` | MSS shortcut — flagged, not the default when a physics form exists | published simulator model | no zero-lift drag, no induced drag; registered as a shortcut | `linear_section` with an identified zero-lift drag |
| `section` | `linear_section` | physics form | published vehicle model | small-angle lift and a fixed zero-lift drag | a lifting-line induced-drag form once its source is read (`lifting_line` refuses until then) |

Sources behind the two tables: the servo settings follow MSS `remus100.m` 109, 113-114 (Fossen, MSS, MIT, `cc07579`; angle limit on the command), Murray-Smith, D. J. (2016), *Inverse simulation methods applied to investigations of actuator nonlinearities in ship steering*, Simulation Notes Europe 26(4), 245-256, pp. 246-247 (steering machine: time constant 3 s, ±35°, ±7 and ±10 °/s) and Sarhadi, P. (2026), *Simple yet effective anti-windup techniques for amplitude and rate saturation: an AUV case study*, arXiv:2601.01302v2, Fig. 4, p. 4 (time constant 0.1 s, ±20°, ±30 °/s; the amplitude limit on the output); the fin forms follow Prestero, T. (2001), *Verification of a six-degree of freedom simulation model for the REMUS autonomous underwater vehicle*, MIT/WHOI MSc thesis, eqs. 4.37 (linear section), 4.40 (rigid-point inflow), 4.41-4.43 (small flow angle), pp. 31-33, and MSS `remus100.m` 234-245 (translational inflow, quadratic-drag section). Each plugin's module docstring carries the full citation and the line of every value.
