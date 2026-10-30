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

| Type | Hull | Reference named part (the plugin itself carries no number) |
|---|---|---|
| `SpheroidAuv` | prolate spheroid, mass from the body density, Lamb added mass | REMUS 100 (MSS) (Prestero 2001; MSS `remus100.m`) |
| `Catamaran` | two pontoons, a hull mass and a point payload, scaled added-mass derivatives | Otter (MSS) (MSS `otter_legacy.m`) |
| `Monohull` | one hull, the same rigid body and added-mass forms as the catamaran | Grethe: document values, **marked estimates**, fitted values and one marked approved guess (the comment of each value says which) |

A plugin declares every vehicle, device or site number as unset (an empty list, `unset()` of `more_transformations`); the numbers come from a named part of `.rppws/parts/` when the composition is built, and a value its selected forms read that the part leaves out refuses the build, naming the plugin and every such value with unit and meaning. Only text selectors, switches and three numerical settings keep a default. `DEFAULTS.md` gives the provenance of each reference part's values.

The graph is the equation of motion of a craft in a current over the children's payloads (`plugins/shared/vehicle_graph.py`): `nu_r = nu - nu_c`, `nu_dot = nu_c_dot + M^-1 (tau - C nu_r)`, `eta_dot = J(eta) nu`, with `tau` the sum of the signed generalized forces of the hydrostatics, the hydrodynamics and every actuator. State `[eta; nu]` after the actuators' states (NED position, ZYX Euler angles, BODY FRD velocity, SI units). A child's input is connected **by the name and size in its description** to a quantity the vehicle computed (`pose`, `velocity`, `relative_velocity`, `mass_matrix`, `water_density`, `gravity`, `length`, ...) or to an output of the hydrostatics (`weight`, `restoring_matrix`, `draft`, ...); a missing or mis-sized one is refused, naming the slot and the child. An actuator input that nothing feeds is a command; the vehicle's inputs are the current values first (when `current_form` is not `none`), then the commands in list order (`getInputDescriptions()`).

### Choosing options

Every option carries its fidelity, source, assumptions and the measurement that raises it; each is documented in its plugin's docstring. Pick each with document values first, after a quay day, and after identification from logs.

| Option | Values | Fidelity | Raised by |
|---|---|---|---|
| `coriolis_form` | `kirchhoff_full` (physics, default); `munk_couplings_removed` (MSS shortcut, `remus100.m` 207-210, a comparison form) | physics form above the shortcut | a frequency-dependent added-mass matrix |
| `current_form` | `none`; `yaw_rate_terms` (MSS `remus100.m` 118-122); `full_rotation_rate` (MSS `otter_legacy.m` 119-124) | none < MSS forms | a measured current profile; the form with the full attitude waits on its port |
| `kinematics_form` (`SpheroidAuv`) | `euler` (default, 12 states: position, ZYX Euler angles, velocity); `unit_quaternion` (13 states: position, unit quaternion, velocity; the output stays position, Euler angles, velocity; the quaternion renormalised after every step of the vehicle's own step; Fossen 2011, eq. 2.69; `remus100.m` 225-227) | Euler < quaternion at a pitch near 90 deg, equal elsewhere | none needed: the form changes the attitude state, not the physics. The `Catamaran` has no such form (`otter_legacy.m` 97 is 12-state) |
| `added_mass_form` | `lamb_spheroid` (`SpheroidAuv`); `scaled_derivatives` (`Catamaran`, `Monohull`); `scaled_derivatives_loaded_surge` (`Catamaran`: the surge term from the loaded mass, hull plus payload) | ideal fluid / printed coefficients | a Capytaine or identified matrix |
| `manoeuvring_damping` (`SurfaceHullLoads`) | `time_constants` (surge from the thrust and the top speed); `time_constants_surge_blend` (surge from a time constant, blended with a quadratic term, MSS `otter_legacy.m` 212-223, 250-257 at `6e78e98`); `linear_coupled`; `linear_coupled_modulus` | MSS calibrations < a coupled law identified from logs | circle and zigzag identification |
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
| REMUS 100 | `SpheroidAuv` with `SubmergedRestoring`, `AuvHullLoads`, `FinPairsDeflectionOnly`, `Propeller` | `remus100` | Prestero (2001) and MSS `remus100.m` @ `ed2ce7b`, the line of every value in the comment beside it in `params/parameters.py` |
| Otter | `Catamaran` with `SurfaceRestoring`, `SurfaceHullLoads` | `otter` | MSS `otter_legacy.m` @ `6e78e98`, the same |

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

Sources behind the two tables: the servo settings follow MSS `remus100.m` 110, 114-115 (Fossen, MSS, MIT, `ed2ce7b`; angle limit on the command), Murray-Smith, D. J. (2016), *Inverse simulation methods applied to investigations of actuator nonlinearities in ship steering*, Simulation Notes Europe 26(4), 245-256, pp. 246-247 (steering machine: time constant 3 s, ±35°, ±7 and ±10 °/s) and Sarhadi, P. (2026), *Simple yet effective anti-windup techniques for amplitude and rate saturation: an AUV case study*, arXiv:2601.01302v2, Fig. 4, p. 4 (time constant 0.1 s, ±20°, ±30 °/s; the amplitude limit on the output); the fin forms follow Prestero, T. (2001), *Verification of a six-degree of freedom simulation model for the REMUS autonomous underwater vehicle*, MIT/WHOI MSc thesis, eqs. 4.37 (linear section), 4.40 (rigid-point inflow), 4.41-4.43 (small flow angle), pp. 31-33, and MSS `remus100.m` 234-245 (translational inflow, quadratic-drag section). Each plugin's module docstring carries the full citation and the line of every value.

## Choosing propulsor options

The `Propulsor` plugin (`plugins/force_producers/propulsor.py`, parts in `models/propulsor/`) is one skeleton with six sockets, each holding one part chosen by a selector; the numbers of the chosen parts are parameters, all unset until a named part gives them (`DEFAULTS.md` lists the reference parts). Several propulsors on one vehicle are several entries of its `actuators` list. The wrench is the thrust along the mounting axis at the thrust point, `[f; r x f]`, plus the shaft torque as a couple scaled by `shaft_torque_scale` (physical: 1), plus the housing's wrench. The drive and the law must pass the same unit (`rev/s` for a shaft speed, `1` for a normalised effort): a pair that does not is refused. The command is `[drive command; steering command]`, the state `[drive state; steering state]`.

| Socket | Selector value | Fidelity | Source | Assumptions | Raised to the next level by |
|---|---|---|---|---|---|
| drive | `speed_command_rpm`, `speed_command_rad_s` | lowest: the shaft follows its command at once, limited | published simulator models | no motor or battery dynamics | a logged step of shaft speed (first-order lag), then the motor envelope |
| drive | `power_limited_speed_command_rpm` | construction of an outboard block, kept as a labelled comparison | outboard block | a power cap `eta P_max` and a thrust cap, no source for either | the vendor's motor data or a logged power at full throttle |
| drive | `speed_command_normalized` | as the plain speed drive, command in `[-1, 1]` | the allocator's signed normalisation | limits ahead and astern are the unit's own | the same logged step |
| drive | `normalized_lag` | comparison: an effort in `[-1, 1]` with separate rise and fall times | a shared jet-nozzle plugin | the thrust follows the effort linearly | logged thrust steps (rise and fall times) |
| steering | `fixed`, `servo_ideal`, `servo_no_limits`, `servo_lag_rate_angle` (default), `servo_output_saturated` | the four settings of the servo part (see Servo above) | see Servo | see Servo | a logged step of the steering angle |
| steering | `jet_nozzle_clamp` | comparison: his nozzle angle with an error clamp | a shared jet-nozzle plugin | equals the lag with rate and angle limits only for symmetric limits, 1 s, command in radians | as the servo |
| mounting | `fixed_euler`, `fixed_direction` | fixed shaft from Euler angles or from a direction | published | the shaft does not move | the measured mounting |
| mounting | `azimuth` (+ `steering_convention`) | a unit turning about the body z axis | Fossen 2011, Table 12.3 | one thrust point that does not move with the angle; the convention names the boat's turn for a unit aft of the CO (a unit at or ahead of it is refused) | a swivel form with the point moving with the angle |
| mounting | `steered_rotation` (default), `jet_nozzle_direction` | the first keeps the thrust length at every trim; the second is the shared plugin's formula, which lengthens it up to 0.7 % at 0.2 rad of trim | geometry; a shared jet-nozzle plugin | a rigid nozzle on a vertical axis | the measured nozzle geometry |
| inflow | `none`, `speed_magnitude`, `surge_component` | none; the speed through the water; its surge component | published vehicle models; an outboard block | wake and thrust deduction are constants | a towed or self-propelled test (wake fraction, thrust deduction) |
| law | `open_water_bollard`, `open_water_linearized`, `open_water_linearized_given`, `open_water_polynomial` | the three open-water models of the propeller block (bollard, linear in advance, full polynomial) | Fossen 2011 eq. 12.265; the Wageningen polynomial | no cavitation, no ventilation, quasi-steady | an open-water or bollard test of the real propeller |
| law | `open_water_signed_advance` | the polynomial at a signed advance number with reverse factors | an outboard block | reverse factors are constructions | a reverse bollard test |
| law | `quadratic_asymmetric` | `T = k n abs(n)` with separate ahead and astern coefficients | Fossen 2011 eq. 12.265 | no advance-speed effect | a bollard test ahead and astern |
| law | `normalized_linear` | comparison: thrust proportional to an effort, no flow effect | a shared jet-nozzle plugin | two numbers of which only the product acts | a bollard thrust at full effort |
| law | `jet_bollard_quadratic`, `jet_bilinear_speed`, `jet_momentum_flux_given`, `jet_momentum_flux_bollard_point`, `jet_pump_system_match` | waterjet laws of rising fidelity: quadratic bollard, bilinear in shaft and boat speed, momentum flux with a flow coefficient (given or from one bollard point), pump line matched to the system head | Tajima et al. 2015; Han et al. 2017; Bulten 2006; Loghis and Xiros 2022 | thrust is zero for a stopped or reversed shaft (no reversing bucket); the momentum-flux form equals the bilinear one exactly when the flow is proportional to the shaft speed | a bollard test (flow coefficient), then the pump's head curve from the maker |
| housing | `none`, `leg_normal_force` | none; the normal force of an outboard leg, proportional to the normal flow and the flow speed | the model file `grethe.m` (a construction: a rudder lift slope adapted to the leg) | the leg area and aspect ratio are estimates | a measured side force on the leg in a tow |

A jet is the composition of a drive, the steering and mounting of the nozzle, a jet law and no housing; no named part of a jet carries numbers yet (the published pump point gives a nozzle diameter, a shaft speed and a torque, not the head curve, the mounting or the limits).
