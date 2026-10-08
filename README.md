# more_dynamics
Repository containing more dynamical models

## Plugin types

`more_dynamics/plugin_types/vehicle_parts.capnp` defines the typed slots of a `MarineCraft6DOF` vehicle: site, current, hull form, rigid body, added mass, added-mass Coriolis term and the section drag law of a cross-flow strip, plus the servo, inflow, flow-angle, interference and section parts of a fin. Every part is a plugin class in the form of `HullVessel`'s children: its `graph()` returns the `CasadyPayload` of `vehicle_model.capnp`, and each input and output is named in the payload's descriptions. Restoring parts are `HydrostaticsModel` plugins, hull loads `HydrodynamicsModel` plugins and thrusters and fins `ForceProducer` plugins, so the plugins of `vehicle_model.capnp` plug into the same vehicle. A matrix crosses the payload as its column-major `vec` (a 6x6 matrix is 36 values).

## Layout: folders nested by level

A folder holds things of one level only; a part of a thing lives inside that thing's folder. `models/` (the CasADi functions), `plugins/` (the plugin classes) and `tests/` use the same tree, so a block, its plugin and its tests are found at the same path in each. Vehicles and force producers are two families of the first level: any force producer plugs into any vehicle's `force_producers` list.

```
more_dynamics/{models,plugins}/            (tests/ mirrors the same tree)
  vehicles/                                level 1: a vehicle
    marine_craft_6dof/                     the equation of motion every vehicle shares
    hull_parts/                            level 2: the vehicle's own parts
      site/  current/  hull_form/
      rigid_body/  added_mass/  added_mass_coriolis/  restoring/
      hydrodynamic_loads/
        damping/  lift_drag/  surge_resistance/
        cross_flow/
          section_drag/                    level 3: inside its parent
  force_producers/                         level 1: plugged into any vehicle
    fin/
      servo/  inflow/  flow_angle/  interference/  section/
    propulsor/
    thrusters/
    shared/                                helpers of the force producers (models)
  shared/                                  helpers of both families: wiring.py (models), payload_io.py (plugins)
plugin_types/
```

Each test folder keeps its frozen reference data beside it in `data/` (the MATLAB generators, their outputs and a `SOURCE.md` naming the source of every file): `tests/vehicles/hull_parts/<part>/data/`, `tests/force_producers/data/`, `tests/force_producers/fin/data/`, and `tests/vehicles/data/` for the vehicle parameter files and the gate trees. References that need MSS skip without `MSS_DIR`; the frozen files run everywhere.

A vehicle is built through it in five steps: (1) pick one plugin per hull part from `plugins/vehicles/hull_parts/<part>/` (a site, a current model, a hull form, a rigid body, an added mass and its Coriolis term, a restoring model, a list of hull loads); (2) pick the force producers from `plugins/force_producers/` (fins, a propeller, a set of them); (3) write them as a tree of part folders with their `parameters.py` under an `.rppws` workspace, the producers as children of the vehicle's `force_producers` slot; (4) build it with rpp's builder; (5) the vehicle's `graph()` is one CasADi graph of all of them, with the commands of the producers as its inputs.

Luka's own files keep their places and are not part of the nested tree: `models/hull_vessel/`, `models/hydrostatics/linear_surface.py`, `models/hydrodynamics/{linear_surface,crossflow_surface}.py`, `plugins/vehicle_models/hull_vessel.py`, `plugins/hydrostatics/`, `plugins/hydrodynamics/`, `plugins/force_producers/{jet_nozzle,thruster}.py`, their tests `tests/{vehicle_models,hydrostatics,hydrodynamics}/`, `scripts/{thruster,vehicle_models}/` and his parts and script descriptions in the root `.rppws/`.

### Swapping a force producer: same vehicle, two sets

The REMUS 100 hull is one tree; only its `force_producers` child differs between `remus100_fin_pairs` (`FinPairsDeflectionOnly`, the rudder and stern-plane pair in the deflection-only form) and `remus100_fin_set` (a `ForceProducerSet` of two `LiftingFin`, each a servo with inflow, flow-angle, interference and section parts). Run from this folder with `more_dynamics` registered:

```python
from pathlib import Path
import numpy as np
from rpp_py.data_manager import DataManager
from rpp_plugin_registrator.library_manager import LibraryManager
from rpp_py.context_builder import ComponentContextBuilder
from more_common.casadi_graph import RppCasadiGraph

script = Path("tests/vehicles/data/.rppws/script_descriptions/vehicles.json")
builder = ComponentContextBuilder(data_manager=DataManager(library_manager=LibraryManager()))

def vehicle(configuration):
    context = builder.build_script_from_description_path(str(script), configuration=configuration)
    context.initialize()
    return context.get_component("vessels")[0]

def wrench(vessel, nu, command):
    graph = RppCasadiGraph(vessel.graph())
    x = np.concatenate([np.zeros(6), nu])        # eta = 0, nu = nu_r (still water)
    u = np.concatenate([np.zeros(3), command])   # no current, then the fin commands
    descriptions = vessel.graph().outputDescription
    out = np.split(np.asarray(graph.output(x, u)).ravel(), np.cumsum([d.size for d in descriptions])[:-1])
    return dict(zip([d.name for d in descriptions], out))["force_producers.0.generated_force"]

nu, command = np.array([1.5, 0.05, 0.02, 0.0, 0.01, 0.03]), np.array([0.1, -0.05])
pairs = wrench(vehicle("remus100_fin_pairs"), nu, command)
fin_set = wrench(vehicle("remus100_fin_set"), nu, command)
print(np.max(np.abs(pairs - fin_set)))           # ~0: the two sets agree on the deflection-only setting
```

The two agree because both are the deflection-only setting; they part the moment one fin part is swapped (a `Servo` with its first-order lag instead of the ideal one, say): change that child's `parameters.py` or plugin in the tree and build again, no code changes.

## Vehicles

`MarineCraft6DOF` takes one part per slot (water, current, hull form, rigid body, added mass, its Coriolis term, restoring) and lists of hull loads and force producers, and connects every part input by name to a vehicle quantity (`pose`, `velocity`, `relative_velocity`, `mass_matrix`) or to an output of an earlier single-slot part with the same name and size; a missing or mis-sized input is refused with the slot and the part named. A force producer's input that nothing feeds is a command, listed by `getInputDescriptions()`. `open_inputs` leaves a part parameter open as a vehicle input (identification, a time-varying current), and `diagnostic_outputs` appends named quantities (`mass_matrix`, `hydrodynamic_loads.0.damping_matrix`, ...) after the twelve states.

A vehicle is a tree of part folders with their `parameters.py` values, built by rpp's builder. The vehicles to use are named parts of the library's own `.rppws/parts/`, beside Luka's hull vessel, and one script description, `.rppws/script_descriptions/vehicle_simulation.json`, holds one configuration per vehicle for `scripts/vehicles/simulate_vehicle.py`:

    python scripts/vehicles/simulate_vehicle.py --configuration remus100 --surge-speed 1.5 --duration 20

| Name | Type | Configuration | Part id | Sources |
|---|---|---|---|---|
| REMUS 100 | `MarineCraft6DOF` (torpedo class) | `remus100` | `more_dynamics__marine_craft6_d_o_f/24d92e9a-fa68-528a-8278-5cdc22cbb17b` | Prestero (2001) and MSS `remus100.m` @ `cc07579`, per value in each part's `SOURCE.md` |
| Otter | `MarineCraft6DOF` (catamaran class) | `otter` | `more_dynamics__marine_craft6_d_o_f/9e3c87b5-9998-5d85-86d8-7bfdf6a80222` | MSS `otter.m` @ `cc07579`, per value in each part's `SOURCE.md` |

The trees the tests gate against are in `tests/vehicles/data/.rppws`. Both sets are written by `tests/vehicles/data/make_trees.py` from the frozen parameter files (it rewrites only its own part folder and script description in the root `.rppws`), and each part folder carries a `SOURCE.md` with the fidelity level of the part and the line of every value. The trees need `more_dynamics` registered in rpp's registry (`rpp library register ./more_dynamics --link`).

## Choosing fin and servo options

A `LiftingFin` composes a servo part with inflow, flow-angle, interference and section parts. Pick each one with document values first; the measurement or identification route raises it to the next fidelity level once that measurement exists.

### Servo (`force_producers/fin/servo/servo.py`)

One device, three switches, not one file per form: `dynamics` (`"none"` / `"first_order_lag"`), `rate_limit` (`False` / `True`, needs the lag's own state), `angle_limit` (`"none"` / `"on_command"` / `"on_output"`). An invalid combination (a rate limit with no lag; an output-side limit with no lag) is refused, naming the switch.

| Setting (`dynamics` / `rate_limit` / `angle_limit`) | Fidelity | Source | Assumptions | Raised to the next level by |
|---|---|---|---|---|
| `none` / `False` / `on_command` | lowest — instantaneous angle limit, no dynamics | published vehicle model | the servo is fast enough to ignore | a step response showing a measurable rate or lag |
| `none` / `False` / `none` | no limits at all | — | a placeholder for an unlimited actuator, not a vehicle's servo | any of the limited settings, once a limit is known |
| `first_order_lag` / `True` / `on_command` | default — amplitude limit on the command, first-order lag, rate limit | published steering-machine model | the lag's time constant is either given directly or derived from a proportional band; the command enters unscaled | a quay-day step test giving the time constant, rate limit and angle limit directly |
| `first_order_lag` / `True` / `on_output` | comparison setting beside the default — amplitude limit on the integrator's own state instead of on the command | published anti-windup case study | equals the default only while the command stays inside the angle limit and the state starts unwound; diverges after a saturating reversal | not a higher-fidelity route on its own — a comparison setting kept beside the default |

Marie's fins and Grethe's steering: `first_order_lag` / `True` / `on_command`, numbers from a quay-day step test and her logs.

### Fin parts (`force_producers/fin/`)

| Socket | Part | Fidelity | Source | Assumptions | Raised to the next level by |
|---|---|---|---|---|---|
| inflow | `rigid_point` | lowest | published vehicle model | the fin sees the vehicle's velocity at one point, uncorrected for its own rotation | the full `translational` form once a rotation rate matters |
| inflow | `translational` | default | published vehicle model | the fin's own lever arm correction is included | — |
| flow angle | `small_angle` | default | published vehicle model | linear ratio, valid away from hover | an `exact` large-angle form once hover-speed manoeuvres are simulated |
| interference | `none` | default (no interference modelled) | — | one fin's flow does not affect another's | the slender-body factor once its source is read |
| section | `quadratic_drag` | MSS shortcut — flagged, not the default when a physics form exists | published simulator model | no zero-lift drag, no induced drag; registered as a shortcut | `linear_section` with an identified zero-lift drag |
| section | `linear_section` | default | published vehicle model | small-angle lift and a fixed zero-lift drag | a lifting-line induced-drag form once its source is read |

Sources behind the two tables: the servo settings follow MSS `remus100.m` 109, 113-114 (Fossen, MSS, MIT, `cc07579`; angle limit on the command), Murray-Smith, D. J. (2016), *Inverse simulation methods applied to investigations of actuator nonlinearities in ship steering*, Simulation Notes Europe 26(4), 245-256, pp. 246-247 (steering machine: time constant 3 s, ±35°, ±7 and ±10 °/s) and Sarhadi, P. (2026), *Simple yet effective anti-windup techniques for amplitude and rate saturation: an AUV case study*, arXiv:2601.01302v2, Fig. 4, p. 4 (time constant 0.1 s, ±20°, ±30 °/s; the amplitude limit on the output); the fin parts follow Prestero, T. (2001), *Verification of a six-degree of freedom simulation model for the REMUS autonomous underwater vehicle*, MIT/WHOI MSc thesis, eqs. 4.37 (linear section), 4.40 (rigid-point inflow), 4.41-4.43 (small flow angle), pp. 31-33, and MSS `remus100.m` 234-245 (translational inflow, quadratic-drag section). Each part's module docstring and the `SOURCE.md` beside each part's parameters carry the full citation and the line of every value.
