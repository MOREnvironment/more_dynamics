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
  wiring.py (models) · payload_io.py (plugins)   helpers shared by both families
plugin_types/
```

A vehicle is built through it in five steps: (1) pick one plugin per hull part from `plugins/vehicles/hull_parts/<part>/` (a site, a current model, a hull form, a rigid body, an added mass and its Coriolis term, a restoring model, a list of hull loads); (2) pick the force producers from `plugins/force_producers/` (fins, a propeller, a set of them); (3) write them as a tree of part folders with their `parameters.py` under an `.rppws` workspace, the producers as children of the vehicle's `force_producers` slot; (4) build it with rpp's builder; (5) the vehicle's `graph()` is one CasADi graph of all of them, with the commands of the producers as its inputs.

Luka's own files keep their places and are not part of the nested tree: `models/hull_vessel/`, `models/hydrostatics/linear_surface.py`, `models/hydrodynamics/linear_surface.py`, `plugins/vehicle_models/hull_vessel.py`, `plugins/hydrostatics/`, `plugins/hydrodynamics/`, `plugins/force_producers/{jet_nozzle,thruster}.py` and the root `.rppws/`.

### Swapping a force producer: same vehicle, two sets

The REMUS 100 hull is one tree; only its `force_producers` child differs between `remus100_fin_pairs` (`FinPairsDeflectionOnly`, the rudder and stern-plane pair in the deflection-only form) and `remus100_fin_set` (a `ForceProducerSet` of two `LiftingFin`, each a servo with inflow, flow-angle, interference and section parts). Run from this folder with `more_dynamics` registered:

```python
from pathlib import Path
import numpy as np
from rpp_py.data_manager import DataManager
from rpp_plugin_registrator.library_manager import LibraryManager
from rpp_py.context_builder import ComponentContextBuilder
from more_common.casadi_graph import RppCasadiGraph

script = Path("tests/data/vehicles/.rppws/script_descriptions/vehicles.json")
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

A vehicle is a tree of part folders with their `parameters.py` values, built by rpp's builder: the vehicles to use are in `scripts/vehicles/.rppws` (REMUS 100 and Otter, one script description each; `scripts/vehicles/simulate_vehicle.py remus100` runs one), the trees the tests gate against in `tests/data/vehicles/.rppws`. Both are written by `tests/data/vehicles/make_trees.py` from the frozen parameter files, and each part folder carries a `SOURCE.md` with the fidelity level of the part and the line of every value. The trees need `more_dynamics` registered in rpp's registry (`rpp library register ./more_dynamics --link`).

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

A composition's choice of parts, and the measurement that would raise each one, are recorded per unit in the corresponding ledger under `agents-more/30_checks/`; this guide carries the content for a collaborator who has only this library.
