# NPS AUV II

A named part of `more_dynamics::DerivativeAuv`: the second AUV type of the library, whose full inertia tensor,
added mass and hull forces are the published nondimensional derivatives of the Naval Postgraduate School AUV II
(Healey and Lienard 1993) as MSS carries them in `CRAFT/AUV/models/npsauv.m` (T. I. Fossen, MIT licence, revision
`cc07579b63630342802a198c34c74ac7e096d5fb`). Every value is read from MSS's instrumented workspace
(`tests/vehicle_models/data/npsauv/nps_auv_ii_parameters.json`), not typed; the paper itself is not read
(`tests/vehicle_models/data/npsauv/SOURCE.md`, needs-access list).

## What this part composes

| Slot | Component | Plugin |
|---|---|---|
| (root) | NPS AUV II | `DerivativeAuv`: full inertia tensor, the 14 acceleration derivatives (added mass) |
| `site` | NPS AUV II site | `Site`, given form: gravity 9.81, density 1025 (`npsauv.m` 105, 109) |
| `hydrostatics` | NPS AUV II hydrostatics | `SubmergedRestoring`, neutral buoyancy (`npsauv.m` 108, 214) |
| `hydrodynamics` | NPS AUV II hull loads | `CoefficientHullLoads`: the 52 hull-only derivatives, the cross-flow strip integral (`npsauv.m` 217-235) and the rigid-body Coriolis-centripetal terms `npsauv.m` writes into the hull force |
| `actuators` | NPS AUV II prescribed wrench | `PrescribedWrench`: an external wrench input, no fins, no propeller |

## What this part does not compose, and why

**No fins and no propeller.** The control forces of `npsauv.m` (`tau_control`, lines 176-212) are an empirical
polynomial in each fin deflection and the vehicle's own relative velocity, scaled by a through-water correction
shared across the stern plane, the bow planes and the roll moment and by the propeller's own sign, and the model
lags the shaft speed as a fifth actuator state. The library's `Fin` is a lifting-surface force and its `Propeller`
carries no lag state; neither reproduces that polynomial, and the slot has no place for a coupling between two
actuators. The 27 actuator-coupled derivatives (`Xqds`, `Ydr`, `Zqn`, ...) are therefore declared by no plugin and
are not in this part. Until an actuator design for them is decided, this part is the **hull-only** vehicle: the
configuration whose dynamics equal `npsauv_hull.m`'s on all 540 frozen cases to 4.4e-15 (maximum difference, measured), and which is driven
through the wrench input.

**No viscosity of its own.** `npsauv.m` has no viscosity. The `Site` plugin declares one, so the part carries
`kinematic_viscosity = 1e-06` m^2/s, the value the REMUS 100 and Otter parts use (`cylinderDrag.m` 78-80). It is
not an NPS value, it is an estimate kept only to satisfy the plugin, and it is read by nothing: changing it to
1e-03 changes no output of the vehicle's dynamics (checked, difference 0.0). Awaiting the owner's approval or a
replacement.

## Verified

* The part builds through rpp's `ComponentContextBuilder` and runs in `scripts/vehicles/vehicle_simulation.py`
  (configuration `npsauv`).
* Its dynamics at 50 random states and wrenches equal those of the hull-only vehicle built by the frozen-data test
  (`tests/vehicle_models/test_derivative_auv.py`, `_build_hull`) with maximum difference 0.0.

## Values

Every parameter carries its source in the comment beside it in `params/parameters.py` of its component; `kind`
is as in the frozen parameter file (`via MSS` = read from MSS, the paper not read).

### NPS AUV II (vehicle)

| Parameter | Value | Source |
|---|---|---|
| `integration_max_step` | 0.05 | numerical setting (plugin default) |
| `length` | 5.3 | npsauv.m:105 L = 5.3 (via MSS, Healey & Lienard 1993 not read) |
| `body_mass` | 5443.425076452599 | npsauv.m:109 mass = W / g (derived from the given weight) (derived, via MSS) |
| `inertia_diagonal` | [2038, 13587, 13587] | npsauv.m:110 [Ix Iy Iz] (via MSS, Healey & Lienard 1993 not read) |
| `inertia_products` | [-13.58, -13.58, -13.58] | npsauv.m:111 [Ixy Iyz Ixz] (full tensor, not diagonal) (via MSS, Healey & Lienard 1993 not read) |
| `body_center_of_gravity` | [0, 0, 0.061] | npsauv.m:106 [xG yG zG] (via MSS, Healey & Lienard 1993 not read) |
| `weight` | 53400 | npsauv.m:108 W = 53400 (via MSS, Healey & Lienard 1993 not read) |
| `buoyancy` | 53400 | npsauv.m:108 B = 53400 (neutral, W = B as given) (via MSS, Healey & Lienard 1993 not read) |
| `center_of_buoyancy` | [0, 0, 0] | npsauv.m:107 [xB yB zB] (via MSS, Healey & Lienard 1993 not read) |
| `derivative_scaling` | 'nondimensional' | npsauv.m 116-124: the published nondimensional form (Fossen 2011, App. D.2 prime-scaling) |
| `current_form` | 'yaw_rate_terms' | npsauv.m 83-90: the MSS form of the current |
| `nondim_Xudot` | -0.0076 | npsauv.m: nondimensional hydrodynamic derivative Xudot (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Yvdot` | -0.055 | npsauv.m: nondimensional hydrodynamic derivative Yvdot (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Ypdot` | 0.00012 | npsauv.m: nondimensional hydrodynamic derivative Ypdot (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Yrdot` | 0.0012 | npsauv.m: nondimensional hydrodynamic derivative Yrdot (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Zwdot` | -0.24 | npsauv.m: nondimensional hydrodynamic derivative Zwdot (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Zqdot` | -0.0068 | npsauv.m: nondimensional hydrodynamic derivative Zqdot (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Kvdot` | 0.00012 | npsauv.m: nondimensional hydrodynamic derivative Kvdot (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Kpdot` | -0.001 | npsauv.m: nondimensional hydrodynamic derivative Kpdot (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Krdot` | -3.4e-05 | npsauv.m: nondimensional hydrodynamic derivative Krdot (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Mwdot` | -0.0068 | npsauv.m: nondimensional hydrodynamic derivative Mwdot (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Mqdot` | -0.017 | npsauv.m: nondimensional hydrodynamic derivative Mqdot (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Nvdot` | 0.0012 | npsauv.m: nondimensional hydrodynamic derivative Nvdot (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Npdot` | -3.4e-05 | npsauv.m: nondimensional hydrodynamic derivative Npdot (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Nrdot` | -0.0034 | npsauv.m: nondimensional hydrodynamic derivative Nrdot (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |

### NPS AUV II site

| Parameter | Value | Source |
|---|---|---|
| `site_form` | 'given' | the given form: gravity and density as npsauv.m states them |
| `latitude` | [] | unset: read by no selected form of this part |
| `gravity` | 9.81 | npsauv.m:105 g = 9.81 (via MSS, Healey & Lienard 1993 not read) |
| `water_density` | 1025 | npsauv.m:109 rho = 1025 (via MSS, Healey & Lienard 1993 not read) |
| `kinematic_viscosity` | 1e-06 | NOT AN NPS VALUE and read by no NPS form: npsauv.m has no viscosity; the Site plugin declares it, so a build needs one; 1e-06 m^2/s is the value cylinderDrag.m 78-80 and the REMUS 100 / Otter parts use; owner to approve or replace |

### NPS AUV II hydrostatics

| Parameter | Value | Source |
|---|---|---|
| `buoyancy_method` | 'neutral' | neutral: B = W as npsauv.m 214 gives it |
| `center_of_buoyancy` | [0, 0, 0] | npsauv.m:107 [xB yB zB] (via MSS, Healey & Lienard 1993 not read) |
| `displaced_volume` | [] | unset: read by no selected form of this part |
| `buoyancy` | 53400 | npsauv.m:108 B = 53400 (neutral, W = B as given) (via MSS, Healey & Lienard 1993 not read) |

### NPS AUV II hull loads

| Parameter | Value | Source |
|---|---|---|
| `nondim_Xpp` | 0.007 | npsauv.m: nondimensional hydrodynamic derivative Xpp (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Xqq` | -0.015 | npsauv.m: nondimensional hydrodynamic derivative Xqq (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Xrr` | 0.004 | npsauv.m: nondimensional hydrodynamic derivative Xrr (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Xpr` | 0.00075 | npsauv.m: nondimensional hydrodynamic derivative Xpr (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Xwq` | -0.2 | npsauv.m: nondimensional hydrodynamic derivative Xwq (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Xvp` | -0.003 | npsauv.m: nondimensional hydrodynamic derivative Xvp (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Xvr` | 0.02 | npsauv.m: nondimensional hydrodynamic derivative Xvr (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Xvv` | 0.053 | npsauv.m: nondimensional hydrodynamic derivative Xvv (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Xww` | 0.17 | npsauv.m: nondimensional hydrodynamic derivative Xww (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Ypq` | 0.004 | npsauv.m: nondimensional hydrodynamic derivative Ypq (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Yqr` | -0.0065 | npsauv.m: nondimensional hydrodynamic derivative Yqr (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Yp` | 0.003 | npsauv.m: nondimensional hydrodynamic derivative Yp (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Yr` | 0.03 | npsauv.m: nondimensional hydrodynamic derivative Yr (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Yvq` | 0.024 | npsauv.m: nondimensional hydrodynamic derivative Yvq (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Ywp` | 0.23 | npsauv.m: nondimensional hydrodynamic derivative Ywp (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Ywr` | -0.019 | npsauv.m: nondimensional hydrodynamic derivative Ywr (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Yv` | -0.1 | npsauv.m: nondimensional hydrodynamic derivative Yv (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Yvw` | 0.068 | npsauv.m: nondimensional hydrodynamic derivative Yvw (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Zpp` | 0.00013 | npsauv.m: nondimensional hydrodynamic derivative Zpp (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Zpr` | 0.0067 | npsauv.m: nondimensional hydrodynamic derivative Zpr (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Zrr` | -0.0074 | npsauv.m: nondimensional hydrodynamic derivative Zrr (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Zq` | -0.14 | npsauv.m: nondimensional hydrodynamic derivative Zq (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Zvp` | -0.048 | npsauv.m: nondimensional hydrodynamic derivative Zvp (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Zvr` | 0.045 | npsauv.m: nondimensional hydrodynamic derivative Zvr (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Zw` | -0.3 | npsauv.m: nondimensional hydrodynamic derivative Zw (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Zvv` | -0.068 | npsauv.m: nondimensional hydrodynamic derivative Zvv (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Kpq` | -6.9e-05 | npsauv.m: nondimensional hydrodynamic derivative Kpq (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Kqr` | 0.017 | npsauv.m: nondimensional hydrodynamic derivative Kqr (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Kp` | -0.011 | npsauv.m: nondimensional hydrodynamic derivative Kp (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Kr` | -0.00084 | npsauv.m: nondimensional hydrodynamic derivative Kr (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Kvq` | -0.0051 | npsauv.m: nondimensional hydrodynamic derivative Kvq (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Kwp` | -0.00013 | npsauv.m: nondimensional hydrodynamic derivative Kwp (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Kwr` | 0.014 | npsauv.m: nondimensional hydrodynamic derivative Kwr (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Kv` | 0.0031 | npsauv.m: nondimensional hydrodynamic derivative Kv (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Kvw` | -0.19 | npsauv.m: nondimensional hydrodynamic derivative Kvw (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Mpp` | 5.3e-05 | npsauv.m: nondimensional hydrodynamic derivative Mpp (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Mpr` | 0.005 | npsauv.m: nondimensional hydrodynamic derivative Mpr (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Mrr` | 0.0029 | npsauv.m: nondimensional hydrodynamic derivative Mrr (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Muq` | -0.068 | npsauv.m: nondimensional hydrodynamic derivative Muq (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Mvp` | 0.0012 | npsauv.m: nondimensional hydrodynamic derivative Mvp (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Mvr` | 0.017 | npsauv.m: nondimensional hydrodynamic derivative Mvr (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Muw` | 0.1 | npsauv.m: nondimensional hydrodynamic derivative Muw (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Mvv` | -0.026 | npsauv.m: nondimensional hydrodynamic derivative Mvv (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Npq` | -0.021 | npsauv.m: nondimensional hydrodynamic derivative Npq (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Nqr` | 0.0027 | npsauv.m: nondimensional hydrodynamic derivative Nqr (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Np` | -0.00084 | npsauv.m: nondimensional hydrodynamic derivative Np (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Nr` | -0.016 | npsauv.m: nondimensional hydrodynamic derivative Nr (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Nvq` | -0.01 | npsauv.m: nondimensional hydrodynamic derivative Nvq (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Nwp` | -0.017 | npsauv.m: nondimensional hydrodynamic derivative Nwp (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Nwr` | 0.0074 | npsauv.m: nondimensional hydrodynamic derivative Nwr (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Nv` | -0.0074 | npsauv.m: nondimensional hydrodynamic derivative Nv (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `nondim_Nvw` | -0.027 | npsauv.m: nondimensional hydrodynamic derivative Nvw (lines 126-156) (via MSS, Healey & Lienard 1993 not read) |
| `cross_flow_drag_coefficients` | [0.5, 0.6] | npsauv.m:112 [Cdy Cdz] (via MSS, Healey & Lienard 1993 not read) |
| `cross_flow_section_dimensions` | [0.53, 0.53] | npsauv.m:220-221 Hx, Bx (average height/width of the submerged body) (via MSS, Healey & Lienard 1993 not read) |
| `cross_flow_sections` | 10 | npsauv.m:218 dxL = L/10 (number of strips; a numerical grid density, not a physical quantity) (numerical device) |
