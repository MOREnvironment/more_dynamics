"""Gates of the higher-fidelity options of ``SpheroidAuv``, ``AuvHullLoads``
and ``Fin``, and of the two REMUS 100 named parts that select them: "REMUS 100
(published reference)" (the reference-model forms: shortcut Coriolis, wing lift
and drag, linear roll damping) and "REMUS 100 (estimated)" (every Coriolis term,
point-force body lift, geometric fin slope, an estimated roll reaction).

The reference is the numpy second implementation in ``data/remus100_fidelity/``
(``reference_torpedo_auv.py``, parameter sets in ``fidelity_parts.py``) and its
frozen turn and stability table (``remus100_fidelity_turn_and_stability.csv``),
which its own generator writes. Published numbers are typed with their page
(Prestero 2001; Fossen 2011). Comparisons with runs of other simulators are not
made here.

Contract (parameter names):

* ``AuvHullLoads``: ``body_lift_form`` (``vertical_plane_wing`` = the lift half
  of the wing block | ``point_force_both_planes``) with
  ``body_lift_slope_method`` (``slender_body_hoerner`` reading
  ``body_lift_slope_per_degree``, refused outside 6.7 <= l/d <= 10 | ``given``
  reading ``body_lift_slope``) and ``pressure_centre_method``
  (``fraction_of_length_from_nose`` reading ``pressure_centre_fraction``, nose at
  +L/2 | ``given`` reading ``pressure_centre``); ``hull_drag_form``
  (``wing_drag`` | ``quadratic_axial`` reading ``hull_frontal_drag_coefficient``);
  ``roll_damping_form`` (``linear`` | ``linear_plus_quadratic`` reading
  ``roll_quadratic_damping``); outputs ``body_lift_force`` and ``hull_drag_force``
  (``lift_drag_force`` stays their sum).
* ``Fin``: ``lift_law`` (``lift_slope`` reading ``lift_slope`` | ``finite_wing``)
  with ``lift_slope_form`` ``low_aspect_ratio`` reading ``fin_span``, ``fin_area``,
  ``aspect_ratio_factor`` and ``section_efficiency``; the sign of a deflection
  is the ``lift_axis`` (no second label); ``FinPairsDeflectionOnly`` declares the
  one ``convention`` its model has.
* ``SpheroidAuv``: ``rigid_body_form`` (``homogeneous_spheroid`` |
  ``given_inertia`` reading ``inertia_diagonal`` and ``inertia_products`` about
  the CO, the names and convention of ``DerivativeAuv``).
* Workspace: configurations ``remus100_published_reference`` and
  ``remus100_estimated`` of ``vehicle_simulation.json``, the named parts with the
  values and selectors of ``remus100_estimated_parameters.json``.

Not contracted here (not derived or not settled): the hull-profile, strip-theory
added-mass and profile cross-flow damping forms, and a part built on Prestero's
own geometry and coefficients.

Tolerances. Derivatives and forces: 1e-12 relative (the same equations in
another summation order). Turn rates and roots: the plugins' CasADi graphs
against the numpy reference over 1200 RK4 steps, measured on the configurations
that build without the new options: states within 2.2e-16 after 2 s, the 10 s
mean yaw rate equal to the last bit, so 1e-9 rad/s is a bound with margin; roots
within 2.1e-7 1/s (the reference's central difference, step 1e-7, against
CasADi's exact Jacobian), so 1e-6. Published numbers: half a unit of the last
printed digit.

References
----------
Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics and Motion
    Control. Wiley. Ch. 3, eq. 3.24, p. 49, eq. 3.34, p. 50; ch. 12,
    eqs. 12.17-12.24, p. 351; eq. 12.226, p. 400.
Prestero, T. (2001). Verification of a six-degree of freedom simulation
    model for the REMUS autonomous underwater vehicle. MSc thesis, MIT/WHOI.
    Tables 2.1-2.2, pp. 15-17; eqs. 4.5-4.6, 4.10, 4.30-4.39, pp. 25-31;
    Fig. 7-7, p. 57; Table 8.1, p. 60; Table A.1, p. 102; Table D.1, p. 112.

Author:    Enio Krizman
"""

import csv
import json
import sys
from pathlib import Path

import casadi as ca
import numpy as np
import pytest

import vehicle_contract as vc

FIDELITY = vc.DATA / "remus100_fidelity"
sys.path.insert(0, str(FIDELITY))
import fidelity_parts as fp  # noqa: E402
import reference_torpedo_auv as ref  # noqa: E402

LIBRARY = Path(__file__).resolve().parents[2]
RPPWS = LIBRARY / ".rppws"
DESCRIPTION = RPPWS / "script_descriptions" / "vehicle_simulation.json"
PUBLISHED_PART, PUBLISHED_CONFIGURATION = "REMUS 100 (published reference)", "remus100_published_reference"
ESTIMATED_PART = fp.part_file()["part"]
ESTIMATED_CONFIGURATION = fp.part_file()["configuration"]
REFERENCE_FILE = FIDELITY / "remus100_fidelity_turn_and_stability.csv"
FORCE_TOLERANCE = 1e-12          # relative, same equations
TURN_TOLERANCE = 1e-9            # rad/s on the 10 s mean yaw rate (measured: equal to the last bit)
ROOT_TOLERANCE = 1e-6            # 1/s (measured 2.1e-7: CasADi Jacobian against central differences, eps 1e-7)
STEP = 0.05                      # s, RK4
DURATION, AVERAGE = 60.0, 10.0
SPEED = 1.54                     # Prestero 2001, Table 8.1, p. 60
RUDDER = np.deg2rad(4.0)         # Prestero 2001, p. 53

HULL_OPTIONS = ("body_lift_form", "body_lift_slope_method", "body_lift_slope_per_degree", "body_lift_slope",
                "pressure_centre_method", "pressure_centre_fraction", "pressure_centre", "hull_drag_form",
                "roll_damping_form", "roll_quadratic_damping")
FIN_OPTIONS = ("lift_law", "lift_slope_form", "fin_span", "aspect_ratio_factor", "section_efficiency")


# ---------------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------------
def _reference(case, fraction=None):
    with open(REFERENCE_FILE, newline="") as handle:
        for row in csv.DictReader(handle):
            if row["case"] == case and row["pressure_centre_fraction"] == ("" if fraction is None else f"{fraction:.2f}"):
                return {k: (v if k == "case" else float(v) if v else None) for k, v in row.items()}
    raise KeyError((case, fraction))


def _require(cls, names, why):
    """Fail (red) naming the options ``cls`` does not declare."""
    missing = [n for n in names if n not in {d.name for d in cls.PARAMETERS}]
    if missing:
        pytest.fail(f"{cls.__name__} does not declare {missing} ({why})")


def _geometric_slope_selectors():
    """The ``Fin`` options that take the lift slope from the fin geometry, from the estimated part's file."""
    record = fp.part_file()["fin"]
    return {"lift_law": record["lift_law"]["value"], "lift_slope_form": record["lift_slope_form"]["value"],
            "fin_span": record["fin_span"]["value"], "aspect_ratio_factor": record["aspect_ratio_factor"]["value"],
            "section_efficiency": record["section_efficiency"]["value"]}


def _fin(fin, selectors=None):
    """A ``Fin`` child from a reference fin set (servo none, angle limit on the command)."""
    p = vc.plugins()
    params = {"fin_position": list(fin["fin_position"]), "chord_axis": list(fin["chord_axis"]),
              "lift_axis": list(fin["lift_axis"]), "fin_area": fin["fin_area"], "max_deflection": fin["max_deflection"],
              "servo_dynamics": "none", "servo_rate_limit": False, "servo_angle_limit": "on_command",
              "inflow": fin["inflow"], "flow_angle": fin["flow_angle"], "interference": "none",
              "section": "lift_times_angle", "lift_slope": fin["lift_slope"],
              "time_constant": 1.0, "max_rate": 1.0, "zero_lift_drag": 0.0}   # read by no selected form
    params.update(selectors or {})
    return (p.Fin, params)


def _fin_pairs_or_fins(part):
    """The reference model's two pairs as ``FinPairsDeflectionOnly`` (two fin entries) or one ``Fin`` per fin."""
    p = vc.plugins()
    if len(part["fins"]) == 2 and part["fins"][0]["flow_angle"] == "none":
        return [(p.FinPairsDeflectionOnly, vc.pick(p.FinPairsDeflectionOnly, part))]
    return [_fin(f) for f in part["fins"]]


def _vehicle(part, hull_options=None, fins=None, vehicle_options=None, actuators=None):
    """A vehicle from a flat reference set: the plugins' names picked from it, the new options handed over."""
    p = vc.plugins()
    hull = {**vc.pick(p.AuvHullLoads, part), **(hull_options or {})}
    vehicle = {**vc.pick(p.SpheroidAuv, part), "coriolis_form": part["coriolis_form"], "current_form": "none",
               **(vehicle_options or {})}
    site = {"site_form": "latitude", "latitude": part["latitude"], "water_density": part["water_density"],
            "kinematic_viscosity": part["kinematic_viscosity"]}
    if actuators is None:
        actuators = (fins if fins is not None else _fin_pairs_or_fins(part)) + [(p.Propeller, vc.pick(p.Propeller, part))]
    return vc.build_context(p.SpheroidAuv, vehicle, site=(p.Site, site),
                            hydrostatics=(p.SubmergedRestoring, vc.pick(p.SubmergedRestoring, part, ("center_of_buoyancy",))),
                            hydrodynamics=(p.AuvHullLoads, hull), actuators=actuators)


def _hull_options(part):
    return {k: part[k] for k in HULL_OPTIONS if k in part}


def _paired_vehicle(part):
    """The part's forms through the plugins."""
    p = vc.plugins()
    _require(p.AuvHullLoads, ("body_lift_form", "body_lift_slope_method", "body_lift_slope_per_degree",
                              "pressure_centre_method", "pressure_centre_fraction", "hull_drag_form", "roll_damping_form"),
             "point-force body lift in both planes, Prestero 2001 eqs. 4.27-4.35")
    if part["fins"][0]["lift_slope"] < 1.0:     # the published lift coefficients: given
        fins = [_fin(f) for f in part["fins"]]
    else:
        _require(p.Fin, FIN_OPTIONS, "fin lift slope from the fin geometry, Prestero 2001 eqs. 4.38-4.39")
        fins = [_fin(f, _geometric_slope_selectors()) for f in part["fins"]]
    return _vehicle(part, _hull_options(part), fins)


def _turn_rate(vessel, command):
    """Mean and spread of r over the last 10 s of a 60 s run from 1.54 m/s straight, input held (RK4, h 0.05)."""
    x = np.zeros(12)
    x[6] = SPEED
    n, keep, rates = int(round(DURATION / STEP)), int(round(AVERAGE / STEP)), []
    for k in range(n):
        x = vessel.step_vector(x, command, STEP)
        if k >= n - keep:
            rates.append(x[11])
    return float(np.mean(rates)), float(np.ptp(rates))


def _for_vehicle(vessel, command):
    """The command of a reference part for ``vessel``: the reference part of the published forms gives one command per
    fin pair (rudder, stern plane) and the shaft speed; a part with one actuator per fin (the rudder pair and the
    stern-plane pair, two fins each) gets each pair's command for both of its fins."""
    command = np.asarray(command, float)
    size = vc.input_size(vessel)
    if len(command) == size:
        return command
    assert len(command) == 3 and size == 5, (len(command), size)
    return np.array([command[0], command[0], command[1], command[1], command[2]])


def _linearisation(vessel, command):
    """Straight-running equilibrium at ``command`` (Newton on the plugin's own dynamics, CasADi Jacobian) and
    the Jacobian over the motion states there."""
    step = vc.graph_of(vessel).step
    nx = step.size1_in(0)
    xs, us = ca.MX.sym("x", nx), ca.MX.sym("u", step.size1_in(1))
    jac = ca.Function("motion_jacobian", [xs, us], [ca.jacobian(step(xs, us), xs)])
    x = np.zeros(nx)
    x[6] = SPEED
    motion = ref.MOTION
    for _ in range(100):
        residual = np.asarray(step(x, command)).ravel()[motion]
        if np.abs(residual).max() < 1e-12:
            break
        j = np.asarray(jac(x, command))[np.ix_(motion, motion)]
        x[motion] -= np.linalg.solve(j, residual)
    return x, np.asarray(jac(x, command))[np.ix_(motion, motion)]


def _roots(vessel, command):
    """Sway-yaw roots, largest real part first (Fossen 2011, eqs. 12.17-12.21, p. 351)."""
    x, a = _linearisation(vessel, command)
    return x, sorted(ref.sway_yaw_roots(a), key=lambda z: (-z.real, z.imag)), np.linalg.eigvals(a).real.max()


def _assert_matches(vessel, part, case, fraction=None):
    expected = _reference(case, fraction)
    x, roots, max_re = _roots(vessel, _for_vehicle(vessel, fp.rudder_command(part, 0.0)))
    assert abs(x[6] - expected["cruise_u"]) <= 1e-9
    assert abs(roots[0].real - expected["sway_yaw_root_1_re"]) <= ROOT_TOLERANCE, roots
    assert abs(roots[1].real - expected["sway_yaw_root_2_re"]) <= ROOT_TOLERANCE, roots
    assert abs(max_re - expected["max_re_motion"]) <= ROOT_TOLERANCE
    mean, _ = _turn_rate(vessel, _for_vehicle(vessel, fp.rudder_command(part, RUDDER)))
    assert abs(mean - expected["turn_rate_mean"]) <= TURN_TOLERANCE, (mean, expected["turn_rate_mean"])
    return roots, mean


def _workspace_vessel(configuration):
    vc.plugins()
    from rpp_py.context_builder import ComponentContextBuilder
    from rpp_py.data_manager import DataManager
    from rpp_plugin_registrator.library_manager import LibraryManager
    if configuration not in json.loads(DESCRIPTION.read_text())["Configurations"]:
        pytest.fail(f"vehicle_simulation.json has no configuration {configuration!r}")
    data_manager = DataManager(library_manager=LibraryManager(), workspace_path=str(LIBRARY))
    context = ComponentContextBuilder(data_manager=data_manager).build_script_from_description_path(
        str(DESCRIPTION), configuration=configuration)
    context.initialize()
    return context.get_component("vessels")[0]


def _part_folder(name):
    found = [d.parent for d in RPPWS.joinpath("parts").rglob("description.json")
             if json.loads(d.read_text()).get("Name") == name]
    if len(found) != 1:
        pytest.fail(f"named part {name!r}: {len(found)} found under .rppws/parts")
    return found[0]


def _part_tree(name):
    """``(root parameters, {slot: [(plugin name, parameters), ...]})`` of a named part, by its description ids."""
    from rpp_orchestrator.component_storage import ComponentParameterStore
    folder = _part_folder(name)
    record = json.loads((folder / "description.json").read_text())
    store = ComponentParameterStore()
    folders = {json.loads(d.read_text()).get("Id"): d.parent for d in RPPWS.joinpath("parts").rglob("description.json")}
    children = {}
    for slot, entry in record["Subcomponents"].items():
        for item in entry if isinstance(entry, list) else [entry]:
            target = folders[item["Id"]]
            linked = json.loads((target / "description.json").read_text()).get("LinkedComponentId")
            if linked is not None:  # a linked child holds its values in the named part it links
                target = folders[linked]
            children.setdefault(slot, []).append((item["PluginName"].split("::")[1], store.load(target)))
    return store.load(folder), children


def _random_nu_r(n, seed):
    g = np.random.default_rng(seed)
    nu = np.column_stack([g.uniform(0.2, 2.5, n), g.uniform(-0.4, 0.4, (n, 2)), g.uniform(-0.3, 0.3, (n, 3))])
    return nu


def _hull_alone(part, hull_options):
    p = vc.plugins()
    return _vehicle(part, hull_options, actuators=[p.PrescribedWrench])


def _signal(vessel, name, nu_r):
    x = np.concatenate([np.zeros(6), nu_r])
    return np.asarray(vessel.signals(x, np.zeros(vc.input_size(vessel)))[name]).ravel()


def _close(a, b, tolerance=FORCE_TOLERANCE):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return np.abs(a - b).max() <= tolerance * max(1.0, np.abs(b).max())


# ---------------------------------------------------------------------------------------------
# the reference itself (guards the frozen data and the published anchors)
# ---------------------------------------------------------------------------------------------
def test_reference_fin_lift_slope_from_geometry_reproduces_prestero_table_2_2():
    fin = fp.part_file()["fin"]
    slope = ref.low_aspect_ratio_lift_slope(fin["fin_span"]["value"], fin["fin_area"]["value"],
                                            fin["section_efficiency"]["value"], fin["aspect_ratio_factor"]["value"])
    aspect = fin["aspect_ratio_factor"]["value"] * fin["fin_span"]["value"] ** 2 / fin["fin_area"]["value"]
    assert abs(slope - 3.12) <= 0.005           # Prestero 2001, Table 2.2, p. 17: c_L_alpha 3.12
    assert abs(aspect - 2.21) <= 0.005          # Prestero 2001, Table 2.2, p. 17: AR_e 2.21


def test_reference_body_lift_reproduces_prestero_table_a1_and_the_sign_of_table_d1():
    """On Prestero's hull (l 1.33 m, d 0.191 m, Table 2.1 p. 15) the slender-body slope is Table A.1's
    1.20 /rad (p. 102); the point force at Table A.1's x_cp = -0.321 m gives M/Z = -x_cp, the ratio
    his linearised Table D.1 prints (M_wl / Z_wl = -11.1 / -34.5, p. 112) -- the sign of eq. 4.36
    as printed is the opposite."""
    slope = ref.slender_body_lift_slope(1.33, 0.191, 0.003)
    assert abs(slope - 1.20) <= 0.005
    tau = ref.point_force_body_lift(np.array([1.5, 0.0, 0.1, 0, 0, 0]), 1030.0, 0.191, slope, -0.321)
    assert abs(tau[4] / tau[2] - (-11.1 / -34.5)) <= 0.002
    tau_v = ref.point_force_body_lift(np.array([1.5, 0.1, 0.0, 0, 0, 0]), 1030.0, 0.191, slope, -0.321)
    assert np.isclose(tau_v[5] / tau_v[1], -0.321) and np.isclose(tau[4] / tau[2], 0.321)   # M_uwl = -N_uvl


def test_reference_paired_set_is_stable_only_with_the_geometric_fin_slope():
    for f in fp.PRESSURE_CENTRE_FRACTIONS:
        assert _reference("paired_geometric_fin_slope", f)["sway_yaw_root_1_re"] < 0
        assert _reference("paired_published_fin_slope", f)["sway_yaw_root_1_re"] > 0
    assert _reference("munk_kept_alone")["sway_yaw_root_1_re"] > 0
    assert _reference("published_forms")["sway_yaw_root_1_re"] < 0


def test_cruise_shaft_speed_of_the_part_file_holds_the_cruise_speed_in_the_reference_model():
    """The part file's cruise shaft speed is the one at which the reference model with the published forms runs straight
    at 1.54 m/s (Prestero 2001, Table 8.1, p. 60)."""
    n = fp.cruise_shaft_speed()
    assert abs(ref.shaft_speed_holding(fp.published_forms(), SPEED) - n) <= 1e-6
    part = fp.published_forms()
    model = ref.TorpedoAuv(part)
    x, residual = ref.straight_running(model.xdot, np.r_[np.zeros(len(part["fins"])), n],
                                       np.r_[np.zeros(6), SPEED, np.zeros(5)])
    assert residual <= 1e-9 and abs(x[6] - SPEED) <= 1e-9


# ---------------------------------------------------------------------------------------------
# the plugins against the reference where the forms are the reference model's own
# ---------------------------------------------------------------------------------------------
@pytest.mark.parametrize("case", ["published_forms", "munk_kept_alone", "geometric_fin_slope_effective_angle"])
def test_turn_tolerance_is_met_where_the_forms_exist(case):
    """Configurations the plugins build with the reference model's hull forms (the fin slope given its computed value)
    equal the reference: this fixes ``TURN_TOLERANCE`` and ``ROOT_TOLERANCE`` from measurement, not from the result
    wanted; the spinning one shows the stability gate can fail."""
    part = dict(fp.cases_by_name()[case])
    _assert_matches(_vehicle(part), part, case)


# ---------------------------------------------------------------------------------------------
# 1. "REMUS 100 (published reference)": the reference model's forms, Munk couplings removed
# ---------------------------------------------------------------------------------------------
def test_published_part_removes_the_munk_couplings():
    root, _ = _part_tree(PUBLISHED_PART)
    assert root["coriolis_form"] == "couplings_removed"


def test_published_part_declares_every_form_by_name():
    root, children = _part_tree(PUBLISHED_PART)
    (_, hull), = children["hydrodynamics"]
    expected = {"body_lift_form": "vertical_plane_wing", "hull_drag_form": "wing_drag", "roll_damping_form": "linear"}
    missing = {k: v for k, v in expected.items() if hull.get(k) != v}
    assert not missing, f"the hull-loads part of {PUBLISHED_PART!r} does not declare {missing}"


def test_published_part_goes_straight_and_turns_as_the_reference():
    vessel = _workspace_vessel(PUBLISHED_CONFIGURATION)
    part = fp.published_forms()
    roots, mean = _assert_matches(vessel, part, "published_forms")
    assert roots[0].real < 0 and mean > 0          # straight at zero fins; a positive rudder turns to starboard


# ---------------------------------------------------------------------------------------------
# the new forms, one at a time
# ---------------------------------------------------------------------------------------------
def test_point_force_body_lift_equals_the_reference_in_any_motion():
    part = fp.estimated()
    _require(vc.plugins().AuvHullLoads, HULL_OPTIONS[:7], "point-force body lift in both planes")
    vessel = _hull_alone(part, _hull_options(part))
    slope = ref.slender_body_lift_slope(part["length"], part["diameter"], part["body_lift_slope_per_degree"])
    x_cp = ref.pressure_centre_from_nose(part["length"] / 2, part["length"], part["pressure_centre_fraction"])
    wing = (part["span"], part["planform_area"], part["parasitic_drag_coefficient"], part["oswald_efficiency"])
    for nu_r in _random_nu_r(50, 1):
        lift = ref.point_force_body_lift(nu_r, part["water_density"], part["diameter"], slope, x_cp)
        assert _close(_signal(vessel, "hydrodynamics.body_lift_force", nu_r), lift)
        assert _close(_signal(vessel, "hydrodynamics.hull_drag_force", nu_r), ref.wing_drag(nu_r, part["water_density"], *wing))


def test_point_force_body_lift_has_the_sign_of_the_cross_product_not_of_the_printed_moment():
    """F_z = -k u w at x_cp < 0 aft of the origin: M = -x_cp F_z = x_cp k u w, N = x_cp F_y = -x_cp k u v with
    k = 1/2 rho d^2 c_ydb > 0 (Fossen 2011, eq. 12.226, p. 400): the moment ratios are M/Z = -x_cp and N/Y = x_cp."""
    part = fp.estimated()
    vessel = _hull_alone(part, _hull_options(part))
    x_cp = part["length"] / 2 - part["pressure_centre_fraction"] * part["length"]
    assert x_cp < 0
    heave = _signal(vessel, "hydrodynamics.body_lift_force", np.array([1.5, 0.0, 0.1, 0, 0, 0]))
    sway = _signal(vessel, "hydrodynamics.body_lift_force", np.array([1.5, 0.1, 0.0, 0, 0, 0]))
    assert heave[2] < 0 and abs(heave[4] / heave[2] + x_cp) <= 1e-12 * abs(x_cp)
    assert sway[1] < 0 and abs(sway[5] / sway[1] - x_cp) <= 1e-12 * abs(x_cp)
    assert np.all(heave[[0, 1, 3, 5]] == 0.0) and np.all(sway[[0, 2, 3, 4]] == 0.0)


def test_slender_body_slope_is_prestero_table_a1_on_his_hull():
    """c_ydb = (l/d) c_yb 180/pi, l 1.33 m, d 0.191 m, c_yb 0.003 per degree: Table A.1 prints 1.20 per rad (p. 102)."""
    from more_dynamics.models.lift_drag.body_lift import slender_body_lift_slope
    function = slender_body_lift_slope()
    slope = float(function(length=1.33, span=0.191, body_lift_slope_per_degree=0.003)["body_lift_slope"])
    assert abs(slope - 1.20) <= 0.005
    assert abs(slope - 1.33 / 0.191 * 0.003 * 180.0 / np.pi) <= 1e-12


@pytest.mark.parametrize("ratio, accepted", [(6.69, False), (6.71, True), (9.99, True), (10.01, False), (14.3, False)])
def test_slender_body_slope_is_refused_outside_its_stated_range(ratio, accepted):
    """The slope per degree is stated for 6.7 <= l/d <= 10 (Prestero 2001, eq. 4.31, p. 30): the vehicle is refused
    outside it, naming the ratio and the way out."""
    p = vc.plugins()
    part = dict(fp.estimated())
    part["length"] = ratio * part["diameter"]
    if accepted:
        _vehicle(part, _hull_options(part))
        return
    with pytest.raises(p.CompositionError, match=r"l/d = .*given"):
        _vehicle(part, _hull_options(part))


def test_a_hull_outside_the_slender_range_takes_the_given_slope():
    part = dict(fp.estimated())
    part["length"] = 14.3 * part["diameter"]
    given = {**_hull_options(part), "body_lift_slope_method": "given", "body_lift_slope": 1.3}
    vessel = _hull_alone(part, given)
    lift = _signal(vessel, "hydrodynamics.body_lift_force", np.array([1.5, 0.1, 0.1, 0, 0, 0]))
    assert np.all(np.isfinite(lift)) and lift[2] < 0


def test_vertical_plane_wing_is_the_lift_drag_block_split_in_two():
    part = fp.published_forms()
    _require(vc.plugins().AuvHullLoads, ("body_lift_form", "hull_drag_form"), "the lift-drag block split")
    vessel = _hull_alone(part, {"body_lift_form": "vertical_plane_wing", "hull_drag_form": "wing_drag"})
    wing = (part["span"], part["planform_area"], part["parasitic_drag_coefficient"], part["oswald_efficiency"])
    for nu_r in _random_nu_r(50, 2):
        lift = _signal(vessel, "hydrodynamics.body_lift_force", nu_r)
        drag = _signal(vessel, "hydrodynamics.hull_drag_force", nu_r)
        assert _close(lift, ref.vertical_plane_wing_lift(nu_r, part["water_density"], *wing))
        assert _close(drag, ref.wing_drag(nu_r, part["water_density"], *wing))
        assert _close(lift + drag, _signal(vessel, "hydrodynamics.lift_drag_force", nu_r))


def test_body_lift_slope_and_pressure_centre_given():
    part = fp.estimated()
    _require(vc.plugins().AuvHullLoads, HULL_OPTIONS[:7], "given body-lift slope and pressure centre")
    given = {**_hull_options(part), "body_lift_slope_method": "given", "body_lift_slope": 1.3,
             "pressure_centre_method": "given", "pressure_centre": -0.2}          # test values, any
    vessel = _hull_alone(part, given)
    for nu_r in _random_nu_r(20, 3):
        lift = ref.point_force_body_lift(nu_r, part["water_density"], part["diameter"], 1.3, -0.2)
        assert _close(_signal(vessel, "hydrodynamics.body_lift_force", nu_r), lift)


def test_quadratic_axial_hull_drag_uses_the_frontal_drag_coefficient():
    """X = -1/2 rho c_d (pi d^2/4) u|u| (Prestero 2001, eqs. 4.5-4.6, p. 25), c_d the declared
    ``hull_frontal_drag_coefficient`` -- no second drag value."""
    part = fp.published_forms()
    _require(vc.plugins().AuvHullLoads, ("hull_drag_form",), "quadratic axial hull drag")
    vessel = _hull_alone(part, {"hull_drag_form": "quadratic_axial"})
    for nu_r in _random_nu_r(20, 4):
        drag = ref.quadratic_axial_drag(nu_r, part["water_density"], part["diameter"], part["hull_frontal_drag_coefficient"])
        assert _close(_signal(vessel, "hydrodynamics.hull_drag_force", nu_r), drag)


def test_quadratic_axial_hull_drag_opposes_the_motion_ahead_and_astern():
    part = fp.published_forms()
    vessel = _hull_alone(part, {"hull_drag_form": "quadratic_axial"})
    ahead = _signal(vessel, "hydrodynamics.hull_drag_force", np.array([1.0, 0, 0, 0, 0, 0]))
    astern = _signal(vessel, "hydrodynamics.hull_drag_force", np.array([-1.0, 0, 0, 0, 0, 0]))
    assert ahead[0] < 0 < astern[0] and abs(ahead[0] + astern[0]) <= 1e-12 * abs(ahead[0])
    assert np.all(ahead[1:] == 0.0)


def test_quadratic_roll_damping_adds_only_k_pp_p_abs_p():
    """K = K_p|p| p|p| added to the linear roll term (Prestero 2001, eq. 4.10, p. 27); zero gives the
    linear form exactly; dissipative for K_p|p| <= 0. The coefficient is a test value: no identified one."""
    part = fp.published_forms()
    _require(vc.plugins().AuvHullLoads, ("roll_damping_form", "roll_quadratic_damping"), "quadratic roll damping")
    linear = _hull_alone(part, {"roll_damping_form": "linear"})
    zero = _hull_alone(part, {"roll_damping_form": "linear_plus_quadratic", "roll_quadratic_damping": 0.0})
    k_pp = -0.13                                                  # test value
    quadratic = _hull_alone(part, {"roll_damping_form": "linear_plus_quadratic", "roll_quadratic_damping": k_pp})
    for nu_r in _random_nu_r(20, 5):
        base = _signal(linear, "hydrodynamics.hydrodynamic_force", nu_r)
        assert np.array_equal(_signal(zero, "hydrodynamics.hydrodynamic_force", nu_r), base)
        delta = _signal(quadratic, "hydrodynamics.hydrodynamic_force", nu_r) - base
        assert _close(delta, ref.quadratic_roll_damping(k_pp, nu_r[3]))
        assert delta[3] * nu_r[3] <= 0


def test_a_quadratic_roll_coefficient_that_adds_energy_is_refused():
    part = fp.published_forms()
    with pytest.raises(ValueError, match="roll_quadratic_damping"):
        _hull_alone(part, {"roll_damping_form": "linear_plus_quadratic", "roll_quadratic_damping": 0.13})


def test_fin_lift_slope_from_geometry_equals_the_computed_value_given():
    p = vc.plugins()
    _require(p.Fin, FIN_OPTIONS, "fin lift slope from the fin geometry, Prestero 2001 eqs. 4.38-4.39")
    part = fp.estimated()
    fin = part["fins"][0]
    computed = _vehicle(part, actuators=[_fin(fin, _geometric_slope_selectors())])
    given = _vehicle(part, actuators=[_fin(fin, {"lift_law": "lift_slope"})])
    g = np.random.default_rng(6)
    for nu_r in _random_nu_r(20, 6):
        x = np.concatenate([np.zeros(6), nu_r])
        u = [g.uniform(-0.3, 0.3)]
        a = np.asarray(computed.signals(x, u)["actuators.0.generated_force"]).ravel()
        b = np.asarray(given.signals(x, u)["actuators.0.generated_force"]).ravel()
        assert _close(a, b)


def test_fin_lift_slope_from_geometry_is_the_closed_form_of_prestero_eqs_4_38_and_4_39():
    """a = [1/(2 pi a~) + 1/(pi AR_e)]^-1, AR_e = k_AR b^2 / S, with limits: a -> 2 pi a~ for a long fin and
    a -> pi AR_e when a~ is large (a slender fin)."""
    from more_dynamics.models.fin.lift_slope import low_aspect_ratio_lift_slope
    function = low_aspect_ratio_lift_slope()

    def slope(span, area, factor, efficiency):
        return float(function(fin_span=span, fin_area=area, aspect_ratio_factor=factor,
                              section_efficiency=efficiency)["lift_slope"])

    assert abs(slope(0.0857, 0.00665, 2.0, 0.9) - 3.1158) <= 5e-5                  # the estimated part's fin
    assert abs(slope(1e4, 1.0, 2.0, 0.9) - 2 * np.pi * 0.9) <= 1e-6                  # AR_e -> infinity
    assert abs(slope(0.1, 1.0, 2.0, 1e9) / (np.pi * 2.0 * 0.01) - 1.0) <= 1e-6      # a~ -> infinity
    assert abs(slope(0.3, 0.1, 1.0, 0.9) - 1.0 / (1 / (2 * np.pi * 0.9) + 1 / (np.pi * 0.9))) <= 1e-12


@pytest.mark.parametrize("form", ["helmbold", "lifting_line"])
def test_the_lift_slope_forms_that_wait_on_their_sources_are_refused(form):
    with pytest.raises(NotImplementedError, match="waits on source"):
        _vehicle(fp.estimated(), actuators=[_fin(fp.estimated()["fins"][0], {**_geometric_slope_selectors(),
                                                                             "lift_slope_form": form})])


def test_given_inertia_with_the_spheroids_own_tensor_equals_the_homogeneous_spheroid():
    """``rigid_body_form = given_inertia`` with the homogeneous spheroid's own tensor about the CO
    (I_o = I_g - m S(r_g)^2, Fossen 2011, eq. 3.34, p. 50) gives the same dynamics."""
    p = vc.plugins()
    _require(p.SpheroidAuv, ("rigid_body_form", "inertia_diagonal", "inertia_products"), "given inertia tensor")
    part = fp.published_forms()
    i_g = ref.spheroid_inertia_about_cg(part["body_mass"], part["length"] / 2, part["diameter"] / 2)
    i_o = ref.inertia_about_co(part["body_mass"], i_g, part["body_center_of_gravity"])
    given = {"rigid_body_form": "given_inertia", "inertia_diagonal": list(np.diag(i_o)),
             "inertia_products": [-i_o[0, 1], -i_o[1, 2], -i_o[0, 2]]}
    a = _vehicle(part, vehicle_options={"rigid_body_form": "homogeneous_spheroid"})
    b = _vehicle(part, vehicle_options=given)
    fa, fb = vc.graph_of(a).step, vc.graph_of(b).step
    g = np.random.default_rng(7)
    for _ in range(20):
        x = np.concatenate([g.uniform(-0.3, 0.3, 6), _random_nu_r(1, int(g.integers(1 << 30)))[0]])
        u = np.r_[g.uniform(-0.3, 0.3, 2), g.uniform(0, 1500)]
        assert _close(fa(x, u), fb(x, u))


def test_given_inertia_places_the_products_of_inertia_with_the_sign_of_the_tensor():
    """I_o = [[Ix, -Ixy, -Ixz], [-Ixy, Iy, -Iyz], [-Ixz, -Iyz, Iz]] (Fossen 2011, eq. 3.19, p. 48) fills the lower block of
    M_RB = [[m I3, -m S(r_g)], [m S(r_g), I_o]] (eq. 3.34, p. 50): a product Ixz = 2 puts -2 at (roll, yaw) and (yaw, roll)."""
    part = fp.published_forms()
    diagonal, products = [3.0, 40.0, 41.0], [0.5, 0.25, 2.0]
    vessel = _vehicle(part, vehicle_options={"rigid_body_form": "given_inertia", "inertia_diagonal": diagonal,
                                             "inertia_products": products})
    m_rb = _signal(vessel, "rigid_body_mass_matrix", np.zeros(6)).reshape(6, 6, order="F")
    r = np.asarray(part["body_center_of_gravity"], float)
    skew = np.array([[0, -r[2], r[1]], [r[2], 0, -r[0]], [-r[1], r[0], 0]])
    expected = np.block([[part["body_mass"] * np.eye(3), -part["body_mass"] * skew],
                         [part["body_mass"] * skew, np.array([[3.0, -0.5, -2.0], [-0.5, 40.0, -0.25], [-2.0, -0.25, 41.0]])]])
    assert np.abs(m_rb - expected).max() <= 1e-12 * np.abs(expected).max()
    assert m_rb[3, 5] == -2.0 and m_rb[5, 3] == -2.0


def test_given_inertia_that_is_not_positive_definite_is_refused():
    part = fp.published_forms()
    with pytest.raises(ValueError):
        _vehicle(part, vehicle_options={"rigid_body_form": "given_inertia", "inertia_diagonal": [1.0, 1.0, 1.0],
                                        "inertia_products": [0.0, 0.0, 5.0]})


# ---------------------------------------------------------------------------------------------
# 2-3. "REMUS 100 (estimated)": straight-line stable, turns as the reference
# ---------------------------------------------------------------------------------------------
@pytest.mark.parametrize("fraction", fp.PRESSURE_CENTRE_FRACTIONS)
def test_paired_set_is_straight_line_stable_and_turns_as_the_reference(fraction):
    part = fp.estimated(fraction)
    roots, mean = _assert_matches(_paired_vehicle(part), part, "paired_geometric_fin_slope", fraction)
    assert roots[0].real < 0 and roots[1].real < 0 and mean > 0


def test_paired_set_with_the_published_fin_coefficients_is_found_unstable():
    """The stability gate must be able to fail: the same pairing with the published reference model's calibrated fin
    coefficients has a sway-yaw root in the right half-plane."""
    part = fp.paired_published_fin_slope(0.65)
    roots, _ = _assert_matches(_paired_vehicle(part), part, "paired_published_fin_slope", 0.65)
    assert roots[0].real > 0


def test_estimated_named_part_holds_its_parameter_file():
    root, children = _part_tree(ESTIMATED_PART)
    record = fp.part_file()
    for name, entry in record["vehicle"].items():
        assert root.get(name) == entry["value"], name
    (_, hull), = children["hydrodynamics"]
    for name, entry in record["hydrodynamics"].items():
        assert hull.get(name) == entry["value"], name
    fins = [params for plugin, params in children["actuators"] if plugin == "Fin"]
    assert len(fins) == len(record["fins"])
    for params, entry in zip(fins, record["fins"]):
        assert list(params["lift_axis"]) == entry["lift_axis"]
        assert list(params["fin_position"]) == [-root["length"] / 2, 0, 0]
        for name in ("fin_area", "lift_law", "lift_slope_form", "fin_span", "aspect_ratio_factor", "section_efficiency",
                     "inflow", "flow_angle", "section", "servo_dynamics", "servo_angle_limit"):
            assert params.get(name) == record["fin"][name]["value"], name
    (propulsor,) = [params for plugin, params in children["actuators"] if plugin == "Propulsor"]
    assert propulsor["shaft_torque_scale"] == record["propeller"]["roll_moment_scale"]["value"]


def test_estimated_named_part_from_the_workspace_is_stable_and_turns_as_the_reference():
    vessel = _workspace_vessel(ESTIMATED_CONFIGURATION)
    part = fp.estimated()
    roots, mean = _assert_matches(vessel, part, "paired_geometric_fin_slope", part["pressure_centre_fraction"])
    assert roots[0].real < 0 and mean > 0


# Every named part that keeps every added-mass Coriolis term passes the straight-line test at its cruise
# state. A part's cruise state is listed here by its name; a new such part without one fails.
CRUISE = {ESTIMATED_PART: (ESTIMATED_CONFIGURATION, fp.estimated)}


def test_every_named_part_keeping_the_munk_moment_is_straight_line_stable():
    description = json.loads(DESCRIPTION.read_text())
    ids = {json.loads(d.read_text()).get("Id"): json.loads(d.read_text()) for d in RPPWS.joinpath("parts").rglob("description.json")}
    checked = []
    for key, configuration in description["Configurations"].items():
        for vessel_entry in configuration["Components"].get("vessels", []):
            record = ids.get(vessel_entry["Id"], {})
            if record.get("PluginName") != "more_dynamics::SpheroidAuv":
                continue
            root, _ = _part_tree(record["Name"])
            if root.get("coriolis_form", "full") != "full":
                continue
            assert record["Name"] in CRUISE, f"{record['Name']!r} keeps the Munk moment: give its cruise state here"
            configuration_key, part = CRUISE[record["Name"]]
            values = part()
            _, roots, _ = _roots(_workspace_vessel(configuration_key), np.array(fp.rudder_command(values, 0.0), float))
            assert roots[0].real < 0 and roots[1].real < 0, (record["Name"], roots)
            checked.append(record["Name"])
    assert ESTIMATED_PART in checked, f"{ESTIMATED_PART!r} is not in the workspace"


# ---------------------------------------------------------------------------------------------
# 4. Prestero's field turn: printed beside ours, never a gate
# ---------------------------------------------------------------------------------------------
def test_informational_turn_rate_beside_prestero_field_record(capsys):
    """Not a gate. Prestero's only horizontal-plane field record (Fig. 7-7, p. 57) is about 10 deg/s
    after a 4-degree rudder (his text), five seconds long, depth under closed-loop control, fin
    misalignment up to 5 degrees (p. 61); his simulator is straight-line unstable. Ours: the reference
    steady rate of the paired set, open loop."""
    for f in fp.PRESSURE_CENTRE_FRACTIONS:
        rate = np.rad2deg(_reference("paired_geometric_fin_slope", f)["turn_rate_mean"])
        print(f"pressure centre {f:.2f} l: {rate:.2f} deg/s (Prestero field record: about 10 deg/s, not comparable)")
    print(f"published forms: {np.rad2deg(_reference('published_forms')['turn_rate_mean']):.2f} deg/s")


# ---------------------------------------------------------------------------------------------
# 5. the sign of a fin deflection: the lift axis, no second label
# ---------------------------------------------------------------------------------------------
def test_fin_declares_no_convention_label():
    """The sign of a deflection is the fin's ``lift_axis`` alone (one source of one sign)."""
    p = vc.plugins()
    assert "convention" not in {d.name for d in p.Fin.PARAMETERS}


@pytest.mark.parametrize("name", [PUBLISHED_PART, ESTIMATED_PART])
def test_every_remus_part_fin_is_starboard_down_positive_by_its_lift_axis(name):
    """A positive rudder deflection lifts along -y (trailing edge to starboard) and a positive stern-plane deflection
    along -z (trailing edge down): the lift axis of every fin has its one non-zero component equal to -1."""
    _, children = _part_tree(name)
    fins = [params for plugin, params in children["actuators"] if plugin == "Fin"]
    assert len(fins) == 4
    for params in fins:
        axis = np.asarray(params["lift_axis"], float)
        assert axis[0] == 0.0 and sorted(np.abs(axis[1:])) == [0.0, 1.0] and axis[1] + axis[2] == -1.0, list(axis)


def test_fin_pairs_deflection_only_declares_the_convention_of_its_model():
    p = vc.plugins()
    declared = {d.name: d.default_value for d in p.FinPairsDeflectionOnly.PARAMETERS}
    assert declared["convention"] == "starboard_down_positive"
    part = fp.published_forms()
    vessel = _vehicle(part, actuators=[(p.FinPairsDeflectionOnly, vc.pick(p.FinPairsDeflectionOnly, part)),
                                       (p.Propeller, vc.pick(p.Propeller, part))])
    state = np.concatenate([np.zeros(6), [1.5, 0, 0, 0, 0, 0]])
    rudder = np.asarray(vessel.signals(state, [0.1, 0.0, 0.0])["actuators.0.generated_force"]).ravel()
    stern = np.asarray(vessel.signals(state, [0.0, 0.1, 0.0])["actuators.0.generated_force"]).ravel()
    # a positive rudder lifts along -y at the tail (x < 0): N = x Y > 0, a turn to starboard; a positive stern plane lifts
    # along -z (upward): M = -x Z < 0, bow down (Fossen 2011, eq. 12.226, p. 400)
    assert rudder[1] < 0 and rudder[5] > 0 and stern[2] < 0 and stern[4] < 0


def test_fin_pairs_deflection_only_refuses_another_convention():
    p = vc.plugins()
    part = fp.published_forms()
    params = {**vc.pick(p.FinPairsDeflectionOnly, part), "convention": "port_up_positive"}
    with pytest.raises(ValueError, match="convention"):
        _vehicle(part, actuators=[(p.FinPairsDeflectionOnly, params), (p.Propeller, vc.pick(p.Propeller, part))])


def test_a_positive_rudder_command_turns_the_published_part_to_starboard():
    vessel = _workspace_vessel(PUBLISHED_CONFIGURATION)
    mean, _ = _turn_rate(vessel, _for_vehicle(vessel, fp.rudder_command(fp.published_forms(), RUDDER)))
    assert mean > 0


# ---------------------------------------------------------------------------------------------
# 6. generality: the same option set on a second torpedo geometry (synthetic, no REMUS number)
# ---------------------------------------------------------------------------------------------
SYNTHETIC = {   # a synthetic torpedo, labelled synthetic: no identified second set with these forms exists
    "length": 2.0, "diameter": 0.24, "body_mass": 88.0, "body_center_of_gravity": [0, 0, 0.015],
    "roll_added_inertia_ratio": 0.2, "latitude": 1.0, "water_density": 1020.0, "kinematic_viscosity": 1.1e-6,
    "center_of_buoyancy": [0, 0, 0], "time_constants": [15, 15, 2], "damping_ratios": [0.4, 0.6],
    "oswald_efficiency": 0.4, "planform_fraction": 0.65, "hull_frontal_drag_coefficient": 0.35,
    "coriolis_form": "full", "body_lift_form": "point_force_both_planes",
    "body_lift_slope_method": "slender_body_hoerner", "body_lift_slope_per_degree": 0.003,   # Hoerner's constant (Prestero 2001, eq. 4.31), not a vehicle number
    "pressure_centre_method": "fraction_of_length_from_nose", "pressure_centre_fraction": 0.62,
    "hull_drag_form": "quadratic_axial", "roll_damping_form": "linear",
}


def test_the_option_set_builds_and_runs_on_a_second_torpedo_geometry():
    p = vc.plugins()
    _require(p.AuvHullLoads, HULL_OPTIONS[:8], "the point-force body lift and axial drag options")
    _require(p.Fin, FIN_OPTIONS, "fin lift slope from the fin geometry")
    s = SYNTHETIC
    hull = {k: s[k] for k in ("time_constants", "damping_ratios", "oswald_efficiency", "planform_fraction",
                              "hull_frontal_drag_coefficient")}
    hull.update({k: s[k] for k in HULL_OPTIONS[:8] if k in s and k not in ("body_lift_slope", "pressure_centre")})
    fins = []
    for axis in ([0, -1, 0], [0, -1, 0], [0, 0, -1], [0, 0, -1]):
        fins.append(_fin({"fin_position": [-s["length"] / 2, 0, 0], "chord_axis": [1, 0, 0], "lift_axis": axis,
                          "fin_area": 0.011, "max_deflection": 0.4, "inflow": "rigid_point", "flow_angle": "small_angle",
                          "lift_slope": 0.0}, {"lift_law": "finite_wing", "lift_slope_form": "low_aspect_ratio",
                                               "fin_span": 0.11, "aspect_ratio_factor": 2.0,
                                               "section_efficiency": 0.85}))
    vehicle = {k: s[k] for k in ("length", "diameter", "body_mass", "body_center_of_gravity", "roll_added_inertia_ratio",
                                 "coriolis_form")}
    vessel = vc.build_context(
        p.SpheroidAuv, {**vehicle, "current_form": "none"},
        site=(p.Site, {"site_form": "latitude", "latitude": s["latitude"], "water_density": s["water_density"],
                       "kinematic_viscosity": s["kinematic_viscosity"]}),
        hydrostatics=(p.SubmergedRestoring, {"center_of_buoyancy": s["center_of_buoyancy"]}),
        hydrodynamics=(p.AuvHullLoads, hull), actuators=fins + [p.PrescribedWrench])
    remus = _paired_vehicle(fp.estimated())
    x = np.zeros(12)
    x[6:12] = [1.2, 0.05, -0.03, 0.02, 0.01, 0.04]
    u = np.r_[np.full(4, 0.05), [30.0, 0, 0, 0, 0, 0]]
    xdot = np.asarray(vc.graph_of(vessel).step(x, u)).ravel()
    assert np.all(np.isfinite(xdot))
    assert np.all(np.isfinite(vessel.step_vector(x, u, 1.0)))
    u_remus = np.r_[np.full(4, 0.05), fp.cruise_shaft_speed()]
    assert not np.allclose(xdot, np.asarray(vc.graph_of(remus).step(x, u_remus)).ravel())
