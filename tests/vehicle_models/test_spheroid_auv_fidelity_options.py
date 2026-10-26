"""Gates of the higher-fidelity options of ``SpheroidAuv``, ``AuvHullLoads``
and ``Fin``, and of the two REMUS 100 named parts that select them.

Written before the options exist: every test that needs a missing option or
part fails with a message naming it. The reference is the numpy second
implementation in ``data/remus100_fidelity/`` (``reference_torpedo_auv.py``,
parameter sets in ``fidelity_parts.py``), whose MSS forms equal the frozen
one-value MSS reference of ``data/remus100/``, and whose frozen turn and
stability table (``remus100_fidelity_turn_and_stability.csv``) its own
generator writes. Published numbers are typed with their page (Prestero
2001; Fossen 2011).

Contract (parameter names; the porter implements them):

* ``AuvHullLoads``: ``body_lift_form`` (``vertical_plane_wing`` = today's
  lift half of the lift-drag block | ``point_force_both_planes``) with
  ``body_lift_slope_method`` (``slender_body_hoerner`` reading
  ``body_lift_slope_per_degree`` | ``given`` reading ``body_lift_slope``) and
  ``pressure_centre_method`` (``fraction_of_length_from_nose`` reading
  ``pressure_centre_fraction``, nose at +L/2 | ``given`` reading
  ``pressure_centre``); ``hull_drag_form`` (``wing_drag`` | ``quadratic_axial``
  reading ``hull_frontal_drag_coefficient``); ``roll_damping_form``
  (``linear`` | ``linear_plus_quadratic`` reading ``roll_quadratic_damping``);
  outputs ``body_lift_force`` and ``hull_drag_force`` (``lift_drag_force``
  stays their sum).
* ``Fin``: ``lift_slope_method`` (``given`` reading ``lift_slope`` |
  ``low_aspect_ratio`` reading ``fin_span``, ``fin_area``,
  ``section_lift_efficiency``); ``convention`` (the four names of
  ``models/fin/fins.py``), also on ``FinPairsDeflectionOnly``.
* ``SpheroidAuv``: ``rigid_body_form`` (``homogeneous_spheroid`` |
  ``given_inertia`` reading ``inertia_diagonal`` and ``inertia_products``
  about the CO, the names and convention of ``DerivativeAuv``).
* Workspace: configuration ``remus100_munk_body_lift`` of
  ``vehicle_simulation.json``, named part "REMUS 100 (Munk + body lift)" with
  the values and selectors of ``remus100_munk_body_lift_parameters.json``.

Not contracted here (not yet derived or settled; stated, never guessed):
the hull-profile, strip-theory added-mass and profile cross-flow damping
forms, and a part built on Prestero's own geometry and coefficients; the
behaviour of ``slender_body_hoerner`` outside 6.7 <= l/d <= 10.

Tolerances. Derivatives and forces: 1e-12 relative (the same equations in
another summation order). Turn rates and roots: the plugins' CasADi graphs
against the numpy reference over 1200 RK4 steps, measured on the three
configurations that build today (2026-10-13): states within 2.2e-16 after
2 s, the 10 s mean yaw rate equal to the last bit, so 1e-9 rad/s (the MSS
gates' derivative tolerance) is a bound with margin; roots within 2.1e-7 1/s
(the reference's central difference, step 1e-7, against CasADi's exact
Jacobian), so 1e-6. Published numbers: half a unit of the last printed
digit.

References
----------
Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics and Motion
    Control. Wiley. Ch. 3, eq. 3.24, p. 49; ch. 12, eqs. 12.17-12.24,
    p. 351; eq. 12.226, p. 400.
Prestero, T. (2001). Verification of a six-degree of freedom simulation
    model for the REMUS autonomous underwater vehicle. MSc thesis, MIT/WHOI.
    Tables 2.1-2.2, pp. 15-17; eqs. 4.5-4.6, 4.10, 4.30-4.39, pp. 25-31;
    Fig. 7-7, p. 57; Table 8.1, p. 60; Table A.1, p. 102; Table D.1, p. 112.
Fossen, T. I. MSS, MIT, @ ed2ce7b: CRAFT/AUV/models/remus100.m 208-211,
    243, 246, 255; CRAFT/AUV/SIMremus100.m 60.

Author:    Enio Krizman
Date:      2026-10-13
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
MSS_PART, MSS_CONFIGURATION = "REMUS 100 (MSS)", "remus100"
PAIRED_PART = fp.part_file()["part"]
PAIRED_CONFIGURATION = fp.part_file()["configuration"]
REFERENCE_FILE = FIDELITY / "remus100_fidelity_turn_and_stability.csv"
FORCE_TOLERANCE = 1e-12          # relative, same equations
TURN_TOLERANCE = 1e-9            # rad/s on the 10 s mean yaw rate (measured: equal to the last bit)
ROOT_TOLERANCE = 1e-6            # 1/s (measured 2.1e-7: CasADi Jacobian against central differences, eps 1e-7)
STEP = 0.05                      # SIMremus100.m 60
DURATION, AVERAGE = 60.0, 10.0
SPEED = 1.54                     # Prestero 2001, Table 8.1, p. 60
RUDDER = np.deg2rad(4.0)         # Prestero 2001, p. 53
CONVENTION_AXES = {"starboard": -1.0, "port": 1.0, "down": -1.0, "up": 1.0}   # lift axis sign (models/fin/fins.py 50-55)

HULL_OPTIONS = ("body_lift_form", "body_lift_slope_method", "body_lift_slope_per_degree", "body_lift_slope",
                "pressure_centre_method", "pressure_centre_fraction", "pressure_centre", "hull_drag_form",
                "roll_damping_form", "roll_quadratic_damping")
FIN_OPTIONS = ("lift_slope_method", "fin_span", "section_lift_efficiency")


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
    """Fail (red) naming the options ``cls`` does not declare yet."""
    missing = [n for n in names if n not in {d.name for d in cls.PARAMETERS}]
    if missing:
        pytest.fail(f"{cls.__name__} does not declare {missing} yet ({why})")


def _fin(fin, selectors=None):
    """A ``Fin`` child from a reference fin set (servo none, angle limit on the command: remus100.m 114-115)."""
    p = vc.plugins()
    params = {"fin_position": list(fin["fin_position"]), "chord_axis": list(fin["chord_axis"]),
              "lift_axis": list(fin["lift_axis"]), "fin_area": fin["fin_area"], "max_deflection": fin["max_deflection"],
              "servo_dynamics": "none", "servo_rate_limit": False, "servo_angle_limit": "on_command",
              "inflow": fin["inflow"], "flow_angle": fin["flow_angle"], "interference": "none",
              "section": "quadratic_drag", "lift_slope": fin["lift_slope"],
              "time_constant": 1.0, "max_rate": 1.0, "zero_lift_drag": 0.0}   # read by no selected form
    params.update(selectors or {})
    return (p.Fin, params)


def _fin_pairs_or_fins(part):
    """MSS's two pairs as ``FinPairsDeflectionOnly`` (two fin entries) or one ``Fin`` per fin."""
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
    """The part's forms through the plugins (red until the options exist)."""
    p = vc.plugins()
    _require(p.AuvHullLoads, ("body_lift_form", "body_lift_slope_method", "body_lift_slope_per_degree",
                              "pressure_centre_method", "pressure_centre_fraction", "hull_drag_form", "roll_damping_form"),
             "point-force body lift in both planes, Prestero 2001 eqs. 4.27-4.35")
    record = fp.part_file()["fin"]
    fins = [_fin(f, {"lift_slope_method": "low_aspect_ratio", "fin_span": record["fin_span"]["value"],
                     "section_lift_efficiency": record["section_lift_efficiency"]["value"]})
            for f in part["fins"]]
    if part["fins"][0]["lift_slope"] < 1.0:     # the calibrated MSS coefficients: given
        fins = [_fin(f) for f in part["fins"]]
    else:
        _require(p.Fin, FIN_OPTIONS, "fin lift slope from the fin geometry, Prestero 2001 eqs. 4.38-4.39")
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
    """The command of a reference part for ``vessel``: the reference part of the MSS forms gives one command per fin
    pair (rudder, stern plane) and the shaft speed; the REMUS 100 (MSS) part has one actuator per fin (the rudder
    pair and the stern-plane pair, two fins each), so each pair's command is given to both of its fins."""
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
        pytest.fail(f"vehicle_simulation.json has no configuration {configuration!r} yet")
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
            target = folders[item.get("LinkedComponentId", item["Id"])]
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
# the reference itself (green: guards the frozen data and the published anchors)
# ---------------------------------------------------------------------------------------------
def test_reference_mss_forms_turn_equals_the_frozen_mss_run():
    """The numpy reference with every MSS form turns as the one-value MSS model: the last 10 s of
    ``rudder_step_4deg_at_1p54ms_consistent.csv`` (rudder from 10 s; both steady after 50 s)."""
    with open(vc.DATA / "remus100" / "remus100_trajectory_rudder_step_4deg_at_1p54ms_consistent.csv", newline="") as h:
        rows = [r for r in csv.DictReader(h) if float(r["t"]) > DURATION - AVERAGE]
    frozen = np.mean([float(r["x_06"]) for r in rows])          # MSS state order [nu; eta]: x_06 = r
    assert abs(_reference("mss_forms")["turn_rate_mean"] - frozen) <= 1e-9


def test_reference_fin_lift_slope_from_geometry_reproduces_prestero_table_2_2():
    fin = fp.part_file()["fin"]
    slope = ref.low_aspect_ratio_lift_slope(fin["fin_span"]["value"], fin["fin_area"]["value"],
                                            fin["section_lift_efficiency"]["value"])
    aspect = 2 * fin["fin_span"]["value"] ** 2 / fin["fin_area"]["value"]
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
        assert _reference("paired_mss_fin_slope", f)["sway_yaw_root_1_re"] > 0
    assert _reference("munk_kept_alone")["sway_yaw_root_1_re"] > 0
    assert _reference("mss_forms")["sway_yaw_root_1_re"] < 0


# ---------------------------------------------------------------------------------------------
# the plugins against the reference where the forms exist today (green: proves the harness)
# ---------------------------------------------------------------------------------------------
@pytest.mark.parametrize("case", ["mss_forms", "munk_kept_alone", "geometric_fin_slope_effective_angle"])
def test_turn_tolerance_is_met_where_the_forms_exist(case):
    """Configurations the plugins can build today (the fin slope given its computed value) equal the
    reference: this fixes ``TURN_TOLERANCE`` and ``ROOT_TOLERANCE`` from measurement, not from the result
    wanted; the spinning one shows the stability gate can fail."""
    part = dict(fp.cases_by_name()[case])
    _assert_matches(_vehicle(part), part, case)


# ---------------------------------------------------------------------------------------------
# 1. "REMUS 100 (MSS)": MSS forms, Munk couplings removed (remus100.m 208-211)
# ---------------------------------------------------------------------------------------------
def test_mss_part_removes_the_munk_couplings():
    root, _ = _part_tree(MSS_PART)
    assert root["coriolis_form"] == "munk_couplings_removed"


def test_mss_part_declares_every_mss_form_by_name():
    root, children = _part_tree(MSS_PART)
    (_, hull), = children["hydrodynamics"]
    expected = {"body_lift_form": "vertical_plane_wing", "hull_drag_form": "wing_drag", "roll_damping_form": "linear"}
    missing = {k: v for k, v in expected.items() if hull.get(k) != v}
    assert not missing, f"the hull-loads part of {MSS_PART!r} does not declare {missing}"


def test_mss_part_goes_straight_and_turns_as_mss():
    vessel = _workspace_vessel(MSS_CONFIGURATION)
    part = fp.mss_forms()
    roots, mean = _assert_matches(vessel, part, "mss_forms")
    assert roots[0].real < 0 and mean > 0          # straight at zero fins; a positive rudder turns to starboard


# ---------------------------------------------------------------------------------------------
# the new forms, one at a time (red until the options exist)
# ---------------------------------------------------------------------------------------------
def test_point_force_body_lift_equals_the_reference_in_any_motion():
    part = fp.munk_body_lift()
    _require(vc.plugins().AuvHullLoads, HULL_OPTIONS[:7], "point-force body lift in both planes")
    vessel = _hull_alone(part, _hull_options(part))
    slope = ref.slender_body_lift_slope(part["length"], part["diameter"], part["body_lift_slope_per_degree"])
    x_cp = ref.pressure_centre_from_nose(part["length"] / 2, part["length"], part["pressure_centre_fraction"])
    wing = (part["span"], part["planform_area"], part["parasitic_drag_coefficient"], part["oswald_efficiency"])
    for nu_r in _random_nu_r(50, 1):
        lift = ref.point_force_body_lift(nu_r, part["water_density"], part["diameter"], slope, x_cp)
        assert _close(_signal(vessel, "hydrodynamics.body_lift_force", nu_r), lift)
        assert _close(_signal(vessel, "hydrodynamics.hull_drag_force", nu_r), ref.wing_drag(nu_r, part["water_density"], *wing))


def test_vertical_plane_wing_is_todays_lift_drag_block_split_in_two():
    part = fp.mss_forms()
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
    part = fp.munk_body_lift()
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
    part = fp.mss_forms()
    _require(vc.plugins().AuvHullLoads, ("hull_drag_form",), "quadratic axial hull drag")
    vessel = _hull_alone(part, {"hull_drag_form": "quadratic_axial"})
    for nu_r in _random_nu_r(20, 4):
        drag = ref.quadratic_axial_drag(nu_r, part["water_density"], part["diameter"], part["hull_frontal_drag_coefficient"])
        assert _close(_signal(vessel, "hydrodynamics.hull_drag_force", nu_r), drag)


def test_quadratic_roll_damping_adds_only_k_pp_p_abs_p():
    """K = K_p|p| p|p| added to the linear roll term (Prestero 2001, eq. 4.10, p. 27); zero gives the
    linear form exactly; dissipative for K_p|p| <= 0. The coefficient is a test value: no identified one."""
    part = fp.mss_forms()
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


def test_fin_lift_slope_from_geometry_equals_the_computed_value_given():
    p = vc.plugins()
    _require(p.Fin, FIN_OPTIONS, "fin lift slope from the fin geometry, Prestero 2001 eqs. 4.38-4.39")
    part = fp.munk_body_lift()
    record = fp.part_file()["fin"]
    fin = part["fins"][0]
    computed = _vehicle(part, actuators=[_fin(fin, {"lift_slope_method": "low_aspect_ratio",
                                                    "fin_span": record["fin_span"]["value"],
                                                    "section_lift_efficiency": record["section_lift_efficiency"]["value"]})])
    given = _vehicle(part, actuators=[_fin(fin, {"lift_slope_method": "given"})])
    g = np.random.default_rng(6)
    for nu_r in _random_nu_r(20, 6):
        x = np.concatenate([np.zeros(6), nu_r])
        u = [g.uniform(-0.3, 0.3)]
        a = np.asarray(computed.signals(x, u)["actuators.0.generated_force"]).ravel()
        b = np.asarray(given.signals(x, u)["actuators.0.generated_force"]).ravel()
        assert _close(a, b)


def test_given_inertia_with_the_spheroids_own_tensor_equals_the_homogeneous_spheroid():
    """``rigid_body_form = given_inertia`` with the homogeneous spheroid's own tensor about the CO
    (I_o = I_g - m S(r_g)^2, Fossen 2011, eq. 3.24, p. 49) gives the same dynamics."""
    p = vc.plugins()
    _require(p.SpheroidAuv, ("rigid_body_form", "inertia_diagonal", "inertia_products"), "given inertia tensor")
    part = fp.mss_forms()
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


# ---------------------------------------------------------------------------------------------
# 2-3. "REMUS 100 (Munk + body lift)": straight-line stable, turns as the reference
# ---------------------------------------------------------------------------------------------
@pytest.mark.parametrize("fraction", fp.PRESSURE_CENTRE_FRACTIONS)
def test_paired_set_is_straight_line_stable_and_turns_as_the_reference(fraction):
    part = fp.munk_body_lift(fraction)
    roots, mean = _assert_matches(_paired_vehicle(part), part, "paired_geometric_fin_slope", fraction)
    assert roots[0].real < 0 and roots[1].real < 0 and mean > 0


def test_paired_set_with_the_calibrated_fin_coefficients_is_found_unstable():
    """The stability gate must be able to fail: the same pairing with MSS's calibrated fin
    coefficients (remus100.m 183, 188) has a sway-yaw root in the right half-plane."""
    part = fp.paired_mss_fin_slope(0.65)
    roots, _ = _assert_matches(_paired_vehicle(part), part, "paired_mss_fin_slope", 0.65)
    assert roots[0].real > 0


def test_munk_body_lift_named_part_holds_its_parameter_file():
    root, children = _part_tree(PAIRED_PART)
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
        for name in ("fin_area", "lift_slope_method", "fin_span", "section_lift_efficiency", "inflow", "flow_angle",
                     "section", "servo_dynamics", "servo_angle_limit", "convention"):
            assert params.get(name) == record["fin"][name]["value"], name


def test_munk_body_lift_named_part_from_the_workspace_is_stable_and_turns_as_the_reference():
    vessel = _workspace_vessel(PAIRED_CONFIGURATION)
    part = fp.munk_body_lift()
    roots, mean = _assert_matches(vessel, part, "paired_geometric_fin_slope", part["pressure_centre_fraction"])
    assert roots[0].real < 0 and mean > 0


# Every named part that keeps every added-mass Coriolis term passes the straight-line test at its cruise
# state. A part's cruise shaft speed is listed here by its name; a new such part without one fails.
CRUISE = {PAIRED_PART: (PAIRED_CONFIGURATION, fp.munk_body_lift)}


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
            if root.get("coriolis_form", "kirchhoff_full") != "kirchhoff_full":
                continue
            assert record["Name"] in CRUISE, f"{record['Name']!r} keeps the Munk moment: give its cruise state here"
            configuration_key, part = CRUISE[record["Name"]]
            values = part()
            _, roots, _ = _roots(_workspace_vessel(configuration_key), np.array(fp.rudder_command(values, 0.0), float))
            assert roots[0].real < 0 and roots[1].real < 0, (record["Name"], roots)
            checked.append(record["Name"])
    assert PAIRED_PART in checked, f"{PAIRED_PART!r} is not in the workspace yet"


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
    print(f"MSS forms: {np.rad2deg(_reference('mss_forms')['turn_rate_mean']):.2f} deg/s")


# ---------------------------------------------------------------------------------------------
# 5. fin convention per part, by name (the four names of models/fin/fins.py)
# ---------------------------------------------------------------------------------------------
@pytest.mark.parametrize("name, expected", [(MSS_PART, "starboard_down_positive"),
                                            (PAIRED_PART, "starboard_down_positive")])
def test_every_remus_part_declares_its_fin_convention(name, expected):
    _, children = _part_tree(name)
    fins = [(plugin, params) for plugin, params in children["actuators"] if plugin in ("Fin", "FinPairsDeflectionOnly")]
    assert fins
    for plugin, params in fins:
        assert params.get("convention") == expected, f"{name!r}: its {plugin} declares no convention {expected!r}"
        if plugin == "Fin":
            rudder_word, plane_word = params["convention"].split("_")[:2]
            axis = np.asarray(params["lift_axis"], float)
            sign = CONVENTION_AXES[rudder_word] if abs(axis[1]) > 0 else CONVENTION_AXES[plane_word]
            assert axis[1] + axis[2] == sign, (params["convention"], list(axis))


def test_a_positive_rudder_command_turns_the_mss_part_to_starboard():
    vessel = _workspace_vessel(MSS_CONFIGURATION)
    mean, _ = _turn_rate(vessel, _for_vehicle(vessel, fp.rudder_command(fp.mss_forms(), RUDDER)))
    assert mean > 0


# ---------------------------------------------------------------------------------------------
# 6. generality: the same option set on a second torpedo geometry (synthetic, no REMUS number)
# ---------------------------------------------------------------------------------------------
SYNTHETIC = {   # a synthetic torpedo, labelled synthetic: no identified second set with these forms exists
    "length": 2.0, "diameter": 0.24, "body_mass": 88.0, "body_center_of_gravity": [0, 0, 0.015],
    "roll_added_inertia_ratio": 0.2, "latitude": 1.0, "water_density": 1020.0, "kinematic_viscosity": 1.1e-6,
    "center_of_buoyancy": [0, 0, 0], "time_constants": [15, 15, 2], "damping_ratios": [0.4, 0.6],
    "oswald_efficiency": 0.4, "planform_fraction": 0.65, "hull_frontal_drag_coefficient": 0.35,
    "coriolis_form": "kirchhoff_full", "body_lift_form": "point_force_both_planes",
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
                          "lift_slope": 0.0}, {"lift_slope_method": "low_aspect_ratio", "fin_span": 0.11,
                                               "section_lift_efficiency": 0.85}))
    vehicle = {k: s[k] for k in ("length", "diameter", "body_mass", "body_center_of_gravity", "roll_added_inertia_ratio",
                                 "coriolis_form")}
    vessel = vc.build_context(
        p.SpheroidAuv, {**vehicle, "current_form": "none"},
        site=(p.Site, {"site_form": "latitude", "latitude": s["latitude"], "water_density": s["water_density"],
                       "kinematic_viscosity": s["kinematic_viscosity"]}),
        hydrostatics=(p.SubmergedRestoring, {"center_of_buoyancy": s["center_of_buoyancy"]}),
        hydrodynamics=(p.AuvHullLoads, hull), actuators=fins + [p.PrescribedWrench])
    remus = _paired_vehicle(fp.munk_body_lift())
    x = np.zeros(12)
    x[6:12] = [1.2, 0.05, -0.03, 0.02, 0.01, 0.04]
    u = np.r_[np.full(4, 0.05), [30.0, 0, 0, 0, 0, 0]]
    xdot = np.asarray(vc.graph_of(vessel).step(x, u)).ravel()
    assert np.all(np.isfinite(xdot))
    assert np.all(np.isfinite(vessel.step_vector(x, u, 1.0)))
    u_remus = np.r_[np.full(4, 0.05), fp.cruise_shaft_speed()]
    assert not np.allclose(xdot, np.asarray(vc.graph_of(remus).step(x, u_remus)).ravel())
