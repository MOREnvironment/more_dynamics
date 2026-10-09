"""The options of the three vehicle types and of the plugins they are built from:
the numerical guard at rest (``smooth_speed``), the derived quantities that
are computed from their primitives by default and given by a named method
(``body_density_method``, ``planform_area_method``,
``parasitic_drag_coefficient_method``, ``fin_positions_method``,
``wetted_surface_method``), the site (``site_form``), the buoyancy of a
positively buoyant vehicle (``buoyancy_method``) and the current form
(``current_form``).

The reference numbers are closed forms of cited equations, written out in each
test, or the frozen REMUS 100 and ITTC reference files under the other test
folders; no value is typed from memory.

Author:    Enio Krizman
Date:      2026-10-09
"""

import csv
from pathlib import Path

import casadi as ca
import numpy as np
import pytest

import test_spheroid_auv as spheroid
import vehicle_contract as vc

DATA_HYDRODYNAMIC_LOADS = Path(__file__).resolve().parents[1] / "hydrodynamic_loads" / "data"
DEG = np.pi / 180.0
REST = np.zeros(12)


def _remus(full=True, **kw):
    return vc.remus(full=full, current="none", **kw)


def _jacobian_at_rest(vessel):
    graph = vc.graph_of(vessel)
    n_input = sum(d.size for d in graph.payload.inputDescription)
    state, command = ca.SX.sym("x", graph.num_states), ca.SX.sym("u", n_input)
    jacobian = ca.Function("J", [state, command], [ca.jacobian(graph.step(state, command), state)])
    return np.asarray(jacobian(np.zeros(graph.num_states), np.zeros(n_input)))


# ---------------------------------------------------------------------------
# The guard at rest
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("full", [False, True], ids=["hull", "full"])
def test_jacobian_at_zero_relative_velocity_is_finite(full):
    """The slender-body lift (``atan2`` of the flow angle) and the speed of the surge-damping fade have 0/0
    derivatives at rest; with the guard (the default) a stiff integrator can start there."""
    assert np.isfinite(_jacobian_at_rest(_remus(full))).all()


def test_the_mss_form_has_no_finite_jacobian_at_rest():
    """The reason for the guard: with ``smooth_speed = False`` (the MSS form, ``remus100.m`` 126-127, 218) eight
    entries of the Jacobian (rows u, w, q; columns u, v, w) are not finite at zero relative velocity."""
    jacobian = _jacobian_at_rest(_remus(values={**vc.remus_values(), "smooth_speed": False}))
    bad = {tuple(int(i) for i in entry) for entry in np.argwhere(~np.isfinite(jacobian))}
    assert bad == {(6, 6), (6, 7), (6, 8), (8, 6), (8, 8), (10, 6), (10, 7), (10, 8)}


@pytest.mark.parametrize("stage", [spheroid.HULL, spheroid.FULL], ids=["hull", "full"])
def test_the_guard_leaves_the_frozen_derivative_rows_unchanged(stage):
    """On every row of the frozen REMUS derivative reference (the rest rows included) the guarded derivative equals
    the unguarded one to 1e-12, far inside the reference's own 1e-9."""
    table = spheroid._read_csv(spheroid.DERIVATIVE_FILE)
    current = np.column_stack([table["Vc"], table["betaVc"], table["w_c"]])
    command = spheroid._columns(table, "ui", 3) if stage == spheroid.FULL else None
    tau = spheroid._columns(table, "tau", 6) if stage == spheroid.HULL else None
    x = spheroid._columns(table, "x", 12)
    values = spheroid._parameter_set(stage)
    guarded = spheroid._evaluate(spheroid._builder(stage, values=values), x, current, tau_ext=tau, u=command)[0]
    exact = spheroid._evaluate(
        spheroid._builder(stage, values={**values, "smooth_speed": False}), x, current, tau_ext=tau, u=command)[0]
    assert np.abs(guarded - exact).max() <= 1e-12


def test_remus_starts_from_rest_in_a_stiff_integrator_and_accelerates_under_thrust():
    """CVODES (as ``more_simulation`` runs a vehicle) from rest, shaft speed commanded: finite, moving forward."""
    graph = vc.graph_of(_remus())
    state, command = graph.step.sx_in(0), graph.step.sx_in(1)
    integrator = ca.integrator("rest", "cvodes", {"x": state, "p": command, "ode": graph.step(state, command)},
                               0, 0.1, {"abstol": 1e-8, "reltol": 1e-6})
    x, u = ca.DM.zeros(graph.num_states), ca.DM([0.0, 0.0, 1000.0])  # fin deflections 0, shaft speed 1000 rpm
    for _ in range(50):
        x = integrator(x0=x, p=u)["xf"]
    x = np.asarray(x).ravel()
    assert np.isfinite(x).all() and x[6] > 0.1  # u > 0.1 m/s after 5 s


# ---------------------------------------------------------------------------
# Derived quantities: computed by default, given by a named method
# ---------------------------------------------------------------------------

def _wrench(vessel, u):
    return np.asarray(vessel.signals(np.r_[REST[:6], 1.5, REST[7:]], u)["actuators.0.generated_force"])


def _lift_drag_force(vessel, x=None):
    x = np.r_[REST[:6], 1.5, 0.0, 0.3, REST[9:]] if x is None else x
    return np.asarray(vessel.signals(x, np.zeros(vc.input_size(vessel)))["hydrodynamics.lift_drag_force"])


def test_the_defaults_compute_density_planform_drag_and_fin_positions_from_the_geometry():
    """The cited REMUS 100 values are what the equations give on the REMUS geometry: body density
    ``m / (4/3 pi a b^2)`` (``spheroid.m`` 36), planform ``0.7 L D`` (``remus100.m`` 133), ``C_D0 = Cd pi (D/2)^2 / S``
    (``remus100.m`` 143-144), fins at ``-a`` (``remus100.m`` 184, 189): the vehicle built on the computed defaults
    equals the vehicle given those values."""
    p = vc.plugins()
    values = vc.remus_values()
    computed = _remus()
    given = vc.remus(full=True, current="none", vehicle_params={"body_density_method": "given"},
                     actuators=[(p.FinPairsDeflectionOnly, {**vc.pick(p.FinPairsDeflectionOnly, values),
                                                            "fin_positions_method": "given"}),
                                (p.Propeller, vc.pick(p.Propeller, values))])
    x = np.r_[1.0, 0.5, 2.0, 0.1, -0.05, 0.3, 1.5, 0.2, 0.3, 0.01, 0.02, 0.03]
    command = [0.1, -0.05, 800.0]
    a, b = (np.asarray(v.signals(x, command)["hydrodynamics.hydrodynamic_force"]) for v in (computed, given))
    assert np.abs(a - b).max() <= 1e-12
    assert np.abs(_wrench(computed, command) - _wrench(given, command)).max() <= 1e-12


def test_a_change_of_the_hull_moves_every_derived_value_and_leaves_the_given_mass():
    """The length and the diameter as a user sets them in the GUI: the mass stays the given 31.9 kg (the density
    follows), the planform and the parasitic drag follow the new geometry, the fin pairs move to ``-L/2``."""
    values = vc.remus_values()
    mass0 = vc.matrix(_remus().signals(REST, [0, 0, 0])["mass_matrix"])
    longer = {**values, "length": 2.4, "diameter": 0.22, "beam": 0.22}
    vessel = _remus(values=longer)
    assert float(vessel.signals(REST, [0, 0, 0])["mass"][0]) == pytest.approx(31.9, rel=1e-14)
    assert np.abs(vc.matrix(vessel.signals(REST, [0, 0, 0])["mass_matrix"]) - mass0).max() > 1.0  # the added mass moved

    # planform and parasitic drag: the same vehicle with the two values given from the equations equals it
    p = vc.plugins()
    planform = 0.7 * 2.4 * 0.22  # remus100.m 133
    parasitic = 0.42 * np.pi * (0.22 / 2) ** 2 / planform  # remus100.m 143-144
    given = {**vc.pick(p.AuvHullLoads, longer), "planform_area_method": "given", "planform_area": planform,
             "parasitic_drag_coefficient_method": "given", "parasitic_drag_coefficient": parasitic}
    reference = vc.build_context(
        p.SpheroidAuv, {**vc.pick(p.SpheroidAuv, vc.remus_values(longer)), "current_form": "none"},
        hydrostatics=(p.SubmergedRestoring, vc.pick(p.SubmergedRestoring, longer, ("center_of_buoyancy",))),
        hydrodynamics=(p.AuvHullLoads, given), actuators=[p.PrescribedWrench])
    computed = vc.remus(full=False, current="none", values=longer)
    assert np.abs(_lift_drag_force(computed) - _lift_drag_force(reference)).max() <= 1e-12
    # a vehicle that keeps the old given planform on the new hull differs: the value is not derived there
    stale = {**given, "planform_area": 0.7 * 1.6 * 0.19, "parasitic_drag_coefficient": 0.05595961914206819}
    old = vc.build_context(
        p.SpheroidAuv, {**vc.pick(p.SpheroidAuv, vc.remus_values(longer)), "current_form": "none"},
        hydrostatics=(p.SubmergedRestoring, vc.pick(p.SubmergedRestoring, longer, ("center_of_buoyancy",))),
        hydrodynamics=(p.AuvHullLoads, stale), actuators=[p.PrescribedWrench])
    assert np.abs(_lift_drag_force(computed) - _lift_drag_force(old)).max() > 1e-3

    # fin pairs at -L/2: the pitch moment of the stern plane over its heave force is the lever arm a = L/2
    wrench = _wrench(vessel, [0.0, 0.1, 800.0])
    assert wrench[4] / wrench[2] == pytest.approx(-(-2.4 / 2), rel=1e-12)  # M = -x_s Z_s (remus100.m 241-245)


def test_the_given_methods_read_their_values_and_the_computed_ones_ignore_them():
    """``given`` takes the typed value; under ``computed`` the typed value is not read."""
    p = vc.plugins()
    values = vc.remus_values()
    base = vc.pick(p.AuvHullLoads, values)
    x = np.r_[REST[:6], 1.5, 0.0, 0.3, REST[9:]]

    def lift(extra):
        v = vc.build_context(p.SpheroidAuv, {**vc.pick(p.SpheroidAuv, values), "current_form": "none"},
                             hydrostatics=(p.SubmergedRestoring, vc.pick(p.SubmergedRestoring, values,
                                                                          ("center_of_buoyancy",))),
                             hydrodynamics=(p.AuvHullLoads, {**base, **extra}), actuators=[p.PrescribedWrench])
        return _lift_drag_force(v, x)

    reference = lift({})
    assert np.array_equal(lift({"planform_area": 9.9, "parasitic_drag_coefficient": 9.9}), reference)  # ignored
    assert np.abs(lift({"planform_area_method": "given", "planform_area": 0.5}) - reference).max() > 1e-3
    assert np.abs(lift({"parasitic_drag_coefficient_method": "given", "parasitic_drag_coefficient": 0.5})
                  - reference).max() > 1e-3
    for name in ("planform_area_method", "parasitic_drag_coefficient_method"):
        with pytest.raises(ValueError, match=name):
            lift({name: "measured"})


def test_the_body_density_method_given_reads_the_density_and_computed_reads_the_mass():
    values = vc.remus_values()
    computed = _remus(vehicle_params=None, values={**values, "body_density": 1.0})  # the density is not read
    assert float(computed.signals(REST, [0, 0, 0])["mass"][0]) == pytest.approx(31.9, rel=1e-14)
    given = vc.remus(full=True, current="none", values={**values, "body_mass": 1.0},  # the mass is not read
                     vehicle_params={"body_density_method": "given"})
    assert float(given.signals(REST, [0, 0, 0])["mass"][0]) == pytest.approx(31.9, rel=1e-12)
    heavier = vc.remus(full=True, current="none", values={**values, "body_density": 1.1 * values["body_density"]},
                       vehicle_params={"body_density_method": "given"})
    assert float(heavier.signals(REST, [0, 0, 0])["mass"][0]) == pytest.approx(1.1 * 31.9, rel=1e-12)


def test_fin_positions_method_refuses_an_unknown_method():
    p = vc.plugins()
    with pytest.raises(ValueError, match="fin_positions_method"):
        _remus(actuators=[(p.FinPairsDeflectionOnly, {"fin_positions_method": "measured"})])


# ---------------------------------------------------------------------------
# The site
# ---------------------------------------------------------------------------

def test_site_form_given_takes_density_and_gravity_as_values():
    """MSS from 2.0.4 takes 1025 kg/m^3 and 9.81 m/s^2 (``mssConstants.m`` 12, 14): the weight is ``m g`` with that
    ``g`` and the water density reaches the signals."""
    vessel = _remus(vehicle_params={"site_form": "given", "gravity": 9.81}, values={**vc.remus_values(),
                                                                                    "water_density": 1025.0})
    signals = vessel.signals(REST, [0, 0, 0])
    assert float(signals["gravity"][0]) == 9.81 and float(signals["water_density"][0]) == 1025.0
    assert float(signals["weight"][0]) == pytest.approx(31.9 * 9.81, rel=1e-14)  # W = m g (Fossen 2011, eq. 4.1, p. 59)


def test_site_form_latitude_is_the_default_of_the_spheroid_and_gives_the_frozen_gravity():
    signals = _remus().signals(REST, [0, 0, 0])
    assert float(signals["gravity"][0]) == pytest.approx(vc.REMUS_FILE["gravity"]["g"], rel=1e-14)


def test_the_surface_types_offer_the_latitude_form_too():
    """Normal gravity at the equator, ``mu = 0``: ``9.7803253359`` m/s^2, the first constant of MSS ``gravity.m`` 11-12
    (``INS/functions``, ``cc07579``)."""
    p = vc.plugins()
    vessel = vc.build_context(p.Catamaran, {**vc.pick(p.Catamaran, vc.OTTER), "current_form": "none",
                                            "site_form": "latitude", "latitude": 0.0},
                              hydrostatics=(p.SurfaceRestoring, {"hull_count": 2}), hydrodynamics=p.SurfaceHullLoads)
    assert float(vessel.signals(REST, [])["gravity"][0]) == pytest.approx(9.7803253359, rel=1e-13)
    assert float(vc.otter().signals(REST, [])["gravity"][0]) == vc.OTTER["gravity"]  # the default stays "given"


# ---------------------------------------------------------------------------
# The wetted surface and the ITTC surge resistance
# ---------------------------------------------------------------------------

def _monohull_of(rho, length, beam, draft, block, surge="ittc"):
    """A monohull whose equilibrium draft is ``draft``: ``m = rho C_b L B T`` (``otter.m`` 121-122)."""
    p = vc.plugins()
    return vc.build_context(
        p.Monohull, {"water_density": rho, "length": length, "beam": beam, "hull_block_coefficient": block,
                     "hull_mass": rho * block * length * beam * draft, "current_form": "none"},
        hydrostatics=(p.SurfaceRestoring, {"hull_count": 1}),
        hydrodynamics=(p.SurfaceHullLoads, {"surge_resistance": surge}))


def _ittc_reference_rows():
    with open(DATA_HYDRODYNAMIC_LOADS / "xuu_ittc_mss_ac77394.csv", newline="") as handle:
        return list(csv.DictReader(handle))


def test_the_computed_wetted_surface_gives_the_ittc_reference_of_every_frozen_row():
    """The five hull sets of the frozen ``XuuITTC`` reference (every speed, the Reynolds floor included): ``Xuu =
    -1/2 rho S (1 + k) C_f`` with ``S`` the vehicle's own ``wetted_surface`` equals the frozen value
    (Fossen 2011, eqs. 6.82-6.85, p. 125; ``XuuITTC.m`` 35-39)."""
    rows = _ittc_reference_rows()
    sets = {}
    for row in rows:
        sets.setdefault(row["set"], []).append(row)
    assert len(sets) == 5
    for set_rows in sets.values():
        first = set_rows[0]
        rho, length, beam, draft, block = (float(first[k]) for k in ("rho", "L", "B", "T", "C_B"))
        vessel = _monohull_of(rho, length, beam, draft, block)
        signals = vessel.signals(REST, [])
        assert float(signals["draft"][0]) == pytest.approx(draft, rel=1e-12)
        surface = float(signals["wetted_surface"][0])
        assert surface == pytest.approx(1.025 * length * (block * beam + 1.7 * draft), rel=1e-13)  # XuuITTC.m 38
        for row in set_rows:
            reynolds = max(length * abs(float(row["u_r"])) / 1e-6, 1e5)  # XuuITTC.m 33-35
            friction = 0.075 / (np.log10(reynolds) - 2.0) ** 2  # (Fossen 2011, eq. 6.83, p. 125)
            xuu = -0.5 * rho * surface * (1 + 0.1) * friction
            assert xuu == pytest.approx(float(row["Xuu"]), rel=1e-12)


def test_the_ittc_surge_option_is_no_longer_refused_and_resists_the_motion():
    """The option builds on every surface craft; at high speed the force is the quadratic law ``Xuu |u| u`` (the
    linear part has faded: ``sigma = 1 - tanh(|u| / u_cross)`` is 6e-7 at 15 m/s, ``forceSurgeDamping.m`` 79 @ ac77394)."""
    rho, length, beam, draft, block = 1025.0, 83.0, 18.0, 5.0, 0.65
    vessel = _monohull_of(rho, length, beam, draft, block)
    speed = 15.0
    x = np.r_[REST[:6], speed, REST[7:]]
    signals = vessel.signals(x, [])
    surge = float(signals["hydrodynamics.surge_resistance_force"][0])
    surface = float(signals["wetted_surface"][0])
    friction = 0.075 / (np.log10(length * speed / 1e-6) - 2.0) ** 2
    assert surge == pytest.approx(-0.5 * rho * surface * 1.1 * friction * speed**2, rel=1e-5)
    p = vc.plugins()
    catamaran = vc.otter(hydrodynamics=(p.SurfaceHullLoads, {**vc.pick(p.SurfaceHullLoads, vc.OTTER),
                                                              "surge_resistance": "ittc"}))
    assert float(catamaran.signals(np.r_[REST[:6], 2.0, REST[7:]], [])["hydrodynamics.surge_resistance_force"][0]) < 0


def test_two_hulls_have_twice_the_wetted_surface_of_one_at_the_same_draft():
    """The Otter set: ``S = n 1.025 L (C_b B + 1.7 T)`` with ``n = 2`` and the pontoon beam, at the equilibrium
    draft of ``otter.m`` 121-122; the given default (1.77 m^2) is that value."""
    signals = vc.otter().signals(REST, [])
    draft = float(signals["draft"][0])
    assert draft == pytest.approx(80.0 / 1025.0 / (2 * 0.4 * 0.25 * 2.0), rel=1e-14)
    assert float(signals["wetted_surface"][0]) == pytest.approx(2 * 1.025 * 2.0 * (0.4 * 0.25 + 1.7 * draft), rel=1e-14)
    assert float(signals["wetted_surface"][0]) == pytest.approx(1.77, rel=1e-14)


def test_the_given_wetted_surface_is_read_and_the_method_refuses_an_unknown_one():
    p = vc.plugins()

    def surface(extra):
        v = vc.otter(hydrostatics=(p.SurfaceRestoring, {**vc.pick(p.SurfaceRestoring, vc.OTTER), "hull_count": 2,
                                                         **extra}))
        return float(v.signals(REST, [])["wetted_surface"][0])

    assert surface({"wetted_surface_method": "given", "wetted_surface": 3.0}) == 3.0
    assert surface({"wetted_surface": 3.0}) == pytest.approx(1.77, rel=1e-14)  # computed ignores the typed value
    with pytest.raises(ValueError, match="wetted_surface_method"):
        surface({"wetted_surface_method": "measured"})


# ---------------------------------------------------------------------------
# A positively buoyant vehicle with the centre of buoyancy off the centre of gravity
# ---------------------------------------------------------------------------

def _marie_like(**restoring):
    p = vc.plugins()
    values = vc.remus_values()
    return vc.build_context(
        p.SpheroidAuv, {**vc.pick(p.SpheroidAuv, values), "current_form": "none"},
        hydrostatics=(p.SubmergedRestoring, restoring), hydrodynamics=(p.AuvHullLoads, vc.pick(p.AuvHullLoads, values)),
        actuators=[p.PrescribedWrench])


def test_given_buoyancy_above_the_weight_gives_a_net_heave_force_and_a_trim_moment():
    """``B = 1.02 W`` given, centre of buoyancy 0.03 m ahead of the origin, centre of gravity at x = 0 (the
    spheroid's): at zero attitude the restoring force is ``[0, 0, W - B, 0, -(x_g W - x_b B), 0]`` (Fossen 2011,
    eqs. 4.5-4.6, p. 60): with z down positive the heave force is negative (the vehicle is pushed up) and the pitch
    moment is ``x_b B`` (a centre of buoyancy ahead of the centre of gravity)."""
    weight = 31.9 * vc.REMUS_FILE["gravity"]["g"]
    buoyancy, x_b = 1.02 * weight, 0.03
    vessel = _marie_like(buoyancy_method="given", buoyancy=buoyancy, center_of_buoyancy=[x_b, 0.0, 0.0])
    signals = vessel.signals(REST, [0, 0, 0, 0, 0, 0])
    force = np.asarray(signals["hydrostatics.restoring_force"])
    x_g = 0.0
    assert float(signals["weight"][0]) == pytest.approx(weight, rel=1e-14)
    assert float(signals["buoyancy"][0]) == pytest.approx(buoyancy, rel=1e-14)
    assert force[2] == pytest.approx(weight - buoyancy, rel=1e-12)  # heave: W - B < 0, upward
    assert force[4] == pytest.approx(-(x_g * weight - x_b * buoyancy), rel=1e-12)  # pitch moment x_b B
    assert force[[0, 1, 3, 5]] == pytest.approx(0.0, abs=1e-12)
    xdot = np.asarray(vc.graph_of(vessel).step(REST, [0, 0, 0, 0, 0, 0])).ravel()
    assert xdot[8] < 0  # accelerates upward (w negative, z down)


def test_buoyancy_from_the_displaced_volume_follows_the_same_equations():
    """``B = rho g nabla`` (Fossen 2011, eq. 4.1, p. 59) with 3 % more volume than the neutral hull."""
    values = vc.remus_values()
    gravity = vc.REMUS_FILE["gravity"]["g"]
    volume = 1.03 * 31.9 / values["water_density"]
    vessel = _marie_like(buoyancy_method="from_volume", displaced_volume=volume, center_of_buoyancy=[0.0, 0.0, -0.01])
    signals = vessel.signals(REST, [0, 0, 0, 0, 0, 0])
    buoyancy = values["water_density"] * gravity * volume
    assert float(signals["buoyancy"][0]) == pytest.approx(buoyancy, rel=1e-14)
    force = np.asarray(signals["hydrostatics.restoring_force"])
    assert force[2] == pytest.approx(31.9 * gravity - buoyancy, rel=1e-12)
    assert abs(force[4]) <= 1e-12  # both centres at x = 0: no trim moment at zero attitude


# ---------------------------------------------------------------------------
# The current form of the vehicle types
# ---------------------------------------------------------------------------

def test_the_default_current_form_is_the_full_attitude_one_and_adds_three_inputs():
    p = vc.plugins()
    for cls, hydrostatics, hydrodynamics in (
            (p.SpheroidAuv, p.SubmergedRestoring, p.AuvHullLoads),
            (p.Catamaran, p.SurfaceRestoring, p.SurfaceHullLoads),
            (p.Monohull, (p.SurfaceRestoring, {"hull_count": 1}), p.SurfaceHullLoads)):
        declared = {d.name: d.default_value for d in cls.PARAMETERS}
        assert declared["current_form"] == "full_attitude", cls.__name__
        vessel = vc.build_context(cls, {}, hydrostatics=hydrostatics, hydrodynamics=hydrodynamics)
        assert vc.input_names(vessel)[:3] == list(vc.CURRENT_NAMES), cls.__name__


def test_the_vehicle_in_a_current_equals_the_mss_forms_level_and_differs_pitched():
    """Level (phi = theta = 0) and without pitch and roll rate the full form gives the REMUS yaw-rate form; pitched
    and pitching it does not: the MSS shortcut drops the heave flow and the rotation-rate terms."""
    full = vc.remus(full=False, current="full_attitude")
    mss = vc.remus(full=False, current="yaw_rate_terms")
    command = [0.5, 0.6, 0.1] + [0.0] * 6  # current speed, set, vertical; prescribed wrench
    level = np.r_[1.0, 2.0, 3.0, 0.0, 0.0, 0.7, 1.5, 0.1, 0.0, 0.0, 0.0, 0.1]
    a = np.asarray(vc.graph_of(full).step(level, command)).ravel()
    b = np.asarray(vc.graph_of(mss).step(level, command)).ravel()
    assert np.abs(a - b).max() > 0  # the vertical current 0.1 m/s enters the full form; the MSS form has no heave flow
    command[2] = 0.0
    a = np.asarray(vc.graph_of(full).step(level, command)).ravel()
    b = np.asarray(vc.graph_of(mss).step(level, command)).ravel()
    assert np.abs(a - b).max() <= 1e-12
    pitched = np.r_[1.0, 2.0, 3.0, 0.0, 20 * DEG, 0.7, 1.5, 0.1, 0.0, 0.0, 0.1, 0.1]
    a = np.asarray(vc.graph_of(full).step(pitched, command)).ravel()
    b = np.asarray(vc.graph_of(mss).step(pitched, command)).ravel()
    assert np.abs(a[6:] - b[6:]).max() > 1e-3
