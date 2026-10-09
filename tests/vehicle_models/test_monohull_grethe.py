"""U8a: the test contract for Grethe's ``Monohull`` and independent references
for the options E-101 chose, written before U8b builds them.

Grethe is not an MSS/Fossen vehicle: there is no MATLAB reference to check
against (unlike REMUS 100 and the Otter). The references in
``data/grethe/`` are a second, independent implementation (numpy, no CasADi,
no ``more_dynamics`` import) of the cited equations — G5 against the published
anchors (Fossen 2011, Radojcic et al. 2014) and G2 against the library's own
formulas where they already match (rigid body, Coriolis, metacentric
restoring, Hoerner cross-flow, ITTC friction). ``Monohull``'s own declared
defaults are already Grethe's (the docstring of ``monohull.py``), the same
convention ``SpheroidAuv``/REMUS 100 and ``Catamaran``/Otter use; this file
reads them the way ``test_spheroid_auv.py`` and ``test_catamaran.py`` read
theirs.

**Sections, in order:** (1) G1 — today's physics against the independent
reference (green); (2) G4 — a perturbed input changes the result (green);
(3) refusals of nonsense on today's existing options (green); (4) the new
options E-101/A-64 Section 3.3 name — not yet built, so these fail today,
naming the missing option or parameter (U8b's work list); (5) sanity bands
and the validity flag (question 1); (6) the 2026 trim check (labelled,
skipped); (7) a named, test-only ``.rppws`` part "Grethe" built through rpp's
own ``ComponentContextBuilder``.

No library, model or plugin file is changed by this job (U8b builds them);
no August data (safeguard rule 20).

References
----------
See ``data/grethe/SOURCE.md`` and ``data/grethe/generate_grethe_references.py``
for every citation; this file states only which named test checks which row.

Author:    Enio Krizman
Date:      2026-10-09
"""

import csv
import json
import math
import shutil
import sys
import uuid
from pathlib import Path

import numpy as np
import pytest

import vehicle_contract as vc

DATA = Path(__file__).resolve().parent / "data" / "grethe"
sys.path.insert(0, str(DATA))  # generate_grethe_references.py (the independent reference, a test-only import)
G1_TOLERANCE = 1e-9
REST = np.zeros(12)

GRETHE_FILE = json.loads((DATA / "grethe_parameters.json").read_text())
GRETHE_RAW = {name: entry["value"] for name, entry in GRETHE_FILE["parameters"].items() if entry["value"] is not None}
# the two selector values E-101/A-64 chose that are not yet valid on today's library
# (U8b's job); "today" uses the valid equivalent so every other test builds cleanly.
GRETHE_TODAY = {**GRETHE_RAW, "surge_resistance": "ittc", "wetted_surface_method": "computed"}


def _read_csv(name):
    with open(DATA / name, newline="") as handle:
        return list(csv.DictReader(handle))


def _row_floats(row, prefix, count):
    return np.array([float(row[f"{prefix}_{i + 1:02d}"]) for i in range(count)])


def _grethe(*, vehicle_params=None, hydrostatics_params=None, hydrodynamics_params=None, values=None,
           actuators=()):
    """Grethe's ``Monohull``: today's valid parameter set (``values`` replaces it, e.g.
    with a perturbation or the not-yet-valid target set), still water, no actuator
    unless given (the same pattern ``vc.otter()`` uses for the Otter)."""
    p = vc.plugins()
    v = GRETHE_TODAY if values is None else values
    vehicle_p = {**vc.pick(p.Monohull, v), "current_form": "none", **(vehicle_params or {})}
    hydrostatics_p = {**vc.pick(p.SurfaceRestoring, v), **(hydrostatics_params or {})}
    hydrodynamics_p = {**vc.pick(p.SurfaceHullLoads, v), **(hydrodynamics_params or {})}
    return vc.build_context(p.Monohull, vehicle_p, hydrostatics=(p.SurfaceRestoring, hydrostatics_p),
                            hydrodynamics=(p.SurfaceHullLoads, hydrodynamics_p), actuators=actuators)


def _mass_matrix(vessel):
    return vc.matrix(vessel.signals(REST, [])["mass_matrix"])


# ---------------------------------------------------------------------------
# (1) G1 — today's physics against the independent reference
# ---------------------------------------------------------------------------

def test_G1_rigid_body_matches_hull_with_point_payload_at_both_payload_masses():
    """``M_RB`` through the rigid body's own ``mass``/``center_of_gravity``/``inertia``
    outputs, bare hull and the owner's +300 kg check case (E-64 row 7)."""
    rows = {float(r["payload_mass"]): r for r in _read_csv("grethe_rigid_body_added_mass.csv")
           if r["added_mass_mass_basis"] == "hull"}
    for payload_mass, row in rows.items():
        vessel = _grethe(vehicle_params={"payload_mass": payload_mass})
        signals = vessel.signals(REST, [])
        assert float(signals["mass"][0]) == pytest.approx(float(row["mass"]), rel=G1_TOLERANCE)
        cg = np.asarray(signals["center_of_gravity"]).ravel()
        assert cg[0] == pytest.approx(float(row["x_g"]), abs=1e-12)
        assert cg[2] == pytest.approx(float(row["z_g"]), rel=G1_TOLERANCE)
        inertia = np.asarray(signals["inertia"]).reshape(3, 3, order="F")
        assert inertia[0, 0] == pytest.approx(float(row["I11"]), rel=G1_TOLERANCE)
        assert inertia[1, 1] == pytest.approx(float(row["I22"]), rel=G1_TOLERANCE)
        assert inertia[2, 2] == pytest.approx(float(row["I33"]), rel=G1_TOLERANCE)


def test_G1_added_mass_today_scales_with_the_hull_mass_not_the_displaced_mass():
    """A-64 Section 3.1 gap 6: today's ``scaled_derivatives`` takes ``hull_mass``, so
    ``A11`` and the two ``m``-scaled entries (surge, sway, heave — indices 0,1,2 of
    ``[A11, m, m, I11, I22, I33]``) are the **same** at payload 0 and +300 kg — this
    is the deviation ``added_mass_mass_basis = "displaced"`` (E-101 Q5 a) is meant to
    fix. The rotational entries still move a little: the payload sits off the hull's
    own CG (``payload_position`` [0,0,0] vs ``hull_center_of_gravity`` [0,0,0.025]),
    which shifts the *rigid-body* inertia the added mass scales, not the mass basis."""
    # test fix (U8b): isolates today's "hull" basis explicitly -- GRETHE_TODAY itself now carries the owner's
    # "displaced" decision (E-101 Q5 a) for every other test, once added_mass_mass_basis exists as a parameter.
    rows = {float(r["payload_mass"]): r for r in _read_csv("grethe_rigid_body_added_mass.csv")
           if r["added_mass_mass_basis"] == "hull"}
    bare = vc.matrix(_grethe(vehicle_params={"payload_mass": 0.0, "added_mass_mass_basis": "hull"})
                     .signals(REST, [])["added_mass_matrix"])
    loaded = vc.matrix(_grethe(vehicle_params={"payload_mass": 300.0, "added_mass_mass_basis": "hull"})
                       .signals(REST, [])["added_mass_matrix"])
    assert np.array_equal(bare[:3, :3], loaded[:3, :3])  # surge/sway/heave unchanged: today's deviation, exact
    assert bare[5, 5] == loaded[5, 5]  # yaw: a pure-z CG shift leaves I33 unchanged (parallel-axis)
    for i in range(6):
        assert bare[i, i] == pytest.approx(float(rows[0.0][f"M_A_{i + 1}{i + 1}"]), rel=G1_TOLERANCE)
        assert loaded[i, i] == pytest.approx(float(rows[300.0][f"M_A_{i + 1}{i + 1}"]), rel=G1_TOLERANCE)


def test_G1_coriolis_matches_kirchhoff_full_and_the_munk_couplings_removed_form():
    rows = _read_csv("grethe_coriolis.csv")
    vessel = _grethe(vehicle_params={"coriolis_form": "kirchhoff_full"})
    vessel_munk = _grethe(vehicle_params={"coriolis_form": "munk_couplings_removed"})
    for row in rows:
        nu = _row_floats(row, "nu", 6)
        state = np.r_[np.zeros(6), nu]
        c_rb = vc.matrix(vessel.signals(state, [])["rigid_body_coriolis_matrix"])
        c_a_full = vc.matrix(vessel.signals(state, [])["added_mass_coriolis_matrix"])
        c_a_munk = vc.matrix(vessel_munk.signals(state, [])["added_mass_coriolis_matrix"])
        assert np.abs(c_rb.ravel(order="F") - _row_floats(row, "C_RB", 36)).max() <= G1_TOLERANCE
        assert np.abs(c_a_full.ravel(order="F") - _row_floats(row, "C_A_full", 36)).max() <= G1_TOLERANCE
        assert np.abs(c_a_munk.ravel(order="F") - _row_floats(row, "C_A_munk_removed", 36)).max() <= G1_TOLERANCE


def test_G1_restoring_matches_the_metacentric_form_at_the_waterline_midpoint_co():
    rows = {float(r["payload_mass"]): r for r in _read_csv("grethe_restoring.csv")}
    for payload_mass, row in rows.items():
        signals = _grethe(vehicle_params={"payload_mass": payload_mass}).signals(REST, [])
        assert float(signals["draft"][0]) == pytest.approx(float(row["draft"]), rel=G1_TOLERANCE)
        assert float(signals["wetted_surface"][0]) == pytest.approx(float(row["wetted_surface_mumford"]), rel=G1_TOLERANCE)
        restoring = vc.matrix(signals["restoring_matrix"])
        assert restoring[2, 2] == pytest.approx(float(row["G_33"]), rel=G1_TOLERANCE)
        assert restoring[3, 3] == pytest.approx(float(row["G_44"]), rel=G1_TOLERANCE)
        assert restoring[4, 4] == pytest.approx(float(row["G_55"]), rel=G1_TOLERANCE)


def test_G1_cross_flow_matches_the_hoerner_strip_integral_on_todays_hull_basis():
    """Today's ``single_hull()`` feeds the strip the overall length/beam (5.2/2.15 m),
    not the waterline dimensions (``tau_overall_*``; ``tau_waterline_*`` is the future
    basis, checked separately below as a red test)."""
    rows = _read_csv("grethe_cross_flow.csv")
    vessel = _grethe()
    for row in rows:
        nu_r = _row_floats(row, "nu_r", 6)
        state = np.r_[np.zeros(6), nu_r]
        force = np.asarray(vessel.signals(state, [])["hydrodynamics.cross_flow_force"]).ravel()
        assert np.abs(force - _row_floats(row, "tau_overall", 6)).max() <= G1_TOLERANCE


def test_G1_ittc_only_matches_the_published_friction_line_on_todays_overall_length():
    """``surge_resistance = "ittc"`` already exists (R9); its reference here is the
    ``L52`` (overall length/beam, today's basis) column of ``grethe_surge_resistance.csv``,
    read in isolation through the named output (the calibrated term is zeroed, as the
    option's own docstring instructs, A-64 Section 3.1 gap 2)."""
    rows = _read_csv("grethe_surge_resistance.csv")
    vessel = _grethe(hydrodynamics_params={"max_forward_thrust": 0.0})
    for row in rows:
        u = float(row["u"])
        state = np.r_[np.zeros(6), u, np.zeros(5)]
        force = float(vessel.signals(state, [])["hydrodynamics.surge_resistance_force"][0])
        assert force == pytest.approx(float(row["X_ittc_only_L52"]), rel=1e-6, abs=1e-9)


def test_G1_linear_coupled_prior_values_match_the_time_constant_formula():
    """``Y_v = -M22/T_sway``, ``N_r = -M66/T_yaw`` on Grethe's own mass matrix, the
    formula the library's ``time_constant_damping_surface`` already uses (otter.m
    203, 207) — the prior ``linear_coupled`` (not yet built) starts from."""
    row = _read_csv("grethe_damping.csv")[0]
    mass = _mass_matrix(_grethe())
    y_v = -mass[1, 1] / GRETHE_TODAY["time_constants"][0]
    n_r = -mass[5, 5] / GRETHE_TODAY["time_constants"][1]
    assert y_v == pytest.approx(float(row["Y_v"]), rel=G1_TOLERANCE)
    assert n_r == pytest.approx(float(row["N_r"]), rel=G1_TOLERANCE)
    assert float(row["Y_r"]) == 0.0 and float(row["N_v"]) == 0.0  # E-101 Q2 a, unidentified


# ---------------------------------------------------------------------------
# (2) G4 — a perturbed input moves the result
# ---------------------------------------------------------------------------

def test_G4_a_heavier_hull_moves_the_draft_and_the_metacentric_height():
    base = _grethe().signals(REST, [])
    heavier = _grethe(values={**GRETHE_TODAY, "hull_mass": 1106.0}).signals(REST, [])  # +300 kg (E-64 row 7)
    assert float(heavier["draft"][0]) > float(base["draft"][0])
    base_g, heavier_g = vc.matrix(base["restoring_matrix"]), vc.matrix(heavier["restoring_matrix"])
    assert abs(heavier_g[4, 4] - base_g[4, 4]) > 1.0  # GM_L moves (N m/rad)


def test_G4_a_wider_waterline_beam_moves_the_cross_flow_force():
    """Perturbing the future ``waterline_beam`` input changes ``tau_waterline`` in the
    frozen reference (confirms the gate is sensitive, independent of the library)."""
    from generate_grethe_references import cross_flow_hoerner
    row = _read_csv("grethe_cross_flow.csv")[0]
    nu_r = _row_floats(row, "nu_r", 6)
    draft = float(_read_csv("grethe_restoring.csv")[0]["draft"])
    reference = _row_floats(row, "tau_waterline", 6)
    wider = cross_flow_hoerner(nu_r, GRETHE_RAW["waterline_length"], GRETHE_RAW["waterline_beam"] * 1.3, draft, 1025.0)
    assert np.abs(wider - reference).max() > 1.0


# ---------------------------------------------------------------------------
# (3) Refusals of nonsense on today's existing options
# ---------------------------------------------------------------------------

def test_an_unknown_surge_resistance_is_refused_by_name():
    with pytest.raises(ValueError, match="surge_resistance"):
        _grethe(hydrodynamics_params={"surge_resistance": "bogus"})


def test_an_unknown_hull_count_is_refused_by_name():
    with pytest.raises(ValueError, match="hull_count"):
        _grethe(hydrostatics_params={"hull_count": 3})


def test_an_unknown_wetted_surface_method_is_refused_by_name():
    with pytest.raises(ValueError, match="wetted_surface_method"):
        _grethe(hydrostatics_params={"wetted_surface_method": "measured"})


# ---------------------------------------------------------------------------
# (4) The new options of E-101 / A-64 Section 3.3 — not yet built (U8b's work list)
# ---------------------------------------------------------------------------

def test_surge_resistance_ittc_residual_is_grethes_default_but_does_not_exist_yet():
    """E-101 Q1 a: ``SurfaceHullLoads.surge_resistance = "ittc_residual"``. Today
    ``SURGE_RESISTANCES = ("none", "ittc")``: this fails, naming the option, until U8b.
    Every other parameter stays at today's valid values, so the failure is isolated
    to this one option."""
    vessel = _grethe(hydrodynamics_params={"surge_resistance": "ittc_residual", "max_forward_thrust": 0.0})
    rows = _read_csv("grethe_surge_resistance.csv")
    for row in rows:
        u = float(row["u"])
        state = np.r_[np.zeros(6), u, np.zeros(5)]
        force = float(vessel.signals(state, [])["hydrodynamics.surge_resistance_force"][0])
        assert force == pytest.approx(float(row["X_ittc_residual_L52"]), rel=1e-6, abs=1e-9)


def test_wetted_surface_method_regression_table_does_not_exist_yet():
    """A-64 Section 3.3 R11: the RA14 wetted-surface option. Today
    ``WETTED_SURFACE_METHODS = ("computed", "given")``: this fails, naming the method."""
    vessel = _grethe(hydrostatics_params={"wetted_surface_method": "regression_table"})
    signals = vessel.signals(REST, [])
    row = _read_csv("grethe_surge_resistance.csv")[0]
    assert float(signals["wetted_surface"][0]) == pytest.approx(float(row["S_RA14_L52"]), rel=1e-6)


def test_added_mass_mass_basis_option_does_not_exist_yet():
    """E-101 Q5 a: added mass should follow the displaced mass (hull + payload), not
    the hull mass alone. Today ``Monohull`` declares no ``added_mass_mass_basis``
    parameter at all: ``vc.build_context`` refuses it, naming it."""
    rows = {(float(r["payload_mass"]), r["added_mass_mass_basis"]): r
           for r in _read_csv("grethe_rigid_body_added_mass.csv")}
    vessel = _grethe(vehicle_params={"payload_mass": 300.0, "added_mass_mass_basis": "displaced"})
    added_mass = vc.matrix(vessel.signals(REST, [])["added_mass_matrix"])
    expected = rows[(300.0, "displaced")]
    for i in range(6):
        assert added_mass[i, i] == pytest.approx(float(expected[f"M_A_{i + 1}{i + 1}"]), rel=G1_TOLERANCE)


def test_manoeuvring_damping_option_does_not_exist_yet():
    """A-64 Section 3.3: ``SurfaceHullLoads.manoeuvring_damping = "linear_coupled"``
    with ``sway_yaw_derivatives``. Neither parameter is declared today."""
    row = _read_csv("grethe_damping.csv")[0]
    sway_yaw_derivatives = [float(row["Y_v"]), float(row["Y_r"]), float(row["N_v"]), float(row["N_r"])]
    vessel = _grethe(hydrodynamics_params={"manoeuvring_damping": "linear_coupled",
                                           "sway_yaw_derivatives": sway_yaw_derivatives})
    nu_r = np.array([0.0, float(row["v_r"]), 0.0, 0.0, 0.0, float(row["r_r"])])
    tau = np.asarray(vessel.signals(np.r_[np.zeros(6), nu_r], [])["hydrodynamics.damping_force"]).ravel()
    assert tau[1] == pytest.approx(float(row["tau_Y"]), rel=G1_TOLERANCE)
    assert tau[5] == pytest.approx(float(row["tau_N"]), rel=G1_TOLERANCE)


def test_cross_flow_switch_does_not_exist_yet():
    """A-64 Section 3.3: a ``cross_flow`` selector (``"hoerner_strips"`` / ``"none"``,
    the second once the modulus terms are identified, A-44 Part 4). Not declared today."""
    _grethe(hydrodynamics_params={"cross_flow": "none"})


def test_waterline_length_and_beam_are_not_separate_monohull_parameters_yet():
    """A-64 Section 3.3: ``waterline_length``/``waterline_beam``, read by restoring,
    cross-flow and resistance instead of ``length``/``beam`` when set. ``Monohull``
    declares neither today: once it does, the vehicle built on them must equal the
    ``tau_waterline_*`` reference above (test_G4 confirms the reference itself moves)."""
    _grethe(vehicle_params={"waterline_length": GRETHE_RAW["waterline_length"],
                            "waterline_beam": GRETHE_RAW["waterline_beam"]})


# ---------------------------------------------------------------------------
# (5) Sanity bands and the validity flag (question 1)
# ---------------------------------------------------------------------------

def test_validity_is_froude_number_on_the_displaced_volume_not_on_the_overall_length():
    """A-64 Section 1.2: at the logged top speed (2.91 m/s) Fn_V is just under 1 (the
    RA14 regime boundary, Section Summary); the "Fn <= 0.35" line written after E-57
    was based on the 5.2 m overall length and is not this validity flag (question 1)."""
    rows = _read_csv("grethe_surge_resistance.csv")
    at_291 = min(rows, key=lambda r: abs(float(r["u"]) - 2.91))
    assert float(at_291["Fn_V_L52"]) < 1.0
    beyond = [r for r in rows if float(r["u"]) > 3.09]
    assert any(float(r["Fn_V_L52"]) > 1.0 for r in beyond)


def _coast_time(vessel, v0, v_target, mass, dt=0.02, t_max=300.0):
    """Explicit-Euler surge-only coast-down: ``m dv/dt = X(v)`` (the vehicle's total
    surge hydrodynamic force at ``nu_r = [v, 0, 0, 0, 0, 0]``), from ``v0`` down to
    ``v_target``. A sanity-band tool, not a frozen reference (no claim on Grethe's
    real coasting distance — that needs the outboard, A-64 Section 5.1, U12)."""
    v, t = v0, 0.0
    while v > v_target and t < t_max:
        state = np.r_[np.zeros(6), v, np.zeros(5)]
        force = float(vessel.signals(state, [])["hydrodynamics.hydrodynamic_force"][0])
        v = v + (force / mass) * dt
        t += dt
    return t


def test_G3_today_the_calibrated_otter_default_coasts_down_much_faster_than_ittc():
    """Sanity band for A-64 Section 6 finding 1 ("the calibrated linear surge law
    loses about four times the logged speed in the first 10 s when coasting", and
    coast-down distance 11 m vs 40-133 m, V1 vs V3-V5): the Otter-default calibrated
    ``Xu`` (``max_forward_thrust``/``max_speed`` are the Otter's, not Grethe's — A-64
    Section 3.1 gap 3/4) must coast down far faster than today's ITTC physics form.
    This only needs the hull, no outboard (A-64's own setup for this comparison)."""
    mass = float(_grethe().signals(REST, [])["mass"][0])
    calibrated = _grethe(hydrodynamics_params={"surge_resistance": "none"})  # today's default: the Otter-default Xu
    ittc = _grethe(hydrodynamics_params={"max_forward_thrust": 0.0})  # isolates the ITTC term (A-64 Section 3.1 gap 2)
    t_calibrated = _coast_time(calibrated, 2.5, 0.3, mass)
    t_ittc = _coast_time(ittc, 2.5, 0.3, mass)
    assert t_ittc > 2.0 * t_calibrated  # A-64 Section 6: about a factor of 4, this just bounds the direction and order


def test_G3_todays_sway_yaw_damping_is_straight_line_unstable_without_the_coupling_terms():
    """Sanity band for A-64 Section 6 finding 3 ("straight-line unstable at every band
    with these priors"): the Jacobian of the real, composed vehicle's sway-yaw
    acceleration at a representative forward speed (today's diagonal
    ``time_constants`` damping, no ``Y_r``/``N_v`` coupling, plus the added-mass Munk
    moment in the Coriolis term) has an eigenvalue with positive real part — the
    reason ``manoeuvring_damping = "linear_coupled"`` (identified ``Y_r``, ``N_v``) is
    needed before closed-loop design (question 2)."""
    u0 = 1.5  # m/s, a representative forward speed (A-45 band)
    graph = vc.graph_of(_grethe())
    step = graph.step

    def rate(v, r):
        state = REST.copy()
        state[6], state[7], state[11] = u0, v, r
        return np.asarray(step(state, [])).ravel()[[7, 11]]  # [v_dot, r_dot]

    h = 1e-4
    jacobian = np.column_stack([(rate(h, 0.0) - rate(-h, 0.0)) / (2 * h),
                               (rate(0.0, h) - rate(0.0, -h)) / (2 * h)])
    eigenvalues = np.linalg.eigvals(jacobian)
    assert np.any(eigenvalues.real > 0.0)


@pytest.mark.skip(reason="A-64 Section 1.1: the 2026 static trim (-4.11 deg, -0.55 deg) is an owner design "
                        "argument, not a measurement, and is used only as a labelled equilibrium check, not a "
                        "gate, until X-9 item G14 (inclinometer + deck reference at the quay) and G2 (drafts "
                        "fore and aft) confirm or refute it (E-101 addendum).")
def test_trim_check_2026_static_pitch_and_roll_from_the_dune_calibration():
    pass


# ---------------------------------------------------------------------------
# (7) A named, test-only `.rppws` part "Grethe"
# ---------------------------------------------------------------------------

def _write_description(part_dir, *, part_id, name, plugin_type, plugin_name, library, spec=None, subcomponents=None):
    (part_dir / "params").mkdir(parents=True, exist_ok=True)
    (part_dir / "description.json").write_text(json.dumps({
        "Id": part_id, "Name": name, "PluginType": plugin_type, "PluginName": plugin_name, "Library": library,
        "SubcomponentSpec": spec or {},
        "Subcomponents": {k: v["entry"] for k, v in (subcomponents or {}).items()},
        "ParentComponentInfo": None,
    }, indent=4))
    (part_dir / "callbacks.py").write_text("from __future__ import annotations\n")


def _write_params(part_dir, values):
    body = "".join(f"    {k} = {v!r}\n" for k, v in values.items()) or "    pass\n"
    (part_dir / "params" / "parameters.py").write_text(
        f"from __future__ import annotations\n\n\nclass ComponentParameters:\n{body}")


def _named_grethe_part(tmp_path, library_root, *, vehicle_params, hydrostatics_params, hydrodynamics_params):
    """A ``.rppws`` tree that exists only for this test (never written beside the
    library's own ``.rppws``): one ``Monohull`` named "Grethe", parametrised the way
    rpp-orchestrator's GUI would save it (E-99: ``.rppws`` vehicles are created and
    parametrised in the GUI, no generator script) — hydrostatics/hydrodynamics
    nested under the vehicle part's own ``subcomponents/``, the shape the real
    ``more_dynamics__catamaran`` Otter part uses."""
    (tmp_path / "scripts" / "vehicles").mkdir(parents=True)
    shutil.copy(library_root / "scripts" / "vehicles" / "vehicle_simulation.py",
               tmp_path / "scripts" / "vehicles" / "vehicle_simulation.py")
    vehicle_id, hydrostatics_id, hydrodynamics_id = (str(uuid.uuid4()) for _ in range(3))
    vehicle_dir = tmp_path / ".rppws" / "parts" / "more_dynamics__monohull" / vehicle_id
    hydrostatics_dir = vehicle_dir / "subcomponents" / hydrostatics_id
    hydrodynamics_dir = vehicle_dir / "subcomponents" / hydrodynamics_id

    _write_description(hydrostatics_dir, part_id=hydrostatics_id, name="Grethe hydrostatics",
                       plugin_type="more_dynamics::HydrostaticsModel", plugin_name="more_dynamics::SurfaceRestoring",
                       library="more_dynamics")
    _write_params(hydrostatics_dir, hydrostatics_params)
    _write_description(hydrodynamics_dir, part_id=hydrodynamics_id, name="Grethe hydrodynamics",
                       plugin_type="more_dynamics::HydrodynamicsModel", plugin_name="more_dynamics::SurfaceHullLoads",
                       library="more_dynamics")
    _write_params(hydrodynamics_dir, hydrodynamics_params)
    _write_description(vehicle_dir, part_id=vehicle_id, name="Grethe", plugin_type="more_dynamics::VehicleModel3D",
                       plugin_name="more_dynamics::Monohull", library="more_dynamics",
                       spec={"actuators": "List[more_dynamics::ForceProducer]",
                            "hydrostatics": "more_dynamics::HydrostaticsModel",
                            "hydrodynamics": "more_dynamics::HydrodynamicsModel"}, subcomponents={
        "hydrostatics": {"entry": {"Id": hydrostatics_id, "PluginType": "more_dynamics::HydrostaticsModel",
                                   "PluginName": "more_dynamics::SurfaceRestoring", "SlotName": "hydrostatics",
                                   "Library": "more_dynamics", "IsLinked": False}},
        "hydrodynamics": {"entry": {"Id": hydrodynamics_id, "PluginType": "more_dynamics::HydrodynamicsModel",
                                    "PluginName": "more_dynamics::SurfaceHullLoads", "SlotName": "hydrodynamics",
                                    "Library": "more_dynamics", "IsLinked": False}},
    })
    _write_params(vehicle_dir, vehicle_params)

    description = {
        "ScriptPath": "scripts/vehicles/vehicle_simulation.py", "Language": "python",
        "Configurations": {"grethe": {"Description": "Grethe, NTNU Mariner 5 USV (test-only part, U8a)",
                                     "Components": {"vessels": [{"Id": vehicle_id, "PluginName": "more_dynamics::Monohull"}]}}},
        "ActiveConfiguration": "grethe", "Spec": {"vessels": "List[more_dynamics::VehicleModel3D]"},
    }
    description_path = tmp_path / ".rppws" / "script_descriptions" / "grethe_vehicle.json"
    description_path.parent.mkdir(parents=True, exist_ok=True)
    description_path.write_text(json.dumps(description, indent=4))
    return description_path


def _build_named_part(tmp_path, description_path):
    try:
        from rpp_py.context_builder import ComponentContextBuilder
        from rpp_py.data_manager import DataManager
        from rpp_plugin_registrator.library_manager import LibraryManager
    except Exception as exc:  # pragma: no cover - environment without rpp
        pytest.skip(f"{vc.REGISTER_MESSAGE} ({type(exc).__name__}: {str(exc)[:80]})")
    vc.plugins()
    data_manager = DataManager(library_manager=LibraryManager(), workspace_path=str(tmp_path))
    context = ComponentContextBuilder(data_manager=data_manager).build_script_from_description_path(
        str(description_path), configuration="grethe")
    context.initialize()
    return context.get_component("vessels")[0]


def test_named_part_grethe_builds_with_todays_valid_options(tmp_path):
    library_root = Path(__file__).resolve().parents[2]
    description_path = _named_grethe_part(
        tmp_path, library_root, vehicle_params={"current_form": "none"},
        hydrostatics_params={"hull_count": 1, "longitudinal_inertia_factor": GRETHE_TODAY["longitudinal_inertia_factor"],
                             "longitudinal_center_of_flotation": GRETHE_TODAY["longitudinal_center_of_flotation"],
                             "reference_point": GRETHE_TODAY["reference_point"]},
        hydrodynamics_params={"surge_resistance": "ittc", "time_constants": GRETHE_TODAY["time_constants"],
                              "damping_ratios": GRETHE_TODAY["damping_ratios"]})
    built = _build_named_part(tmp_path, description_path)
    assert type(built).__name__ == "Monohull"
    reference = _grethe()
    assert np.array_equal(np.asarray(vc.graph_of(built).step(REST, [])),
                         np.asarray(vc.graph_of(reference).step(REST, [])))


def test_named_part_grethe_with_the_full_g0_target_set_does_not_build_yet(tmp_path):
    """The option values E-101/A-64 Section 3.4 name as Grethe's "G0 documents" set
    (``surge_resistance = "ittc_residual"``, ``wetted_surface_method =
    "regression_table"``): the test-only part is written the same way, and fails to
    initialize until U8b lands those options."""
    library_root = Path(__file__).resolve().parents[2]
    description_path = _named_grethe_part(
        tmp_path, library_root, vehicle_params={"current_form": "none"},
        hydrostatics_params={"hull_count": 1, "wetted_surface_method": "regression_table",
                             "reference_point": GRETHE_TODAY["reference_point"]},
        hydrodynamics_params={"surge_resistance": "ittc_residual", "time_constants": GRETHE_TODAY["time_constants"]})
    _build_named_part(tmp_path, description_path)
