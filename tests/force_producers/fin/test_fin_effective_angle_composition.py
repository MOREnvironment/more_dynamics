"""Gate tests of the effective-angle fin composition: ``LiftingFin`` with
servo ``ideal``, inflow ``rigid_point``, flow angle ``small_angle``,
interference ``none`` and section ``linear_section`` (lift only) or
``quadratic_drag``.

Written 2026-10-08, before the parts exist (contract: ``fin_parts_contract``).

Gates:

* G5 (form): a rudder lifting along +y and a stern plane lifting along -z,
  both at ``(x_fin, 0, 0)``, with ``linear_section`` and no fin drag, equal
  Prestero 2001 eq. 4.43 (p. 33) in any motion:
  ``Y_r = 1/2 rho c_L_alpha S_fin [u^2 delta_r - u v - x_fin u r]``,
  ``Z_s = -1/2 rho c_L_alpha S_fin [u^2 delta_s + u w - x_fin u q]``,
  ``M_s = -x_fin Z_s``, ``N_r = x_fin Y_r``, ``X = K = 0``; numbers from his
  Tables A.1 and A.5 (pp. 102-103).
* G2: off the centre line and with any perpendicular chord and lift axes,
  the fin equals the chain eq. 4.40 (all terms) -> eq. 4.42 -> eq. 4.41 ->
  eq. 4.37 -> ``tau = [f; r x f]`` (Fossen 2011, eq. 12.226, p. 400) written
  out here.
* G5 (table anchor): four fins (a rudder pair, a stern-plane pair, one
  command per pair) give the twelve linear coefficients of Table C.10
  (p. 111) in the ratios printed there. The table's scale (``Y_uudr =
  9.64``) does **not** follow from eq. 4.44 with Tables A.1 and A.5
  (``rho c_L_alpha S_fin = 1030 * 3.12 * 0.00665 = 21.4``, with the pair's
  two halves); the composition follows the equation, and the printed scale
  is pinned below as an inconsistency of the source, so the test that
  anchors to the table uses the ratios only.
* The gradient of ``tau`` with respect to an open part parameter
  (``section.lift_slope``) in one call equals a central difference.

References
----------
[Fossen 2011] Fossen, T. I. (2011). *Handbook of Marine Craft Hydrodynamics
    and Motion Control*, 1st ed. John Wiley & Sons, Chichester. Ch. 12,
    eq. 12.226, p. 400.
[Prestero 2001] Prestero, T. (2001). *Verification of a six-degree of
    freedom simulation model for the REMUS autonomous underwater vehicle*.
    MSc thesis, MIT/WHOI. Ch. 4, eqs. 4.37, 4.40-4.45, pp. 31-33; App. A,
    Tables A.1, A.5, pp. 102-103; App. C, Table C.10, p. 111.

Author:    Enio Krizman
Date:      2026-10-08
"""

import casadi as ca
import numpy as np
import pytest

from fin_parts_contract import (
    E_X,
    E_Y,
    G2_TOLERANCE,
    N_RANDOM,
    NEG_Y,
    NEG_Z,
    call,
    fin_from_forms,
    geometry,
    max_diff,
    module,
    part,
    prestero,
    prestero_value,
    producer_set,
    random_nu_r,
    rng,
    skeleton,
    tau,
)

EFFECTIVE_ANGLE = {"servo": "ideal", "inflow": "rigid_point", "flow_angle": "small_angle",
                   "interference": "none", "section": "linear_section"}
MAX_DEFLECTION = 0.4   # rad, a test construction; commands stay inside it


def prestero_fin(position, lift_axis, chord_axis=E_X, section="linear_section"):
    forms = {**EFFECTIVE_ANGLE, "section": section}
    values = {"lift_slope": prestero_value("lift_slope")}
    if section == "linear_section":
        values["zero_lift_drag"] = 0.0
    return fin_from_forms(forms, {"servo": {"max_deflection": MAX_DEFLECTION}, "section": values},
                          geometry(position, lift_axis, prestero_value("fin_area"), chord_axis))


def _cases(n, generator):
    nu_r = random_nu_r(generator, n)
    nu_r[:, 0] = np.where(np.abs(nu_r[:, 0]) < 0.05, 0.05, nu_r[:, 0])   # away from u = 0 (small angles)
    return generator.uniform(-0.3, 0.3, (n, 2)), nu_r


def test_G5_rudder_and_stern_plane_equal_prestero_eq_4_43():
    x_fin, rho = prestero_value("fin_position_x"), prestero_value("water_density")
    k = 0.5 * rho * prestero_value("lift_slope") * prestero_value("fin_area")
    fin_set = producer_set([prestero_fin([x_fin, 0, 0], E_Y), prestero_fin([x_fin, 0, 0], NEG_Z)], np.eye(2))
    commands, nus = _cases(N_RANDOM, rng())
    for (delta_r, delta_s), nu_r in zip(commands, nus):
        u, v, w, _, q, r = nu_r
        lateral = u ** 2 * delta_r - u * v - x_fin * u * r          # (Prestero 2001, eq. 4.43, p. 33)
        vertical = u ** 2 * delta_s + u * w - x_fin * u * q
        expected = [0.0, k * lateral, -k * vertical, 0.0, k * x_fin * vertical, k * x_fin * lateral]
        assert max_diff(tau(fin_set, [delta_r, delta_s], nu_r, rho), expected) <= G2_TOLERANCE, nu_r


def test_G2_fin_anywhere_equals_the_chain_of_prestero_equations():
    g = rng()
    rho, area, slope = 1025.0, 0.013, 2.9
    for _ in range(100):
        q, _ = np.linalg.qr(g.normal(size=(3, 3)))
        chord, lift, r = q[:, 0], q[:, 1], g.uniform(-1, 1, 3)
        one = fin_from_forms(EFFECTIVE_ANGLE, {"servo": {"max_deflection": MAX_DEFLECTION},
                                               "section": {"lift_slope": slope, "zero_lift_drag": 0.0}},
                             geometry(r, lift, area, chord))
        delta, nu_r = g.uniform(-0.3, 0.3), random_nu_r(g, 1)[0]
        v_fin = nu_r[0:3] + np.cross(nu_r[3:6], r)                     # (eq. 4.40, all terms)
        u_c, u_n = v_fin @ chord, v_fin @ lift
        if abs(u_c) < 0.05:
            continue
        alpha = delta - u_n / u_c                                     # (eqs. 4.41-4.42)
        f = 0.5 * rho * area * u_c ** 2 * slope * alpha * lift        # (eq. 4.37; lift on u_c^2, eq. 4.43)
        expected = np.concatenate([f, np.cross(r, f)])                # (Fossen 2011, eq. 12.226)
        out = call(one, command=delta, state=np.zeros(0), nu_r=nu_r, water_density=rho)
        assert max_diff(out["tau"], expected) <= G2_TOLERANCE
        assert abs(out["angle_of_attack"][0] - alpha) <= G2_TOLERANCE * max(1.0, abs(alpha))
        assert abs(out["deflection"][0] - delta) <= G2_TOLERANCE


# --------------------------------------------------------------------------
# Table C.10 (table anchor)
# --------------------------------------------------------------------------
TABLE_ENTRIES = {   # name -> (output row, "delta" or "nu", column)
    "Y_uudr": (1, "delta", 0), "Z_uuds": (2, "delta", 1), "M_uuds": (4, "delta", 1), "N_uudr": (5, "delta", 0),
    "Y_uvf": (1, "nu", 1), "Z_uwf": (2, "nu", 2), "Y_urf": (1, "nu", 5), "Z_uqf": (2, "nu", 4),
    "M_uwf": (4, "nu", 2), "N_uvf": (5, "nu", 1), "M_uqf": (4, "nu", 4), "N_urf": (5, "nu", 5),
}


def _remus_coefficients(x_fin):
    """The twelve linear coefficients of four fins at x_fin (two lifting
    along +y on one command, two along -z on the other), at u = 1 m/s."""
    fins = [prestero_fin([x_fin, 0, 0], E_Y), prestero_fin([x_fin, 0, 0], E_Y),
            prestero_fin([x_fin, 0, 0], NEG_Z), prestero_fin([x_fin, 0, 0], NEG_Z)]
    fin_set = producer_set(fins, [[1, 0], [1, 0], [0, 1], [0, 1]])
    command, nu_r = ca.SX.sym("command", 2), ca.SX.sym("nu_r", 6)
    out = fin_set(command=command, state=ca.SX(0, 1), nu_r=nu_r, water_density=prestero_value("water_density"))["tau"]
    jac = ca.Function("coefficients", [command, nu_r], [ca.jacobian(out, command), ca.jacobian(out, nu_r)])
    j_delta, j_nu = (np.array(m, dtype=float) for m in jac([0, 0], [1, 0, 0, 0, 0, 0]))
    return {name: (j_delta if kind == "delta" else j_nu)[row, col] for name, (row, kind, col) in TABLE_ENTRIES.items()}


def test_G5_table_C10_ratios_follow_from_the_composition():
    table = prestero()["control_fin_coefficients"]
    half = prestero()["printed_half_unit"]
    x_fin, dx = prestero_value("fin_position_x"), half["Table A.5 x_finpost"]
    scale, d_scale = table["Y_uudr"]["value"], half["Table C.10 9.64"]
    ratios = [{k: v / c["Y_uudr"] for k, v in c.items()} for c in (_remus_coefficients(x) for x in (x_fin - dx, x_fin + dx))]
    for name in TABLE_ENTRIES:
        printed = table[name]["value"]
        candidates = [s * r[name] for s in (scale - d_scale, scale + d_scale) for r in ratios]
        low, high = min(candidates) - 0.005, max(candidates) + 0.005       # half a printed unit
        assert low <= printed <= high, (name, printed, low, high)


def test_composition_scale_is_eq_4_44_and_table_C10_scale_is_not():
    """Y_uudr of the pair = rho c_L_alpha S_fin (eq. 4.44, p. 33) = 21.4
    kg/(m rad) on Tables A.1 and A.5; Table C.10 prints 9.64 (an
    inconsistency of the source, pinned here)."""
    rho, slope, area = (prestero_value(n) for n in ("water_density", "lift_slope", "fin_area"))
    got = _remus_coefficients(prestero_value("fin_position_x"))["Y_uudr"]
    assert abs(got - rho * slope * area) <= G2_TOLERANCE * rho * slope * area
    assert got / prestero()["control_fin_coefficients"]["Y_uudr"]["value"] > 2.0


# --------------------------------------------------------------------------
# One parameter open: gradient in one call
# --------------------------------------------------------------------------
@pytest.mark.parametrize("section", ["linear_section", "quadratic_drag"])
def test_gradient_with_respect_to_an_open_section_parameter(section):
    s = skeleton()
    from more_transformations.more_casadi_transformations import freeze

    values = s.check_lifting_fin_values(geometry([-0.75, 0.02, -0.03], NEG_Y, 0.013, E_X))
    parts = {"servo": part("servo", "ideal", {"max_deflection": MAX_DEFLECTION}),
             "inflow": part("inflow", "rigid_point"), "flow_angle": part("flow_angle", "small_angle"),
             "interference": part("interference", "none"),
             "section": part("section", section)}   # every section parameter open
    open_fin = freeze(s.lifting_fin_casadi(parts), s.lifting_fin_parameters(), values)
    assert "section.lift_slope" in open_fin.name_in()
    x0 = np.array([0.2, 1.4, 0.2, -0.1, 0.05, 0.1, 0.3])
    slope = ca.SX.sym("slope")
    extra = {"section.zero_lift_drag": 0.01} if section == "linear_section" else {}
    out = open_fin(command=x0[0], state=ca.SX(0, 1), nu_r=x0[1:7], water_density=1026.0,
                   **{"section.lift_slope": slope}, **extra)["tau"]
    f = ca.Function("f", [slope], [out, ca.jacobian(out, slope)])
    gradient = np.array(f(2.9)[1], dtype=float).reshape(-1)
    h = 1e-6
    central = (np.array(f(2.9 + h)[0], dtype=float) - np.array(f(2.9 - h)[0], dtype=float)).reshape(-1) / (2 * h)
    assert np.max(np.abs(gradient)) > 0.1
    assert max_diff(gradient, central) <= 1e-6


def test_open_part_parameters_are_named_by_slot_in_a_set():
    """A set exposes producer i's open input '<slot>.<name>' as
    'producers[i].<slot>.<name>'."""
    s = skeleton()
    from more_transformations.more_casadi_transformations import freeze

    parts = {"servo": part("servo", "ideal", {"max_deflection": MAX_DEFLECTION}),
             "inflow": part("inflow", "translational"), "flow_angle": part("flow_angle", "none"),
             "interference": part("interference", "none"), "section": part("section", "quadratic_drag")}
    open_fin = freeze(s.lifting_fin_casadi(parts), s.lifting_fin_parameters(),
                      s.check_lifting_fin_values(geometry([-0.75, 0, 0], NEG_Y, 0.013, E_X)))
    m = module("force_producer_set")
    block = m.force_producer_set_casadi([open_fin, open_fin], 2)
    assert "producers[0].section.lift_slope" in block.name_in()
    assert "producers[1].section.lift_slope" in block.name_in()
