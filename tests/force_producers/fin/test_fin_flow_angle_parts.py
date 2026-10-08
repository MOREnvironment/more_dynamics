"""Gate tests of the fin flow-angle parts (water velocity at the fin -> flow
angle and speed squared): ``none``, ``small_angle``.

Written 2026-10-08, before the parts exist (contract: ``fin_parts_contract``).
Coupling inputs ``fin_velocity`` (3, m/s, BODY), ``chord_axis`` and
``lift_axis`` (3, unit, BODY; the lift axis is the direction of the lift for
a positive deflection); outputs ``flow_angle`` (1, rad) and
``speed_squared`` (1, m^2/s^2). With ``u_c = fin_velocity . chord_axis`` and
``u_n = fin_velocity . lift_axis``:

* ``none``: ``flow_angle = 0``, ``speed_squared = u_c^2 + u_n^2`` — the speed
  in the fin's section plane (MSS ``remus100.m`` 234-235: ``U_rh^2 = u_r^2 +
  v_r^2`` for the rudder, ``U_rv^2 = u_r^2 + w_r^2`` for the stern planes).
  G1 against MATLAB.
* ``small_angle``: ``flow_angle = u_n / u_c``, ``speed_squared = u_c^2``
  (Prestero 2001, eqs. 4.41-4.43, pp. 32-33: ``beta_re = v_fin / u_fin``,
  ``beta_se = w_fin / u_fin``, lift on ``u^2``). The sign: the flow angle is
  positive when the fin moves along ``+lift_axis``, so the skeleton's
  ``alpha = delta - flow_angle`` gives Prestero's ``delta_re = delta_r -
  beta_re`` for a rudder lifting along +y and ``delta_se = delta_s +
  beta_se`` for stern planes lifting along -z (eq. 4.41).

Not settled by the source and stated here: the small-angle form is not
defined at ``u_c = 0`` (Prestero assumes small angles, p. 32); the part
must stay finite there (no NaN), its value at ``u_c = 0`` is the part's own
guard and is not gated beyond ``speed_squared = 0``.

References
----------
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2
    with the fixes of 2026-10-07. https://github.com/cybergalactic/MSS, MIT
    licence, revision ``cc07579``: ``CRAFT/AUV/models/remus100.m`` 234-235.
[Prestero 2001] Prestero, T. (2001). *Verification of a six-degree of
    freedom simulation model for the REMUS autonomous underwater vehicle*.
    MSc thesis, MIT/WHOI. Ch. 4, eqs. 4.40-4.43, pp. 32-33, Figs. 4-1, 4-2.

Author:    Enio Krizman
Date:      2026-10-08
"""

import numpy as np
import pytest

from fin_parts_contract import (
    E_X,
    E_Y,
    G1_TOLERANCE,
    G2_TOLERANCE,
    N_RANDOM,
    NEG_Y,
    NEG_Z,
    call,
    mss_columns,
    mss_fins,
    part,
    part_declared,
    rng,
)


def _orthonormal_pairs(generator, n):
    """Random perpendicular unit vectors (chord, lift)."""
    pairs = []
    for _ in range(n):
        q, _ = np.linalg.qr(generator.normal(size=(3, 3)))
        pairs.append((q[:, 0], q[:, 1]))
    return pairs


@pytest.mark.parametrize("form", ["none", "small_angle"])
def test_flow_angle_declares_no_parameter(form):
    assert tuple(part_declared("flow_angle", form)) == ()


def test_G1_no_flow_angle_speed_equals_mss_section_plane_speeds():
    """Rudder (lift -y) -> MATLAB U_rh^2; stern planes (lift -z) -> U_rv^2
    (remus100.m 234-235), fin velocity = nu_r[0:3]."""
    ref = mss_fins()
    part_ = part("flow_angle", "none")
    nu_r = mss_columns(ref, "nu_r", 6)
    for lift, column in ((NEG_Y, "U_rh"), (NEG_Z, "U_rv")):
        for k, v in enumerate(nu_r[:, 0:3]):
            out = call(part_, fin_velocity=v, chord_axis=E_X, lift_axis=lift)
            assert out["flow_angle"][0] == 0.0
            assert abs(out["speed_squared"][0] - ref[column][k] ** 2) <= G1_TOLERANCE, (column, k)


def test_G2_no_flow_angle_in_any_section_plane():
    part_ = part("flow_angle", "none")
    g = rng()
    for (c, n), v in zip(_orthonormal_pairs(g, 200), g.uniform(-3, 3, (200, 3))):
        out = call(part_, fin_velocity=v, chord_axis=c, lift_axis=n)
        assert abs(out["speed_squared"][0] - (v @ c) ** 2 - (v @ n) ** 2) <= G2_TOLERANCE


def test_G2_small_angle_equals_prestero_eq_4_42():
    """Rudder lifting along +y: flow angle = v_fin / u_fin (beta_re); stern
    planes lifting along -z: flow angle = -w_fin / u_fin (= -beta_se); both
    with speed squared u_fin^2 (eqs. 4.41-4.43)."""
    part_ = part("flow_angle", "small_angle")
    g = rng()
    v = g.uniform(-3, 3, (N_RANDOM, 3))
    v[:, 0] = np.where(np.abs(v[:, 0]) < 0.05, 0.05, v[:, 0])       # away from u_fin = 0
    for u_fin, v_fin, w_fin in v:
        rudder = call(part_, fin_velocity=[u_fin, v_fin, w_fin], chord_axis=E_X, lift_axis=E_Y)
        stern = call(part_, fin_velocity=[u_fin, v_fin, w_fin], chord_axis=E_X, lift_axis=NEG_Z)
        assert abs(rudder["flow_angle"][0] - v_fin / u_fin) <= G2_TOLERANCE * max(1.0, abs(v_fin / u_fin))
        assert abs(stern["flow_angle"][0] + w_fin / u_fin) <= G2_TOLERANCE * max(1.0, abs(w_fin / u_fin))
        for out in (rudder, stern):
            assert abs(out["speed_squared"][0] - u_fin ** 2) <= G2_TOLERANCE


def test_G2_small_angle_in_any_section_plane():
    part_ = part("flow_angle", "small_angle")
    g = rng()
    for (c, n), v in zip(_orthonormal_pairs(g, 200), g.uniform(-3, 3, (200, 3))):
        u_c, u_n = v @ c, v @ n
        if abs(u_c) < 0.05:
            continue
        out = call(part_, fin_velocity=v, chord_axis=c, lift_axis=n)
        assert abs(out["flow_angle"][0] - u_n / u_c) <= G2_TOLERANCE * max(1.0, abs(u_n / u_c))
        assert abs(out["speed_squared"][0] - u_c ** 2) <= G2_TOLERANCE


@pytest.mark.parametrize("v", [[0.0, 0.0, 0.0], [0.0, 0.4, 0.0], [0.0, 0.0, -0.3], [0.0, 0.2, 0.1]])
def test_small_angle_stays_finite_where_the_chordwise_speed_is_zero(v):
    import casadi as ca

    part_ = part("flow_angle", "small_angle")
    out = call(part_, fin_velocity=v, chord_axis=E_X, lift_axis=E_Y)
    assert np.all(np.isfinite(out["flow_angle"])) and out["speed_squared"][0] == 0.0
    x = ca.SX.sym("v", 3)
    jac = ca.Function("j", [x], [ca.jacobian(part_(fin_velocity=x, chord_axis=E_X, lift_axis=E_Y)["speed_squared"], x)])
    assert np.all(np.isfinite(np.array(jac(v), dtype=float)))
