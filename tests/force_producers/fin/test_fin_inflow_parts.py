"""Gate tests of the fin inflow parts (vessel motion -> water velocity at
the fin): ``translational``, ``rigid_point``.

Written 2026-10-08, before the parts exist (contract: ``fin_parts_contract``).
Coupling inputs ``nu_r`` (6, BODY FRD, relative to the water) and
``fin_position`` (3, m, from the CO); output ``fin_velocity`` (3, m/s, BODY):
the velocity of the fin relative to the water.

* ``translational``: ``fin_velocity = nu_r[0:3]`` — the fin sees the
  vehicle's translation only (MSS ``remus100.m`` 234-235 use ``nu_r(1:3)``
  alone).
* ``rigid_point``: ``fin_velocity = nu_r[0:3] + omega x r`` with every term
  kept: ``u_fin = u + z_fin q - y_fin r``, ``v_fin = v + x_fin r - z_fin p``,
  ``w_fin = w + y_fin p - x_fin q`` (Prestero 2001, eq. 4.40, p. 32;
  Prestero drops the ``y_fin``, ``z_fin`` terms for REMUS, this part keeps
  them).

References
----------
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2
    with the fixes of 2026-10-07. https://github.com/cybergalactic/MSS, MIT
    licence, revision ``cc07579``: ``CRAFT/AUV/models/remus100.m`` 234-235.
[Prestero 2001] Prestero, T. (2001). *Verification of a six-degree of
    freedom simulation model for the REMUS autonomous underwater vehicle*.
    MSc thesis, MIT/WHOI. Ch. 4, eq. 4.40, p. 32.

Author:    Enio Krizman
Date:      2026-10-08
"""

import numpy as np
import pytest

from fin_parts_contract import G2_TOLERANCE, N_RANDOM, call, part, part_declared, random_nu_r, rng


@pytest.mark.parametrize("form", ["translational", "rigid_point"])
def test_inflow_declares_no_parameter(form):
    assert tuple(part_declared("inflow", form)) == ()


def test_G2_translational_inflow_is_the_vehicle_translation():
    inflow = part("inflow", "translational")
    g = rng()
    for nu_r, r in zip(random_nu_r(g, N_RANDOM), g.uniform(-1, 1, (N_RANDOM, 3))):
        got = call(inflow, nu_r=nu_r, fin_position=r)["fin_velocity"]
        assert np.max(np.abs(got - nu_r[0:3])) <= G2_TOLERANCE


def test_G2_rigid_point_inflow_equals_prestero_eq_4_40_all_terms():
    inflow = part("inflow", "rigid_point")
    g = rng()
    for nu_r, r in zip(random_nu_r(g, N_RANDOM), g.uniform(-1, 1, (N_RANDOM, 3))):
        u, v, w, p, q, rr = nu_r
        x, y, z = r
        expected = [u + z * q - y * rr,           # (Prestero 2001, eq. 4.40, p. 32)
                    v + x * rr - z * p,
                    w + y * p - x * q]
        got = call(inflow, nu_r=nu_r, fin_position=r)["fin_velocity"]
        assert np.max(np.abs(got - expected)) <= G2_TOLERANCE, (nu_r, r)


def test_rigid_point_inflow_equals_translational_without_rotation():
    rigid, translational = part("inflow", "rigid_point"), part("inflow", "translational")
    g = rng()
    for nu_r, r in zip(random_nu_r(g, 100), g.uniform(-1, 1, (100, 3))):
        nu_r[3:6] = 0.0
        a = call(rigid, nu_r=nu_r, fin_position=r)["fin_velocity"]
        b = call(translational, nu_r=nu_r, fin_position=r)["fin_velocity"]
        assert np.max(np.abs(a - b)) <= G2_TOLERANCE
