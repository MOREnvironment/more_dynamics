"""Gate tests of DUNE's simulator fins as one-fin leaves in a
``ForceProducerSet``.

Written 2026-10-08, before the parts exist (contract: ``fin_parts_contract``).
DUNE's simulated fin is not a ``LiftingFin``: its moments are not ``r x f``
but its own arm law, so a one-fin leaf ``proportional_force_fin`` (a force
producer with no sockets) carries it, and a vehicle's fins are a set of
leaves with an identity command map (DUNE commands each fin; the roll,
pitch and yaw allocation stays in the control library).

Gates:

* G2: one leaf equals today's block ``vsim_fins_casadi(fin_count=1)`` on the
  same numbers.
* G2: the set of four leaves, identity map, equals ``vsim_fins_casadi(
  fin_count=4)`` on the LAUV fin forces and positions of DUNE's simulator
  configuration (``tests/force_producers/data/parameter_sets.json``,
  ``vsim_fins``: ``etc/common/vsim-models.ini`` 54-62), including speeds
  below the block's zero-speed guard.
* The leaf takes ``water_density`` (the producer interface) and its force
  does not depend on it: DUNE's force gain holds the density it was made
  with.

References
----------
[DUNE] LSTS, Universidade do Porto. *DUNE: Unified Navigation Environment*,
    EUPL v1.1, ``src/Simulators/VSIM/VSIM/Fin.cpp`` 56-102 (B. Terra,
    J. Braga), ``etc/common/vsim-models.ini`` 54-62, revision ``555ef0b`` of
    the vehicle's DUNE source (read through the existing block, not copied).

Author:    Enio Krizman
Date:      2026-10-08
"""

import numpy as np

from fin_parts_contract import (
    G2_TOLERANCE,
    LEAF_DECLARED,
    LEAF_MODULE,
    N_RANDOM,
    leaf,
    max_diff,
    module,
    parameter_set,
    producer_set,
    random_nu_r,
    rng,
    tau,
)


def _vsim_block(fin_count, values):
    from more_dynamics.models.fin.vsim_fins import vsim_fins_casadi, vsim_fins_parameters
    from more_transformations.more_casadi_transformations import freeze

    return freeze(vsim_fins_casadi(fin_count=fin_count), vsim_fins_parameters(fin_count=fin_count), values)


def _leaves(p):
    return [leaf({"max_force": p["max_forces"][i], "fin_position": p["positions"][i],
                  "max_deflection": p["max_deflection"]}) for i in range(len(p["max_forces"]))]


def _cases(n_fins):
    g = rng()
    nu_r = random_nu_r(g, N_RANDOM)
    nu_r[:20, :3] = g.uniform(-1e-7, 1e-7, (20, 3))          # below the zero-speed guard
    return g.uniform(-0.6, 0.6, (N_RANDOM, n_fins)), nu_r


def test_leaf_declares_its_parameters():
    declared = module(LEAF_MODULE).proportional_force_fin_parameters()
    assert [d.name for d in declared] == list(LEAF_DECLARED)
    for d in declared:
        shape, unit, minimum, exclusive = LEAF_DECLARED[d.name]
        assert (d.shape, d.unit, d.minimum, d.minimum_exclusive) == (shape, unit, minimum, exclusive), d
        assert d.meaning


def test_G2_one_leaf_equals_the_one_fin_block():
    p = parameter_set("vsim_fins")
    one = {"max_forces": [p["max_forces"][0]], "positions": [p["positions"][0]], "max_deflection": p["max_deflection"]}
    block, single = _vsim_block(1, one), _leaves(p)[0]
    for delta, nu_r in zip(*_cases(1)):
        expected = np.array(block(delta=delta, nu_r=nu_r)["tau"], dtype=float).reshape(-1)
        assert max_diff(tau(single, delta, nu_r, 1026.0), expected) <= G2_TOLERANCE


def test_G2_set_of_four_leaves_equals_the_four_fin_block():
    p = parameter_set("vsim_fins")
    block = _vsim_block(4, p)
    fin_set = producer_set(_leaves(p), np.eye(4))
    for delta, nu_r in zip(*_cases(4)):
        expected = np.array(block(delta=delta, nu_r=nu_r)["tau"], dtype=float).reshape(-1)
        assert max_diff(tau(fin_set, delta, nu_r, 1026.0), expected) <= G2_TOLERANCE


def test_leaf_force_does_not_depend_on_the_fed_density():
    single = _leaves(parameter_set("vsim_fins"))[1]
    assert "water_density" in single.name_in()
    for delta, nu_r in zip(*_cases(1)):
        assert max_diff(tau(single, delta, nu_r, 1000.0), tau(single, delta, nu_r, 1030.0)) == 0.0
