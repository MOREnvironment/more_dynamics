"""Deflection-only fin pairs as a leaf of the ``force_producers`` list slot, fed
``nu_r`` and ``water_density`` by the vehicle rather than its own copies.

References
----------
[MSS] Fossen, T. I. MSS, MIT, @ cc07579: CRAFT/AUV/models/remus100.m
    228-245 (fin pairs).

Author:    Enio Krizman
Date:      2026-10-08
"""

import casadi as ca

from more_dynamics.models.fin.fins import fins_casadi, fins_parameters
from more_dynamics.models.shared.as_producer import as_producer
from more_dynamics.models.shared.wiring import function_from


def fin_pairs_deflection_only():
    """Deflection-only fin pairs (remus100.m 228-245)."""
    return as_producer(fins_casadi(convention="starboard_down_positive"), "delta", "fin_pairs_deflection_only")


def fin_pairs_deflection_only_parameters(*, positions_given=True):
    """Every declared parameter but ``water_density`` (fed by the ``site``
    coupling, not frozen here); with ``positions_given=False`` the two fin
    positions are not declared either: ``fin_pair_positions`` computes them."""
    own = tuple(d for d in fins_parameters() if d.name != "water_density")
    if positions_given:
        return own
    return tuple(d for d in own if d.name not in ("rudder_position", "stern_plane_position"))



def fin_pair_positions():
    """``(semi_major_axis) -> (rudder_position, stern_plane_position)``: both fin pairs at the tail of the body,
    ``x = -a`` with ``a`` the body's half length (``remus100.m`` 184, 189)."""
    semi_major_axis = ca.SX.sym("semi_major_axis")
    return function_from("fin_pair_positions", {"semi_major_axis": semi_major_axis},
                         {"rudder_position": -semi_major_axis,  # x_r = -a (remus100.m 184)
                          "stern_plane_position": -semi_major_axis})  # x_s = -a (remus100.m 189)
