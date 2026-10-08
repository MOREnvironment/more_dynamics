"""Deflection-only fin pairs as a leaf of the ``force_producers`` list slot, fed
``nu_r`` and ``water_density`` by the vehicle rather than its own copies.

References
----------
[MSS] Fossen, T. I. MSS, MIT, @ cc07579: CRAFT/AUV/models/remus100.m
    228-245 (fin pairs).

Author:    Enio Krizman
Date:      2026-10-08
"""

from more_dynamics.models.force_producers.fins import fins_casadi, fins_parameters
from more_dynamics.models.force_producers.shared.as_producer import as_producer


def fin_pairs_deflection_only():
    """Deflection-only fin pairs (remus100.m 228-245)."""
    return as_producer(fins_casadi(convention="starboard_down_positive"), "delta", "fin_pairs_deflection_only")


def fin_pairs_deflection_only_parameters():
    """Every declared parameter but ``water_density`` (fed by the ``site``
    coupling, not frozen here)."""
    return tuple(d for d in fins_parameters() if d.name != "water_density")

