"""Linearised open-water propeller with given coefficients as a leaf of the
``force_producers`` list slot, fed ``nu_r`` and ``water_density`` by the
vehicle rather than its own copies.

References
----------
[MSS] Fossen, T. I. MSS, MIT, @ cc07579: CRAFT/AUV/models/remus100.m
    148-176 (propeller).

Author:    Enio Krizman
Date:      2026-10-08
"""

from more_dynamics.models.propeller.propeller import propeller_casadi, propeller_parameters
from more_dynamics.models.shared.as_producer import as_producer


def propeller_linearized_open_water():
    """Linearised open-water propeller, given coefficients (remus100.m
    148-176)."""
    return as_producer(
        propeller_casadi(open_water_model="linearized", open_water_coefficients="given"),
        "n", "propeller_linearized_open_water")


def propeller_linearized_open_water_parameters():
    """Every declared parameter but ``water_density`` (fed by the ``site``
    coupling, not frozen here)."""
    return tuple(d for d in propeller_parameters(open_water_coefficients="given") if d.name != "water_density")
