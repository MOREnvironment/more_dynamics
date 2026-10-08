"""Members of the ``force_producers`` list slot: a prescribed wrench (stage
1, no actuator model), and today's fin and propeller blocks as leaves — fed
``nu_r`` and ``water_density`` by the vehicle rather than their own copies.
Each block's own command input is renamed ``command``; the plugin layer
gives it a payload name.

References
----------
[MSS] Fossen, T. I. MSS, MIT, @ cc07579: CRAFT/AUV/models/remus100.m
    228-245 (fin pairs); 148-176 (propeller).

Author:    Enio Krizman
Date:      2026-10-08
"""

import casadi as ca

from more_dynamics.models.force_producers.fins import fins_casadi, fins_parameters
from more_dynamics.models.force_producers.propulsor.propeller import propeller_casadi, propeller_parameters

from more_dynamics.models.wiring import function_from


def prescribed_wrench():
    """Stage 1: the command is the wrench on the vehicle (no actuator
    model)."""
    command = ca.SX.sym("command", 6)
    return function_from("prescribed_wrench", {"command": command}, {"tau": command})


PRESCRIBED_WRENCH_PARAMETERS = ()


def _as_producer(block, command_name, name, keep=("tau",)):
    """``block`` as a ``ForceProducer`` leaf: its own command input renamed
    ``command``."""
    ins = {("command" if n == command_name else n): block.sx_in(i) for i, n in enumerate(block.name_in())}
    out = block.call({n: ins["command" if n == command_name else n] for n in block.name_in()})
    return function_from(name, ins, {k: out[k] for k in keep})


def fin_pairs_deflection_only():
    """Deflection-only fin pairs (remus100.m 228-245)."""
    return _as_producer(fins_casadi(convention="starboard_down_positive"), "delta", "fin_pairs_deflection_only")


def fin_pairs_deflection_only_parameters():
    """Every declared parameter but ``water_density`` (fed by the ``site``
    coupling, not frozen here)."""
    return tuple(d for d in fins_parameters() if d.name != "water_density")


def propeller_linearized_open_water():
    """Linearised open-water propeller, given coefficients (remus100.m
    148-176)."""
    return _as_producer(
        propeller_casadi(open_water_model="linearized", open_water_coefficients="given"),
        "n", "propeller_linearized_open_water")


def propeller_linearized_open_water_parameters():
    """Every declared parameter but ``water_density`` (fed by the ``site``
    coupling, not frozen here)."""
    return tuple(d for d in propeller_parameters(open_water_coefficients="given") if d.name != "water_density")
