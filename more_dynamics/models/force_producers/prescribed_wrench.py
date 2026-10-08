"""Stage-1 force producer: the commanded wrench acts on the vehicle unchanged
(no actuator model), the external-wrench boundary of ``hydroVessel.m``.

Author:    Enio Krizman
Date:      2026-10-08
"""

import casadi as ca

from more_dynamics.models.shared.wiring import function_from


def prescribed_wrench():
    """Stage 1: the command is the wrench on the vehicle (no actuator
    model)."""
    command = ca.SX.sym("command", 6)
    return function_from("prescribed_wrench", {"command": command}, {"tau": command})


PRESCRIBED_WRENCH_PARAMETERS = ()
