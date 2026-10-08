"""Force producers: thrusters, propellers, fins and an outboard motor as all-CasADi blocks.

Author:    Enio Krizman
Date:      2026-10-05
"""

from .thrusters.differential_thruster import (
    check_differential_thruster_values,
    differential_thruster_casadi,
    differential_thruster_parameters,
)
from .fins import fins_casadi, fins_parameters
from .thrusters.outboard_motor import (
    outboard_motor_rpm_casadi,
    outboard_motor_rpm_parameters,
    outboard_motor_throttle_casadi,
    outboard_motor_throttle_parameters,
)
from .propulsor.propeller import propeller_casadi, propeller_parameters
from .vsim_fins import vsim_fins_casadi, vsim_fins_parameters
from .propulsor.wageningen_kt_kq import wageningen_casadi, wageningen_expression, wageningen_parameters

__all__ = [
    "check_differential_thruster_values",
    "differential_thruster_casadi",
    "differential_thruster_parameters",
    "fins_casadi",
    "fins_parameters",
    "outboard_motor_rpm_casadi",
    "outboard_motor_rpm_parameters",
    "outboard_motor_throttle_casadi",
    "outboard_motor_throttle_parameters",
    "propeller_casadi",
    "propeller_parameters",
    "vsim_fins_casadi",
    "vsim_fins_parameters",
    "wageningen_casadi",
    "wageningen_expression",
    "wageningen_parameters",
]
