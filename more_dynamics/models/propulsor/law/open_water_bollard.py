"""Bollard law: thrust and torque of the open-water coefficients at zero advance, ``T = rho D^4 K_T0 n |n|`` (MSS ``remus100.m`` 176-177).
Declared numbers and the models are in ``open_water``.

References
----------
[MSS] Fossen, T. I. (2026). Marine Systems Simulator, MIT, CRAFT/AUV/models/remus100.m 149-178
    @ 6e78e98.

Author:    Enio Krizman
Date:      2026-10-23
"""
from more_dynamics.models.propulsor.law.open_water import LAW_INPUT_UNIT, open_water_law, open_water_parameters  # noqa: F401


def open_water_bollard_parameters(*, open_water_coefficients="polynomial"):
    return open_water_parameters(open_water_coefficients=open_water_coefficients)


def open_water_bollard_casadi(*, open_water_coefficients="polynomial"):
    return open_water_law("open_water_bollard", "bollard", open_water_coefficients)
