"""Linearised law: for ``n > 0`` the open-water coefficients move linearly from their value at ``J = 0`` to ``J = Ja_max`` with the advance speed, the bollard law for ``n <= 0`` (MSS ``remus100.m`` 168-178).
Declared numbers and the models are in ``open_water``.

References
----------
[MSS] Fossen, T. I. (2026). Marine Systems Simulator, MIT, CRAFT/AUV/models/remus100.m 149-178
    @ 99fdf4d.

Author:    Enio Krizman
Date:      2026-10-23
"""
from more_dynamics.models.propulsor.law.open_water import LAW_INPUT_UNIT, open_water_law, open_water_parameters  # noqa: F401


def open_water_linearized_parameters(*, open_water_coefficients="polynomial"):
    return open_water_parameters(open_water_coefficients=open_water_coefficients)


def open_water_linearized_casadi(*, open_water_coefficients="polynomial"):
    return open_water_law("open_water_linearized", "linearized", open_water_coefficients)
