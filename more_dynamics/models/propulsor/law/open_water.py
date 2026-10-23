"""Open-water propeller laws of the law socket, the three of the propeller block as shaft-speed laws:
the shaft speed (rev/s) and the advance speed give the thrust and the shaft torque.

``K_T`` and ``K_Q`` of the Wageningen polynomial; ``T = rho D^4 K_T n |n|``, ``Q = rho D^5 K_Q n |n|``
(Fossen 2011, eq. 12.265, p. 411). The models are those of ``models.propeller.propeller``
(``"bollard"``, ``"linearized"``, ``"full"``, MSS ``remus100.m`` 149-178); the open-water
coefficients ``K_T0``, ``K_Q0``, ``K_Tmax``, ``K_Qmax`` come from the polynomial at ``J = 0`` and
``J = Ja_max`` (``open_water_coefficients="polynomial"``) or are declared as
``thrust_torque_coefficients`` (``"given"``, the four numbers the model file prints).

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics and Motion Control.
    Wiley. Ch. 12, eq. 12.265, p. 411.
[MSS] Fossen, T. I. (2026). Marine Systems Simulator, MIT, CRAFT/AUV/models/remus100.m 149-178
    @ 99fdf4d.

Author:    Enio Krizman
Date:      2026-10-23
"""
import casadi as ca
from more_transformations.more_casadi_transformations import Parameter, symbols

from more_dynamics.models.propeller.propeller import OPEN_WATER_COEFFICIENTS, open_water_loads
from more_dynamics.models.propeller.wageningen_kt_kq import wageningen_expression, wageningen_parameters
from more_dynamics.models.propulsor.declarations import part_function, positive

LAW_INPUT_UNIT = "rev/s"


def _check(open_water_coefficients):
    if open_water_coefficients not in OPEN_WATER_COEFFICIENTS:
        raise ValueError(f"open_water_coefficients must be one of {OPEN_WATER_COEFFICIENTS}, "
                         f"got {open_water_coefficients!r}")


def open_water_parameters(*, open_water_coefficients="polynomial"):
    _check(open_water_coefficients)
    declared = (positive("diameter", "m", "propeller diameter D"), *wageningen_parameters(),
                positive("max_advance_number", "1", "Ja_max = (1 - w) U_max / (D n_max), J >= 0 on the forward branch"))
    if open_water_coefficients == "given":
        declared += (Parameter("thrust_torque_coefficients", (4, 1), "1", "[K_T0, K_Q0, K_Tmax, K_Qmax]"),)
    return declared


def open_water_law(name, model, open_water_coefficients):
    """The law function of ``model`` (``"bollard"``, ``"linearized"``, ``"full"``)."""
    declared = open_water_parameters(open_water_coefficients=open_water_coefficients)
    p = symbols(declared)
    n, advance, density = ca.SX.sym("drive_output"), ca.SX.sym("advance_speed"), ca.SX.sym("water_density")
    geometry = (p["pitch_diameter_ratio"], p["blade_area_ratio"], p["blade_count"])
    if open_water_coefficients == "given":
        coefficients = p["thrust_torque_coefficients"]
    else:
        kt0, kq0 = wageningen_expression(0.0, *geometry)  # at J = 0 (remus100.m 157)
        kt_max, kq_max = wageningen_expression(p["max_advance_number"], *geometry)  # at J = Ja_max (remus100.m 160)
        coefficients = ca.vertcat(kt0, kq0, kt_max, kq_max)
    block = {"water_density": density, "propeller_diameter": p["diameter"], "max_advance_number": p["max_advance_number"],
             "pitch_diameter_ratio": p["pitch_diameter_ratio"], "blade_area_ratio": p["blade_area_ratio"],
             "blade_count": p["blade_count"]}
    thrust, torque = open_water_loads(block, coefficients, model, False, n, advance)
    return part_function(name, {"drive_output": n, "advance_speed": advance, "water_density": density},
                         {"thrust": thrust, "torque": torque}, p, declared)
