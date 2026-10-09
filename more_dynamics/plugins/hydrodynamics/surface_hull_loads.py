"""The hydrodynamic loads on the hull of a surface craft, as Luka's
``HydrodynamicsModel``: linear damping, either calibrated from the top speed
and thrust with a quadratic yaw term (MSS ``otter.m`` 195-240) or, by
``manoeuvring_damping``, a coupled linear sway/yaw law
(``models/damping/hydrodynamic_load_parts.py``); the cross-flow strip
integral with the rectangular Hoerner section (``crossFlowDrag.m`` 54-69,
``Hoerner.m`` 47-51), switched off by ``cross_flow = "none"``; and optionally
the ITTC-1957 surge resistance, plain or with the RA14 residual term
(Fossen 2011, eqs. 6.82-6.85, p. 125; Radojcic et al. 2014). The first output
is the sum, the signed generalized force the vehicle adds
(``vehicle_model.capnp`` 18-21); the loads follow as ``damping_force``,
``cross_flow_force`` and ``surge_resistance_force``, the damping matrix as
``damping_matrix``.

Inputs by name, from the vehicle and the hydrostatics: ``relative_velocity``,
``mass_matrix``, ``restoring_matrix_at_flotation``, ``length``, ``draft``,
``water_density``, ``section_beam``; with ``surge_resistance`` in
``("ittc", "ittc_residual")`` also ``mass``, ``wetted_surface`` and
``kinematic_viscosity``. A vehicle that gives no ``wetted_surface`` is
refused with the missing name. The ITTC/residual resistance is added to the
calibrated damping, which already contains a surge term; whether the two are
used together is the composition's choice, not a default (``"none"``).

Fidelity: the time-constant calibration is an MSS calibration (raise:
identification from logs); ``linear_coupled`` is a given or identified 2x2
sway/yaw law, ``Y_r = N_v = 0`` until identified (raise: circle-test
identification); the Hoerner section is a published table (raise: a
tank-measured section drag); the ITTC line is a physics form with a printed
form factor (raise: a measured resistance curve); the RA14 residual term is a
published regression keyed on an unmeasured length (raise: an identified
coast-down curve).

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics
    and Motion Control. Wiley. Eqs. 6.82-6.85, p. 125; eqs. 6.91-6.92, p. 127.
[MSS] Fossen, T. I. MSS, MIT, @ cc07579: CRAFT/USV/models/otter.m 99-101,
    195-240; LIBRARY/modeling/crossFlowDrag.m 54-69; Hoerner.m 47-51.
[RA14] Radojcic, D., Zgradic, A., Kalajdzic, M., Simic, A. (2014). Resistance
    prediction for hard chine hulls in the pre-planing regime. Polish
    Maritime Research 21(2):9-26. Appendix 1, p. 24 (through
    ``models/surge_resistance/ra14_residual_resistance.py``).

Author:    Enio Krizman
Date:      2026-10-09
"""
from rpp_plugin_types.more_dynamics import HydrodynamicsModel
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_dynamics.models.cross_flow.cross_flow_strip import cross_flow_strip_rectangular_section_hoerner
from more_dynamics.models.damping.hydrodynamic_load_parts import (
    coupled_manoeuvring_damping_surface, coupled_manoeuvring_damping_surface_parameters,
    surge_resistance_ittc, surge_resistance_ittc_parameters, surge_resistance_ittc_residual,
    surge_resistance_ittc_residual_parameters, time_constant_damping_surface,
    time_constant_damping_surface_parameters)
from more_dynamics.plugins.shared.payload_io import PayloadBuilder, frozen_block, payload_name

SURGE_RESISTANCES = ("none", "ittc", "ittc_residual")
MANOEUVRING_DAMPINGS = ("time_constants", "linear_coupled", "linear_coupled_modulus")
CROSS_FLOW_FORMS = ("hoerner_strips", "none")

_SURGE_RESISTANCE_FUNCTIONS = {"ittc": surge_resistance_ittc, "ittc_residual": surge_resistance_ittc_residual}
_SURGE_RESISTANCE_PARAMETER_FUNCTIONS = {
    "ittc": surge_resistance_ittc_parameters, "ittc_residual": surge_resistance_ittc_residual_parameters}


class SurfaceHullLoads(HydrodynamicsModel):
    PARAMETERS = [
        ParameterDescription("max_forward_thrust", 239.364),  # otter.m, 24.4 kgf = 24.4 * 9.81 N
        ParameterDescription("max_speed", 3.0864),  # otter.m, 6 knots = 6 * 0.5144 m/s
        ParameterDescription("time_constants", [1.0, 1.0]),  # otter.m 99-101, 202-240
        ParameterDescription("damping_ratios", [0.3, 0.2, 0.4]),  # otter.m 99-101, 202-240
        ParameterDescription("yaw_damping_nonlinearity", 10.0),  # otter.m 99-101, 202-240
        ParameterDescription("manoeuvring_damping", "time_constants"),
        ParameterDescription("sway_yaw_derivatives", [0.0, 0.0, 0.0, 0.0]),  # read with "linear_coupled"/"_modulus"
        ParameterDescription("modulus_coefficients", [0.0, 0.0]),  # read with "linear_coupled_modulus"; zero (module docstring)
        ParameterDescription("cross_flow", "hoerner_strips"),
        ParameterDescription("surge_resistance", "none"),
        ParameterDescription("time_constant", 100.0),  # osv.m 128, vessel.T1 = 100 s (read with "ittc"/"ittc_residual")
        ParameterDescription("form_factor", 0.1),  # XuuITTC.m 33, k = 0.1; Fossen 2011, p. 125 (read with "ittc"/"ittc_residual")
        ParameterDescription("crossover_speed", 2.0),  # forceSurgeDamping.m 57 @ ac77394, u_cross = 2 m/s (read with "ittc"/"ittc_residual")
        ParameterDescription("surge_added_mass_factor", 2.7),  # addedMassSurge.m 33, A11 = 2.7 rho nabla^(5/3) / L^2 (read with "ittc"/"ittc_residual")
        ParameterDescription("ittc_reynolds_floor", 100000.0),  # XuuITTC.m 35, Re_min = 1e5 (read with "ittc"/"ittc_residual")
    ]

    def __init__(self) -> None:
        self._damping = self._cross_flow = self._surge = None

    def initialize(self, context: ComponentContext) -> None:
        surge = context.get_parameter("surge_resistance")
        if surge not in SURGE_RESISTANCES:
            raise ValueError(f"surge_resistance must be one of {SURGE_RESISTANCES}, got {surge!r}")
        manoeuvring = context.get_parameter("manoeuvring_damping")
        if manoeuvring not in MANOEUVRING_DAMPINGS:
            raise ValueError(f"manoeuvring_damping must be one of {MANOEUVRING_DAMPINGS}, got {manoeuvring!r}")
        cross_flow_form = context.get_parameter("cross_flow")
        if cross_flow_form not in CROSS_FLOW_FORMS:
            raise ValueError(f"cross_flow must be one of {CROSS_FLOW_FORMS}, got {cross_flow_form!r}")
        if manoeuvring == "time_constants":
            self._damping = frozen_block(context, time_constant_damping_surface(),
                                         time_constant_damping_surface_parameters())
        else:
            modulus = manoeuvring == "linear_coupled_modulus"
            self._damping = frozen_block(context, coupled_manoeuvring_damping_surface(modulus=modulus),
                                         coupled_manoeuvring_damping_surface_parameters(modulus=modulus))
        self._cross_flow = (cross_flow_strip_rectangular_section_hoerner()
                            if cross_flow_form == "hoerner_strips" else None)
        self._surge = None
        if surge in _SURGE_RESISTANCE_FUNCTIONS:
            self._surge = frozen_block(context, _SURGE_RESISTANCE_FUNCTIONS[surge](),
                                       _SURGE_RESISTANCE_PARAMETER_FUNCTIONS[surge]())

    def graph(self) -> HydrodynamicsModel.CasadyPayload:
        if self._damping is None:
            raise RuntimeError("SurfaceHullLoads must be initialized before graph()")
        io = PayloadBuilder(HydrodynamicsModel.CasadyPayload())
        damping = io.call(self._damping, rename={"restoring_matrix": "restoring_matrix_at_flotation"})
        total = damping["tau"]
        cross_flow = None
        if self._cross_flow is not None:
            cross_flow = io.call(self._cross_flow, rename={"beam": "section_beam"})
            total = total + cross_flow["tau"]
        surge = None
        if self._surge is not None:
            surge = io.call(self._surge)
            total = total + surge["tau"]
        io.output("hydrodynamic_force", total, "signed generalized force added by the vehicle, BODY")
        io.output("damping_force", damping["tau"])
        if cross_flow is not None:
            io.output("cross_flow_force", cross_flow["tau"])
        if surge is not None:
            io.output("surge_resistance_force", surge["tau"])
        for name in ("D", "damping_derivatives"):
            io.output(payload_name(name), damping[name])
        return io.payload()
