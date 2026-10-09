"""The hydrodynamic loads on the hull of a slender submerged vehicle, as
Luka's ``HydrodynamicsModel``: time-constant linear damping (``Dmtrx.m``;
``remus100.m`` 217-218), slender-body lift and induced drag
(``forceLiftDrag.m``) and the cross-flow strip integral with the
circular-cylinder section (``crossFlowDrag.m`` 54-69, ``cylinderDrag.m``
78-110). The first output is the sum, the signed generalized force the
vehicle adds (``vehicle_model.capnp`` 18-21); the three loads follow as
``damping_force``, ``lift_drag_force`` and ``cross_flow_force``, and the
damping matrix as ``damping_matrix``.

Inputs by name, from the vehicle and the hydrostatics: ``relative_velocity``,
the two mass matrices, ``weight``, ``center_of_gravity``,
``center_of_buoyancy``, ``span``, ``water_density``, ``length``, ``draft``,
``section_beam``, ``kinematic_viscosity``. Own parameters: the damping time
constants and ratios, and the planform, parasitic drag and Oswald factor of
the lift.

Fidelity: the time-constant damping is an MSS calibration (raise: damping
identified from decay tests); the lift and the cross-flow are physics forms
with printed coefficients (raise: towing-tank or log identification).

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics
    and Motion Control. Wiley. Eqs. 6.91-6.92, p. 127 (strip integral).
[MSS] Fossen, T. I. MSS, MIT, @ cc07579: LIBRARY/modeling/Dmtrx.m;
    CRAFT/AUV/models/remus100.m 133, 143-144, 192-196, 217-221;
    LIBRARY/modeling/forceLiftDrag.m; coeffLiftDrag.m 54;
    LIBRARY/modeling/crossFlowDrag.m 54-69; cylinderDrag.m 78-110.

Author:    Enio Krizman
Date:      2026-10-09
"""
from rpp_plugin_types.more_dynamics import HydrodynamicsModel
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_dynamics.models.cross_flow.cross_flow_strip import cross_flow_strip_circular_cylinder_reynolds
from more_dynamics.models.damping.hydrodynamic_load_parts import (
    hull_lift_drag, hull_lift_drag_parameters, time_constant_damping_submerged,
    time_constant_damping_submerged_parameters)
from more_dynamics.plugins.shared.payload_io import PayloadBuilder, frozen_block, payload_name


class AuvHullLoads(HydrodynamicsModel):
    PARAMETERS = [
        ParameterDescription("time_constants", [20, 20, 1]),  # remus100.m:192, 193, 196 [T1 T2 T6] (Dmtrx.m call at :217)
        ParameterDescription("damping_ratios", [0.3, 0.8]),  # remus100.m:194, 195 [zeta4 zeta5]
        ParameterDescription("planform_area", 0.2128),  # remus100.m:133 S = 0.7 * L_auv * D_auv
        ParameterDescription("parasitic_drag_coefficient", 0.05595961914206819),  # remus100.m:143-144 CD_0 = Cd * pi * b^2 / S, b = D_auv / 2 (one geometry)
        ParameterDescription("oswald_efficiency", 0.3),  # coeffLiftDrag.m:54 e = 0.3
    ]

    def __init__(self) -> None:
        self._damping = self._lift_drag = self._cross_flow = None

    def initialize(self, context: ComponentContext) -> None:
        self._damping = frozen_block(context, time_constant_damping_submerged(),
                                     time_constant_damping_submerged_parameters())
        self._lift_drag = frozen_block(context, hull_lift_drag(), hull_lift_drag_parameters())
        self._cross_flow = cross_flow_strip_circular_cylinder_reynolds()

    def graph(self) -> HydrodynamicsModel.CasadyPayload:
        if self._damping is None:
            raise RuntimeError("AuvHullLoads must be initialized before graph()")
        io = PayloadBuilder(HydrodynamicsModel.CasadyPayload())
        rename = {"beam": "section_beam"}
        damping = io.call(self._damping)
        lift_drag = io.call(self._lift_drag)
        cross_flow = io.call(self._cross_flow, rename=rename)
        io.output("hydrodynamic_force", damping["tau"] + lift_drag["tau"] + cross_flow["tau"],
                  "signed generalized force added by the vehicle, BODY")
        io.output("damping_force", damping["tau"])
        io.output("lift_drag_force", lift_drag["tau"])
        io.output("cross_flow_force", cross_flow["tau"])
        for name in ("D", "damping_coefficients"):
            io.output(payload_name(name), damping[name])
        return io.payload()
