"""One lifting fin as Luka's ``ForceProducer``, assembled from five child parts
(servo, inflow, flow angle, interference, section) in the way ``HullVessel``
takes its children: each child's ``graph()`` payload is read through its
descriptions and handed to the model-layer fin. The fin's wrench is the first
output (``generated_force``, BODY, N and N m about the CO); the servo's state
is the fin's state.

Inputs: ``fin_angle_command`` (rad, from outside), ``relative_velocity`` and
``water_density`` (fed by the vehicle by name). Outputs: ``generated_force``,
``deflection``, ``angle_of_attack``.

The force is resolved along fixed BODY axes, the small-angle form of the
cited sources, not a large-angle flow-axis model: f = 1/2 rho A U^2 (C_L
lift_axis - C_D chord_axis), tau = [f; r x f] (Prestero 2001, eqs. 4.37,
4.41-4.43, pp. 31-33; Fossen 2011, eq. 12.226, p. 400). A parameter of a
child part left open is refused: inside a fin every child is frozen.

Defaults: one REMUS 100 fin at the fin post, a rudder lifting along +y (the
sign of the lift axis is a convention of the composition).

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics
    and Motion Control. Wiley. Eq. 12.226, p. 400.
[Prestero 2001] Prestero, T. (2001). Verification of a six-degree of freedom
    simulation model for the REMUS autonomous underwater vehicle. MIT/WHOI
    MSc thesis. Eqs. 4.37, 4.41-4.43, pp. 31-33; Table A.5, p. 103.
[MSS] Fossen, T. I. MSS, MIT, CRAFT/AUV/models/remus100.m 238-254 @ cc07579.

Author:    Enio Krizman
Date:      2026-10-08
"""
import casadi as ca
from rpp_plugin_types.more_dynamics import ForceProducer
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_transformations.more_casadi_transformations import freeze
from more_dynamics.models.force_producers.fin.lifting_fin import (check_lifting_fin_values, lifting_fin_casadi,
                                                                       lifting_fin_parameters)
from more_dynamics.plugins.payload_io import MODEL_NAMES, PayloadBuilder, named_function

FED_INPUTS = ["command", "state", "nu_r", "water_density"]


class LiftingFin(ForceProducer):
    COMPONENTS = {"servo": "more_dynamics::ActuatorServo", "inflow": "more_dynamics::FinInflow",
                  "flow_angle": "more_dynamics::FinFlowAngle", "interference": "more_dynamics::FinInterference",
                  "section": "more_dynamics::FinSection"}
    PARAMETERS = [
        ParameterDescription("fin_position", [-0.638, 0.0, 0.0]),  # m, Prestero 2001, Table A.5, p. 103 (x_finpost)
        ParameterDescription("chord_axis", [1.0, 0.0, 0.0]),  # unit chord direction, BODY FRD (convention)
        ParameterDescription("lift_axis", [0.0, 1.0, 0.0]),  # unit positive-lift direction, BODY FRD (convention)
        ParameterDescription("fin_area", 0.00665),  # m^2, Prestero 2001, Table A.5, p. 103 (S_fin)
    ]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        parts = {}
        for slot in self.COMPONENTS:
            child = context.get_component(slot)
            if child == []:
                raise ValueError(f"LiftingFin: slot {slot!r} is empty: a fin needs a part in every slot")
            parts[slot] = named_function(child.graph(), slot, rename=MODEL_NAMES,
                                         state_name="servo_state" if slot == "servo" else None)
        declared = lifting_fin_parameters()
        values = check_lifting_fin_values({d.name: context.get_parameter(d.name) for d in declared})
        model = freeze(lifting_fin_casadi(parts), declared, values)
        left_open = [name for name in model.name_in() if name not in FED_INPUTS]
        if left_open:
            raise ValueError(f"LiftingFin: parameters {left_open} of its child parts are open; "
                             "a fin takes every child parameter from the composition")
        self._model = model

    def graph(self) -> ForceProducer.CasadyPayload:
        if self._model is None:
            raise RuntimeError("LiftingFin must be initialized before graph()")
        io = PayloadBuilder(ForceProducer.CasadyPayload())
        n = self._model.size1_in("state")
        state = io.state("servo_state", n, [0.0] * n, "servo state") if n else ca.SX.sym("state", 0)
        out = io.call(self._model, rename={"command": "fin_angle_command"}, given={"state": state})
        io.output("generated_force", out["tau"], "fin wrench about the CO, BODY, N and N m")
        io.output("deflection", out["deflection"], "fin deflection, rad")
        io.output("angle_of_attack", out["angle_of_attack"], "fin angle of attack, rad")
        if n:
            io.state_dot(out["state_dot"])
        return io.payload()
