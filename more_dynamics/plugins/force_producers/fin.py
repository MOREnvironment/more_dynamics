"""One lifting fin as Luka's ``ForceProducer``: the servo, the inflow at the
fin, the flow angle, the interference and the section are model functions of
``models/{servo,fin}`` chosen by selectors, their numbers declared here as
parameters. The fin's wrench is the first output (``generated_force``, BODY,
N and N m about the CO); the servo's state, when it has one, is the fin's
state.

Inputs: ``fin_angle_command`` (rad, from outside), ``relative_velocity`` and
``water_density`` (fed by the vehicle by name). Outputs: ``generated_force``,
``deflection``, ``angle_of_attack``.

Servo (one device whose features are switched on or off, Murray-Smith 2016,
Fig. 1, p. 246): ``servo_dynamics`` ``"none"`` or ``"first_order_lag"`` (needs
``time_constant``), ``servo_rate_limit`` ``False`` or ``True`` (needs the lag
and ``max_rate``), ``servo_angle_limit`` ``"none"``, ``"on_command"`` or
``"on_output"`` (needs ``max_deflection``). Defaults: the lag with rate and
angle limit on the command, the numbers of the AUV fin servo of Sarhadi
(2026): 0.1 s, +-20 deg, +-30 deg/s (Fig. 4 and text, p. 4; a demanding
choice of that study, not a measurement; raise: identify the time constant
and the rate limit of the real servo from a logged step in the fin command).

Forms, each with its fidelity (the next level needs the measurement named):

* ``inflow``: ``"translational"`` (the fin sees the vehicle's translation,
  ``remus100.m`` 234-235) or ``"rigid_point"`` (the velocity at a point fixed
  to the body, Prestero 2001, eq. 4.40, p. 32). Raise: the flow measured at
  the fin.
* ``flow_angle``: ``"none"`` (``remus100.m`` 234-245) or ``"small_angle"``
  (``u_n / u_c``, Prestero 2001, eqs. 4.42-4.43, pp. 32-33).
* ``interference``: ``"none"`` (both factors one, Prestero 2001, eq. 4.41,
  p. 32); ``"slender_body"`` (Pitts 1957) waits on its source and refuses.
* ``section``: ``"quadratic_drag"`` (lift slope, drag proportional to the
  squared angle, ``remus100.m`` 238-245) or ``"linear_section"`` (lift slope,
  constant zero-lift drag, Prestero 2001, eq. 4.37, p. 31); ``"lifting_line"``
  waits on its source and refuses.

The force is resolved along fixed BODY axes, the small-angle form of the
cited sources, not a large-angle flow-axis model: f = 1/2 rho A U^2 (C_L
lift_axis - C_D chord_axis), tau = [f; r x f] (Prestero 2001, eqs. 4.37,
4.41-4.43, pp. 31-33; Fossen 2011, eq. 12.226, p. 400).

Defaults: one fin at the fin post, a rudder lifting along +y (the sign of
the lift axis is a convention of the composition); provenance of every
default number is in ``DEFAULTS.md``.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics
    and Motion Control. Wiley. Eq. 12.226, p. 400.
[Prestero 2001] Prestero, T. (2001). Verification of a six-degree of freedom
    simulation model for the REMUS autonomous underwater vehicle. MIT/WHOI
    MSc thesis. Eqs. 4.37, 4.40-4.43, pp. 31-33; Table A.5, p. 103.
[Murray-Smith 2016] Murray-Smith, D. J. (2016). Inverse simulation methods
    applied to investigations of actuator nonlinearities in ship steering.
    Simulation Notes Europe 26(4), 245-254. Fig. 1 and text, p. 246.
[Sarhadi 2026] Sarhadi, P. (2026). Simple yet effective anti-windup
    techniques for amplitude and rate saturation: an AUV case study.
    arXiv:2601.01302v2, Fig. 4, p. 4.
[MSS] Fossen, T. I. MSS, MIT, CRAFT/AUV/models/remus100.m 109, 113-114,
    234-254 @ cc07579.

Author:    Enio Krizman
Date:      2026-10-09
"""
import casadi as ca
from rpp_plugin_types.more_dynamics import ForceProducer
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_transformations.more_casadi_transformations import freeze
from more_dynamics.models.fin.flow_angle import none as flow_angle_none, small_angle
from more_dynamics.models.fin.inflow import rigid_point, translational
from more_dynamics.models.fin.interference import none as interference_none
from more_dynamics.models.fin.lifting_fin import (check_lifting_fin_values, lifting_fin_casadi,
                                                  lifting_fin_parameters)
from more_dynamics.models.fin.section import linear_section, quadratic_drag
from more_dynamics.models.servo.servo import servo_casadi, servo_parameters
from more_dynamics.plugins.shared.payload_io import PayloadBuilder, frozen_block

FED_INPUTS = ["command", "state", "nu_r", "water_density"]
INFLOWS = {"translational": (translational.translational_casadi, translational.translational_parameters),
           "rigid_point": (rigid_point.rigid_point_casadi, rigid_point.rigid_point_parameters)}
FLOW_ANGLES = {"none": (flow_angle_none.none_casadi, flow_angle_none.none_parameters),
               "small_angle": (small_angle.small_angle_casadi, small_angle.small_angle_parameters)}
INTERFERENCES = {"none": (interference_none.none_casadi, interference_none.none_parameters)}
SECTIONS = {"quadratic_drag": (quadratic_drag.quadratic_drag_casadi, quadratic_drag.quadratic_drag_parameters),
            "linear_section": (linear_section.linear_section_casadi, linear_section.linear_section_parameters)}
WAITING = {("interference", "slender_body"): "Pitts 1957 slender-body interference factors",
           ("section", "lifting_line"): "lifting-line induced drag"}


def _form(slot, name, forms):
    if (slot, name) in WAITING:
        raise NotImplementedError(f"Fin: {slot} {name!r} waits on source: {WAITING[(slot, name)]}")
    if name not in forms:
        raise ValueError(f"Fin: {slot} must be one of {sorted(forms) + [w[1] for w in WAITING if w[0] == slot]}, "
                         f"got {name!r}")
    return forms[name]


class Fin(ForceProducer):
    PARAMETERS = [
        ParameterDescription("fin_position", [-0.638, 0.0, 0.0]),  # m, Prestero 2001, Table A.5, p. 103 (x_finpost)
        ParameterDescription("chord_axis", [1.0, 0.0, 0.0]),  # unit chord direction, BODY FRD (convention)
        ParameterDescription("lift_axis", [0.0, 1.0, 0.0]),  # unit positive-lift direction, BODY FRD (convention)
        ParameterDescription("fin_area", 0.00665),  # m^2, Prestero 2001, Table A.5, p. 103 (S_fin)
        ParameterDescription("servo_dynamics", "first_order_lag"),  # Murray-Smith 2016, Fig. 1, p. 246
        ParameterDescription("servo_rate_limit", True),
        ParameterDescription("servo_angle_limit", "on_command"),
        ParameterDescription("time_constant", 0.1),  # s, Sarhadi 2026, Fig. 4, p. 4
        ParameterDescription("max_deflection", 0.3490658503988659),  # rad = 20 deg, Sarhadi 2026, Fig. 4, p. 4
        ParameterDescription("max_rate", 0.5235987755982988),  # rad/s = 30 deg/s, Sarhadi 2026, Fig. 4, p. 4
        ParameterDescription("inflow", "rigid_point"),  # Prestero 2001, eq. 4.40, p. 32
        ParameterDescription("flow_angle", "small_angle"),  # Prestero 2001, eqs. 4.42-4.43, pp. 32-33
        ParameterDescription("interference", "none"),  # Prestero 2001, eq. 4.41, p. 32 (both factors one)
        ParameterDescription("section", "quadratic_drag"),  # remus100.m 238-245
        ParameterDescription("lift_slope", 3.12),  # 1/rad, Prestero 2001, Table A.5, p. 103 (c_L_alpha)
        ParameterDescription("zero_lift_drag", 0.0),  # Prestero 2001, eq. 4.37, p. 31 carries no fin drag (read by "linear_section")
    ]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        selectors = {name: context.get_parameter(name)
                     for name in ("servo_dynamics", "servo_rate_limit", "servo_angle_limit")}
        servo = {"dynamics": selectors["servo_dynamics"], "rate_limit": selectors["servo_rate_limit"],
                 "angle_limit": selectors["servo_angle_limit"]}
        parts = {"servo": frozen_block(context, servo_casadi(**servo), servo_parameters(**servo))}
        for slot, forms in (("inflow", INFLOWS), ("flow_angle", FLOW_ANGLES), ("interference", INTERFERENCES),
                            ("section", SECTIONS)):
            casadi_function, parameters = _form(slot, context.get_parameter(slot), forms)
            parts[slot] = frozen_block(context, casadi_function(), parameters())
        declared = lifting_fin_parameters()
        values = check_lifting_fin_values({d.name: context.get_parameter(d.name) for d in declared})
        model = freeze(lifting_fin_casadi(parts), declared, values)
        left_open = [name for name in model.name_in() if name not in FED_INPUTS]
        if left_open:
            raise ValueError(f"Fin: parameters {left_open} of its parts are open; "
                             "a fin takes every number from its parameters")
        self._model = model

    def graph(self) -> ForceProducer.CasadyPayload:
        if self._model is None:
            raise RuntimeError("Fin must be initialized before graph()")
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
