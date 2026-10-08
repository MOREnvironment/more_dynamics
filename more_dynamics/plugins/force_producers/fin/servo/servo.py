"""Fin and steering servo as an ``ActuatorServo`` plugin: one device whose lag,
rate limit and angle limit are switched on or off by build-time selectors,
the numbers declared parameters. The command enters as ``command``, the
output is ``deflection``; with the lag there is one state, ``servo_state``.

Selectors (the features of one device, not different physics):

* ``dynamics``: ``"none"`` (the output is the command; MSS's fin servo) or
  ``"first_order_lag"`` (one state; Murray-Smith 2016, Fig. 1, p. 246).
* ``rate_limit``: ``False`` or ``True`` (needs the lag; refused otherwise).
* ``angle_limit``: ``"none"``, ``"on_command"`` (limit before the lag;
  Murray-Smith 2016, Fig. 1, p. 246; MSS ``remus100.m`` 113-114 without the
  lag) or ``"on_output"`` (limit on the lag's own state; Sarhadi 2026, Fig. 4,
  p. 4).

Defaults are the lag with rate and angle limit on the command and the
numbers of the AUV fin servo of Sarhadi (2026): time constant 0.1 s, +-20 deg,
+-30 deg/s (Fig. 4 and text, p. 4; a deliberately demanding choice of that
study, not a measurement). Raising the level: identify the time constant and
the rate limit of the real servo from a logged step in the fin command.

References
----------
[Murray-Smith 2016] Murray-Smith, D. J. (2016). Inverse simulation methods
    applied to investigations of actuator nonlinearities in ship steering.
    Simulation Notes Europe 26(4), 245-254. Fig. 1 and text, p. 246.
[Sarhadi 2026] Sarhadi, P. (2026). Simple yet effective anti-windup
    techniques for amplitude and rate saturation: an AUV case study.
    arXiv:2601.01302v2, Fig. 4, p. 4.
[MSS] Fossen, T. I. MSS, MIT, CRAFT/AUV/models/remus100.m 109, 113-114
    @ cc07579.

Author:    Enio Krizman
Date:      2026-10-08
"""
import casadi as ca
from rpp_plugin_types.more_dynamics import ActuatorServo
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_dynamics.models.force_producers.fin.servo.servo import servo_casadi, servo_parameters
from more_dynamics.plugins.shared.payload_io import OPEN_PARAMETERS, PayloadBuilder, frozen_block


class Servo(ActuatorServo):
    PARAMETERS = [
        ParameterDescription("dynamics", "first_order_lag"),
        ParameterDescription("rate_limit", True),
        ParameterDescription("angle_limit", "on_command"),
        ParameterDescription("time_constant", 0.1),  # s, Sarhadi 2026, Fig. 4, p. 4
        ParameterDescription("max_deflection", 0.3490658503988659),  # rad = 20 deg, Sarhadi 2026, Fig. 4, p. 4
        ParameterDescription("max_rate", 0.5235987755982988),  # rad/s = 30 deg/s, Sarhadi 2026, Fig. 4, p. 4
        OPEN_PARAMETERS,
    ]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        selectors = {name: context.get_parameter(name) for name in ("dynamics", "rate_limit", "angle_limit")}
        self._model = frozen_block(context, servo_casadi(**selectors), servo_parameters(**selectors))

    def graph(self) -> ActuatorServo.CasadyPayload:
        if self._model is None:
            raise RuntimeError("Servo must be initialized before graph()")
        io = PayloadBuilder(ActuatorServo.CasadyPayload())
        n = self._model.size1_in("servo_state")
        state = io.state("servo_state", n, [0.0] * n, "servo lag state, rad") if n else ca.SX.sym("servo_state", 0)
        out = io.call(self._model, given={"servo_state": state})
        io.output("deflection", out["deflection"], "fin deflection, rad")
        if n:
            io.state_dot(out["servo_state_dot"])
        return io.payload()
