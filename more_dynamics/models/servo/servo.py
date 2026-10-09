"""Fin and propulsor-steering servo: one physical device with switches for
which features are present, not one file per form.

A servo has at most one lag state (``dynamics``), may or may not rate-limit
its approach to the target (``rate_limit``, meaningless and refused without
the lag's state), and places its amplitude limit, if any, either on the
command before the lag or on the lag's own output after it
(``angle_limit``). These are features of one device switched on or off, not
different physics (rule 17): a quadratic-drag section and a lifting-line
section are separate files because they are different lift/drag laws; a fin
servo with its lag turned off is the same law with a state removed.

* ``dynamics="none"``: the output is the (optionally limited) command
  directly, no state. This is MSS's own fin servo (``remus100.m`` 109,
  113-114): instantaneous saturation, no dynamics at all.
* ``dynamics="first_order_lag"``, ``angle_limit="on_command"``: the servo's
  target is bounded before the lag (Murray-Smith 2016, Fig. 1 and text,
  p. 246: "If the required rudder deflection delta_c is the variable
  subjected to limiting, the gain factor G_a is unity" — the command enters
  the limiter unscaled). This is the default fin and steering servo.
* ``dynamics="first_order_lag"``, ``angle_limit="on_output"``: the limit
  sits on the integrator's own state instead (Sarhadi 2026, Fig. 4, p. 4:
  lag -> rate limit -> integrator -> amplitude limit on the output). The
  integrator's state is not itself bounded, so a command held beyond the
  limit winds the state up past the stop; after a reversal the output
  leaves the stop late. Equals the command-side form only while the command
  never exceeds the limit and the state starts inside it. Kept as a
  comparison part beside the default, not a higher-fidelity route of its
  own.

A rate limit needs the lag's own state (there is nothing to rate-limit
without one) and ``angle_limit="on_output"`` needs it too (there is no
state distinct from the output to wind up without one); both selector
combinations are refused, naming the selector.

References
----------
[MSS] Fossen, T. I. (2026). Marine Systems Simulator, MIT licence,
    CRAFT/AUV/models/remus100.m 109, 113-114 @ cc07579.
[Murray-Smith 2016] Murray-Smith, D. J. (2016). Inverse simulation methods
    applied to investigations of actuator nonlinearities in ship steering.
    Simulation Notes Europe 26(4), 245-254. Fig. 1 and text, p. 246.
[Sarhadi 2026] Sarhadi, P. (2026). Simple yet effective anti-windup
    techniques for amplitude and rate saturation: an AUV case study.
    arXiv:2601.01302v2, Fig. 4, p. 4.

Author:    Enio Krizman
Date:      2026-10-08
"""
import casadi as ca
from more_transformations.more_casadi_transformations import Parameter, symbols
from more_dynamics.models.shared.force_producer_common import saturate

DYNAMICS = ("none", "first_order_lag")
ANGLE_LIMIT = ("none", "on_command", "on_output")


def _check_selectors(dynamics, rate_limit, angle_limit):
    if dynamics not in DYNAMICS:
        raise ValueError(f"dynamics: must be one of {DYNAMICS}, got {dynamics!r}")
    if angle_limit not in ANGLE_LIMIT:
        raise ValueError(f"angle_limit: must be one of {ANGLE_LIMIT}, got {angle_limit!r}")
    if rate_limit not in (False, True):
        raise ValueError(f"rate_limit: must be False or True, got {rate_limit!r}")
    if rate_limit and dynamics == "none":
        raise ValueError("rate_limit: a rate limit needs the lag's state (dynamics='first_order_lag')")
    if angle_limit == "on_output" and dynamics == "none":
        raise ValueError("angle_limit: 'on_output' needs the lag's own state (dynamics='first_order_lag'); "
                          "with no dynamics the output is the command and 'on_command' already saturates it")


def _angle_limit_meaning(angle_limit):
    if angle_limit == "on_command":
        return ("maximum commanded angle, limiting the target before the lag "
                "(MSS remus100.m 113-114 with no lag; Murray-Smith 2016 Fig. 1 p. 246 with the lag)")
    return "maximum output angle, limiting the integrator's own state after the lag (Sarhadi 2026 Fig. 4 p. 4)"


def servo_parameters(*, dynamics, rate_limit, angle_limit):
    _check_selectors(dynamics, rate_limit, angle_limit)
    declared = []
    if angle_limit != "none":
        declared.append(Parameter("max_deflection", (1, 1), "rad",
                                   _angle_limit_meaning(angle_limit), 0.0, minimum_exclusive=True))
    if rate_limit:
        declared.append(Parameter("max_rate", (1, 1), "rad/s",
                                   "maximum angular rate (Murray-Smith 2016 Fig. 1 p. 246; Sarhadi 2026 Fig. 4 p. 4)",
                                   0.0, minimum_exclusive=True))
    if dynamics == "first_order_lag":
        declared.append(Parameter("time_constant", (1, 1), "s",
                                   "first-order lag time constant (Murray-Smith 2016 Fig. 1 p. 246)",
                                   0.0, minimum_exclusive=True))
    return tuple(declared)


def servo_casadi(*, dynamics, rate_limit, angle_limit):
    _check_selectors(dynamics, rate_limit, angle_limit)
    declared = servo_parameters(dynamics=dynamics, rate_limit=rate_limit, angle_limit=angle_limit)
    p = symbols(declared)
    names = [d.name for d in declared]
    command = ca.SX.sym("command")

    if dynamics == "none":
        state = ca.SX.sym("servo_state", 0)
        if angle_limit == "on_command":
            delta = saturate(command, -p["max_deflection"], p["max_deflection"])  # (MSS remus100.m 113-114)
        else:
            delta = command
        return ca.Function("servo", [command, state, *[p[n] for n in names]], [delta, ca.SX.zeros(0, 1)],
                           ["command", "servo_state", *names], ["deflection", "servo_state_dot"])

    state = ca.SX.sym("servo_state")
    if angle_limit == "on_command":
        target = saturate(command, -p["max_deflection"], p["max_deflection"])  # (Murray-Smith 2016, Fig. 1, p. 246)
    else:
        target = command  # (Sarhadi 2026, Fig. 4, p. 4: the limit is on the output, not on the command)
    rate = (target - state) / p["time_constant"]
    if rate_limit:
        rate = saturate(rate, -p["max_rate"], p["max_rate"])
    if angle_limit == "on_output":
        delta = saturate(state, -p["max_deflection"], p["max_deflection"])  # (Sarhadi 2026, Fig. 4, p. 4)
    else:
        delta = state
    return ca.Function("servo", [command, state, *[p[n] for n in names]], [delta, rate],
                       ["command", "servo_state", *names], ["deflection", "servo_state_dot"])
