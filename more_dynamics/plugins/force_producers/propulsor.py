"""One propulsor as Luka's ``ForceProducer``: six typed sockets (drive, steering, mounting, inflow, law,
housing), each holding one part chosen by a selector, the numbers of the chosen parts declared here as
parameters. The first output is the wrench on the vehicle (``generated_force``, BODY, N and N m about the
CO), then the thrust along the axis and the steering angle; several propulsors on one vehicle are several
entries of the vehicle's own actuator list. The model layer is ``models/propulsor``.

Inputs: ``drive_command`` then ``steering_command`` (a command of size zero is absent: a fixed unit has no
steering command), from outside; ``relative_velocity`` and ``water_density`` are fed by the vehicle by
name. States: ``drive_state`` then ``servo_state`` (each absent when its part has none).

Selectors (the choices are the keys of ``FORMS``): ``drive`` (a shaft-speed command in rpm or rad/s, the
same with a shaft-power cap, a normalised command mapped to a shaft speed, or a normalised effort with a
first-order lag), ``steering`` (fixed, the servo with its three switches set as four named settings, or
the clamped nozzle angle), ``mounting`` (a fixed shaft by Euler angles or by a direction, an azimuth about
the body z axis, a steered rotation, or the direction formula of the shared jet nozzle), ``inflow`` (none,
the speed magnitude, the surge component), ``law`` (the open-water propeller laws, the quadratic law of
two thrusters, the waterjet laws, the normalised linear law), ``housing`` (none, or the normal force of an
outboard leg). The drive and the law must pass the same unit (rev/s for a shaft speed, 1 for a normalised
effort); a pair that does not is refused when the plugin is initialised, naming both.

``steering_convention`` (``"port_positive"`` or ``"starboard_positive"``) is read by the azimuth mounting
only: the name says which way the boat turns for a positive angle when the unit is aft of the CO.

Every number of every part is unset in this plugin and given by the composition's named part; the
selectors default to a fixed-thrust-vector jet (the lag with rate and angle limits for the nozzle, a
unit-length rotation for its direction, a normalised lag and linear law for its thrust); the numbers in
``DEFAULTS.md`` come from the reference parts.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics and Motion Control. Wiley.
    Ch. 12, eq. 12.226, p. 400.
[Mandic] Mandic, L. rpp_control, ``plugins/controllers/cascade_controller.hpp`` and more_dynamics,
    ``plugins/force_producers/fin.py`` (the form: one plugin, parts as selectors).
The part laws and their sources are cited in ``models/propulsor``.

Author:    Enio Krizman
Date:      2026-10-23
"""
import casadi as ca
from rpp_plugin_types.more_dynamics import ForceProducer
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_transformations.more_casadi_transformations import check_values, freeze, unset
from more_dynamics.models.propulsor import propulsor as skeleton
from more_dynamics.models.propulsor.drive import (normalized_lag, power_limited_speed_command, speed_command,
                                                  speed_command_normalized)
from more_dynamics.models.propulsor.housing import leg_normal_force, none as housing_none
from more_dynamics.models.propulsor.inflow import none as inflow_none, speed_magnitude, surge_component
from more_dynamics.models.propulsor.law import (jet_bilinear_speed, jet_bollard_quadratic, jet_momentum_flux,
                                                jet_pump_system_match, normalized_linear, open_water_bollard,
                                                open_water_linearized, open_water_polynomial,
                                                open_water_signed_advance, quadratic_asymmetric)
from more_dynamics.models.propulsor.mounting import (azimuth, fixed_direction, fixed_euler, jet_nozzle_direction,
                                                     steered_rotation)
from more_dynamics.models.propulsor.steering import fixed, jet_nozzle_clamp
from more_dynamics.models.servo import servo
from more_dynamics.plugins.shared.payload_io import PayloadBuilder, read_declared, refuses_unset

FED_INPUTS = ["command", "state", "nu_r", "water_density"]
SERVO_SETTINGS = {  # the three switches of the servo, as four named settings
    "servo_ideal": {"dynamics": "none", "rate_limit": False, "angle_limit": "on_command"},
    "servo_no_limits": {"dynamics": "none", "rate_limit": False, "angle_limit": "none"},
    "servo_lag_rate_angle": {"dynamics": "first_order_lag", "rate_limit": True, "angle_limit": "on_command"},
    "servo_output_saturated": {"dynamics": "first_order_lag", "rate_limit": True, "angle_limit": "on_output"},
}


def _form(module, stem=None, **kwargs):
    return module, stem or module.__name__.rsplit(".", 1)[1], kwargs


# slot -> {selector value: (module, function stem, keywords that build the part)}
FORMS = {
    "drive": {"speed_command_rpm": _form(speed_command, command_unit="rpm"),
              "speed_command_rad_s": _form(speed_command, command_unit="rad/s"),
              "power_limited_speed_command_rpm": _form(power_limited_speed_command, command_unit="rpm"),
              "speed_command_normalized": _form(speed_command_normalized),
              "normalized_lag": _form(normalized_lag)},
    "steering": {"fixed": _form(fixed),
                 **{name: _form(servo, "servo", **settings) for name, settings in SERVO_SETTINGS.items()},
                 "jet_nozzle_clamp": _form(jet_nozzle_clamp)},
    "mounting": {"fixed_euler": _form(fixed_euler), "fixed_direction": _form(fixed_direction),
                 "azimuth": _form(azimuth), "jet_nozzle_direction": _form(jet_nozzle_direction),
                 "steered_rotation": _form(steered_rotation)},
    "inflow": {"none": _form(inflow_none), "speed_magnitude": _form(speed_magnitude),
               "surge_component": _form(surge_component)},
    "law": {"open_water_bollard": _form(open_water_bollard, open_water_coefficients="polynomial"),
            "open_water_linearized": _form(open_water_linearized, open_water_coefficients="polynomial"),
            "open_water_linearized_given": _form(open_water_linearized, open_water_coefficients="given"),
            "open_water_polynomial": _form(open_water_polynomial, open_water_coefficients="polynomial"),
            "open_water_signed_advance": _form(open_water_signed_advance),
            "quadratic_asymmetric": _form(quadratic_asymmetric), "normalized_linear": _form(normalized_linear),
            "jet_bollard_quadratic": _form(jet_bollard_quadratic), "jet_bilinear_speed": _form(jet_bilinear_speed),
            "jet_momentum_flux_given": _form(jet_momentum_flux, flow_coefficient_method="given"),
            "jet_momentum_flux_bollard_point": _form(jet_momentum_flux, flow_coefficient_method="from_bollard_point"),
            "jet_pump_system_match": _form(jet_pump_system_match)},
    "housing": {"none": _form(housing_none), "leg_normal_force": _form(leg_normal_force)},
}


def _declared(slot, name, **extra):
    module, stem, kwargs = FORMS[slot][name]
    return getattr(module, f"{stem}_parameters")(**kwargs, **extra)


def _number_names():
    """Every number of every variant, once, in the order of the sockets."""
    names = []
    for slot in skeleton.SLOTS:
        for name in FORMS[slot]:
            extra = {"steering_convention": "port_positive"} if (slot, name) == ("mounting", "azimuth") else {}
            for d in _declared(slot, name, **extra):
                if d.name not in names:
                    names.append(d.name)
    return names


def _select(slot, value):
    if value not in FORMS[slot]:
        raise ValueError(f"Propulsor: {slot} must be one of {sorted(FORMS[slot])}, got {value!r}")
    return value


class Propulsor(ForceProducer):
    PARAMETERS = [
        ParameterDescription("drive", "normalized_lag"),  # selector: the drive part
        ParameterDescription("steering", "servo_lag_rate_angle"),  # selector: the steering part
        ParameterDescription("mounting", "steered_rotation"),  # selector: the mounting part
        ParameterDescription("inflow", "none"),  # selector: the inflow part
        ParameterDescription("law", "normalized_linear"),  # selector: the law part
        ParameterDescription("housing", "none"),  # selector: the housing part
        ParameterDescription("steering_convention", "port_positive"),  # selector, azimuth mounting only: the boat's turn for a positive angle
        ParameterDescription("shaft_torque_scale", unset()),  # 1: share of the shaft torque that reacts on the body
        ParameterDescription("min_command", unset()),  # rpm or rad/s as the drive says: lowest shaft-speed command
        ParameterDescription("max_command", unset()),  # rpm or rad/s as the drive says: highest shaft-speed command
        ParameterDescription("max_power", unset()),  # W: electrical input power limit P_max
        ParameterDescription("efficiency", unset()),  # 1: drive efficiency eta: the shaft power cap is eta P_max
        ParameterDescription("thrust_saturation", unset()),  # N: thrust limit of the unit
        ParameterDescription("max_shaft_speed", unset()),  # rev/s: shaft speed at a command of +1
        ParameterDescription("min_shaft_speed", unset()),  # rev/s: shaft speed at a command of -1
        ParameterDescription("thrust_rise_time", unset()),  # s: time constant when the effort rises
        ParameterDescription("thrust_fall_time", unset()),  # s: time constant when the effort falls
        ParameterDescription("max_deflection", unset()),  # rad: maximum commanded angle, limiting the target before the lag
        ParameterDescription("max_rate", unset()),  # rad/s: maximum angular rate
        ParameterDescription("time_constant", unset()),  # s: first-order lag time constant
        ParameterDescription("min_angle", unset()),  # rad: lowest nozzle angle
        ParameterDescription("max_angle", unset()),  # rad: highest nozzle angle
        ParameterDescription("nozzle_velocity", unset()),  # rad/s: bound of the angle error that drives the nozzle
        ParameterDescription("position", unset()),  # m: CO to the thrust point, BODY
        ParameterDescription("orientation", unset()),  # rad: shaft orientation
        ParameterDescription("direction", unset()),  # 1: thrust direction, BODY
        ParameterDescription("trim_angle", unset()),  # rad: nozzle trim
        ParameterDescription("yaw_bias", unset()),  # rad: yaw bias of the nozzle
        ParameterDescription("thrust_deduction", unset()),  # 1: thrust deduction t of the thrust (1 - t) T
        ParameterDescription("wake_fraction", unset()),  # 1: wake fraction w of the advance speed (1 - w) U
        ParameterDescription("advance_speed_factor", unset()),  # 1: factor c_Va of the advance speed
        ParameterDescription("diameter", unset()),  # m: propeller diameter D
        ParameterDescription("pitch_diameter_ratio", unset()),  # 1: pitch ratio P/D
        ParameterDescription("blade_area_ratio", unset()),  # 1: expanded blade-area ratio A_E/A_O
        ParameterDescription("blade_count", unset()),  # 1: number of blades z
        ParameterDescription("max_advance_number", unset()),  # 1: advance number at the top speed and top shaft speed
        ParameterDescription("thrust_torque_coefficients", unset()),  # 1: [K_T0, K_Q0, K_Tmax, K_Qmax]
        ParameterDescription("reverse_thrust_factor", unset()),  # 1: factor on K_T for n < 0
        ParameterDescription("reverse_torque_factor", unset()),  # 1: factor on K_Q for n < 0
        ParameterDescription("positive_thrust_coefficient", unset()),  # N*s^2/(rev*m): b2 of T = b1 n^2 + b2 n Varad^2: k_pos: T = k_pos omega |omega| for omega > 0
        ParameterDescription("negative_thrust_coefficient", unset()),  # N*s^2/(rev*m): b2 of T = b1 n^2 + b2 n Varad^2: k_neg: T = k_neg omega |omega| for omega <= 0
        ParameterDescription("max_thrust", unset()),  # N: thrust at an effort of 1
        ParameterDescription("thrust_coefficient", unset()),  # 1: scale of the thrust
        ParameterDescription("bollard_thrust_coefficient", unset()),  # N*s^2/(rev*m): b2 of T = b1 n^2 + b2 n Varev^2: k of T = k n^2 at the bollard
        ParameterDescription("speed_squared_coefficient", unset()),  # N*s^2/(rev*m): b2 of T = b1 n^2 + b2 n Varev^2: b1 of T = b1 n^2 + b2 n Va
        ParameterDescription("speed_advance_coefficient", unset()),  # N*s^2/(rev*m): b2 of T = b1 n^2 + b2 n Va
        ParameterDescription("nozzle_area", unset()),  # m^2: nozzle outlet area A_n
        ParameterDescription("flow_coefficient", unset()),  # m^3/rev: c_Q of Q = c_Q n
        ParameterDescription("bollard_thrust", unset()),  # N: thrust T_b measured at the bollard
        ParameterDescription("bollard_shaft_speed", unset()),  # rev/s: shaft speed n_b of that bollard point
        ParameterDescription("impeller_diameter", unset()),  # m: impeller diameter D
        ParameterDescription("pump_head_coefficient_at_zero_flow", unset()),  # 1: psi_0 of the pump line psi = psi_0 - s phi
        ParameterDescription("pump_head_slope", unset()),  # 1: s of the pump line psi = psi_0 - s phi
        ParameterDescription("nozzle_loss_coefficient", unset()),  # 1: phi_n: nozzle loss as a fraction of the jet head
        ParameterDescription("inlet_loss_coefficient", unset()),  # 1: eps: share of the intake ram head lost
        ParameterDescription("pump_efficiency", unset()),  # 1: eta_p of the shaft torque
        ParameterDescription("leg_area", unset()),  # m^2: leg and fin area A_leg
        ParameterDescription("leg_aspect_ratio", unset()),  # 1: leg aspect ratio AR
    ]

    def __init__(self) -> None:
        self._model = None
        self._sizes = None

    @refuses_unset
    def initialize(self, context: ComponentContext) -> None:
        names = {slot: _select(slot, context.get_parameter(slot)) for slot in skeleton.SLOTS}
        convention = context.get_parameter("steering_convention")
        skeleton.check_units({"drive": FORMS["drive"][names["drive"]][0].DRIVE_OUTPUT_UNIT,
                              "law": FORMS["law"][names["law"]][0].LAW_INPUT_UNIT})
        parts = {}
        for slot in skeleton.SLOTS:
            parts[slot] = self._frozen_part(context, slot, names[slot], convention)
        units = {"drive": FORMS["drive"][names["drive"]][0].DRIVE_OUTPUT_UNIT,
                 "law": FORMS["law"][names["law"]][0].LAW_INPUT_UNIT}
        declared = skeleton.propulsor_parameters()
        values = skeleton.check_propulsor_values(read_declared(context, declared))
        model = freeze(skeleton.propulsor_casadi(parts, units), declared, values)
        left_open = [name for name in model.name_in() if name not in FED_INPUTS]
        if left_open:
            raise ValueError(f"Propulsor: parameters {left_open} of its parts are open; "
                             "a propulsor takes every number from its parameters")
        self._model = model
        self._sizes = {"drive_command": parts["drive"].size1_in("command"),
                       "steering_command": parts["steering"].size1_in("command"),
                       "drive_state": parts["drive"].size1_in("drive_state"),
                       "servo_state": parts["steering"].size1_in("servo_state")}

    @staticmethod
    def _frozen_part(context, slot, name, convention):
        module, stem, kwargs = FORMS[slot][name]
        if (slot, name) == ("mounting", "azimuth"):
            kwargs = {**kwargs, "steering_convention": convention}
        block = getattr(module, f"{stem}_casadi")(**kwargs)
        declared = getattr(module, f"{stem}_parameters")(**kwargs)
        if not declared:
            return block
        values = read_declared(context, declared)
        check = getattr(module, f"check_{stem}_values", None)
        numbers = check(values, **({"steering_convention": convention} if stem == "azimuth" else {})) \
            if check else check_values(declared, values)
        return freeze(block, declared, numbers)

    def graph(self) -> ForceProducer.CasadyPayload:
        if self._model is None:
            raise RuntimeError("Propulsor must be initialized before graph()")
        io = PayloadBuilder(ForceProducer.CasadyPayload())
        sizes = self._sizes

        def command(name, description):
            return io.input(name, sizes[name], description) if sizes[name] else ca.SX.sym(name, 0)

        def state(name, description):
            n = sizes[name]
            return io.state(name, n, [0.0] * n, description) if n else ca.SX.sym(name, 0)

        drive_command = command("drive_command", "command of the drive (shaft speed, or a normalised effort)")
        steering_command = command("steering_command", "command of the steering (angle, or a normalised angle)")
        drive_state = state("drive_state", "state of the drive")
        servo_state = state("servo_state", "state of the steering")
        out = io.call(self._model, given={"command": ca.vertcat(drive_command, steering_command),
                                          "state": ca.vertcat(drive_state, servo_state)})
        io.output("generated_force", out["tau"], "propulsor wrench about the CO, BODY, N and N m")
        io.output("thrust", out["thrust"], "thrust along the axis after deduction, power cap and limit, N")
        io.output("steering_angle", out["steering_angle"], "deflection of the steering, rad")
        if sizes["drive_state"] + sizes["servo_state"]:
            io.state_dot(out["state_dot"])
        return io.payload()
