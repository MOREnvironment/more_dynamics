@0xb94c85a190380959;

using Anot = import "rpp_common/anot.capnp";
using Vm = import "more_dynamics/vehicle_model.capnp";
using Msgs = import "rpp_common/msgs.capnp";

# Parts of a MarineCraft6DOF vehicle. Each is a stateless graph in the shape of
# HydrostaticsModel and HydrodynamicsModel (vehicle_model.capnp): its inputs and
# outputs are named in the payload descriptions; the vehicle connects an input
# to the vehicle quantity or the earlier part output of the same name and size.
# Vehicle quantities: pose (6, NED position and ZYX Euler angles), velocity
# (6, BODY FRD), relative_velocity (6, velocity minus the current, BODY),
# mass_matrix (36). A matrix travels as its column-major vec (6x6 -> 36).
# SI units throughout. Restoring parts use HydrostaticsModel, hull loads
# HydrodynamicsModel, actuators ForceProducer (vehicle_model.capnp).

# Link-only. The registrator links only the files a signature or struct field names, and vehicle_model's
# generated code needs rpp_common/msgs.capnp.c++ (Odometry3D, Command); -Wl,-z,defs then fails the link.
interface LinkOnlyMsgs {
    link @0 (odometry :Msgs.Odometry3D, command :Msgs.Command) -> ();
}

# () -> water_density (1), gravity (1), kinematic_viscosity (1)
interface SiteModel $Anot.plugin("SiteModel") {
    graph @0 () -> (graph :Vm.CasadyPayload);
}

# pose, velocity -> current_velocity (6), current_acceleration (6), BODY
interface CurrentModel $Anot.plugin("CurrentModel") {
    graph @0 () -> (graph :Vm.CasadyPayload);
}

# () -> length (1) and the form's dimensions (beam, draft, span, section_beam, ...)
interface HullForm $Anot.plugin("HullForm") {
    graph @0 () -> (graph :Vm.CasadyPayload);
}

# velocity, dimensions -> rigid_body_mass_matrix (36), rigid_body_coriolis_matrix (36),
# mass (1), center_of_gravity (3), inertia (9)
interface RigidBodyModel $Anot.plugin("RigidBodyModel") {
    graph @0 () -> (graph :Vm.CasadyPayload);
}

# dimensions, water_density, mass properties -> added_mass_matrix (36)
interface AddedMassModel $Anot.plugin("AddedMassModel") {
    graph @0 () -> (graph :Vm.CasadyPayload);
}

# added_mass_matrix, relative_velocity -> added_mass_coriolis_matrix (36)
interface AddedMassCoriolisModel $Anot.plugin("AddedMassCoriolisModel") {
    graph @0 () -> (graph :Vm.CasadyPayload);
}

# The 2-D section drag law of a cross-flow strip: -> section_drag_coefficient (1)
interface SectionDragModel $Anot.plugin("SectionDragModel") {
    graph @0 () -> (graph :Vm.CasadyPayload);
}

# Parts of a LiftingFin (a ForceProducer), same shape. The servo may declare states.
interface ActuatorServo $Anot.plugin("ActuatorServo") {      # command -> deflection
    graph @0 () -> (graph :Vm.CasadyPayload);
}
interface FinInflow $Anot.plugin("FinInflow") {              # relative_velocity, fin_position -> fin_velocity
    graph @0 () -> (graph :Vm.CasadyPayload);
}
interface FinFlowAngle $Anot.plugin("FinFlowAngle") {        # fin_velocity, chord_axis, lift_axis -> flow_angle, speed_squared
    graph @0 () -> (graph :Vm.CasadyPayload);
}
interface FinInterference $Anot.plugin("FinInterference") {  # -> deflection_factor, flow_angle_factor
    graph @0 () -> (graph :Vm.CasadyPayload);
}
interface FinSection $Anot.plugin("FinSection") {            # angle_of_attack -> lift_coefficient, drag_coefficient
    graph @0 () -> (graph :Vm.CasadyPayload);
}
