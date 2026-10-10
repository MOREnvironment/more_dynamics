@0x932ab4a8c784ba9c;

using Anot = import "rpp_common/anot.capnp";
using Vm = import "more_dynamics/vehicle_model.capnp";
using Msgs = import "rpp_common/msgs.capnp";

# The water and the gravity a vehicle floats in, in the shape of HydrostaticsModel and HydrodynamicsModel
# (vehicle_model.capnp): a stateless graph with no input. Its outputs are named in the payload descriptions:
# water_density (1, kg/m^3), gravity (1, m/s^2) and kinematic_viscosity (1, m^2/s), each of size 1; a vehicle
# connects every input of its other children to the output of the same name and size. One Site component linked
# into the `site` slot of several vehicles puts them in the same water; each vehicle may hold its own.
#
# Link-only. The registrator links only the files a signature or struct field names, and vehicle_model's
# generated code needs rpp_common/msgs.capnp.c++ (Odometry3D, Command); -Wl,-z,defs then fails the link.
interface LinkOnlyMsgs {
    link @0 (odometry :Msgs.Odometry3D, command :Msgs.Command) -> ();
}

# () -> water_density (1), gravity (1), kinematic_viscosity (1)
interface SiteModel $Anot.plugin("SiteModel") {
    graph @0 () -> (graph :Vm.CasadyPayload);
}
