@0xaaaabbbbccccdddd;

using Msgs = import "rpp_common/msgs.capnp";
using Anot = import "rpp_common/anot.capnp";

interface ForceProducer $Anot.plugin("ForceProducer") {
    step @0 (state :Msgs.Odometry3D, command :Msgs.Command, t :Float64, dt :Float64) -> (force :Msgs.Wrench3D);
    graph @1 () -> (graph : CasadyPayload);
}

# A stateless hydrostatics graph maps vessel pose to a signed generalized
# force, displaced volume, and wetted surface area. Vehicle dynamics add the
# returned generalized force to the other applied forces.
interface HydrostaticsModel $Anot.plugin("HydrostaticsModel") {
    graph @0 () -> (graph :CasadyPayload);
}

# A stateless hydrodynamics graph maps vessel velocity to a signed generalized
# force. Vehicle dynamics add the returned force to the other applied forces.
interface HydrodynamicsModel $Anot.plugin("HydrodynamicsModel") {
    graph @0 () -> (graph :CasadyPayload);
}

interface DynamicsModel(State, CommandType) {
    step @0 (state :State, command :CommandType, t :Float64, dt :Float64) -> (new_state :State);
    graph @1 () -> (graph :CasadyPayload);
}

interface VehicleModel3D extends(DynamicsModel(Msgs.Odometry3D, List(Msgs.Command)))
$Anot.plugin("VehicleModel3D") {
    getInputDescriptions @0 () -> (inputs :IODescription);
}

struct IODescription {
    size @0 :UInt32;
    min @1 :List(Float64);
    max @2 :List(Float64);
    name @3 :Text;
    description @4 :Text;
}

struct StateDescription {
    size @0 :UInt32;
    min @1 :List(Float64);
    max @2 :List(Float64);
    ic @3 :List(Float64);
    name @4 :Text;
    description @5 :Text;
}

struct CasadyPayload {
    inputDescription @0 :List(IODescription);
    outputDescription @1 :List(IODescription);
    stateDescription @2 :List(StateDescription);
    # Empty when the component has no state dynamics.
    dynamics @3 :Data;
    output @4 :Data;
}
