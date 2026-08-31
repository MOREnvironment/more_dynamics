from typing import List
from rpp_plugin_types.more_dynamics import ForceProducer
from rpp_py.context import ComponentContext


class Thruster(ForceProducer):

    COMPONENTS = {
        "sensors" : "List[more_sensors::Sensor]",
        "actuators" : "List[more_dynamics::ForceProducer]"
    }

    PARAMETERS = [
    ]

    def __init__(self):
        pass

    def initialize(self, context: ComponentContext):
        pass
