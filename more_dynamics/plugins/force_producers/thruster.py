from typing import List
from rpp_plugin_types.more_dynamics import ForceProducer
from rpp_py.context import ComponentContext
from rpp_py.lifecycle import CallbackResult


class Thruster(ForceProducer):

    COMPONENTS = {
        "sensors" : "List[more_sensors::Sensor]",
        "actuators" : "List[more_dynamics::ForceProducer]"
    }

    PARAMETERS = [
    ]

    def __init__(self):
        pass

    def on_configure(self, context: ComponentContext) -> CallbackResult:
        pass
        return CallbackResult()

    def on_cleanup(self) -> CallbackResult:
        self.__init__()
        return CallbackResult()
