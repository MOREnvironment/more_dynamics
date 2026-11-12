import casadi as ca
import rclpy
from rclpy.node import Node
from nav_msgs.msg import Odometry
from rosgraph_msgs.msg import Clock
from std_msgs.msg import Float64MultiArray
from rpp_py.context_builder import ComponentRuntimeBuilder
from rpp_py.data_manager import DataManager

_DATA_MANAGER = DataManager()

from rpp_plugin_types.more_dynamics import VehicleModel3D  # noqa: E402


class SimulationRos(Node):

    COMPONENTS = {
        "vessel": "more_dynamics::VehicleModel3D"
    }

    def __init__(self):
        super().__init__("simulation_node")
        self.context_builder = ComponentRuntimeBuilder(
            data_manager=_DATA_MANAGER
        )
        self.rpp_context = self.context_builder.build_script_from_path(__file__)
        self.rpp_context.configure().require_success("configure simulation")
        self.rpp_context.activate().require_success("activate simulation")
        context = self.rpp_context.get_context()
        context.activate_subcomponents().require_success("activate vessel")
        self.vessel: VehicleModel3D = context.get_subcomponent_instance("vessel")

        self.clock_pub = self.create_publisher(Clock, "/clock", 10)
        self.odom_pub = self.create_publisher(Odometry, "/odom", 10)


        self.delta_t = self.declare_parameter("delta_t", 0.1).value

        self.clock_msg = Clock()
        self.time = 0.0

        self.cmd_sub = self.create_subscription(
            Float64MultiArray,
            "/vessel/cmd",
            self.cmd_callback,
            10
        )

    def cmd_callback(self, msg: Float64MultiArray):
        # Update simulation

        # Update clock
        self.time += self.delta_t
        self.publish_clock()

    def publish_clock(self):
        self.clock_msg.clock.sec = int(self.time)
        self.clock_msg.clock.nanosec = int((self.time - int(self.time)) * 1e9)
        self.clock_pub.publish(self.clock_msg)


def main():
    rclpy.init()
    simulation = Simulation()
    rclpy.spin(simulation)
    rclpy.shutdown()







if __name__ == "__main__":
    main()
