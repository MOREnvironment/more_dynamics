import casadi as ca
import rclpy
from rclpy.node import Node
from nav_msgs.msg import Odometry
from rosgraph_msgs.msg import Clock
from std_msgs.msg import Float64MultiArray
from rpp_plugin_types.more_dynamics import VehicleModel3D
from rpp_py.context_builder import ComponentContextBuilder



class SimulationRos(Node):

    COMPONENTS = {
        "vessel": "more_dynamics::VehicleModel3D"
    }

    def __init__(self):
        super().__init__("simulation_node")
        self.context_builder = ComponentContextBuilder()
        self.rpp_context = self.context_builder.build_from_script(__file__)
        self.rpp_context.initialize()
        self.vessel: VehicleModel3D = self.rpp_context.get_component("vessel")

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