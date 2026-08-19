import casadi as ca
import numpy as np
from rpp_plugin_types.more_dynamics import VehicleModel3D
from more_dynamics.plugins.casadi_helpers import RppCasadiGraph
from rpp_py.context_builder import ComponentContextBuilder



class Simulation:

    COMPONENTS = {
        "vessel": "more_dynamics::VehicleModel3D"
    }

    def __init__(self):
        self.context_builder = ComponentContextBuilder()
        self.rpp_context = self.context_builder.build_from_script(__file__)
        self.rpp_context.initialize()
        self.vessel: VehicleModel3D = self.rpp_context.get_component("vessel")


class SimulationResult:

    def __init__(self, time, states, outputs):
        self.time = time
        self.states = states
        self.outputs = outputs

def simulate(graph: RppCasadiGraph, inputs, delta_t: float, num_steps: int):
    ic = extract_initial_conditions_from_graph(graph.payload)
    opts = {
        'abstol': 1e-8,
        'reltol': 1e-6,
    }
    x_sim = graph.step.sx_in(0)
    p_sim = graph.step.sx_in(1)
    ode_dict = {'x': x_sim,
                'p': p_sim,
                'ode': graph.step(x_sim, p_sim)
    }
    simulator = ca.integrator('sim', 'cvodes', ode_dict, 0, delta_t, opts)
    state = ca.DM(ic)  # Initial state vector
    time = np.zeros(num_steps+1)
    states = np.zeros((num_steps+1, ic.shape[0]))
    outputs = np.zeros((num_steps+1, graph.num_outputs))
    time[0] = 0.0
    states[0, :] = state.full().flatten()
    outputs[0, :] = graph.output(state, inputs[0, :]).full().flatten()
    for step in range(num_steps):
        # Here you would typically compute the control inputs based on your control strategy.
        # For this example, we'll just use zero inputs.
        u = inputs[step, :]
        # Integrate the system
        result = simulator(x0=state, p=u)
        state = result["xf"].full().flatten()  # Update state for the next step
        output = graph.output(state, u).full().flatten()  # Get the output for the current state and input
        states[step+1, :] = state
        outputs[step+1, :] = output
        time[step+1] = (step+1) * delta_t
    return SimulationResult(time, states, outputs)


def extract_initial_conditions_from_graph(graph: VehicleModel3D.CasadyPayload):
    initial_conditions = []
    for state_desc in graph.stateDescription:
        if state_desc.ic:
            initial_conditions.extend(state_desc.ic)
        else:
            initial_conditions.extend([0.0] * state_desc.size)
    return ca.DM(initial_conditions)

def main():
    simulation = Simulation()

    delta_t = 0.1  # Time step for the simulation
    t_sim = 35.0

    num_steps = int(t_sim / delta_t)
    vessel_graph = RppCasadiGraph(simulation.vessel.graph())
    inputs = ca.DM.zeros(num_steps, vessel_graph.num_inputs)
    simulate(vessel_graph, inputs, delta_t, num_steps)







if __name__ == "__main__":
    main()