class MolecularDynamics:

    def standard(self, simulation, params):
        """

        :param simulation:
        :return:
        """
        n_steps = int(params.get("simulation time (ns)")*1e6\
                  / params.get("timestep (fs)"))
        simulation.step(n_steps)

