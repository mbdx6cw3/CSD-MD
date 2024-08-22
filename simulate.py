from openmm import unit

class MolecularDynamics:

    def standard(self, md_params, system, simulation):
        """

        :param simulation:
        :return:
        """

        dt = md_params.get("timestep (fs)") * unit.femtoseconds
        time = md_params.get("simulation time (ns)") * unit.nanoseconds
        n_steps = time / dt

        if md_params.get("pair-net"):
            pass
        else:
            simulation.step(n_steps)
        return None
