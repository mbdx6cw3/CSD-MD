from openmm import unit
import numpy as np
import tensorflow as tf
from network import Network
import os

class MolecularDynamics:

    def standard(self, md_params, simulation, force):
        """

        :param simulation:
        :return:
        """

        dt = md_params.get("timestep (fs)") * unit.femtoseconds
        time = md_params.get("simulation time (ns)") * unit.nanoseconds
        n_steps = int(time / dt)

        if md_params.get("pair-net"):
            model, atoms = load_pairnet()
            n_atoms = len(atoms)
            for i in range(n_steps):
                coords = simulation.context.getState(getPositions=True). \
                    getPositions(asNumpy=True).in_units_of(unit.angstrom)
                if (i % 1000) == 0:
                    tf.keras.backend.clear_session()
                prediction = model.predict_on_batch([np.reshape(coords[:n_atoms]
                    / unit.angstrom, (1, -1, 3)), np.reshape(atoms, (1, -1))])
                forces = prediction[0]*unit.kilocalories_per_mole/unit.angstrom
                forces = np.reshape(forces, (-1, 3))
                # TODO: reorder forces here???
                for j in range(n_atoms):
                    force.setParticleParameters(j, j, forces[j])
                force.updateParametersInContext(simulation.context)
                # TODO: predict potential energy???
                simulation.step(1)
        else:
            simulation.step(n_steps)
        return None


def load_pairnet():
    input_dir = "trained_model"
    isExist = os.path.exists(input_dir)
    if not isExist:
        print("ERROR. Previously trained model could not be located.")
        exit()

    print("Loading a previously trained model...")
    prescale = np.loadtxt(f"./{input_dir}/prescale.txt", dtype=np.float32).reshape(-1)
    atoms = np.loadtxt(f"./{input_dir}/atoms.txt", dtype=np.float32).reshape(-1)
    ann_params = Network().read_params(f"./{input_dir}/ann_params.txt")
    model = Network().build(len(atoms), ann_params, prescale)
    model.summary()
    model.load_weights(f"./{input_dir}/best_ever_model")
    return model, atoms