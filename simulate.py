from openmm import unit
import numpy as np
from network import Network
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import tensorflow as tf

class MolecularDynamics:

    def standard(self, md_params, simulation, ml_force):
        """

        :param simulation:
        :return:
        """
        tf.get_logger().setLevel('ERROR')

        dt = md_params.get("timestep (fs)") * unit.femtoseconds
        time = md_params.get("simulation time (ns)") * unit.nanoseconds
        n_steps = int(time / dt)

        if md_params.get("pair-net"):
            model, atoms = load_pairnet()
            n_atoms = len(atoms)
            # this step is necessary because pairnet predicts forces for an arbitrary atom ordering
            # TODO: Ideally would not need a manually created mapping file.
            # TODO: For now, can remap atoms in this way because we only have models for a few molecules anyway.
            # TODO: This can be made more efficient if either:
            # TODO: 1) atom ordering of CCDC outputs can be changed
            # TODO: 2) mapping of topology/structure prior to starting simulation
            mapping = np.loadtxt(f"./atom_mapping.dat", dtype=int)
            csd2pairnet = mapping[:, 0]
            pairnet2csd = mapping[:, 1]
            for i in range(n_steps):
                coords = simulation.context.getState(getPositions=True). \
                    getPositions(asNumpy=True).value_in_unit(unit.angstrom)
                # map CSD to pairnet atom order
                coords = coords[csd2pairnet]
                if (i % 1000) == 0:
                    tf.keras.backend.clear_session()
                prediction = model.predict_on_batch([np.reshape(coords[:n_atoms]
                    / unit.angstrom, (1, -1, 3)), np.reshape(atoms, (1, -1))])
                ML_forces = prediction[0]*unit.kilocalories_per_mole/unit.angstrom
                ML_forces = np.reshape(ML_forces, (-1, 3))
                # map pairnet back to CSD
                ML_forces = ML_forces[pairnet2csd]
                #print(ML_forces)
                for j in range(n_atoms):
                    ml_force.setParticleParameters(j, j, ML_forces[j])
                ml_force.updateParametersInContext(simulation.context)

                # get total forces
                forces = simulation.context.getState(getForces=True). \
                    getForces(asNumpy=True).in_units_of(
                    unit.kilocalories_per_mole / unit.angstrom)

                # check MM contribution to forces (should be 0 for pure ML simulation)
                if md_params.get("pair-net"):
                    MM_forces = forces[:n_atoms] - ML_forces

                # TODO: these ML forces appear to be correct and are very similar to those
                # TODO: predicted by PairNetOps-2 using the same structure.
                # TODO: however, the ligand blows up over the first few hundred timesteps
                # TODO: not clear why the potential energy is 600 kJ/mol when forces
                # TODO: predicted from Pair-Net. It should be zero?
                # TODO: predict potential energy as well???
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
    atoms = np.loadtxt(f"./{input_dir}/atoms.txt", dtype=np.float32).reshape(-1)
    ann_params = Network.read_params(f"{input_dir}/ann_params.txt")
    prescale = np.loadtxt(f"./{input_dir}/prescale.txt", dtype=np.float64).reshape(-1)
    network = Network()
    model = network.build(len(atoms), ann_params, prescale)
    model.summary()
    model.load_weights(f"./{input_dir}/best_ever_model").expect_partial()
    return model, atoms

