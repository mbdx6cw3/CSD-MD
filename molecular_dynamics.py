class MolecularDynamics():

    import warnings
    warnings.filterwarnings("ignore", category=DeprecationWarning)

    def __init__(self):
        import sys
        self.resetopenflags = sys.getdlopenflags()
        print(self.resetopenflags)
        pass


    def read_inputs(self):
        '''

        :returns md_params: contains all molecular dynamics parameters required to
                            run a simulation.
        '''
        import yaml, argparse, os, shutil

        # parse input arguments
        parser = argparse.ArgumentParser()
        parser.add_argument("--md_params", default="input.yaml")
        args = parser.parse_args()
        input_file = args.md_params
        isExist = os.path.exists(input_file)
        if not isExist:
            print("Error: input file does not exist.")
        else:
            with open(input_file, 'r') as input:
                md_params = yaml.full_load(input)

        system_type = md_params.get("system type")
        if system_type != "ligand" and system_type != "protein" and \
            system_type != "ligand-protein":
            print("ERROR - system type not allowed.")
            print("Allowed system types: ligand, protein or ligand-protein")
            exit()

        if system_type == "ligand" or system_type == "ligand-protein":
            self.ligand = True
        else:
            self.ligand = False

        if system_type == "protein" or system_type == "ligand-protein":
            self.protein = True
        else:
            self.protein = False

        self.CSD = md_params.get("CSD identifier")
        self.PDB = md_params.get("PDB identifier")
        self.solvate = md_params.get("solvate system")
        self.temp = md_params.get("temperature (K)")
        self.ensemble = md_params.get("ensemble")
        self.dt = md_params.get("timestep (fs)")
        self.pairnet_path = md_params.get("pair-net model path")
        self.time = md_params.get("simulation time (ns)")
        self.simulation_type = md_params.get("simulation type")
        self.analyse_torsions = md_params.get("analyse torsions")
        self.charges = md_params.get("partial charges")

        self.input_dir = "md_input"
        isExist = os.path.exists(self.input_dir)
        if isExist:
            shutil.rmtree(self.input_dir)
        os.makedirs(self.input_dir)

        self.output_dir = "md_data"
        isExist = os.path.exists(self.output_dir)
        if isExist:
            shutil.rmtree(self.output_dir)
        os.makedirs(self.output_dir)

        return md_params


    def setup(self):
        """Set up an MD simulation.
        ...
        """
        from openmm import LangevinMiddleIntegrator, app, unit
        from openmmforcefields.generators import SystemGenerator
        from openff.toolkit.topology import Molecule
        import warnings, sys
        warnings.filterwarnings("ignore", category=DeprecationWarning)

        # this is required to reset paths after loading ccdc modules
        sys.setdlopenflags(self.resetopenflags)

        # set name of input file for this system
        if self.ligand and not self.protein:
            input_file = f"{self.input_dir}/ligand_0.pdb"

        if self.protein and not self.ligand:
            input_file = f"{self.input_dir}/protein.pdb"

        if self.ligand and self.protein:
            input_file = f"{self.input_dir}/protein-ligand.pdb"

        # retrieve pdb file
        pdb = self.get_pdb(input_file)

        # define force field
        std_ff = ["amber/ff14SB.xml", "amber/tip3p_standard.xml"]
        small_mol_ff = "gaff-2.11"

        # get structure and topology from pdb file
        self.topology = pdb.topology
        self.positions = pdb.positions

        # create topology for non-standard residue from SMILES
        if self.ligand:
            print("Creating topology...")
            ligand = Molecule.from_smiles(self.smiles, allow_undefined_stereo=True)
            self.ligand_n_atom = ligand.n_atoms
            molecules = ligand
            if not self.solvate:
                self.topology.setUnitCellDimensions([3.0]*3)
        else:
            molecules = None

        # add water to system using PDBfixer
        if self.solvate:
            solvated_pdb = self.solvate_system(input_file)
            n_solvent = len(solvated_pdb.positions) - len(pdb.positions)
            self.topology = solvated_pdb.topology
            self.positions = solvated_pdb.positions
            print(f"Adding {int(n_solvent / 3)} solvent molecules...")

        # construct OpenMM system object using topology and forcefield
        system_generator = SystemGenerator(forcefields=std_ff,
            small_molecule_forcefield=small_mol_ff,
            forcefield_kwargs={"rigidWater": True}, cache="db.json")
        system = system_generator.create_system(self.topology, molecules=molecules)

        '''
        if self.simulation_type == "enhanced":
            from openmmplumed import PlumedForce
            plumed_file = open(f"plumed.dat", "r")
            plumed_script = plumed_file.read()
            system.addForce(PlumedForce(plumed_script))
            plumed_file.close()
        '''

        # if using pair-net must set all intramolecular ligand interactions to zero
        if self.ligand:
            if self.pairnet_path != "none":
                print("Using ML force field (PairNet) for ligand...")
                self.create_MLP(system)
            else:
                print("Using MM force field (GAFF2) for ligand...")
                self.ml_force = None
                self.assign_charges(system)

        # setup integrator
        self.temp = self.temp*unit.kelvin
        self.dt = self.dt*unit.femtoseconds
        temp_coupling = 1/unit.picosecond
        if self.ensemble == "NVT":
            integrator = LangevinMiddleIntegrator(self.temp, temp_coupling, self.dt)

        # setup simulation and output
        # Simulation object ties together topology, system and integrator
        # and maintains list of reporter objects that record or analyse data
        self.simulation = app.Simulation(self.topology, system, integrator)
        self.simulation.context.setPositions(self.positions)
        self.simulation.reporters.append(app.PDBReporter("output.pdb", 1000,
            enforcePeriodicBox=True))
        # TODO: replace below with CSDDataReporter?
        '''
        self.simulation.reporters.append(app.StateDataReporter(stdout, 100,
            step=True, potentialEnergy=True, temperature=True))
        '''
        return None


    def simulate(self):
        """

        :param:
        :return:
        """
        from openmm import unit
        import os
        os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
        import tensorflow as tf
        import numpy as np

        tf.get_logger().setLevel('ERROR')

        time = self.time*unit.nanoseconds
        n_steps = int(time/self.dt)+1

        if self.ligand:

            if self.pairnet_path != "none":
                input_dir = f"{self.pairnet_path}trained_model"
                print(f"Using pairnet model: {input_dir}")
                isExist = os.path.exists(input_dir)
                if not isExist:
                    print("ERROR. Previously trained model could not be located.")
                    exit()
                model, atoms = self.load_pairnet(input_dir)
                map = False
                isExist = os.path.exists(f"{self.pairnet_path}atom_mapping.dat")
                if isExist:
                    map = True
                    # this step is necessary because pairnet models have a specific (arbitrary) atom ordering
                    # TODO: Ideally would not need a manually created mapping file.
                    # TODO: For now, can remap atoms in this way because we only have models for a few molecules anyway.
                    # TODO: This can be made more efficient if either:
                    # TODO: 1) atom ordering of CCDC outputs can be changed
                    # TODO: 2) mapping of topology/structure prior to starting simulation
                    mapping = np.loadtxt(f"{self.pairnet_path}/atom_mapping.dat", dtype=int)
                    csd2pairnet = mapping[:, 0]
                    pairnet2csd = mapping[:, 1]

            # open files for storing PairNetOps compatible datasets
            f1 = open(f"{self.output_dir}/coords.txt", "w")
            f2 = open(f"{self.output_dir}/forces.txt", "w")
            f3 = open(f"{self.output_dir}/energies.txt", "w")
            f4 = open(f"{self.output_dir}/charges.txt", "w")

            # open files for storing PairNetOps compatible datasets
            data_files = [f1, f2, f3, f4]

        # loop over conformers
        for i_conf in range(self.n_conf):

            if self.ligand:
                # for first conformer use initial coords and vels defined in setup
                # for all other conformers reset coords/vels including solvent
                print(f"Conformer number: {i_conf+1}")
                if i_conf > 0:
                    input_file = f"{self.input_dir}/ligand_{i_conf}.pdb"
                    pdb = self.get_pdb(input_file)
                    if self.solvate:
                        print("ERROR - cannot yet simulate multiple solvated conformers.")
                        print("Try single solvated conformer or multi-conformer in gas phase instead.")
                        exit()
                    self.simulation.context.setPositions(pdb.positions)
                    self.simulation.context.setVelocitiesToTemperature(self.temp)

            if self.protein:
                self.simulation.context.setPositions(self.positions)
                print("Minimising initial protein structure...")
                self.simulation.minimizeEnergy()

            if self.ensemble == "NVT":
                self.simulation.context.setVelocitiesToTemperature(self.temp)

            print("Performing MD simulation...")
            print("Time (ps) | PE (kcal/mol)")
            for i in range(n_steps):

                if self.ligand:

                    if self.pairnet_path != "none":

                        # this stops us running out of memory
                        if (i % 1000) == 0:
                            tf.keras.backend.clear_session()

                        coords = self.simulation.context.getState(getPositions=True). \
                            getPositions(asNumpy=True).value_in_unit(unit.angstrom)

                        if map:
                            coords = coords[csd2pairnet] # map CSD to pairnet atom order

                        prediction = self.get_pairnet_prediction(model, atoms, coords)
                        ML_forces = prediction[0] * unit.kilocalories_per_mole / unit.angstrom
                        ML_forces = np.reshape(ML_forces, (-1, 3))
                        energy = prediction[2][0][0]

                        if map:
                            ML_forces = ML_forces[pairnet2csd] # map pairnet back to CSD

                        for j in range(self.ligand_n_atom):
                            self.ml_force.setParticleParameters(j, j, ML_forces[j])
                        self.ml_force.updateParametersInContext(self.simulation.context)

                # every 1000 steps save data for PairNetOps compatible dataset
                if self.ligand:
                    if (i % 1000) == 0:

                        coords = self.simulation.context.getState(getPositions=True). \
                            getPositions(asNumpy=True).value_in_unit(unit.angstrom)
                        forces = self.simulation.context.getState(getForces=True). \
                            getForces(asNumpy=True).in_units_of(unit.kilocalories_per_mole / unit.angstrom)

                        if self.pairnet_path == "none":
                            state = self.simulation.context.getState(getEnergy=True)
                            energy = state.getPotentialEnergy() / unit.kilocalories_per_mole
                            charges = self.mm_charges

                        self.write_dataset(coords, forces, energy, charges, data_files)

                # TODO: CSDDataReporter output?
                if ((i+1) % 1000) == 0:
                    if self.ligand and self.pairnet_path != "none":
                        energy = prediction[2][0][0]
                    else:
                        state = self.simulation.context.getState(getEnergy=True)
                        energy = state.getPotentialEnergy() / unit.kilocalories_per_mole
                    time = self.simulation.context.getState().getTime(). \
                        value_in_unit(unit.picoseconds)
                    print(f"{time:9.1f} | {energy:13.2f}")

                # advance trajectory one timestep
                self.simulation.step(1)

        print("MD simulation has completed.")
        if self.ligand:
            f1.close()
            f2.close()
            f3.close()
            f4.close()
        return None


    def load_pairnet(self, input_dir):
        """

        :param input_dir: path to pairnet model
        :returns: model - pairnet model
                  atoms - list of atoms
        """
        from network import Network
        import numpy as np
        print("Loading a previously trained model...")
        atoms = np.loadtxt(f"{input_dir}/atoms.txt", dtype=np.float32).reshape(-1)
        ann_params = Network.read_params(f"{input_dir}/ann_params.txt")
        prescale = np.loadtxt(f"{input_dir}/prescale.txt", dtype=np.float64).reshape(-1)
        network = Network()
        model = network.build(len(atoms), ann_params, prescale)
        model.summary()
        model.load_weights(f"{input_dir}/best_ever_model").expect_partial()
        return model, atoms


    def get_pairnet_prediction(self, model, atoms, coords):
        import numpy as np
        from openmm import unit

        prediction = model.predict_on_batch([np.reshape(coords[:self.ligand_n_atom]
            / unit.angstrom, (1, -1, 3)), np.reshape(atoms, (1, -1))])

        return prediction


    def get_pdb(self, filename):
        from openmm import app
        # TODO: ***CONVERSION FUNCTIONS***
        # TODO: This is where the coordinates are read in.
        pdb = app.PDBFile(filename)
        return pdb


    def solvate_system(self, filename):
        from pdbfixer import PDBFixer
        from openmm import Vec3
        fixer = PDBFixer(filename=filename)
        boxSize = 3.0 * Vec3(1, 1, 1)
        fixer.addSolvent(boxSize=boxSize)
        return fixer


    def create_MLP(self, system):
        from openmm import HarmonicBondForce, HarmonicAngleForce, PeriodicTorsionForce
        from openmm import NonbondedForce, CustomExternalForce

        # exclude all non-bonded interactions
        nb = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]
        for i in range(self.ligand_n_atom):
            for j in range(i):
                nb.addException(i, j, 0, 1, 0, replace=True)

        # set all bond force constants to zero
        bond_force = [f for f in system.getForces() if isinstance(f, HarmonicBondForce)][0]
        for i in range(bond_force.getNumBonds()):
            p1, p2, length, k = bond_force.getBondParameters(i)
            bond_force.setBondParameters(i, p1, p2, length, 0)

        # set all angles force constants to zero
        angle_force = [f for f in system.getForces() if isinstance(f, HarmonicAngleForce)][0]
        for i in range(angle_force.getNumAngles()):
            p1, p2, p3, angle, k = angle_force.getAngleParameters(i)
            angle_force.setAngleParameters(i, p1, p2, p3, angle, 0)

        # set all dihedral force constants/barrier heights to zero
        torsion_force = [f for f in system.getForces() if isinstance(f, PeriodicTorsionForce)][0]
        for i in range(torsion_force.getNumTorsions()):
            p1, p2, p3, p4, n, phase, k = torsion_force.getTorsionParameters(i)
            torsion_force.setTorsionParameters(i, p1, p2, p3, p4, n, phase, 0)

        # create custom force for PairNet predictions
        self.ml_force = CustomExternalForce("-fx*x-fy*y-fz*z")
        system.addForce(self.ml_force)
        self.ml_force.addPerParticleParameter("fx")
        self.ml_force.addPerParticleParameter("fy")
        self.ml_force.addPerParticleParameter("fz")

        for j in range(self.ligand_n_atom):
            self.ml_force.addParticle(j, (0, 0, 0))

        return None


    def assign_charges(self, system):
        from openmm import NonbondedForce, unit
        import numpy as np
        self.mm_charges = np.zeros(self.ligand_n_atom)
        if self.charges == "am1-bcc":
            nb = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]
            for i in range(system.getNumParticles()):
                charge, sigma, epsilon = nb.getParticleParameters(i)
                self.mm_charges[i] = charge.value_in_unit(unit.elementary_charge)
        elif self.charges == "from_file":
            self.mm_charges = np.loadtxt("charges.txt")
        return self.mm_charges


    def write_dataset(self, coords, forces, energy, charges, data_files):
        import numpy as np
        np.savetxt(data_files[0], coords[:self.ligand_n_atom])
        np.savetxt(data_files[1], forces[:self.ligand_n_atom])
        data_files[2].write(f"{energy}\n")
        np.savetxt(data_files[3], charges[:self.ligand_n_atom])
        return None

