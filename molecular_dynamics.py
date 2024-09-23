class MolecularDynamics():

    def __init__(self):
        pass


    def read_inputs(self):
        '''

        :returns md_params: contains all molecular dynamics parameters required to
                            run a simulation.
        '''
        import yaml, argparse, os

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

        self.system_type = md_params.get("system type")
        if self.system_type != "ligand" and self.system_type != "protein" and \
            self.system_type != "ligand-protein":
            print("ERROR - system type not allowed.")
            print("Allowed system types: ligand, protein or ligand-protein")
            exit()

        self.CSD = md_params.get("CSD identifier")
        self.PDB = md_params.get("PDB identifier")
        self.solvate = md_params.get("solvate system")
        self.temp = md_params.get("temperature (K)")
        self.ensemble = md_params.get("ensemble")
        self.dt = md_params.get("timestep (fs)")
        self.pairnet = md_params.get("pair-net model")
        self.time = md_params.get("simulation time (ns)")

        return md_params


    def setup(self):
        """Set up an MD simulation.

        :param pdb: The input coordinates.

        :returns: Simulation object
        """
        from sys import stdout
        from openmm import LangevinMiddleIntegrator, CustomExternalForce, app, unit
        from openmm import NonbondedForce, HarmonicBondForce, HarmonicAngleForce, PeriodicTorsionForce
        from openmmforcefields.generators import GAFFTemplateGenerator
        from openff.toolkit.topology import Molecule

        print("Reading structure from PDB file...")
        pdb = app.PDBFile("input.pdb")

        # setup force field
        forcefield = app.ForceField("amber14-all.xml", "amber14/tip3p.xml")

        # non-standard residue needs to generate a force field template
        if self.system_type == "ligand":
            # TODO: ***CONVERSION FUNCTIONS***
            # TODO: This is where the coordinates are read in.
            # TODO: http://docs.openmm.org/latest/api-python/generated/openmm.app.pdbfile.PDBFile.html
            # TODO: PDBFile is a class from the OpenMM Python application layer.
            # TODO: This is where the molecule connectivity is read in.
            # TODO: https://docs.openforcefield.org/projects/toolkit/en/stable/api/generated/openff.toolkit.topology.Molecule.html
            ligand = Molecule.from_file("ligand.sdf")
            # TODO: Molecule is a class from OpenFF-toolkits.
            # TODO: ***CONVERSION FUNCTIONS***
            topology = ligand.to_topology().to_openmm()
            gaff = GAFFTemplateGenerator(molecules=ligand)
            forcefield.registerTemplateGenerator(gaff.generator)
            topology.setUnitCellDimensions([2.5]*3)
        elif self.system_type == "protein":
            # TODO: ***CONVERSION FUNCTIONS***
            # TODO: This is where the protein coordinates are read in.
            # TODO: http://docs.openmm.org/latest/api-python/generated/openmm.app.pdbfile.PDBFile.html
            # TODO: this is where the protein topology is defined.
            topology = pdb.topology
            # TODO: ***CONVERSION FUNCTIONS***
        elif self.system_type == "ligand-protein":
            pass

        modeller = app.Modeller(topology, pdb.positions)
        cutoff = 1.0*unit.nanometer

        if self.solvate:
            modeller.addSolvent(forcefield)
            n_solvent = len(modeller.getPositions()) - len(pdb.positions)
            print(f"{n_solvent} solvent molecules added...")

        # construct OpenMM system object using topology and other MD run parameters
        if self.system_type == "protein":
            system = forcefield.createSystem(modeller.topology,
                nonbondedCutoff=cutoff, nonbondedMethod=app.PME)

        elif self.system_type == "ligand":
            # charges assigned using am1bcc
            # bug in OpenFF Toolkit means this doesn't always work first time, so in a loop for now
            # TODO: sorting OpenEye license will probably avoid this problem
            while True:
                try:
                    # construct OpenMM system object using topology and other MD run parameters
                    system = forcefield.createSystem(modeller.topology, nonbondedCutoff=cutoff,
                                                     nonbondedMethod=app.PME)
                    print("Charge assignment succeeded...")
                    break
                except:
                    print("Charge assignment failed...")
                    pass

        # if using pair-net must set all intramolecular ligand interactions to zero
        if self.pairnet != "none":
            # exclude all non-bonded interactions
            nb = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]
            for i in range(ligand.n_atoms):
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
            for j in range(ligand.n_atoms):
                self.ml_force.addParticle(j, (0, 0, 0))
        else:
            self.ml_force = None

        # setup integrator
        self.temp = self.temp*unit.kelvin
        self.dt = self.dt*unit.femtoseconds
        temp_coupling = 1/unit.picosecond
        if self.ensemble == "NVT":
            integrator = LangevinMiddleIntegrator(self.temp, temp_coupling, self.dt)

        # setup simulation and output
        # Simulation object ties together topology, system and integrator
        # and maintains list of reporter objects that record or analyse data
        self.simulation = app.Simulation(modeller.topology, system, integrator)
        self.simulation.context.setPositions(modeller.positions)
        if self.ensemble == "NVT":
            self.simulation.context.setVelocitiesToTemperature(self.temp)
        self.simulation.reporters.append(app.PDBReporter("output.pdb", 10, enforcePeriodicBox=True))
        self.simulation.reporters.append(app.StateDataReporter(stdout, 10, step=True,
            potentialEnergy=True, temperature=True))
        if self.system_type == "protein":
            print("Minimising initial protein structure...")
            self.simulation.minimizeEnergy()

        return None


    def standard(self):
        """

        :param simulation:
        :return:
        """
        from openmm import unit
        import os
        os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
        import tensorflow as tf
        import numpy as np

        tf.get_logger().setLevel('ERROR')

        self.time = self.time*unit.nanoseconds
        n_steps = int(self.time/self.dt)

        if self.pairnet != "none":
            input_dir = f"./pair-net_models/{self.pairnet}"
            isExist = os.path.exists(input_dir)
            if not isExist:
                print("ERROR. Previously trained model could not be located.")
                exit()
            model, atoms = self.load_pairnet(input_dir)
            n_atoms = len(atoms)
            # this step is necessary because pairnet predicts forces for an arbitrary atom ordering
            # TODO: Ideally would not need a manually created mapping file.
            # TODO: For now, can remap atoms in this way because we only have models for a few molecules anyway.
            # TODO: This can be made more efficient if either:
            # TODO: 1) atom ordering of CCDC outputs can be changed
            # TODO: 2) mapping of topology/structure prior to starting simulation
            mapping = np.loadtxt(f"{input_dir}/atom_mapping.dat", dtype=int)
            csd2pairnet = mapping[:, 0]
            pairnet2csd = mapping[:, 1]
            for i in range(n_steps):
                coords = self.simulation.context.getState(getPositions=True). \
                    getPositions(asNumpy=True).value_in_unit(unit.angstrom)
                coords = coords[csd2pairnet] # map CSD to pairnet atom order
                if (i % 1000) == 0:
                    tf.keras.backend.clear_session()
                prediction = model.predict_on_batch([np.reshape(coords[:n_atoms]
                    / unit.angstrom, (1, -1, 3)), np.reshape(atoms, (1, -1))])
                ML_forces = prediction[0]*unit.kilocalories_per_mole/unit.angstrom
                ML_forces = np.reshape(ML_forces, (-1, 3))
                ML_forces = ML_forces[pairnet2csd] # map pairnet back to CSD
                for j in range(n_atoms):
                    self.ml_force.setParticleParameters(j, j, ML_forces[j])
                self.ml_force.updateParametersInContext(self.simulation.context)

                self.simulation.step(1)
        else:
            self.simulation.step(n_steps)

        return None


    def load_pairnet(self, input_dir):
        """

        :param simulation:
        :return:
        """
        from network import Network
        import numpy as np
        print("Loading a previously trained model...")
        atoms = np.loadtxt(f"./{input_dir}/trained_model/atoms.txt", dtype=np.float32).reshape(-1)
        ann_params = Network.read_params(f"{input_dir}/trained_model/ann_params.txt")
        prescale = np.loadtxt(f"./{input_dir}/trained_model/prescale.txt", dtype=np.float64).reshape(-1)
        network = Network()
        model = network.build(len(atoms), ann_params, prescale)
        model.summary()
        model.load_weights(f"./{input_dir}/trained_model/best_ever_model").expect_partial()
        return model, atoms

