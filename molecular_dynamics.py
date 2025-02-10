import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import tensorflow as tf
tf.get_logger().setLevel('ERROR')
import ccdc_convertor

class MolecularDynamics():

    import warnings
    warnings.filterwarnings("ignore", category=DeprecationWarning)

    def __init__(self):
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

        self.CSD = md_params.get("CSD identifier")
        self.PDB = md_params.get("PDB identifier")
        self.solvate = md_params.get("solvate system")
        self.temp = md_params.get("temperature (K)")
        self.ensemble = md_params.get("ensemble")
        self.dt = md_params.get("timestep (fs)")
        self.pairnet_path = md_params.get("pair-net model path")
        self.time = md_params.get("simulation time (ns)")
        self.type = md_params.get("simulation type")
        if md_params.get("analyse torsions"):
            self.analyse_torsions = True
        else:
            self.analyse_torsions = False

        self.partial_charges = md_params.get("partial charges")

        # run some checks
        system_type = md_params.get("system type")
        if system_type != "ligand" and system_type != "protein" and \
            system_type != "ligand-protein":
            print("ERROR - system type not allowed.")
            print("Allowed system types: ligand, protein or ligand-protein")
            exit()

        if system_type== "ligand":
            if self.partial_charges != "am1-bcc" and self.partial_charges != "from_file" \
                    and self.partial_charges != "predicted":
                print("ERROR - charge scheme not recognised.")
                print("Charge scheme must either be 'am1-bcc', 'from_file' or 'predicted'")
                print("Switching to 'am1-bcc'")
                self.partial_charges == "am1-bcc"

        if system_type == "ligand" or system_type == "ligand-protein":
            self.ligand = True
        else:
            self.ligand = False

        if system_type == "protein" or system_type == "ligand-protein":
            self.protein = True
        else:
            self.protein = False

        if self.partial_charges == "predicted" and self.pairnet_path == "none":
            print("ERROR - cannot predict partial charges without a pairnet model.")
            print("Provide path to pairnet model or change charge scheme.")
            exit()
        if system_type == "protein":
            if self.pairnet_path != "none":
                print("WARNING - PairNet cannot be used for protein simulations.")
                print("Switching to MM potential.")
                self.pairnet_path = "none"
            if self.partial_charges != "none":
                print("WARNING - must use standard partial charges for standard residues.")
                print("Resetting charge scheme...")
                self.partial_charges = "none"

        self.output_dir = "md_data/"
        isExist = os.path.exists(self.output_dir)
        if isExist:
            shutil.rmtree(self.output_dir)
        os.makedirs(self.output_dir)

        return md_params


    def setup(self):
        """Set up an MD simulation.
        ...
        """
        from openmm import app, unit
        from openmmforcefields.generators import SystemGenerator
        from openff.toolkit.topology import Molecule
        import warnings
        import numpy as np
        warnings.filterwarnings("ignore", category=DeprecationWarning)

        if self.CSD != "from_gro":

            # protein only simulation, positions from downloaded PDB file
            if self.protein and not self.ligand:
                self.topology, self.positions = get_pdb(f"{self.input_dir}protein.pdb")
                # TODO - can't yet use CCDC->OpenMM converter for proteins due to incorrect charge state of CCDC protein:
                # self.topology, self.positions = ccdc_convertor.openmm_topology_and_positions_from_ccdc_molecule(self.protein)

            # ligand only simulation, positions from first conformer
            if self.ligand and not self.protein:
                self.topology, self.positions = ccdc_convertor.\
                    openmm_topology_and_positions_from_ccdc_molecule(self.conformers[0])
                if not self.solvate:
                    self.topology.setUnitCellDimensions([3.0]*3)

            # ligand-protein simulation, ligand positions from first docked pose
            if self.ligand and self.protein:
                # TODO - can't yet use CCDC->OpenMM converter for proteins due to incorrect charge state of CCDC protein:
                # self.topology, self.positions = ccdc_convertor.openmm_topology_and_positions_from_ccdc_molecule(self.protein)
                protein_topology, protein_positions = get_pdb(f"{self.input_dir}protein.pdb")
                ligand_topology, ligand_positions = ccdc_convertor. \
                    openmm_topology_and_positions_from_ccdc_molecule(self.conformers[0])
                self.topology, self.positions = merge_topology(ligand_topology,
                    ligand_positions, protein_topology, protein_positions)

             # define force field
            std_ff = ["amber/ff14SB.xml", "amber/tip3p_standard.xml"]
            small_mol_ff = "gaff-2.11"

            # create topology for non-standard residue from SMILES
            if self.ligand:
                print("Creating topology...")
                ligand = Molecule.from_smiles(self.smiles, allow_undefined_stereo=True)
                self.ligand_n_atom = ligand.n_atoms
                if self.partial_charges == "from_file":
                    print("Reading partial charges from file...")
                    ligand.partial_charges = np.loadtxt("charges.txt") * unit.elementary_charges
                molecules = ligand
            else:
                molecules = None

            # add water to system using PDBfixer
            if self.solvate:
                input_file = f"{self.input_dir}ligand.pdb"
                solvated_pdb = solvate_system_pdbfixer(input_file)
                n_solvent = len(solvated_pdb.positions) - len(self.positions)
                self.topology = solvated_pdb.topology
                self.positions = solvated_pdb.positions
                print(f"Adding {int(n_solvent / 3)} solvent molecules...")

            # construct M system using topology and forcefield
            print("Creating system...")
            system_generator = SystemGenerator(forcefields=std_ff, small_molecule_forcefield=small_mol_ff,
                forcefield_kwargs={"constraints": None, "rigidWater": True}, cache="db.json")
            self.system = system_generator.create_system(self.topology, molecules=molecules)

        else:
            self.input_dir = "md_input_gro/"
            self.n_conf = 1
            gro = app.GromacsGroFile(f"{self.input_dir}input.gro")
            top = app.GromacsTopFile(f"{self.input_dir}input.top",
                                 periodicBoxVectors=gro.getPeriodicBoxVectors())
            residues = list(top.topology.residues())
            self.ligand_n_atom = len(list(residues[0].atoms()))
            self.system = top.createSystem(constraints=None, removeCMMotion=True,
                                      rigidWater=True)
            self.positions = gro.positions
            self.topology = top.topology

        if self.ligand:

            self.ml_force = None
            self.fixed_charges = assign_fixed_charges(self.system, self.partial_charges, self.ligand_n_atom)

            # add metadynamics force to system
            if self.type == "enhanced":
                from openmmplumed import PlumedForce
                print("Using metadynamics bias...")
                print()
                if self.CSD != "from_gro":
                    plumed_script = write_plumed_script(self.unique_torsions)
                else:
                    plumed_file = open(f"{self.input_dir}/plumed.dat", "r")
                    plumed_script = plumed_file.read()
                print("PLUMED input...")
                print(plumed_script)
                self.system.addForce(PlumedForce(plumed_script))

        # setup simulation and output tying together topology, system and integrator
        print("Creating simulation...")
        self.integrator = get_integrator(self.temp, self.dt, self.ensemble, self.pairnet_path)
        self.simulation = app.Simulation(self.topology, self.system, self.integrator)
        self.simulation.context.setPositions(self.positions)

        print("Minimising initial structure...")
        self.simulation.minimizeEnergy()

        # if using pair-net create new system with ML force and minised coords
        if self.pairnet_path != "none":
            print("Using ML force field (PairNet) for ligand...")
            minimised_coords = self.simulation.context.getState(
                getPositions=True).getPositions(asNumpy=True)
            self.system = turn_off_MM(self.system, self.ligand_n_atom)
            self.system, self.ml_force = create_MLP(self.system, self.ligand_n_atom)
            self.integrator = get_integrator(self.temp, self.dt, self.ensemble, self.pairnet_path)
            self.simulation = app.Simulation(self.topology, self.system, self.integrator)
            self.simulation.context.setPositions(minimised_coords)

        # fix initial velocity seed for now
        if self.ensemble == "NVT":
            self.simulation.context.setVelocitiesToTemperature(self.temp)

        self.simulation.reporters.append(app.PDBReporter("output.pdb", 1000,
                                    enforcePeriodicBox=True))

        return None


    def simulate(self):
        """

        :param:
        :return:
        """
        from openmm import unit, app
        import numpy as np
        import time as timer
        start_time = timer.time()

        self.time = self.time*unit.nanoseconds
        self.dt = self.dt*unit.femtoseconds
        n_steps = int(self.time/self.dt)+2

        if self.ligand:

            if self.pairnet_path != "none":

                input_dir = f"{self.pairnet_path}trained_model/"
                print(f"Using pairnet model: {input_dir}")
                isExist = os.path.exists(input_dir)
                if not isExist:
                    print("ERROR. Previously trained model could not be located.")
                    exit()
                model, atoms = load_pairnet_v2(input_dir, self.ligand_n_atom)

            charges = self.fixed_charges

            # open files for storing PairNetOps compatible datasets
            f1 = open(f"{self.output_dir}coords.txt", "w")
            f2 = open(f"{self.output_dir}forces.txt", "w")
            f3 = open(f"{self.output_dir}energies.txt", "w")
            f4 = open(f"{self.output_dir}charges.txt", "w")

            # open files for storing PairNetOps compatible datasets
            data_files = [f1, f2, f3, f4]

        print("Performing MD simulation...")
        print("Time (ps) | PE (kcal/mol)")
        stable = True
        for i in range(n_steps):

            if self.pairnet_path != "none":

                # this stops us running out of memory
                if (i % 1000) == 0:
                    tf.keras.backend.clear_session()

                coords = self.simulation.context.getState(getPositions=True). \
                    getPositions(asNumpy=True).value_in_unit(unit.angstrom)[:self.ligand_n_atom]

                prediction = get_pairnet_prediction(model, atoms, coords)
                stable = check_stability(i, prediction[0])
                if not stable:
                    break
                ML_forces = prediction[0] * unit.kilocalories_per_mole / unit.angstrom

                ML_forces = np.reshape(ML_forces, (-1, 3))

                if self.partial_charges == "predicted":
                    charges = prediction[2].T
                    if i == 0:
                        net_charge = sum(charges[:,0])
                        net_charge = round(net_charge)
                    charges = set_charges(charges, self.system, self.simulation,
                        self.ligand_n_atom, net_charge)

                for j in range(self.ligand_n_atom):
                    self.ml_force.setParticleParameters(j, j, ML_forces[j])
                self.ml_force.updateParametersInContext(self.simulation.context)

            # every 1000 steps save data for PairNetOps compatible dataset
            if (i % 100) == 0:

                state = self.simulation.context.getState(getEnergy=True)
                energy = state.getPotentialEnergy() / unit.kilocalories_per_mole

                if self.ligand:
                    coords = self.simulation.context.getState(getPositions=True). \
                        getPositions(asNumpy=True).value_in_unit(unit.angstrom)
                    forces = self.simulation.context.getState(getForces=True). \
                        getForces(asNumpy=True).in_units_of(unit.kilocalories_per_mole / unit.angstrom)

                    if self.pairnet_path != "none":
                        energy = prediction[1][0][0]

                    write_dataset(self.ligand_n_atom, coords, forces, energy, charges, data_files)

                time = self.simulation.context.getState().getTime(). \
                    value_in_unit(unit.picoseconds)
                print(f"{time:9.1f} | {energy:13.2f}")

            # advance trajectory one timestep
            self.simulation.step(1)

        end_time = timer.time()
        run_time = end_time - start_time

        if not stable:
            positions = self.simulation.context.getState(getPositions=True).getPositions()
            app.PDBFile.writeFile(self.simulation.topology, positions, open("final.pdb", "w"))
            print(f"MD simulation ended prematurely in {timer.strftime('%H:%M:%S', timer.gmtime(run_time))}.")
        else:
            print(f"MD simulation has completed in {timer.strftime('%H:%M:%S', timer.gmtime(run_time))}.")

        if self.ligand:
            f1.close()
            f2.close()
            f3.close()
            f4.close()

        return None


def load_pairnet_v2(input_dir, ligand_n_atom):
    """

    :param input_dir: path to pairnet model
    :returns: model - pairnet model
              atoms - list of atoms
    """
    from network_v2 import Network
    print("Loading a previously trained model...")
    network = Network()
    model = network.load(input_dir)
    model.summary()
    model.load_weights(f"{input_dir}best_ever_model").expect_partial()
    atoms = network.atoms
    if len(atoms) != ligand_n_atom:
        print("ERROR - number of atoms in trained network is incompatible with number of atoms in topology")
        exit()
    return model, atoms


def get_pairnet_prediction(model, atoms, coords):
    import numpy as np
    from openmm import unit

    prediction = model.predict_on_batch([np.reshape(coords[:len(atoms)]
        / unit.angstrom, (1, -1, 3)), np.reshape(atoms, (1, -1))])

    return prediction


def get_pdb(filename):
    from openmm import app
    pdb = app.PDBFile(filename)
    return pdb.topology, pdb.positions


def solvate_system_pdbfixer(filename):
    from pdbfixer import PDBFixer
    from openmm import Vec3
    fixer = PDBFixer(filename=filename)
    boxSize = 3.0 * Vec3(1, 1, 1)
    fixer.addSolvent(boxSize=boxSize)
    return fixer

def solvate_system_modeller(topology, positions, std_ff):
    from openmm import app, Vec3
    modeller = app.Modeller(topology, positions)
    boxSize = 3.0 * Vec3(1, 1, 1)
    forcefield = app.ForceField(std_ff[0], std_ff[1])
    modeller.addSolvent(forcefield, boxSize=boxSize)
    return modeller


def create_MLP(system, n_atom):

    from openmm import CustomExternalForce

    # create custom force for PairNet predictions
    ml_force = CustomExternalForce("-fx*x-fy*y-fz*z")
    system.addForce(ml_force)
    ml_force.addPerParticleParameter("fx")
    ml_force.addPerParticleParameter("fy")
    ml_force.addPerParticleParameter("fz")

    for j in range(n_atom):
        ml_force.addParticle(j, (0, 0, 0))

    return system, ml_force


def turn_off_MM(system, n_atom):

    from openmm import NonbondedForce, HarmonicBondForce, HarmonicAngleForce
    from openmm import PeriodicTorsionForce

    # exclude all non-bonded interactions
    nb = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]
    for i in range(n_atom):
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

    return system


def assign_fixed_charges(system, charge_scheme, ligand_n_atom):
    from openmm import NonbondedForce, unit
    import numpy as np
    fixed_charges = np.zeros(ligand_n_atom)
    if charge_scheme == "from_file":
        file_charge = np.loadtxt("charges.txt")
    nb = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]
    for i in range(ligand_n_atom):
        [charge, sigma, epsilon] = nb.getParticleParameters(i)
        if charge_scheme == "from_file":
            fixed_charges[i] = file_charge[i]
        else:
            fixed_charges[i] = charge / unit.elementary_charge
        nb.setParticleParameters(i, fixed_charges[i] * unit.elementary_charge, sigma, epsilon)
    return fixed_charges


def write_dataset(n_atom, coords, forces, energy, charges, data_files):
    import numpy as np
    np.savetxt(data_files[0], coords[:n_atom])
    np.savetxt(data_files[1], forces[:n_atom])
    data_files[2].write(f"{energy}\n")
    np.savetxt(data_files[3], charges[:n_atom])
    return None


def write_plumed_script(torsions):

    # each torsion will be allocated to a CV
    torsion_CV = [None]*len(torsions)

    # start CV counter
    CV_count = 1

    # first torsion will always be allocated to new CV
    torsion_CV[0] = CV_count

    # loop through all torsions to identify all CVs
    for i_tors in range (1, len(torsions)):

        # get central atom pair for this torsion
        pair_1 = torsions[i_tors][1:-1]

        # assume we will create a new CV
        new_CV = True

        # loop over other torsions
        for j_tors in range(i_tors):

            # get central atom pair for this torsion
            pair_2 = torsions[j_tors][1:-1]

            # check to see if there are any shared atom indices in the two torsions
            adjacent_torsions = not set(pair_1).isdisjoint(pair_2)

            # if the two pairs do not share any atoms it may be a new CV
            if adjacent_torsions:
                new_CV = False
                torsion_CV[i_tors] = torsion_CV[j_tors]
                break

        # if no adjacent torsions were found assign this torsion to a new CV
        if new_CV:
            CV_count += 1
            torsion_CV[i_tors] = CV_count

    bias_pace = 500 # frequency of depositing Gaussians
    bias_height = 1.0 # Gaussian height in kJ/mol
    bias_sigma = 0.35 # Gaussian width

    template_text = [": TORSION ATOMS=", ": METAD ARG=", f" PACE={bias_pace}"
                    f" HEIGHT={bias_height} SIGMA=", " FILE=./HILLS_",
                     "PRINT STRIDE=10 ARG=", ".bias FILE=./COLVAR_"]

    plumed_text = ""
    # write all torsion angles for biasing
    for i_tors in range(len(torsions)):
        indices = ",".join(str(x+1) for x in torsions[i_tors])
        text = f"phi{i_tors+1}" + template_text[0] + indices + "\n"
        plumed_text = plumed_text + text
    plumed_text = plumed_text + "\n"

    # for each CV bias determine which torsions are in it
    for i_CV in range(CV_count):
        text = f"metad{i_CV + 1}" + template_text[1]
        torsion_list = [i for i, n in enumerate(torsion_CV) if n == (i_CV+1)]
        for i_tors in range(len(torsion_list)):
            text = text + f"phi{torsion_list[i_tors]+1}"
            if i_tors != len(torsion_list)-1:
                text = text + ","
        text = text + template_text[2]
        for i_tors in range(len(torsion_list)):
            text = text + f"{bias_sigma}"
            if i_tors != len(torsion_list)-1:
                text = text + ","
        text = text + template_text[3] + f"{i_CV+1}" + "\n"
        plumed_text = plumed_text + text
    plumed_text = plumed_text + "\n"

    for i_CV in range(CV_count):
        text = template_text[4]
        torsion_list = [i for i, n in enumerate(torsion_CV) if n == (i_CV+1)]
        for i_tors in range(len(torsion_list)):
            text = text + f"phi{torsion_list[i_tors]+1}" + ","
        text = text + f"metad{i_CV+1}" + template_text[5] + f"{i_CV+1}" + "\n"
        plumed_text = plumed_text + text

    return plumed_text


def get_integrator(temp, dt, ensemble, path):
    from openmm import unit
    from openmm import LangevinMiddleIntegrator, NoseHooverIntegrator
    temp = temp * unit.kelvin
    dt = dt * unit.femtoseconds
    temp_coupling = 1 / unit.picosecond
    if ensemble == "NVT":
        if path != "none":
            print("Using Langevin thermostat...")
            integrator = LangevinMiddleIntegrator(temp, temp_coupling, dt)
        else:
            print("Using Nose-Hoover thermostat...")
            integrator = NoseHooverIntegrator(temp, temp_coupling, dt)
    return integrator


def set_charges(predicted_charges, system, simulation, n_atom, net_charge):
    from openmm import NonbondedForce
    from openmm import unit

    # adjust to give correct net charge
    corr = (sum(predicted_charges[:,0]) - net_charge) / n_atom
    corrected_charges = (predicted_charges - corr) * unit.elementary_charge

    nbforce = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]
    for j in range(n_atom):
        [old_charge, sigma, epsilon] = nbforce.getParticleParameters(j)
        nbforce.setParticleParameters(j, corrected_charges[j][0], sigma, epsilon)
    nbforce.updateParametersInContext(simulation.context)

    return corrected_charges


def merge_topology(topology_1, positions_1, topology_2, positions_2):
    from openmm import app
    cell_dims = topology_2.getUnitCellDimensions()
    modeller = app.Modeller(topology_1, positions_1)
    modeller.add(topology_2, positions_2)
    modeller.topology.setUnitCellDimensions(cell_dims)
    return modeller.topology, modeller.positions


def check_stability(i, prediction):
    import numpy as np
    # exit simulation if any predicted forces are greater than threshold
    stable = True
    if np.any(prediction >= 1000.0):
        stable = False
        print(f"Error - predicted force exceed stability threshold in step {i}")
        print("Final forces:")
        print(prediction)
        print()
    return stable

