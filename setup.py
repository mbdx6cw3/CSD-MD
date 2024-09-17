from sys import stdout
from openmm import LangevinMiddleIntegrator, CustomExternalForce, app
from openmm import NonbondedForce, HarmonicBondForce, HarmonicAngleForce, PeriodicTorsionForce
from openmm import unit
from openmmforcefields.generators import GAFFTemplateGenerator
from openff.toolkit.topology import Molecule

class MolecularDynamics:

    def openMM(self, md_params):
        """Set up an MD simulation.

        :param pdb: The input coordinates.

        :returns: Simulation object
        """

        # setup force field
        forcefield = app.ForceField("amber14-all.xml", "amber14/tip3p.xml")
        # TODO: ***CONVERSION FUNCTIONS***
        # TODO: This is where the coordinates are read in.
        # TODO: http://docs.openmm.org/latest/api-python/generated/openmm.app.pdbfile.PDBFile.html
        # TODO: https://simtk.org/api_docs/openmm/api6_0/python/classsimtk_1_1openmm_1_1app_1_1pdbfile_1_1PDBFile.html#a5e8a38af13069a0cc3fff9aae26892e4
        print("Reading structure from PDB file...")
        pdb = app.PDBFile('input.pdb')
        # TODO: PDBFile is a class from the OpenMM Python application layer.
        # TODO: ***CONVERSION FUNCTIONS***

        # non-standard residue needs to generate a force field template
        if md_params.get("system type") == "ligand":
            # TODO: ***CONVERSION FUNCTIONS***
            # TODO: This is where the molecule connectivity is read in.
            # TODO: https://docs.openforcefield.org/projects/toolkit/en/stable/api/generated/openff.toolkit.topology.Molecule.html
            ligand = Molecule.from_file("input.sdf")
            # TODO: Molecule is a class from OpenFF-toolkits.
            # TODO: ***CONVERSION FUNCTIONS***
            topology = ligand.to_topology().to_openmm()
            gaff = GAFFTemplateGenerator(molecules=ligand)
            forcefield.registerTemplateGenerator(gaff.generator)
            topology.setUnitCellDimensions([2.5]*3)
        elif md_params.get("system type") == "protein":
            topology = pdb.topology

        modeller = app.Modeller(topology, pdb.positions)
        cutoff = 1.0*unit.nanometer

        if md_params.get("solvate system"):
            modeller.addSolvent(forcefield)
            n_solvent = len(modeller.getPositions()) - ligand.n_atoms
            print(f"{n_solvent} solvent molecules added...")

        # construct OpenMM system object using topology and other MD run parameters
        if md_params.get("system type") == "protein":
            system = forcefield.createSystem(modeller.topology,
                nonbondedCutoff=cutoff, nonbondedMethod=app.PME)

        elif md_params.get("system type") == "ligand":
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
        if md_params.get("pair-net model"):

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
            ml_force = CustomExternalForce("-fx*x-fy*y-fz*z")
            system.addForce(ml_force)
            ml_force.addPerParticleParameter("fx")
            ml_force.addPerParticleParameter("fy")
            ml_force.addPerParticleParameter("fz")
            for j in range(ligand.n_atoms):
                ml_force.addParticle(j, (0, 0, 0))
        else:
            ml_force = None

        # setup integrator
        temp = md_params.get("temperature (K)")*unit.kelvin
        dt = md_params.get("timestep (fs)")*unit.femtoseconds
        temp_coupling = 1/unit.picosecond
        if md_params.get("ensemble") == "NVT":
            integrator = LangevinMiddleIntegrator(temp, temp_coupling, dt)

        # setup simulation and output
        # Simulation object ties together topology, system and integrator
        # and maintains list of reporter objects that record or analyse data
        simulation = app.Simulation(modeller.topology, system, integrator)
        simulation.context.setPositions(modeller.positions)
        if md_params.get("ensemble") == "NVT":
            simulation.context.setVelocitiesToTemperature(temp)
        simulation.reporters.append(app.PDBReporter("output.pdb", 1, enforcePeriodicBox=True))
        simulation.reporters.append(app.StateDataReporter(stdout, 1, step=True,
            potentialEnergy=True, temperature=True))

        if not md_params.get("pair-net model"):
            print("Minimising initial structure...")
            simulation.minimizeEnergy()

        return simulation, ml_force

