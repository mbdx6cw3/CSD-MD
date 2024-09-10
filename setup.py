from sys import stdout
from openmm import LangevinMiddleIntegrator, NonbondedForce, CustomExternalForce, app
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
        # TODO: *** CONVERSION FUNCTIONS ***

        # non-standard residue needs to generate a force field template
        if md_params.get("system type") == "ligand":
            # TODO: *** CONVERSION FUNCTIONS ***
            # TODO: This is where the molecule connectivity is read in.
            # TODO: https://docs.openforcefield.org/projects/toolkit/en/stable/api/generated/openff.toolkit.topology.Molecule.html
            ligand = Molecule.from_file("input.sdf")
            # TODO: *** CONVERSION FUNCTIONS ***
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

        if md_params.get("system type") == "protein":
            system = forcefield.createSystem(modeller.topology,
                                             nonbondedCutoff=cutoff,
                                             nonbondedMethod=app.PME)
        elif md_params.get("system type") == "ligand":
            # charges assigned using am1bcc
            # bug in OpenFF Toolkit means this doesn't always work so in a loop for now
            # TODO: sorting OpenEye license will probably avoid this problem
            while True:
                try:
                    system = forcefield.createSystem(modeller.topology, nonbondedCutoff=cutoff,
                                                     nonbondedMethod=app.PME)
                    print("Charge assignment succeeded...")
                    break
                except:
                    print("Charge assignment failed...")
                    pass

        # if using pair-net load a model and set all intramolecular ligand interactions to zero
        if md_params.get("pair-net"):
            nb = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]
            for i in range(ligand.n_atoms):
                for j in range(i):
                    nb.addException(i, j, 0, 1, 0, replace=True)

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
        simulation = app.Simulation(modeller.topology, system, integrator)
        simulation.context.setPositions(modeller.positions)
        if md_params.get("ensemble") == "NVT":
            simulation.context.setVelocitiesToTemperature(temp)
        simulation.reporters.append(app.PDBReporter("output.pdb", 10, enforcePeriodicBox=True))
        simulation.reporters.append(app.StateDataReporter(stdout, 10, step=True,
            potentialEnergy=True, temperature=True))

        return simulation, ml_force

