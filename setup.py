from sys import stdout
from openmm import LangevinMiddleIntegrator, NonbondedForce, CustomExternalForce, app
from openmm import unit
from openmmforcefields.generators import GAFFTemplateGenerator
from openff.toolkit.topology import Molecule

class MolecularDynamics:

    def openMM(self, pdb, md_params):
        """Set up an MD simulation.

        :param pdb: The input coordinates.

        :returns: Simulation object
        """

        # setup force field
        forcefield = app.ForceField("amber14-all.xml", "amber14/tip3p.xml")

        # non-standard residue needs to generate a force field template
        if md_params.get("system type") == "ligand":
            ligand = Molecule.from_file("input.sdf")
            topology = ligand.to_topology().to_openmm()
            topology.setUnitCellDimensions([3.0]*3)
            gaff = GAFFTemplateGenerator(molecules=ligand)
            forcefield.registerTemplateGenerator(gaff.generator)
        elif md_params.get("system type") == "protein":
            topology = pdb.topology

        modeller = app.Modeller(topology, pdb.positions)
        cutoff = 1.0 * unit.nanometer

        # charges assigned using am1bcc - bug in OpenFF Toolkit means this doesn't always work.
        while True:
            try:
                if md_params.get("solvate system"):
                    modeller.addSolvent(forcefield)
                system = forcefield.createSystem(modeller.topology, nonbondedCutoff=cutoff)
                print("Charge assignment succeeded...")
                break
            except:
                print("Charge assignment failed...")
                pass

        # if using pair-net set all intramolecular ligand interactions to zero
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

        # setup integrator
        temp = md_params.get("temperature (K)")*unit.kelvin
        dt = md_params.get("timestep (fs)")*unit.femtoseconds
        temp_coupling = 1/unit.picosecond
        if md_params.get("ensemble") == "NVT":
            integrator = LangevinMiddleIntegrator(temp, temp_coupling, dt)

        # setup simulation and output
        simulation = app.Simulation(modeller.topology, system, integrator)
        simulation.context.setPositions(modeller.positions)
        simulation.context.setVelocitiesToTemperature(temp)
        simulation.reporters.append(app.PDBReporter("output.pdb", 10))
        simulation.reporters.append(app.StateDataReporter(stdout, 10, step=True,
            potentialEnergy=True, temperature=True))

        return system, simulation

