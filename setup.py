from sys import stdout
from openmm import LangevinMiddleIntegrator, NonbondedForce, app
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
        if md_params.get("force field") == "amber":
            solute_FF = "amber14-all.xml"
            solvent_FF = "amber14/tip3p.xml"
        forcefield = app.ForceField(solute_FF, solvent_FF)

        # non-standard residue needs to generate a force field template
        if md_params.get("system type") == "ligand":
            molecule = Molecule.from_file("input.sdf")
            topology = molecule.to_topology().to_openmm()
            topology.setUnitCellDimensions([3.0]*3)
            gaff = GAFFTemplateGenerator(molecules=molecule)
            forcefield.registerTemplateGenerator(gaff.generator)

        modeller = app.Modeller(topology, pdb.positions)
        cutoff = 1.0 * unit.nanometer

        # charges assigned using am1bcc - bug in OpenFF Toolkit means this doesn't always work.
        while True:
            try:
                if md_params.get("solvated" == "yes"):
                    modeller.addSolvent(forcefield)
                system = forcefield.createSystem(modeller.topology, nonbondedCutoff=cutoff)
                print("Charge assignment succeeded...")
                break
            except:
                print("Charge assignment failed...")
                pass

        if md_params.get("system type") == "ligand":
            nonbonded = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]
            for i in range(system.getNumParticles()):
                charge, sigma, epsilon = nonbonded.getParticleParameters(i)

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

        return simulation

