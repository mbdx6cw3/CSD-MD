from sys import stdout
from openmm.app import *
from openmm import *
from openmm.unit import *
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
        forces = ForceField(solute_FF, solvent_FF)

        # TODO: if non-standard molecule generate a force field template
        if md_params.get("system type") == "small_molecule":
            molecule = Molecule.from_file("input.sdf")
            gaff = GAFFTemplateGenerator(molecules=molecule, forcefield="gaff-2.11")
            forces.registerTemplateGenerator(gaff.generator)

        system = forces.createSystem(pdb.topology, nonbondedMethod=PME, 
            nonbondedCutoff=1*nanometer, constraints=HBonds)
        system.setDefaultPeriodicBoxVectors(Vec3(3.0, 0, 0), Vec3(0, 3.0, 0), Vec3(0, 0, 3.0))

        # setup integrator
        temp = md_params.get("temperature (K)")*kelvin
        dt = md_params.get("timestep (fs)")*femtoseconds
        temp_coupling = 1/picosecond
        if md_params.get("ensemble") == "NVT":
            integrator = LangevinMiddleIntegrator(temp, temp_coupling, dt)

        # setup simulation and output
        simulation = Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        simulation.context.setVelocitiesToTemperature(temp)
        simulation.reporters.append(PDBReporter("output.pdb", 10))
        simulation.reporters.append(StateDataReporter(stdout, 10, step=True,
            potentialEnergy=True, temperature=True))

        return simulation

