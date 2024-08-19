from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
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

        if md_params.get("system type") == "small_molecule":
            molecule = Molecule.from_file("input.sdf")
            gaff = GAFFTemplateGenerator(molecules=molecule, forcefield="gaff-2.11")
            forces = ForceField(solute_FF, solvent_FF)
            forces.registerTemplateGenerator(gaff.generator)
        elif md_params.get("system type") == "protein":
            forces = ForceField(solute_FF, solvent_FF)

        system = forces.createSystem(pdb.topology, nonbondedMethod=PME, 
            nonbondedCutoff=1*nanometer, constraints=HBonds)

        # setup integrator
        temp = md_params.get("temperature (K)")*kelvin
        dt = md_params.get("timestep (fs)")*femtoseconds
        temp_coupling = 1/picosecond
        print(temp, dt, temp_coupling)
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

