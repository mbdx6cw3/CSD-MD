#!/usr/bin/env python3
import read_inputs, setup, simulate, get_structure

# TODO: apply to all modules
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

def main():
    """

    :return:
    """

    print("Reading input parameters...")
    md_params = read_inputs.csdMD()

    print("Simulation:", md_params.get("name"))
    print()

    if md_params.get("identifier") != "n/a":
        print("Retrieving entry from CSD...")
        entry = get_structure.CSDDatabase()
        if md_params.get("system type") == "ligand":
            entry.ligand(md_params.get("identifier"))
        elif md_params.get("system type") == "protein":
            entry.protein(md_params.get("identifier"))

    # TODO: protein-ligand docking.
    # TODO:  atoms

    print("Reading structure from PDB file...")
    input_structure = get_structure.pdb()

    print("Setting up MD simulation...")
    simulation, force = setup.MolecularDynamics().openMM(input_structure, md_params)

    print("Minimising initial structure...")
    simulation.minimizeEnergy()

    print("Performing MD simulation...")
    simulate.MolecularDynamics().standard(md_params, simulation, force)

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()
