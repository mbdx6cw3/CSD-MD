#!/usr/bin/env python3
import read_inputs, get_structure, setup, simulate

def main():
    """

    :return:
    """

    print("Reading input parameters...")
    md_params = read_inputs.csdMD()

    print(md_params.get("name"))
    print()

    # retrieve initial structure from CSD or directly from PDB
    if md_params.get("identifier") != "n/a":
        print("Getting structure from CSD...")
        entry = get_structure.CSDDatabase()
        if md_params.get("small_molecule"):
            input_structure = entry.smallmolecule()
    else:
        print("Reading structure from PDB file...")
        input_structure = get_structure.pdb()

    # TODO: prepare structure for simulation (e.g. solvation, box size, add Hs)

    print("Initialising MD simulation...")
    simulation_setup = setup.MolecularDynamics().openMM(input_structure, md_params)

    print("Minimising initial structure...")
    simulation_setup.minimizeEnergy()

    print("Performing MD simulation...")
    simulate.MolecularDynamics().standard(simulation_setup)

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()
