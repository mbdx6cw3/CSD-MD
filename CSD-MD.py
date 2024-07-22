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

    '''
    identifier = md_params.get("identifier")
    if identifier != "n/a":
        print("Retrieving entry from CSD...")
        entry = get_structure.CSDDatabase()
        if md_params.get("system type") == "small_molecule":
            entry.smallmolecule(identifier)
    '''

    print("Reading structure from PDB file...")
    input_structure = get_structure.pdb()

    # TODO: prepare structure for simulation (e.g. solvation, box size, add Hs)

    print("Initialising MD simulation...")
    simulation_setup = setup.MolecularDynamics().openMM(input_structure, md_params)

    print("Minimising initial structure...")
    simulation_setup.minimizeEnergy()

    print("Performing a", md_params.get("simulation time (ns)"), "ns MD simulation...")
    simulate.MolecularDynamics().standard(simulation_setup, md_params)

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()
