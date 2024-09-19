#!/usr/bin/env python3
__author__ = ['Christopher D Williams']
__credits__ = [...]
__license__ = '...'
__maintainer__ = 'Christopher D Williams'
__email__ = 'christopher.williams@manchester.ac.uk'
__status__ = 'Development'


def main():
    """

    :return:
    """
    import read_inputs, setup, simulate, get_structure

    # TODO: apply to all modules
    import warnings
    warnings.filterwarnings("ignore", category=DeprecationWarning)

    print("Reading input parameters...")
    md_params = read_inputs.csdMD()

    print("Simulation:", md_params.get("name"))
    print()


    if md_params.get("system type") == "ligand":
        print("Retrieving CSD entry...")
        entry = get_structure.CSDDatabase()
        entry.ligand(md_params.get("CSD identifier"))
    elif md_params.get("system type") == "protein":
        print("Retrieving PDB entry...")
        entry = get_structure.PDBDatabase()
        entry.protein(md_params.get("PDB identifier"))
    else:
        pass

    # TODO: protein-ligand docking option here.

    print("Setting up MD simulation...")
    simulation, force = setup.MolecularDynamics().openMM(md_params)

    print("Performing MD simulation...")
    simulate.MolecularDynamics().standard(md_params, simulation, force)

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()
