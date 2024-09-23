#!/usr/bin/env python3
__author__ = ['Christopher D Williams']
__credits__ = ['Kepa Burusco-Goni, Simon Cottrell, Austin Lloyd, '
               'Bojana Popovic, Richard Bryce ']
__license__ = '...'
__maintainer__ = 'Christopher D Williams'
__email__ = 'christopher.williams@manchester.ac.uk'
__status__ = 'Development'


def main():
    """

    :return:
    """
    import read_inputs, setup, simulate, get_structure
    import warnings
    warnings.filterwarnings("ignore", category=DeprecationWarning)

    print("Reading input parameters...")
    md_params = read_inputs.csdMD()

    print("Simulation:", md_params.get("name"))
    print()

    if md_params.get("system type") == "ligand":
        ligand = True
        protein = False
    elif md_params.get("system type") == "protein":
        ligand = True
        protein = False
    elif md_params.get("system type") == "ligand-protein":
        ligand = True
        protein = True
    else:
        print("ERROR - system type not allowed.")
        print("Allowed system types: ligand, protein or ligand-protein")
        exit()

    if ligand:
        print("Retrieving CSD entry...")
        entry = get_structure.CSDDatabase()
        entry.ligand(md_params.get("CSD identifier"))

    if protein:
        print("Retrieving PDB entry...")
        entry = get_structure.PDBDatabase()
        entry.protein(md_params.get("PDB identifier"))

    # TODO: protein-ligand docking option here.

    print("Setting up MD simulation...")
    simulation, force = setup.MolecularDynamics().openMM(md_params)

    print("Performing MD simulation...")
    simulate.MolecularDynamics().standard(md_params, simulation, force)

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()
