#!/usr/bin/env python3
__author__ = ['Christopher D Williams']
__credits__ = ['Kepa Burusco-Goni, Simon Cottrell, Austin Lloyd, '
               'Bojana Popovic, Richard Bryce ']
__license__ = '...'
__maintainer__ = 'Christopher D Williams'
__email__ = 'christopher.williams@manchester.ac.uk'
__status__ = 'Development'
__funding__ = "Funded by UKRI IAA"


def main():
    """

    :return:
    """
    from molecular_dynamics import MolecularDynamics
    import os, shutil
    import warnings, get_structure
    warnings.filterwarnings("ignore", category=DeprecationWarning)

    print("Reading input parameters...")
    simulation = MolecularDynamics()
    # TODO: move read_inputs to init method?
    simulation.read_inputs()

    if simulation.CSD != "from_gro":
        simulation.input_dir = "md_input/"
        isExist = os.path.exists(simulation.input_dir)
        if isExist:
            shutil.rmtree(simulation.input_dir)
        os.makedirs(simulation.input_dir)

        if simulation.ligand and not simulation.protein:
            print(f"Retrieving CSD entry for {simulation.CSD}...")
            get_structure.ligand(simulation.CSD, simulation)
            print(f"SMILES notation: {simulation.smiles} ")
        elif simulation.protein and not simulation.ligand:
            print("Retrieving PDB...")
            get_structure.get_protein(simulation)
            print("Fixing PDB...")
            get_structure.fix_protein_ccdcfixer(simulation)
        else:
            get_structure.docking(simulation)

    print("Setting up MD simulation...")
    simulation.setup()

    simulation.simulate()

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()

