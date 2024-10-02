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
    from get_structure import GetStructure
    from molecular_dynamics import MolecularDynamics
    import warnings
    warnings.filterwarnings("ignore", category=DeprecationWarning)

    simulation = MolecularDynamics()
    structure = GetStructure()

    print("Reading input parameters...")
    simulation.read_inputs()

    if simulation.system_type != "protein":
        print("Retrieving CSD entry...")
        n_conf = structure.ligand(simulation.CSD, simulation)
        print(f"Using {n_conf} conformer(s)...")

    if simulation.system_type != "ligand":
        print("Retrieving PDB entry...")
        structure.protein(simulation.PDB)

    if simulation.system_type == "ligand-protein":
        # TODO: protein-ligand docking option here.
        pass

    print("Setting up MD simulation...")
    simulation.setup()

    print("Performing MD simulation...")
    simulation.simulate(n_conf)

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()

