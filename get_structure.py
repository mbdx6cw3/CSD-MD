# TODO: this is needed on the CSF to reset environment variables after loading CSD module
# TODO: remove when CSF environment working
#import sys
#sys.path.insert(1, "/mnt/iusers01/rb01/mbdx6cw3/mambaforge/envs/ccdc_new/lib/python37.zip")
#sys.path.insert(1, "/mnt/iusers01/rb01/mbdx6cw3/mambaforge/envs/ccdc_new/lib/python3.7")
#sys.path.insert(1, "/mnt/iusers01/rb01/mbdx6cw3/mambaforge/envs/ccdc_new/lib/python3.7/lib-dynload")
#sys.path.insert(1, "/mnt/iusers01/rb01/mbdx6cw3/.local/lib/python3.7/site-packages")
#sys.path.insert(1, "/mnt/iusers01/rb01/mbdx6cw3/mambaforge/envs/ccdc_new/lib/python3.7/site-packages")


def ligand(identifier):
    """

    :return:
    """
    from ccdc.conformer import ConformerGenerator
    from ccdc.io import MoleculeWriter
    from ccdc.io import EntryReader

    csd_reader = EntryReader("CSD")
    entry = csd_reader.entry(identifier)
    ligand = entry.molecule
    conformer_generator = ConformerGenerator()
    conformers = conformer_generator.generate(ligand)

    with MoleculeWriter("input.pdb") as mol_writer:
        mol_writer.write(conformers[0].molecule)

    with MoleculeWriter("ligand.sdf") as mol_writer:
        mol_writer.write(conformers[0].molecule)

    return None


def protein(identifier):
    """

    :return:
    """
    from pdbfixer import PDBFixer
    from openmm import app

    fixer = PDBFixer(pdbid=identifier)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(True)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    app.PDBFile.writeFile(fixer.topology, fixer.positions, open('input.pdb', 'w'))
    return None

