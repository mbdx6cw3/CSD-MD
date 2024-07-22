from openmm.app import *
#from ccdc.conformer import ConformerGenerator
#from ccdc.io import MoleculeWriter
#from ccdc.io import EntryReader

class CSDDatabase:

    def smallmolecule(self, identifier):
        """

        :return:
        """

        '''
        csd_reader = EntryReader("CSD")
        entry = csd_reader.entry("ACSALA")
        molecule = entry.molecule
        conformer_generator = ConformerGenerator()
        conformers = conformer_generator.generate(molecule)
        mol_writer = MoleculeWriter("conformers.pdb")
        mol_writer.write(conformers[0].molecule)
        '''

        return


def pdb():
    """

    :return:
    """
    pdb = PDBFile('input.pdb')
    return pdb