class GetStructure():
    def __init__(self):
        pass

    def ligand(self, identifier):
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


    def protein(self, identifier):
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


