class GetStructure():
    def __init__(self):
        pass

    def ligand(self, identifier, simulation):
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
        if simulation.simulation_type != "multi-conformer":
            conformer_generator.settings.max_conformers = 1
        conformers = conformer_generator.generate(ligand)

        # TODO: CCDC MoleculeWriter does not print MODEL or ENDMDL
        # TODO: these keywords are used to look at different structures of the same molecule
        # TODO: at present PDBFile sees only one structure with n_conf aspirins
        # TODO: this problem will go away with conversion functions
        i_conf = 0
        for c in conformers:
            with MoleculeWriter(f"{simulation.input_dir}/ligand_{i_conf}.pdb") as mol_writer:
                mol_writer.write(c.molecule)
            i_conf = i_conf + 1
        n_conf = i_conf

        #with MoleculeWriter("input.pdb") as mol_writer:
        #    mol_writer.write(conformers[0].molecule)

        with MoleculeWriter(f"{simulation.input_dir}/ligand.sdf") as mol_writer:
            mol_writer.write(conformers[0].molecule)

        return n_conf


    def protein(self, identifier, simulation):
        """

        :return:
        """
        from pdbfixer import PDBFixer
        from openmm import app
        fixer = PDBFixer(pdbid=identifier)
        print("Fixing PDB file...")
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.removeHeterogens(False)
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.0)
        # TODO: add solvent box here instead of using OpenMM modeller?
        app.PDBFile.writeFile(fixer.topology, fixer.positions,
                              open(f"{simulation.input_dir}/protein.pdb", "w"))
        return None


    def docking(self):
        return None



