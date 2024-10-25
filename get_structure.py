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
        simulation.smiles = ligand.smiles
        conformer_generator = ConformerGenerator()
        if simulation.simulation_type != "multi-conformer":
            conformer_generator.settings.max_conformers = 1
        # TODO: deprecated???
        # conformer_generator.settings.normalised_score_threshold = 0.1
        # print(conformer_generator.settings.normalised_score_threshold)

        print("Generating conformers...")
        conformers = conformer_generator.generate(ligand)

        # TODO: CCDC MoleculeWriter does not print MODEL or ENDMDL
        # TODO: these keywords are used to look at different structures of the same molecule
        # TODO: at present PDBFile sees only one structure with n_conf aspirins
        # TODO: this problem will go away with conversion functions
        i_conf = 0
        print("Conformer | Probability Score | RMSD wrt Input")
        print("----------------------------------------------")
        for c in conformers:
            with MoleculeWriter(f"{simulation.input_dir}/ligand_{i_conf}.pdb") as mol_writer:
                mol_writer.write(c.molecule)
                print(f"{i_conf:9d} | {c.normalised_score:17.3f} | {c.rmsd():14.3f}")
            i_conf = i_conf + 1
        print("----------------------------------------------")
        simulation.n_conf = i_conf

        return None


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



