class GetStructure():
    def __init__(self):
        pass

    def ligand(self, identifier, simulation):
        """

        :return:
        """
        from ccdc.conformer import ConformerGenerator, GeometryAnalyser
        from ccdc.io import MoleculeWriter
        from ccdc.io import EntryReader

        csd_reader = EntryReader("CSD")
        entry = csd_reader.entry(identifier)
        ligand = entry.molecule
        simulation.smiles = ligand.smiles
        conformer_generator = ConformerGenerator()
        if simulation.simulation_type != "multi-conformer":
            conformer_generator.settings.max_conformers = 1

        print("Generating conformers...")
        conformers = conformer_generator.generate(ligand)
        simulation.n_conf = len(conformers)
        for i_conf, c in enumerate(conformers):
            with MoleculeWriter(f"{simulation.input_dir}/ligand_{i_conf}.pdb") as mol_writer:
                mol_writer.write(c.molecule)


        # TODO: CCDC MoleculeWriter does not print MODEL or ENDMDL
        # TODO: these keywords are used to look at different structures of the same molecule
        # TODO: at present PDBFile sees only one structure with n_conf aspirins
        # TODO: this problem will go away with conversion functions
        if simulation.analyse_torsions:
            print("Torsion Angle | Atom Indices")
            engine = GeometryAnalyser()
            geometry_analysed_mol = engine.analyse_molecule(conformers[0].molecule)
            for count, tor in enumerate(geometry_analysed_mol.analysed_torsions):
                print(f"{count:12}: {(', '.join(label for label in tor.atom_labels))}")
            print()
            print("Conformer | Probability Score | RMSD wrt Input | Initial Torsion Angles")
            print("-----------------------------------------------------------------------")

            for i_conf, c in enumerate(conformers):
                engine = GeometryAnalyser()
                all_torsions = engine.analyse_molecule(c.molecule).analysed_torsions
                torsions = []
                for count, tor in enumerate(all_torsions):
                    torsions.append(round(tor.value, 1))
                print(f"{i_conf:9d} | {c.normalised_score:17.3f} | {c.rmsd():14.3f} | {(' '.join(str(x) for x in torsions))}")
            print("----------------------------------------------")

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

