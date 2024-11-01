def ligand(identifier, simulation):
    """

    :return:
    """
    import sys
    resetopenflags = sys.getdlopenflags()
    from ccdc.conformer import ConformerGenerator
    from ccdc.io import MoleculeWriter
    from ccdc.io import EntryReader
    # this is required to reset paths after loading ccdc modules
    sys.setdlopenflags(resetopenflags)

    csd_reader = EntryReader("CSD")
    entry = csd_reader.entry(identifier)
    ligand = entry.molecule
    simulation.smiles = ligand.smiles
    conformer_generator = ConformerGenerator()
    if simulation.type != "multi-conformer":
        conformer_generator.settings.max_conformers = 1

    print("Generating conformers...")
    conformers = conformer_generator.generate(ligand)
    simulation.n_conf = len(conformers)
    for i_conf, c in enumerate(conformers):
        with MoleculeWriter(f"{simulation.input_dir}/ligand_{i_conf}.pdb") as mol_writer:
            mol_writer.write(c.molecule)

    if simulation.type == "enhanced":
        geometry_analysed_mol = get_torsions(conformers[0].molecule)
        simulation.unique_torsions = get_unique_torsions(geometry_analysed_mol)

    # TODO: CCDC MoleculeWriter does not print MODEL or ENDMDL
    # TODO: these keywords are used to look at different structures of the same molecule
    # TODO: at present PDBFile sees only one structure with n_conf aspirins
    # TODO: this problem will go away with conversion functions
    if simulation.analyse_torsions:
        print("Torsion Angle | Atom Indices")
        geometry_analysed_mol = get_torsions(conformers[0].molecule)
        for count, tor in enumerate(geometry_analysed_mol.analysed_torsions):
            print(f"{count:12}  | {tor.atom_indices}")
        print()
        print("Conformer | Probability Score | RMSD wrt Input | Initial Torsion Angles")
        print("-----------------------------------------------------------------------")

        for i_conf, c in enumerate(conformers):
            all_torsions = get_torsions(c.molecule).analysed_torsions
            torsions = []
            for count, tor in enumerate(all_torsions):
                torsions.append(round(tor.value, 1))
            print(f"{i_conf:9d} | {c.normalised_score:17.3f} | {c.rmsd():14.3f} | {(' '.join(str(x) for x in torsions))}")
        print("----------------------------------------------")

    return None


def protein(identifier, simulation):
    from openmm import app
    from pdbfixer import PDBFixer
    if identifier == "from_file":
        fixer = PDBFixer(filename="input.pdb")
    else:
        fixer = PDBFixer(pdbid=identifier)
    # will also need unfixed structure for docking if protein-ligand simulation
    if simulation.ligand:
        app.PDBFile.writeFile(fixer.topology, fixer.positions,
            open(f"{simulation.input_dir}/protein-unfixed.pdb", "w"))
    print("Fixing PDB file...")
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    app.PDBFile.writeFile(fixer.topology, fixer.positions,
        open(f"{simulation.input_dir}/protein.pdb", "w"))
    return None


def docking():
    # TODO: read in protein-unfixed and ligand pdbs from md_input
    # TODO: here we will do some protein-ligand docking to generate initial structure
    # TODO: will need to use the unfixed pdb file because we need the native ligand to define binding site
    # TODO: then tidy up protein structure and recombine with ligand
    # TODO: write new PDB file for the ligand and protein
    return None


def fix_pdb(identifier, simulation):
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
    app.PDBFile.writeFile(fixer.topology, fixer.positions,
        open(f"{simulation.input_dir}/protein.pdb", "w"))
    return None


def get_torsions(molecule):
    from ccdc.conformer import GeometryAnalyser
    engine = GeometryAnalyser()
    geometry_analysed_mol = engine.analyse_molecule(molecule)
    return geometry_analysed_mol


def get_unique_torsions(molecule):
    all_torsions = molecule.analysed_torsions
    unique_torsions = []
    for count, tor in enumerate(all_torsions):
        if count == 0:
            unique_torsions.append(all_torsions[count].atom_indices)
            continue
        delete_torsion = False
        i_pair = all_torsions[count].atom_indices[1:-1]
        for i_tors in range(len(unique_torsions)):
            j_pair = unique_torsions[i_tors][1:-1]
            if (i_pair == j_pair):
                delete_torsion = True
        if delete_torsion:
            continue
        else:
            unique_torsions.append(all_torsions[count].atom_indices)
    return unique_torsions

