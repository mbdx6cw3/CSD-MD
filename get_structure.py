def ligand(identifier, simulation):
    """

    :return:
    """
    import sys
    resetopenflags = sys.getdlopenflags()
    from ccdc.conformer import ConformerGenerator
    from ccdc.io import EntryReader, MoleculeWriter
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
    simulation.conformers = [c.molecule for c in conformers]

    # TODO: solvated systems required ionised ligands
    # TODO: for now need to print pdb, ionise using OB and then write and fix pdb
    if simulation.solvate:
        for i_conf, c in enumerate(conformers):
            with MoleculeWriter(f"{simulation.input_dir}ligand.pdb") as mol_writer:
                mol_writer.write(c.molecule)
        simulation.smiles = ionize_smiles(simulation.smiles)
        ionize_pdb(f"{simulation.input_dir}ligand.pdb")

    if simulation.type == "enhanced":
        geometry_analysed_mol = get_torsions(conformers[0].molecule)
        simulation.unique_torsions = get_unique_torsions(geometry_analysed_mol)

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


def ionize_smiles(smiles):
    from openbabel import openbabel
    pH = 7.4
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("smi", "smi")
    mol = openbabel.OBMol()
    obConversion.ReadString(mol, smiles)
    mol.CorrectForPH(pH)
    smiles = obConversion.WriteString(mol)
    return smiles


def ionize_pdb(filename):
    from openbabel import openbabel
    pH = 7.4
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("pdb", "pdb")
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, filename)
    mol.CorrectForPH(pH)
    mol.AddHydrogens()
    obConversion.WriteFile(mol, filename)
    fix_ob_output(filename)
    return None


def fix_ob_output(filename):
    '''
    fix atom labelling
    ---------
    '''
    with open(filename) as infile:
        # count element in a dict
        ele_count = {}
        text = ""
        for pdb_line in infile:
            line = pdb_line.strip().split()
            if line[0] == "HETATM":
                ele = line[2]
                try:
                    # rename if more than one (add count)
                    ele_count[ele] += 1
                except KeyError:
                    ele_count[ele] = 1
                atom_name = ele + str(ele_count[ele])
                pdb_line = pdb_line.replace(f"{ele}  ", f"{atom_name:<3}")
            text = text + pdb_line
    with open(filename, "w") as outfile:
        outfile.write(text)
    return None