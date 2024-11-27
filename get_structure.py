import os

def ligand(identifier, simulation):
    """

    :return:
    """
    import sys
    resetopenflags = sys.getdlopenflags()
    from ccdc.conformer import ConformerGenerator
    from ccdc.io import EntryReader, MoleculeWriter
    # this is required to reset paths due to bug in CSD-Python-API
    sys.setdlopenflags(resetopenflags)

    csd_reader = EntryReader("CSD")
    entry = csd_reader.entry(identifier)
    ligand = entry.molecule
    simulation.smiles = ligand.smiles
    conformer_generator = ConformerGenerator()
    '''
    if simulation.type != "multi-conformer":
        conformer_generator.settings.max_conformers = 1
    '''

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


def get_protein(simulation):
    import sys
    import urllib.request
    resetopenflags = sys.getdlopenflags()
    sys.setdlopenflags(resetopenflags)
    URL = f"https://files.rcsb.org/download/{simulation.PDB}.pdb"
    file_path = f"{simulation.input_dir}protein-unsanitised.pdb"
    urllib.request.urlretrieve(URL, file_path)
    return None


def fix_protein_pdbfixer(simulation):
    from openmm import app
    from pdbfixer import PDBFixer
    file_path = f"{simulation.input_dir}protein-unsanitised.pdb"
    fixer = PDBFixer(filename=file_path)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    app.PDBFile.writeFile(fixer.topology, fixer.positions,
        open(f"{simulation.input_dir}protein.pdb", "w"))
    return None


def fix_protein_ccdcfixer(simulation):
    import sys
    # this is required to reset paths due to bug in CSD-Python-API
    resetopenflags = sys.getdlopenflags()
    from ccdc.protein import Protein
    from ccdc.io import MoleculeWriter
    sys.setdlopenflags(resetopenflags)
    simulation.protein = Protein.from_file(f"{simulation.input_dir}protein-unsanitised.pdb")
    simulation.protein.remove_all_metals()
    simulation.protein.remove_all_waters()
    for ligand in simulation.protein.ligands:
        simulation.protein.remove_ligand(ligand.identifier)
    simulation.protein.add_hydrogens(mode="all", rules_file=None)
    simulation.protein.sort_atoms_by_residue()
    with MoleculeWriter(f"{simulation.input_dir}protein.pdb") as w:
        w.write(simulation.protein)
    return None


def docking(simulation):
    import os, sys, shutil
    from pathlib import Path
    resetopenflags = sys.getdlopenflags()
    from ccdc.docking import Docker
    from ccdc.io import MoleculeReader, EntryWriter
    sys.setdlopenflags(resetopenflags)

    docker = Docker()
    settings = docker.settings

    input_dir = Path("docking_input").absolute()

    # get the protein structure and load into docker settings
    protein_file = input_dir / "protein.pdb"
    settings.add_protein_file(str(protein_file))

    # load the native ligand and use it to define the binding site
    native_ligand_file = input_dir / "native_ligand.pdb"
    native_ligand = MoleculeReader(str(native_ligand_file))[0]
    protein = settings.proteins[0]

    # this defines a binding site as within a radius of 6A from the ligand
    settings.binding_site = settings.BindingSiteFromLigand(protein, native_ligand)

    # define docking parameters including fitness function
    settings.fitness_function = "plp"
    settings.autoscale = 10.
    settings.early_termination = False

    settings.output_file = "docked_ligands.mol2"

    # select ligands to dock and number of docking runs per ligand
    ligand_file = input_dir / "ligand.pdb"
    settings.add_ligand_file(str(ligand_file), 10)

    # create output directory
    output_dir = Path("docking_output")
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir()
    os.chdir(output_dir)

    # launch the docker
    results = docker.dock()
    print(results.return_code)
    batch_conf_file = settings.conf_file

    settings = Docker.Settings.from_file(batch_conf_file)
    results = Docker.Results(settings)

    ligands = results.ligands
    first_dock = ligands[0]
    "Gold.PLP.Fitness" in first_dock.attributes

    print(first_dock.fitness())
    print(first_dock.fitness(settings.fitness_function))
    print(first_dock.scoring_term())
    print(first_dock.scoring_term("plp", "chemscore", "hbond"))
    print(first_dock.hbonds())
    complexed = results.make_complex(ligands[0])
    complexed.remove_unknown_atoms()
    print(len(protein.ligands))

    with EntryWriter("complexed.pdb") as writer:
        writer.write(complexed)

    # TODO: CCDC complex - to be converted to OpenMM object
    simulation.complexed = complexed
    os.chdir("..")

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