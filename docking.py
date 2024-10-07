from ccdc.io import MoleculeReader, EntryWriter
from ccdc.docking import Docker
from ccdc.diagram import DiagramGenerator

protein_file = "protein.pdb"
ligand_file = "ligand.pdb"
input_file = "ligand.sdf"

radius = 6
ndocks = 20
fitness_function = 'plp'
autoscale = 30
write_options = ['NO_LINK_FILES', 'NO_RNK_FILES', 'NO_GOLD_LIGAND_MOL2_FILE']
save_lone_pairs = False

diagram_generator = DiagramGenerator()

diagram_generator.settings.return_type = 'SVG'
diagram_generator.settings.explicit_polar_hydrogens = False
diagram_generator.settings.shrink_symbols = False

settings = Docker.Settings()
settings.add_protein_file(str(protein_file))
native_ligand = MoleculeReader(str(ligand_file))[0]
settings.binding_site = settings.BindingSiteFromLigand(settings.proteins[0], native_ligand, radius)
settings.add_ligand_file(str(input_file), ndocks=ndocks)

#settings.output_format = output_format
settings.fitness_function = fitness_function
settings.autoscale = autoscale
settings.write_options = write_options
settings.save_lone_pairs = save_lone_pairs

docker = Docker(settings=settings)
results = docker.dock(mode='foreground', file_name='api_gold.conf')
assert results.return_code == 0, "Error! GOLD did not run successfully."
soln = results.ligands[0]
soln.fitness()
soln.scoring_term()

export_format = 'pdb'  # File format in which to export protein-ligand complexes

complexed = results.make_complex(soln)
complexed.remove_unknown_atoms()  # Remove lone pairs for export

file_path = f'complexed.{export_format}'

with EntryWriter(file_path) as writer:
    writer.write(complexed)

