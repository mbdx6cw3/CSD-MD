CSD-MD is a Python package that enables the user to setup and run a molecular
dynamics simulation using CSD entries and tools.

Installation:
module load apps/binapps/anaconda3/2023.09
conda install -n base conda-libmamba-solver
conda config --set solver libmamba
conda create -n csd-md python=3.9
conda activate csd-md
conda install openmm openmmforcefields pyyaml pdbfixer openmm-plumed openbabel

...then on CSF (Red Hat linux)
conda install -c /opt/apps/apps/binapps/ccdc-csds/2024.1/ccdc_conda_channel_py39 Pillow six lxml numpy matplotlib
conda install -c /opt/apps/apps/binapps/ccdc-csds/2024.1/ccdc_conda_channel_py39 csd-python-api
export CCDC_LICENSING_CONFIGURATION='lf-server;http://login1:8090'

...or on local machine (OS(X)):
python -m pip install --extra-index-url https://pip.ccdc.cam.ac.uk/ csd-python-api

Running MD simulations:
A single .yaml input file is required. This contains all information
required to retrieve structures, construct topologies and run MD simulations.

Example usage:
    "python CSD-MD.py --md_params input.yaml > md.log"

Input options:
    name:                   name of the simulation
    system type:            "ligand", "protein" or "protein-ligand"
    CSD identifier:         identifier associated with the ligand
    PDB identifier:         four letter identifier associated with protein,
    pair-net model path:    path to trained PairNet model or "none"
    solvate system:         "yes" or "no"
    simulation type:        "standard" or "enhanced"
    simulation time (ns):   total simulation time in nanoseconds
    timestep (fs):          integration timestep in femtoseconds
    temperature (K):        temperature in Kelvin
    ensemble:               "NVT"

--------------------------------------------------------------------------------
System types:

"Ligand" will retrieve a ligand from a CSD entry and generate the initial
structure using CCDC conformer generator.

"Protein" will retrieve a protein from RCSB Protein Data Bank and
generate the initial (sanitised) structure using PDBFixer.

"Ligand-protein" will retrieve a ligand from a CSD entry and a protein from
RCSB Protein Data Bank, and then generate the initial structure by docking the
ligand to the protein, defining the binding site using a native ligand in the
unsanitised protein structure.

--------------------------------------------------------------------------------
Solvate system:
"yes" (ligand only) adds water to the system and ionises functional groups
appropriate for pH 7.4.
"no" will perform a gas phase simulation
Note that since PairNet has a fixed number of input descriptors, the number of
atoms in the ligand must match the number of atoms in the PairNet model.

--------------------------------------------------------------------------------
Simulation types:
"standard" will perform an MD simulation
"enhanced" will perform a metadynamics simulation with sampling enhanced
with respect to the rotatable bonds identified using the CCDC conformer
generator

--------------------------------------------------------------------------------
PairNet Model Library (in "models" directory):
models/aspirin/neutral/MD-300K/
models/aspirin/neutral/Meta-300K/
models/aspirin/ionised/MD-300K/
models/aspirin/ionised/Meta-300K/
models/ibuprofen/neutral/MD-300K/
models/ibuprofen/neutral/Meta-300K/
models/ibuprofen/ionised/MD-300K/
models/ibuprofen/ionised/Meta-300K/

Note that "none" will use an MM potential (GAFF2) instead of PairNet. Water is
modelled using TIP3P.

--------------------------------------------------------------------------------
Example Library (in "examples" directory):
asp-gas-MM.yaml:                MD simulation of aspirin using an MM potential
asp-solution-MM.yaml:           MD simulation of aspirin in water using an MM potential
asp-gas-MM-enhanced.yaml:       Metaydnamics simulation of aspirin using an MM potential
asp-gas-ML.yaml:                MD simulation of aspirin using a PairNet potential
asp-solution-ML.yaml:           MD simulation of aspirin in water using a PairNet potential
ibu-gas-ML.yaml:                MD simulation of ibuprofen using a PairNet potential
ibu-solution-MM-enhanced.yaml:  Metadynamics simulation of ibuprofen in water using an MM potential
4ph9-protein-MM.yaml:           MD simulation of cyclooxygenase-2 using an MM potential
asp-4ph9-MM.yaml:               MD simulation of cyclooxygenase-2 bound aspirin using an MM potential
asp-4ph9-ML.yaml:               MD simulation of cyclooxygenase-2 bound aspirin using a PairNet potential
ibu-4ph9-MM.yaml:               MD simulation of cyclooxygenase-2 bound ibuprofen using a PairNet potential

 -------------------------------------------------------------------------------
References:
[1] CR Groom, IJ Bruno, MP Lightfoot and SC Ward, The Cambridge Structural
    Database, 2016, Acta Cryst. B72: 171-179.

[2] JC Cole, O Korb, P McCabe, MG Read, R Taylor, Knowledge-Based Conformer
    Generation Using the Cambridge Structural Database, 2018, J. Chem. Inf.
    Model. 58: 615-629.

[3] G Jones, P Willett, RC Glen, AR Leach, R Taylor, Development and
    Validation of a Genetic Algorithm for Flexible Docking, 1997, J. Mol.
    Bio. 267: 727-748.

[4] P Eastman, J Swails, JD Chodera, RT McGibbon, Y Zhao, KA Beauchamp,
    LP Wang, AC Simmonett, MP Harrigan, CD Stern, RP Wiewiora, BR Brooks,
    VS Pande, OpenMM 7: Rapid Development of High Performance Algorithms for
    Molecular Dynamics, 2017, PLOS Comp. Biol. 13(7): e1005659.

[5] CD Williams, J Kalayan, NA Burton, RA Bryce, Stable and Accurate
    Atomistic Simulations of Flexible Molecules using Conformationally
    Generalisable Machine Learned Potentials, 2024, Chem. Sci., 15: 12780-12795.

[6] RA Sykes, NT Johnson, CJ Kingsbury et al, What Has Scripting Ever Done For Us?
    The CSD Python Application Programming Interface (API), J. Appl. Cryst., 2024,
    57, 1235-1250.
