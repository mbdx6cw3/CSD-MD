CSD-MD is a Python package that enables the user to setup and run a molecular
dynamics simulation from an entry in the Cambridge Structural Database.

Installation:
module load apps/binapps/anaconda3/2023.09
conda install -n base conda-libmamba-solver
conda config --set solver libmamba
conda create -n csd-md python=3.9
conda activate csd-md
conda install openmm openmmforcefields pyyaml pdbfixer openmm-plumed

...then on CSF (Red Hat linux)
conda install -c /opt/apps/apps/binapps/ccdc-csds/2024.1/ccdc_conda_channel_py39 Pillow six lxml numpy matplotlib
conda install -c /opt/apps/apps/binapps/ccdc-csds/2024.1/ccdc_conda_channel_py39 csd-python-api
export CCDC_LICENSING_CONFIGURATION='lf-server;http://login1:8090'

...or on local machine (OS(X)):
python -m pip install --extra-index-url https://pip.ccdc.cam.ac.uk/ csd-python-api

Notes:
 Problems finding the CSD Database so had to follow instructions here to set
 up mamba and save environment variables in ./etc/conda/activate.d/env_vars.sh:
 https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#set-env-vars

Running MD simulations:
A single .yaml input file is required. This contains all information
required to retrieve a structure from CSD and run an MD simulation.

Example usage:
    "python CSD-MD.py --md_params input.yaml > md.log"

Example .yaml input file (input.yaml by default):
    name: demonstration
    system type: ligand
    CSD identifier: ACSALA
    PDB identifier: 3I40
    analyse torsions: False
    simulation type: standard
    pair-net model path: none
    charges: am1-bcc
    solvate system: no
    ensemble: NVT
    temperature (K): 300.0
    simulation time (ns): 0.001
    timestep (fs): 1.0

Input options:
    name:                   name of the simulation

    system type:            "ligand", "protein" or "protein-ligand"

    CSD identifier:         identifier associated with CSD ligand entry

    PDB identifier:         four letter identifier associated with PDB entry

    pair-net model path:    name of trained pair-net model to search for in
                            "pair-net_models" directory. If set to "none"
                             will default to MM simulation using the Amber14 potential.
                             Water is always TIP3P.

    solvate system:         "yes" will fill box with explicit water molecules,
                            modelled by default using TIP3P

    simulation type:        "standard" or "multi-conformer" for now

    simulation time (ns):   total simulation time (per conformer) in nanoseconds

    timestep (fs):          integration timestep in femtoseconds

    temperature (K):        temperature in Kelvin

    ensemble:               "NVT"

Notes:
For a simulation involving a ligand OpenFF is used to generate a non-standard
residue template with AM1-BCC charges assigned.

Only neutral molecules?

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

