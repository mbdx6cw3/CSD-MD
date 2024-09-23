CSD-MD is a Python package that enables the user to setup and run a molecular
dynamics simulation from an entry in the Cambridge Structural Database.

Installation:

Note, CSD System must first be installed:
https://www.ccdc.cam.ac.uk/support-and-resources/csdsdownloads/

Quicker, easier and more robust using mamba instead of conda:
Mamba installation instructions: https://github.com/conda-forge/miniforge
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh ---> yes to initialise at the end (add executables to the path)

Example environment setup:
mamba create -n ccdc python=3.9
mamba activate ccdc
mamba install openmm openmmforcefields pyyaml tensorflow=2.12.0 pdbfixer
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
    pair-net model: none
    solvate system: yes
    simulation type: standard
    simulation time (ns): 0.01
    timestep (fs): 1.0
    temperature (K): 300.0
    ensemble: NVT

Input options:
    name:                   name of the simulation

    system type:            "ligand", "protein" or "protein-ligand"

    CSD identifier:         identifier associated with CSD ligand entry

    PDB identifier:         four letter identifier associated with PDB entry

    pair-net model:         name of trained pair-net model to search for in
                            "pair-net_models" directory. "none" will default
                            to Amber potential. More MM potentials to be added

    solvate system:         "yes" will fill box with explicit waters, modelled
                            by default using TIP3P

    simulation type:        "standard" is the only option for now

    simulation time (ns):   total simulation time in ns

    timestep (fs):          integration timestep in femtoseconds

    temperature (K):        temperature in Kelvin

    ensemble:               "NVT"

Notes:
For a simulation involving a ligand OpenFF is used to generate a non-standard
residue template with AM1-BCC charges assigned.

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

