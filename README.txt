CSD-MD is a Python package that enables the user to setup and run a molecular
dynamics simulation from an entry in the Cambridge Structural Database.

Environment setup:

Note, CSD System must first be installed:
https://www.ccdc.cam.ac.uk/support-and-resources/csdsdownloads/

Quicker, easier and more robust using mamba instead of conda:
Mamba installation instructions: https://github.com/conda-forge/miniforge
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh ---> yes to initialise at the end (add executables to the path)

Example:
mamba create -n ccdc python=3.7.12
mamba activate ccdc
mamba install openmm openmmforcefields pyyaml
python -m pip install --extra-index-url https://pip.ccdc.cam.ac.uk/ csd-python-api

Running MD simulations:

A single .yaml input file is required. This contains all information
required to retrieve a structure from CSD and run an MD simulation.

For a small molecule simulation openmmforcefields is used to generate a
non-standard residue template. Charges assigned using AM1-BCC from Antechamber /
Ambertools.

Input options:
...to follow

Available force fields:
...to follow

System type:
...to follow

Simulation type:
...to follow

Example usage:
"python CSD-MD.py --md_params input.yaml > md.log"

Notes:
 - Problems finding the CSD Database so had to follow instructions here to set
 up mamba and save environment variables in ./etc/conda/activate.d/env_vars.sh:
 https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#set-env-vars

