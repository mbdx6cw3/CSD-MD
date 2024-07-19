CSD-MD enables the user to setup and run a molecular dynamics simulation from
an entry in the Cambridge Structural Database.

Required packages:
    OpenMM (e.g. conda install conda-forge::openmm)
    CCDC

A single input file is required in .yaml format. This contains all information
required to run an MD simulation.

Input options:
...to follow

Available force fields:
...to follow

Example usage:
"python CSD-MD.py --md_params input.yaml > md.log"
