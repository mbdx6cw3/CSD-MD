import sys
resetopenflags = sys.getdlopenflags()
from ccdc.molecule import Molecule as CCDCMolecule
sys.setdlopenflags(resetopenflags)

from openmm import Vec3
from openmm.app.topology import Topology as OpenMMTopology
from openmm.app.element import Element as OpenMMElement

from typing import Any


ANGSTROM_TO_NANOMETRE = 0.1


def openmm_topology_and_positions_from_ccdc_molecule(ccdc_mol: CCDCMolecule) -> tuple[OpenMMTopology, list[Vec3]]:
    """Converts a Molecule from the CSD Python API to OpenMM types.
    
    :param ccdc_mol: the CSD Python API Molecule to convert
    :return: a tuple containg an OpenMM Topology, and a list of Vec3 with the atomic positions (in nm)
    """
    topology = OpenMMTopology()
    positions = []

    chain = topology.addChain()
    residue = topology.addResidue('ligand', chain)

    ccdc_to_openmm_atom_map = dict()

    for ccdc_atom in ccdc_mol.atoms:
        ccdc_to_openmm_atom_map[ccdc_atom] = topology.addAtom(
            ccdc_atom.label, OpenMMElement.getByAtomicNumber(ccdc_atom.atomic_number), residue)
        positions.append(Vec3(*(coord * ANGSTROM_TO_NANOMETRE for coord in ccdc_atom.coordinates)))

    for ccdc_bond in ccdc_mol.bonds:
        topology.addBond(
            ccdc_to_openmm_atom_map[ccdc_bond.atoms[0]],
            ccdc_to_openmm_atom_map[ccdc_bond.atoms[1]],
        )

    return topology, positions
