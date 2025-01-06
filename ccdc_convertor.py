import sys
resetopenflags = sys.getdlopenflags()
from ccdc.molecule import Molecule as CCDCMolecule, Atom as CCDCAtom, Bond as CCDCBond
sys.setdlopenflags(resetopenflags)

from openmm import Vec3
from openmm.app.topology import Topology as OpenMMTopology
from openmm.app.element import Element as OpenMMElement

from openbabel import OBMol as OpenBabelMolecule

from typing import Any


ANGSTROM_TO_NANOMETRE = 0.1


def openmm_topology_and_positions_from_ccdc_molecule(ccdc_mol: CCDCMolecule) -> tuple[OpenMMTopology, list[Vec3]]:
    """Converts a Molecule from the CSD Python API to OpenMM types.

    Retains protein residue information if present.

    :param ccdc_mol: the CSD Python API Molecule to convert
    :return: a tuple containg an OpenMM Topology, and a list of Vec3 with the atomic positions (in nm)
    """
    topology = OpenMMTopology()
    positions = []

    chains = dict()
    residues_by_chain = dict()

    ccdc_to_openmm_atom_map = dict()

    for ccdc_atom in ccdc_mol.atoms:
        chain_id = ccdc_atom.chain_label or ''
        if chain_id not in chains:
            chains[chain_id] = topology.addChain(chain_id)
            residues_by_chain[chain_id] = dict()
        chain = chains[chain_id]
        
        residue_id = ccdc_atom.residue_label or 'UNK1'
        if residue_id not in residues_by_chain[chain_id]:
            residue_name = residue_id[:3]
            residues_by_chain[chain_id][residue_id] = topology.addResidue(residue_name, chain)              
        residue = residues_by_chain[chain_id][residue_id]


        ccdc_to_openmm_atom_map[ccdc_atom] = topology.addAtom(
            ccdc_atom.label, OpenMMElement.getByAtomicNumber(ccdc_atom.atomic_number), residue)
        positions.append(Vec3(*(coord * ANGSTROM_TO_NANOMETRE for coord in ccdc_atom.coordinates)))

    for ccdc_bond in ccdc_mol.bonds:
        topology.addBond(
            ccdc_to_openmm_atom_map[ccdc_bond.atoms[0]],
            ccdc_to_openmm_atom_map[ccdc_bond.atoms[1]],
        )

    return topology, positions


def openbabel_molecule_from_ccdc_molecule(ccdc_mol: CCDCMolecule) -> OpenBabelMolecule:
    """Converts a molecule from CSD Python API to OpenBabel types.

    :param ccdc_mol: the CSD Python API Molecule to convert
    :return: an openbabel.OBMol object
    """
    ob_mol = OpenBabelMolecule()
    ob_mol.SetAromaticPerceived()   # tell openbabel we are defining aromaticity ourselves

    ccdc_to_ob_atom_map = dict()

    for ccdc_atom in ccdc_mol.atoms:
        ob_atom = ob_mol.NewAtom()
        ob_atom.SetTitle(ccdc_atom.label)
        ob_atom.SetAtomicNum(ccdc_atom.atomic_number)
        ob_atom.SetVector(*ccdc_atom.coordinates)
        ccdc_to_ob_atom_map[ccdc_atom] = ob_atom

    for ccdc_bond in ccdc_mol.bonds:
        ob_bond = ob_mol.NewBond()
        begin_atom = ccdc_to_ob_atom_map[ccdc_bond.atoms[0]]
        end_atom = ccdc_to_ob_atom_map[ccdc_bond.atoms[1]]
        ob_bond.SetBegin(begin_atom)
        ob_bond.SetEnd(end_atom)
        begin_atom.AddBond(ob_bond)
        end_atom.AddBond(ob_bond)

        if ccdc_bond.bond_type == 'Single':
            ob_bond.SetBondOrder(1)
        elif ccdc_bond.bond_type == 'Double':
            ob_bond.SetBondOrder(2)
        elif ccdc_bond.bond_type == 'Triple':
            ob_bond.SetBondOrder(3)
        elif ccdc_bond.bond_type == 'Aromatic':
            ob_bond.SetBondOrder(5)
            ob_bond.SetAromatic()

    return ob_mol


def ccdc_molecule_from_openbabel_molecule(ob_mol: OpenBabelMolecule) -> CCDCMolecule:
    """Converts a molecule from OpenBabel to CSD Python API types.

    :param ccdc_mol: the openbabel.OBMol object to convert
    :return: a CSD Python API Molecule
    """
    ccdc_mol = CCDCMolecule()

    for i in range(ob_mol.NumAtoms()):
        ob_atom  = ob_mol.GetAtom(i + 1)
        ccdc_atom = CCDCAtom(
            atomic_number=ob_atom.GetAtomicNum(),
            coordinates=(ob_atom.x(), ob_atom.y(), ob_atom.z())
        )
        ccdc_mol.add_atom(ccdc_atom)

    for i in range(ob_mol.NumBonds()):
        ob_bond  = ob_mol.GetBond(i)

        ccdc_bond_type = 'Unknown'
        if ob_bond.IsAromatic():
            ccdc_bond_type = 'Aromatic'
        elif ob_bond.GetBondOrder() == 1:
            ccdc_bond_type = 'Single'
        elif ob_bond.GetBondOrder() == 2:
            ccdc_bond_type = 'Double'
        elif ob_bond.GetBondOrder() == 3:
            ccdc_bond_type = 'Triple'

        ccdc_mol.add_bond(bond_type=CCDCBond.BondType(ccdc_bond_type), atom1=ccdc_mol.atoms[ob_bond.GetBeginAtomIdx() - 1], atom2=ccdc_mol.atoms[ob_bond.GetEndAtomIdx() - 1])
    
    return ccdc_mol
