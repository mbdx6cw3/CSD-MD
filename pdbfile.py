"""
pdbfile.py: Used for loading PDB files.

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012-2021 Stanford University and the Authors.
Authors: Peter Eastman
Contributors:

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
from __future__ import print_function, division, absolute_import
__author__ = "Peter Eastman"
__version__ = "1.0"

from openmm import Vec3
from openmm.app.internal.pdbstructure import PdbStructure
from openmm.app import Topology
from openmm.unit import nanometers
from . import element as elem

class PDBFile(object):
    """PDBFile parses a Protein Data Bank (PDB) file and constructs a Topology and a set of atom positions from it.

    This class also provides methods for creating PDB files.  To write a file containing a single model, call
    writeFile().  You also can create files that contain multiple models.  To do this, first call writeHeader(),
    then writeModel() once for each model in the file, and finally writeFooter() to complete the file."""

    _residueNameReplacements = {}
    _atomNameReplacements = {}
    _standardResidues = ['ALA', 'ASN', 'CYS', 'GLU', 'HIS', 'LEU', 'MET', 'PRO', 'THR', 'TYR',
                         'ARG', 'ASP', 'GLN', 'GLY', 'ILE', 'LYS', 'PHE', 'SER', 'TRP', 'VAL',
                         'A', 'G', 'C', 'U', 'I', 'DA', 'DG', 'DC', 'DT', 'DI', 'HOH']

    def __init__(self, file, extraParticleIdentifier='EP'):
        """Load a PDB file.

        The atom positions and Topology can be retrieved by calling getPositions() and getTopology().

        Parameters
        ----------
        file : string or file
            the name of the file to load.  Alternatively you can pass an open file object.
        extraParticleIdentifier : string='EP'
            if this value appears in the element column for an ATOM record, the Atom's element will be set to None to mark it as an extra particle
        """
        
        metalElements = ['Al','As','Ba','Ca','Cd','Ce','Co','Cs','Cu','Dy','Fe','Gd','Hg','Ho','In','Ir','K','Li','Mg',
        'Mn','Mo','Na','Ni','Pb','Pd','Pt','Rb','Rh','Sm','Sr','Te','Tl','V','W','Yb','Zn']
        
        top = Topology()
        ## The Topology read from the PDB file
        self.topology = top

        # Load the PDB file

        if isinstance(file, PdbStructure):
            pdb = file
        else:
            inputfile = file
            own_handle = False
            if isinstance(file, str):
                inputfile = open(file)
                own_handle = True
            pdb = PdbStructure(inputfile, load_all_models=True, extraParticleIdentifier=extraParticleIdentifier)
            if own_handle:
                inputfile.close()
        PDBFile._loadNameReplacementTables()

        # Build the topology
        atomByNumber = {}
        for chain in pdb.iter_chains():
            c = top.addChain(chain.chain_id)
            for residue in chain.iter_residues():
                resName = residue.get_name()
                if resName in PDBFile._residueNameReplacements:
                    resName = PDBFile._residueNameReplacements[resName]
                r = top.addResidue(resName, c, str(residue.number), residue.insertion_code)
                if resName in PDBFile._atomNameReplacements:
                    atomReplacements = PDBFile._atomNameReplacements[resName]
                else:
                    atomReplacements = {}
                processedAtomNames = set()
                for atom in residue.atoms_by_name.values():
                    atomName = atom.get_name()
                    if atomName in processedAtomNames or atom.residue_name != residue.get_name():
                        continue
                    processedAtomNames.add(atomName)
                    if atomName in atomReplacements:
                        atomName = atomReplacements[atomName]
                    atomName = atomName.strip()
                    element = atom.element
                    if element == 'EP':
                        element = None
                    elif element is None:
                        # Try to guess the element.

                        upper = atomName.upper()
                        while len(upper) > 1 and upper[0].isdigit():
                            upper = upper[1:]
                        if upper.startswith('CL'):
                            element = elem.chlorine
                        elif upper.startswith('NA'):
                            element = elem.sodium
                        elif upper.startswith('MG'):
                            element = elem.magnesium
                        elif upper.startswith('BE'):
                            element = elem.beryllium
                        elif upper.startswith('LI'):
                            element = elem.lithium
                        elif upper.startswith('K'):
                            element = elem.potassium
                        elif upper.startswith('ZN'):
                            element = elem.zinc
                        elif len(residue) == 1 and upper.startswith('CA'):
                            element = elem.calcium
                        elif upper.startswith('D') and any(a.name == atomName[1:] for a in residue.iter_atoms()):
                            pass # A Drude particle
                        else:
                            try:
                                element = elem.get_by_symbol(upper[0])
                            except KeyError:
                                pass
                    newAtom = top.addAtom(atomName, element, r, str(atom.serial_number))
                    atomByNumber[atom.serial_number] = newAtom
        self._positions = []
        for model in pdb.iter_models(True):
            coords = []
            for chain in model.iter_chains():
                for residue in chain.iter_residues():
                    processedAtomNames = set()
                    for atom in residue.atoms_by_name.values():
                        if atom.get_name() in processedAtomNames or atom.residue_name != residue.get_name():
                            continue
                        processedAtomNames.add(atom.get_name())
                        pos = atom.get_position().value_in_unit(nanometers)
                        coords.append(Vec3(pos[0], pos[1], pos[2]))
            self._positions.append(coords*nanometers)
        ## The atom positions read from the PDB file.  If the file contains multiple frames, these are the positions in the first frame.
        self.positions = self._positions[0]
        self.topology.setPeriodicBoxVectors(pdb.get_periodic_box_vectors())
        self.topology.createStandardBonds()
        self.topology.createDisulfideBonds(self.positions)
        self._numpyPositions = None

        # Add bonds based on CONECT records. Bonds between metals of elements specified in metalElements and residues in standardResidues are not added.

        connectBonds = []
        for connect in pdb.models[-1].connects:
            i = connect[0]
            for j in connect[1:]:
                if i in atomByNumber and j in atomByNumber:    
                    if atomByNumber[i].element is not None and atomByNumber[j].element is not None:
                        if atomByNumber[i].element.symbol not in metalElements and atomByNumber[j].element.symbol not in metalElements:
                            connectBonds.append((atomByNumber[i], atomByNumber[j])) 
                        elif atomByNumber[i].element.symbol in metalElements and atomByNumber[j].residue.name not in PDBFile._standardResidues:
                            connectBonds.append((atomByNumber[i], atomByNumber[j])) 
                        elif atomByNumber[j].element.symbol in metalElements and atomByNumber[i].residue.name not in PDBFile._standardResidues:
                            connectBonds.append((atomByNumber[i], atomByNumber[j]))     
                    else:
                        connectBonds.append((atomByNumber[i], atomByNumber[j]))         
        if len(connectBonds) > 0:
            # Only add bonds that don't already exist.
            existingBonds = set(top.bonds())
            for bond in connectBonds:
                if bond not in existingBonds and (bond[1], bond[0]) not in existingBonds:
                    top.addBond(bond[0], bond[1])
                    existingBonds.add(bond)

