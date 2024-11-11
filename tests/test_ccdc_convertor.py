import unittest

# Add the location of the script being tested to the path
import sys
from os.path import dirname, abspath
sys.path.append(dirname(dirname(abspath(__file__))) )

import ccdc_convertor

from ccdc.io import EntryReader
from openmm.app.element import Element as OpenMMElement


class TestCCDCConvertor(unittest.TestCase):
    def test_openmm_from_ccdc_molecule(self):

        reader = EntryReader('csd')
        entry = reader.entry('ACSALA')
        ccdc_molecule = entry.molecule

        openmm_topology, openmm_positions = ccdc_convertor.openmm_topology_and_positions_from_ccdc_molecule(ccdc_molecule)

        self.assertEqual(21, openmm_topology.getNumAtoms())
        self.assertEqual(21, openmm_topology.getNumBonds())

        atoms = [atom for atom in openmm_topology.atoms()]
        bonds = [bond for bond in openmm_topology.bonds()]

        self.assertEqual(atoms[0].name, 'C1')
        self.assertEqual(atoms[0].element, OpenMMElement.getBySymbol('C'))
        self.assertAlmostEqual(0.1678, openmm_positions[0][0], places=4)
        self.assertAlmostEqual(0.2879, openmm_positions[0][1], places=4)
        self.assertAlmostEqual(0.0764, openmm_positions[0][2], places=4)

        self.assertEqual(atoms[17].name, 'O1')
        self.assertEqual(atoms[17].element, OpenMMElement.getBySymbol('O'))
        self.assertAlmostEqual(0.0015, openmm_positions[17][0], places=4)
        self.assertAlmostEqual(0.1237, openmm_positions[17][1], places=4)
        self.assertAlmostEqual(0.1095, openmm_positions[17][2], places=4)

        self.assertEqual((atoms[0], atoms[1]), bonds[0])    # C1-C2
        self.assertEqual((atoms[17], atoms[6]), bonds[16])  # O1-C7
