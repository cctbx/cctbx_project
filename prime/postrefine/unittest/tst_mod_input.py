"""Unit test for mod_input.py"""
from __future__ import division
from prime.postrefine import mod_input
import unittest

class readBadInput(unittest.TestCase):
  def testMissingData(self):
    """mod_input should fail with invalid or missing data"""
    self.assertRaises(mod_input.InvalidData, mod_input.process_input, \
        ['dummy=dummier'], False)
  def testInvalidCrystalSystem(self):
    """mod_input should fail with invalid crystal system"""
    self.assertRaises(mod_input.InvalidCrystalSystem, mod_input.process_input, \
        ['data=test_data','target_crystal_system=dummy'], False)
  def testMissingNumberOfResidues(self):
    """mod_input should fail with invalid or missing n_residues"""
    self.assertRaises(mod_input.InvalidNumberOfResidues, mod_input.process_input, \
        ['data=test_data'], False)
  def testInvalidPixelSize(self):
    """mod_input should fail with invalid or missing pixel size"""
    self.assertRaises(mod_input.InvalidPixelSize, mod_input.process_input, \
        ['data=test_data','n_residues=100'], False)
  def testBadRunNo(self):
    """mod_input should fail with invalid run no."""
    self.assertRaises(mod_input.InvalidRunNo, mod_input.process_input, \
        ['data=test_data','n_residues=100','pixel_size_mm=0.07934','run_no=dummy_run_no'], True)
  def testNoData(self):
    """mod_input should fail when data is empty."""
    self.assertRaises(mod_input.InvalidData, mod_input.read_pickles, ['empty_data'])


if __name__ == "__main__":
    unittest.main()
