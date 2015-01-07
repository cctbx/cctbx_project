from __future__ import division
from cctbx.geometry_restraints.clash_score import check_and_add_hydrogen
from libtbx.utils import null_out
from iotbx.pdb import fetch
import iotbx.pdb
import unittest
import os

__author__ = 'Youval'


class MyTestCase(unittest.TestCase):

  def setUp(self):
    self.file_to_delete = []
    # import files used in tests
    fn = '1a18'
    # fetch pdb file
    self.file_name = fetch.get_pdb (fn,'pdb',mirror='rcsb',log=null_out())
    self.file_to_delete.append(self.file_name)

  def test_identifying_and_addition_of_hydrogen(self):
    """ test identifying and addition of hydrogen """
    pdb_inp = iotbx.pdb.input(file_name=self.file_name)
    pdb_hierarchy = pdb_inp.construct_hierarchy()
    elements = pdb_hierarchy.atoms().extract_element()
    h_count_0 = elements.count(' H') + elements.count(' D')
    new_pdb_str,_ = check_and_add_hydrogen(
      file_name=self.file_name,
      verbose=False)

    pdb_inp = iotbx.pdb.input(source_info=None, lines=new_pdb_str)
    pdb_hierarchy = pdb_inp.construct_hierarchy()
    elements = pdb_hierarchy.atoms().extract_element()
    h_count_1 = elements.count(' H') + elements.count(' D')

    self.assertEqual(h_count_0,0)
    self.assertTrue(h_count_1>0)

  def tearDown(self):
    """ delete files created in during testing"""
    if self.file_to_delete:
      for fn in self.file_to_delete:
        if os.path.isfile(fn): os.remove(fn)

def run_selected_tests():
  """  Run selected tests

  1) List in "tests" the names of the particular test you want to run
  2) Comment out unittest.main()
  3) Un-comment unittest.TextTestRunner().run(run_selected_tests())
  """
  tests = ['test_something']
  suite = unittest.TestSuite(map(MyTestCase, tests))
  return suite


if __name__ == '__main__':
  # use for individual tests
  # unittest.TextTestRunner().run(run_selected_tests())

  # Use to run all tests
  unittest.main(verbosity=0)
