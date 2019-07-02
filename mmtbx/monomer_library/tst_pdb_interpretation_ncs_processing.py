from __future__ import absolute_import, division, print_function
import mmtbx.monomer_library.pdb_interpretation as pdb_interpretation
from libtbx.utils import null_out
from datetime import datetime
import iotbx.phil
import unittest
import shutil
import os
from six.moves import map

__author__ = 'Youval Dar'

# control if temp folder will be deleted after test
DEBUG_MODE = False

test_pdb_str  = '''\
CRYST1  577.812  448.715  468.790  90.00  90.00  90.00 P 1
SCALE1      0.001731  0.000000  0.000000        0.00000
SCALE2      0.000000  0.002229  0.000000        0.00000
SCALE3      0.000000  0.000000  0.002133        0.00000
ATOM      1  CA  LYS A 151      10.766   9.333  12.905  1.00 44.22           C
ATOM      2  CA  LYS A 152      10.117   9.159  11.610  1.00 49.42           C
ATOM      3  CA  LYS A 153       9.099   8.000  11.562  1.00 46.15           C
ATOM      4  CA  LYS A 154       8.000   8.202  11.065  1.00 52.97           C
ATOM      5  CA  LYS A 155      11.146   9.065  10.474  1.00 41.68           C
ATOM      6  CA  LYS A 156      10.547   9.007   9.084  1.00 55.55           C
ATOM      7  CA  LYS A 157      11.545   9.413   8.000  1.00 72.27           C
ATOM      8  CA  LYS A 158      12.277  10.718   8.343  1.00 75.78           C
ATOM      9  CA  LYS A 159      11.349  11.791   8.809  1.00 75.88           C
TER
ATOM    222  CA  LEU X  40      94.618  -5.253  91.582  1.00 87.10           C
ATOM    223  CA  ARG X  41      62.395  51.344  80.786  1.00107.25           C
ATOM    224  CA  ARG X  42      62.395  41.344  80.786  1.00107.25           C
TER
ATOM      1  CA  THR D   1       8.111  11.080  10.645  1.00 20.00           C
ATOM      2  CA  THR D   2       8.000   9.722  10.125  1.00 20.00           C
ATOM      3  CA  THR D   3       8.075   8.694  11.249  1.00 20.00           C
ATOM      4  CA  THR D   4       8.890   8.818  12.163  1.00 20.00           C
ATOM      5  CA  THR D   5       9.101   9.421   9.092  1.00 20.00           C
ATOM      6  CA  THR D   6       9.001  10.343   8.000  1.00 20.00           C
ATOM      7  CA  THR D   7       8.964   8.000   8.565  1.00 20.00           C
TER
ATOM      1  CA  LYS B 151       6.855   8.667  15.730  1.00 44.22           C
ATOM      2  CA  LYS B 152       5.891   8.459  14.655  1.00 49.42           C
ATOM      3  CA  LYS B 153       6.103   7.155  13.858  1.00 46.15           C
ATOM      4  CA  LYS B 154       5.138   6.438  13.633  1.00 52.97           C
ATOM      5  CA  LYS B 155       5.801   9.685  13.736  1.00 41.68           C
ATOM      6  CA  LYS B 156       4.731   9.594  12.667  1.00 55.55           C
ATOM      7  CA  LYS B 157       4.334  10.965  12.119  1.00 72.27           C
ATOM      8  CA  LYS B 158       4.057  11.980  13.238  1.00 75.78           C
ATOM      9  CA  LYS B 159       3.177  11.427  14.310  1.00 75.88           C
TER
ATOM      1  CA  LYS C 151       6.987   4.106  17.432  1.00 44.22           C
ATOM      2  CA  LYS C 152       6.017   3.539  16.502  1.00 49.42           C
ATOM      3  CA  LYS C 153       6.497   3.492  15.036  1.00 46.15           C
ATOM      4  CA  LYS C 154       6.348   2.458  14.400  1.00 52.97           C
ATOM      5  CA  LYS C 155       4.647   4.221  16.634  1.00 41.68           C
ATOM      6  CA  LYS C 156       3.552   3.605  15.788  1.00 55.55           C
ATOM      7  CA  LYS C 157       2.154   3.953  16.298  1.00 72.27           C
ATOM      8  CA  LYS C 158       2.014   3.732  17.811  1.00 75.78           C
ATOM      9  CA  LYS C 159       2.558   2.413  18.250  1.00 75.88           C
TER
ATOM    222  CA  LEU Y  40     194.618   5.253  81.582  1.00 87.10           C
ATOM    223  CA  ARG Y  41     162.395  41.344  70.786  1.00107.25           C
ATOM    224  CA  ARG Y  42     162.395  31.344  70.786  1.00107.25           C
TER
ATOM      1  CA  THR E   1       8.111 -10.645  11.080  1.00 20.00           C
ATOM      2  CA  THR E   2       8.000 -10.125   9.722  1.00 20.00           C
ATOM      3  CA  THR E   3       8.075 -11.249   8.694  1.00 20.00           C
ATOM      4  CA  THR E   4       8.890 -12.163   8.818  1.00 20.00           C
ATOM      5  CA  THR E   5       9.101  -9.092   9.421  1.00 20.00           C
ATOM      6  CA  THR E   6       9.001  -8.000  10.343  1.00 20.00           C
ATOM      7  CA  THR E   7       8.964  -8.565   8.000  1.00 20.00           C
TER
'''

master_params = '''\
    include scope mmtbx.monomer_library.pdb_interpretation.master_params_str
'''


class TestPDBinterpretationNCSSearch(unittest.TestCase):

  def setUp(self):
    """ Create temporary folder for temp files produced during test """
    self.current_dir = os.getcwd()
    now = datetime.now().strftime("%I%p_%m_%d_%Y")
    test_name = self.__class__.__name__
    self.tempdir = '{}_{}'.format(test_name, now)
    if not os.path.isdir(self.tempdir):
      os.mkdir(self.tempdir)
    os.chdir(self.tempdir)
    #
    open('test.pdb','w').write(test_pdb_str)

  def test_calling_pdb_interpretation(self):
    """ Make sure can create NCS object and change search parameters """
    params_phil = iotbx.phil.parse(
      input_string=master_params,
      process_includes=True)
    params = params_phil.extract()
    # Turn on NCS search and adjust one of the parameters
    params.ncs_search.enabled = True
    pdb_processed_file = pdb_interpretation.run(
    args=['test.pdb'],
    params=params,
    log=null_out())
    ncs_obj = pdb_processed_file.ncs_obj
    self.assertEqual(ncs_obj.number_of_ncs_groups,2)

  def tearDown(self):
    """
    remove temp files and folder
    When DEBUG_MODE = True, the temporary folder will not be deleted
    """
    os.chdir(self.current_dir)
    if not DEBUG_MODE:
      try:
        shutil.rmtree(self.tempdir)
      except Exception:
        pass

def run_selected_tests():
  """  Run selected tests

  1) List in "tests" the names of the particular test you want to run
  2) Comment out unittest.main()
  3) Un-comment unittest.TextTestRunner().run(run_selected_tests())
  """
  tests = ['some_test_name']
  suite = unittest.TestSuite(list(map(TestPDBinterpretationNCSSearch, tests)))
  return suite


if __name__ == '__main__':
  # use for individual tests
  # DEBUG_MODE = True
  # unittest.TextTestRunner().run(run_selected_tests())

  # Use to run all tests
  unittest.main(verbosity=0)
