from __future__ import division
from iotbx import pdb
from libtbx.test_utils import approx_equal
from  iotbx.pdb.multimer_reconstruction import multimer
import os

'''
@author: Youval Dar (LBL)
'''


#def exercise_mtrix(regression_pdb):
  #pdb_inp = pdb.input(file_name=op.join(regression_pdb, "pdb1a1q.ent"))
  #mtrix_info = pdb_inp.process_mtrix_records()
  #assert len(mtrix_info) == 2
  #assert approx_equal(mtrix_info[0].values, [[
    #-0.952060, 0.305556,-0.014748,
    #-0.305663,-0.952124, 0.005574,
    #-0.012339, 0.009815, 0.999876], [
      #-22.67001, 73.03197, 0.78307]])
  #assert mtrix_info[0].coordinates_present
  #assert mtrix_info[1].serial_number == "  2"


def exercise_multimer_reconstruction():
  '''
  Test multimer reconstruction

  - create a file, 'multimer_test_data.pdb', using the pdb_test_data
  >>> v = multimer('multimer_test_data.pdb')
  >>> v.assembled_multimer
  '''
  pdb_test_data = '''\
REMARK   0 Test molecule with BIOMOLECULE: 1
REMARK   0
REMARK   0 The test will generate the biomolecule (the multimer assembly)
REMARK   0 from the transformation matrices writen below
REMARK   0 and then compare the results to the calculated expected one
REMARK 350 CRYSTALLOGRAPHIC OPERATIONS ARE GIVEN.
REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.000000
REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.000000
REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.000000
REMARK 350   BIOMT1   2  1.000000  0.000000  0.000000        0.000000
REMARK 350   BIOMT2   2  0.000000  0.000000 -1.000000        0.000000
REMARK 350   BIOMT3   2  0.000000  1.000000  0.000000        0.000000
REMARK 350   BIOMT1   3  0.000000  0.000000  1.000000        0.000000
REMARK 350   BIOMT2   3  0.000000  1.000000  0.000000        0.000000
REMARK 350   BIOMT3   3 -1.000000  0.000000  0.000000        0.000000
REMARK 350   BIOMT1   4  0.000000 -1.000000  0.000000        0.000000
REMARK 350   BIOMT2   4  1.000000  0.000000  0.000000        0.000000
REMARK 350   BIOMT3   4  0.000000  0.000000  1.000000        0.000000
REMARK 350   BIOMT1   5  0.000000  0.000000  1.000000        0.000000
REMARK 350   BIOMT2   5  0.000000  1.000000  0.000000        0.000000
REMARK 350   BIOMT3   5 -1.000000  0.000000  0.000000        0.000000
REMARK 350   BIOMT1   6  0.000000 -1.000000  0.000000        0.000000
REMARK 350   BIOMT2   6  0.000000  0.000000  1.000000        0.000000
REMARK 350   BIOMT3   6 -1.000000  0.000000  0.000000        0.000000
REMARK 350   BIOMT1   7  0.500000 -0.866025  0.000000        0.000000
REMARK 350   BIOMT2   7  0.866025  0.500000  0.000000        0.000000
REMARK 350   BIOMT3   7  0.000000  0.000000  1.000000        0.000000
REMARK 350   BIOMT1   8 -0.500000 -0.866025  0.000000        0.000000
REMARK 350   BIOMT2   8  0.866025 -0.500000  0.000000        0.000000
REMARK 350   BIOMT3   8  0.000000  0.000000  1.000000        0.000000
REMARK 350   BIOMT1   9  1.000000  0.000000  0.000000        0.000000
REMARK 350   BIOMT2   9  0.000000  1.000000  0.000000        0.500000
REMARK 350   BIOMT3   9  0.000000  0.000000  1.000000        0.000000
REMARK 350   BIOMT1  10 -0.500000 -0.866025  0.000000        0.000000
REMARK 350   BIOMT2  10  0.866025 -0.500000  0.000000        0.000000
REMARK 350   BIOMT3  10  0.000000  0.000000  1.000000       -1.000000
MTRIX1   1  1.000000  0.000000  0.000000        0.00000    1
MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1
MTRIX3   1  0.000000  0.000000  1.000000        0.00000    1
MTRIX1   2  1.000000  0.000000  0.000000        0.00000
MTRIX2   2  0.000000  0.000000 -1.000000        0.00000
MTRIX3   2  0.000000  1.000000  0.000000        0.00000
MTRIX1   3  0.500000 -0.866025  0.000000        0.00000
MTRIX2   3  0.866025  0.500000  0.000000        0.00000
MTRIX3   3  0.000000  0.000000  1.000000        0.00000
MTRIX1   4 -0.500000 -0.866025  0.000000        0.00000
MTRIX2   4  0.866025 -0.500000  0.000000        0.00000
MTRIX3   4  0.000000  0.000000  1.000000        0.00000
MTRIX1   5  1.000000  0.000000  0.000000        0.00000
MTRIX2   5  0.000000  1.000000  0.000000        0.50000
MTRIX3   5  0.000000  0.000000  1.000000        0.00000
MTRIX1   6 -0.500000 -0.866025  0.000000        0.00000
MTRIX2   6  0.866025 -0.500000  0.000000        0.00000
MTRIX3   6  0.000000  0.000000  1.000000       -1.00000
ATOM      1   N  ILE A  40       1.000   1.000   1.000  1.00162.33           C
ATOM      2  CA  LEU A  40      94.618  -5.253  91.582  1.00 87.10           C
ATOM      3   C  ARG B  40      62.395  51.344  80.786  1.00107.25           C
HETATM    4  C1  EDO A  40      39.954  51.526  72.372  0.33 60.93           C
'''
  ba_expected_results  = [[1.0, 1.0, 1.0], [1.0, -1.0, 1.0], [1.0, 1.0, -1.0], [-1.0, 1.0, 1.0], [1.0, 1.0, -1.0],
                     [-1.0, 1.0, -1.0], [-0.366025, 1.366025, 1.0], [-1.366025, 0.366025, 1.0], [1.0, 1.5, 1.0],
                     [-1.366025, 0.366025, 0.0], [94.618, -5.253, 91.582], [94.618, -91.582, -5.253],
                     [91.582, -5.253, -94.618], [5.253, 94.618, 91.582], [91.582, -5.253, -94.618],
                     [5.253, 91.582, -94.618], [51.858229, 79.315053, 91.582], [-42.759771, 84.568053, 91.582],
                     [94.618, -4.753, 91.582], [-42.759771, 84.568053, 90.582], [62.395, 51.344, 80.786],
                     [62.395, -80.786, 51.344], [80.786, 51.344, -62.395], [-51.344, 62.395, 80.786],
                     [80.786, 51.344, -62.395], [-51.344, 80.786, -62.395], [-13.267688, 79.70763, 80.786],
                     [-75.662688, 28.36363, 80.786], [62.395, 51.844, 80.786], [-75.662688, 28.36363, 79.786],
                     [39.954, 51.526, 72.372], [39.954, -72.372, 51.526], [72.372, 51.526, -39.954],
                     [-51.526, 39.954, 72.372], [72.372, 51.526, -39.954], [-51.526, 72.372, -39.954],
                     [-24.645804, 60.364163, 72.372], [-64.599804, 8.838163, 72.372], [39.954, 52.026, 72.372],
                     [-64.599804, 8.838163, 71.372]]
  cau_expected_results  = [[1.0, 1.0, 1.0], [1.0, -1.0, 1.0], [-0.366025, 1.366025, 1.0], [-1.366025, 0.366025, 1.0],
                           [1.0, 1.5, 1.0], [-1.366025, 0.366025, 0.0], [94.618, -5.253, 91.582],
                           [94.618, -91.582, -5.253], [51.858229, 79.315053, 91.582], [-42.759771, 84.568053, 91.582],
                           [94.618, -4.753, 91.582], [-42.759771, 84.568053, 90.582], [62.395, 51.344, 80.786],
                           [62.395, -80.786, 51.344], [-13.267688, 79.70763, 80.786], [-75.662688, 28.36363, 80.786],
                           [62.395, 51.844, 80.786], [-75.662688, 28.36363, 79.786], [39.954, 51.526, 72.372],
                           [39.954, -72.372, 51.526], [-24.645804, 60.364163, 72.372], [-64.599804, 8.838163, 72.372],
                           [39.954, 52.026, 72.372], [-64.599804, 8.838163, 71.372]]

  # Write the test data into a file
  #print 'current working directory is: {}'.format(os.getcwd())
  f = open('multimer_test_data.pdb', 'w')
  f.write(pdb_test_data)
  f.close()

  # use BIOMT data
  ba_multimer_data = multimer('multimer_test_data.pdb','ba')
  ba_multimer_xyz = ba_multimer_data.get_xyz()
  # use MTRIX data
  cau_multimer_data = multimer('multimer_test_data.pdb','cau')
  cau_multimer_xyz = cau_multimer_data.get_xyz()

  # The multimer processing is my chain, so we need to sort lists before we compare them
  ba_multimer_xyz.sort()
  ba_expected_results.sort()
  cau_multimer_xyz.sort()
  cau_expected_results.sort()

  assert approx_equal(ba_expected_results,ba_multimer_xyz,eps=0.001)
  assert approx_equal(cau_expected_results,cau_multimer_xyz,eps=0.001)

  # remove the file created for the test
  os.remove('multimer_test_data.pdb')


def run():
  exercise_multimer_reconstruction()


if (__name__ == "__main__"):
  run()
