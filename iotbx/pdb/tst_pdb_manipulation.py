from __future__ import division
from  iotbx.pdb.multimer_reconstruction import format_num_as_str
from  iotbx.pdb.multimer_reconstruction import ncs_group_object
from  iotbx.pdb.multimer_reconstruction import multimer
from libtbx.test_utils import approx_equal
from libtbx.phil import parse
from iotbx import pdb
import string
import unittest
import tempfile
import shutil
import os


class TestMultimerReconstruction(unittest.TestCase):
  '''Test multimer reconstruction from BIOMT and MTRIX records'''


  def setUp(self):
    '''
    Create temporary file: multimer_test_data.pdb
    '''
    self.currnet_dir = os.getcwd()
    self.tempdir = tempfile.mkdtemp('tempdir')
    os.chdir(self.tempdir)
    # Write the test data into a file
    open('multimer_test_data.pdb', 'w').write(pdb_test_data)
    open('multimer_test_data2.pdb', 'w').write(pdb_test_data2)

    self.user_phil = parse("""
      ncs_refinement {
        dont_apply_when_coordinates_present = False
        apply_to_all_chains = False
        ncs_group {
          transform {
            rotation = (1.0,1.0,1.0,0.2,0.5,0.6,0.7,0.8,0.9)
            translation = (1,2,3)
            }
          transform {
            rotation = 0.1,1.0,1.0,0.2,0.5,0.6,0.7,0.8,0.9
            translation = (0,2,1)
            }
          transform {
            rotation = 1.0,0.2,1.0,0.2,0.5,0.6,0.7,0.8,0.9
            translation = (-1,3,-2)
            coordinates_present = True
            }
          apply_to_selection = chain A
        }
        ncs_group {
          transform {
            rotation = 0.2,1.0,1.0,0.2,0.5,0.6,0.7,0.8,0.9
            translation = (0,0,0)
          }
          apply_to_selection = chain A or chain C
          apply_to_selection = and name ca
        }
      }
      """)

  def test_MTRIX(self):
    '''Test MTRIX record processing'''

    cau_expected_results  = [
    [1.0, 1.0, 1.0], [1.0, -1.0, 1.0], [-0.366025, 1.366025, 1.0], [-1.366025, 0.366025, 1.0],
    [1.0, 1.5, 1.0], [94.618, -5.253, 91.582],
    [94.618, -91.582, -5.253], [51.858229, 79.315053, 91.582], [-42.759771, 84.568053, 91.582],
    [94.618, -4.753, 91.582], [62.395, 51.344, 80.786],
    [62.395, -80.786, 51.344], [-13.267688, 79.70763, 80.786], [-75.662688, 28.36363, 80.786],
    [62.395, 51.844, 80.786], [39.954, 51.526, 72.372],
    [39.954, -72.372, 51.526], [-24.645804, 60.364163, 72.372], [-64.599804, 8.838163, 72.372],
    [39.954, 52.026, 72.372]]

    # use MTRIX data
    cau_multimer_data = multimer(
      file_name='multimer_test_data.pdb',
      reconstruction_type='cau')
    cau_multimer_xyz = list(cau_multimer_data.sites_cart())

    cau_multimer_xyz.sort()
    cau_expected_results.sort()
    assert approx_equal(cau_expected_results,cau_multimer_xyz,eps=0.001)

    # Test that the number of MTRIX record to process is correct
    self.assertEqual(cau_multimer_data.number_of_transforms,4)

    # Test getting non-rounded ASU
    transforms_obj = cau_multimer_data.transforms_obj
    source_xyz = cau_multimer_data.get_source_hierarchy().atoms().extract_xyz()
    xyz = transforms_obj.apply_transforms(
      ncs_coordinates = source_xyz,
      round_coordinates=False)
    cau_multimer_xyz = list(xyz)
    cau_multimer_xyz.sort()
    assert approx_equal(cau_expected_results,cau_multimer_xyz,eps=0.00001)

    # Test multimer without rounding
    cau_multimer_data = multimer(
      file_name='multimer_test_data.pdb',
      round_coordinates=False,
      reconstruction_type='cau')
    cau_multimer_xyz = list(cau_multimer_data.sites_cart())
    cau_multimer_xyz.sort()
    cau_expected_results.sort()
    assert approx_equal(cau_expected_results,cau_multimer_xyz,eps=0.00001)




  def test_BIOMT(self):
    '''Test MTRIX record processing'''

    ba_expected_results  = [
    [1.0, 1.0, 1.0], [1.0, -1.0, 1.0], [1.0, 1.0, -1.0], [-1.0, 1.0, 1.0], [1.0, 1.0, -1.0],
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

    # use BIOMT data
    ba_multimer_data = multimer(
      file_name='multimer_test_data.pdb',
      reconstruction_type='ba')
    ba_multimer_xyz = list(ba_multimer_data.sites_cart())

    # The multimer processing is my chain, so we need to sort lists before we compare them
    ba_multimer_xyz.sort()
    ba_expected_results.sort()

    assert approx_equal(ba_expected_results,ba_multimer_xyz,eps=0.001)
    self.assertEqual(ba_multimer_data.number_of_transforms,9)

  def test_transform_count(self):
    '''
    Test correct MTRIX transform counting and
    new chains naming'''
    # use MTRIX data
    m = multimer(
      file_name='multimer_test_data2.pdb',reconstruction_type='cau')
    self.assertEqual(m.number_of_transforms, 1)

  def test_trasformation_application_order(self):
    """
    Test that we build the new assembly by applying each transformation
    to all chains by iterating over the chains first. To be constant with
    NCS/ASU transformations that are done during refinement.

    Note that if chains are broken, their copy will be concatenated
    """
    m = multimer(pdb_str=pdb_test_data3,reconstruction_type='cau')
    self.assertEquals(m.ncs_unique_chains_ids, ('A', 'B'))
    self.assertEquals(pdb_test_data3_expected_results,
                      m.assembled_multimer.as_pdb_string())


  def test_selection_naming(self):
    """
    Verify naming of selection done as expected
    """
    transforms_obj = ncs_group_object()
    result = []
    for i in range(26*2):
      result.append(transforms_obj.produce_selection_name())
    chr_list = string.ascii_uppercase
    expected = ['S'+a+b for a in chr_list[:2] for b in chr_list]
    self.assertEqual(result,expected)

  def test_ncs_group_object(self):
      """ Verify that phil parameters are properly read   """
      transforms_obj = ncs_group_object()
      transforms_obj.populate_ncs_group_object(
        ncs_refinement_params = self.user_phil)
      gr1 = transforms_obj.ncs_refinement_groups[0]
      gr2 = transforms_obj.ncs_refinement_groups[1]
      result = gr1.transform[1].rotation
      expected = [0.1, 1.0, 1.0, 0.2, 0.5, 0.6, 0.7, 0.8, 0.9]
      self.assertEqual(result,expected)
      result = gr1.transform[1].translation
      expected = [0.0, 2.0, 1.0]
      self.assertEqual(result,expected)
      result = gr1.transform[2].coordinates_present
      expected = True
      self.assertEqual(result,expected)
      result = gr2.apply_to_selection
      expected = ['chain A or chain C', 'and name ca']
      self.assertEqual(result,expected)
      result = transforms_obj.dont_apply_when_coordinates_present
      self.assertEqual(result,False)
      result = transforms_obj.apply_to_all_chains
      self.assertEqual(result,False)

  def test_ncs_copies_naming(self):
    transforms_obj = ncs_group_object()
    result =  transforms_obj.make_chains_names((1,2),('A','B'))
    expected = {'A_001': 'C', 'A_002': 'D', 'B_001': 'E', 'B_002': 'F'}
    self.assertEqual(result,expected)

  def test_coordinate_mapping(self):
    """
    Verify that atoms mapping to chain-transformation ID
    """
    transforms_obj = ncs_group_object()
    pdb_inp = pdb.input(source_info=None, lines=pdb_test_data2)
    pdb_obj = pdb.hierarchy.input(pdb_string=pdb_test_data2)
    transform_info = pdb_inp.process_mtrix_records()
    transforms_obj.populate_ncs_group_object(
      transform_info = transform_info,
      pdb_hierarchy_inp = pdb_obj)

    result = {k:list(v) for (k,v) in transforms_obj.asu_to_ncs_map.iteritems()}
    expected = {'A': [0, 1, 2, 5, 6], 'B': [3, 4]}
    self.assertEqual(result,expected)
    result =  transforms_obj.ncs_to_asu_map
    expected = {'A_003': [7, 12], 'B_003': [12, 14]}
    self.assertEqual(result,expected)

    result = transforms_obj.ncs_copies_chains_names
    expected = {'B_003': 'D', 'A_003': 'C'}
    self.assertEqual(result,expected)

  def test_num_to_str(self):
    self.assertEqual(format_num_as_str(946),'946')
    self.assertEqual(format_num_as_str(46),'046')
    self.assertEqual(format_num_as_str(6),'006')
    self.assertEqual(format_num_as_str(0),'000')

  def tearDown(self):
    '''remove temp files and folder'''
    os.chdir(self.currnet_dir)
    shutil.rmtree(self.tempdir)


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
MTRIX1   1  1.000000  0.000000  0.000000        0.00000
MTRIX2   1  0.000000  1.000000  0.000000        0.00000
MTRIX3   1  0.000000  0.000000  1.000000        0.00000
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
MTRIX1   6 -0.500000 -0.866025  0.000000        0.00000    1
MTRIX2   6  0.866025 -0.500000  0.000000        0.00000    1
MTRIX3   6  0.000000  0.000000  1.000000       -1.00000    1
ATOM      1   N  ILE A  40       1.000   1.000   1.000  1.00162.33           C
ATOM      2  CA  LEU A  40      94.618  -5.253  91.582  1.00 87.10           C
ATOM      3   C  ARG B  40      62.395  51.344  80.786  1.00107.25           C
HETATM    4  C1  EDO A  40      39.954  51.526  72.372  0.33 60.93           C
'''

pdb_test_data2="""\
MTRIX1   1  1.000000  0.000000  0.000000        0.00000    1
MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1
MTRIX3   1  0.000000  0.000000  1.000000        0.00000    1
MTRIX1   2  0.496590 -0.643597  0.582393        0.00000    1
MTRIX2   2  0.867925  0.376088 -0.324443        0.00000    1
MTRIX3   2 -0.010221  0.666588  0.745356        0.00000    1
MTRIX1   3 -0.317946 -0.173437  0.932111        0.00000
MTRIX2   3  0.760735 -0.633422  0.141629        0.00000
MTRIX3   3  0.565855  0.754120  0.333333        0.00000
ATOM      1  N   THR A   1       9.670  10.289  11.135  1.00 20.00           N
ATOM      2  CA  THR A   1       9.559   8.931  10.615  1.00 20.00           C
ATOM      3  C   THR A   1       9.634   7.903  11.739  1.00 20.00           C
ATOM      4  O   THR B   1      10.449   8.027  12.653  1.00 20.00           O
ATOM      5  CB  THR B   1      10.660   8.630   9.582  1.00 20.00           C
ATOM      6  OG1 THR A   1      10.560   9.552   8.490  1.00 20.00           O
ATOM      7  CG2 THR A   1      10.523   7.209   9.055  1.00 20.00           C
TER
"""

pdb_test_data3="""\
MTRIX1   1  1.000000  0.000000  0.000000        0.00000    1
MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1
MTRIX3   1  0.000000  0.000000  1.000000        0.00000    1
MTRIX1   2  0.309017 -0.809017  0.500000        0.00000
MTRIX2   2  0.809017  0.500000  0.309017        0.00000
MTRIX3   2 -0.500000  0.309017  0.809017        0.00000
MTRIX1   3 -0.809017 -0.500000  0.309017        0.00000
MTRIX2   3  0.500000 -0.309017  0.809017        0.00000
MTRIX3   3 -0.309017  0.809017  0.500000        0.00000
ATOM    749  O   UNK A  90      28.392  67.262  97.682  1.00  0.00           O
ATOM    750  N   UNK A  91      30.420  66.924  98.358  1.00  0.00           N
TER
ATOM   1495  N   UNK B  67      33.124   2.704 114.920  1.00  0.00           N
END
"""

pdb_test_data3_expected_results = """\
ATOM    749  O   UNK A  90      28.392  67.262  97.682  1.00  0.00           O
ATOM    750  N   UNK A  91      30.420  66.924  98.358  1.00  0.00           N
TER
ATOM   1495  N   UNK B  67      33.124   2.704 114.920  1.00  0.00           N
TER
ATOM    749  O   UNK C  90       3.199  86.786  85.616  1.00  0.00           O
ATOM    750  N   UNK C  91       4.437  88.467  85.044  1.00  0.00           N
TER
ATOM   1495  N   UNK E  67      65.508  63.662  77.246  1.00  0.00           N
TER
ATOM    749  O   UNK D  90     -26.415  72.437  94.483  1.00  0.00           O
ATOM    750  N   UNK D  91     -27.678  74.103  93.921  1.00  0.00           N
TER
ATOM   1495  N   UNK F  67       7.362 108.699  49.412  1.00  0.00           N
TER
"""

if __name__ == "__main__":
  unittest.main()
