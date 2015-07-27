from __future__ import division
from  iotbx.pdb.multimer_reconstruction import multimer
from iotbx.ncs.ncs_preprocess import format_num_as_str
from iotbx.ncs.ncs_preprocess import ncs_group_object
from mmtbx.ncs.ncs_search import is_same_transform
from mmtbx.ncs.ncs_utils import apply_transforms
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
from scitbx import matrix
import iotbx.ncs as ncs
from iotbx import pdb
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

  # @unittest.SkipTest
  def test_MTRIX(self):
    '''Test MTRIX record processing'''
    # print sys._getframe().f_code.co_name
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
    multimer_data = multimer(
      file_name='multimer_test_data.pdb',
      reconstruction_type='cau')
    cau_multimer_xyz = list(multimer_data.sites_cart())

    cau_multimer_xyz.sort()
    cau_expected_results.sort()
    assert approx_equal(cau_expected_results,cau_multimer_xyz,eps=0.001)

    # Test that the number of MTRIX record to process is correct
    self.assertEqual(multimer_data.number_of_transforms,4)

    # Test getting non-rounded ASU
    source_xyz = multimer_data.get_ncs_hierarchy().atoms().extract_xyz()
    xyz = apply_transforms(
      ncs_coordinates = source_xyz,
      ncs_restraints_group_list = multimer_data.get_ncs_restraints_group_list(),
      total_asu_length = multimer_data.total_asu_length(),
      extended_ncs_selection = flex.size_t_range(source_xyz.size()),
      round_coordinates=False)
    cau_multimer_xyz = list(xyz)
    cau_multimer_xyz.sort()
    assert approx_equal(cau_expected_results,cau_multimer_xyz,eps=0.00001)

    # Test multimer without rounding
    multimer_data = multimer(
      file_name='multimer_test_data.pdb',
      round_coordinates=False,
      reconstruction_type='cau')
    cau_multimer_xyz = list(multimer_data.sites_cart())
    cau_multimer_xyz.sort()
    cau_expected_results.sort()
    assert approx_equal(cau_expected_results,cau_multimer_xyz,eps=0.00001)

  # @unittest.SkipTest
  def test_BIOMT(self):
    '''Test MTRIX record processing'''
    # print sys._getframe().f_code.co_name
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
    assert ba_multimer_data.new_annotation is None

    # The multimer processing is my chain, so we need to sort lists before we compare them
    ba_multimer_xyz.sort()
    ba_expected_results.sort()

    assert approx_equal(ba_expected_results,ba_multimer_xyz,eps=0.001)
    self.assertEqual(ba_multimer_data.number_of_transforms,9)

  # @unittest.SkipTest
  def test_ncs_copies_naming(self):
    # print sys._getframe().f_code.co_name
    transforms_obj = ncs_group_object()
    result =  transforms_obj.make_chains_names(
      ['chain A_001','chain B_001','chain A_002','chain B_002'],('A','B'))
    expected = {'chain A_001': 'C', 'chain B_001': 'D', 'chain A_002': 'E',
                'chain B_002': 'F'}
    self.assertEqual(result,expected)

  def test_identity_tranform_insertion(self):
    """
    Verify that insertion and reordering of the identity transform is done
    properly
    """
    # print sys._getframe().f_code.co_name
    for pdb_str in [pdb_test_data5,pdb_test_data6]:
      ncs_inp = ncs.input(pdb_string=pdb_test_data5)
      transform_info = ncs_inp.build_MTRIX_object()
      self.assertEqual(len(transform_info.r),3)
      self.assertEqual(len(transform_info.t),3)
      self.assertEqual(transform_info.r[0].is_r3_identity_matrix(),True)
      self.assertEqual(transform_info.t[0].is_col_zero(),True)
      sn = [int(x) for x in transform_info.serial_number]
      self.assertEqual(sn,[1,2,3])


  # @unittest.SkipTest
  def test_adding_transforms_directly(self):
    """
    Verify that processing of transforms provided manually is done properly """
    # print sys._getframe().f_code.co_name
    pdb_obj = pdb.hierarchy.input(pdb_string=pdb_test_data4)
    r = [matrix.sqr([0.1,1.0,1.0,0.2,0.5,0.6,0.7,0.8,0.9])]
    r.append(matrix.sqr([1.0,0.2,1.0,0.2,0.5,0.6,0.7,0.8,0.4]))
    t = [matrix.col([0,2,1])]
    t.append(matrix.col([-1,3,-2]))
    transforms_obj = ncs.input(
      pdb_hierarchy_inp = pdb_obj,
      rotations=r,
      translations=t)

    result = transforms_obj.transform_to_ncs
    expected = {'002': ['chain A_002', 'chain B_002'],
                '003': ['chain A_003', 'chain B_003']}
    self.assertEqual(result,expected)
    result = transforms_obj.ncs_selection_str
    expected = 'chain A or chain B'
    self.assertEqual(result,expected)
    # check that if transforms are provided MTRIX record ignored
    pdb_obj = pdb.hierarchy.input(pdb_string=pdb_test_data2)
    transforms_obj = ncs.input(
      pdb_hierarchy_inp = pdb_obj,
      rotations=r,
      translations=t)
    result = transforms_obj.transform_to_ncs
    expected = {'002': ['chain A_002', 'chain B_002'],
                '003': ['chain A_003', 'chain B_003']}
    self.assertEqual(result,expected)
    result = transforms_obj.ncs_selection_str
    expected = 'chain A or chain B'
    self.assertEqual(result,expected)
    # transforms that are not present
    result = transforms_obj.transform_to_ncs.keys()
    expected = ['003', '002']
    self.assertEqual(result,expected)
    # all transforms
    result = transforms_obj.ncs_transform.keys()
    expected = ['001', '002', '003']
    result.sort()
    self.assertEqual(result,expected)

  # @unittest.SkipTest
  def test_transform_application_order(self):
    """
    Verify that transform order is kept even when chain selection is complex
    """
    # print sys._getframe().f_code.co_name
    pdb_inp = pdb.input(source_info=None, lines=pdb_test_data2)
    pdb_obj = pdb.hierarchy.input(pdb_string=pdb_test_data2)
    transform_info = pdb_inp.process_mtrix_records()
    transforms_obj = ncs.input(
      transform_info=transform_info,
      pdb_hierarchy_inp=pdb_obj)

    expected = ['chain A_002', 'chain B_002', 'chain A_003', 'chain B_003']
    self.assertEqual(transforms_obj.transform_chain_assignment,expected)

    expected = {
      'chain A_002': 'C','chain A_003': 'E','chain A_001': 'A',
      'chain B_002': 'D','chain B_003': 'F','chain B_001': 'B'}
    self.assertEqual(transforms_obj.ncs_copies_chains_names,expected)

    expected = [0, 1, 2, 5, 6]
    results = list(transforms_obj.asu_to_ncs_map['chain A'])
    self.assertEqual(results,expected)

    expected = [3, 4]
    results = list(transforms_obj.asu_to_ncs_map['chain B'])
    self.assertEqual(results,expected)

    expected = [7, 8, 9, 12, 13]
    results = list(transforms_obj.ncs_to_asu_map['chain A_002'])
    self.assertEqual(results,expected)

    expected = [17, 18]
    results = list(transforms_obj.ncs_to_asu_map['chain B_003'])
    self.assertEqual(results,expected)

    self.assertEqual(len(transforms_obj.ncs_atom_selection),21)
    self.assertEqual(transforms_obj.ncs_atom_selection.count(True),7)

  # @unittest.SkipTest
  def test_num_to_str(self):
    # print sys._getframe().f_code.co_name
    self.assertEqual(format_num_as_str(946),'946')
    self.assertEqual(format_num_as_str(46),'046')
    self.assertEqual(format_num_as_str(6),'006')
    self.assertEqual(format_num_as_str(0),'000')

  def test_pdb_writing(self):
    # print sys._getframe().f_code.co_name
    """
    Verify that there are no errors processing the write command
    No inception of the output is done.
    To view the output change the write=False to ,write=True
    """
    transforms_obj = ncs_group_object()
    pdb_inp = pdb.input(source_info=None, lines=pdb_test_data2)
    pdb_obj = pdb.hierarchy.input(pdb_string=pdb_test_data2)
    transform_info = pdb_inp.process_mtrix_records()
    transforms_obj.preprocess_ncs_obj(
      transform_info=transform_info,
      pdb_hierarchy_inp=pdb_obj)

    multimer_data = multimer(
      pdb_str=pdb_test_data2,
      reconstruction_type='cau')

    pdb_hierarchy_asu = multimer_data.assembled_multimer

    # print '--- using ASU hierarchy ---'
    pdbstr = transforms_obj.get_transform_records(
      ncs_only=True,
      pdb_hierarchy=pdb_hierarchy_asu,
      write=False)
    # print pdbstr
    # print '='*50

    pdbstr = transforms_obj.get_transform_records(
      ncs_only=False,
      pdb_hierarchy=pdb_hierarchy_asu,
      write=False)
    # print pdbstr

    # print '--- using the hierarchy of only the master NCS ---'
    pdbstr = transforms_obj.get_transform_records(
      pdb_hierarchy=pdb_obj.hierarchy,
      biomt=True,
      write=False)
    # print pdbstr

    # print '--- from xray structure ---'
    xrs = pdb_hierarchy_asu.extract_xray_structure()
    pdbstr = transforms_obj.get_transform_records(
      # xrs=pdb_obj.xray_structure_simple(),
      xrs=xrs,
      biomt=True,
      write=False)
    # print pdbstr
    # print '='*50

  def test_spec_file_format(self):
    """ Verify that spec object are produced properly """
    # print sys._getframe().f_code.co_name

    multimer_data = multimer(
      pdb_str=pdb_test_data2,
      reconstruction_type='cau')

    trans_obj = ncs.input(pdb_string=pdb_test_data2)

    pdb_hierarchy_asu = multimer_data.assembled_multimer
    spec_output = trans_obj.get_ncs_info_as_spec(
      pdb_hierarchy_asu=pdb_hierarchy_asu,write=False)

    trans_obj2 = ncs.input(spec_ncs_groups=spec_output)

    t1 = trans_obj.ncs_transform['002'].r
    t2 = trans_obj2.ncs_transform['002'].r
    self.assertEqual(t1,t2)
    self.assertEqual(len(trans_obj.ncs_transform),len(trans_obj2.ncs_transform))

    t1 = trans_obj.ncs_to_asu_selection
    t1_expected = {'chain A or chain B':
                     ['chain E or chain F', 'chain C or chain D']}
    self.assertEqual(t1,t1_expected)
    t2 = trans_obj2.ncs_to_asu_selection
    t2_expected = {
      'chain A and (resseq 1:3 or resseq 6:7)':
        ['chain C and (resseq 1:3 or resseq 6:7)',
         'chain E and (resseq 1:3 or resseq 6:7)'],
      'chain B and (resseq 4:5)':
        ['chain D and (resseq 4:5)', 'chain F and (resseq 4:5)']}
    self.assertEqual(t2,t2_expected)

    t1 = trans_obj.tr_id_to_selection['chain A_003']
    t1_expected = ('chain A',
                   'chain E')
    self.assertEqual(t1,t1_expected)

    t2 = trans_obj2.tr_id_to_selection['chain A_003']
    t2_expected = ('chain A and (resseq 1:3 or resseq 6:7)',
                   'chain E and (resseq 1:3 or resseq 6:7)')
    self.assertEqual(t2,t2_expected)

  # @unittest.SkipTest
  def test_writing_spec_file(self):
    # print sys._getframe().f_code.co_name
    """
    Verify that there are no errors processing the write command
    No inception of the output is done. Just making sure it does not break
    To view the output change the write=False to ,write=True
    """
    transforms_obj = ncs_group_object()
    pdb_inp = pdb.input(source_info=None, lines=pdb_test_data2)
    pdb_obj = pdb.hierarchy.input(pdb_string=pdb_test_data2)
    transform_info = pdb_inp.process_mtrix_records()
    transforms_obj.preprocess_ncs_obj(
      # transform_info=transform_info,
      # pdb_hierarchy_inp=pdb_obj,
      pdb_inp=pdb_inp)

    asu = multimer(
      pdb_str=pdb_test_data2,
      reconstruction_type='cau')
    transforms_obj.get_ncs_info_as_spec(
      pdb_hierarchy_asu=asu.assembled_multimer,write=False)

  # @unittest.SkipTest
  def test_proper_biomat_application(self):
    """ Test that when building bio-molecule and then finding NCS relatin
    from it, we get the same rotation and translation"""
    pdb_strings = [pdb_test_data7,pdb_test_data8]
    methods = ['ba','cau']
    for method,pdb_string in zip(methods,pdb_strings):
      pdb_inp = pdb.input(source_info=None, lines=pdb_string)
      crystal_symmetry = pdb_inp.crystal_symmetry()
      m = multimer(
        pdb_str=pdb_string,
        round_coordinates=False,
        reconstruction_type=method,
        error_handle=True,eps=1e-2)
      # The exact transforms from pdb_string
      r1_expected = matrix.sqr(
        [0.309017, -0.951057, 0.0,0.951057, 0.309017,-0.0,0.0,0.0,1.0])
      r2_expected = matrix.sqr(
        [-0.809017,-0.587785,0.0,0.587785,-0.809017,-0.0,0.0,0.0,1.0])
      t1_expected = matrix.col([0,0,7])
      t2_expected = matrix.col([0,0,0])
      # Look at biomt records retrieved from PDB file
      if method == 'ba':
        rec = pdb_inp.process_BIOMT_records()
      else:
        rec = pdb_inp.process_mtrix_records()
      r1 = rec.r[1]
      r2 = rec.r[2]
      t1 = rec.t[1]
      t2 = rec.t[2]
      (the_same, transpose) = is_same_transform(r1,t1,r1_expected,t1_expected)
      self.assertTrue(the_same)
      (the_same, transpose)= is_same_transform(r2,t2,r2_expected,t2_expected)
      self.assertTrue(the_same)
      # Look at the rotation and translation found by the NCS search
      s = m.assembled_multimer.as_pdb_string(crystal_symmetry=crystal_symmetry)
      ncs_obj = ncs.input(pdb_string=s)
      r1 = ncs_obj.ncs_transform['002'].r
      t1 = ncs_obj.ncs_transform['002'].t
      r2 = ncs_obj.ncs_transform['003'].r
      t2 = ncs_obj.ncs_transform['003'].t
      (the_same, transpose) = is_same_transform(r1,t1,r1_expected,t1_expected)
      self.assertTrue(the_same)
      (the_same, transpose)= is_same_transform(r2,t2,r2_expected,t2_expected)
      self.assertTrue(the_same)
      self.assertEqual(ncs_obj.number_of_ncs_groups,1)

  # @unittest.SkipTest
  def test_ignoring_mtrix_rec(self):
    """
    Test ignoring MTRIX record when copies already present in file
    """
    pdb_inp = pdb.input(source_info=None, lines=test_pdb_9)
    m = multimer(
        pdb_str=test_pdb_9,
        round_coordinates=False,
        reconstruction_type='cau',
        error_handle=True,eps=1e-2)
    n1 = m.assembled_multimer.atoms().size()
    n2 = pdb_inp.atoms().size()
    self.assertEqual(n1,n2)


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
ATOM      1   N  ILE A  40       1.000   1.000   1.000  1.00162.33           C
ATOM      2  CA  LEU A  40      94.618  -5.253  91.582  1.00 87.10           C
ATOM      3   C  ARG B  40      62.395  51.344  80.786  1.00107.25           C
HETATM    4  C1  EDO A  40      39.954  51.526  72.372  0.33 60.93           C
'''

pdb_test_data2="""\
MTRIX1   1  1.000000  0.000000  0.000000        0.00000    1
MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1
MTRIX3   1  0.000000  0.000000  1.000000        0.00000    1
MTRIX1   2  0.496590 -0.643597  0.582393        0.00000
MTRIX2   2  0.867925  0.376088 -0.324443        0.00000
MTRIX3   2 -0.010221  0.666588  0.745356        0.00000
MTRIX1   3 -0.317946 -0.173437  0.932111        0.00000
MTRIX2   3  0.760735 -0.633422  0.141629        0.00000
MTRIX3   3  0.565855  0.754120  0.333333        0.00000
ATOM      1  N   THR A   1       9.670  10.289  11.135  1.00 20.00           N
ATOM      2  CA  THR A   2       9.559   8.931  10.615  1.00 20.00           C
ATOM      3  C   THR A   3       9.634   7.903  11.739  1.00 20.00           C
ATOM      4  O   THR B   4      10.449   8.027  12.653  1.00 20.00           O
ATOM      5  CB  THR B   5      10.660   8.630   9.582  1.00 20.00           C
ATOM      6  OG1 THR A   6      10.560   9.552   8.490  1.00 20.00           O
ATOM      7  CG2 THR A   7      10.523   7.209   9.055  1.00 20.00           C
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

pdb_test_data4="""\
ATOM      1  N   THR A   1       9.670  10.289  11.135  1.00 20.00           N
ATOM      2  CA  THR A   2       9.559   8.931  10.615  1.00 20.00           C
ATOM      3  C   THR A   3       9.634   7.903  11.739  1.00 20.00           C
ATOM      4  O   THR B   4      10.449   8.027  12.653  1.00 20.00           O
ATOM      5  CB  THR B   5      10.660   8.630   9.582  1.00 20.00           C
ATOM      6  OG1 THR A   6      10.560   9.552   8.490  1.00 20.00           O
ATOM      7  CG2 THR A   7      10.523   7.209   9.055  1.00 20.00           C
"""

pdb_test_data5="""\
MTRIX1   1  0.309017 -0.809017  0.500000        0.00000
MTRIX2   1  0.809017  0.500000  0.309017        0.00000
MTRIX3   1 -0.500000  0.309017  0.809017        0.00000
MTRIX1   2 -0.809017 -0.500000  0.309017        0.00000
MTRIX2   2  0.500000 -0.309017  0.809017        0.00000
MTRIX3   2 -0.309017  0.809017  0.500000        0.00000
MTRIX1   3  1.000000  0.000000  0.000000        0.00000    1
MTRIX2   3  0.000000  1.000000  0.000000        0.00000    1
MTRIX3   3  0.000000  0.000000  1.000000        0.00000    1
ATOM    749  O   UNK A  90      28.392  67.262  97.682  1.00  0.00           O
ATOM    750  N   UNK A  91      30.420  66.924  98.358  1.00  0.00           N
TER
ATOM   1495  N   UNK B  67      33.124   2.704 114.920  1.00  0.00           N
END
"""

pdb_test_data6="""\
MTRIX1   1  0.309017 -0.809017  0.500000        0.00000
MTRIX2   1  0.809017  0.500000  0.309017        0.00000
MTRIX3   1 -0.500000  0.309017  0.809017        0.00000
MTRIX1   2 -0.809017 -0.500000  0.309017        0.00000
MTRIX2   2  0.500000 -0.309017  0.809017        0.00000
MTRIX3   2 -0.309017  0.809017  0.500000        0.00000
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
ATOM   1495  N   UNK D  67      65.508  63.662  77.246  1.00  0.00           N
TER
ATOM    749  O   UNK E  90     -26.415  72.437  94.483  1.00  0.00           O
ATOM    750  N   UNK E  91     -27.678  74.103  93.921  1.00  0.00           N
TER
ATOM   1495  N   UNK F  67       7.362 108.699  49.412  1.00  0.00           N
TER
"""

pdb_test_data7 = '''\
REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000
REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000
REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000
REMARK 350   BIOMT1   2  0.309017 -0.951057  0.000000        0.00000
REMARK 350   BIOMT2   2  0.951057  0.309017 -0.000000        0.00000
REMARK 350   BIOMT3   2  0.000000  0.000000  1.000000        7.00000
REMARK 350   BIOMT1   3 -0.809017 -0.587785  0.000000        0.00000
REMARK 350   BIOMT2   3  0.587785 -0.809017 -0.000000        0.00000
REMARK 350   BIOMT3   3  0.000000  0.000000  1.000000        0.00000
CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1
ATOM      1  N   ALA A   2      64.807-112.186 260.746  1.00160.99           N
ATOM      2  CA  ALA A   2      64.727-111.450 262.002  1.00159.36           C
ATOM      3  C   ALA A   2      63.960-110.148 261.805  1.00154.38           C
ATOM      4  O   ALA A   2      62.935-109.914 262.452  1.00149.47           O
ATOM      5  CB  ALA A   2      66.123-111.175 262.542  1.00156.98           C
ATOM      6  N   SER A   3      64.474-109.323 260.896  1.00135.75           N
ATOM      7  CA  SER A   3      63.887-108.040 260.510  1.00131.97           C
ATOM      8  C   SER A   3      64.863-107.340 259.575  1.00140.51           C
ATOM      9  O   SER A   3      65.864-107.925 259.165  1.00148.46           O
ATOM     10  CB  SER A   3      63.641-107.147 261.726  1.00126.01           C
ATOM     11  OG  SER A   3      64.002-105.804 261.453  1.00119.04           O
END
'''

pdb_test_data8 = '''\
MTRIX1   1  1.000000  0.000000  0.000000        0.00000    1
MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1
MTRIX3   1  0.000000  0.000000  1.000000        0.00000    1
MTRIX1   2  0.309017 -0.951057  0.000000        0.00000
MTRIX2   2  0.951057  0.309017 -0.000000        0.00000
MTRIX3   2  0.000000  0.000000  1.000000        7.00000
MTRIX1   3 -0.809017 -0.587785  0.000000        0.00000
MTRIX2   3  0.587785 -0.809017 -0.000000        0.00000
MTRIX3   3  0.000000  0.000000  1.000000        0.00000
CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1
ATOM    757  N   ASP A 247      16.068 -20.882 -28.984  1.00 35.93           N
ATOM    758  CA  ASP A 247      15.914 -22.265 -28.600  1.00 47.90           C
ATOM    759  C   ASP A 247      17.130 -23.042 -29.116  1.00 42.32           C
ATOM    760  O   ASP A 247      17.461 -22.986 -30.301  1.00 47.25           O
ATOM    761  CB  ASP A 247      14.621 -22.814 -29.198  1.00 47.22           C
ATOM    762  CG  ASP A 247      14.068 -23.974 -28.412  1.00 61.15           C
ATOM    763  OD1 ASP A 247      14.359 -24.061 -27.196  1.00 63.66           O
ATOM    764  OD2 ASP A 247      13.341 -24.798 -29.012  1.00 77.01           O
ATOM    765  N   VAL A 248      17.808 -23.746 -28.218  1.00 44.08           N
ATOM    766  CA  VAL A 248      19.008 -24.503 -28.584  1.00 46.18           C
ATOM    767  C   VAL A 248      18.668 -25.988 -28.583  1.00 53.97           C
ATOM    768  O   VAL A 248      18.049 -26.478 -27.638  1.00 51.48           O
ATOM    769  CB  VAL A 248      20.185 -24.226 -27.608  1.00 47.55           C
ATOM    770  CG1 VAL A 248      21.414 -25.015 -28.012  1.00 41.43           C
ATOM    771  CG2 VAL A 248      20.513 -22.743 -27.567  1.00 41.64           C
ATOM    772  N   VAL A 249      19.057 -26.697 -29.641  1.00 54.29           N
ATOM    773  CA  VAL A 249      18.662 -28.097 -29.810  1.00 60.17           C
ATOM    774  C   VAL A 249      19.859 -29.041 -29.982  1.00 57.98           C
ATOM    775  O   VAL A 249      20.731 -28.827 -30.828  1.00 58.31           O
ATOM    776  CB  VAL A 249      17.671 -28.280 -30.997  1.00 60.85           C
ATOM    777  CG1 VAL A 249      16.500 -27.300 -30.884  1.00 48.00           C
ATOM    778  CG2 VAL A 249      18.386 -28.110 -32.337  1.00 59.99           C
TER
ATOM    780  N   LYS B 151       4.045  -6.858 -32.823  1.00 45.22           N
ATOM    781  CA  LYS B 151       4.686  -6.715 -34.123  1.00 50.40           C
ATOM    782  C   LYS B 151       5.707  -5.554 -34.172  1.00 47.13           C
ATOM    783  O   LYS B 151       6.820  -5.764 -34.625  1.00 52.91           O
ATOM    784  CB  LYS B 151       3.657  -6.646 -35.268  1.00 40.73           C
ATOM    785  CG  LYS B 151       4.264  -6.627 -36.661  1.00 55.98           C
ATOM    786  CD  LYS B 151       3.272  -7.051 -37.745  1.00 72.14           C
ATOM    787  CE  LYS B 151       2.529  -8.338 -37.375  1.00 75.11           C
ATOM    788  NZ  LYS B 151       3.451  -9.400 -36.884  1.00 75.46           N
ATOM    789  N   ARG B 152       5.369  -4.349 -33.709  1.00 42.01           N
ATOM    790  CA  ARG B 152       6.399  -3.290 -33.702  1.00 40.51           C
ATOM    791  C   ARG B 152       6.155  -2.002 -32.909  1.00 34.21           C
ATOM    792  O   ARG B 152       5.015  -1.605 -32.636  1.00 33.77           O
ATOM    793  CB  ARG B 152       6.845  -2.937 -35.130  1.00 40.62           C
ATOM    794  CG  ARG B 152       5.842  -2.126 -35.925  1.00 45.94           C
ATOM    795  CD  ARG B 152       6.341  -1.926 -37.341  1.00 42.75           C
ATOM    796  NE  ARG B 152       7.478  -1.006 -37.404  1.00 45.27           N
ATOM    797  CZ  ARG B 152       8.177  -0.763 -38.509  1.00 49.68           C
ATOM    798  NH1 ARG B 152       7.860  -1.382 -39.644  1.00 47.81           N
ATOM    799  NH2 ARG B 152       9.192   0.096 -38.482  1.00 48.06           N
END
'''

test_pdb_9 = """\
MTRIX1   1  1.000000  0.000000  0.000000        0.00000    1
MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1
MTRIX3   1  0.000000  0.000000  1.000000        0.00000    1
MTRIX1   2  0.496590 -0.643597  0.582393        0.00000    1
MTRIX2   2  0.867925  0.376088 -0.324443        0.00000    1
MTRIX3   2 -0.010221  0.666588  0.745356        0.00000    1
MTRIX1   3 -0.317946 -0.173437  0.932111        0.00000    1
MTRIX2   3  0.760735 -0.633422  0.141629        0.00000    1
MTRIX3   3  0.565855  0.754120  0.333333        0.00000    1
ATOM      1  N   THR A   1       9.670  10.289  11.135  1.00 20.00           N
ATOM      2  CA  THR A   2       9.559   8.931  10.615  1.00 20.00           C
ATOM      3  C   THR A   3       9.634   7.903  11.739  1.00 20.00           C
ATOM      4  O   THR B   4      10.449   8.027  12.653  1.00 20.00           O
ATOM      5  CB  THR B   5      10.660   8.630   9.582  1.00 20.00           C
ATOM      6  OG1 THR A   6      10.560   9.552   8.490  1.00 20.00           O
ATOM      7  CG2 THR A   7      10.523   7.209   9.055  1.00 20.00           C
END
"""

def exercise_ss_multiplication():
  pdb_str = """\
REMARK 350   BIOMT1   1  1.000000 -0.000000  0.000000        0.00000
REMARK 350   BIOMT2   1 -0.000000  1.000000 -0.000000        0.00000
REMARK 350   BIOMT3   1  0.000000 -0.000000  1.000000        0.00000
REMARK 350   BIOMT1   2 -0.809017 -0.500000  0.309017        0.00000
REMARK 350   BIOMT2   2 -0.500000  0.309017 -0.809017        0.00000
REMARK 350   BIOMT3   2  0.309017 -0.809017 -0.500000        0.00000
REMARK 350   BIOMT1   3 -0.000000  1.000000  0.000000        0.00000
REMARK 350   BIOMT2   3 -0.000000 -0.000000 -1.000000        0.00000
REMARK 350   BIOMT3   3 -1.000000  0.000000 -0.000000        0.00000
REMARK 350   BIOMT1   4  0.809017 -0.500000 -0.309017        0.00000
REMARK 350   BIOMT2   4 -0.500000 -0.309017 -0.809017        0.00000
REMARK 350   BIOMT3   4  0.309017  0.809017 -0.500000       -0.00000
REMARK 350   BIOMT1   5  0.500000  0.309017 -0.809017        0.00000
REMARK 350   BIOMT2   5 -0.309017 -0.809017 -0.500000        0.00000
REMARK 350   BIOMT3   5 -0.809017  0.500000 -0.309017        0.00000
REMARK 350   BIOMT1   6 -0.309017 -0.809017 -0.500000        0.00000
REMARK 350   BIOMT2   6 -0.809017  0.500000 -0.309017        0.00000
REMARK 350   BIOMT3   6  0.500000  0.309017 -0.809017       -0.00000
REMARK 350   BIOMT1   7 -0.809017  0.500000 -0.309017        0.00000
REMARK 350   BIOMT2   7  0.500000  0.309017 -0.809017        0.00000
REMARK 350   BIOMT3   7 -0.309017 -0.809017 -0.500000        0.00000
REMARK 350   BIOMT1   8 -0.809017 -0.500000 -0.309017        0.00000
REMARK 350   BIOMT2   8  0.500000 -0.309017 -0.809017        0.00000
REMARK 350   BIOMT3   8  0.309017 -0.809017  0.500000        0.00000
REMARK 350   BIOMT1   9 -0.309017  0.809017 -0.500000        0.00000
REMARK 350   BIOMT2   9 -0.809017 -0.500000 -0.309017        0.00000
REMARK 350   BIOMT3   9 -0.500000  0.309017  0.809017        0.00000
HELIX    1   1 PRO A   28  LEU A   33  1                                   6
HELIX    7   7 ASP B   74  ALA B   79  1                                   6
HELIX   12  12 LYS C  151  TYR C  159  1                                   9
CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1
ATOM     62  N   PRO A  28      33.884  22.417 110.984 14.01 14.01           N
ATOM     63  CA  PRO A  28      34.083  23.080 112.268 12.01 12.01           C
ATOM     64  C   PRO A  28      32.908  24.004 112.604 12.01 12.01           C
ATOM     65  O   PRO A  28      32.254  24.556 111.723 16.00 16.00           O
ATOM     66  CB  PRO A  28      35.381  23.873 112.096 12.01 12.01           C
ATOM     67  CG  PRO A  28      35.369  24.249 110.611 12.01 12.01           C
ATOM     68  CD  PRO A  28      34.692  23.045 109.949 12.01 12.01           C
ATOM     69  N   ASN A  29      32.694  24.235 113.899 14.01 14.01           N
ATOM     70  CA  ASN A  29      31.700  25.147 114.485 12.01 12.01           C
ATOM     71  C   ASN A  29      30.215  24.812 114.276 12.01 12.01           C
ATOM     72  O   ASN A  29      29.406  25.211 115.111 16.00 16.00           O
ATOM     73  CB  ASN A  29      32.013  26.616 114.122 12.01 12.01           C
ATOM     74  CG  ASN A  29      33.381  27.082 114.575 12.01 12.01           C
ATOM     75  OD1 ASN A  29      34.061  26.468 115.379 16.00 16.00           O
ATOM     76  ND2 ASN A  29      33.841  28.204 114.073 14.01 14.01           N
ATOM     77  N   GLU A  30      30.550  21.845 114.346 14.01 14.01           N
ATOM     78  CA  GLU A  30      29.386  21.544 115.199 12.01 12.01           C
ATOM     79  C   GLU A  30      29.802  20.720 116.428 12.01 12.01           C
ATOM     80  O   GLU A  30      29.237  19.675 116.751 16.00 16.00           O
ATOM     81  CB  GLU A  30      28.306  20.825 114.365 12.01 12.01           C
ATOM     82  CG  GLU A  30      27.795  21.680 113.200 12.01 12.01           C
ATOM     83  CD  GLU A  30      26.582  21.020 112.532 12.01 12.01           C
ATOM     84  OE1 GLU A  30      26.801  20.115 111.689 16.00 16.00           O
ATOM     85  OE2 GLU A  30      25.450  21.412 112.878 16.00 16.00           O
ATOM     86  N   LEU A  31      29.828  23.005 116.792 14.01 14.01           N
ATOM     87  CA  LEU A  31      29.299  22.073 115.816 12.01 12.01           C
ATOM     88  C   LEU A  31      28.149  21.279 116.426 12.01 12.01           C
ATOM     89  O   LEU A  31      28.305  20.133 116.853 16.00 16.00           O
ATOM     90  CB  LEU A  31      30.413  21.167 115.268 12.01 12.01           C
ATOM     91  CG  LEU A  31      31.519  21.894 114.483 12.01 12.01           C
ATOM     92  CD1 LEU A  31      32.517  20.866 113.958 12.01 12.01           C
ATOM     93  CD2 LEU A  31      30.995  22.696 113.286 12.01 12.01           C
ATOM     94  N   GLY A  32      26.947  21.837 116.303 14.01 14.01           N
ATOM     95  CA  GLY A  32      25.679  21.115 116.378 12.01 12.01           C
ATOM     96  C   GLY A  32      25.506  20.004 115.323 12.01 12.01           C
ATOM     97  O   GLY A  32      24.382  19.638 114.991 16.00 16.00           O
ATOM     98  N   LEU A  33      26.624  19.470 114.815 14.01 14.01           N
ATOM     99  CA  LEU A  33      26.791  18.378 113.866 12.01 12.01           C
ATOM    100  C   LEU A  33      27.494  17.166 114.509 12.01 12.01           C
ATOM    101  O   LEU A  33      27.549  16.110 113.890 16.00 16.00           O
ATOM    102  CB  LEU A  33      27.595  18.877 112.644 12.01 12.01           C
ATOM    103  CG  LEU A  33      27.192  20.239 112.037 12.01 12.01           C
ATOM    104  CD1 LEU A  33      28.098  20.578 110.853 12.01 12.01           C
ATOM    105  CD2 LEU A  33      25.745  20.281 111.549 12.01 12.01           C
ATOM   2090  N   ASP B  74      21.209  34.482 113.075 14.01 14.01           N
ATOM   2091  CA  ASP B  74      20.506  33.939 114.233 12.01 12.01           C
ATOM   2092  C   ASP B  74      21.375  33.841 115.474 12.01 12.01           C
ATOM   2093  O   ASP B  74      22.182  32.923 115.611 16.00 16.00           O
ATOM   2094  CB  ASP B  74      19.906  32.571 113.907 12.01 12.01           C
ATOM   2095  CG  ASP B  74      18.979  32.068 114.997 12.01 12.01           C
ATOM   2096  OD1 ASP B  74      19.110  32.533 116.149 16.00 16.00           O
ATOM   2097  OD2 ASP B  74      18.121  31.208 114.708 16.00 16.00           O
ATOM   2098  N   MET B  75      21.194  34.811 116.364 14.01 14.01           N
ATOM   2099  CA  MET B  75      21.934  34.884 117.610 12.01 12.01           C
ATOM   2100  C   MET B  75      22.089  33.512 118.239 12.01 12.01           C
ATOM   2101  O   MET B  75      23.123  33.213 118.837 16.00 16.00           O
ATOM   2102  CB  MET B  75      21.234  35.824 118.583 12.01 12.01           C
ATOM   2103  CG  MET B  75      21.525  37.288 118.322 12.01 12.01           C
ATOM   2104  SD  MET B  75      23.242  37.703 118.659 32.07 32.07           S
ATOM   2105  CE  MET B  75      23.058  38.659 120.162 12.01 12.01           C
ATOM   2106  N   ILE B  76      21.067  32.672 118.091 14.01 14.01           N
ATOM   2107  CA  ILE B  76      21.116  31.324 118.645 12.01 12.01           C
ATOM   2108  C   ILE B  76      22.422  30.683 118.212 12.01 12.01           C
ATOM   2109  O   ILE B  76      22.994  29.858 118.922 16.00 16.00           O
ATOM   2110  CB  ILE B  76      19.929  30.457 118.191 12.01 12.01           C
ATOM   2111  CG1 ILE B  76      20.068  30.070 116.720 12.01 12.01           C
ATOM   2112  CG2 ILE B  76      18.618  31.182 118.436 12.01 12.01           C
ATOM   2113  CD1 ILE B  76      18.979  29.144 116.230 12.01 12.01           C
ATOM   2114  N   LYS B  77      22.902  31.098 117.044 14.01 14.01           N
ATOM   2115  CA  LYS B  77      24.160  30.591 116.525 12.01 12.01           C
ATOM   2116  C   LYS B  77      25.209  30.734 117.613 12.01 12.01           C
ATOM   2117  O   LYS B  77      26.219  30.038 117.612 16.00 16.00           O
ATOM   2118  CB  LYS B  77      24.577  31.362 115.281 12.01 12.01           C
ATOM   2119  CG  LYS B  77      23.624  31.194 114.112 12.01 12.01           C
ATOM   2120  CD  LYS B  77      24.208  30.282 113.048 12.01 12.01           C
ATOM   2121  CE  LYS B  77      25.656  30.627 112.741 12.01 12.01           C
ATOM   2122  NZ  LYS B  77      25.986  32.047 113.045 14.01 14.01           N
ATOM   2123  N   ILE B  78      24.954  31.644 118.547 14.01 14.01           N
ATOM   2124  CA  ILE B  78      25.852  31.863 119.664 12.01 12.01           C
ATOM   2125  C   ILE B  78      25.432  30.891 120.752 12.01 12.01           C
ATOM   2126  O   ILE B  78      26.264  30.367 121.489 16.00 16.00           O
ATOM   2127  CB  ILE B  78      25.788  33.314 120.171 12.01 12.01           C
ATOM   2128  CG1 ILE B  78      26.541  34.227 119.215 12.01 12.01           C
ATOM   2129  CG2 ILE B  78      26.405  33.447 121.552 12.01 12.01           C
ATOM   2130  CD1 ILE B  78      27.879  33.670 118.785 12.01 12.01           C
ATOM   2131  N   ALA B  79      24.131  30.644 120.844 14.01 14.01           N
ATOM   2132  CA  ALA B  79      23.615  29.721 121.839 12.01 12.01           C
ATOM   2133  C   ALA B  79      24.401  28.426 121.758 12.01 12.01           C
ATOM   2134  O   ALA B  79      25.260  28.159 122.591 16.00 16.00           O
ATOM   2135  CB  ALA B  79      22.140  29.461 121.605 12.01 12.01           C
ATOM   4569  N   LYS C 151      65.128  35.635 115.156 14.01 14.01           N
ATOM   4570  CA  LYS C 151      65.530  35.389 113.777 12.01 12.01           C
ATOM   4571  C   LYS C 151      64.991  36.467 112.836 12.01 12.01           C
ATOM   4572  O   LYS C 151      65.493  37.600 112.882 16.00 16.00           O
ATOM   4573  CB  LYS C 151      65.144  33.970 113.358 12.01 12.01           C
ATOM   4574  CG  LYS C 151      65.912  32.877 114.088 12.01 12.01           C
ATOM   4575  CD  LYS C 151      67.128  32.399 113.304 12.01 12.01           C
ATOM   4576  CE  LYS C 151      67.638  31.055 113.817 12.01 12.01           C
ATOM   4577  NZ  LYS C 151      68.536  30.367 112.845 14.01 14.01           N
ATOM   4578  N   PRO C 152      63.957  36.119 111.972 14.01 14.01           N
ATOM   4579  CA  PRO C 152      63.508  37.232 111.112 12.01 12.01           C
ATOM   4580  C   PRO C 152      63.109  38.448 111.934 12.01 12.01           C
ATOM   4581  O   PRO C 152      63.440  39.574 111.564 16.00 16.00           O
ATOM   4582  CB  PRO C 152      62.297  36.661 110.371 12.01 12.01           C
ATOM   4583  CG  PRO C 152      62.611  35.226 110.237 12.01 12.01           C
ATOM   4584  CD  PRO C 152      63.111  34.901 111.606 12.01 12.01           C
ATOM   4585  N   ALA C 153      62.412  38.211 113.040 14.01 14.01           N
ATOM   4586  CA  ALA C 153      61.963  39.291 113.912 12.01 12.01           C
ATOM   4587  C   ALA C 153      63.098  40.266 114.150 12.01 12.01           C
ATOM   4588  O   ALA C 153      63.069  41.394 113.658 16.00 16.00           O
ATOM   4589  CB  ALA C 153      61.447  38.743 115.234 12.01 12.01           C
ATOM   4590  N   TYR C 154      64.115  39.823 114.883 14.01 14.01           N
ATOM   4591  CA  TYR C 154      65.267  40.665 115.165 12.01 12.01           C
ATOM   4592  C   TYR C 154      65.756  41.261 113.864 12.01 12.01           C
ATOM   4593  O   TYR C 154      66.400  42.302 113.847 16.00 16.00           O
ATOM   4594  CB  TYR C 154      66.376  39.866 115.831 12.01 12.01           C
ATOM   4595  CG  TYR C 154      67.016  40.606 116.973 12.01 12.01           C
ATOM   4596  CD2 TYR C 154      67.159  40.015 118.220 12.01 12.01           C
ATOM   4597  CD1 TYR C 154      67.468  41.902 116.804 12.01 12.01           C
ATOM   4598  CE2 TYR C 154      67.743  40.697 119.266 12.01 12.01           C
ATOM   4599  CE1 TYR C 154      68.055  42.593 117.842 12.01 12.01           C
ATOM   4600  CZ  TYR C 154      68.190  41.987 119.072 12.01 12.01           C
ATOM   4601  OH  TYR C 154      68.775  42.675 120.111 16.00 16.00           O
ATOM   4602  N   GLY C 155      65.432  40.583 112.771 14.01 14.01           N
ATOM   4603  CA  GLY C 155      65.788  41.039 111.449 12.01 12.01           C
ATOM   4604  C   GLY C 155      65.033  42.277 111.006 12.01 12.01           C
ATOM   4605  O   GLY C 155      64.881  43.241 111.751 16.00 16.00           O
ATOM   4606  N   GLN C 156      64.550  42.233 109.773 14.01 14.01           N
ATOM   4607  CA  GLN C 156      63.853  43.353 109.162 12.01 12.01           C
ATOM   4608  C   GLN C 156      62.850  43.968 110.123 12.01 12.01           C
ATOM   4609  O   GLN C 156      62.506  45.141 110.012 16.00 16.00           O
ATOM   4610  CB  GLN C 156      63.138  42.899 107.894 12.01 12.01           C
ATOM   4611  CG  GLN C 156      64.021  42.114 106.945 12.01 12.01           C
ATOM   4612  CD  GLN C 156      63.662  40.645 106.916 12.01 12.01           C
ATOM   4613  OE1 GLN C 156      62.516  40.271 107.160 16.00 16.00           O
ATOM   4614  NE2 GLN C 156      64.642  39.802 106.616 14.01 14.01           N
ATOM   4615  N   ALA C 157      62.387  43.173 111.075 14.01 14.01           N
ATOM   4616  CA  ALA C 157      61.412  43.650 112.046 12.01 12.01           C
ATOM   4617  C   ALA C 157      61.892  44.877 112.817 12.01 12.01           C
ATOM   4618  O   ALA C 157      61.110  45.516 113.514 16.00 16.00           O
ATOM   4619  CB  ALA C 157      61.039  42.537 113.014 12.01 12.01           C
ATOM   4620  N   LYS C 158      63.168  45.211 112.703 14.01 14.01           N
ATOM   4621  CA  LYS C 158      63.684  46.365 113.420 12.01 12.01           C
ATOM   4622  C   LYS C 158      63.183  47.677 112.829 12.01 12.01           C
ATOM   4623  O   LYS C 158      62.833  48.598 113.561 16.00 16.00           O
ATOM   4624  CB  LYS C 158      65.215  46.350 113.441 12.01 12.01           C
ATOM   4625  CG  LYS C 158      65.885  46.986 112.230 12.01 12.01           C
ATOM   4626  CD  LYS C 158      67.364  47.232 112.484 12.01 12.01           C
ATOM   4627  CE  LYS C 158      67.585  48.021 113.764 12.01 12.01           C
ATOM   4628  NZ  LYS C 158      67.696  47.123 114.945 14.01 14.01           N
ATOM   4629  N   TYR C 159      63.215  47.748 111.504 14.01 14.01           N
ATOM   4630  CA  TYR C 159      62.848  48.954 110.776 12.01 12.01           C
ATOM   4631  C   TYR C 159      61.362  49.227 110.625 12.01 12.01           C
ATOM   4632  O   TYR C 159      60.808  50.102 111.288 16.00 16.00           O
ATOM   4633  CB  TYR C 159      63.486  48.914 109.383 12.01 12.01           C
ATOM   4634  CG  TYR C 159      64.927  48.468 109.382 12.01 12.01           C
ATOM   4635  CD2 TYR C 159      65.676  48.472 110.547 12.01 12.01           C
ATOM   4636  CD1 TYR C 159      65.540  48.048 108.213 12.01 12.01           C
ATOM   4637  CE2 TYR C 159      66.996  48.072 110.545 12.01 12.01           C
ATOM   4638  CE1 TYR C 159      66.858  47.646 108.202 12.01 12.01           C
ATOM   4639  CZ  TYR C 159      67.578  47.660 109.370 12.01 12.01           C
ATOM   4640  OH  TYR C 159      68.888  47.260 109.362 16.00 16.00           O
END
"""
  m = multimer(reconstruction_type = 'ba', pdb_str=pdb_str)
  assert m.new_annotation is not None
  assert m.new_annotation.get_n_helices() == 24, \
      "expecing 24 helices, got %d" % m.new_annotation.get_n_helices()
  assert m.new_annotation.get_n_sheets() == 0, \
      "expecing 0 sheets, got %d" % m.new_annotation.get_n_helices()

def run_selected_tests():
  """  Run selected tests

  1) List in "tests" the names of the particular test you want to run
  2) Comment out unittest.main()
  3) Un-comment unittest.TextTestRunner().run(run_selected_tests())
  """
  tests = ['test_proper_biomat_application']
  suite = unittest.TestSuite(map(TestMultimerReconstruction,tests))
  return suite

if __name__=='__main__':
  exercise_ss_multiplication()
  # use for individual tests
  # unittest.TextTestRunner().run(run_selected_tests())

  # Use to run all tests
  unittest.main(verbosity=0)
