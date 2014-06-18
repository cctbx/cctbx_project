from __future__ import division
from  iotbx.pdb.multimer_reconstruction import ncs_group_object
from libtbx.utils import null_out
from scitbx import matrix
from iotbx import pdb
import unittest
import tempfile
import shutil
import os
import sys


class TestNcsGroupPreprocessing(unittest.TestCase):
  """ Test spec and phil, NCS and ASU information processing"""

  def setUp(self):
    """ Create temporary folder for temp files produced during test """
    self.currnet_dir = os.getcwd()
    self.tempdir = tempfile.mkdtemp('tempdir')
    os.chdir(self.tempdir)

  def test_phil_param_read(self):
    """ Verify that phil parameters are properly read   """
    print 'Running ',sys._getframe().f_code.co_name
    # check correctness
    expected_ncs_selection =['(chain A)','(chain A) or (chain B)']
    expected_ncs_to_asu = [
      {'chain A': ['chain B', 'chain C']},
      {'chain A': ['chain C', 'chain E'], 'chain B': ['chain D', 'chain F']}]
    expected_ncs_chains = [['chain A'],['chain A', 'chain B']]
    for i,phil_case in enumerate([user_phil1,user_phil2]):
      trans_obj = ncs_group_object()
      trans_obj.build_ncs_obj_from_phil(
        ncs_selection_params = phil_case)
      self.assertEqual(trans_obj.ncs_selection_str,expected_ncs_selection[i])
      self.assertEqual(trans_obj.ncs_to_asu_selection,expected_ncs_to_asu[i])
      self.assertEqual(trans_obj.ncs_chain_selection,expected_ncs_chains[i])
    # error reporting
    for pc in [user_phil3,user_phil4,user_phil5]:
      trans_obj = ncs_group_object()
      self.assertRaises(
        IOError,trans_obj.build_ncs_obj_from_phil,ncs_selection_params=pc)

  def test_phil_processing(self):
    """ Verify that phil parameters are properly processed   """
    print 'Running ',sys._getframe().f_code.co_name
    # read file and create pdb object
    pdb_obj = pdb.hierarchy.input(pdb_string=pdb_test_data2)
    trans_obj = ncs_group_object()
    trans_obj.build_ncs_obj_from_phil(
        ncs_selection_params = pdb_test_data2_phil,
        pdb_hierarchy_inp=pdb_obj)

    expected = '(chain A) or (chain B or chain C)'
    self.assertEqual(trans_obj.ncs_selection_str,expected)

    expected = {'chain A': ['chain D', 'chain G'],
                'chain B or chain C': ['chain E or chain F',
                                       'chain H or chain I']}
    self.assertEqual(trans_obj.ncs_to_asu_selection,expected)
    # check ncs_transform
    group_ids = [x.ncs_group_id for x in trans_obj.ncs_transform.itervalues()]
    tran_sn = {x.serial_num for x in trans_obj.ncs_transform.itervalues()}
    group_keys = {x for x in trans_obj.ncs_transform.iterkeys()}
    #
    self.assertEqual(len(group_ids),6)
    self.assertEqual(set(group_ids),{1,2})
    self.assertEqual(tran_sn,{1,2,3,4,5,6})
    self.assertEqual(group_keys,{'s005','s004','s006','s001','s003','s002'})
    self.assertEqual(trans_obj.ncs_atom_selection.count(True),4)

  def test_superpos_pdb(self):
    """  verify creation of transformations using superpose_pdb   """
    print 'Running ',sys._getframe().f_code.co_name
    # read file and create pdb object
    pdb_obj = pdb.hierarchy.input(pdb_string=pdb_test_data1)
    trans_obj = ncs_group_object()
    trans_obj.build_ncs_obj_from_phil(
        ncs_selection_params = pdb_test_data1_phil,
        pdb_hierarchy_inp=pdb_obj)

    self.assertEqual(trans_obj.ncs_selection_str,'(chain A) or (chain B)')
    expected = {'chain A': ['chain C', 'chain E'],
                'chain B': ['chain D', 'chain F']}
    self.assertEqual(trans_obj.ncs_to_asu_selection,expected)
    # check ncs_transform
    group_ids = [x.ncs_group_id for x in trans_obj.ncs_transform.itervalues()]
    tran_sn = {x.serial_num for x in trans_obj.ncs_transform.itervalues()}
    group_keys = {x for x in trans_obj.ncs_transform.iterkeys()}
    r1 = trans_obj.ncs_transform['s004'].r
    r2 = trans_obj.ncs_transform['s002'].r
    #
    self.assertEqual(len(group_ids),6)
    self.assertEqual(set(group_ids),{1,2})
    self.assertEqual(tran_sn,{1,2,3,4,5,6})
    self.assertEqual(group_keys,{'s005','s004','s006','s001','s003','s002'})
    #
    self.assertTrue(r1.is_r3_identity_matrix())
    expected_r = matrix.sqr(
      [0.309017,-0.809017,0.5,0.809017,0.5,0.309017,-0.5,0.309017,0.809017])
    d = r2 - expected_r
    d = map(abs,d)
    self.assertTrue(max(d)<0.01)

    # test that ncs_asu does not contain the identity transforms
    expected = {'chain A_002', 'chain A_003', 'chain B_005', 'chain B_006'}
    self.assertEqual(expected,set(trans_obj.ncs_to_asu_map.keys()))

    # test mapping of the different selection in the NCS
    self.assertEqual(list(trans_obj.asu_to_ncs_map['chain A']),[0,1])
    self.assertEqual(list(trans_obj.asu_to_ncs_map['chain B']),[2])

    # test that transform_chain_assignment contains all transforms
    expected = {'chain A_002', 'chain A_003', 'chain B_005', 'chain B_006'}
    self.assertEqual(expected,set(trans_obj.transform_chain_assignment))


  def test_spec_reading(self):
    """ verify creating and processing spec """
    print 'Running ',sys._getframe().f_code.co_name
    # creating a spec file
    from phenix.command_line.simple_ncs_from_pdb import simple_ncs_from_pdb
    from phenix.command_line.simple_ncs_from_pdb import ncs_master_params
    params = ncs_master_params.extract()
    params.simple_ncs_from_pdb.min_length = 1
    pdb_inp = pdb.input(source_info=None, lines=test_pdb_ncs_spec)
    ph = pdb_inp.construct_hierarchy()
    xrs = pdb_inp.xray_structure_simple()
    xrs_unit_cell = xrs.orthorhombic_unit_cell_around_centered_scatterers(
      buffer_size=8)
    ph.adopt_xray_structure(xrs_unit_cell)
    of = open("test_ncs_spec.pdb", "w")
    print >> of, ph.as_pdb_string(crystal_symmetry=xrs.crystal_symmetry())
    of.close()
    ncs_from_pdb=simple_ncs_from_pdb(
      pdb_file="test_ncs_spec.pdb",
      quiet=True,
      log=null_out(),
      params=params)

    # reading and processing the spec file
    pdb_obj = pdb.hierarchy.input(pdb_string=test_pdb_ncs_spec)
    trans_obj = ncs_group_object()
    trans_obj.build_ncs_obj_from_spec_file(
      file_name="simple_ncs_from_pdb.ncs_spec",
      # file_str=test_ncs_spec,  # use output string directly
      pdb_hierarchy_inp = pdb_obj)

    # test created object
    self.assertEqual(len(trans_obj.ncs_group),3)
    expected = '(chain A and (resseq 151:159)) or (chain D and (resseq 1:7))'
    self.assertEqual(trans_obj.ncs_selection_str,expected)
    # check that static parts are included in NCS and ASU
    self.assertEqual(len(trans_obj.ncs_atom_selection),3*9+2*7+3+3)
    self.assertEqual(trans_obj.ncs_atom_selection.count(True),9+7+3+3)
    #
    expected = {
      'chain A and (resseq 151:159)':
        ['chain B and (resseq 151:159)','chain C and (resseq 151:159)'],
      'chain D and (resseq 1:7)':
        ['chain E and (resseq 1:7)']}
    self.assertEqual(trans_obj.ncs_to_asu_selection,expected)

    # check ncs_transform
    group_ids = [x.ncs_group_id for x in trans_obj.ncs_transform.itervalues()]
    tran_sn = {x.serial_num for x in trans_obj.ncs_transform.itervalues()}
    group_keys = {x for x in trans_obj.ncs_transform.iterkeys()}
    r1 = trans_obj.ncs_transform['s004'].r
    r2 = trans_obj.ncs_transform['s002'].r
    #
    self.assertEqual(len(group_ids),5)
    self.assertEqual(set(group_ids),{1,2})
    self.assertEqual(tran_sn,{1,2,3,4,5})
    self.assertEqual(group_keys,{'s001', 's002', 's003', 's004', 's005'})
    #
    self.assertTrue(r1.is_r3_identity_matrix())
    expected_r = matrix.sqr(
      [0.4966,0.8679,-0.0102,-0.6436,0.3761,0.6666,0.5824,-0.3245,0.7453])
    d = r2 - expected_r
    d = map(abs,d)
    self.assertTrue(max(d)<0.01)


  def test_mmcif_reading(self):
    print 'Running ',sys._getframe().f_code.co_name
    # Todo: test_mmcif_reading
    pass

  def test_processing_of_asu_with_mtrix(self):
    print 'Running ',sys._getframe().f_code.co_name
    # Todo: test_processing_of_asu_with_mtrix
    pass

  def test_rotaion_translation_input(self):
    """
    Verify correct processing
    :return:
    """
  # def test_build_pdb(self):
  #   """ produce test pdb file """
  #   from iotbx import pdb
  #   os.chdir(r'C:\Phenix\Dev\Work\work\junk')
  #   pdb_inp = pdb.input(source_info=None, lines=test_ncs_spec)
  #   ph = pdb_inp.construct_hierarchy()
  #   xrs = pdb_inp.xray_structure_simple()
  #   xrs_one_ncs = xrs.orthorhombic_unit_cell_around_centered_scatterers(
  #     buffer_size=8)
  #   ph.adopt_xray_structure(xrs)
  #   of = open("test_ncs_spec.pdb", "w")
  #   print >> of, ph.as_pdb_string(crystal_symmetry=xrs.crystal_symmetry())
  #   of.close()

  def tearDown(self):
    """ remove temp files and folder """
    os.chdir(self.currnet_dir)
    shutil.rmtree(self.tempdir)


pdb_test_data1="""\
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

pdb_test_data1_phil = '''\
ncs_group_selection {
  ncs_group {
    master_ncs_selection = 'chain A'
    selection_copy = 'chain C'
    selection_copy = 'chain E'
  }
  ncs_group {
    master_ncs_selection = 'chain B'
    selection_copy = 'chain D'
    selection_copy = 'chain F'
  }
}
'''

pdb_test_data2="""\
ATOM    749  O   UNK A  90      28.392  67.262  97.682  1.00  0.00           O
ATOM    750  N   UNK B  91      30.420  66.924  98.358  1.00  0.00           N
ATOM    750  N   UNK X  93      38.420  76.924  58.358  1.00  0.00           N
ATOM   1495  N   UNK C  67      33.124   2.704 114.920  1.00  0.00           N
TER
ATOM    749  O   UNK D  90       3.199  86.786  85.616  1.00  0.00           O
ATOM    750  N   UNK E  91       4.437  88.467  85.044  1.00  0.00           N
ATOM   1495  N   UNK F  67      65.508  63.662  77.246  1.00  0.00           N
TER
ATOM    749  O   UNK G  90     -26.415  72.437  94.483  1.00  0.00           O
ATOM    750  N   UNK H  91     -27.678  74.103  93.921  1.00  0.00           N
ATOM   1495  N   UNK I  67       7.362 108.699  49.412  1.00  0.00           N
TER
"""

pdb_test_data2_phil = '''\
ncs_group_selection {
  ncs_group {
    master_ncs_selection = 'chain A'
    selection_copy = 'chain D'
    selection_copy = 'chain G'
  }
  ncs_group {
    master_ncs_selection = 'chain B or chain C'
    selection_copy = 'chain E or chain F'
    selection_copy = 'chain H or chain I'
  }
}
'''

user_phil1 = '''\
ncs_group_selection {
  ncs_group {
    master_ncs_selection = 'chain A'
    selection_copy = 'chain B'
    selection_copy = 'chain C'
  }
}
'''

user_phil2 = '''\
ncs_group_selection {
  ncs_group {
    master_ncs_selection = 'chain A'
    selection_copy = 'chain C'
    selection_copy = 'chain E'
  }
  ncs_group {
    master_ncs_selection = 'chain B'
    selection_copy = 'chain D'
    selection_copy = 'chain F'
  }
}
'''

user_phil3 = '''\
ncs_group_selection {
  ncs_group {
    master_ncs_selection = 'chain A'
    selection_copy = 'chain B'
    selection_copy = 'chain C'
  }
  ncs_group {
    master_ncs_selection = 'chain B'
    selection_copy = 'chain D'
  }
}
'''

user_phil4 = '''\
ncs_group_selection {
  ncs_group {
    master_ncs_selection = 'chain A'
    selection_copy = 'chain C'
    selection_copy = 'chain D'
  }
  ncs_group {
    master_ncs_selection = 'chain B'
    selection_copy = 'chain D'
  }
}
'''

user_phil5 = '''\
ncs_group_selection {
  ncs_group {
    master_ncs_selection = 'chain A'
    selection_copy = 'chain C'
    selection_copy = 'chain D'
  }
  ncs_group {
    master_ncs_selection = 'chain C'
    selection_copy = 'chain E'
  }
}
'''

test_pdb_ncs_spec = '''\
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
test_ncs_spec = '''\

Summary of NCS information
Fri Jun 13 13:18:12 2014
C:\Phenix\Dev\Work\work\junk

new_ncs_group
new_operator

rota_matrix    1.0000    0.0000    0.0000
rota_matrix    0.0000    1.0000    0.0000
rota_matrix    0.0000    0.0000    1.0000
tran_orth     0.0000    0.0000    0.0000

center_orth   10.5384    9.4098   10.2058
CHAIN A
RMSD 0.0
MATCHING 9
  RESSEQ 151:159

new_operator

rota_matrix    0.4966    0.8679   -0.0102
rota_matrix   -0.6436    0.3761    0.6666
rota_matrix    0.5824   -0.3245    0.7453
tran_orth    -0.0003   -0.0002    0.0003

center_orth    5.1208    9.3744   13.7718
CHAIN B
RMSD 0.000490318808456
MATCHING 9
  RESSEQ 151:159

new_operator

rota_matrix   -0.3180    0.7607    0.5659
rota_matrix   -0.1734   -0.6334    0.7541
rota_matrix    0.9321    0.1416    0.3334
tran_orth     0.0002    0.0004   -0.0006

center_orth    4.5304    3.5021   16.4612
CHAIN C
RMSD 0.000475964539946
MATCHING 9
  RESSEQ 151:159

new_ncs_group
new_operator

rota_matrix    1.0000    0.0000    0.0000
rota_matrix    0.0000    1.0000    0.0000
rota_matrix    0.0000    0.0000    1.0000
tran_orth     0.0000    0.0000    0.0000

center_orth    8.5917    9.4397    9.9770
CHAIN D
RMSD 0.0
MATCHING 7
  RESSEQ 1:7

new_operator

rota_matrix    1.0000    0.0000    0.0000
rota_matrix    0.0000   -0.0000    1.0000
rota_matrix    0.0000   -1.0000   -0.0000
tran_orth     0.0000   -0.0000    0.0000

center_orth    8.5917   -9.9770    9.4397
CHAIN E
RMSD 7.50647770368e-16
MATCHING 7
  RESSEQ 1:7




'''

if __name__ == "__main__":
  unittest.main()







