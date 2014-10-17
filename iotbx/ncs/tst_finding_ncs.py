from __future__ import division
from mmtbx.ncs import ncs_search
import iotbx.ncs
import iotbx.pdb
import unittest
import sys

class test_find_ncs_operators(unittest.TestCase):

  def test_two_chain_master(self):
    """
    Testing one ncs group and master ncs with two chains of different length
    and chains are not in order
    """
    print sys._getframe().f_code.co_name
    trans_obj = iotbx.ncs.input(pdb_string=test_pdb1)
    t = trans_obj.ncs_to_asu_selection
    self.assertEqual(t.keys(),['chain A or chain B'])
    self.assertEqual(t.values(),[['chain C or chain D', 'chain E or chain F']])

  def test_three_chain_master(self):
    """
    Test case of:
    Two groups, different number of chain in each group.
    two chains that are NOT ncs related
    """
    print sys._getframe().f_code.co_name
    trans_obj = iotbx.ncs.input(pdb_string=test_pdb2)
    t = trans_obj.ncs_to_asu_selection
    self.assertEqual(t.keys(),['chain D', 'chain A'])
    self.assertEqual(t.values(),[['chain E'], ['chain B', 'chain C']])
    #
    self.assertEqual(trans_obj.ncs_group_map[1][0],{'chain A'})
    self.assertEqual(trans_obj.ncs_group_map[2][0],{'chain D'})

  def test_master_build_from_two_related_chains(self):
    """
    Test that minimal number of chains in master ncs are selected (not the
    minimal number of transformations)
    """
    print sys._getframe().f_code.co_name
    trans_obj = iotbx.ncs.input(pdb_string=test_pdb3)
    t = trans_obj.ncs_to_asu_selection
    self.assertEqual(t.keys(),['chain A'])
    self.assertEqual(t.values(),
                     [['chain B', 'chain C', 'chain D', 'chain E',
                       'chain F', 'chain G', 'chain H', 'chain I']])

  def test_largest_common_ncs(self):
    """
    Test that minimal number of transformations in
    master ncs are selected (not the minimal number of chains)
    """
    print sys._getframe().f_code.co_name
    trans_obj = iotbx.ncs.input(
      pdb_string=test_pdb3,
      use_simple_ncs_from_pdb=False,
      use_minimal_master_ncs=False)
    t = trans_obj.ncs_to_asu_selection
    self.assertEqual(t.keys(),['chain A or chain B or chain C'])
    self.assertEqual(t.values(),
                     [['chain D or chain E or chain F',
                       'chain G or chain H or chain I']])

  def test_group_chains_by_domains_minimal_master_ncs_grouping(self):
    """ """
    print sys._getframe().f_code.co_name
    #
    pdb_inp = iotbx.pdb.input(source_info=None, lines=test_pdb3)
    ph = pdb_inp.construct_hierarchy()
    #
    chain_match_list = ncs_search.search_ncs_relations(
      ph=ph,min_contig_length=0,min_fraction_domain=0)
    # make sure that all possible chains compared
    model  = ph.models()[0]
    chain_ids = list({x.id for x in model.chains()})
    chain_ids = sorted(chain_ids)
    n_chains = len(chain_ids)
    self.assertEqual(n_chains,9)
    self.assertEqual(len(chain_match_list),36)
    #
    match_dict = ncs_search.clean_chain_matching(chain_match_list,ph)
    # Test minimal master NCS
    transform_to_group,match_dict = ncs_search.minimal_master_ncs_grouping(
    match_dict=match_dict)
    group_dict = ncs_search.build_group_dict(
      transform_to_group,match_dict)
    self.assertEqual(len(group_dict),1)
    gr_obj = group_dict[('A',)]
    self.assertEqual(len(gr_obj.transforms),len(gr_obj.copies))
    self.assertEqual(len(gr_obj.iselections),len(gr_obj.copies))
    expected = [['A'], ['B'], ['C'], ['D'], ['E'], ['F'], ['G'], ['H'], ['I']]
    self.assertEqual(gr_obj.copies,expected)
    tr = gr_obj.transforms[0]
    self.assertTrue(tr.r.is_r3_identity_matrix())
    self.assertTrue(tr.t.is_col_zero())
    tr = gr_obj.transforms[1]
    self.assertFalse(tr.r.is_r3_identity_matrix())
    self.assertFalse(tr.t.is_col_zero())



  def test_group_chains_by_domains_minimal_NCS_operators(self):
    """ Test minimal NCS operators """
    print sys._getframe().f_code.co_name

    pdb_inp = iotbx.pdb.input(source_info=None, lines=test_pdb3)
    ph = pdb_inp.construct_hierarchy()
    #
    chain_match_list = ncs_search.search_ncs_relations(
      ph=ph,min_contig_length=0,min_fraction_domain=0)
    #
    match_dict = ncs_search.clean_chain_matching(chain_match_list,ph)
    #
    transform_to_group,match_dict = ncs_search.minimal_ncs_operators_grouping(
    match_dict=match_dict)
    group_dict = ncs_search.build_group_dict(
      transform_to_group,match_dict)
    self.assertEqual(len(group_dict),1)
    gr_obj = group_dict[('A', 'B', 'C')]
    self.assertEqual(len(gr_obj.transforms),len(gr_obj.copies))
    self.assertEqual(len(gr_obj.iselections),len(gr_obj.copies))
    expected = [['A', 'B', 'C'], ['D', 'E', 'F'], ['G', 'H', 'I']]
    self.assertEqual(gr_obj.copies,expected)
    tr = gr_obj.transforms[0]
    self.assertTrue(tr.r.is_r3_identity_matrix())
    self.assertTrue(tr.t.is_col_zero())
    tr = gr_obj.transforms[1]
    self.assertFalse(tr.r.is_r3_identity_matrix())
    self.assertFalse(tr.t.is_col_zero())

test_pdb1 = '''\
CRYST1  577.812  448.715  468.790  90.00  90.00  90.00 P 1
ATOM      1  CA  LYS A 151      10.766   9.333  12.905  1.00 44.22           C
ATOM      2  CA  LYS A 152      10.117   9.159  11.610  1.00 49.42           C
ATOM      3  CA  LYS A 153       9.099   8.000  11.562  1.00 46.15           C
ATOM      4  CA  LYS A 154       8.000   8.202  11.065  1.00 52.97           C
ATOM      5  CA  LYS A 155      11.146   9.065  10.474  1.00 41.68           C
ATOM      6  CA  LYS A 156      10.547   9.007   9.084  1.00 55.55           C
TER
ATOM      7  CA  LYS B 157      11.545   9.413   8.000  1.00 72.27           C
ATOM      8  CA  LYS B 158      12.277  10.718   8.343  1.00 75.78           C
ATOM      9  CA  LYS B 159      11.349  11.791   8.809  1.00 75.88           C
TER
ATOM      7  CA  LYS F 157       2.154   3.953  16.298  1.00 72.27           C
ATOM      8  CA  LYS F 158       2.014   3.732  17.811  1.00 75.78           C
ATOM      9  CA  LYS F 159       2.558   2.413  18.250  1.00 75.88           C
TER
ATOM      7  CA  LYS D 157       4.334  10.965  12.119  1.00 72.27           C
ATOM      8  CA  LYS D 158       4.057  11.980  13.238  1.00 75.78           C
ATOM      9  CA  LYS D 159       3.177  11.427  14.310  1.00 75.88           C
TER
ATOM      1  CA  LYS C 151       6.855   8.667  15.730  1.00 44.22           C
ATOM      2  CA  LYS C 152       5.891   8.459  14.655  1.00 49.42           C
ATOM      3  CA  LYS C 153       6.103   7.155  13.858  1.00 46.15           C
ATOM      4  CA  LYS C 154       5.138   6.438  13.633  1.00 52.97           C
ATOM      5  CA  LYS C 155       5.801   9.685  13.736  1.00 41.68           C
ATOM      6  CA  LYS C 156       4.731   9.594  12.667  1.00 55.55           C
TER
ATOM      1  CA  LYS E 151       6.987   4.106  17.432  1.00 44.22           C
ATOM      2  CA  LYS E 152       6.017   3.539  16.502  1.00 49.42           C
ATOM      3  CA  LYS E 153       6.497   3.492  15.036  1.00 46.15           C
ATOM      4  CA  LYS E 154       6.348   2.458  14.400  1.00 52.97           C
ATOM      5  CA  LYS E 155       4.647   4.221  16.634  1.00 41.68           C
ATOM      6  CA  LYS E 156       3.552   3.605  15.788  1.00 55.55           C
TER
'''

test_pdb2 = '''\
CRYST1  577.812  448.715  468.790  90.00  90.00  90.00 P 1
ATOM      1  CA  LYS A 151      10.766   9.333  12.905  1.00 44.22           C
ATOM      2  CA  LYS A 152      10.117   9.159  11.610  1.00 49.42           C
ATOM      3  CA  LYS A 153       9.099   8.000  11.562  1.00 46.15           C
ATOM      4  CA  LYS A 154       8.000   8.202  11.065  1.00 52.97           C
ATOM      5  CA  LYS A 155      11.146   9.065  10.474  1.00 41.68           C
TER
ATOM    222  CA  LEU X  40      94.618  -5.253  91.582  1.00 87.10           C
ATOM    223  CA  ARG X  41      62.395  51.344  80.786  1.00107.25           C
ATOM    224  CA  ARG X  42      62.395  41.344  80.786  1.00107.25           C
TER
ATOM      1  CA  THR D   1       8.111  11.080  10.645  1.00 20.00           C
ATOM      2  CA  THR D   2       8.000   9.722  10.125  1.00 20.00           C
ATOM      3  CA  THR D   3       8.075   8.694  11.249  1.00 20.00           C
ATOM      4  CA  THR D   4       8.890   8.818  12.163  1.00 20.00           C
TER
ATOM      1  CA  LYS B 151       6.855   8.667  15.730  1.00 44.22           C
ATOM      2  CA  LYS B 152       5.891   8.459  14.655  1.00 49.42           C
ATOM      3  CA  LYS B 153       6.103   7.155  13.858  1.00 46.15           C
ATOM      4  CA  LYS B 154       5.138   6.438  13.633  1.00 52.97           C
ATOM      5  CA  LYS B 155       5.801   9.685  13.736  1.00 41.68           C
TER
ATOM      1  CA  LYS C 151       6.987   4.106  17.432  1.00 44.22           C
ATOM      2  CA  LYS C 152       6.017   3.539  16.502  1.00 49.42           C
ATOM      3  CA  LYS C 153       6.497   3.492  15.036  1.00 46.15           C
ATOM      4  CA  LYS C 154       6.348   2.458  14.400  1.00 52.97           C
ATOM      5  CA  LYS C 155       4.647   4.221  16.634  1.00 41.68           C
TER
ATOM    222  CA  LEU Y  40     194.618   5.253  81.582  1.00 87.10           C
ATOM    223  CA  ARG Y  41     162.395  41.344  70.786  1.00107.25           C
ATOM    224  CA  ARG Y  42     162.395  31.344  70.786  1.00107.25           C
TER
ATOM      1  CA  THR E   1       8.111 -10.645  11.080  1.00 20.00           C
ATOM      2  CA  THR E   2       8.000 -10.125   9.722  1.00 20.00           C
ATOM      3  CA  THR E   3       8.075 -11.249   8.694  1.00 20.00           C
ATOM      4  CA  THR E   4       8.890 -12.163   8.818  1.00 20.00           C
TER
'''

test_pdb3 = '''\
ATOM      1  N   THR A   1      13.014  18.419   8.520  1.00 20.00           N
ATOM      2  CA  THR A   1      12.903  17.061   8.000  1.00 20.00           C
ATOM      3  C   THR A   1      12.978  16.033   9.124  1.00 20.00           C
ATOM      4  O   THR A   1      13.793  16.157  10.038  1.00 20.00           O
TER
ATOM      1  N   THR C   1      10.325   8.000  14.368  1.00 20.00           N
ATOM      2  CA  THR C   1      10.111   8.702  13.108  1.00 20.00           C
ATOM      3  C   THR C   1      11.313   9.570  12.750  1.00 20.00           C
ATOM      4  O   THR C   1      11.885  10.241  13.609  1.00 20.00           O
TER
ATOM      1  N   THR D   1      -0.430  15.458  18.495  1.00 20.00           N
ATOM      2  CA  THR D   1       0.086  15.020  17.204  1.00 20.00           C
ATOM      3  C   THR D   1       1.440  14.334  17.355  1.00 20.00           C
ATOM      4  O   THR D   1       2.297  14.791  18.111  1.00 20.00           O
TER
ATOM      1  N   THR E   1       1.895   8.366  19.752  1.00 20.00           N
ATOM      2  CA  THR E   1       1.682   9.068  18.491  1.00 20.00           C
ATOM      3  C   THR E   1       2.884   9.935  18.133  1.00 20.00           C
ATOM      4  O   THR E   1       3.455  10.606  18.993  1.00 20.00           O
TER
ATOM      1  N   THR F   1       8.346   7.308  15.936  1.00 20.00           N
ATOM      2  CA  THR F   1       7.054   7.796  15.467  1.00 20.00           C
ATOM      3  C   THR F   1       6.884   9.281  15.767  1.00 20.00           C
ATOM      4  O   THR F   1       7.237   9.751  16.849  1.00 20.00           O
TER
ATOM      1  N   THR G   1       0.609  -0.560  24.094  1.00 20.00           N
ATOM      2  CA  THR G   1       0.395   0.142  22.834  1.00 20.00           C
ATOM      3  C   THR G   1       1.598   1.009  22.476  1.00 20.00           C
ATOM      4  O   THR G   1       2.169   1.680  23.335  1.00 20.00           O
TER
ATOM      1  N   THR H   1       7.061  -1.617  20.279  1.00 20.00           N
ATOM      2  CA  THR H   1       5.768  -1.130  19.810  1.00 20.00           C
ATOM      3  C   THR H   1       5.599   0.356  20.109  1.00 20.00           C
ATOM      4  O   THR H   1       5.950   0.825  21.191  1.00 20.00           O
TER
ATOM      1  N   THR I   1       8.722   4.822  16.665  1.00 20.00           N
ATOM      2  CA  THR I   1       7.494   4.036  16.653  1.00 20.00           C
ATOM      3  C   THR I   1       6.628   4.350  17.868  1.00 20.00           C
ATOM      4  O   THR I   1       7.130   4.482  18.984  1.00 20.00           O
TER
ATOM      1  N   THR B   1       8.000  15.093  13.112  1.00 20.00           N
ATOM      2  CA  THR B   1       8.516  14.654  11.820  1.00 20.00           C
ATOM      3  C   THR B   1       9.870  13.968  11.972  1.00 20.00           C
ATOM      4  O   THR B   1      10.727  14.426  12.727  1.00 20.00           O
TER
'''

def run_selected_tests():
  """  Run selected tests

  1) List in "tests" the names of the particular test you want to run
  2) Comment out unittest.main()
  3) Un-comment unittest.TextTestRunner().run(run_selected_tests())
  """
  tests = ['test_three_chain_master']
  suite = unittest.TestSuite(map(test_find_ncs_operators,tests))
  return suite

if __name__=='__main__':
  # use for individual tests
  # unittest.TextTestRunner().run(run_selected_tests())

  # Use to run all tests
  unittest.main()
