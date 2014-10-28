from __future__ import division
from mmtbx.ncs.ncs_search import res_alignment
from scitbx.array_family import flex
from mmtbx.ncs import ncs_search
from libtbx.utils import Sorry
from iotbx import pdb
import unittest
import sys


class TestSimpleAlignment(unittest.TestCase):
  """ Test alignment of closely similar sequences  """

  def setUp(self):
    a = "ssssAESSADKFKRQxxxHMDTEGPSKSSPTYCNQMMKRQGMTKGSCKPVNTFVHEPLEDVQ" \
        "NGRNNCHKSSSTLRITDCRLKGSSKYPNCDYTTTDSQKHIIIACDGNPYVPVHFDASV"
    b = "AESSADKFKRQHMDTEGPSKSSPTYCNQMMKRQGMTKGSCKPVNTFVHEPLEDVQ" \
        "NGRNNCHKSSSTLRITDCRLKGSSKYPNCDYTTTDSQkhIIIACDGNPYVPVHFDASVtttt"
    self.seq_a = list(a)
    self.seq_b = list(b)
    self.length_diff = len(b)/len(a)
    # Gaps needed for the inserted, not aligned, letters
    self.gaps_needed = 4 + 3 + 4
    #
    pdb_obj = pdb.hierarchy.input(pdb_string=pdb_str)
    self.ph = pdb_obj.hierarchy
    cache = self.ph.atom_selection_cache()
    self.chain_a = self.ph.models()[0].chains()[0]
    self.chain_b = self.ph.models()[0].chains()[1]
    self.hierarchy_a = self.ph.select(cache.selection('chain A'))
    self.hierarchy_b = self.ph.select(cache.selection('chain B'))

  def test_1(self):
    print sys._getframe().f_code.co_name
    a = 'abcfadx'
    b = 'cabfa'
    seq_a = list(a)
    seq_b = list(b)
    # Test that aligned segments are at least min_contig_length long
    sel_a, sel_b, similarity = res_alignment(
      seq_a=seq_a, seq_b=seq_b,
      min_percent =0.1,min_contig_length=3)
    self.assertEqual([],list(sel_a))
    self.assertEqual([],list(sel_b))
    # Without length limitation
    sel_a, sel_b,similarity = res_alignment(
      seq_a=seq_a, seq_b=seq_b,
      min_percent =0.1,min_contig_length=2)
    self.assertEqual([0,1,3,4],list(sel_a))
    self.assertEqual([1,2,3,4],list(sel_b))

  def test_2(self):
    print sys._getframe().f_code.co_name
    sel_a, sel_b, similarity = res_alignment(
      self.seq_a, self.seq_b, min_percent=0.9,min_contig_length=0)
    expected_1 = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 18, 19,
                  20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
                  32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43,
                  44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55,
                  56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67,
                  68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
                  80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91,
                  92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102,
                  103, 104, 105, 106, 107, 108, 109, 110, 111, 112,
                  113, 114, 115, 116, 117, 118, 119]
    expected_2 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
                  14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
                  26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37,
                  38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
                  50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61,
                  62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73,
                  74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85,
                  86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97,
                  98, 99, 100, 101, 102, 103, 104, 105, 106, 107,
                  108, 109, 110, 111, 112]
    self.assertEqual(expected_1,list(sel_a))
    self.assertEqual(expected_2,list(sel_b))

  def test_3(self):
    print sys._getframe().f_code.co_name
    sel_a, sel_b, similarity = res_alignment(
      self.seq_a,self.seq_b,min_percent=0.95,min_contig_length=0)
    # difference is to large
    expected_1 = []
    expected_2 = []
    self.assertEqual(expected_1,list(sel_a))
    self.assertEqual(expected_2,list(sel_b))

  def test_search_ncs_relations(self):
    print sys._getframe().f_code.co_name
    chain_match_list = ncs_search.search_ncs_relations(
      ph=self.ph,
      min_percent=0.70,
      min_contig_length=3,
      check_atom_order=True)
    [chain_a_id,chain_b_id,sel_a,sel_b,r1,r2,_] = chain_match_list[0]
    #
    self.assertEqual(chain_a_id,'A')
    self.assertEqual(chain_b_id,'B')
    #
    atoms_in_A = self.hierarchy_a.atoms().size()
    atoms_in_B = self.hierarchy_b.atoms().size()
    # atom count including  water
    self.assertEqual(32,atoms_in_A)
    self.assertEqual(44,atoms_in_B)
    #
    self.assertEqual(sel_a.size(),25)
    self.assertEqual(sel_b.size(),25)
    # atom count including  water
    self.assertEqual(32,atoms_in_A)
    self.assertEqual(44,atoms_in_B)
    #
    # number of waters is different as well
    expected_sel_a = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                      16, 17, 21, 22, 23, 24, 25, 26, 27]
    expected_sel_b = [28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41,
                      42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52]
    self.assertEqual(list(sel_a), expected_sel_a)
    self.assertEqual(list(sel_b), expected_sel_b)
    #
    self.assertEqual(r1,[0, 1, 2, 3])
    self.assertEqual(r2,[0, 1, 2, 3])

  def test_get_chains_info(self):
    print sys._getframe().f_code.co_name
    pdb_inp = pdb.hierarchy.input(pdb_string=test_pdb_1)
    ph = pdb_inp.hierarchy
    # test without selection
    chains_info1 = ncs_search.get_chains_info(ph)
    self.assertEqual(sorted(chains_info1),['A', 'D', 'X'])

    # test with selection strings: compare to no selection
    selection_list = ['chain A','chain x','chain D']
    chains_info2 = ncs_search.get_chains_info(ph,selection_list)
    self.assertEqual(sorted(chains_info2),['chain A', 'chain D', 'chain x'])

    self.assertEqual(
      chains_info1['A'].atom_names,chains_info2['chain A'].atom_names)
    self.assertEqual(
      chains_info1['A'].atom_selection,chains_info2['chain A'].atom_selection)
    self.assertEqual(
      chains_info1['A'].chains_atom_number,
      chains_info2['chain A'].chains_atom_number)
    self.assertEqual(
      chains_info1['X'].res_names,chains_info2['chain x'].res_names)
    self.assertEqual(
      chains_info1['X'].resid,chains_info2['chain x'].resid)

    # test with selection strings: multiple chain selection
    selection_list = [
      '(chain A and resseq 151:156) or (chain x and resname LEU)','chain D']
    chains_info = ncs_search.get_chains_info(ph,selection_list)
    self.assertEqual(
      sorted(chains_info),
      ['(chain A and resseq 151:156) or (chain x and resname LEU)', 'chain D'])

    key = '(chain A and resseq 151:156) or (chain x and resname LEU)'
    expected = [[' CA '],[' CA '],[' CA '],[' CA '],[' CA '],[' CA '],[' CA ']]
    self.assertEqual(chains_info[key].atom_names,expected)
    self.assertEqual(
      chains_info[key].atom_selection,[[0],[1],[2],[3],[4],[5],[9]])
    self.assertEqual(chains_info[key].chains_atom_number,7)
    self.assertEqual(
      chains_info[key].res_names,['LYS','LYS','LYS','LYS','LYS','LYS','LEU'])
    self.assertEqual(
      chains_info[key].resid,
      [' 151 ',' 152 ',' 153 ',' 154 ',' 155 ',' 156 ','  40 '])

    # test with selection strings: overlapping selection
    selection_list = ['chain A and resseq 151:156','chain A or chain D']
    self.assertRaises(Sorry,ncs_search.get_chains_info,
                      ph=ph,
                      selection_list=selection_list)

    # test with selection strings: empty selection
    selection_list = ['chain A and resseq 1:6']
    self.assertRaises(Sorry,ncs_search.get_chains_info,
                      ph=ph,
                      selection_list=selection_list)

  def test_remove_masters_if_appear_in_copies(self):
    print sys._getframe().f_code.co_name
    transform_to_group = {1:[['A','B','E'],['B','C','F']]}
    tg = ncs_search.remove_masters_if_appear_in_copies(transform_to_group)
    self.assertEqual(tg.values(),[[['A','E'],['B','F']]])
    #
    transform_to_group = {
    1:[['A'],['B']],2:[['B'],['C']],3:[['A'],['D']],4:[['A'],['E']],
    5:[['A'],['C']],6:[['B'],['D']],7:[['G'],['J']],8:[['C'],['H']],
    9:[['A','1'],['X','2']]}
    tg = ncs_search.remove_masters_if_appear_in_copies(transform_to_group)
    expected = [
      [['A'], ['B']], [['A'], ['C']], [['A'], ['D']], [['A'], ['E']],
      [['A', '1'], ['X', '2']], [['G'], ['J']]]
    self.assertEqual(sorted(tg.values()),expected)

  def test_remove_far_atoms(self):
    print sys._getframe().f_code.co_name
    pdb_inp = pdb.hierarchy.input(pdb_string=AB112_116)
    ph = pdb_inp.hierarchy
    pdb_inp_b = pdb.hierarchy.input(pdb_string=B112_116_fitted)
    ph_b = pdb_inp_b.hierarchy
    info_ab =  ncs_search.get_chains_info(ph)
    # matching residues
    res_n_a = [0, 1, 2, 3, 4, 5, 6, 7, 8]
    res_n_b = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    sel_a = flex.size_t([x for xi in info_ab['A'].atom_selection for x in xi])
    other_sites = ph.select(sel_a).atoms().extract_xyz()
    ref_sites = ph_b.atoms().extract_xyz()
    self.assertEqual(other_sites.size(),ref_sites.size())
    self.assertEqual(other_sites.size(),56)
    #
    sel_a,sel_b,res_n_a,res_n_b,ref_sites,other_sites = \
      ncs_search.remove_far_atoms(
        'A','B',
        res_n_a,res_n_b,
        info_ab,
        ref_sites,other_sites,
        max_dist_diff=4.0)
    #
    self.assertEqual([0, 1, 5, 7, 8],res_n_a)
    self.assertEqual([1, 2, 6, 8, 9],res_n_b)

    pass

test_pdb_1 = '''\
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
END
'''

pdb_str = '''\
CRYST1   37.760   43.710  107.440  90.00 108.54  90.00 P 1 21 1      4
ATOM      1  N   GLY A   1      10.100  14.506  -6.813  1.00 11.02           N
ATOM      2  CA  GLY A   1       9.618  15.911  -6.685  1.00  9.34           C
ATOM      3  C   GLY A   1      10.369  16.727  -7.707  1.00  7.85           C
ATOM      4  O   GLY A   1      11.481  16.388  -8.042  1.00  8.86           O
ATOM      5  N   LYS A   2       9.783  17.791  -8.213  1.00  7.53           N
ATOM      6  CA  LYS A   2      10.481  18.573  -9.209  1.00  7.43           C
ATOM      7  C   LYS A   2      10.199  20.063  -9.068  1.00  6.50           C
ATOM      8  O   LYS A   2       9.065  20.460  -8.837  1.00  4.58           O
ATOM      9  CB  LYS A   2      10.064  18.096 -10.596  1.00 10.22           C
ATOM     10  CG  LYS A   2      10.871  18.687 -11.722  1.00 11.46           C
ATOM     11  CD  LYS A   2      10.092  18.615 -13.014  1.00 19.76           C
ATOM     12  CE  LYS A   2       9.725  17.195 -13.381  1.00 20.40           C
ATOM     13  NZ  LYS A   2      10.922  16.324 -13.540  1.00 26.54           N
ATOM     14  N   ILE A   3      11.226  20.875  -9.282  1.00  6.47           N
ATOM     15  CA  ILE A   3      11.128  22.321  -9.186  1.00  4.49           C
ATOM     16  C   ILE A   3      12.072  22.938 -10.224  1.00  4.38           C
ATOM     17  O   ILE A   3      13.107  22.372 -10.517  1.00  3.75           O
ATOM     18  CB  ILE A   3      11.550  22.813  -7.771  1.00  2.00           C
ATOM     19  CG1 ILE A   3      11.272  24.308  -7.617  1.00  2.44           C
ATOM     20  CG2 ILE A   3      13.048  22.554  -7.541  1.00  2.00           C
ATOM     21  CD1 ILE A   3      11.435  24.816  -6.208  1.00  2.18           C
ATOM     22  N   THR A   4      11.699  24.079 -10.788  1.00  4.32           N
ATOM     23  CA  THR A   4      12.543  24.760 -11.751  1.00  7.18           C
ATOM     24  C   THR A   4      12.633  26.205 -11.328  1.00  5.50           C
ATOM     25  O   THR A   4      11.602  26.824 -11.045  1.00  5.90           O
ATOM     26  CB  THR A   4      11.950  24.736 -13.155  1.00  8.96           C
ATOM     27  OG1 THR A   4      11.629  23.389 -13.508  1.00 14.82           O
ATOM     28  CG2 THR A   4      12.954  25.288 -14.146  1.00 10.81           C
TER
ATOM   1492  N   GLY B   3      30.298  10.660  57.402  1.00  8.41           N
ATOM   1493  CA  GLY B   3      30.111   9.211  57.319  1.00  6.72           C
ATOM   1494  C   GLY B   3      30.875   8.547  58.431  1.00  6.27           C
ATOM   1495  O   GLY B   3      31.907   9.064  58.860  1.00  6.95           O
ATOM   1496  N   LYS B   4      30.353   7.443  58.946  1.00  6.16           N
ATOM   1497  CA  LYS B   4      31.044   6.737  60.009  1.00  9.57           C
ATOM   1498  C   LYS B   4      30.699   5.272  59.908  1.00  8.82           C
ATOM   1499  O   LYS B   4      29.544   4.909  59.671  1.00 10.99           O
ATOM   1500  CB  LYS B   4      30.681   7.309  61.385  1.00 15.37           C
ATOM   1501  CG  LYS B   4      31.286   6.548  62.561  1.00 17.25           C
ATOM   1502  CD  LYS B   4      30.583   6.873  63.875  1.00 19.93           C
ATOM   1503  CE  LYS B   4      30.944   8.264  64.374  1.00 25.63           C
ATOM   1504  NZ  LYS B   4      32.411   8.433  64.713  1.00 32.30           N
ATOM   1505  N   ILE B   5      31.719   4.436  60.038  1.00  6.78           N
ATOM   1506  CA  ILE B   5      31.566   3.002  59.948  1.00  2.83           C
ATOM   1507  C   ILE B   5      32.501   2.424  60.999  1.00  3.23           C
ATOM   1508  O   ILE B   5      33.516   3.037  61.326  1.00  5.50           O
ATOM   1509  CB  ILE B   5      31.947   2.522  58.541  1.00  2.00           C
ATOM   1513  N   THR B   6      32.140   1.279  61.559  1.00  2.00           N
ATOM   1514  CA  THR B   6      32.944   0.619  62.593  1.00  2.91           C
ATOM   1515  C   THR B   6      33.142  -0.836  62.182  1.00  2.66           C
ATOM   1516  O   THR B   6      32.180  -1.560  61.883  1.00  2.00           O
ATOM   1517  CB  THR B   6      32.226   0.656  63.970  1.00  5.04           C
ATOM   1518  OG1 THR B   6      31.894   2.006  64.303  1.00  7.93           O
ATOM   1519  CG2 THR B   6      33.099   0.080  65.064  1.00  5.59           C
ATOM   1520  N   PHE B   7      34.393  -1.258  62.142  1.00  3.13           N
ATOM   1521  CA  PHE B   7      34.709  -2.610  61.743  1.00  4.16           C
ATOM   1522  C   PHE B   7      35.004  -3.435  62.989  1.00  7.18           C
ATOM   1523  O   PHE B   7      35.698  -2.951  63.889  1.00  8.29           O
ATOM   1524  CB  PHE B   7      35.936  -2.591  60.832  1.00  4.23           C
ATOM   1525  CG  PHE B   7      35.889  -1.521  59.766  1.00  5.16           C
ATOM   1526  CD1 PHE B   7      35.279  -1.760  58.545  1.00  4.71           C
ATOM   1527  CD2 PHE B   7      36.474  -0.268  59.990  1.00  4.91           C
ATOM   1528  CE1 PHE B   7      35.248  -0.782  57.557  1.00  5.29           C
ATOM   1529  CE2 PHE B   7      36.447   0.709  59.016  1.00  3.36           C
ATOM   1530  CZ  PHE B   7      35.826   0.443  57.791  1.00  8.22           C
TER
HETATM 2976  O   HOH A 200      -5.015  35.539  13.500  1.00 13.39           O
HETATM 2977  O   HOH A 201      27.397  33.353  -7.193  1.00 42.92           O
HETATM 2978  O   HOH A 202       9.999  12.737  -8.741  1.00 22.09           O
HETATM 2979  O   HOH A 203      20.918  31.300   9.008  1.00 11.57           O
HETATM 3128  O   HOH B 303      38.892  -4.641  42.351  1.00 28.05           O
HETATM 3129  O   HOH B 304      23.863   8.189  55.409  1.00  9.33           O
HETATM 3130  O   HOH B 305      37.020   0.527  42.707  1.00  6.39           O
HETATM 3131  O   HOH B 306      29.889  -7.894  21.701  1.00 17.26           O
HETATM 3132  O   HOH B 307      36.038  -0.748  40.367  1.00  8.62           O
HETATM 3133  O   HOH B 308      31.840 -18.209  58.265  1.00 34.09           O
HETATM 3134  O   HOH B 309      21.378  -6.597  55.543  1.00 21.44           O
HETATM 3135  O   HOH B 310      39.908 -13.756  54.831  1.00 11.05           O
END
'''

AB112_116 = '''\
CRYST1  435.970  183.000  225.390  90.00 108.99  90.00 C 1 2 1      48
SCALE1      0.002294  0.000000  0.000789        0.00000
SCALE2      0.000000  0.005464  0.000000        0.00000
SCALE3      0.000000  0.000000  0.004692        0.00000
ATOM    642  N   GLY A 112     -74.574 -54.455  36.175  1.00 79.40           N
ATOM    643  CA  GLY A 112     -74.465 -54.666  37.607  1.00 83.00           C
ATOM    644  C   GLY A 112     -73.722 -55.943  37.951  1.00 84.28           C
ATOM    645  O   GLY A 112     -73.439 -56.211  39.118  1.00 87.64           O
ATOM    646  N   ALA A 113     -73.414 -56.738  36.932  1.00 83.50           N
ATOM    647  CA  ALA A 113     -72.596 -57.931  37.110  1.00 82.93           C
ATOM    648  C   ALA A 113     -71.163 -57.534  37.451  1.00 80.54           C
ATOM    649  O   ALA A 113     -70.557 -56.721  36.754  1.00 78.62           O
ATOM    650  CB  ALA A 113     -72.628 -58.785  35.855  1.00 83.57           C
ATOM    651  N   PRO A 114     -70.616 -58.108  38.532  1.00 83.82           N
ATOM    652  CA  PRO A 114     -69.275 -57.762  39.016  1.00 84.57           C
ATOM    653  C   PRO A 114     -68.151 -58.247  38.106  1.00 85.91           C
ATOM    654  O   PRO A 114     -67.833 -59.436  38.098  1.00 88.66           O
ATOM    655  CB  PRO A 114     -69.190 -58.486  40.366  1.00 87.53           C
ATOM    656  CG  PRO A 114     -70.605 -58.786  40.742  1.00 88.24           C
ATOM    657  CD  PRO A 114     -71.307 -59.029  39.447  1.00 86.78           C
ATOM    658  N   ILE A 115     -67.558 -57.328  37.350  1.00 89.78           N
ATOM    659  CA  ILE A 115     -66.334 -57.615  36.614  1.00 92.13           C
ATOM    660  C   ILE A 115     -65.174 -57.526  37.600  1.00 98.02           C
ATOM    661  O   ILE A 115     -64.008 -57.698  37.241  1.00 96.36           O
ATOM    662  CB  ILE A 115     -66.127 -56.622  35.451  1.00 89.89           C
ATOM    663  CG1 ILE A 115     -64.961 -57.062  34.562  1.00 92.90           C
ATOM    664  CG2 ILE A 115     -65.917 -55.210  35.982  1.00 87.17           C
ATOM    665  CD1 ILE A 115     -64.727 -56.159  33.371  1.00 91.28           C
ATOM    666  N   ASP A 116     -65.525 -57.262  38.855  1.00104.50           N
ATOM    667  CA  ASP A 116     -64.564 -57.096  39.937  1.00112.90           C
ATOM    668  C   ASP A 116     -63.761 -58.370  40.179  1.00119.75           C
ATOM    669  O   ASP A 116     -62.562 -58.318  40.455  1.00117.63           O
ATOM    670  CB  ASP A 116     -65.301 -56.701  41.219  1.00115.98           C
ATOM    671  CG  ASP A 116     -64.408 -55.988  42.212  1.00120.07           C
ATOM    672  OD1 ASP A 116     -63.170 -56.093  42.083  1.00122.69           O
ATOM    673  OD2 ASP A 116     -64.944 -55.323  43.123  1.00120.16           O
ATOM    674  N   GLY A 117     -64.431 -59.513  40.073  1.00119.73           N
ATOM    675  CA  GLY A 117     -63.801 -60.795  40.328  1.00125.54           C
ATOM    676  C   GLY A 117     -63.734 -61.135  41.806  1.00127.67           C
ATOM    677  O   GLY A 117     -63.286 -62.219  42.179  1.00131.50           O
ATOM    678  N   LYS A 118     -64.179 -60.210  42.650  1.00123.44           N
ATOM    679  CA  LYS A 118     -64.151 -60.418  44.095  1.00122.70           C
ATOM    680  C   LYS A 118     -65.232 -61.380  44.595  1.00122.68           C
ATOM    681  O   LYS A 118     -65.086 -61.981  45.659  1.00129.18           O
ATOM    682  CB  LYS A 118     -64.228 -59.082  44.840  1.00119.98           C
ATOM    683  CG  LYS A 118     -62.940 -58.274  44.781  1.00118.22           C
ATOM    684  CD  LYS A 118     -63.051 -56.980  45.569  1.00118.70           C
ATOM    685  CE  LYS A 118     -61.757 -56.180  45.501  1.00119.82           C
ATOM    686  NZ  LYS A 118     -61.886 -54.842  46.145  1.00119.54           N
ATOM    687  N   GLY A 119     -66.313 -61.522  43.834  1.00119.90           N
ATOM    688  CA  GLY A 119     -67.364 -62.458  44.194  1.00119.27           C
ATOM    689  C   GLY A 119     -68.748 -62.047  43.729  1.00116.99           C
ATOM    690  O   GLY A 119     -68.897 -61.427  42.676  1.00113.26           O
ATOM    691  N   PRO A 120     -69.777 -62.412  44.510  1.00119.66           N
ATOM    692  CA  PRO A 120     -71.160 -61.997  44.260  1.00115.62           C
ATOM    693  C   PRO A 120     -71.347 -60.537  44.650  1.00108.77           C
ATOM    694  O   PRO A 120     -70.388 -59.894  45.078  1.00114.46           O
ATOM    695  CB  PRO A 120     -71.961 -62.899  45.196  1.00122.11           C
ATOM    696  CG  PRO A 120     -71.037 -63.146  46.333  1.00125.70           C
ATOM    697  CD  PRO A 120     -69.662 -63.237  45.726  1.00124.48           C
TER
ATOM   4292  N   LEU B 111     -74.170 -53.490  32.745  1.00118.81           N
ATOM   4293  CA  LEU B 111     -75.497 -54.028  33.028  1.00124.23           C
ATOM   4294  C   LEU B 111     -75.744 -54.116  34.532  1.00120.20           C
ATOM   4295  O   LEU B 111     -76.714 -53.552  35.020  1.00124.85           O
ATOM   4296  CB  LEU B 111     -75.771 -55.355  32.297  1.00128.26           C
ATOM   4297  CG  LEU B 111     -75.161 -56.698  32.702  1.00130.01           C
ATOM   4298  CD1 LEU B 111     -75.848 -57.275  33.931  1.00135.24           C
ATOM   4299  CD2 LEU B 111     -75.278 -57.670  31.533  1.00128.53           C
ATOM   4300  N   GLY B 112    -111.059 -18.236  -2.457  1.00108.88           N
ATOM   4301  CA  GLY B 112    -111.950 -17.386  -3.228  1.00103.58           C
ATOM   4302  C   GLY B 112    -112.916 -16.511  -2.453  1.00100.28           C
ATOM   4303  O   GLY B 112    -112.527 -15.793  -1.532  1.00 94.05           O
ATOM   4304  N   ALA B 113    -114.190 -16.591  -2.829  1.00113.45           N
ATOM   4305  CA  ALA B 113    -115.241 -15.733  -2.282  1.00128.46           C
ATOM   4306  C   ALA B 113    -115.192 -15.490  -0.769  1.00141.15           C
ATOM   4307  O   ALA B 113    -115.505 -14.387  -0.326  1.00134.98           O
ATOM   4308  CB  ALA B 113    -116.623 -16.246  -2.698  1.00132.53           C
ATOM   4309  N   PRO B 114    -114.804 -16.509   0.026  1.00158.95           N
ATOM   4310  CA  PRO B 114    -114.735 -16.318   1.477  1.00164.22           C
ATOM   4311  C   PRO B 114    -114.476 -14.863   1.850  1.00156.52           C
ATOM   4312  O   PRO B 114    -113.387 -14.327   1.638  1.00159.74           O
ATOM   4313  CB  PRO B 114    -113.570 -17.220   1.866  1.00170.89           C
ATOM   4314  CG  PRO B 114    -113.739 -18.401   0.932  1.00174.42           C
ATOM   4315  CD  PRO B 114    -114.398 -17.878  -0.341  1.00170.85           C
ATOM   4316  N   ILE B 115    -115.508 -14.239   2.409  1.00146.97           N
ATOM   4317  CA  ILE B 115    -115.596 -12.788   2.496  1.00134.85           C
ATOM   4318  C   ILE B 115    -115.260 -12.189   3.859  1.00129.85           C
ATOM   4319  O   ILE B 115    -115.262 -12.879   4.879  1.00127.28           O
ATOM   4320  CB  ILE B 115    -117.002 -12.313   2.086  1.00133.00           C
ATOM   4321  CG1 ILE B 115    -118.058 -13.325   2.542  1.00135.78           C
ATOM   4322  CG2 ILE B 115    -117.086 -12.147   0.582  1.00131.22           C
ATOM   4323  CD1 ILE B 115    -118.178 -13.482   4.047  1.00135.38           C
ATOM   4324  N   ASP B 116    -114.985 -10.886   3.848  1.00124.18           N
ATOM   4325  CA  ASP B 116    -114.709 -10.108   5.049  1.00115.99           C
ATOM   4326  C   ASP B 116    -114.924  -8.631   4.722  1.00103.80           C
ATOM   4327  O   ASP B 116    -114.195  -8.062   3.910  1.00107.79           O
ATOM   4328  CB  ASP B 116    -113.263 -10.311   5.511  1.00113.78           C
ATOM   4329  CG  ASP B 116    -112.825 -11.765   5.469  1.00109.78           C
ATOM   4330  OD1 ASP B 116    -112.961 -12.458   6.499  1.00108.02           O
ATOM   4331  OD2 ASP B 116    -112.344 -12.215   4.406  1.00108.69           O
ATOM   4332  N   GLY B 117    -115.923  -8.013   5.345  1.00101.31           N
ATOM   4333  CA  GLY B 117    -116.795  -8.697   6.279  1.00 98.56           C
ATOM   4334  C   GLY B 117    -118.204  -8.856   5.743  1.00116.01           C
ATOM   4335  O   GLY B 117    -119.000  -9.613   6.299  1.00110.03           O
ATOM   4336  N   LYS B 118    -118.520  -8.138   4.666  1.00130.29           N
ATOM   4337  CA  LYS B 118    -119.836  -8.247   4.040  1.00154.94           C
ATOM   4338  C   LYS B 118    -120.063  -9.655   3.495  1.00156.32           C
ATOM   4339  O   LYS B 118    -119.113 -10.404   3.274  1.00163.40           O
ATOM   4340  CB  LYS B 118    -120.015  -7.202   2.933  1.00168.24           C
ATOM   4341  CG  LYS B 118    -119.996  -5.763   3.429  1.00171.89           C
ATOM   4342  CD  LYS B 118    -121.137  -4.941   2.834  1.00177.02           C
ATOM   4343  CE  LYS B 118    -120.967  -4.729   1.336  1.00176.66           C
ATOM   4344  NZ  LYS B 118    -122.043  -3.862   0.774  1.00180.46           N
ATOM   4345  N   GLY B 119    -121.326 -10.010   3.283  1.00157.13           N
ATOM   4346  CA  GLY B 119    -121.693 -11.366   2.912  1.00152.24           C
ATOM   4347  C   GLY B 119    -121.499 -11.774   1.461  1.00155.42           C
ATOM   4348  O   GLY B 119    -120.887 -12.806   1.186  1.00145.45           O
ATOM   4349  N   PRO B 120    -122.036 -10.981   0.521  1.00169.33           N
ATOM   4350  CA  PRO B 120    -122.021 -11.349  -0.899  1.00176.41           C
ATOM   4351  C   PRO B 120    -120.755 -10.917  -1.634  1.00171.02           C
ATOM   4352  O   PRO B 120    -119.877 -10.280  -1.053  1.00173.27           O
ATOM   4353  CB  PRO B 120    -123.215 -10.576  -1.453  1.00187.63           C
ATOM   4354  CG  PRO B 120    -123.225  -9.326  -0.639  1.00188.14           C
ATOM   4355  CD  PRO B 120    -122.778  -9.727   0.750  1.00180.96           C
END
'''

B112_116_fitted = '''\
CRYST1  435.970  183.000  225.390  90.00 108.99  90.00 C 1 2 1
SCALE1      0.002294  0.000000  0.000789        0.00000
SCALE2      0.000000  0.005464  0.000000        0.00000
SCALE3      0.000000  0.000000  0.004692        0.00000
ATOM   4300  N   GLY B 112     -74.275 -54.368  36.100  1.00108.88           N
ATOM   4301  CA  GLY B 112     -74.350 -55.041  37.385  1.00103.58           C
ATOM   4302  C   GLY B 112     -73.130 -55.823  37.831  1.00100.28           C
ATOM   4303  O   GLY B 112     -72.008 -55.318  37.814  1.00 94.05           O
ATOM   4304  N   ALA B 113     -73.362 -57.075  38.218  1.00113.45           N
ATOM   4305  CA  ALA B 113     -72.332 -57.938  38.794  1.00128.46           C
ATOM   4306  C   ALA B 113     -70.954 -57.881  38.125  1.00141.15           C
ATOM   4307  O   ALA B 113     -69.943 -57.969  38.819  1.00134.98           O
ATOM   4308  CB  ALA B 113     -72.828 -59.385  38.876  1.00132.53           C
ATOM   4309  N   PRO B 114     -70.904 -57.736  36.784  1.00158.95           N
ATOM   4310  CA  PRO B 114     -69.608 -57.668  36.104  1.00164.22           C
ATOM   4311  C   PRO B 114     -68.514 -57.121  37.014  1.00156.52           C
ATOM   4312  O   PRO B 114     -68.494 -55.939  37.360  1.00159.74           O
ATOM   4313  CB  PRO B 114     -69.903 -56.728  34.942  1.00170.89           C
ATOM   4314  CG  PRO B 114     -71.315 -57.115  34.553  1.00174.42           C
ATOM   4315  CD  PRO B 114     -72.009 -57.615  35.816  1.00170.85           C
ATOM   4316  N   ILE B 115     -67.606 -58.015  37.392  1.00146.97           N
ATOM   4317  CA  ILE B 115     -66.716 -57.801  38.525  1.00134.85           C
ATOM   4318  C   ILE B 115     -65.291 -57.384  38.172  1.00129.85           C
ATOM   4319  O   ILE B 115     -64.835 -57.558  37.042  1.00127.28           O
ATOM   4320  CB  ILE B 115     -66.654 -59.065  39.401  1.00133.00           C
ATOM   4321  CG1 ILE B 115     -66.741 -60.321  38.529  1.00135.78           C
ATOM   4322  CG2 ILE B 115     -67.793 -59.071  40.401  1.00131.22           C
ATOM   4323  CD1 ILE B 115     -65.576 -60.512  37.574  1.00135.38           C
ATOM   4324  N   ASP B 116     -65.771 -57.644  39.567  1.00124.18           N
ATOM   4325  CA  ASP B 116     -64.325 -57.578  39.395  1.00115.99           C
ATOM   4326  C   ASP B 116     -63.713 -57.054  40.694  1.00103.80           C
ATOM   4327  O   ASP B 116     -63.949 -55.907  41.073  1.00107.79           O
ATOM   4328  CB  ASP B 116     -63.954 -56.631  38.251  1.00113.78           C
ATOM   4329  CG  ASP B 116     -64.825 -56.816  37.020  1.00109.78           C
ATOM   4330  OD1 ASP B 116     -64.438 -57.601  36.129  1.00108.02           O
ATOM   4331  OD2 ASP B 116     -65.896 -56.176  36.943  1.00108.69           O
ATOM   4332  N   GLY B 117     -62.935 -57.889  41.375  1.00101.31           N
ATOM   4333  CA  GLY B 117     -62.684 -59.247  40.936  1.00 98.56           C
ATOM   4334  C   GLY B 117     -63.330 -60.276  41.843  1.00116.01           C
ATOM   4335  O   GLY B 117     -63.416 -61.452  41.489  1.00110.03           O
ATOM   4336  N   LYS B 118     -63.779 -59.840  43.020  1.00130.29           N
ATOM   4337  CA  LYS B 118     -64.457 -60.736  43.953  1.00154.94           C
ATOM   4338  C   LYS B 118     -65.758 -61.263  43.352  1.00156.32           C
ATOM   4339  O   LYS B 118     -66.301 -60.676  42.418  1.00163.40           O
ATOM   4340  CB  LYS B 118     -64.721 -60.047  45.297  1.00168.24           C
ATOM   4341  CG  LYS B 118     -63.458 -59.663  46.055  1.00171.89           C
ATOM   4342  CD  LYS B 118     -63.535 -60.063  47.527  1.00177.02           C
ATOM   4343  CE  LYS B 118     -64.580 -59.254  48.284  1.00176.66           C
ATOM   4344  NZ  LYS B 118     -64.598 -59.595  49.736  1.00180.46           N
ATOM   4345  N   GLY B 119     -66.252 -62.372  43.892  1.00157.13           N
ATOM   4346  CA  GLY B 119     -67.397 -63.062  43.322  1.00152.24           C
ATOM   4347  C   GLY B 119     -68.776 -62.494  43.612  1.00155.42           C
ATOM   4348  O   GLY B 119     -69.562 -62.279  42.690  1.00145.45           O
ATOM   4349  N   PRO B 120     -69.090 -62.267  44.896  1.00169.33           N
ATOM   4350  CA  PRO B 120     -70.437 -61.844  45.298  1.00176.41           C
ATOM   4351  C   PRO B 120     -70.647 -60.333  45.254  1.00171.02           C
ATOM   4352  O   PRO B 120     -69.724 -59.579  44.947  1.00173.27           O
ATOM   4353  CB  PRO B 120     -70.516 -62.323  46.745  1.00187.63           C
ATOM   4354  CG  PRO B 120     -69.117 -62.161  47.239  1.00188.14           C
ATOM   4355  CD  PRO B 120     -68.218 -62.486  46.066  1.00180.96           C
END
'''

def run_selected_tests():
  """  Run selected tests

  1) List in "tests" the names of the particular test you want to run
  2) Comment out unittest.main()
  3) Un-comment unittest.TextTestRunner().run(run_selected_tests())
  """
  tests = ['test_remove_far_atoms']
  suite = unittest.TestSuite(map(TestSimpleAlignment,tests))
  return suite

if __name__=='__main__':
  # use for individual tests
  # unittest.TextTestRunner().run(run_selected_tests())

  # Use to run all tests
  unittest.main()
