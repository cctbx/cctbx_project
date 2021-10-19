from __future__ import absolute_import, division, print_function
from mmtbx.ncs.ncs_search import mmtbx_res_alignment
from scitbx.array_family import flex
from mmtbx.ncs import ncs_search
import iotbx.ncs as ncs
import iotbx.pdb
import unittest
from iotbx.pdb.amino_acid_codes import three_letter_given_one_letter
from six.moves import map

__author__ = 'Youval'

class TestSimpleAlignment(unittest.TestCase):
  """ Test alignment of closely similar sequences  """

  def setUp(self):
    a = "ssssAESSADKFKRQxxxHMDTEGPSKSSPTYCNQMMKRQGMTKGSCKPVNTFVHEPLEDVQ" \
        "NGRNNCHKSSSTLRITDCRLKGSSKYPNCDYTTTDSQKHIIIACDGNPYVPVHFDASV"
    b = "AESSADKFKRQHMDTEGPSKSSPTYCNQMMKRQGMTKGSCKPVNTFVHEPLEDVQ" \
        "NGRNNCHKSSSTLRITDCRLKGSSKYPNCDYTTTDSQkhIIIACDGNPYVPVHFDASVtttt"
    self.seq_a = [three_letter_given_one_letter.get(x, 'X') for x in a]
    self.seq_b = [three_letter_given_one_letter.get(x, 'X') for x in b]
    self.length_diff = len(b)/len(a)
    # Gaps needed for the inserted, not aligned, letters
    self.gaps_needed = 4 + 3 + 4
    #
    pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
    self.ph = pdb_inp.construct_hierarchy()
    cache = self.ph.atom_selection_cache()
    self.chain_a = self.ph.models()[0].chains()[0]
    self.chain_b = self.ph.models()[0].chains()[1]
    self.hierarchy_a = self.ph.select(cache.selection('chain A'))
    self.hierarchy_b = self.ph.select(cache.selection('chain B'))
    self.p = ncs.input.get_default_params()

  def test_1(self):
    # print sys._getframe().f_code.co_name
    seq_aa = "AAKDVKFGVNVLADAVKVTVASKANDAAGDGTTTATVLAQAI"
    seq_bb = "AAKDVKFGNDARVKMLRGVNVLADAVKVTVASTATVLAQAI"
    seq_a = [three_letter_given_one_letter.get(x, 'X') for x in seq_aa]
    seq_b = [three_letter_given_one_letter.get(x, 'X') for x in seq_bb]
    # Test that aligned segments are at least min_contig_length long
    sel_a, sel_b, similarity = mmtbx_res_alignment(
      seq_a=seq_a, seq_b=seq_b,
      min_percent =0.1)
    self.assertEqual([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                16, 17, 18, 19, 20, 21, 33, 34, 35, 36, 37, 38, 39, 40, 41],
                    list(sel_a))
    self.assertEqual([0, 1, 2, 3, 4, 5, 6, 7, 18, 19, 20, 21, 22, 23, 24,
            25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40],
                    list(sel_b))
    # Without length limitation
    sel_a, sel_b,similarity = mmtbx_res_alignment(
      seq_a=seq_a, seq_b=seq_b,
      min_percent =0.8)
    self.assertEqual([],list(sel_a))
    self.assertEqual([],list(sel_b))

  def test_2(self):
    # print sys._getframe().f_code.co_name
    seq_a = ["CB", "NB", "OB", "SB", "CB1", "NB1", "OB1", "SB1"]
    seq_b = ["CB", "NB", "OB", "SB",        "NB1", "OB1", "SB1"]
    sel_a, sel_b, similarity = mmtbx_res_alignment(
      seq_a, seq_b, min_percent=0.6, atomnames=True)
    expected_1 = [0,1,2,3,5,6,7]
    expected_2 = [0,1,2,3,4,5,6]
    # print "sela", list(sel_a)
    self.assertEqual(similarity, 0.875)
    self.assertEqual(expected_1,list(sel_a))
    self.assertEqual(expected_2,list(sel_b))

  def test_3(self):
    # print sys._getframe().f_code.co_name
    sel_a, sel_b, similarity = mmtbx_res_alignment(
      self.seq_a,self.seq_b,min_percent=0.95)
    # difference is to large
    expected_1 = []
    expected_2 = []
    self.assertEqual(expected_1,list(sel_a))
    self.assertEqual(expected_2,list(sel_b))

  def test_search_ncs_relations(self):
    # print sys._getframe().f_code.co_name
    asc = self.ph.atom_selection_cache()
    no_water_h = self.ph.select(asc.selection("not water"))
    match_dict = ncs_search.search_ncs_relations(
      ph=no_water_h,
      chain_similarity_threshold=0.7)

    # FIXME indexing keys breaks py2/3 compatibility, this could be broken if more than 1 key
    chain_a_id, chain_b_id = list(match_dict.keys())[0]
    sel_a,sel_b,r1,r2,_,_,_ = match_dict[chain_a_id, chain_b_id]
    #
    self.assertEqual(chain_a_id,'A')
    self.assertEqual(chain_b_id,'B')
    #
    atoms_in_A = self.hierarchy_a.atoms_size()
    atoms_in_B = self.hierarchy_b.atoms_size()
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
    # print sys._getframe().f_code.co_name
    pdb_inp = iotbx.pdb.input(source_info=None, lines=test_pdb_1)
    ph = pdb_inp.construct_hierarchy()
    # test without selection
    chains_info1 = ncs_search.get_chains_info(ph)
    self.assertEqual(sorted(chains_info1),['A', 'D', 'X'])

  def test_remove_far_atoms(self):
    # print sys._getframe().f_code.co_name
    pdb_inp = iotbx.pdb.input(source_info=None, lines=AB112_116)
    ph = pdb_inp.construct_hierarchy()
    pdb_inp_b = iotbx.pdb.input(source_info=None, lines=B112_116_fitted)
    ph_b = pdb_inp_b.construct_hierarchy()
    info_ab =  ncs_search.get_chains_info(ph)
    # matching residues
    res_n_a = [0, 1, 2, 3, 4, 5, 6, 7, 8]
    res_n_b = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    list_a = [flex.size_t(x) for x in info_ab['A'].atom_selection]
    list_b = [flex.size_t(x) for x in info_ab['B'].atom_selection]
    xxx = [x for y in list_a for x in y]
    xxx.sort()
    sel_a = flex.size_t(xxx)
    other_sites = ph.select(sel_a).atoms().extract_xyz()
    ref_sites = ph_b.atoms().extract_xyz()
    self.assertEqual(other_sites.size(),ref_sites.size())
    self.assertEqual(other_sites.size(),56)
    #
    sel_a,sel_b,res_n_a,res_n_b,ref_sites,other_sites = \
      ncs_search.remove_far_atoms(
        list_a, list_b,
        res_n_a,res_n_b,
        ref_sites,other_sites,
        residue_match_radius=4.0)
    #
    self.assertEqual([0, 1, 5, 7, 8],res_n_a)
    self.assertEqual([1, 2, 6, 8, 9],res_n_b)

  def test_update_atom_selections(self):
    # print sys._getframe().f_code.co_name
    pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str_3)
    ph = pdb_inp.construct_hierarchy(sort_atoms=False)
    chains_info = ncs_search.get_chains_info(ph)
    match_dict = ncs_search.search_ncs_relations(
        ph=ph,
        chains_info=chains_info,
        chain_similarity_threshold=0.10,
        chain_max_rmsd=2.3)
    nrgl = ncs_search.ncs_grouping_and_group_dict(match_dict, ph)
    # print group_dict
    nrgl.update_str_selections_if_needed(hierarchy=ph)
    # nrgl._show()
    strsels = nrgl.get_array_of_str_selections()
    # print strsels
    assert strsels == [[
        "(chain 'A' and resid 2 through 216)",
        "chain 'B'",
        "chain 'C'",
        "(chain 'D' and (resid 2 through 171 or resid 186 through 216))"]]
    assert nrgl.get_n_groups() == 1
    isels = nrgl[0].get_iselections_list()
    isels_list = [list(a) for a in isels]
    # print isels_list
    assert isels_list == [
        [9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
         27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44,
         45, 46, 47, 48, 49, 50, 51],
        [52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69,
         70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87,
         88, 89, 90, 91, 92, 93, 94],
        [95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109,
         110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123,
         124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137],

        [138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151,
         152, 153, 154, 155, 156, 157, 158, 159, 160, 165, 166, 167, 168, 169,
         170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184]]

  def test_split_chain_with_altloc(self):
    # print sys._getframe().f_code.co_name
    pdb_inp = iotbx.pdb.input(source_info=None, lines=test_pdb_4)
    chains_info = ncs_search.get_chains_info(pdb_inp.construct_hierarchy())
    ch_A = chains_info['A']
    self.assertEqual(len(ch_A.res_names),len(ch_A.no_altloc))
    #
    pdb_inp = iotbx.pdb.input(source_info=None, lines=test_pdb_5)
    chains_info = ncs_search.get_chains_info(pdb_inp.construct_hierarchy())
    ch_A = chains_info['A']
    self.assertEqual(len(ch_A.res_names),len(ch_A.no_altloc))

  def test_groups_with_chains_of_different_size(self):
    # print sys._getframe().f_code.co_name
    pdb_inp = iotbx.pdb.input(source_info=None, lines=test_pdb_6)
    ph=pdb_inp.construct_hierarchy()
    nrgl = ncs_search.find_ncs_in_hierarchy(ph=ph)
    nrgl.update_str_selections_if_needed(hierarchy=ph)
    str_sels = nrgl.get_array_of_str_selections()
    self.assertEqual(str_sels,
        [["(chain 'H' and (resid 5 through 6 or (resid 640 and (name C1 ))))",
          "chain 'J'", "(chain 'L' and (resid 5 through 6 or (resid 646 and (name C1 ))))"],
        ["chain 'I'", "chain 'K'", "chain 'M'"]])


  def test_iselection_is_in_correct_order(self):
    """ Make sure can calling get_ncs_restraints_group_list does not raise
    error when atom selection in NCS group are out of order"""
    h = iotbx.pdb.input(source_info=None, lines=test_pdb_6).construct_hierarchy()
    ncs_obj = ncs.input(hierarchy=h)
    ncs_obj.get_ncs_restraints_group_list()
    pass


  def test_correct_transform_selection(self):
    """
    When [A,B] can be the master for [D,E], [F,G]
    but the transform from B to G is the same as the inverse of A to D,
    make sure the non-transform match is selected
    """
    # print sys._getframe().f_code.co_name
    pdb_inp = iotbx.pdb.input(source_info=None, lines=test_pdb_7)
    ph = pdb_inp.construct_hierarchy()
    chains_info = ncs_search.get_chains_info(ph)
    match_dict = ncs_search.search_ncs_relations(ph=ph,
      chains_info=chains_info,chain_similarity_threshold=0.10)
    nrgl = ncs_search.ncs_grouping_and_group_dict(match_dict, ph)
    nrgl.update_str_selections_if_needed(hierarchy=ph)
    str_sels = nrgl.get_array_of_str_selections()
    self.assertEqual(str_sels,
        [["chain 'A'", "chain 'D'", "chain 'F'"],
        ["chain 'B'", "chain 'E'", "chain 'G'"]])
    # nrgl._show()
    #
    # self.assertEqual(group_dict.keys(),[('A',), ('B',)])

  def test_split_groups_to_spec(self):
    # print sys._getframe().f_code.co_name
    pdb_str = test_pdb_7 + test_pdb_7_addition
    h = iotbx.pdb.input(source_info=None, lines=pdb_str).construct_hierarchy()
    ncs_obj = ncs.input(hierarchy=h)
    spec = ncs_obj.get_ncs_info_as_spec()
    gr = spec.ncs_groups()
    self.assertEqual(len(gr),2)
    self.assertEqual(gr[0].chain_residue_id()[0],["A", "D", "F"])
    self.assertEqual(gr[1].chain_residue_id()[0],["B", "E", "G", "X"])

  def test_processing_ncs_with_hierarchy_input(self):
    pdb_inp = iotbx.pdb.input(lines=test_pdb_6,source_info=None)
    ph = pdb_inp.construct_hierarchy()
    ncs_obj = ncs.input(hierarchy=ph)
    nrgl = ncs_obj.get_ncs_restraints_group_list()
    self.assertEqual(nrgl.get_n_groups(), 2)
    self.assertEqual(nrgl[0].get_number_of_copies(), 2)
    self.assertEqual(nrgl[1].get_number_of_copies(), 2)

  def test_min_contig_length(self):
    """ Test correct handling of min_contig_length and min_percent
    alignment with 0.73 idenitiy is really exist:
    referenc      AAKDVKFG----------VNVLADAVKVTVASKANDAAGDGTTTATVLAQAI
                  ||||||||          ||||||||||||||           |||||||||
    query         AAKDVKFGNDARVKMLRGVNVLADAVKVTVAS-----------TATVLAQAI
    """
    h = iotbx.pdb.input(source_info=None, lines=test_pdb_8).construct_hierarchy()
    self.p.ncs_search.chain_max_rmsd=10
    self.p.ncs_search.residue_match_radius=1000
    self.p.ncs_search.chain_similarity_threshold=0.5
    ncs_obj = ncs.input(
        hierarchy=h,
        params=self.p.ncs_search)
    self.assertEqual(ncs_obj.number_of_ncs_groups,1)
    ncs_groups = ncs_obj.get_ncs_restraints_group_list()
    assert len(ncs_groups) == 1
    assert ncs_groups[0].master_iselection.size() == 31
    assert ncs_groups[0].copies[0].iselection.size() == 31

  def test_file_formats(self):
    """ Make sure .ncs file format is as needed """
    # fixme: finish test when we decide how to change .ncs files
    # create a pdb string for the test and flip to atoms order
    pdb_str = test_pdb_7 + test_pdb_7_addition
    pdb_lines = pdb_str.splitlines()
    pdb_lines[2:4] = [
      'ATOM      3  O   ASP A   5      89.553 -29.857  72.225  1.00 77.56           O',
      'ATOM      4  C   ASP A   5      90.136 -30.762  71.617  1.00 77.70           C'
    ]
    h = iotbx.pdb.input(source_info=None, lines='\n'.join(pdb_lines)).construct_hierarchy()
    self.p.ncs_search.chain_max_rmsd=10
    self.p.ncs_search.chain_similarity_threshold=0.20
    ncs_obj = ncs.input(
        hierarchy=h,
        params=self.p.ncs_search)
    # x = ncs_obj.get_ncs_info_as_spec(write=True,write_ncs_phil=True)
    # x = ncs_obj.show(format='spec')
    # check another pdb string
    h = iotbx.pdb.input(source_info=None, lines=test_pdb_6).construct_hierarchy()
    ncs_obj = ncs.input(hierarchy=h)
    ncs_obj.get_ncs_restraints_group_list()
    # x = ncs_obj.get_ncs_info_as_spec(write=True,write_ncs_phil=True)

  def test_overlapping_chains(self):
    '''  Test processing of overlapping identical chains  '''
    h = iotbx.pdb.input(source_info=None, lines=pdb_str10).construct_hierarchy()
    ncs_inp = ncs.input(hierarchy=h)
    nrgl = ncs_inp.get_ncs_restraints_group_list()
    # nrgl._show()
    self.assertEqual(nrgl[0].master_str_selection, "chain 'A '")


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

pdb_str_2 = '''\
CRYST1   53.930  108.200  147.160  90.00  90.00  90.00 P 21 21 21   16
SCALE1      0.018543  0.000000  0.000000        0.00000
SCALE2      0.000000  0.009242  0.000000        0.00000
SCALE3      0.000000  0.000000  0.006795        0.00000
ATOM      1  N   GLN A  -1     -11.627  10.458 -69.449  1.00 66.71           N
ATOM      2  CA  GLN A  -1     -12.125  11.752 -68.882  1.00 71.10           C
ATOM      3  C   GLN A  -1     -11.010  12.528 -68.080  1.00 70.59           C
ATOM      4  O   GLN A  -1      -9.956  12.796 -68.658  1.00 72.03           O
ATOM      5  CB  GLN A  -1     -13.372  11.454 -68.047  1.00 72.25           C
ATOM      6  CG  GLN A  -1     -14.331  12.612 -67.869  1.00 72.60           C
ATOM      7  CD  GLN A  -1     -15.171  12.445 -66.641  1.00 68.34           C
ATOM      8  OE1 GLN A  -1     -14.698  12.648 -65.520  1.00 64.02           O
ATOM      9  NE2 GLN A  -1     -16.418  12.080 -66.838  1.00 67.72           N
ATOM     10  N   SER A   0     -11.274  12.959 -66.830  1.00 66.92           N
ATOM     11  CA  SER A   0     -10.269  13.265 -65.790  1.00 65.35           C
ATOM     12  C   SER A   0     -10.382  12.220 -64.655  1.00 63.87           C
ATOM     13  O   SER A   0      -9.628  12.224 -63.681  1.00 69.02           O
ATOM     14  CB  SER A   0     -10.492  14.676 -65.206  1.00 69.65           C
ATOM     15  OG  SER A   0     -11.321  14.650 -64.034  1.00 62.15           O
ATOM     16  N   MET A   1     -11.405  11.381 -64.772  1.00 60.87           N
ATOM     17  CA  MET A   1     -11.551  10.116 -64.057  1.00 54.12           C
ATOM     18  C   MET A   1     -10.548   9.045 -64.525  1.00 48.89           C
ATOM     19  O   MET A   1     -10.553   7.938 -63.982  1.00 43.78           O
ATOM     20  CB  MET A   1     -12.984   9.610 -64.327  1.00 54.20           C
ATOM     21  CG  MET A   1     -14.056  10.134 -63.367  1.00 53.02           C
ATOM     22  SD  MET A   1     -13.776   9.668 -61.636  1.00 51.21           S
ATOM     23  CE  MET A   1     -14.325   7.960 -61.693  1.00 52.78           C
ATOM     24  N   SER A   2      -9.665   9.382 -65.470  1.00 40.80           N
ATOM     25  CA  SER A   2      -8.906   8.392 -66.236  1.00 37.26           C
ATOM     26  C   SER A   2      -7.689   7.905 -65.466  1.00 34.89           C
ATOM     27  O   SER A   2      -6.972   8.684 -64.878  1.00 37.04           O
ATOM     28  CB  SER A   2      -8.477   8.972 -67.586  1.00 35.56           C
ATOM     29  OG  SER A   2      -7.350   8.293 -68.131  1.00 34.81           O
ATOM     30  N   LEU A   3      -7.454   6.604 -65.494  1.00 35.57           N
ATOM     31  CA  LEU A   3      -6.272   6.004 -64.855  1.00 36.13           C
ATOM     32  C   LEU A   3      -5.327   5.402 -65.901  1.00 37.73           C
ATOM     33  O   LEU A   3      -4.604   4.461 -65.615  1.00 34.22           O
ATOM     34  CB  LEU A   3      -6.723   4.909 -63.909  1.00 35.14           C
ATOM     35  CG  LEU A   3      -7.869   5.303 -62.975  1.00 34.10           C
ATOM     36  CD1 LEU A   3      -8.258   4.166 -62.074  1.00 33.01           C
ATOM     37  CD2 LEU A   3      -7.489   6.520 -62.156  1.00 37.25           C
ATOM     38  N   GLN A   4      -5.331   5.970 -67.110  1.00 43.39           N
ATOM     39  CA  GLN A   4      -4.570   5.405 -68.225  1.00 46.98           C
ATOM     40  C   GLN A   4      -3.121   5.461 -67.776  1.00 41.88           C
ATOM     41  O   GLN A   4      -2.709   6.446 -67.224  1.00 41.47           O
ATOM     42  CB  GLN A   4      -4.792   6.221 -69.516  1.00 52.77           C
ATOM     43  CG  GLN A   4      -4.258   5.608 -70.825  1.00 62.32           C
ATOM     44  CD  GLN A   4      -5.145   4.498 -71.405  1.00 70.16           C
ATOM     45  OE1 GLN A   4      -6.073   4.004 -70.746  1.00 75.80           O
ATOM     46  NE2 GLN A   4      -4.848   4.090 -72.643  1.00 67.99           N
ATOM     47  N   GLY A   5      -2.382   4.380 -67.922  1.00 42.31           N
ATOM     48  CA  GLY A   5      -0.969   4.392 -67.606  1.00 43.10           C
ATOM     49  C   GLY A   5      -0.693   4.036 -66.168  1.00 45.85           C
ATOM     50  O   GLY A   5       0.469   4.059 -65.725  1.00 42.79           O
ATOM     51  N   LYS A   6      -1.733   3.711 -65.403  1.00 42.68           N
ATOM     52  CA  LYS A   6      -1.459   3.217 -64.065  1.00 41.28           C
ATOM     53  C   LYS A   6      -1.774   1.775 -63.808  1.00 34.91           C
ATOM     54  O   LYS A   6      -2.720   1.196 -64.355  1.00 38.16           O
ATOM     55  CB  LYS A   6      -2.156   4.078 -63.032  1.00 44.16           C
ATOM     56  CG  LYS A   6      -1.870   5.543 -63.255  1.00 41.61           C
ATOM     57  CD  LYS A   6      -1.933   6.314 -61.986  1.00 41.03           C
ATOM     58  CE  LYS A   6      -1.903   7.786 -62.314  1.00 42.29           C
ATOM     59  NZ  LYS A   6      -0.581   8.346 -62.089  1.00 47.26           N
ATOM     60  N   VAL A   7      -0.958   1.197 -62.943  1.00 28.08           N
ATOM     61  CA  VAL A   7      -1.094  -0.189 -62.571  1.00 25.74           C
ATOM     62  C   VAL A   7      -1.829  -0.205 -61.232  1.00 24.26           C
ATOM     63  O   VAL A   7      -1.364   0.414 -60.252  1.00 23.35           O
ATOM     64  CB  VAL A   7       0.295  -0.849 -62.426  1.00 25.45           C
ATOM     65  CG1 VAL A   7       0.193  -2.277 -61.911  1.00 24.74           C
ATOM     66  CG2 VAL A   7       1.029  -0.815 -63.765  1.00 27.28           C
ATOM     67  N   ALA A   8      -2.910  -0.967 -61.176  1.00 22.55           N
ATOM     68  CA  ALA A   8      -3.650  -1.120 -59.945  1.00 22.99           C
ATOM     69  C   ALA A   8      -3.615  -2.568 -59.493  1.00 23.34           C
ATOM     70  O   ALA A   8      -3.778  -3.464 -60.296  1.00 23.15           O
ATOM     71  CB  ALA A   8      -5.078  -0.656 -60.119  1.00 21.71           C
TER    1758      SER A 247
ATOM   1760  CA  ASN B  -5     -32.763   5.274   2.167  1.00 76.86           C
ATOM   1761  C   ASN B  -5     -32.157   6.017   0.968  1.00 73.12           C
ATOM   1762  O   ASN B  -5     -31.225   5.512   0.325  1.00 66.84           O
ATOM   1763  CB  ASN B  -5     -34.276   5.540   2.248  1.00 85.21           C
ATOM   1764  CG  ASN B  -5     -34.902   5.011   3.533  1.00 92.73           C
ATOM   1765  OD1 ASN B  -5     -34.335   4.150   4.222  1.00 96.14           O
ATOM   1766  ND2 ASN B  -5     -36.076   5.537   3.868  1.00 93.99           N
ATOM   1767  N   LEU B  -4     -32.681   7.215   0.680  1.00 69.04           N
ATOM   1768  CA  LEU B  -4     -32.100   8.095  -0.347  1.00 66.62           C
ATOM   1769  C   LEU B  -4     -30.656   8.416   0.031  1.00 68.55           C
ATOM   1770  O   LEU B  -4     -29.765   8.408  -0.816  1.00 63.84           O
ATOM   1771  CB  LEU B  -4     -32.866   9.429  -0.472  1.00 64.15           C
ATOM   1772  CG  LEU B  -4     -32.175  10.515  -1.327  1.00 62.37           C
ATOM   1773  CD1 LEU B  -4     -31.835  10.043  -2.744  1.00 60.90           C
ATOM   1774  CD2 LEU B  -4     -33.008  11.780  -1.389  1.00 59.79           C
ATOM   1775  N   TYR B  -3     -30.440   8.711   1.306  1.00 64.98           N
ATOM   1776  CA  TYR B  -3     -29.117   9.171   1.763  1.00 67.34           C
ATOM   1777  C   TYR B  -3     -28.007   8.102   1.965  1.00 65.79           C
ATOM   1778  O   TYR B  -3     -26.830   8.428   1.845  1.00 58.29           O
ATOM   1779  CB  TYR B  -3     -29.262  10.057   3.008  1.00 67.25           C
ATOM   1780  CG  TYR B  -3     -29.855  11.415   2.711  1.00 70.44           C
ATOM   1781  CD1 TYR B  -3     -29.079  12.412   2.112  1.00 68.62           C
ATOM   1782  CD2 TYR B  -3     -31.185  11.716   3.026  1.00 71.73           C
ATOM   1783  CE1 TYR B  -3     -29.603  13.668   1.843  1.00 68.12           C
ATOM   1784  CE2 TYR B  -3     -31.716  12.971   2.753  1.00 72.09           C
ATOM   1785  CZ  TYR B  -3     -30.919  13.946   2.160  1.00 69.15           C
ATOM   1786  OH  TYR B  -3     -31.432  15.200   1.886  1.00 61.93           O
ATOM   1787  N   PHE B  -2     -28.385   6.858   2.262  1.00 75.41           N
ATOM   1788  CA  PHE B  -2     -27.460   5.743   2.479  1.00 79.50           C
ATOM   1789  C   PHE B  -2     -27.673   4.666   1.409  1.00 75.48           C
ATOM   1790  O   PHE B  -2     -28.637   3.899   1.427  1.00 70.21           O
ATOM   1791  CB  PHE B  -2     -27.647   5.195   3.892  1.00 81.95           C
ATOM   1792  CG  PHE B  -2     -27.604   6.251   4.964  1.00 85.79           C
ATOM   1793  CD1 PHE B  -2     -26.400   6.668   5.504  1.00 86.27           C
ATOM   1794  CD2 PHE B  -2     -28.781   6.820   5.446  1.00 93.44           C
ATOM   1795  CE1 PHE B  -2     -26.369   7.628   6.502  1.00 93.03           C
ATOM   1796  CE2 PHE B  -2     -28.751   7.783   6.449  1.00 99.13           C
ATOM   1797  CZ  PHE B  -2     -27.542   8.181   6.982  1.00 96.60           C
ATOM   1798  N   GLN B  -1     -26.754   4.667   0.452  1.00 78.86           N
ATOM   1799  CA  GLN B  -1     -26.814   3.855  -0.754  1.00 84.42           C
ATOM   1800  C   GLN B  -1     -25.455   3.220  -1.060  1.00 82.15           C
ATOM   1801  O   GLN B  -1     -25.354   2.182  -1.756  1.00 70.29           O
ATOM   1802  CB  GLN B  -1     -27.137   4.792  -1.922  1.00 83.03           C
ATOM   1803  CG  GLN B  -1     -28.574   4.845  -2.404  1.00 83.73           C
ATOM   1804  CD  GLN B  -1     -28.770   5.981  -3.394  1.00 85.99           C
ATOM   1805  OE1 GLN B  -1     -27.881   6.261  -4.192  1.00 78.88           O
ATOM   1806  NE2 GLN B  -1     -29.922   6.656  -3.332  1.00 87.79           N
ATOM   1807  N   SER B   0     -24.442   3.833  -0.441  1.00 84.33           N
ATOM   1808  CA  SER B   0     -23.223   4.266  -1.163  1.00 83.58           C
ATOM   1809  C   SER B   0     -22.542   3.354  -2.160  1.00 81.31           C
ATOM   1810  O   SER B   0     -22.481   2.138  -2.031  1.00 64.12           O
ATOM   1811  CB  SER B   0     -22.129   4.807  -0.226  1.00 81.41           C
ATOM   1812  OG  SER B   0     -22.202   6.219  -0.147  1.00 78.43           O
ATOM   1813  N   MET B   1     -21.988   4.060  -3.135  1.00 84.64           N
ATOM   1814  CA  MET B   1     -21.048   3.590  -4.104  1.00 91.39           C
ATOM   1815  C   MET B   1     -19.622   3.737  -3.534  1.00 84.93           C
ATOM   1816  O   MET B   1     -18.750   4.408  -4.100  1.00 83.07           O
ATOM   1817  CB  MET B   1     -21.251   4.465  -5.329  1.00 95.18           C
ATOM   1818  CG  MET B   1     -22.638   4.345  -5.917  1.00 90.71           C
ATOM   1819  SD  MET B   1     -23.168   5.933  -6.562  1.00 95.39           S
ATOM   1820  CE  MET B   1     -21.782   6.648  -7.464  1.00 77.69           C
ATOM   1821  N   SER B   2     -19.435   3.129  -2.367  1.00 68.33           N
ATOM   1822  CA  SER B   2     -18.122   2.823  -1.841  1.00 56.23           C
ATOM   1823  C   SER B   2     -17.667   1.496  -2.480  1.00 49.42           C
ATOM   1824  O   SER B   2     -18.495   0.790  -3.086  1.00 47.15           O
ATOM   1825  CB  SER B   2     -18.223   2.683  -0.315  1.00 61.95           C
ATOM   1826  OG  SER B   2     -18.227   1.337   0.111  1.00 58.21           O
ATOM   1827  N   LEU B   3     -16.383   1.136  -2.352  1.00 39.00           N
ATOM   1828  CA  LEU B   3     -15.872  -0.124  -2.895  1.00 36.22           C
ATOM   1829  C   LEU B   3     -15.419  -1.064  -1.800  1.00 34.80           C
ATOM   1830  O   LEU B   3     -14.566  -1.939  -2.035  1.00 34.78           O
ATOM   1831  CB  LEU B   3     -14.700   0.131  -3.826  1.00 35.39           C
ATOM   1832  CG  LEU B   3     -14.944   1.249  -4.837  1.00 33.93           C
ATOM   1833  CD1 LEU B   3     -13.728   1.448  -5.727  1.00 34.16           C
ATOM   1834  CD2 LEU B   3     -16.173   0.974  -5.697  1.00 35.81           C
ATOM   1835  N   GLN B   4     -16.041  -0.944  -0.627  1.00 35.17           N
ATOM   1836  CA  GLN B   4     -15.615  -1.715   0.531  1.00 38.05           C
ATOM   1837  C   GLN B   4     -15.840  -3.138   0.162  1.00 36.52           C
ATOM   1838  O   GLN B   4     -16.852  -3.448  -0.431  1.00 38.48           O
ATOM   1839  CB  GLN B   4     -16.403  -1.379   1.809  1.00 43.87           C
ATOM   1840  CG  GLN B   4     -16.314   0.095   2.254  1.00 47.31           C
ATOM   1841  CD  GLN B   4     -14.943   0.501   2.804  1.00 53.62           C
ATOM   1842  OE1 GLN B   4     -14.059  -0.348   3.027  1.00 59.13           O
ATOM   1843  NE2 GLN B   4     -14.767   1.802   3.059  1.00 54.06           N
ATOM   1844  N   GLY B   5     -14.866  -3.998   0.411  1.00 38.17           N
ATOM   1845  CA  GLY B   5     -15.029  -5.426   0.135  1.00 39.12           C
ATOM   1846  C   GLY B   5     -14.637  -5.835  -1.267  1.00 37.64           C
ATOM   1847  O   GLY B   5     -14.778  -7.006  -1.639  1.00 45.69           O
ATOM   1848  N   LYS B   6     -14.153  -4.897  -2.071  1.00 34.00           N
ATOM   1849  CA  LYS B   6     -13.638  -5.295  -3.377  1.00 32.13           C
ATOM   1850  C   LYS B   6     -12.151  -5.184  -3.592  1.00 28.44           C
ATOM   1851  O   LYS B   6     -11.483  -4.302  -3.055  1.00 25.94           O
ATOM   1852  CB  LYS B   6     -14.372  -4.554  -4.473  1.00 33.55           C
ATOM   1853  CG  LYS B   6     -15.874  -4.652  -4.285  1.00 37.63           C
ATOM   1854  CD  LYS B   6     -16.616  -4.577  -5.591  1.00 39.25           C
ATOM   1855  CE  LYS B   6     -18.085  -4.933  -5.433  1.00 38.87           C
ATOM   1856  NZ  LYS B   6     -18.911  -3.745  -5.162  1.00 39.20           N
ATOM   1857  N   VAL B   7     -11.657  -6.109  -4.407  1.00 26.29           N
ATOM   1858  CA  VAL B   7     -10.281  -6.166  -4.745  1.00 26.12           C
ATOM   1859  C   VAL B   7     -10.144  -5.498  -6.104  1.00 26.13           C
ATOM   1860  O   VAL B   7     -10.774  -5.928  -7.072  1.00 25.12           O
ATOM   1861  CB  VAL B   7      -9.771  -7.633  -4.787  1.00 28.15           C
ATOM   1862  CG1 VAL B   7      -8.322  -7.707  -5.258  1.00 28.49           C
ATOM   1863  CG2 VAL B   7      -9.933  -8.301  -3.425  1.00 28.37           C
END
'''

pdb_str_3 = '''\
CRYST1   53.930  108.200  147.160  90.00  90.00  90.00 P 21 21 21   16
SCALE1      0.018543  0.000000  0.000000        0.00000
SCALE2      0.000000  0.009242  0.000000        0.00000
SCALE3      0.000000  0.000000  0.006795        0.00000
ATOM      1  N   GLN A  -1     -11.627  10.458 -69.449  1.00 66.71           N
ATOM      2  CA  GLN A  -1     -12.125  11.752 -68.882  1.00 71.10           C
ATOM      3  C   GLN A  -1     -11.010  12.528 -68.080  1.00 70.59           C
ATOM      4  O   GLN A  -1      -9.956  12.796 -68.658  1.00 72.03           O
ATOM      5  CB  GLN A  -1     -13.372  11.454 -68.047  1.00 72.25           C
ATOM      6  CG  GLN A  -1     -14.331  12.612 -67.869  1.00 72.60           C
ATOM      7  CD  GLN A  -1     -15.171  12.445 -66.641  1.00 68.34           C
ATOM      8  OE1 GLN A  -1     -14.698  12.648 -65.520  1.00 64.02           O
ATOM      9  NE2 GLN A  -1     -16.418  12.080 -66.838  1.00 67.72           N
ATOM     24  N   SER A   2      -9.665   9.382 -65.470  1.00 40.80           N
ATOM     25  CA  SER A   2      -8.906   8.392 -66.236  1.00 37.26           C
ATOM     26  C   SER A   2      -7.689   7.905 -65.466  1.00 34.89           C
ATOM     27  O   SER A   2      -6.972   8.684 -64.878  1.00 37.04           O
ATOM     28  CB  SER A   2      -8.477   8.972 -67.586  1.00 35.56           C
ATOM     29  OG  SER A   2      -7.350   8.293 -68.131  1.00 34.81           O
ATOM    498  N   LEU A  73       4.762  -9.035 -55.093  1.00 54.92           N
ATOM    499  CA  LEU A  73       5.149  -7.682 -55.481  1.00 56.67           C
ATOM    500  C   LEU A  73       6.330  -7.693 -56.456  1.00 56.93           C
ATOM    501  O   LEU A  73       6.296  -7.006 -57.477  1.00 53.60           O
ATOM    502  CB  LEU A  73       5.425  -6.806 -54.257  1.00 58.11           C
ATOM    503  CG  LEU A  73       4.139  -6.316 -53.586  1.00 63.07           C
ATOM    504  CD1 LEU A  73       4.494  -5.448 -52.389  1.00 68.99           C
ATOM    505  CD2 LEU A  73       3.221  -5.562 -54.542  1.00 62.63           C
ATOM   1189  N   GLU A 171       0.727   9.279 -41.839  1.00 33.10           N
ATOM   1190  C   GLU A 171       2.659   9.149 -43.370  1.00 32.44           C
ATOM   1191  O   GLU A 171       3.759   9.624 -43.525  1.00 31.72           O
ATOM   1192  CA  GLU A 171       2.123   8.890 -41.964  0.50 33.32           C
ATOM   1193  CB  GLU A 171       2.301   7.408 -41.595  0.50 34.92           C
ATOM   1194  CG  GLU A 171       3.770   6.950 -41.526  0.50 37.34           C
ATOM   1195  CD  GLU A 171       3.984   5.453 -41.479  0.50 37.95           C
ATOM   1196  OE1 GLU A 171       3.024   4.749 -41.141  0.50 39.46           O
ATOM   1197  OE2 GLU A 171       5.118   4.979 -41.745  0.50 38.28           O
ATOM   1295  N   PHE A 186     -25.235  -3.045 -46.568  1.00 25.89           N
ATOM   1296  CA  PHE A 186     -25.800  -4.292 -47.005  1.00 26.70           C
ATOM   1297  C   PHE A 186     -25.763  -4.404 -48.510  1.00 25.41           C
ATOM   1298  O   PHE A 186     -26.508  -3.741 -49.214  1.00 23.13           O
ATOM   1299  CB  PHE A 186     -27.228  -4.359 -46.511  1.00 28.93           C
ATOM   1300  CG  PHE A 186     -27.843  -5.732 -46.569  1.00 33.84           C
ATOM   1301  CD1 PHE A 186     -27.121  -6.860 -46.230  1.00 37.15           C
ATOM   1302  CD2 PHE A 186     -29.175  -5.885 -46.906  1.00 38.18           C
ATOM   1303  CE1 PHE A 186     -27.698  -8.116 -46.239  1.00 39.95           C
ATOM   1304  CE2 PHE A 186     -29.770  -7.141 -46.913  1.00 42.87           C
ATOM   1305  CZ  PHE A 186     -29.026  -8.255 -46.588  1.00 42.49           C
ATOM   1528  N   GLU A 216     -21.253  -2.962 -57.909  1.00 23.63           N
ATOM   1529  C   GLU A 216     -20.899  -0.498 -58.040  1.00 23.99           C
ATOM   1530  O   GLU A 216     -20.144   0.434 -58.386  1.00 22.97           O
ATOM   1531  CA  GLU A 216     -20.839  -1.836 -58.751  0.40 23.33           C
ATOM   1532  CB  GLU A 216     -21.728  -1.708 -59.978  0.40 24.06           C
ATOM   1533  CG  GLU A 216     -21.277  -0.602 -60.946  0.40 23.98           C
ATOM   1534  CD  GLU A 216     -22.119  -0.501 -62.201  0.40 23.91           C
ATOM   1535  OE1 GLU A 216     -23.358  -0.345 -62.095  0.40 23.91           O
ATOM   1536  OE2 GLU A 216     -21.522  -0.553 -63.291  0.40 23.68           O
TER
ATOM   1821  N   SER B   2     -19.435   3.129  -2.367  1.00 68.33           N
ATOM   1822  CA  SER B   2     -18.122   2.823  -1.841  1.00 56.23           C
ATOM   1823  C   SER B   2     -17.667   1.496  -2.480  1.00 49.42           C
ATOM   1824  O   SER B   2     -18.495   0.790  -3.086  1.00 47.15           O
ATOM   1825  CB  SER B   2     -18.223   2.683  -0.315  1.00 61.95           C
ATOM   1826  OG  SER B   2     -18.227   1.337   0.111  1.00 58.21           O
ATOM   2295  N   LEU B  73      -1.962 -13.414 -11.627  1.00 39.07           N
ATOM   2296  C   LEU B  73      -3.519 -14.714 -10.237  1.00 43.61           C
ATOM   2297  O   LEU B  73      -4.239 -14.550  -9.256  1.00 44.31           O
ATOM   2298  CA  LEU B  73      -3.364 -13.601 -11.270  0.50 42.12           C
ATOM   2299  CB  LEU B  73      -4.240 -13.842 -12.506  0.50 41.67           C
ATOM   2300  CG  LEU B  73      -4.567 -12.548 -13.254  0.50 41.29           C
ATOM   2301  CD1 LEU B  73      -5.465 -12.886 -14.437  0.50 42.75           C
ATOM   2302  CD2 LEU B  73      -5.234 -11.495 -12.371  0.50 40.84           C
ATOM   2960  N   GLU B 171     -19.507  -8.000 -25.568  1.00 33.23           N
ATOM   2961  CA  GLU B 171     -19.305  -9.408 -25.319  1.00 35.30           C
ATOM   2962  C   GLU B 171     -19.653  -9.810 -23.897  1.00 31.64           C
ATOM   2963  O   GLU B 171     -20.228 -10.853 -23.684  1.00 36.45           O
ATOM   2964  CB  GLU B 171     -17.847  -9.790 -25.646  1.00 39.36           C
ATOM   2965  CG  GLU B 171     -17.549 -11.290 -25.624  1.00 42.46           C
ATOM   2966  CD  GLU B 171     -16.057 -11.580 -25.682  1.00 46.99           C
ATOM   2967  OE1 GLU B 171     -15.316 -10.762 -26.273  1.00 43.99           O
ATOM   2968  OE2 GLU B 171     -15.619 -12.618 -25.115  1.00 57.62           O
ATOM   3060  N   PHE B 186      -4.138  16.584 -21.892  1.00 30.38           N
ATOM   3061  C   PHE B 186      -2.849  17.033 -19.855  1.00 29.34           C
ATOM   3062  O   PHE B 186      -3.439  17.908 -19.219  1.00 28.85           O
ATOM   3063  CA  PHE B 186      -2.844  16.990 -21.351  0.50 30.79           C
ATOM   3064  CB  PHE B 186      -2.438  18.324 -21.877  0.50 34.30           C
ATOM   3065  CG  PHE B 186      -2.165  18.302 -23.332  0.50 36.69           C
ATOM   3066  CD1 PHE B 186      -0.909  17.954 -23.766  0.50 38.63           C
ATOM   3067  CD2 PHE B 186      -3.157  18.554 -24.258  0.50 38.07           C
ATOM   3068  CE1 PHE B 186      -0.618  17.909 -25.095  0.50 40.04           C
ATOM   3069  CE2 PHE B 186      -2.861  18.507 -25.599  0.50 40.28           C
ATOM   3070  CZ  PHE B 186      -1.595  18.190 -26.021  0.50 40.65           C
ATOM   3301  N   GLU B 216      -5.083  13.281 -10.276  1.00 26.06           N
ATOM   3302  C   GLU B 216      -7.573  13.217 -10.183  1.00 24.45           C
ATOM   3303  O   GLU B 216      -8.578  12.594  -9.847  1.00 22.48           O
ATOM   3304  CA  GLU B 216      -6.250  13.054  -9.419  0.50 24.68           C
ATOM   3305  CB  GLU B 216      -6.214  14.036  -8.239  0.50 24.66           C
ATOM   3306  CG  GLU B 216      -5.125  13.775  -7.196  0.50 25.28           C
ATOM   3307  CD  GLU B 216      -3.695  14.109  -7.630  0.50 25.48           C
ATOM   3308  OE1 GLU B 216      -3.457  14.775  -8.675  0.50 22.77           O
ATOM   3309  OE2 GLU B 216      -2.793  13.688  -6.871  0.50 27.08           O
TER
ATOM   3656  N   SER C   2     -21.208  16.849  -0.522  1.00 46.87           N
ATOM   3657  CA  SER C   2     -20.000  17.531  -0.933  1.00 49.27           C
ATOM   3658  C   SER C   2     -20.315  18.781  -1.786  1.00 48.14           C
ATOM   3659  O   SER C   2     -19.417  19.458  -2.317  1.00 47.91           O
ATOM   3660  CB  SER C   2     -19.248  16.576  -1.820  1.00 49.15           C
ATOM   3661  OG  SER C   2     -20.044  16.287  -2.951  1.00 48.27           O
ATOM   4130  N   LEU C  73     -35.100  33.259 -13.297  1.00 30.37           N
ATOM   4131  C   LEU C  73     -33.669  34.512 -11.749  1.00 31.67           C
ATOM   4132  O   LEU C  73     -33.066  34.339 -10.704  1.00 30.92           O
ATOM   4133  CA  LEU C  73     -33.737  33.424 -12.784  0.50 32.11           C
ATOM   4134  CB  LEU C  73     -32.737  33.755 -13.889  0.50 32.72           C
ATOM   4135  CG  LEU C  73     -31.281  33.568 -13.444  0.50 32.58           C
ATOM   4136  CD1 LEU C  73     -31.123  32.091 -13.133  0.50 33.06           C
ATOM   4137  CD2 LEU C  73     -30.263  34.044 -14.475  0.50 31.36           C
ATOM   4790  N   GLU C 171     -16.383  27.806 -25.479  1.00 26.63           N
ATOM   4791  CA  GLU C 171     -16.634  29.227 -25.280  1.00 27.35           C
ATOM   4792  C   GLU C 171     -16.472  29.632 -23.830  1.00 26.45           C
ATOM   4793  O   GLU C 171     -15.858  30.669 -23.552  1.00 27.99           O
ATOM   4794  CB  GLU C 171     -18.033  29.579 -25.749  1.00 29.87           C
ATOM   4795  CG  GLU C 171     -18.392  31.049 -25.714  1.00 32.88           C
ATOM   4796  CD  GLU C 171     -19.897  31.278 -25.925  1.00 33.29           C
ATOM   4797  OE1 GLU C 171     -20.491  30.491 -26.671  1.00 36.00           O
ATOM   4798  OE2 GLU C 171     -20.463  32.254 -25.419  1.00 29.64           O
ATOM   4890  N   PHE C 186     -32.084   3.244 -23.359  1.00 28.88           N
ATOM   4891  CA  PHE C 186     -33.431   2.830 -22.971  1.00 33.70           C
ATOM   4892  C   PHE C 186     -33.595   2.761 -21.478  1.00 33.49           C
ATOM   4893  O   PHE C 186     -33.026   1.892 -20.818  1.00 32.50           O
ATOM   4894  CB  PHE C 186     -33.768   1.469 -23.537  1.00 38.02           C
ATOM   4895  CG  PHE C 186     -33.888   1.497 -24.994  1.00 40.36           C
ATOM   4896  CD1 PHE C 186     -35.038   1.954 -25.573  1.00 43.21           C
ATOM   4897  CD2 PHE C 186     -32.804   1.161 -25.780  1.00 47.75           C
ATOM   4898  CE1 PHE C 186     -35.164   2.004 -26.943  1.00 46.94           C
ATOM   4899  CE2 PHE C 186     -32.890   1.230 -27.155  1.00 53.40           C
ATOM   4900  CZ  PHE C 186     -34.085   1.642 -27.741  1.00 53.46           C
ATOM   5123  N   GLU C 216     -32.293   6.544 -11.770  1.00 29.23           N
ATOM   5124  CA  GLU C 216     -31.230   6.774 -10.803  1.00 29.17           C
ATOM   5125  C   GLU C 216     -29.843   6.616 -11.430  1.00 28.37           C
ATOM   5126  O   GLU C 216     -28.876   7.239 -10.988  1.00 25.66           O
ATOM   5127  CB  GLU C 216     -31.385   5.803  -9.618  1.00 30.29           C
ATOM   5128  CG  GLU C 216     -32.587   6.059  -8.702  1.00 35.54           C
ATOM   5129  CD  GLU C 216     -33.973   5.699  -9.284  1.00 38.79           C
ATOM   5130  OE1 GLU C 216     -34.110   5.042 -10.354  1.00 35.53           O
ATOM   5131  OE2 GLU C 216     -34.965   6.108  -8.644  1.00 47.25           O
TER
ATOM   5402  N   SER D   2     -22.275  11.190 -66.206  1.00 42.42           N
ATOM   5403  CA  SER D   2     -22.944  12.197 -67.038  1.00 43.54           C
ATOM   5404  C   SER D   2     -24.228  12.691 -66.401  1.00 43.32           C
ATOM   5405  O   SER D   2     -25.028  11.902 -65.893  1.00 39.72           O
ATOM   5406  CB  SER D   2     -23.245  11.636 -68.426  1.00 44.54           C
ATOM   5407  OG  SER D   2     -24.332  12.329 -69.043  1.00 41.22           O
ATOM   5876  N   LEU D  73     -37.667  29.195 -57.139  1.00 46.87           N
ATOM   5877  CA  LEU D  73     -38.031  27.855 -57.596  1.00 53.31           C
ATOM   5878  C   LEU D  73     -39.085  27.881 -58.710  1.00 58.14           C
ATOM   5879  O   LEU D  73     -38.933  27.212 -59.728  1.00 54.43           O
ATOM   5880  CB  LEU D  73     -38.453  26.947 -56.428  1.00 51.85           C
ATOM   5881  CG  LEU D  73     -37.263  26.459 -55.589  1.00 57.14           C
ATOM   5882  CD1 LEU D  73     -37.740  25.563 -54.450  1.00 56.76           C
ATOM   5883  CD2 LEU D  73     -36.203  25.763 -56.430  1.00 54.73           C
ATOM   6589  N   GLU D 171     -34.939  10.760 -43.719  1.00 36.58           N
ATOM   6590  CA  GLU D 171     -36.316  11.130 -44.031  1.00 38.45           C
ATOM   6591  C   GLU D 171     -36.674  10.889 -45.492  1.00 34.39           C
ATOM   6592  O   GLU D 171     -37.743  10.403 -45.796  1.00 36.15           O
ATOM   6593  CB  GLU D 171     -36.636  12.593 -43.597  1.00 46.32           C
ATOM   6594  CG  GLU D 171     -37.237  12.647 -42.205  1.00 51.40           C
ATOM   6595  CD  GLU D 171     -38.750  12.455 -42.174  1.00 56.11           C
ATOM   6596  OE1 GLU D 171     -39.470  13.018 -43.048  1.00 55.38           O
ATOM   6597  OE2 GLU D 171     -39.225  11.781 -41.234  1.00 57.72           O
ATOM   6685  N   GLY D 185     -11.457  21.336 -45.861  1.00 26.90           N
ATOM   6686  CA  GLY D 185     -10.452  21.743 -44.885  1.00 28.00           C
ATOM   6687  C   GLY D 185      -9.985  23.139 -45.235  1.00 28.89           C
ATOM   6688  O   GLY D 185     -10.812  24.017 -45.460  1.00 30.81           O
ATOM   6689  N   PHE D 186      -8.671  23.306 -45.345  1.00 28.34           N
ATOM   6690  C   PHE D 186      -7.994  24.656 -47.258  1.00 29.13           C
ATOM   6691  O   PHE D 186      -7.170  24.013 -47.905  1.00 28.20           O
ATOM   6692  CA  PHE D 186      -8.041  24.559 -45.733  0.50 28.01           C
ATOM   6693  CB  PHE D 186      -6.620  24.597 -45.114  0.50 28.01           C
ATOM   6694  CG  PHE D 186      -5.978  25.966 -45.078  0.50 28.93           C
ATOM   6695  CD1 PHE D 186      -6.708  27.075 -44.715  0.50 28.34           C
ATOM   6696  CD2 PHE D 186      -4.620  26.126 -45.367  0.50 29.31           C
ATOM   6697  CE1 PHE D 186      -6.132  28.317 -44.669  0.50 28.18           C
ATOM   6698  CE2 PHE D 186      -4.030  27.371 -45.319  0.50 29.22           C
ATOM   6699  CZ  PHE D 186      -4.794  28.470 -44.974  0.50 29.31           C
ATOM   6930  N   GLU D 216     -11.467  23.461 -57.189  1.00 35.15           N
ATOM   6931  CA  GLU D 216     -11.751  22.379 -58.122  1.00 35.59           C
ATOM   6932  C   GLU D 216     -11.772  20.999 -57.440  1.00 31.90           C
ATOM   6933  O   GLU D 216     -12.469  20.088 -57.895  1.00 31.94           O
ATOM   6934  CB  GLU D 216     -10.688  22.369 -59.250  1.00 38.37           C
ATOM   6935  CG  GLU D 216     -10.759  23.561 -60.212  1.00 42.56           C
ATOM   6936  CD  GLU D 216     -10.288  24.920 -59.644  1.00 46.88           C
ATOM   6937  OE1 GLU D 216      -9.676  25.012 -58.542  1.00 45.82           O
ATOM   6938  OE2 GLU D 216     -10.525  25.929 -60.342  1.00 49.35           O
END
'''

test_pdb_4 = '''\
CRYST1  129.340  136.710  168.120  90.00  90.00  90.00 I 2 2 2      32
SCALE1      0.007732  0.000000  0.000000        0.00000
SCALE2      0.000000  0.007315  0.000000        0.00000
SCALE3      0.000000  0.000000  0.005948        0.00000
ATOM      1  N   GLY A   0      20.178   9.597  30.380  1.00 32.39           N
ATOM      2  CA  GLY A   0      19.786   9.406  31.801  1.00 38.52           C
ATOM      3  C   GLY A   0      20.415  10.525  32.584  1.00 38.95           C
ATOM      4  O   GLY A   0      20.229  10.681  33.787  1.00 44.05           O
ATOM    345  N  AASP A  43      43.208  62.734  38.734  0.50 21.83           N
ATOM    346  N  BASP A  43      43.214  62.729  38.785  0.50 22.10           N
ATOM    347  CA AASP A  43      44.141  61.954  37.935  0.50 23.37           C
ATOM    348  CA BASP A  43      44.182  62.012  37.961  0.50 23.35           C
ATOM    349  C  AASP A  43      44.271  60.570  38.560  0.50 22.51           C
ATOM    350  C  BASP A  43      44.528  60.667  38.595  0.50 23.04           C
ATOM    351  O  AASP A  43      44.292  59.538  37.875  0.50 20.85           O
ATOM    352  O  BASP A  43      45.032  59.760  37.945  0.50 22.31           O
ATOM    353  CB AASP A  43      45.503  62.651  37.781  0.50 22.81           C
ATOM    354  CB BASP A  43      45.450  62.838  37.738  0.50 23.04           C
ATOM    355  CG AASP A  43      46.251  62.823  39.095  0.50 27.19           C
ATOM    356  CG BASP A  43      45.207  64.072  36.867  0.50 23.04           C
ATOM    357  OD1AASP A  43      45.629  63.099  40.150  0.50 31.88           O
ATOM    358  OD1BASP A  43      44.045  64.332  36.464  0.50 31.20           O
ATOM    359  OD2AASP A  43      47.496  62.695  39.064  0.50 36.29           O
ATOM    360  OD2BASP A  43      46.189  64.797  36.603  0.50 28.48           O
TER    8500      GLU D 264
HETATM 8501 NA    NA A 265      38.174  35.784  51.193  1.00 18.09          NA
HETATM 8502 NA    NA A 266      37.726  40.289  60.423  1.00 22.78          NA
HETATM 8503 NA    NA A 267      13.212  41.425  36.424  1.00 30.14          NA
HETATM 8504 CL    CL A 268      44.904  37.770  59.641  1.00 29.60          CL
'''

test_pdb_5 = '''\
CRYST1  129.340  136.710  168.120  90.00  90.00  90.00 I 2 2 2      32
SCALE1      0.007732  0.000000  0.000000        0.00000
SCALE2      0.000000  0.007315  0.000000        0.00000
SCALE3      0.000000  0.000000  0.005948        0.00000
ATOM      1  N   GLY A   0      20.178   9.597  30.380  1.00 32.39           N
ATOM      2  CA  GLY A   0      19.786   9.406  31.801  1.00 38.52           C
ATOM      3  C   GLY A   0      20.415  10.525  32.584  1.00 38.95           C
ATOM      4  O   GLY A   0      20.229  10.681  33.787  1.00 44.05           O
TER    8500      GLU D 264
ATOM    345  N  AASP A  43      43.208  62.734  38.734  0.50 21.83           N
ATOM    346  N  BASP A  43      43.214  62.729  38.785  0.50 22.10           N
ATOM    347  CA AASP A  43      44.141  61.954  37.935  0.50 23.37           C
ATOM    348  CA BASP A  43      44.182  62.012  37.961  0.50 23.35           C
ATOM    349  C  AASP A  43      44.271  60.570  38.560  0.50 22.51           C
ATOM    350  C  BASP A  43      44.528  60.667  38.595  0.50 23.04           C
ATOM    351  O  AASP A  43      44.292  59.538  37.875  0.50 20.85           O
ATOM    352  O  BASP A  43      45.032  59.760  37.945  0.50 22.31           O
ATOM    353  CB AASP A  43      45.503  62.651  37.781  0.50 22.81           C
ATOM    354  CB BASP A  43      45.450  62.838  37.738  0.50 23.04           C
ATOM    355  CG AASP A  43      46.251  62.823  39.095  0.50 27.19           C
ATOM    356  CG BASP A  43      45.207  64.072  36.867  0.50 23.04           C
ATOM    357  OD1AASP A  43      45.629  63.099  40.150  0.50 31.88           O
ATOM    358  OD1BASP A  43      44.045  64.332  36.464  0.50 31.20           O
ATOM    359  OD2AASP A  43      47.496  62.695  39.064  0.50 36.29           O
ATOM    360  OD2BASP A  43      46.189  64.797  36.603  0.50 28.48           O
HETATM 8501 NA    NA A 265      38.174  35.784  51.193  1.00 18.09          NA
HETATM 8502 NA    NA A 266      37.726  40.289  60.423  1.00 22.78          NA
HETATM 8503 NA    NA A 267      13.212  41.425  36.424  1.00 30.14          NA
HETATM 8504 CL    CL A 268      44.904  37.770  59.641  1.00 29.60          CL
'''

test_pdb_6 = '''\
CRYST1  203.106   83.279  178.234  90.00 106.67  90.00 C 1 2 1      12
SCALE1      0.004924  0.000000  0.001474        0.00000
SCALE2      0.000000  0.012008  0.000000        0.00000
SCALE3      0.000000  0.000000  0.005857        0.00000
ATOM      1  N   ASP H   5      91.286 -31.834  73.572  1.00 77.83           N
ATOM      2  CA  ASP H   5      90.511 -32.072  72.317  1.00 78.04           C
ATOM      3  C   ASP H   5      90.136 -30.762  71.617  1.00 77.70           C
ATOM      4  O   ASP H   5      89.553 -29.857  72.225  1.00 77.56           O
ATOM      9  N   THR H   6      91.286 -31.834  73.572  1.00 77.83           N
ATOM     10  CA  THR H   6      90.511 -32.072  72.317  1.00 78.04           C
TER
ATOM   2517  N   GLY I 501      91.286 -31.834  73.572  1.00 77.83           N
ATOM   2518  CA  GLY I 501      90.511 -32.072  72.317  1.00 78.04           C
ATOM   2519  C   GLY I 501      90.136 -30.762  71.617  1.00 77.70           C
ATOM   2520  O   GLY I 501      89.553 -29.857  72.225  1.00 77.56           O
TER
ATOM   3802  N   ASP J   5      92.487   3.543  81.144  1.00 70.91           N
ATOM   3803  CA  ASP J   5      93.100   3.556  79.781  1.00 70.52           C
ATOM   3804  C   ASP J   5      92.161   2.961  78.728  1.00 70.38           C
ATOM   3805  O   ASP J   5      91.661   1.839  78.880  1.00 69.56           O
ATOM   3810  N   THR J   6      92.487   3.543  81.144  1.00 70.91           N
ATOM   3811  CA  THR J   6      93.100   3.556  79.781  1.00 70.52           C
TER
ATOM   6318  N   GLY K 501      92.487   3.543  81.144  1.00 70.91           N
ATOM   6319  CA  GLY K 501      93.100   3.556  79.781  1.00 70.52           C
ATOM   6320  C   GLY K 501      92.161   2.961  78.728  1.00 70.38           C
ATOM   6321  O   GLY K 501      91.661   1.839  78.880  1.00 69.56           O
TER
ATOM   7603  N   ASP L   5      61.028 -14.273  81.262  1.00 69.80           N
ATOM   7604  CA  ASP L   5      60.761 -13.218  80.242  1.00 70.31           C
ATOM   7605  C   ASP L   5      61.755 -13.281  79.080  1.00 70.96           C
ATOM   7606  O   ASP L   5      62.973 -13.267  79.280  1.00 70.38           O
ATOM   7611  N   THR L   6      61.028 -14.273  81.262  1.00 69.80           N
ATOM   7612  CA  THR L   6      60.761 -13.218  80.242  1.00 70.31           C
TER
ATOM  10119  N   GLY M 501      61.028 -14.273  81.262  1.00 69.80           N
ATOM  10120  CA  GLY M 501      60.761 -13.218  80.242  1.00 70.31           C
ATOM  10121  C   GLY M 501      61.755 -13.281  79.080  1.00 70.96           C
ATOM  10122  O   GLY M 501      62.973 -13.267  79.280  1.00 70.38           O
TER
HETATM11404  C1  NDG H 640      91.286 -31.834  73.572  1.00 77.83           C
HETATM11405  C2  NDG H 640      91.286 -31.834  73.572  1.00 77.83           C
HETATM11449  C1  NDG J 643      92.487   3.543  81.144  1.00 70.91           C
HETATM11494  C1  NDG L 646      61.028 -14.273  81.262  1.00 69.80           C
HETATM11495  C2  NDG L 646      61.028 -14.273  81.262  1.00 69.80           C
END
'''

test_pdb_7 = '''\
ATOM      1  N   ASP A   5      91.286 -31.834  73.572  1.00 77.83           N
ATOM      2  CA  ASP A   5      90.511 -32.072  72.317  1.00 78.04           C
ATOM      3  C   ASP A   5      90.136 -30.762  71.617  1.00 77.70           C
ATOM      4  O   ASP A   5      89.553 -29.857  72.225  1.00 77.56           O
ATOM      9  N   THR A   6      91.286 -31.834  73.572  1.00 77.83           N
ATOM     10  CA  THR A   6      90.511 -32.072  72.317  1.00 78.04           C
TER
ATOM   2517  N   GLY B 501      91.286 -31.834  73.572  1.00 77.83           N
ATOM   2518  CA  GLY B 501      90.511 -32.072  72.317  1.00 78.04           C
ATOM   2519  C   GLY B 501      90.136 -30.762  71.617  1.00 77.70           C
ATOM   2520  O   GLY B 501      89.553 -29.857  72.225  1.00 77.56           O
TER
ATOM   3802  N   ASP D   5      92.487   3.543  81.144  1.00 70.91           N
ATOM   3803  CA  ASP D   5      93.100   3.556  79.781  1.00 70.52           C
ATOM   3804  C   ASP D   5      92.161   2.961  78.728  1.00 70.38           C
ATOM   3805  O   ASP D   5      91.661   1.839  78.880  1.00 69.56           O
ATOM   3810  N   THR D   6      92.487   3.543  81.144  1.00 70.91           N
ATOM   3811  CA  THR D   6      93.100   3.556  79.781  1.00 70.52           C
TER
ATOM   6318  N   GLY E 501      92.487   3.543  81.144  1.00 70.91           N
ATOM   6319  CA  GLY E 501      93.100   3.556  79.781  1.00 70.52           C
ATOM   6320  C   GLY E 501      92.161   2.961  78.728  1.00 70.38           C
ATOM   6321  O   GLY E 501      91.661   1.839  78.880  1.00 69.56           O
TER
ATOM   7603  N   ASP F   5      61.217 -13.228  81.321  1.00 69.80           N
ATOM   7604  CA  ASP F   5      60.975 -12.176  80.288  1.00 70.31           C
ATOM   7605  C   ASP F   5      61.953 -12.288  79.115  1.00 70.96           C
ATOM   7606  O   ASP F   5      63.174 -12.323  79.302  1.00 70.38           O
ATOM   7611  N   THR F   6      61.217 -13.228  81.321  1.00 69.80           N
ATOM   7612  CA  THR F   6      60.975 -12.176  80.288  1.00 70.31           C
TER
ATOM  10119  N   GLY G 501      61.217 -13.228  81.321  1.00 69.80           N
ATOM  10120  CA  GLY G 501      60.975 -12.176  80.288  1.00 70.31           C
ATOM  10121  C   GLY G 501      61.953 -12.288  79.115  1.00 70.96           C
ATOM  10122  O   GLY G 501      63.174 -12.323  79.302  1.00 70.38           O
'''

test_pdb_7_addition = '''\
TER
ATOM  10119  N   GLY X 501     -44.932 -80.765  78.851  1.00 69.80           N
ATOM  10120  CA  GLY X 501     -45.174 -79.714  77.818  1.00 70.31           C
ATOM  10121  C   GLY X 501     -44.196 -79.826  76.644  1.00 70.96           C
ATOM  10122  O   GLY X 501     -42.975 -79.861  76.832  1.00 70.38           O
'''

test_pdb_8 = '''\
CRYST1  115.986  115.986   44.151  90.00  90.00 120.00 P 6
ATOM      2  CA  ALA A   1      64.128 -67.160 -16.901  1.00  2.19           C
ATOM      7  CA  ALA A   2      63.524 -68.826 -13.574  1.00  2.00           C
ATOM     12  CA  LYS A   3      59.889 -68.274 -12.586  1.00  2.00           C
ATOM     21  CA  ASP A   4      58.199 -67.416  -9.309  1.00  2.22           C
ATOM     29  CA  VAL A   5      55.099 -69.499  -8.580  1.00  2.00           C
ATOM     36  CA  LYS A   6      52.462 -68.681  -5.994  1.00  2.00           C
ATOM     45  CA  PHE A   7      49.510 -70.712  -4.872  1.00  2.00           C
ATOM     56  CA  GLY A   8      46.007 -70.315  -3.583  1.00  2.00           C
ATOM    139  CA  VAL A   9      43.283 -55.029  -9.675  1.00  2.35           C
ATOM    146  CA  ASN A  10      42.676 -52.469  -6.956  1.00  2.36           C
ATOM    154  CA  VAL A  11      45.554 -50.168  -7.937  1.00  2.29           C
ATOM    161  CA  LEU A  12      44.244 -50.043 -11.502  1.00  2.00           C
ATOM    169  CA  ALA A  13      40.561 -49.956 -10.406  1.00  2.35           C
ATOM    174  CA  ASP A  14      41.229 -47.208  -7.865  1.00  3.03           C
ATOM    182  CA  ALA A  15      43.030 -45.028 -10.390  1.00  2.69           C
ATOM    187  CA  VAL A  16      40.208 -45.442 -12.860  1.00  3.02           C
ATOM    194  CA  LYS A  17      37.331 -45.055 -10.390  1.00  2.71           C
ATOM    203  CA  VAL A  18      38.172 -41.474  -9.488  1.00  2.13           C
ATOM    210  CA  THR A  19      36.904 -40.296 -12.812  1.00  2.00           C
ATOM    566  CA  VAL A  20      37.819 -55.002 -19.732  1.00  2.05           C
ATOM    573  CA  ALA A  21      38.383 -51.225 -19.995  1.00  2.13           C
ATOM    578  CA  SER A  22      40.244 -51.871 -23.171  1.00  2.00           C
ATOM    584  CA  LYS A  23      37.585 -54.070 -24.767  1.00  2.49           C
ATOM    593  CA  ALA A  24      34.690 -51.803 -23.726  1.00  2.00           C
ATOM    598  CA  ASN A  25      36.487 -49.018 -25.568  1.00  2.00           C
ATOM    606  CA  ASP A  26      36.910 -50.949 -28.863  1.00  3.88           C
ATOM    614  CA  ALA A  27      33.274 -51.972 -29.032  1.00  3.99           C
ATOM    619  CA  ALA A  28      31.290 -49.213 -27.336  1.00  4.20           C
ATOM    624  CA  GLY A  29      33.537 -46.101 -27.308  1.00  3.15           C
ATOM    628  CA  ASP A  30      32.467 -45.237 -23.735  1.00  2.42           C
ATOM    636  CA  GLY A  31      31.652 -46.933 -20.435  1.00  2.00           C
ATOM    640  CA  THR A  32      35.290 -48.038 -19.921  1.00  2.00           C
ATOM    647  CA  THR A  33      35.265 -46.726 -16.371  1.00  2.06           C
ATOM    654  CA  THR A  34      31.920 -48.406 -15.620  1.00  2.56           C
ATOM    661  CA  ALA A  35      33.225 -51.616 -17.116  1.00  2.54           C
ATOM    666  CA  THR A  36      36.291 -51.335 -14.856  1.00  2.42           C
ATOM    673  CA  VAL A  37      34.269 -50.635 -11.747  1.00  3.40           C
ATOM    680  CA  LEU A  38      31.961 -53.605 -12.305  1.00  3.15           C
ATOM    688  CA  ALA A  39      34.947 -55.852 -12.956  1.00  2.38           C
ATOM    693  CA  GLN A  40      36.596 -54.839  -9.731  1.00  2.11           C
ATOM    702  CA  ALA A  41      33.351 -55.531  -7.927  1.00  2.78           C
ATOM    707  CA  ILE A  42      32.667 -58.960  -9.253  1.00  2.08           C
ATOM   3858  CA  ALA B   1      81.555 -73.470 -17.124  1.00  2.07           C
ATOM   3863  CA  ALA B   2      82.510 -74.911 -13.782  1.00  2.00           C
ATOM   3868  CA  LYS B   3      79.826 -77.479 -12.785  1.00  2.00           C
ATOM   3877  CA  ASP B   4      78.102 -78.275  -9.522  1.00  2.22           C
ATOM   3885  CA  VAL B   5      77.738 -81.976  -8.764  1.00  2.03           C
ATOM   3892  CA  LYS B   6      75.488 -83.519  -6.170  1.00  2.00           C
ATOM   3901  CA  PHE B   7      75.245 -87.087  -5.013  1.00  2.00           C
ATOM   3912  CA  GLY B   8      72.723 -89.591  -3.785  1.00  2.00           C
ATOM   3916  CA  ASN B   9      70.099 -88.361  -1.408  1.00  2.52           C
ATOM   3924  CA  ASP B  10      70.962 -84.672  -1.863  1.00  3.75           C
ATOM   3932  CA  ALA B  11      70.712 -84.850  -5.605  1.00  2.57           C
ATOM   3937  CA  ARG B  12      67.366 -86.601  -5.201  1.00  2.33           C
ATOM   3948  CA  VAL B  13      65.695 -84.360  -2.631  1.00  2.00           C
ATOM   3955  CA  LYS B  14      66.474 -81.440  -4.956  1.00  2.05           C
ATOM   3964  CA  MET B  15      65.273 -83.306  -8.068  1.00  2.49           C
ATOM   3972  CA  LEU B  16      62.014 -84.040  -6.264  1.00  2.32           C
ATOM   3980  CA  ARG B  17      61.348 -80.533  -5.048  1.00  2.45           C
ATOM   3991  CA  GLY B  18      61.641 -79.632  -8.714  1.00  2.12           C
ATOM   3995  CA  VAL B  19      59.091 -82.211  -9.831  1.00  2.22           C
ATOM   4002  CA  ASN B  20      56.721 -81.085  -7.152  1.00  2.44           C
ATOM   4010  CA  VAL B  21      56.739 -77.364  -8.133  1.00  2.42           C
ATOM   4017  CA  LEU B  22      55.793 -78.315 -11.709  1.00  2.00           C
ATOM   4025  CA  ALA B  23      53.474 -81.149 -10.626  1.00  2.36           C
ATOM   4030  CA  ASP B  24      51.680 -78.946  -8.069  1.00  2.96           C
ATOM   4038  CA  ALA B  25      51.101 -76.134 -10.576  1.00  2.68           C
ATOM   4043  CA  VAL B  26      49.693 -78.571 -13.084  1.00  2.94           C
ATOM   4050  CA  LYS B  27      47.585 -80.601 -10.605  1.00  2.85           C
ATOM   4059  CA  VAL B  28      45.313 -77.718  -9.655  1.00  2.15           C
ATOM   4066  CA  THR B  29      43.588 -77.953 -13.009  1.00  2.15           C
ATOM   4422  CA  VAL B  30      55.643 -86.411 -19.883  1.00  2.06           C
ATOM   4429  CA  ALA B  31      53.071 -83.639 -20.174  1.00  2.02           C
ATOM   4434  CA  SER B  32      54.722 -82.603 -23.360  1.00  2.04           C
ATOM   4510  CA  THR B  33      46.796 -86.923 -15.766  1.00  2.45           C
ATOM   4517  CA  ALA B  34      50.129 -87.956 -17.268  1.00  2.48           C
ATOM   4522  CA  THR B  35      51.827 -85.331 -15.045  1.00  2.35           C
ATOM   4529  CA  VAL B  36      50.007 -86.467 -11.913  1.00  3.42           C
ATOM   4536  CA  LEU B  37      50.897 -90.094 -12.505  1.00  3.16           C
ATOM   4544  CA  ALA B  38      54.531 -89.280 -13.126  1.00  2.34           C
ATOM   4549  CA  GLN B  39      54.795 -87.269  -9.909  1.00  2.02           C
ATOM   4558  CA  ALA B  40      53.272 -90.229  -8.089  1.00  2.64           C
ATOM   4563  CA  ILE B  41      55.515 -92.921  -9.376  1.00  2.00           C
'''

test_pdb_9 = '''\
ATOM      1  CA  ASP A   1       0.000   0.000   0.000  1.00 78.04           C
TER
ATOM      1  CA  ASP B   1      10.000  10.000  10.000  1.00 78.04           C
TER
ATOM      1  CA  ASP C   1      11.000  11.000  11.000  1.00 78.04           C
TER
ATOM      1  CA  ASP E   1      20.000  20.000  20.000  1.00 78.04           C
TER
ATOM      1  CA  ASP F   1      21.000  21.000  21.000  1.00 78.04           C
TER
ATOM      1  CA  ASP D   1       1.000   1.000   1.000  1.00 78.04           C
TER
'''

pdb_str10 = """\n
CRYST1   25.287   40.217   65.471  90.00  90.00  90.00 P 21 21 21
SCALE1      0.039546  0.000000  0.000000        0.00000
SCALE2      0.000000  0.024865  0.000000        0.00000
SCALE3      0.000000  0.000000  0.015274        0.00000
ATOM      1  P   DA A   5      -8.062  -5.965 -15.755  1.00 42.17           P
ATOM      2  OP1 DA A   5      -8.426  -7.228 -16.405  1.00 50.61           O
ATOM      3  OP2 DA A   5      -8.689  -5.557 -14.457  1.00 51.75           O
ATOM      4  O5  DA A   5      -6.496  -5.961 -15.638  1.00 34.89           O
ATOM      5  C5  DA A   5      -5.791  -6.321 -16.790  1.00 30.71           C
ATOM      6  C4  DA A   5      -4.355  -5.917 -16.600  1.00 34.43           C
ATOM      7  O4  DA A   5      -4.303  -4.509 -16.239  1.00 33.96           O
ATOM      8  C3  DA A   5      -3.630  -6.687 -15.491  1.00 35.56           C
ATOM      9  O3  DA A   5      -2.407  -7.257 -16.020  1.00 33.08           O
ATOM     10  C2  DA A   5      -3.531  -5.654 -14.384  1.00 32.41           C
ATOM     11  C1  DA A   5      -3.435  -4.334 -15.130  1.00 28.44           C
ATOM     12  N9  DA A   5      -3.904  -3.143 -14.449  1.00 28.37           N
ATOM     13  C8  DA A   5      -5.187  -2.933 -14.022  1.00 27.53           C
ATOM     14  N7  DA A   5      -5.401  -1.724 -13.565  1.00 29.33           N
ATOM     15  C5  DA A   5      -4.187  -1.082 -13.747  1.00 23.78           C
ATOM     16  C6  DA A   5      -3.761   0.226 -13.474  1.00 25.22           C
ATOM     17  N6  DA A   5      -4.519   1.150 -12.896  1.00 25.69           N
ATOM     18  N1  DA A   5      -2.485   0.535 -13.749  1.00 24.39           N
ATOM     19  C2  DA A   5      -1.712  -0.389 -14.320  1.00 24.89           C
ATOM     20  N3  DA A   5      -2.001  -1.641 -14.653  1.00 28.33           N
ATOM     21  C4  DA A   5      -3.268  -1.935 -14.326  1.00 27.45           C
TER
ATOM     83  P   DA B  17     -10.220   1.260  -1.207  1.00 27.94           P
ATOM     84  OP1 DA B  17     -11.370   2.143  -0.856  1.00 34.83           O
ATOM     85  OP2 DA B  17     -10.221   0.599  -2.553  1.00 31.17           O
ATOM     86  O5  DA B  17      -8.842   2.020  -1.098  1.00 26.12           O
ATOM     87  C5  DA B  17      -8.558   2.683   0.094  1.00 25.41           C
ATOM     88  C4  DA B  17      -7.407   3.619  -0.107  1.00 26.38           C
ATOM     89  O4  DA B  17      -6.208   2.886  -0.440  1.00 24.41           O
ATOM     90  C3  DA B  17      -7.600   4.631  -1.214  1.00 27.57           C
ATOM     91  O3  DA B  17      -6.972   5.834  -0.764  1.00 29.89           O
ATOM     92  C2  DA B  17      -6.902   3.980  -2.406  1.00 26.29           C
ATOM     93  C1  DA B  17      -5.771   3.225  -1.781  1.00 23.13           C
ATOM     94  N9  DA B  17      -5.444   1.986  -2.460  1.00 22.66           N
ATOM     95  C8  DA B  17      -6.295   0.942  -2.750  1.00 23.38           C
ATOM     96  N7  DA B  17      -5.700  -0.094  -3.288  1.00 20.62           N
ATOM     97  C5  DA B  17      -4.344   0.242  -3.234  1.00 20.59           C
ATOM     98  C6  DA B  17      -3.178  -0.447  -3.603  1.00 17.89           C
ATOM     99  N6  DA B  17      -3.184  -1.685  -4.072  1.00 20.22           N
ATOM    100  N1  DA B  17      -1.995   0.205  -3.497  1.00 19.61           N
ATOM    101  C2  DA B  17      -1.992   1.465  -3.030  1.00 20.38           C
ATOM    102  N3  DA B  17      -3.021   2.207  -2.621  1.00 20.80           N
ATOM    103  C4  DA B  17      -4.182   1.540  -2.774  1.00 19.17           C
TER
ATOM    247  P   DA D  17     -10.220   1.260  -1.207  1.00 27.94           P
ATOM    248  OP1 DA D  17     -11.370   2.143  -0.856  1.00 34.83           O
ATOM    249  OP2 DA D  17     -10.221   0.599  -2.553  1.00 31.17           O
ATOM    250  O5  DA D  17      -8.842   2.020  -1.098  1.00 26.12           O
ATOM    251  C5  DA D  17      -8.558   2.683   0.094  1.00 25.41           C
ATOM    252  C4  DA D  17      -7.407   3.619  -0.107  1.00 26.38           C
ATOM    253  O4  DA D  17      -6.208   2.886  -0.440  1.00 24.41           O
ATOM    254  C3  DA D  17      -7.600   4.631  -1.214  1.00 27.57           C
ATOM    255  O3  DA D  17      -6.972   5.834  -0.764  1.00 29.89           O
ATOM    256  C2  DA D  17      -6.902   3.980  -2.406  1.00 26.29           C
ATOM    257  C1  DA D  17      -5.771   3.225  -1.781  1.00 23.13           C
ATOM    258  N9  DA D  17      -5.444   1.986  -2.460  1.00 22.66           N
ATOM    259  C8  DA D  17      -6.295   0.942  -2.750  1.00 23.38           C
ATOM    260  N7  DA D  17      -5.700  -0.094  -3.288  1.00 20.62           N
ATOM    261  C5  DA D  17      -4.344   0.242  -3.234  1.00 20.59           C
ATOM    262  C6  DA D  17      -3.178  -0.447  -3.603  1.00 17.89           C
ATOM    263  N6  DA D  17      -3.184  -1.685  -4.072  1.00 20.22           N
ATOM    264  N1  DA D  17      -1.995   0.205  -3.497  1.00 19.61           N
ATOM    265  C2  DA D  17      -1.992   1.465  -3.030  1.00 20.38           C
ATOM    266  N3  DA D  17      -3.021   2.207  -2.621  1.00 20.80           N
ATOM    267  C4  DA D  17      -4.182   1.540  -2.774  1.00 19.17           C
END
"""


def run_selected_tests():
  """  Run selected tests

  1) List in "tests" the names of the particular test you want to run
  2) Comment out unittest.main()
  3) Un-comment unittest.TextTestRunner().run(run_selected_tests())
  """
  tests = ['test_correct_transform_selection']
  suite = unittest.TestSuite(list(map(TestSimpleAlignment,tests)))
  return suite

if __name__=='__main__':
  # use for individual tests
  # unittest.TextTestRunner().run(run_selected_tests())

  # Use to run all tests
  unittest.main(verbosity=0)
