from __future__ import absolute_import, division, print_function

from scitbx.array_family import flex
from scitbx import matrix
import iotbx.ncs as ncs
import iotbx.pdb
from iotbx.ncs import ncs_group_master_phil
import iotbx.phil
from time import time
from libtbx.test_utils import approx_equal
from six.moves import range


test_pdb_str = '''\
ATOM      1  N   THR A   1       9.670  10.289  11.135  1.00 20.00           N
ATOM      2  CA  THR A   1       9.559   8.931  10.615  1.00 20.00           C
ATOM      3  C   THR A   1       9.634   7.903  11.739  1.00 20.00           C
ATOM      4  O   THR B   1      10.449   8.027  12.653  1.00 20.00           O
ATOM      5  CB  THR B   1      10.660   8.630   9.582  1.00 20.00           C
ATOM      6  OG1 THR A   1      10.560   9.552   8.490  1.00 20.00           O
ATOM      7  CG2 THR A   1      10.523   7.209   9.055  1.00 20.00
'''

pdb_answer_0 = """\
CRYST1   18.415   14.419   12.493  90.00  90.00  90.00 P 1
ATOM      1  N   THR A   1      11.782  12.419   4.645  1.00 10.00           N
ATOM      2  CA  THR A   1      11.671  11.061   4.125  1.00 10.00           C
ATOM      3  C   THR A   1      11.746  10.033   5.249  1.00 10.00           C
ATOM      4  O   THR A   1      12.561  10.157   6.163  1.00 10.00           O
ATOM      5  CB  THR A   1      12.772  10.760   3.092  1.00 10.00           C
ATOM      6  OG1 THR A   1      12.672  11.682   2.000  1.00 10.00           O
ATOM      7  CG2 THR A   1      12.635   9.340   2.565  1.00 10.00           C
TER
ATOM      1  N   THR B   1       6.768   9.093   9.237  1.00 10.00           N
ATOM      2  CA  THR B   1       7.284   8.654   7.945  1.00 10.00           C
ATOM      3  C   THR B   1       8.638   7.968   8.097  1.00 10.00           C
ATOM      4  O   THR B   1       9.495   8.426   8.852  1.00 10.00           O
ATOM      5  CB  THR B   1       7.423   9.832   6.963  1.00 10.00           C
ATOM      6  OG1 THR B   1       6.144  10.446   6.765  1.00 10.00           O
ATOM      7  CG2 THR B   1       7.962   9.350   5.625  1.00 10.00           C
TER
ATOM      1  N   THR C   1       9.093   2.000  10.493  1.00 10.00           N
ATOM      2  CA  THR C   1       8.879   2.702   9.233  1.00 10.00           C
ATOM      3  C   THR C   1      10.081   3.570   8.875  1.00 10.00           C
ATOM      4  O   THR C   1      10.652   4.241   9.734  1.00 10.00           O
ATOM      5  CB  THR C   1       7.618   3.584   9.284  1.00 10.00           C
ATOM      6  OG1 THR C   1       6.472   2.770   9.559  1.00 10.00           O
ATOM      7  CG2 THR C   1       7.417   4.305   7.960  1.00 10.00           C
END
"""

test_pdb_str_2 = '''\
ATOM     45  N   PHEAa   6     219.693 144.930 112.416  1.00 50.00           N
ATOM     46  CA  PHEAa   6     218.871 146.020 112.886  1.00 50.00           C
ATOM     47  C   PHEAa   6     217.413 145.628 112.926  1.00 50.00           C
ATOM     48  O   PHEAa   6     216.730 145.905 113.908  1.00 50.00           O
TER
ATOM   1244  N   ARGAb   6     303.367 160.705 103.239  1.00 50.00           N
ATOM   1245  CA  ARGAb   6     302.396 160.991 104.331  1.00 50.00           C
ATOM   1246  C   ARGAb   6     302.285 162.473 104.586  1.00 50.00           C
ATOM   1247  O   ARGAb   6     302.837 163.292 103.851  1.00 50.00           O
TER
ATOM   2754  N   PHEAc   6     242.472 151.067 115.352  1.00 50.00           N
ATOM   2755  CA  PHEAc   6     241.314 151.789 115.823  1.00 50.00           C
ATOM   2756  C   PHEAc   6     240.094 150.900 115.864  1.00 50.00           C
ATOM   2757  O   PHEAc   6     239.358 150.912 116.847  1.00 50.00           O
TER
ATOM   3953  N   ARGAd   6     314.882 195.854 106.123  1.00 50.00           N
ATOM   3954  CA  ARGAd   6     313.875 195.773 107.215  1.00 50.00           C
ATOM   3955  C   ARGAd   6     313.239 197.116 107.471  1.00 50.00           C
ATOM   3956  O   ARGAd   6     313.460 198.078 106.736  1.00 50.00           O
TER
ATOM   5463  N   PHEAe   6     261.525 165.024 118.275  1.00 50.00           N
ATOM   5464  CA  PHEAe   6     260.185 165.283 118.746  1.00 50.00           C
ATOM   5465  C   PHEAe   6     259.365 164.014 118.787  1.00 50.00           C
ATOM   5466  O   PHEAe   6     258.673 163.762 119.769  1.00 50.00           O
TER
ATOM   6662  N   ARGAf   6     313.035 232.818 109.051  1.00 50.00           N
ATOM   6663  CA  ARGAf   6     312.124 232.379 110.143  1.00 50.00           C
ATOM   6664  C   ARGAf   6     311.048 233.405 110.399  1.00 50.00           C
ATOM   6665  O   ARGAf   6     310.909 234.383 109.665  1.00 50.00           O
END
'''

phil_str = '''\
ncs_group {
  reference = chain 'Aa'
  selection = chain 'Ac'
  selection = chain 'Ae'
}

ncs_group {
  reference = chain 'Ab'
  selection = chain 'Ad'
  selection = chain 'Af'
}
'''

phil_str2 = '''\
ncs_group {
  reference = chain 'Aa' or chain 'Ab'
  selection = chain 'Ae' or chain 'Af'
}
'''

def test_transform_update():
  """ Test update of rotation and translation using selection """
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_answer_0)
  ncs_obj = ncs.input(hierarchy=pdb_inp.construct_hierarchy())
  pdb_inp = iotbx.pdb.input(lines=pdb_answer_0,source_info=None)
  nrgl = ncs_obj.get_ncs_restraints_group_list()
  asu_site_cart = pdb_inp.atoms().extract_xyz()
  # reference matrices
  r1 = nrgl[0].copies[0].r
  t1 = nrgl[0].copies[0].t
  r2 = nrgl[0].copies[1].r
  t2 = nrgl[0].copies[1].t
  # modify matrices in the ncs group list
  nrgl[0].copies[0].r = r1 + r2
  nrgl[0].copies[0].t = t1 + t2
  nrgl[0].copies[1].r = r1 + r2
  nrgl[0].copies[1].t = t1 + t2
  nrgl.recalculate_ncs_transforms(asu_site_cart)
  # Get the updated values
  r1_n = nrgl[0].copies[0].r
  t1_n = nrgl[0].copies[0].t
  r2_n = nrgl[0].copies[1].r
  t2_n = nrgl[0].copies[1].t
  #
  assert approx_equal(r1, r1_n, eps=0.001)
  assert approx_equal(t1, t1_n, eps=0.1)
  assert approx_equal(r2, r2_n, eps=0.001)
  assert approx_equal(t2, t2_n, eps=0.1)

def test_check_for_max_rmsd():
  """ Test that ncs_restraints_group_list test is working properly """
  phil_groups = ncs_group_master_phil.fetch(
      iotbx.phil.parse(phil_str)).extract()
  pdb_inp = iotbx.pdb.input(source_info=None, lines=test_pdb_str_2)
  ncs_obj_phil = ncs.input(
      hierarchy=pdb_inp.construct_hierarchy(),
      ncs_phil_groups=phil_groups.ncs_group)
  nrgl = ncs_obj_phil.get_ncs_restraints_group_list()
  pdb_inp = iotbx.pdb.input(lines=test_pdb_str_2,source_info=None)
  ph = pdb_inp.construct_hierarchy()
  # passing test
  assert nrgl.check_for_max_rmsd(ph.atoms().extract_xyz() ,chain_max_rmsd=1)
  # make sure test fails when it suppose to
  nrgl[0].copies[1].t = matrix.col([100, -89.7668, 5.8996])
  assert not nrgl.check_for_max_rmsd(ph.atoms().extract_xyz(),chain_max_rmsd=1)

def test_center_of_coordinates_shift():
  """
  test shifting translation to and from the center of coordinates of the
  master ncs copy
  """
  # print sys._getframe().f_code.co_name
  # c = commons()

  pdb_inp = iotbx.pdb.input(source_info=None, lines=test_pdb_str_2)
  ncs_obj_phil = ncs.input(
      hierarchy=pdb_inp.construct_hierarchy())
  ncs_restraints_group_list = ncs_obj_phil.get_ncs_restraints_group_list()

  # ncs_restraints_group_list = c.ncs_restraints_group_list
  xrs = pdb_inp.xray_structure_simple()
  shifts = ncs_restraints_group_list.get_ncs_groups_centers(
    sites_cart = xrs.sites_cart())

  xyz = pdb_inp.atoms().extract_xyz()
  center_of_coor = (flex.vec3_double([xyz.sum()]) * (1/xyz.size())).round(8)
  # test shifts
  t1 = shifts[0].round(8)
  t2 = shifts[1].round(8)
  d1 = flex.sqrt((center_of_coor-t1).dot()).min_max_mean().as_tuple()
  d2 = flex.sqrt((center_of_coor-t2).dot()).min_max_mean().as_tuple()
  assert (d1 == d2)

  # test shift to center
  new_nrg = ncs_restraints_group_list.shift_translation_to_center(shifts = shifts)
  expected = (22.63275, 5.54625, 2.9375)
  assert (new_nrg[0].copies[0].t.round(5)).elems == expected
  # back to original coordinates system
  old_nrg = new_nrg.shift_translation_back_to_place(shifts=shifts)
  expected = (old_nrg[0].copies[0].t.round(5)).elems
  result = (ncs_restraints_group_list[0].copies[0].t.round(5)).elems
  assert result == expected

def test_ncs_selection():
  """
  verify that extended_ncs_selection, which include the master ncs copy and
  the portion of the protein we want to refine.
  """
  pdb_inp = iotbx.pdb.input(source_info=None, lines=test_pdb_str_2)
  ncs_obj_phil = ncs.input(
      hierarchy=pdb_inp.construct_hierarchy())
  ncs_restraints_group_list = ncs_obj_phil.get_ncs_restraints_group_list()
  # ncs_restraints_group_list._show()
  refine_selection = flex.size_t(range(30))
  result = ncs_restraints_group_list.get_extended_ncs_selection(
      refine_selection=refine_selection)
  # print list(result)
  expected = [0, 1, 2, 3, 4, 5, 6, 7, 24, 25, 26, 27, 28, 29]
  assert list(result) == expected

def test_whole_group_iselection():
  """ selection of a complete NCS group """
  phil_groups = ncs_group_master_phil.fetch(
      iotbx.phil.parse(phil_str)).extract()
  pdb_inp = iotbx.pdb.input(source_info=None, lines=test_pdb_str_2)
  ncs_obj = ncs.input(hierarchy=pdb_inp.construct_hierarchy(),
      ncs_phil_groups=phil_groups.ncs_group)
  nrgl = ncs_obj.get_ncs_restraints_group_list()
  assert len(nrgl) == nrgl.get_n_groups() == 2
  isel = nrgl[1].whole_group_iselection()
  expected = [4, 5, 6, 7, 12, 13, 14, 15, 20, 21, 22, 23]
  assert list(isel) == expected
  isel = nrgl[0].whole_group_iselection()
  expected = [0, 1, 2, 3, 8, 9, 10, 11, 16, 17, 18, 19]
  assert list(isel) == expected

def test_selection():
  """
  test that a atom selection propagates correctly to ncs_restraints_group_list
  """
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_answer_0)
  ncs_obj_phil = ncs.input(
      hierarchy=pdb_inp.construct_hierarchy())
  nrg = ncs_obj_phil.get_ncs_restraints_group_list()

  m1 = list(nrg[0].master_iselection)
  c1 = list(nrg[0].copies[0].iselection)
  c2 = list(nrg[0].copies[1].iselection)

  assert len(m1) == len(c1)  #                                           renumbering
  assert m1 == [0,   1,  2,  3,  4,  5,  6] #   0,  1, X,  3, X,  5, X | 0, 1, 3
  assert c1 == [7,   8,  9, 10, 11, 12, 13] #   7,  8, 9,  X, X, 12, X | 4, 5, 7
  assert c2 == [14, 15, 16, 17, 18, 19, 20] #  14, 15, X, 17, X, 19, X | 8, 9, 11

  selection1 = flex.size_t([0,1,5,3,100,101])
  selection2 = flex.size_t([0,1,5,3,7,8,9,12,100,101])
  selection3 = flex.size_t([0,1,5,3,7,8,9,12,14,15,19,17,100,101])
  # gone iseqs for selection3: 2,4,6,10,11,13,16,18,20-99

  new_nrg = nrg.select(flex.bool(102, selection1))
  # only atoms in master are selected
  mt = list(new_nrg[0].master_iselection)
  c1t = list(new_nrg[0].copies[0].iselection)

  assert mt == []
  assert c1t == []

  # atoms selected in both master and copies
  new_nrg = nrg.select(flex.bool(102, selection2))
  # only atoms in master are selected
  mt = list(new_nrg[0].master_iselection)
  c1t = list(new_nrg[0].copies[0].iselection)

  assert mt == []
  assert c1t == []

  new_nrg = nrg.select(flex.bool(102, selection3))
  # only atoms in master are selected
  mt = list(new_nrg[0].master_iselection)
  c1t = list(new_nrg[0].copies[0].iselection)
  c2t = list(new_nrg[0].copies[1].iselection)

  assert mt == [0, 1, 3], list(mt)
  assert c1t == [4, 5, 7], list(c1t)
  assert c2t == [8, 9, 11], list(c2t)

def test_split_by_chain():
  phil_groups = ncs_group_master_phil.fetch(
      iotbx.phil.parse(phil_str2)).extract()
  pdb_inp = iotbx.pdb.input(source_info=None, lines=test_pdb_str_2)
  pars = ncs.input.get_default_params()
  pars.ncs_search.chain_max_rmsd = 100
  h = pdb_inp.construct_hierarchy()
  ncs_obj_phil = ncs.input(
      hierarchy=h,
      ncs_phil_groups=phil_groups.ncs_group,
      params=pars.ncs_search)
  nrgl = ncs_obj_phil.get_ncs_restraints_group_list()
  assert nrgl.get_array_of_str_selections() == \
      [["chain 'Aa' or chain 'Ab'", "chain 'Ae' or chain 'Af'"]]
  splitted_nrgl = nrgl.split_by_chains(hierarchy=h)
  assert splitted_nrgl.get_n_groups() == 2
  for g in splitted_nrgl:
    assert g.get_number_of_copies() == 1
    assert approx_equal(g.copies[0].r, nrgl[0].copies[0].r)
    assert approx_equal(g.copies[0].t, nrgl[0].copies[0].t)
  splitted_nrgl.update_str_selections_if_needed(hierarchy=h)
  assert splitted_nrgl.get_array_of_str_selections() == \
     [["chain 'Aa'", "chain 'Ae'"], ["chain 'Ab'", "chain 'Af'"]]

def test_split_by_chain_2():
  """ Test of split_by_chains where chains are splitted, e.g.
  protein and ligand are separated by TER
  """
  pdb_str = """\
CRYST1  174.866  238.091  243.466  90.00  90.00  90.00 P 21 21 21
SCALE1      0.005719  0.000000  0.000000        0.00000
SCALE2      0.000000  0.004200  0.000000        0.00000
SCALE3      0.000000  0.000000  0.004107        0.00000
ATOM  15295  N   SER D  74      59.475 -46.773 121.838  1.00383.47           N
ATOM  15296  CA  SER D  74      59.914 -46.196 120.574  1.00383.24           C
ATOM  15297  C   SER D  74      59.605 -44.703 120.512  1.00379.51           C
ATOM  15298  O   SER D  74      58.830 -44.258 119.666  1.00374.31           O
ATOM  15299  CB  SER D  74      59.247 -46.915 119.397  1.00380.53           C
ATOM  15300  OG  SER D  74      59.544 -48.301 119.400  1.00383.68           O
ATOM  15301  N   ASN D  75      60.241 -43.943 121.402  1.00382.45           N
ATOM  15302  CA  ASN D  75      60.129 -42.482 121.464  1.00380.05           C
ATOM  15303  C   ASN D  75      58.725 -41.896 121.293  1.00372.68           C
ATOM  15304  O   ASN D  75      58.581 -40.746 120.881  1.00370.19           O
ATOM  15305  CB  ASN D  75      61.046 -41.855 120.418  1.00382.20           C
ATOM  15306  CG  ASN D  75      62.508 -42.175 120.656  1.00389.96           C
ATOM  15307  OD1 ASN D  75      62.958 -42.273 121.796  1.00394.06           O
ATOM  15308  ND2 ASN D  75      63.263 -42.322 119.575  1.00392.78           N
ATOM  15309  N   GLY D  76      57.703 -42.688 121.603  1.00369.47           N
ATOM  15310  CA  GLY D  76      56.321 -42.306 121.367  1.00362.43           C
ATOM  15311  C   GLY D  76      55.844 -41.019 122.011  1.00358.49           C
ATOM  15312  O   GLY D  76      56.116 -40.749 123.182  1.00360.25           O
ATOM  15313  N   ARG D  77      55.124 -40.219 121.235  1.00353.03           N
ATOM  15314  CA  ARG D  77      54.452 -39.045 121.771  1.00347.98           C
ATOM  15315  C   ARG D  77      53.192 -38.773 120.956  1.00340.73           C
ATOM  15316  O   ARG D  77      53.132 -39.072 119.764  1.00339.65           O
ATOM  15317  CB  ARG D  77      55.383 -37.824 121.801  1.00349.75           C
ATOM  15318  CG  ARG D  77      55.931 -37.385 120.451  1.00350.35           C
ATOM  15319  CD  ARG D  77      57.293 -37.995 120.136  1.00358.77           C
ATOM  15320  NE  ARG D  77      58.029 -37.121 119.227  1.00359.99           N
ATOM  15321  CZ  ARG D  77      59.247 -37.355 118.749  1.00366.18           C
ATOM  15322  NH1 ARG D  77      59.921 -38.437 119.111  1.00371.93           N
ATOM  15323  NH2 ARG D  77      59.809 -36.473 117.933  1.00366.55           N
TER
ATOM  15324  N   TYR D  78      52.166 -38.264 121.625  1.00335.62           N
ATOM  15325  CA  TYR D  78      50.919 -37.914 120.961  1.00328.66           C
ATOM  15326  C   TYR D  78      50.974 -36.476 120.456  1.00325.47           C
ATOM  15327  O   TYR D  78      50.504 -35.558 121.125  1.00322.51           O
ATOM  15328  CB  TYR D  78      49.734 -38.100 121.913  1.00324.55           C
ATOM  15329  CG  TYR D  78      49.238 -39.526 122.016  1.00325.38           C
ATOM  15330  CD1 TYR D  78      49.026 -40.289 120.877  1.00324.97           C
ATOM  15331  CD2 TYR D  78      48.973 -40.106 123.251  1.00326.27           C
ATOM  15332  CE1 TYR D  78      48.567 -41.589 120.962  1.00325.32           C
ATOM  15333  CE2 TYR D  78      48.514 -41.408 123.345  1.00326.09           C
ATOM  15334  CZ  TYR D  78      48.314 -42.144 122.197  1.00325.64           C
ATOM  15335  OH  TYR D  78      47.859 -43.439 122.280  1.00325.77           O
ATOM  19064  N   ASN E  75      29.075   5.259 120.437  1.00232.69           N
ATOM  19065  CA  ASN E  75      28.096   5.658 121.456  1.00227.54           C
ATOM  19066  C   ASN E  75      27.246   4.544 122.068  1.00223.14           C
ATOM  19067  O   ASN E  75      26.146   4.804 122.555  1.00219.10           O
ATOM  19068  CB  ASN E  75      27.159   6.709 120.876  1.00226.97           C
ATOM  19069  CG  ASN E  75      27.882   7.970 120.467  1.00230.84           C
ATOM  19070  OD1 ASN E  75      28.845   8.385 121.111  1.00232.56           O
ATOM  19071  ND2 ASN E  75      27.417   8.595 119.393  1.00233.24           N
ATOM  19072  N   GLY E  76      27.759   3.318 122.053  1.00224.32           N
ATOM  19073  CA  GLY E  76      26.998   2.150 122.461  1.00220.68           C
ATOM  19074  C   GLY E  76      26.384   2.164 123.847  1.00215.47           C
ATOM  19075  O   GLY E  76      27.013   2.573 124.823  1.00215.60           O
ATOM  19076  N   ARG E  77      25.135   1.722 123.922  1.00211.38           N
ATOM  19077  CA  ARG E  77      24.469   1.517 125.200  1.00206.91           C
ATOM  19078  C   ARG E  77      23.459   0.384 125.076  1.00203.78           C
ATOM  19079  O   ARG E  77      22.895   0.146 124.008  1.00203.81           O
ATOM  19080  CB  ARG E  77      23.800   2.801 125.707  1.00204.73           C
ATOM  19081  CG  ARG E  77      22.666   3.330 124.849  1.00204.05           C
ATOM  19082  CD  ARG E  77      23.143   4.321 123.799  1.00208.96           C
ATOM  19083  NE  ARG E  77      22.118   5.331 123.560  1.00208.57           N
ATOM  19084  CZ  ARG E  77      22.230   6.335 122.698  1.00211.78           C
ATOM  19085  NH1 ARG E  77      23.334   6.480 121.983  1.00215.57           N
ATOM  19086  NH2 ARG E  77      21.241   7.205 122.563  1.00211.19           N
TER
ATOM  19087  N   TYR E  78      23.286  -0.346 126.170  1.00201.65           N
ATOM  19088  CA  TYR E  78      22.334  -1.441 126.223  1.00198.95           C
ATOM  19089  C   TYR E  78      20.973  -0.908 126.642  1.00195.03           C
ATOM  19090  O   TYR E  78      20.616  -0.967 127.815  1.00193.71           O
ATOM  19091  CB  TYR E  78      22.813  -2.517 127.201  1.00199.92           C
ATOM  19092  CG  TYR E  78      23.854  -3.449 126.623  1.00204.23           C
ATOM  19093  CD1 TYR E  78      23.681  -4.014 125.369  1.00205.06           C
ATOM  19094  CD2 TYR E  78      25.014  -3.757 127.326  1.00208.09           C
ATOM  19095  CE1 TYR E  78      24.627  -4.863 124.832  1.00209.52           C
ATOM  19096  CE2 TYR E  78      25.969  -4.606 126.795  1.00212.29           C
ATOM  19097  CZ  TYR E  78      25.769  -5.155 125.547  1.00213.03           C
ATOM  19098  OH  TYR E  78      26.711  -6.001 125.009  1.00218.33           O
TER
  """
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  h = pdb_inp.construct_hierarchy()
  pars = ncs.input.get_default_params()
  pars.ncs_search.chain_similarity_threshold = 0.5

  ncs_obj = ncs.input(hierarchy=h, params=pars.ncs_search)
  spec = ncs_obj.get_ncs_info_as_spec()
  # ncs_obj.show(format='spec')
  assert len(spec.ncs_groups()) ==1

def run_tests():
  test_transform_update()
  test_check_for_max_rmsd()
  test_center_of_coordinates_shift()
  test_ncs_selection()
  test_whole_group_iselection()
  test_selection()
  test_split_by_chain()
  test_split_by_chain_2()

if __name__=='__main__':
  t0 = time()
  run_tests()
  print("OK. Time: %8.3f"%(time()-t0))
