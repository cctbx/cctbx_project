from __future__ import division, print_function
import time
import mmtbx.model
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal
import iotbx.pdb
import cctbx.geometry_restraints.process_nonbonded_proxies as pnp


#-------------------------------------------------------------------------------

def obtain_model(raw_records):
  '''
  Helper function to obtain model class from raw records
  '''
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.pdb_interpretation.allow_polymer_cross_special_position=True
  pdb_inp = iotbx.pdb.input(lines=raw_records.split("\n"), source_info=None)
  model = mmtbx.model.manager(
    model_input = pdb_inp,
    pdb_interpretation_params = params,
    build_grm   = True,
    log         = null_out())
  return model

#-------------------------------------------------------------------------------

def test_inline_angle():
  '''
  Test cos_vec(u,v)
  '''
  u = (40.196,3.261,-48.474)
  v = (39.466,4.279,-52.202)
  w = (37.407,5.077,-51.025)
  result = pnp.cos_vec(u, v, w)
  expected = 0.5072
  msg = 'The difference is: {}'.format(result - expected)
  assert approx_equal(result, expected, eps=0.001), msg


#-------------------------------------------------------------------------------

def test_1_5_overlaps():
  '''
  Test if 1_5 overlaps are correctly identified
  '''
  model = obtain_model(raw_records_0)
  hd_sel = model.get_hd_selection()
  grm = model.get_restraints_manager().geometry
  full_connectivity_table = grm.shell_sym_tables[0].full_simple_connectivity()

  outstring = '1-5 Interaction test error. {}'
  # check that direction of function calling does not matter
  tst = pnp.check_if_1_5_interaction(21, 33,hd_sel,full_connectivity_table)
  msg = 'Test results depend on atoms order'
  assert(tst), msg
  tst = pnp.check_if_1_5_interaction(33, 21,hd_sel,full_connectivity_table)
  assert(tst), msg
  # check 1-4 interaction
  tst = pnp.check_if_1_5_interaction(33, 20,hd_sel,full_connectivity_table)
  msg = 'Test fails on 1-4 interaction'
  assert(not tst), msg
  # check 1-6 interaction
  tst = pnp.check_if_1_5_interaction(33, 38,hd_sel,full_connectivity_table)
  msg = 'Test fails on 1-6 interaction'
  assert(not tst), msg
  # 1-5 interaction of atoms other than hydrogen
  tst = pnp.check_if_1_5_interaction(38, 25,hd_sel,full_connectivity_table)
  msg = 'Test fails on 1-5 non hydrogen interaction'
  assert(not tst), msg
  # 1-5 interaction of two hydrogens
  tst = pnp.check_if_1_5_interaction(33, 31,hd_sel,full_connectivity_table)
  msg = 'Test fails on 1-5 two hydrogen interaction'
  assert(not tst), msg

#-------------------------------------------------------------------------------

raw_records_0 = """\
CRYST1   80.020   97.150   49.850  90.00  90.00  90.00 C 2 2 21
ATOM   1271  N   ILE A  83      31.347   4.310 -43.960  1.00  9.97           N
ATOM   1272  CA  ILE A  83      32.076   3.503 -44.918  1.00 19.49           C
ATOM   1273  C   ILE A  83      32.062   4.261 -46.243  1.00 15.26           C
ATOM   1274  O   ILE A  83      31.006   4.660 -46.740  1.00 16.12           O
ATOM   1275  CB  ILE A  83      31.486   2.080 -45.048  1.00 18.05           C
ATOM   1276  CG1 ILE A  83      32.328   1.263 -46.027  1.00 17.98           C
ATOM   1277  CG2 ILE A  83      30.020   2.134 -45.469  1.00 28.62           C
ATOM   1278  CD1 ILE A  83      32.130  -0.226 -45.888  1.00 38.88           C
ATOM   1279  H   ILE A  83      30.735   4.798 -44.316  1.00  9.97           H
ATOM   1280  HA  ILE A  83      32.999   3.421 -44.630  1.00 19.49           H
ATOM   1281  HB  ILE A  83      31.534   1.655 -44.178  1.00 18.05           H
ATOM   1282 HG12 ILE A  83      32.087   1.512 -46.933  1.00 17.98           H
ATOM   1283 HG13 ILE A  83      33.266   1.454 -45.871  1.00 17.98           H
ATOM   1284 HG21 ILE A  83      29.627   1.256 -45.347  1.00 28.62           H
ATOM   1285 HG22 ILE A  83      29.554   2.783 -44.919  1.00 28.62           H
ATOM   1286 HG23 ILE A  83      29.967   2.392 -46.402  1.00 28.62           H
ATOM   1287 HD11 ILE A  83      33.003  -0.644 -45.858  1.00 38.88           H
ATOM   1288 HD12 ILE A  83      31.654  -0.404 -45.062  1.00 38.88           H
ATOM   1289 HD13 ILE A  83      31.608  -0.541 -46.640  1.00 38.88           H
ATOM   1290  N   ILE A  84      33.252   4.542 -46.758  1.00 11.36           N
ATOM   1291  CA  ILE A  84      33.381   5.290 -47.993  1.00  9.27           C
ATOM   1292  C   ILE A  84      33.906   4.417 -49.101  1.00 11.06           C
ATOM   1293  O   ILE A  84      34.881   3.694 -48.921  1.00 14.52           O
ATOM   1294  CB  ILE A  84      34.308   6.508 -47.813  1.00 11.89           C
ATOM   1295  CG1 ILE A  84      33.730   7.436 -46.746  1.00  9.03           C
ATOM   1296  CG2 ILE A  84      34.471   7.257 -49.132  1.00 10.52           C
ATOM   1297  CD1 ILE A  84      34.436   8.761 -46.632  1.00 12.96           C
ATOM   1298  H   ILE A  84      34.002   4.308 -46.408  1.00 11.36           H
ATOM   1299  HA  ILE A  84      32.515   5.625 -48.269  1.00  9.27           H
ATOM   1300  HB  ILE A  84      35.179   6.198 -47.520  1.00 11.89           H
ATOM   1301 HG12 ILE A  84      32.800   7.587 -46.942  1.00  9.03           H
ATOM   1302 HG13 ILE A  84      33.808   6.995 -45.886  1.00  9.03           H
ATOM   1303 HG21 ILE A  84      35.193   7.898 -49.045  1.00 10.52           H
ATOM   1304 HG22 ILE A  84      34.688   6.643 -49.848  1.00 10.52           H
ATOM   1305 HG23 ILE A  84      33.642   7.718 -49.335  1.00 10.52           H
ATOM   1306 HD11 ILE A  84      33.911   9.335 -46.052  1.00 12.96           H
ATOM   1307 HD12 ILE A  84      35.300   8.606 -46.229  1.00 12.96           H
ATOM   1308 HD13 ILE A  84      34.530   9.185 -47.497  1.00 12.96           H
ATOM   1309  N   THR A  85      33.228   4.470 -50.242  1.00 17.41           N
ATOM   1310  CA  THR A  85      33.615   3.692 -51.411  1.00 14.87           C
ATOM   1311  C   THR A  85      33.743   4.584 -52.633  1.00 15.68           C
ATOM   1312  O   THR A  85      33.232   5.702 -52.660  1.00 10.37           O
ATOM   1313  CB  THR A  85      32.581   2.585 -51.737  1.00 13.04           C
ATOM   1314  OG1 THR A  85      31.285   3.169 -51.914  1.00 19.63           O
ATOM   1315  CG2 THR A  85      32.524   1.546 -50.626  1.00 14.50           C
ATOM   1316  H   THR A  85      32.530   4.957 -50.365  1.00 17.41           H
ATOM   1317  HA  THR A  85      34.473   3.268 -51.259  1.00 14.87           H
ATOM   1318  HB  THR A  85      32.872   2.178 -52.553  1.00 13.04           H
ATOM   1319  HG1 THR A  85      31.049   3.564 -51.211  1.00 19.63           H
ATOM   1320 HG21 THR A  85      32.000   0.780 -50.908  1.00 14.50           H
ATOM   1321 HG22 THR A  85      33.420   1.251 -50.402  1.00 14.50           H
ATOM   1322 HG23 THR A  85      32.116   1.932 -49.835  1.00 14.50           H
"""


if (__name__ == "__main__"):
  t0 = time.time()
  test_1_5_overlaps()
  test_inline_angle()
  print("OK. Time: %8.3f"%(time.time()-t0))
