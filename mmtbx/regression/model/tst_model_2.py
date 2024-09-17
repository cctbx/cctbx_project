from __future__ import absolute_import, division, print_function
import mmtbx.model
import libtbx.load_env
import iotbx.pdb
import os
import tempfile
import time
from libtbx.utils import null_out
from cctbx import crystal, adptbx
from six.moves import range
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex

pdb_str_1 = """
CRYST1   79.110   79.110   37.465  90.00  90.00  90.00 P 43 21 2
ATOM      1  N   GLY A  67      11.351   9.426  29.699  1.00 16.57      A    N
ATOM      2  CA  GLY A  67      12.344   8.654  30.419  1.00 16.65      A    C
ATOM      3  C   GLY A  67      13.703   9.318  30.525  1.00 17.27      A    C
ATOM      4  O   GLY A  67      14.613   8.754  31.138  1.00 18.12      A    O
ATOM      5  HA2 GLY A  67      11.164   8.290  32.128  1.00 16.65      A    H
ATOM      6  HA3 GLY A  67      12.899   8.562  30.893  1.00 16.65      A    H
ATOM      7  N   ARG A  68      13.850  10.509  29.946  1.00 17.04      A    N
ATOM      8  CA  ARG A  68      15.122  11.232  30.010  1.00 18.21      A    C
ATOM      9  C   ARG A  68      14.955  12.712  30.337  1.00 19.28      A    C
ATOM     10  O   ARG A  68      15.856  13.506  30.060  1.00 18.39      A    O
ATOM     11  CB  ARG A  68      15.888  11.102  28.683  1.00 17.72      A    C
ATOM     12  CG  ARG A  68      15.297  11.878  27.495  1.00 17.98      A    C
ATOM     13  CD  ARG A  68      16.259  11.919  26.349  1.00 18.60      A    C
ATOM     14  NE  ARG A  68      15.794  12.838  25.306  1.00 17.24      A    N
ATOM     15  CZ  ARG A  68      16.143  12.775  24.023  1.00 19.10      A    C
ATOM     16  NH1 ARG A  68      16.951  11.815  23.586  1.00 18.15      A    N
ATOM     17  NH2 ARG A  68      15.667  13.670  23.169  1.00 19.17      A    N
ATOM     18  HA  ARG A  68      15.977  11.764  31.624  1.00 18.21      A    H
ATOM     19  HB2 ARG A  68      17.228  11.790  28.507  1.00 17.72      A    H
ATOM     20  HB3 ARG A  68      16.065   9.301  28.834  1.00 17.72      A    H
ATOM     21  HG2 ARG A  68      14.632  11.422  26.483  1.00 17.98      A    H
ATOM     22  HG3 ARG A  68      15.102  11.944  27.956  1.00 17.98      A    H
ATOM     23  HD2 ARG A  68      17.713  11.903  27.452  1.00 18.60      A    H
ATOM     24  HD3 ARG A  68      17.031  10.479  25.184  1.00 18.60      A    H
ATOM     25  H  AARG A  68      13.747  10.274  28.969  0.50 17.04      A    H
ATOM     26  HE AARG A  68      14.652  14.091  26.079  0.50 17.24      A    H
ATOM     27 HH11AARG A  68      16.935  11.771  24.831  0.50 18.15      A    H
ATOM     28 HH12AARG A  68      17.117  12.550  21.627  0.50 18.15      A    H
ATOM     29 HH21AARG A  68      14.846  14.096  22.818  0.50 19.17      A    H
ATOM     30 HH22AARG A  68      15.127  14.548  23.046  0.50 19.17      A    H
ATOM     31  D  BARG A  68      12.950  10.336  30.200  0.50 17.04      A    D
ATOM     32  DE BARG A  68      14.587  13.771  24.772  0.50 17.24      A    D
ATOM     33 DH11BARG A  68      16.347  11.268  24.478  0.50 18.15      A    D
ATOM     34 DH12BARG A  68      16.772  11.748  23.372  0.50 18.15      A    D
ATOM     35 DH21BARG A  68      14.983  15.028  23.897  0.50 19.17      A    D
ATOM     36 DH22BARG A  68      15.963  14.351  21.579  0.50 19.17      A    D
ATOM     37  N   THR A  69      13.811  13.085  30.907  1.00 18.60      A    N
ATOM     38  CA  THR A  69      13.578  14.468  31.314  1.00 21.27      A    C
ATOM     39  C   THR A  69      13.418  14.522  32.825  1.00 23.37      A    C
ATOM     40  O   THR A  69      12.315  14.330  33.336  1.00 22.36      A    O
ATOM     41  CB  THR A  69      12.332  15.076  30.647  1.00 18.85      A    C
ATOM     42  OG1 THR A  69      12.366  14.826  29.239  1.00 16.82      A    O
ATOM     43  CG2 THR A  69      12.281  16.578  30.889  1.00 22.56      A    C
ATOM     44  HA  THR A  69      14.647  15.056  30.363  1.00 21.27      A    H
ATOM     45  HB  THR A  69      11.388  15.356  31.910  1.00 18.85      A    H
ATOM     46 HG21 THR A  69      11.173  17.433  30.334  1.00 22.56      A    H
ATOM     47 HG22 THR A  69      12.747  16.477  32.550  1.00 22.56      A    H
ATOM     48 HG23 THR A  69      12.759  16.559  30.395  1.00 22.56      A    H
ATOM     49  H  ATHR A  69      13.914  11.773  30.712  0.50 18.60      A    H
ATOM     50  HG1ATHR A  69      10.933  14.249  27.986  0.50 16.82      A    H
ATOM     51  D  BTHR A  69      13.543  11.891  30.715  0.50 18.60      A    D
ATOM     52  DG1BTHR A  69      12.128  14.870  28.337  0.50 16.82      A    D
TER
END
"""

pdb_str_2 = """
CRYST1   34.238   35.096   43.858  90.00  90.00  90.00 P 21 21 21    0
ATOM      1  N   LYS A  45       6.154   2.754   1.212  1.00 12.39           N
ATOM      2  C   LYS A  45       7.533   2.537  -0.815  1.00  8.18           C
ATOM      3  O   LYS A  45       8.546   2.217  -1.437  1.00  8.18           O
ATOM      4  CA ALYS A  45       7.469   2.388   0.702  0.18  8.95           C
ATOM      5  CB ALYS A  45       7.820   0.954   1.105  0.18 15.56           C
ATOM      6  CG ALYS A  45       7.880   0.729   2.607  0.18  7.12           C
ATOM      7  CD ALYS A  45       8.227  -0.714   2.935  0.18 22.74           C
ATOM      8  CE ALYS A  45       8.402  -0.914   4.432  0.18 43.83           C
ATOM      9  NZ ALYS A  45       7.134  -0.676   5.175  0.18 58.88           N
ATOM     23  CA BLYS A  45       7.396   2.217   0.670  0.82  9.42           C
ATOM     24  CB BLYS A  45       7.467   0.704   0.891  0.82 15.14           C
ATOM     25  CG BLYS A  45       7.428   0.288   2.352  0.82  9.76           C
ATOM     26  CD BLYS A  45       7.539  -1.221   2.500  0.82 22.13           C
ATOM     27  CE BLYS A  45       6.396  -1.931   1.793  0.82 22.20           C
ATOM     28  NZ BLYS A  45       6.495  -3.411   1.927  0.82 22.09           N
ATOM      0  HA ALYS A  45       8.198   3.068   1.142  0.18  9.42           H   new
ATOM      0  HA BLYS A  45       8.223   2.691   1.198  0.82  9.42           H   new
ATOM      0  HB2ALYS A  45       7.082   0.276   0.676  0.18 15.14           H   new
ATOM      0  HB2BLYS A  45       6.636   0.232   0.368  0.82 15.14           H   new
ATOM      0  HB3ALYS A  45       8.784   0.691   0.670  0.18 15.14           H   new
ATOM      0  HB3BLYS A  45       8.384   0.324   0.441  0.82 15.14           H   new
ATOM      0  HG2ALYS A  45       8.624   1.393   3.047  0.18  9.76           H   new
ATOM      0  HG2BLYS A  45       8.244   0.769   2.891  0.82  9.76           H   new
ATOM      0  HG3ALYS A  45       6.919   0.985   3.054  0.18  9.76           H   new
ATOM      0  HG3BLYS A  45       6.499   0.633   2.806  0.82  9.76           H   new
ATOM      0  HD2ALYS A  45       7.439  -1.372   2.568  0.18 22.13           H   new
ATOM      0  HD2BLYS A  45       8.490  -1.560   2.089  0.82 22.13           H   new
ATOM      0  HD3ALYS A  45       9.145  -0.996   2.419  0.18 22.13           H   new
ATOM      0  HD3BLYS A  45       7.536  -1.486   3.557  0.82 22.13           H   new
ATOM      0  HE2ALYS A  45       8.751  -1.929   4.625  0.18 22.20           H   new
ATOM      0  HE2BLYS A  45       5.446  -1.592   2.207  0.82 22.20           H   new
ATOM      0  HE3ALYS A  45       9.172  -0.237   4.801  0.18 22.20           H   new
ATOM      0  HE3BLYS A  45       6.398  -1.661   0.737  0.82 22.20           H   new
ATOM      0  HZ1ALYS A  45       7.263  -0.936   6.174  0.18 22.09           H   new
ATOM      0  HZ1BLYS A  45       5.698  -3.859   1.432  0.82 22.09           H   new
ATOM      0  HZ2ALYS A  45       6.877   0.330   5.110  0.18 22.09           H   new
ATOM      0  HZ2BLYS A  45       7.390  -3.738   1.510  0.82 22.09           H   new
ATOM      0  HZ3ALYS A  45       6.375  -1.255   4.761  0.18 22.09           H   new
ATOM      0  HZ3BLYS A  45       6.467  -3.671   2.934  0.82 22.09           H   new
TER
END
"""

def run():
  if (not libtbx.env.has_module("reduce")):
    print("Reduce not installed, needed for model.idealize_h_minimization(). skipping")
    return
  for pdb_str in [pdb_str_1, pdb_str_2]:
    for use_neutron_distances in [True, False]:
      print("use_neutron_distances:", use_neutron_distances, "*"*30)
      params = mmtbx.model.manager.get_default_pdb_interpretation_params()
      params.pdb_interpretation.use_neutron_distances = use_neutron_distances
      inp = iotbx.pdb.input(lines=pdb_str, source_info=None)
      m = mmtbx.model.manager(
        model_input = inp,
        log         = null_out())
      m.process(pdb_interpretation_params=params, make_restraints=True)
      r1 = m.geometry_statistics()
      #m.setup_riding_h_manager()
      m.idealize_h_minimization(show=False)
      r2 = m.geometry_statistics()
      print("%6.3f %6.3f %6.3f %6.3f"%(
        r1.angle().mean,r1.bond().mean, r2.angle().mean,r2.bond().mean))
      assert r2.angle().mean < 1.0, "assertion %f < 1.0" % r2.angle().mean
      assert r2.bond().mean < 0.01, "assertion %f < 0.01" % r2.bond().mean

def exercise_2():
  inp = iotbx.pdb.input(lines=pdb_str_1, source_info=None)
  model = mmtbx.model.manager(
      model_input = inp)
  model.setup_scattering_dictionaries(scattering_table="wk1995")
  g_stat = model.geometry_statistics()

def exercise_3():
  """
  Make sure model can do deep_copy() straight away.
  """
  inp = iotbx.pdb.input(lines=pdb_str_1, source_info=None)
  model = mmtbx.model.manager(
      model_input = inp)
  new_model = model.deep_copy()

pdb_str_3 = """CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
SCALE1      0.050000  0.000000  0.000000        0.00000
SCALE2      0.000000  0.050000  0.000000        0.00000
SCALE3      0.000000  0.000000  0.050000        0.00000
ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00 30.00           C
ATOM      2  CA  GLY A   2       1.000   1.000   1.000  1.00 30.00           C
ATOM      3  CA  GLY A   3       2.000   2.000   2.000  1.00 30.00           C
ATOM      4  CA  GLY A   4       3.000   3.000   3.000  1.00 30.00           C
ATOM      5  CA  GLY A   5       4.000   4.000   4.000  1.00 30.00           C
ATOM      6  CA  GLY A   6       5.000   5.000   5.000  1.00 30.00           C
ATOM      7  CA  GLY A   7       6.000   6.000   6.000  1.00 30.00           C
ATOM      8  CA  GLY A   8       7.000   7.000   7.000  1.00 30.00           C
ATOM      9  CA  GLY A   9       8.000   8.000   8.000  1.00 30.00           C
ATOM     10  CA  GLY A  10       9.000   9.000   9.000  1.00 30.00           C
TER
END
"""

def exercise_convert_to_isotropic():
  pdb_str = """
ATOM    104  N   SER A   8      14.526  19.060  18.223  1.00  5.10           N
ANISOU  104  N   SER A   8      500    632    808   -107     58    104       N
ATOM    105  CA  SER A   8      14.099  17.792  18.758  1.00  5.95           C
"""
  cs = crystal.symmetry((5.01, 5.01, 5.47, 90, 90, 120), "P 62 2 2")
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  model = mmtbx.model.manager(
    model_input      = pdb_inp,
    crystal_symmetry = cs)
  #
  m = model.deep_copy()
  m.convert_to_isotropic()
  assert approx_equal(
    m.get_hierarchy().atoms().extract_b(), [5.105875, 5.95])
  assert approx_equal(
    m.get_hierarchy().atoms().extract_uij(),
    [(-1.0, -1.0, -1.0, -1.0, -1.0, -1.0), (-1.0, -1.0, -1.0, -1.0, -1.0, -1.0)])
  #
  m = model.deep_copy()
  m.convert_to_isotropic(selection=flex.bool([False, True]))
  assert approx_equal(
    m.get_hierarchy().atoms().extract_b(),
    model.get_hierarchy().atoms().extract_b())
  assert approx_equal(
    m.get_hierarchy().atoms().extract_uij(),
    model.get_hierarchy().atoms().extract_uij())
  #
  m1 = model.deep_copy()
  m2 = model.deep_copy()
  m1.convert_to_isotropic(selection=flex.bool([True, False]))
  m2.convert_to_isotropic()
  assert approx_equal(
    m1.get_hierarchy().atoms().extract_b(),
    m2.get_hierarchy().atoms().extract_b())
  assert approx_equal(
    m1.get_hierarchy().atoms().extract_uij(),
    m2.get_hierarchy().atoms().extract_uij())

def exercise_set_b_iso():
  pdb_str = """
ATOM    104  N   SER A   8      14.526  19.060  18.223  1.00  5.10           N
ANISOU  104  N   SER A   8      500    632    808   -107     58    104       N
ATOM    105  CA  SER A   8      14.099  17.792  18.758  1.00  5.95           C
"""
  cs = crystal.symmetry((5.01, 5.01, 5.47, 90, 90, 120), "P 62 2 2")
  uc = cs.unit_cell()
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  model = mmtbx.model.manager(
    model_input      = pdb_inp,
    crystal_symmetry = cs)
  #
  m = model.deep_copy()
  bs = flex.double([10,20])
  m.set_b_iso(values = bs)
  assert approx_equal(m.get_hierarchy().atoms().extract_b(), [-1.0, 20.0])
  us = m.get_hierarchy().atoms().extract_uij()
  assert approx_equal(us[0], model.get_hierarchy().atoms().extract_uij()[0])
  assert approx_equal(us[1], (-1.0, -1.0, -1.0, -1.0, -1.0, -1.0))
  u_iso_1 = m.get_xray_structure().scatterers().extract_u_iso()
  b_iso_2 = m.get_hierarchy().atoms().extract_b()
  assert approx_equal(u_iso_1[0], -1)
  assert approx_equal(b_iso_2[0], -1)
  assert approx_equal(adptbx.u_as_b(u_iso_1[1]), 20)
  assert approx_equal(b_iso_2[1], 20)
  #
  m = model.deep_copy()
  m.convert_to_isotropic()
  bs = flex.double([10,20])
  sel = flex.bool([True, False])
  m.set_b_iso(values=bs, selection=sel)
  assert approx_equal(m.get_b_iso(), [10, 5.95])

def exercise_from_sites_cart():
  from cctbx import crystal
  from scitbx.matrix import col
  from scitbx.array_family import flex
  sites_cart=flex.vec3_double()
  crystal_symmetry=crystal.symmetry(
      unit_cell=(20,20,20,90,90,90),
      space_group_symbol="P 1")
  for i in range(10):
    sites_cart.append(col((i,i,i)))
  model=mmtbx.model.manager.from_sites_cart(sites_cart=sites_cart,
      crystal_symmetry=crystal_symmetry)
  # print model.model_as_pdb()
  for i, a in enumerate(model.get_hierarchy().atoms()):
    assert a.i_seq == i, "iseqs were not set properly"
  assert model.model_as_pdb()==pdb_str_3

def exercise_has_hd():
  inp = iotbx.pdb.input(lines=pdb_str_1, source_info=None)
  m = mmtbx.model.manager(
      model_input               = inp,
      log                       = null_out())
  assert m.has_hd()
  inp = iotbx.pdb.input(lines=pdb_str_3, source_info=None)
  m = mmtbx.model.manager(
      model_input               = inp,
      log                       = null_out())
  assert not m.has_hd()

def exercise_ss_creation_crash():
  pdb_str = """
CRYST1  145.350  135.090  157.320  90.00  90.00  90.00 P 1
SCALE1      0.006880  0.000000  0.000000        0.00000
SCALE2      0.000000  0.007402  0.000000        0.00000
SCALE3      0.000000  0.000000  0.006356        0.00000
ATOM      1  N   ASN A   1      47.095 160.279  31.220  1.00 30.00           N
ATOM      2  CA  ASN A   1      65.985 120.233  34.727  1.00 30.00           C
ATOM      3  C   ASN A   1      56.657 138.700  33.374  1.00 30.00           C
ATOM      4  O   ASN A   1      56.353 138.977  34.561  1.00 30.00           O
ATOM      5  CB  ASN A   1      65.238 120.133  36.068  1.00 30.00           C
ATOM      6  CG  ASN A   1      66.087 119.360  37.057  1.00 30.00           C
ATOM      7  OD1 ASN A   1      65.746 118.217  37.441  1.00 30.00           O
ATOM      8  ND2 ASN A   1      67.240 119.920  37.395  1.00 30.00           N
ATOM      9  N   ASN A   2      56.939 137.441  33.021  1.00 30.00           N
ATOM     10  CA  ASN A   2      67.135 117.384  35.354  1.00 30.00           C
ATOM     11  C   ASN A   2      74.935 104.398  35.546  1.00 30.00           C
ATOM     12  O   ASN A   2      74.423 104.166  34.444  1.00 30.00           O
ATOM     13  CB  ASN A   2      65.828 116.703  35.809  1.00 30.00           C
ATOM     14  CG  ASN A   2      66.092 115.518  36.718  1.00 30.00           C
ATOM     15  OD1 ASN A   2      66.641 114.515  36.266  1.00 30.00           O
ATOM     16  ND2 ASN A   2      65.744 115.556  38.000  1.00 30.00           N
ATOM     17  N   ASN A   3      76.102 103.886  35.920  1.00 30.00           N
ATOM     18  CA  ASN A   3      68.960 115.076  35.163  1.00 30.00           C
ATOM     19  C   ASN A   3      86.047  90.376  35.591  1.00 30.00           C
ATOM     20  O   ASN A   3      87.134  90.903  35.535  1.00 30.00           O
ATOM     21  CB  ASN A   3      70.251 115.882  34.903  1.00 30.00           C
ATOM     22  CG  ASN A   3      71.023 116.208  36.192  1.00 30.00           C
ATOM     23  OD1 ASN A   3      70.637 117.096  36.957  1.00 30.00           O
ATOM     24  ND2 ASN A   3      72.106 115.481  36.436  1.00 30.00           N
ATOM     25  OXT ASN A   3      85.912  89.104  36.045  1.00 30.00           O
TER
END


"""
  with open("exercise_ss_creation_crash_model.pdb","w") as fo:
    fo.write(pdb_str)
  from iotbx.data_manager import DataManager
  dm=DataManager()
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.pdb_interpretation.secondary_structure.enabled=True
  model = dm.get_model('exercise_ss_creation_crash_model.pdb')
  model.process(make_restraints=True,
    pdb_interpretation_params=params)

def exercise_flip_nqh():
  """
  Expected flips:
     C 190  ASN
     D 212  ASN
     A 190  ASN
     B 212  ASN
  """
  pdb_str = """
CRYST1   19.391   14.196   16.627  90.00  90.00  90.00 P 1
SCALE1      0.051570  0.000000  0.000000        0.00000
SCALE2      0.000000  0.070442  0.000000        0.00000
SCALE3      0.000000  0.000000  0.060143        0.00000
ATOM      1  N   ASN A 190      14.391   7.875  10.959  1.00 25.08           N
ATOM      2  CA  ASN A 190      12.965   7.729  10.710  1.00 25.31           C
ATOM      3  C   ASN A 190      12.031   8.507  11.627  1.00 25.45           C
ATOM      4  O   ASN A 190      11.129   9.196  11.149  1.00 24.83           O
ATOM      5  CB  ASN A 190      12.578   6.250  10.760  1.00 29.15           C
ATOM      6  CG  ASN A 190      11.141   6.007  10.323  1.00 31.49           C
ATOM      7  OD1 ASN A 190      10.458   5.135  10.858  1.00 30.98           O
ATOM      8  ND2 ASN A 190      10.682   6.773   9.339  1.00 34.25           N
TER
ATOM      9  N   ASN B 212       7.583   7.891   7.004  1.00 63.62           N
ATOM     10  CA  ASN B 212       6.962   6.727   6.387  1.00 70.84           C
ATOM     11  C   ASN B 212       5.525   7.117   6.075  1.00 76.13           C
ATOM     12  O   ASN B 212       5.000   6.819   5.000  1.00 78.13           O
ATOM     13  CB  ASN B 212       6.978   5.532   7.341  1.00 71.04           C
ATOM     14  CG  ASN B 212       8.376   5.163   7.790  1.00 72.56           C
ATOM     15  OD1 ASN B 212       9.289   5.017   6.974  1.00 70.76           O
ATOM     16  ND2 ASN B 212       8.551   5.000   9.096  1.00 73.07           N
TER
ATOM      1  N   ASN C 190      24.391  17.875  20.959  1.00 25.08           N
ATOM      2  CA  ASN C 190      22.965  17.729  20.710  1.00 25.31           C
ATOM      3  C   ASN C 190      22.031  18.507  21.627  1.00 25.45           C
ATOM      4  O   ASN C 190      21.129  19.196  21.149  1.00 24.83           O
ATOM      5  CB  ASN C 190      22.578  16.250  20.760  1.00 29.15           C
ATOM      6  CG  ASN C 190      21.141  16.007  20.323  1.00 31.49           C
ATOM      7  OD1 ASN C 190      20.458  15.135  20.858  1.00 30.98           O
ATOM      8  ND2 ASN C 190      20.682  16.773  19.339  1.00 34.25           N
TER
ATOM      9  N   ASN D 212      17.583  17.891  17.004  1.00 63.62           N
ATOM     10  CA  ASN D 212      16.962  16.727  16.387  1.00 70.84           C
ATOM     11  C   ASN D 212      15.525  17.117  16.075  1.00 76.13           C
ATOM     12  O   ASN D 212      15.000  16.819  15.000  1.00 78.13           O
ATOM     13  CB  ASN D 212      16.978  15.532  17.341  1.00 71.04           C
ATOM     14  CG  ASN D 212      18.376  15.163  17.790  1.00 72.56           C
ATOM     15  OD1 ASN D 212      19.289  15.017  16.974  1.00 70.76           O
ATOM     16  ND2 ASN D 212      18.551  15.000  19.096  1.00 73.07           N
END
"""
  pdb_str_answer = """
CRYST1   19.391   14.196   16.627  90.00  90.00  90.00 P 1
SCALE1      0.051570  0.000000  0.000000        0.00000
SCALE2      0.000000  0.070442  0.000000        0.00000
SCALE3      0.000000  0.000000  0.060143        0.00000
ATOM      1  N   ASN A 190      14.391   7.875  10.959  1.00 25.08           N
ATOM      2  CA  ASN A 190      12.965   7.729  10.710  1.00 25.31           C
ATOM      3  C   ASN A 190      12.031   8.507  11.627  1.00 25.45           C
ATOM      4  O   ASN A 190      11.129   9.196  11.149  1.00 24.83           O
ATOM      5  CB  ASN A 190      12.578   6.250  10.760  1.00 29.15           C
ATOM      6  CG  ASN A 190      11.141   6.007  10.323  1.00 31.49           C
ATOM      7  OD1 ASN A 190      10.633   6.678   9.426  1.00 30.98           O
ATOM      8  ND2 ASN A 190      10.478   5.051  10.966  1.00 34.25           N
TER
ATOM      9  N   ASN B 212       7.583   7.891   7.004  1.00 63.62           N
ATOM     10  CA  ASN B 212       6.962   6.727   6.387  1.00 70.84           C
ATOM     11  C   ASN B 212       5.525   7.117   6.075  1.00 76.13           C
ATOM     12  O   ASN B 212       5.000   6.819   5.000  1.00 78.13           O
ATOM     13  CB  ASN B 212       6.978   5.532   7.341  1.00 71.04           C
ATOM     14  CG  ASN B 212       8.376   5.163   7.790  1.00 72.56           C
ATOM     15  OD1 ASN B 212       8.639   4.999   8.984  1.00 70.76           O
ATOM     16  ND2 ASN B 212       9.288   5.039   6.833  1.00 73.07           N
TER
ATOM     17  N   ASN C 190      24.391  17.875  20.959  1.00 25.08           N
ATOM     18  CA  ASN C 190      22.965  17.729  20.710  1.00 25.31           C
ATOM     19  C   ASN C 190      22.031  18.507  21.627  1.00 25.45           C
ATOM     20  O   ASN C 190      21.129  19.196  21.149  1.00 24.83           O
ATOM     21  CB  ASN C 190      22.578  16.250  20.760  1.00 29.15           C
ATOM     22  CG  ASN C 190      21.141  16.007  20.323  1.00 31.49           C
ATOM     23  OD1 ASN C 190      20.633  16.678  19.426  1.00 30.98           O
ATOM     24  ND2 ASN C 190      20.478  15.051  20.966  1.00 34.25           N
TER
ATOM     25  N   ASN D 212      17.583  17.891  17.004  1.00 63.62           N
ATOM     26  CA  ASN D 212      16.962  16.727  16.387  1.00 70.84           C
ATOM     27  C   ASN D 212      15.525  17.117  16.075  1.00 76.13           C
ATOM     28  O   ASN D 212      15.000  16.819  15.000  1.00 78.13           O
ATOM     29  CB  ASN D 212      16.978  15.532  17.341  1.00 71.04           C
ATOM     30  CG  ASN D 212      18.376  15.163  17.790  1.00 72.56           C
ATOM     31  OD1 ASN D 212      18.639  14.999  18.984  1.00 70.76           O
ATOM     32  ND2 ASN D 212      19.288  15.039  16.833  1.00 73.07           N
TER
END
"""
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  m = mmtbx.model.manager(model_input = pdb_inp)
  m.flip_nqh()
  sc_flipped = m.get_sites_cart()
  sc_orig    = iotbx.pdb.input(source_info=None, lines=pdb_str).atoms().extract_xyz()
  sc_answer  = iotbx.pdb.input(source_info=None, lines=pdb_str_answer).atoms().extract_xyz()
  assert flex.max(flex.sqrt((sc_flipped - sc_answer).dot())) < 0.001
  assert flex.max(flex.sqrt((sc_flipped - sc_orig).dot())) > 2.3

def exercise_macromolecule_plus_hetatms_by_chain_selections():
  pdb_str = """
CRYST1   29.071   35.014   35.025  90.00  90.00  90.00 P 1
ATOM      0  N   GLN A 595      15.775   6.909   5.000  1.00 77.72           N
ATOM      1  CA  GLN A 595      15.250   6.223   6.187  1.00 77.72           C
ATOM      2  C   GLN A 595      14.070   5.304   5.864  1.00 77.72           C
ATOM      3  O   GLN A 595      13.250   5.000   6.737  1.00 77.72           O
ATOM      4  CB  GLN A 595      14.885   7.217   7.303  1.00 77.72           C
ATOM      5  CG  GLN A 595      13.621   8.046   7.104  1.00 77.72           C
ATOM      6  CD  GLN A 595      13.921   9.467   6.669  1.00 77.72           C
ATOM      7  OE1 GLN A 595      15.074   9.827   6.433  1.00 77.72           O
ATOM      8  NE2 GLN A 595      12.881  10.285   6.561  1.00 77.72           N
TER
HETATM    9 MG    MG A4002     136.594 129.025 155.295  1.00 76.36          Mg
TER
HETATM   10  P   64T T  21       9.610  14.398   8.714  1.00 78.87           P
HETATM   11  OP1 64T T  21      10.939  13.757   9.035  1.00 78.87           O
HETATM   12  OP2 64T T  21       9.433  15.199   7.448  1.00 78.87           O
HETATM   13  O5' 64T T  21       8.476  13.255   8.731  1.00 78.87           O
HETATM   14  C5' 64T T  21       8.535  12.150   7.834  1.00 78.87           C
HETATM   15  C4' 64T T  21       9.475  11.087   8.384  1.00 78.87           C
HETATM   16  O4' 64T T  21       9.311  10.945   9.798  1.00 78.87           O
HETATM   17  C3' 64T T  21       9.192   9.737   7.756  1.00 78.87           C
HETATM   18  O3' 64T T  21      10.364   9.248   7.103  1.00 78.87           O
HETATM   19  C2' 64T T  21       8.794   8.819   8.892  1.00 78.87           C
HETATM   20  C1' 64T T  21       9.189   9.565  10.156  1.00 78.87           C
HETATM   21  N1  64T T  21       8.104   9.357  11.105  1.00 78.87           N
HETATM   22  C2  64T T  21       8.020   8.233  11.809  1.00 78.87           C
HETATM   23  O2  64T T  21       8.605   7.219  11.466  1.00 78.87           O
HETATM   24  N3  64T T  21       7.275   8.245  12.907  1.00 78.87           N
HETATM   25  C4  64T T  21       6.769   9.387  13.359  1.00 78.87           C
HETATM   26  O4  64T T  21       6.720   9.660  14.546  1.00 78.87           O
HETATM   27  C5  64T T  21       6.246  10.366  12.344  1.00 78.87           C
HETATM   28  C6  64T T  21       7.210  10.501  11.171  1.00 78.87           C
HETATM   29  C5M 64T T  21       6.024  11.715  13.011  1.00 78.87           C
HETATM   30  O5  64T T  21       5.000   9.880  11.839  1.00 78.87           O
TER
END
"""
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  m = mmtbx.model.manager(model_input = pdb_inp)
  sels = m.macromolecule_plus_hetatms_by_chain_selections()
  sels = [list(s.iselection()) for s in sels]
  assert len(sels)==3
  assert sels[0] == [0, 1, 2, 3, 4, 5, 6, 7, 8]
  assert sels[1] == [10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,
                     29, 30]
  assert sels[2] == [9]
  #
  sels = m.macromolecule_plus_hetatms_by_chain_selections(max_radius=20)
  sels = [list(s.iselection()) for s in sels]
  assert len(sels)==2
  assert sels[0] == [0, 1, 2, 3, 4, 5, 6, 7, 8]
  assert sels[1] == [9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,
                     29, 30]

def exercise_model_filename():
  for suffix in ['.pdb', '.cif']:
    with tempfile.NamedTemporaryFile(dir=os.getcwd(), suffix=suffix) as f:
      filename = f.name
      inp = iotbx.pdb.input(file_name=os.path.basename(filename))
      m = mmtbx.model.manager(model_input=inp)
      assert os.path.abspath(filename) == m.get_source_filename()
      assert os.path.basename(filename) == m.get_source_filename(full_path=False)

  inp = iotbx.pdb.input(source_info=None, lines=pdb_str_1)
  m = mmtbx.model.manager(model_input=inp)
  assert m.get_source_filename() is None

def exercise_altlocs_present():
  pdb_str1="""
CRYST1   13.135   16.057   12.855  90.00  90.00  90.00 P 1
ATOM      1  N   HIS A   3       8.135   7.414   7.855  1.00 22.05           N
ATOM      2  CA  HIS A   3       7.032   6.873   7.032  1.00 18.37           C
ATOM      3  C   HIS A   3       7.233   5.393   6.739  1.00 13.79           C
ATOM      4  O   HIS A   3       8.085   5.000   5.991  1.00 12.07           O
ATOM      5  CB  HIS A   3       6.816   7.572   5.708  1.00 21.10           C
ATOM      6  CG  HIS A   3       6.273   8.956   5.759  1.00 24.64           C
ATOM      7  ND1 HIS A   3       5.237   9.382   6.577  1.00 21.60           N
ATOM      8  CD2 HIS A   3       6.625  10.039   5.000  1.00 33.23           C
ATOM      9  CE1 HIS A   3       5.000  10.626   6.356  1.00 23.45           C
ATOM     10  NE2 HIS A   3       5.834  11.057   5.402  1.00 24.68           N
ATOM     11  OXT HIS A   3       6.509   4.555   7.276  1.00 13.79           O
ATOM     12  H   HIS A   3       7.879   8.235   8.404  1.00 22.05           H
ATOM     13  H2  HIS A   3       8.446   6.714   8.498  1.00 22.05           H
ATOM     14  H3  HIS A   3       8.894   7.681   7.262  1.00 22.05           H
ATOM     15  HA  HIS A   3       6.135   7.043   7.628  1.00 18.37           H
ATOM     16  HB2 HIS A   3       7.776   7.626   5.195  1.00 21.10           H
ATOM     17  HB3 HIS A   3       6.114   6.977   5.123  1.00 21.10           H
ATOM     18  HD2 HIS A   3       7.384  10.073   4.233  1.00 33.23           H
ATOM     19  HE1 HIS A   3       4.253  11.224   6.857  1.00 23.45           H

ATOM     20  HD1AHIS A   3       4.740   8.798   7.250  0.50 21.60           H
ATOM     21  HE2AHIS A   3       5.867  12.008   5.035  0.50 24.68           H

ATOM     22  DD1BHIS A   3       4.740   8.798   7.250  0.50 21.60           D
ATOM     23  DE2BHIS A   3       5.867  12.008   5.035  0.50 24.68           D
TER
END
  """
  pdb_str2="""
CRYST1   13.135   16.057   12.855  90.00  90.00  90.00 P 1
ATOM      1  N   HIS A   3       8.135   7.414   7.855  1.00 22.05           N
ATOM      2  CA AHIS A   3       7.032   6.873   7.032  1.00 18.37           C
ATOM      2  CA BHIS A   3       7.032   6.873   7.032  1.00 18.37           C
ATOM      3  C   HIS A   3       7.233   5.393   6.739  1.00 13.79           C
ATOM      4  O   HIS A   3       8.085   5.000   5.991  1.00 12.07           O
ATOM      5  CB  HIS A   3       6.816   7.572   5.708  1.00 21.10           C
ATOM      6  CG  HIS A   3       6.273   8.956   5.759  1.00 24.64           C
ATOM      7  ND1 HIS A   3       5.237   9.382   6.577  1.00 21.60           N
ATOM      8  CD2 HIS A   3       6.625  10.039   5.000  1.00 33.23           C
ATOM      9  CE1 HIS A   3       5.000  10.626   6.356  1.00 23.45           C
ATOM     10  NE2 HIS A   3       5.834  11.057   5.402  1.00 24.68           N
ATOM     11  OXT HIS A   3       6.509   4.555   7.276  1.00 13.79           O
ATOM     12  H   HIS A   3       7.879   8.235   8.404  1.00 22.05           H
ATOM     13  H2  HIS A   3       8.446   6.714   8.498  1.00 22.05           H
ATOM     14  H3  HIS A   3       8.894   7.681   7.262  1.00 22.05           H
ATOM     15  HA  HIS A   3       6.135   7.043   7.628  1.00 18.37           H
ATOM     16  HB2 HIS A   3       7.776   7.626   5.195  1.00 21.10           H
ATOM     17  HB3 HIS A   3       6.114   6.977   5.123  1.00 21.10           H
ATOM     18  HD2 HIS A   3       7.384  10.073   4.233  1.00 33.23           H
ATOM     19  HE1 HIS A   3       4.253  11.224   6.857  1.00 23.45           H

ATOM     20  HD1AHIS A   3       4.740   8.798   7.250  0.50 21.60           H
ATOM     21  HE2AHIS A   3       5.867  12.008   5.035  0.50 24.68           H

ATOM     22  DD1BHIS A   3       4.740   8.798   7.250  0.50 21.60           D
ATOM     23  DE2BHIS A   3       5.867  12.008   5.035  0.50 24.68           D
TER
END
  """
  pdb_str3="""
CRYST1   13.135   16.057   12.855  90.00  90.00  90.00 P 1
ATOM      1  N   HIS A   3       8.135   7.414   7.855  1.00 22.05           N
ATOM      2  CA AHIS A   3       7.032   6.873   7.032  1.00 18.37           C
ATOM      2  CA BHIS A   3       7.032   6.873   7.032  1.00 18.37           C
ATOM      3  C   HIS A   3       7.233   5.393   6.739  1.00 13.79           C
ATOM      4  O   HIS A   3       8.085   5.000   5.991  1.00 12.07           O
ATOM      5  CB  HIS A   3       6.816   7.572   5.708  1.00 21.10           C
ATOM      6  CG  HIS A   3       6.273   8.956   5.759  1.00 24.64           C
ATOM      7  ND1 HIS A   3       5.237   9.382   6.577  1.00 21.60           N
ATOM      8  CD2 HIS A   3       6.625  10.039   5.000  1.00 33.23           C
ATOM      9  CE1 HIS A   3       5.000  10.626   6.356  1.00 23.45           C
ATOM     10  NE2 HIS A   3       5.834  11.057   5.402  1.00 24.68           N
ATOM     11  OXT HIS A   3       6.509   4.555   7.276  1.00 13.79           O
ATOM     12  H   HIS A   3       7.879   8.235   8.404  1.00 22.05           H
ATOM     13  H2  HIS A   3       8.446   6.714   8.498  1.00 22.05           H
ATOM     14  H3  HIS A   3       8.894   7.681   7.262  1.00 22.05           H
ATOM     15  HA  HIS A   3       6.135   7.043   7.628  1.00 18.37           H
ATOM     16  HB2 HIS A   3       7.776   7.626   5.195  1.00 21.10           H
ATOM     17  HB3 HIS A   3       6.114   6.977   5.123  1.00 21.10           H
ATOM     18  HD2 HIS A   3       7.384  10.073   4.233  1.00 33.23           H
ATOM     19  HE1 HIS A   3       4.253  11.224   6.857  1.00 23.45           H
TER
END
  """
  pdb_str4="""
CRYST1   13.135   16.057   12.855  90.00  90.00  90.00 P 1
ATOM      1  N   HIS A   3       8.135   7.414   7.855  1.00 22.05           N
ATOM      2  CA  HIS A   3       7.032   6.873   7.032  1.00 18.37           C
ATOM      3  C   HIS A   3       7.233   5.393   6.739  1.00 13.79           C
ATOM      4  O   HIS A   3       8.085   5.000   5.991  1.00 12.07           O
ATOM      5  CB  HIS A   3       6.816   7.572   5.708  1.00 21.10           C
ATOM      6  CG  HIS A   3       6.273   8.956   5.759  1.00 24.64           C
ATOM      7  ND1 HIS A   3       5.237   9.382   6.577  1.00 21.60           N
ATOM      8  CD2 HIS A   3       6.625  10.039   5.000  1.00 33.23           C
ATOM      9  CE1 HIS A   3       5.000  10.626   6.356  1.00 23.45           C
ATOM     10  NE2 HIS A   3       5.834  11.057   5.402  1.00 24.68           N
ATOM     11  OXT HIS A   3       6.509   4.555   7.276  1.00 13.79           O
ATOM     12  H   HIS A   3       7.879   8.235   8.404  1.00 22.05           H
ATOM     13  H2  HIS A   3       8.446   6.714   8.498  1.00 22.05           H
ATOM     14  H3  HIS A   3       8.894   7.681   7.262  1.00 22.05           H
ATOM     15  HA  HIS A   3       6.135   7.043   7.628  1.00 18.37           H
ATOM     16  HB2 HIS A   3       7.776   7.626   5.195  1.00 21.10           H
ATOM     17  HB3 HIS A   3       6.114   6.977   5.123  1.00 21.10           H
ATOM     18  HD2 HIS A   3       7.384  10.073   4.233  1.00 33.23           H
ATOM     19  HE1 HIS A   3       4.253  11.224   6.857  1.00 23.45           H
TER
END
  """
  pdb_inp1 = iotbx.pdb.input(source_info=None, lines=pdb_str1)
  m1 = mmtbx.model.manager(model_input = pdb_inp1)
  #
  pdb_inp2 = iotbx.pdb.input(source_info=None, lines=pdb_str2)
  m2 = mmtbx.model.manager(model_input = pdb_inp2)
  #
  pdb_inp3 = iotbx.pdb.input(source_info=None, lines=pdb_str3)
  m3 = mmtbx.model.manager(model_input = pdb_inp3)
  #
  pdb_inp4 = iotbx.pdb.input(source_info=None, lines=pdb_str4)
  m4 = mmtbx.model.manager(model_input = pdb_inp4)
  #
  assert m1.altlocs_present()
  assert m1.altlocs_present_only_hd()
  #
  assert m2.altlocs_present()
  assert not m2.altlocs_present_only_hd()
  #
  assert m3.altlocs_present()
  assert not m3.altlocs_present_only_hd()
  #
  assert not m4.altlocs_present()
  assert not m4.altlocs_present_only_hd()


if (__name__ == "__main__"):
  t0 = time.time()
  exercise_altlocs_present()
  exercise_macromolecule_plus_hetatms_by_chain_selections()
  exercise_ss_creation_crash()
  exercise_set_b_iso()
  exercise_convert_to_isotropic()
  run()
  exercise_2()
  exercise_3()
  exercise_from_sites_cart()
  exercise_has_hd()
  exercise_flip_nqh()
  exercise_model_filename()
  print("Time: %6.3f"%(time.time()-t0))
