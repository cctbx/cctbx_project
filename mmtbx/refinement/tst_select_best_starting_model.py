from __future__ import absolute_import, division, print_function

from libtbx.test_utils import Exception_expected
from libtbx.utils import null_out, Sorry
from libtbx import easy_run
from scitbx.array_family import flex
from cctbx import uctbx, sgtbx
import iotbx.pdb
from iotbx import file_reader
from mmtbx.refinement import select_best_starting_model
import random
from six.moves import range

def exercise_main():
  unit_cell = (24.937, 8.866, 25.477, 90.00, 107.08, 90.00)
  space_group = "P21"
  pdb_base = """\
CRYST1   24.937    8.866   25.477  90.00 107.08  90.00 P 1 21 1
SCALE1      0.040101  0.000000  0.012321        0.00000
SCALE2      0.000000  0.112790  0.000000        0.00000
SCALE3      0.000000  0.000000  0.041062        0.00000
ATOM      1  N   GLY A   1       8.992   0.474  -6.096  1.00 16.23           N
ATOM      2  CA  GLY A   1       9.033   0.047  -4.707  1.00 16.20           C
ATOM      3  C   GLY A   1       7.998  -1.029  -4.448  1.00 15.91           C
ATOM      4  O   GLY A   1       7.548  -1.689  -5.385  1.00 16.11           O
ATOM      5  N   ASN A   2       7.625  -1.218  -3.185  1.00 15.02           N
ATOM      6  CA  ASN A   2       6.523  -2.113  -2.848  1.00 13.92           C
ATOM      7  C   ASN A   2       5.220  -1.618  -3.428  1.00 12.24           C
ATOM      8  O   ASN A   2       4.955  -0.418  -3.432  1.00 11.42           O
ATOM      9  CB  ASN A   2       6.376  -2.261  -1.340  1.00 14.42           C
ATOM     10  CG  ASN A   2       7.620  -2.786  -0.697  1.00 13.92           C
ATOM     11  OD1 ASN A   2       8.042  -3.915  -0.978  1.00 14.39           O
ATOM     12  ND2 ASN A   2       8.232  -1.975   0.168  1.00 12.78           N
ATOM     13  N   ASN A   3       4.406  -2.553  -3.904  1.00 12.20           N
ATOM     14  CA  ASN A   3       3.164  -2.226  -4.594  1.00 11.81           C
ATOM     15  C   ASN A   3       1.925  -2.790  -3.910  1.00 10.59           C
ATOM     16  O   ASN A   3       1.838  -3.991  -3.653  1.00 10.32           O
ATOM     17  CB  ASN A   3       3.231  -2.727  -6.046  1.00 12.51           C
ATOM     18  CG  ASN A   3       1.973  -2.405  -6.848  1.00 12.59           C
ATOM     19  OD1 ASN A   3       1.662  -1.239  -7.106  1.00 13.64           O
ATOM     20  ND2 ASN A   3       1.260  -3.443  -7.268  1.00 12.39           N
ATOM     21  N   GLN A   4       0.973  -1.913  -3.608  1.00 10.34           N
ATOM     22  CA  GLN A   4      -0.366  -2.335  -3.208  1.00 10.00           C
ATOM     23  C   GLN A   4      -1.402  -1.637  -4.085  1.00 10.21           C
ATOM     24  O   GLN A   4      -1.514  -0.414  -4.070  1.00  8.99           O
ATOM     25  CB  GLN A   4      -0.656  -2.027  -1.736  1.00 10.00           C
ATOM     26  CG  GLN A   4      -1.927  -2.705  -1.229  1.00 10.50           C
ATOM     27  CD  GLN A   4      -2.482  -2.102   0.060  1.00 11.36           C
ATOM     28  OE1 GLN A   4      -2.744  -0.900   0.151  1.00 12.29           O
ATOM     29  NE2 GLN A   4      -2.684  -2.951   1.055  1.00 10.43           N
ATOM     30  N   GLN A   5      -2.154  -2.406  -4.857  1.00 10.48           N
ATOM     31  CA  GLN A   5      -3.247  -1.829  -5.630  1.00 11.24           C
ATOM     32  C   GLN A   5      -4.591  -2.382  -5.178  1.00 11.40           C
ATOM     33  O   GLN A   5      -4.789  -3.599  -5.092  1.00 11.94           O
ATOM     34  CB  GLN A   5      -3.024  -2.023  -7.129  1.00 11.14           C
ATOM     35  CG  GLN A   5      -1.852  -1.222  -7.653  1.00 10.65           C
ATOM     36  CD  GLN A   5      -1.338  -1.748  -8.965  1.00 10.73           C
ATOM     37  OE1 GLN A   5      -0.794  -2.845  -9.028  1.00 10.14           O
ATOM     38  NE2 GLN A   5      -1.511  -0.968 -10.027  1.00 11.31           N
ATOM     39  N   ASN A   6      -5.504  -1.471  -4.872  1.00 11.56           N
ATOM     40  CA  ASN A   6      -6.809  -1.838  -4.359  1.00 12.07           C
ATOM     41  C   ASN A   6      -7.856  -1.407  -5.353  1.00 13.18           C
ATOM     42  O   ASN A   6      -8.257  -0.251  -5.362  1.00 13.64           O
ATOM     43  CB  ASN A   6      -7.053  -1.149  -3.017  1.00 12.12           C
ATOM     44  CG  ASN A   6      -5.966  -1.446  -1.998  1.00 12.31           C
ATOM     45  OD1 ASN A   6      -5.833  -2.579  -1.517  1.00 13.43           O
ATOM     46  ND2 ASN A   6      -5.198  -0.423  -1.645  1.00 11.88           N
ATOM     47  N   TYR A   7      -8.298  -2.332  -6.193  1.00 14.34           N
ATOM     48  CA  TYR A   7      -9.162  -1.980  -7.317  1.00 15.00           C
ATOM     49  C   TYR A   7     -10.603  -1.792  -6.893  1.00 15.64           C
ATOM     50  O   TYR A   7     -11.013  -2.278  -5.838  1.00 15.68           O
ATOM     51  CB  TYR A   7      -9.064  -3.041  -8.412  1.00 15.31           C
ATOM     52  CG  TYR A   7      -7.657  -3.197  -8.931  1.00 15.06           C
ATOM     53  CD1 TYR A   7      -6.785  -4.118  -8.368  1.00 15.24           C
ATOM     54  CD2 TYR A   7      -7.193  -2.400  -9.960  1.00 14.96           C
ATOM     55  CE1 TYR A   7      -5.489  -4.253  -8.830  1.00 14.94           C
ATOM     56  CE2 TYR A   7      -5.905  -2.526 -10.429  1.00 15.13           C
ATOM     57  CZ  TYR A   7      -5.055  -3.451  -9.861  1.00 14.97           C
ATOM     58  OH  TYR A   7      -3.768  -3.572 -10.335  1.00 14.93           O
ATOM     59  OXT TYR A   7     -11.378  -1.149  -7.601  1.00 15.89           O
TER
"""
  pdb_base_water = """\
HETATM   64  O   HOH S   1     -10.466  -2.347  -3.168  1.00 17.57           O
HETATM   65  O   HOH S   2       6.469   1.081  -7.070  1.00 21.27           O
HETATM   66  O   HOH S   3     -11.809   0.108  -9.956  1.00 27.52           O
HETATM   67  O   HOH S   4       1.580  -3.455 -11.035  1.00 44.76           O
END
"""
  with open("tst_start_model_base.pdb", "w") as f:
    f.write(pdb_base+pdb_base_water)
  params = """
    high_resolution = 1.75
    add_sigmas = True
    output {
      label = F
      type = *real complex
      file_name = tst_start_model_base.mtz
    }
  """
  with open("tst_start_model_fmodel.eff", "w") as f:
    f.write(params)
  assert (easy_run.fully_buffered(
    "phenix.fmodel tst_start_model_base.pdb tst_start_model_fmodel.eff"
  ).raise_if_errors().return_code == 0)
  mtz_in = file_reader.any_file("tst_start_model_base.mtz")
  f_obs = mtz_in.file_server.miller_arrays[0]
  symm = f_obs.crystal_symmetry().customized_copy(
    space_group_info=sgtbx.space_group_info("P2"))
  f_obs = f_obs.customized_copy(crystal_symmetry=symm)
  random.seed(12345) # XXX makes results more predictable
  flags = f_obs.generate_r_free_flags(fraction=0.1)
  mtz_data = f_obs.as_mtz_dataset(
    column_root_label="F")
  mtz_data.add_miller_array(flags,
    column_root_label="FreeR_flag")
  mtz_data.mtz_object().write("tst_start_model.mtz")
  pdb_in = iotbx.pdb.input("tst_start_model_base.pdb")
  hierarchy_in = pdb_in.construct_hierarchy()
  xrs_in = pdb_in.xray_structure_simple()
  selection = hierarchy_in.atom_selection_cache().selection
  # Model 1: very few changes, but shifted by (1,0,0.5)
  symm2 = xrs_in.crystal_symmetry().customized_copy(
    unit_cell=uctbx.unit_cell((24.932, 8.841, 25.501, 90.00, 107.5, 90.00)))
  #u_iso = xrs_in.extract_u_iso_or_u_equiv()
  #xrs_out = xrs_in.deep_copy_scatterers().set_u_iso(
  #  selection=flex.bool(u_iso.size(), True), values=u_iso*1.1)
  xrs_out = xrs_in.deep_copy_scatterers()
  sites_cart = xrs_out.sites_cart()
  sites_cart += flex.vec3_double(sites_cart.size(), (1.0, 0.0, 0.5))
  xrs_out.set_sites_cart(sites_cart)
  hierarchy_out = hierarchy_in.deep_copy()
  hierarchy_out.adopt_xray_structure(xrs_out)
  with open("tst_start_model_1.pdb", "w") as f:
    f.write(hierarchy_out.as_pdb_string(crystal_symmetry=symm2))
  # Model 2: no sidechains
  mc_sele = selection("(name N or name C or name O or name CA or name CB)")
  hierarchy_out = hierarchy_in.select(mc_sele)
  xrs_out = xrs_in.select(mc_sele)
  with open("tst_start_model_2.pdb", "w") as f:
    f.write(hierarchy_out.as_pdb_string(crystal_symmetry=xrs_out))
  # Model 3: P1 symmetry
  symm3 = xrs_in.crystal_symmetry().customized_copy(
    space_group_info=sgtbx.space_group_info("P1"))
  with open("tst_start_model_3.pdb", "w") as f:
    f.write(hierarchy_out.as_pdb_string(crystal_symmetry=symm3))
  # Model 4: shaken coordinates and ADPs
  def random_double(size, factor=1):
    d = flex.double()
    for x in range(size):
      d.append(random.random() * factor)
    return d
  xrs_out = xrs_in.customized_copy()
  xrs_out.shake_sites_in_place(0.3, random_double=random_double)
  xrs_out.shake_adp()
  hierarchy_out = hierarchy_in.deep_copy()
  hierarchy_out.adopt_xray_structure(xrs_out)
  with open("tst_start_model_4.pdb", "w") as f:
    f.write(hierarchy_out.as_pdb_string(crystal_symmetry=xrs_out))
  # Model 5: perfect, but missing CRYST1
  with open("tst_start_model_5.pdb", "w") as f:
    f.write(hierarchy_in.as_pdb_string())
  # run method
  params = select_best_starting_model.master_phil.extract()
  params.rigid_body_refine = False
  model_names = [
    "tst_start_model_1.pdb",
    "tst_start_model_2.pdb",
    "tst_start_model_3.pdb",
    "tst_start_model_4.pdb",
    "tst_start_model_5.pdb",
  ]
  result = select_best_starting_model.select_model(
    model_names=model_names,
    model_data=None,
    f_obs=f_obs,
    r_free_flags=flags,
    params=params,
    skip_twin_detection=True,
    log=null_out())
  # result.show(verbose=True)
  assert (result.best_model_name == "tst_start_model_4.pdb"), result.best_model_name
  params.rigid_body_refine = True
  result = select_best_starting_model.select_model(
    model_names=model_names,
    model_data=None,
    f_obs=f_obs,
    r_free_flags=flags,
    params=params,
    skip_twin_detection=True,
    log=null_out())
  # result.show(verbose=True)
  assert (result.best_model_name == "tst_start_model_1.pdb"), result.best_model_name

def exercise_misc():
  pdb_str = """\
REMARK this is a remark record!
CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1
SCALE1      0.045585  0.000000  0.014006        0.00000
SCALE2      0.000000  0.205508  0.000000        0.00000
SCALE3      0.000000  0.000000  0.044560        0.00000
ATOM      1  N   GLY A   1      -9.009   4.612   6.102  0.00 16.77           N
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  0.00 16.57           C
ATOM      3  C   GLY A   1      -8.015   3.140   4.419  0.00 16.16           C
ATOM      4  O   GLY A   1      -7.523   2.521   5.381  0.00 16.78           O
ATOM      5  H1  GLY A   1      -9.802   4.938   6.343  0.00 16.77           H
ATOM      6  H2  GLY A   1      -8.816   3.902   6.603  0.00 16.77           H
ATOM      7  H3  GLY A   1      -8.385   5.236   6.218  0.00 16.77           H
ATOM      8  HA2 GLY A   1      -9.928   3.856   4.426  0.00 16.57           H
ATOM      9  HA3 GLY A   1      -8.858   4.970   4.084  0.00 16.57           H
ATOM     10  N   ASN A   2      -7.656   2.923   3.155  1.00 15.02           N
ATOM     11  CA  ASN A   2      -6.522   2.038   2.831  1.00 14.10           C
ATOM     12  C   ASN A   2      -5.241   2.537   3.427  1.00 13.13           C
ATOM     13  O   ASN A   2      -4.978   3.742   3.426  1.00 11.91           O
ATOM     14  CB  ASN A   2      -6.346   1.881   1.341  1.00 15.38           C
ATOM     15  CG  ASN A   2      -7.584   1.342   0.692  1.00 14.08           C
ATOM     16  OD1 ASN A   2      -8.025   0.227   1.016  1.00 17.46           O
ATOM     17  ND2 ASN A   2      -8.204   2.155  -0.169  1.00 11.72           N
ATOM     18  H   ASN A   2      -8.044   3.269   2.470  1.00 15.02           H
ATOM     19  HA  ASN A   2      -6.698   1.159   3.202  1.00 14.10           H
ATOM     20  HB2 ASN A   2      -6.150   2.746   0.949  1.00 15.38           H
ATOM     21  HB3 ASN A   2      -5.619   1.262   1.169  1.00 15.38           H
ATOM     22 HD21 ASN A   2      -8.919   1.893  -0.569  1.00 11.72           H
ATOM     23 HD22 ASN A   2      -7.888   2.940  -0.323  1.00 11.72           H
ATOM     24  N  AASN A   3      -4.438   1.590   3.905  0.50 12.26           N
ATOM     25  CA AASN A   3      -3.193   1.904   4.589  0.50 11.74           C
ATOM     26  C  AASN A   3      -1.955   1.332   3.895  0.50 11.10           C
ATOM     27  O  AASN A   3      -1.872   0.119   3.648  0.50 10.42           O
ATOM     28  CB AASN A   3      -3.259   1.378   6.042  0.50 12.15           C
ATOM     29  CG AASN A   3      -2.006   1.739   6.861  0.50 12.82           C
ATOM     30  OD1AASN A   3      -1.702   2.925   7.072  0.50 15.05           O
ATOM     31  ND2AASN A   3      -1.271   0.715   7.306  0.50 13.48           N
ATOM     32  H  AASN A   3      -4.597   0.747   3.843  0.50 12.26           H
ATOM     33  HA AASN A   3      -3.083   2.866   4.633  0.50 11.74           H
ATOM     34  HB2AASN A   3      -4.029   1.767   6.485  0.50 12.15           H
ATOM     35  HB3AASN A   3      -3.339   0.411   6.025  0.50 12.15           H
ATOM     36 HD21AASN A   3      -0.560   0.864   7.766  0.50 13.48           H
ATOM     37 HD22AASN A   3      -1.509  -0.093   7.132  0.50 13.48           H
ATOM     38  N  BASN A   3      -4.338   1.590   3.905  0.50 12.26           N
ATOM     39  CA BASN A   3      -3.093   1.904   4.589  0.50 11.74           C
ATOM     40  C  BASN A   3      -1.855   1.332   3.895  0.50 11.10           C
ATOM     41  O  BASN A   3      -1.772   0.119   3.648  0.50 10.42           O
ATOM     42  CB BASN A   3      -3.159   1.378   6.042  0.50 12.15           C
ATOM     43  CG BASN A   3      -4.127   2.189   6.923  0.50 12.82           C
ATOM     44  OD1BASN A   3      -4.598   1.573   8.012  0.50 15.05           O
ATOM     45  ND2BASN A   3      -4.430   3.358   6.630  0.50 13.48           N
ATOM     46  H  BASN A   3      -4.464   0.743   3.823  0.50 12.26           H
ATOM     47  HA BASN A   3      -2.983   2.866   4.621  0.50 11.74           H
ATOM     48  HB2BASN A   3      -3.463   0.457   6.031  0.50 12.15           H
ATOM     49  HB3BASN A   3      -2.275   1.432   6.438  0.50 12.15           H
ATOM     50 HD21BASN A   3      -3.922   3.821   6.115  0.50 13.48           H
ATOM     51 HD22BASN A   3      -5.147   3.708   6.952  0.50 13.48           H
ATOM     52  N   GLN A   4      -1.005   2.228   3.598  1.00 10.29           N
ATOM     53  CA  GLN A   4       0.384   1.888   3.199  1.00 10.53           C
ATOM     54  C   GLN A   4       1.435   2.606   4.088  1.00 10.24           C
ATOM     55  O   GLN A   4       1.547   3.843   4.115  1.00  8.86           O
ATOM     56  CB  GLN A   4       0.656   2.148   1.711  1.00  9.80           C
ATOM     57  CG  GLN A   4       1.944   1.458   1.213  1.00 10.25           C
ATOM     58  CD  GLN A   4       2.504   2.044  -0.089  1.00 12.43           C
ATOM     59  OE1 GLN A   4       2.744   3.268  -0.190  1.00 14.62           O
ATOM     60  NE2 GLN A   4       2.750   1.161  -1.091  1.00  9.05           N
ATOM     61  HA  GLN A   4       0.509   0.941   3.335  1.00 10.53           H
ATOM     62  HB2 GLN A   4      -0.088   1.808   1.189  1.00  9.80           H
ATOM     63  HB3 GLN A   4       0.752   3.103   1.569  1.00  9.80           H
ATOM     64  HG2 GLN A   4       2.630   1.544   1.893  1.00 10.25           H
ATOM     65  HG3 GLN A   4       1.753   0.520   1.057  1.00 10.25           H
ATOM     66 HE21 GLN A   4       2.592   0.324  -0.975  1.00  9.05           H
ATOM     67 HE22 GLN A   4       3.064   1.439  -1.842  1.00  9.05           H
ATOM     68  H  AGLN A   4      -1.142   3.077   3.620  0.50 10.29           H
ATOM     69  H  BGLN A   4      -1.168   3.072   3.606  0.50 10.29           H
ATOM     70  N   GLN A   5       2.154   1.821   4.871  1.00 10.38           N
ATOM     71  CA  GLN A   5       3.270   2.361   5.640  1.00 11.39           C
ATOM     72  C   GLN A   5       4.594   1.768   5.172  1.00 11.52           C
ATOM     73  O   GLN A   5       4.768   0.546   5.054  1.00 12.05           O
ATOM     74  CB  GLN A   5       3.056   2.183   7.147  1.00 11.96           C
ATOM     75  CG  GLN A   5       1.829   2.950   7.647  1.00 10.81           C
ATOM     76  CD  GLN A   5       1.344   2.414   8.954  1.00 13.10           C
ATOM     77  OE1 GLN A   5       0.774   1.325   9.002  1.00 10.65           O
ATOM     78  NE2 GLN A   5       1.549   3.187  10.039  1.00 12.30           N
ATOM     79  H   GLN A   5       2.021   0.978   4.977  1.00 10.38           H
ATOM     80  HA  GLN A   5       3.333   3.315   5.487  1.00 11.39           H
ATOM     81  HB2 GLN A   5       2.924   1.242   7.341  1.00 11.96           H
ATOM     82  HB3 GLN A   5       3.834   2.516   7.621  1.00 11.96           H
ATOM     83  HG2 GLN A   5       2.063   3.883   7.769  1.00 10.81           H
ATOM     84  HG3 GLN A   5       1.111   2.866   7.000  1.00 10.81           H
ATOM     85 HE21 GLN A   5       1.943   3.947   9.957  1.00 12.30           H
ATOM     86 HE22 GLN A   5       1.286   2.920  10.813  1.00 12.30           H
ATOM     87  N   ASN A   6       5.514   2.664   4.856  1.00 11.99           N
ATOM     88  CA  ASN A   6       6.831   2.310   4.318  1.00 12.30           C
ATOM     89  C   ASN A   6       7.854   2.761   5.324  1.00 13.40           C
ATOM     90  O   ASN A   6       8.219   3.943   5.374  1.00 13.92           O
ATOM     91  CB  ASN A   6       7.065   3.016   2.993  1.00 12.13           C
ATOM     92  CG  ASN A   6       5.961   2.735   2.003  1.00 12.77           C
ATOM     93  OD1 ASN A   6       5.798   1.604   1.551  1.00 14.27           O
ATOM     94  ND2 ASN A   6       5.195   3.747   1.679  1.00 10.07           N
ATOM     95  H   ASN A   6       5.401   3.512   4.945  1.00 11.99           H
ATOM     96  HA  ASN A   6       6.909   1.352   4.189  1.00 12.30           H
ATOM     97  HB2 ASN A   6       7.098   3.974   3.140  1.00 12.13           H
ATOM     98  HB3 ASN A   6       7.901   2.707   2.610  1.00 12.13           H
ATOM     99 HD21 ASN A   6       4.553   3.638   1.117  1.00 10.07           H
ATOM    100 HD22 ASN A   6       5.334   4.521   2.028  1.00 10.07           H
ATOM    101  N   MSE A   7       8.292   1.817   6.147  1.00 14.70           N
ATOM    102  CA  MSE A   7       9.159   2.144   7.299  1.00 15.18           C
ATOM    103  C   MSE A   7      10.603   2.331   6.885  1.00 15.91           C
ATOM    104  O   MSE A   7      11.041   1.811   5.855  1.00 15.76           O
ATOM    105  CB  MSE A   7       9.071   1.027   8.337  1.00 20.00           C
ATOM    106  CG  MSE A   7       7.703   0.904   8.987  1.00 20.00           C
ATOM    107 SE   MSE A   7       7.114   2.484   9.902  1.00 20.00          Se
ATOM    108  CE  MSE A   7       8.400   2.453  11.367  1.00 20.00           C
ATOM    109  H   MSE A   7       8.109   0.980   6.070  1.00 14.70           H
ATOM    110  HA  MSE A   7       8.841   2.961   7.708  1.00 15.18           H
ATOM    111  HB2 MSE A   7       9.270   0.186   7.900  1.00 20.00           H
ATOM    112  HB3 MSE A   7       9.724   1.195   9.033  1.00 20.00           H
ATOM    113  HG2 MSE A   7       7.042   0.539   8.377  1.00 20.00           H
ATOM    114  HG3 MSE A   7       8.123   0.264   9.563  1.00 20.00           H
ATOM    115  HE1 MSE A   7       8.976   3.219  11.304  1.00 20.00           H
ATOM    116  HE2 MSE A   7       8.333   2.133  12.198  1.00 20.00           H
ATOM    117  HE3 MSE A   7       9.010   1.805  10.983  1.00 20.00           H
TER
HETATM  118 CL    CL     1      -6.471   5.227   7.124  1.00 30.00          Cl
TER
HETATM  119  O   HOH A   9      10.431   1.858   3.216  1.00 19.71           O
HETATM  120  O   HOH A  10     -11.286   1.756  -1.468  1.00 17.08           O
HETATM  121  O   HOH A  11      11.808   4.179   9.970  1.00 23.99           O
HETATM  122  O   HOH A  12      13.605   1.327   9.198  1.00 26.17           O
HETATM  123  O   HOH A  13      -2.749   3.429  10.024  1.00 39.15           O
HETATM  124  O   HOH A  14      -1.500   0.682  10.967  1.00 43.49           O
TER
END
"""
  pdb_file = "tst_start_model_misc.pdb"
  with open(pdb_file, "w") as f:
    f.write(pdb_str)
  pdb_in = iotbx.pdb.input(file_name=pdb_file)
  pdb_hierarchy = pdb_in.construct_hierarchy()
  xray_structure = pdb_in.xray_structure_simple()
  hd_sel = xray_structure.hd_selection()
  assert hd_sel.count(True) == 55
  xray_structure.convert_to_anisotropic(selection=~hd_sel)
  pdb_hierarchy.adopt_xray_structure(xray_structure)
  hierarchy_new, xrs_new = select_best_starting_model.strip_model(
    pdb_hierarchy=pdb_hierarchy,
    xray_structure=xray_structure,
    output_file="tst_start_model_misc_new.pdb",
    log=null_out())
  # basic changes
  hd_sel = xray_structure.hd_selection()
  assert hd_sel.count(True) == 55 # should be left alone
  assert xrs_new.hd_selection().count(True) == 0
  sel = hierarchy_new.atom_selection_cache().selection
  mse_sel = sel("resname MSE")
  assert mse_sel.count(True) == 0
  hoh_sel = sel("resname HOH")
  assert hoh_sel.count(True) == 0
  alt_sel = sel("altloc 'B'")
  assert alt_sel.count(True) == 0
  cl_sel = sel("element CL")
  assert cl_sel.count(True) == 1
  occ_xrs = xrs_new.scatterers().extract_occupancies()
  occ_pdb = hierarchy_new.atoms().extract_occ()
  assert occ_xrs.all_eq(occ_pdb) and occ_xrs.all_eq(1.0)
  assert xrs_new.use_u_aniso().all_eq(False)
  # remove all ligands
  assert xray_structure.hd_selection().count(True) == 55
  hierarchy_new, xrs_new = select_best_starting_model.strip_model(
    pdb_hierarchy=pdb_hierarchy,
    xray_structure=xray_structure,
    remove_ligands=True,
    remove_hydrogens=False,
    output_file="tst_start_model_misc_new2.pdb",
    log=null_out())
  sel = hierarchy_new.atom_selection_cache().selection
  assert sel("resname CL").count(True) == 0
  assert xrs_new.hd_selection().count(True) == 48
  # multi-model hierarchy
  model2 = pdb_hierarchy.models()[0].detached_copy()
  pdb_hierarchy.append_model(model2)
  try :
    hierarchy_new, xrs_new = select_best_starting_model.strip_model(
      pdb_hierarchy=pdb_hierarchy,
      xray_structure=xray_structure,
      log=null_out())
  except Sorry :
    pass
  else :
    raise Exception_expected
  # both structure and file name passed
  try :
    hierarchy_new, xrs_new = select_best_starting_model.strip_model(
      pdb_hierarchy=pdb_hierarchy,
      xray_structure=xray_structure,
      file_name=pdb_file,
      log=null_out())
  except AssertionError :
    pass
  else :
    raise Exception_expected
  # preserve REMARK
  hierarchy_new, xrs_new = select_best_starting_model.strip_model(
    file_name=pdb_file,
    preserve_remarks=True,
    output_file="tst_start_model_misc_new.pdb",
    log=null_out())
  pdb_new = iotbx.pdb.input(file_name="tst_start_model_misc_new.pdb")
  assert len(pdb_new.remark_section()) == 1
  # reset HETATM
  hierarchy_new, xrs_new = select_best_starting_model.strip_model(
    file_name=pdb_file,
    preserve_remarks=True,
    reset_hetatm_flag=True,
    output_file="tst_start_model_misc_new.pdb",
    log=null_out())
  with open("tst_start_model_misc_new.pdb") as f:
    lines = f.readlines()
  for line in lines :
    assert (not line.startswith("HETATM"))

if (__name__ == "__main__"):
  exercise_misc()
  exercise_main()
  print("OK")
