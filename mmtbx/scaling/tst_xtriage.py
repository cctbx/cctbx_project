
from __future__ import absolute_import, division, print_function
from mmtbx.scaling import data_statistics as ds
from mmtbx.scaling import xtriage
from mmtbx.command_line import fmodel
from iotbx import file_reader
import iotbx.pdb
from cctbx import crystal
from cctbx import miller
from cctbx import sgtbx
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal, Exception_expected, show_diff
from libtbx.development import show_pickle_sizes
from libtbx.easy_pickle import dumps, loads
from libtbx.utils import null_out, Sorry
import libtbx.load_env
from six.moves import cStringIO as StringIO
import warnings
import os.path
import sys
from six.moves import range

# synthetic data
def exercise_1():
  pdb_raw = """\
CRYST1   23.000    6.666   25.000  90.00 107.08  90.00 P 1 21 1      2
ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      3  C   GLY A   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.523   2.521   5.381  1.00 16.78           O
ATOM      5  N   ASN A   2      -7.656   2.923   3.155  1.00 15.02           N
ATOM      6  CA  ASN A   2      -6.522   2.038   2.831  1.00 14.10           C
ATOM      7  C   ASN A   2      -5.241   2.537   3.427  1.00 13.13           C
ATOM      8  O   ASN A   2      -4.978   3.742   3.426  1.00 11.91           O
ATOM      9  CB  ASN A   2      -6.346   1.881   1.341  1.00 15.38           C
ATOM     10  CG  ASN A   2      -7.584   1.342   0.692  1.00 14.08           C
ATOM     11  OD1 ASN A   2      -8.025   0.227   1.016  1.00 17.46           O
ATOM     12  ND2 ASN A   2      -8.204   2.155  -0.169  1.00 11.72           N
ATOM     13  N   ASN A   3      -4.438   1.590   3.905  1.00 12.26           N
ATOM     14  CA  ASN A   3      -3.193   1.904   4.589  1.00 11.74           C
ATOM     15  C   ASN A   3      -1.955   1.332   3.895  1.00 11.10           C
ATOM     16  O   ASN A   3      -1.872   0.119   3.648  1.00 10.42           O
ATOM     17  CB  ASN A   3      -3.259   1.378   6.042  1.00 12.15           C
ATOM     18  CG  ASN A   3      -2.006   1.739   6.861  1.00 12.82           C
ATOM     19  OD1 ASN A   3      -1.702   2.925   7.072  1.00 15.05           O
ATOM     20  ND2 ASN A   3      -1.271   0.715   7.306  1.00 13.48           N
ATOM     21  N   MET A   4      -1.005   2.228   3.598  1.00 10.29           N
ATOM     22  CA  MET A   4       0.384   1.888   3.199  1.00 10.53           C
ATOM     23  C   MET A   4       1.435   2.606   4.088  1.00 10.24           C
ATOM     24  O   MET A   4       1.547   3.843   4.115  1.00  8.86           O
ATOM     25  CB  MET A   4       0.616   2.241   1.729  1.00 20.00           C
ATOM     26  CG  MET A   4      -0.207   1.416   0.754  1.00 20.00           C
ATOM     27  SD  MET A   4       0.132  -0.349   0.876  1.00 20.00           S
ATOM     28  CE  MET A   4       1.822  -0.411   0.285  1.00 20.00           C
ATOM     29  N   GLN A   5       2.154   1.821   4.871  1.00 10.38           N
ATOM     30  CA  GLN A   5       3.270   2.361   5.640  1.00 11.39           C
ATOM     31  C   GLN A   5       4.594   1.768   5.172  1.00 11.52           C
ATOM     32  O   GLN A   5       4.768   0.546   5.054  1.00 12.05           O
ATOM     33  CB  GLN A   5       3.056   2.183   7.147  1.00 11.96           C
ATOM     34  CG  GLN A   5       1.829   2.950   7.647  1.00 10.81           C
ATOM     35  CD  GLN A   5       1.344   2.414   8.954  1.00 13.10           C
ATOM     36  OE1 GLN A   5       0.774   1.325   9.002  1.00 10.65           O
ATOM     37  NE2 GLN A   5       1.549   3.187  10.039  1.00 12.30           N
ATOM     38  N   ASN A   6       5.514   2.664   4.856  1.00 11.99           N
ATOM     39  CA  ASN A   6       6.831   2.310   4.318  1.00 12.30           C
ATOM     40  C   ASN A   6       7.854   2.761   5.324  1.00 13.40           C
ATOM     41  O   ASN A   6       8.219   3.943   5.374  1.00 13.92           O
ATOM     42  CB  ASN A   6       7.065   3.016   2.993  1.00 12.13           C
ATOM     43  CG  ASN A   6       5.961   2.735   2.003  1.00 12.77           C
ATOM     44  OD1 ASN A   6       5.798   1.604   1.551  1.00 14.27           O
ATOM     45  ND2 ASN A   6       5.195   3.747   1.679  1.00 10.07           N
ATOM     46  N   TYR A   7       8.292   1.817   6.147  1.00 14.70           N
ATOM     47  CA  TYR A   7       9.159   2.144   7.299  1.00 15.18           C
ATOM     48  C   TYR A   7      10.603   2.331   6.885  1.00 15.91           C
ATOM     49  O   TYR A   7      11.041   1.811   5.855  1.00 15.76           O
ATOM     50  CB  TYR A   7       9.061   1.065   8.369  1.00 15.35           C
ATOM     51  CG  TYR A   7       7.665   0.929   8.902  1.00 14.45           C
ATOM     52  CD1 TYR A   7       6.771   0.021   8.327  1.00 15.68           C
ATOM     53  CD2 TYR A   7       7.210   1.756   9.920  1.00 14.80           C
ATOM     54  CE1 TYR A   7       5.480  -0.094   8.796  1.00 13.46           C
ATOM     55  CE2 TYR A   7       5.904   1.649  10.416  1.00 14.33           C
ATOM     56  CZ  TYR A   7       5.047   0.729   9.831  1.00 15.09           C
ATOM     57  OH  TYR A   7       3.766   0.589  10.291  1.00 14.39           O
ATOM     58  OXT TYR A   7      11.358   2.999   7.612  1.00 17.49           O
TER      59      TYR A   7
HETATM    1 CA    CA A   8      10.431   1.858   3.216  1.00 30.00          CA
HETATM   60  O   HOH A   9      -6.471   5.227   7.124  1.00 22.62           O
HETATM   62  O   HOH A  10     -11.286   1.756  -1.468  1.00 17.08           O
HETATM   63  O   HOH A  11      11.808   4.179   9.970  1.00 23.99           O
HETATM   64  O   HOH A  12      13.605   1.327   9.198  1.00 26.17           O
HETATM   65  O   HOH A  13      -2.749   3.429  10.024  1.00 39.15           O
HETATM   66  O   HOH A  14      -1.500   0.682  10.967  1.00 43.49           O
END
"""
  pdb_file = "tst_xtriage_in.pdb"
  with open(pdb_file, "w") as f:
    f.write(pdb_raw)
  fmodel_args = [
    pdb_file,
    "high_resolution=1.5",
    "k_sol=0.35",
    "b_sol=20",
    "wavelength=1.54",
    "add_random_error_to_amplitudes_percent=3",
    "random_seed=12345",
    "output.type=real",
    "output.label=F",
    "output.file_name=tst_xtriage_fmodel.mtz",
  ]

  #  read it instead so python3 will be the same
  #  fmodel.run(args=fmodel_args, log=null_out())
  hkl_file = libtbx.env.find_in_repositories(
    relative_path="mmtbx/regression/mtz/tst_xtriage_fmodel.mtz",
    test=os.path.isfile)
  mtz_in = file_reader.any_file(
    hkl_file).assert_file_type("hkl")
  f_obs = mtz_in.file_server.miller_arrays[0].remove_cone(0.1)
  data = f_obs.data()
  # add some outliers
  #data[17] = 20
  #data[334] = 26
  #data[1908] = 13
  # and sigmas
  sigf = flex.double(f_obs.size(), 0.1) + (f_obs.data() * 0.03)
  f_obs = f_obs.customized_copy(sigmas=sigf)
  mtz_file = "tst_xtriage_in.mtz"
  f_obs.as_mtz_dataset(column_root_label="F").mtz_object().write(mtz_file)
  seq_file = "tst_xtriage_in.fa"
  with open(seq_file, "w") as f:
    f.write("> tst_xtriage\nGNNMQNY")

  # check with completeness_as_non_anomalous=True

  xtriage_args = [
    mtz_file,
    pdb_file,
    seq_file,
    "log=tst_xtriage_1.log",
    "l_test_dhkl=2,2,2",
    "completeness_as_non_anomalous=True",
  ]
  result = xtriage.run(args=xtriage_args, out=null_out())
  test_pickle_consistency_and_size(result)
  assert (result.matthews.n_copies == 1)
  assert (str(result.matthews.table) == """\
Solvent content analysis
Copies             Solvent content    Matthews coeff.    P(solvent content)
1                  0.472              2.33               1.000
""")
  data_strength = result.data_strength_and_completeness
  assert approx_equal(data_strength.data_strength.resolution_cut, 1.5351,
    eps=0.001)
  out1 = data_strength.low_resolution_completeness.format()
  assert (out1 == """\
---------------------------------------------------------
| Resolution range  | N(obs)/N(possible) | Completeness |
---------------------------------------------------------
| 21.9858 - 10.4368 | [6/7]              | 0.857        |
| 10.4368 -  8.4369 | [3/3]              | 1.000        |
|  8.4369 -  7.4172 | [3/4]              | 0.750        |
|  7.4172 -  6.7606 | [4/4]              | 1.000        |
|  6.7606 -  6.2882 | [5/5]              | 1.000        |
|  6.2882 -  5.9252 | [3/4]              | 0.750        |
|  5.9252 -  5.6337 | [7/7]              | 1.000        |
|  5.6337 -  5.3922 | [5/5]              | 1.000        |
|  5.3922 -  5.1874 | [4/4]              | 1.000        |
|  5.1874 -  5.0106 | [4/4]              | 1.000        |
---------------------------------------------------------"""), out1
  # ANOMALOUS SIGNAL
  a_meas = result.anomalous_info.measurability
  #assert approx_equal(a_meas.high_d_cut, 4.7636, eps=0.0001) # Why it's None?
  assert approx_equal(a_meas.low_d_cut, 2.3566, eps=0.0001)
  # ABSOLUTE SCALING
  ws = result.wilson_scaling
  assert ("%.2f" % ws.iso_p_scale) == "0.65",ws.iso_p_scale
  assert ("%.2f" % ws.iso_b_wilson) == "14.42", ws.iso_b_wilson
  # FIXME these may need to be adjusted for different hardware/OS
  assert approx_equal(ws.aniso_p_scale, 0.64723, eps=0.001)
  assert approx_equal(ws.aniso_u_star, [0.00034229, 0.00475982, 0.000285989,
                                        -0.0, 8.95386085999e-05, 0.0])
  assert approx_equal(ws.aniso_b_cart, (13.218423, 16.840142, 12.948426,
    1.0354e-15, -0.0685311, -7.92862e-16), 0.3)
  # convenience methods for GUI
  assert approx_equal(result.aniso_b_min, 12.895580)
  assert approx_equal(result.aniso_range_of_b, 3.804215)
  #
  assert approx_equal(ws.outlier_shell_table.data[0], # d_spacing
    [9.865131, 8.369653, 4.648634])
  assert approx_equal(ws.outlier_shell_table.data[1], # z_score
    [5.306713, 18.068284, 5.319230])
  assert (len(ws.outliers.acentric_outliers_table.data[0]) == 2)
  assert (ws.outliers.acentric_outliers_table.data[1] == [(0,-1,-1), (0,1,1)])
  assert approx_equal(ws.outliers.acentric_outliers_table.data[2],
    [3.507247, 3.315550])
  assert (ws.outliers.centric_outliers_table.data is None)
  assert (len(ws.ice_rings.table._rows) == 10)
  assert (ws.ice_rings.table._rows[0] ==
          ['    3.897', '     1.000', '   0.76', '   1.00']), \
          ws.ice_rings.table._rows[0]
  tw = result.twin_results
  wm = tw.wilson_moments
  out = StringIO()
  wm.show(out)
  assert not show_diff(out.getvalue(), """
                  ----------Wilson ratio and moments----------

Acentric reflections:


   <I^2>/<I>^2    :2.063   (untwinned: 2.000; perfect twin 1.500)
   <F>^2/<F^2>    :0.778   (untwinned: 0.785; perfect twin 0.885)
   <|E^2 - 1|>    :0.745   (untwinned: 0.736; perfect twin 0.541)

Centric reflections:


   <I^2>/<I>^2    :3.076   (untwinned: 3.000; perfect twin 2.000)
   <F>^2/<F^2>    :0.628   (untwinned: 0.637; perfect twin 0.785)
   <|E^2 - 1|>    :0.999   (untwinned: 0.968; perfect twin 0.736)

""")
  # XXX PDB validation server
  assert approx_equal(result.iso_b_wilson, 14.51, eps=0.1)
  assert approx_equal(result.aniso_b_ratio, 0.271, eps=0.1)
  assert (result.number_of_wilson_outliers == 2)
  assert approx_equal(result.l_test_mean_l, 0.481, eps=0.1)
  assert approx_equal(result.l_test_mean_l_squared, 0.322, eps=0.1)
  assert approx_equal(result.i_over_sigma_outer_shell, 10.71, eps=0.01)
  assert ("indicating pseudo-translationa" in result.patterson_verdict)
  # check relative Wilson
  # FIXME
  #result.relative_wilson.show()
  #assert (result.relative_wilson.n_outliers() == 0)
  #show_pickled_object_sizes(result)
  #

  # check with completeness_as_non_anomalous=False

  xtriage_args = [
    mtz_file,
    pdb_file,
    seq_file,
    "log=tst_xtriage_1.log",
    "l_test_dhkl=2,2,2",
    "completeness_as_non_anomalous=False",
  ]
  result = xtriage.run(args=xtriage_args, out=null_out())
  test_pickle_consistency_and_size(result)
  assert (result.matthews.n_copies == 1)
  assert (str(result.matthews.table) == """\
Solvent content analysis
Copies             Solvent content    Matthews coeff.    P(solvent content)
1                  0.472              2.33               1.000
""")
  data_strength = result.data_strength_and_completeness
  assert approx_equal(data_strength.data_strength.resolution_cut, 1.5351,
    eps=0.001)
  out1 = data_strength.low_resolution_completeness.format()
  assert (out1 == """\
---------------------------------------------------------
| Resolution range  | N(obs)/N(possible) | Completeness |
---------------------------------------------------------
| 21.9858 - 10.4368 | [ 6/7 ]            | 0.857        |
| 10.4368 -  8.4369 | [ 3/3 ]            | 1.000        |
|  8.4369 -  7.4172 | [ 3/4 ]            | 0.750        |
|  7.4172 -  6.7606 | [ 4/4 ]            | 1.000        |
|  6.7606 -  6.2882 | [ 8/8 ]            | 1.000        |
|  6.2882 -  5.9252 | [ 4/5 ]            | 0.800        |
|  5.9252 -  5.6337 | [11/11]            | 1.000        |
|  5.6337 -  5.3922 | [ 7/7 ]            | 1.000        |
|  5.3922 -  5.1874 | [ 6/6 ]            | 1.000        |
|  5.1874 -  5.0106 | [ 7/7 ]            | 1.000        |
---------------------------------------------------------"""), out1
  # ANOMALOUS SIGNAL
  a_meas = result.anomalous_info.measurability
  #assert approx_equal(a_meas.high_d_cut, 4.7636, eps=0.0001) # Why?
  assert approx_equal(a_meas.low_d_cut, 2.3565, eps=0.0001)
  # ABSOLUTE SCALING
  ws = result.wilson_scaling
  assert ("%.2f" % ws.iso_p_scale) == "0.65", ws.iso_p_scale
  assert ("%.2f" % ws.iso_b_wilson) == "14.42", ws.iso_b_wilson
  # FIXME these may need to be adjusted for different hardware/OS
  assert approx_equal(ws.aniso_p_scale, 0.64723, eps=0.001)
  assert approx_equal(ws.aniso_u_star, [0.00034473, 0.00479983, 0.000287162,
                                        -0.0, 9.00962e-05, 0.0], 6.e-5)
  assert approx_equal(ws.aniso_b_cart, [13.12, 16.69, 12.89,
    0, -0.08, 0], 0.01)
  # convenience methods for GUI
  assert approx_equal(result.aniso_b_min, 12.9, 0.1)
  assert approx_equal(result.aniso_range_of_b, 3.8, 0.1)
  #
  assert approx_equal(ws.outlier_shell_table.data[0], # d_spacing
    [9.86, 8.36, 4.64], 0.02)
  assert approx_equal(ws.outlier_shell_table.data[1], # z_score
    [5.30, 18.06, 5.31], 0.01)
  assert (len(ws.outliers.acentric_outliers_table.data[0]) == 2)
  assert (ws.outliers.acentric_outliers_table.data[1] == [(0,-1,-1), (0,1,1)])
  assert approx_equal(ws.outliers.acentric_outliers_table.data[2],
    [3.5, 3.3], 0.1)
  assert (ws.outliers.centric_outliers_table.data is None)
  assert (len(ws.ice_rings.table._rows) == 10)
  assert (ws.ice_rings.table._rows[0] ==
          ['    3.897', '     1.000', '   0.76', '   1.00']), \
          ws.ice_rings.table._rows[0]
  tw = result.twin_results
  wm = tw.wilson_moments
  out = StringIO()
  wm.show(out)
  assert not show_diff(out.getvalue(), """
                  ----------Wilson ratio and moments----------

Acentric reflections:


   <I^2>/<I>^2    :2.063   (untwinned: 2.000; perfect twin 1.500)
   <F>^2/<F^2>    :0.778   (untwinned: 0.785; perfect twin 0.885)
   <|E^2 - 1|>    :0.745   (untwinned: 0.736; perfect twin 0.541)

Centric reflections:


   <I^2>/<I>^2    :3.076   (untwinned: 3.000; perfect twin 2.000)
   <F>^2/<F^2>    :0.628   (untwinned: 0.637; perfect twin 0.785)
   <|E^2 - 1|>    :0.999   (untwinned: 0.968; perfect twin 0.736)

""")
  # XXX PDB validation server
  assert approx_equal(result.iso_b_wilson, 14.51, eps=0.1)
  assert approx_equal(result.aniso_b_ratio, 0.271, eps=0.1)
  assert (result.number_of_wilson_outliers == 2)
  assert approx_equal(result.l_test_mean_l, 0.481, eps=0.1)
  assert approx_equal(result.l_test_mean_l_squared, 0.322, eps=0.1)
  assert approx_equal(result.i_over_sigma_outer_shell, 10.71, eps=0.01)
  assert ("indicating pseudo-translationa" in result.patterson_verdict)
  # check relative Wilson
  # FIXME
  #result.relative_wilson.show()
  #assert (result.relative_wilson.n_outliers() == 0)
  #show_pickled_object_sizes(result)
  #
  # test without sigmas
  f_obs_2 = f_obs.customized_copy(sigmas=None)
  mtz_file = "tst_xtriage_in_2.mtz"
  f_obs_2.as_mtz_dataset(column_root_label="F").mtz_object().write(mtz_file)
  xtriage_args = [
    mtz_file,
    pdb_file,
    seq_file,
    "log=tst_xtriage_1.log",
  ]


  result = xtriage.run(args=xtriage_args, out=null_out())
  result.summarize_issues()
  # test in lower symmetry
  f_obs_3 = f_obs.expand_to_p1()
  mtz_file = "tst_xtriage_in_3.mtz"
  f_obs_3.as_mtz_dataset(column_root_label="F").mtz_object().write(mtz_file)
  xtriage_args = [
    mtz_file,
    seq_file,
    "log=tst_xtriage_2.log",
  ]
  result = xtriage.run(args=xtriage_args, out=null_out())
  assert ((1, 'One or more symmetry operators suggest that the data has a higher crystallographic symmetry (P 2 1 1).', 'Point group and R-factor analysis') in result.summarize_issues()._issues)
  # test with elliptical truncation
  f_obs_3 = f_obs.customized_copy(
    crystal_symmetry=crystal.symmetry((23,5,20,90,107.8,90), "P 21"))
  f_obs_3 = f_obs_3.resolution_filter(d_min=1.5)
  f_obs_3 = f_obs_3.customized_copy(crystal_symmetry=f_obs.crystal_symmetry())
  reso = ds.analyze_resolution_limits(f_obs_3)
  out = StringIO()
  reso.show(out=out)
  assert ("max. difference between axes = 0.652" in out.getvalue()), \
    out.getvalue()
  assert ("elliptically truncated" in out.getvalue())
  # make sure the elliptical truncation detection still works in higher space
  # groups - we only need a miller.set for this
  miller_set = miller.build_set(
    crystal_symmetry=crystal.symmetry((20,20,20,90,90,90), "P422"),
    d_min=1.5,
    anomalous_flag=False)
  reso = ds.analyze_resolution_limits(miller_set)
  out = StringIO()
  reso.show(out=out)
  assert ("Resolution limits are within expected tolerances" in out.getvalue())
  # log binning
  out = StringIO()
  log_binned = ds.log_binned_completeness(f_obs_3)
  log_binned.show(out=out)
  assert ("""| 1.9724 - 1.5094  | 368/1230    | 29.9%        |""" in
          out.getvalue()), out.getvalue()
  # test with no acentrics
  cf = f_obs.centric_flags().data()
  centrics = f_obs.select(cf)
  acentrics = f_obs.select(~cf)
  mtz_file = "tst_xtriage_in_3.mtz"
  centrics.as_mtz_dataset(column_root_label="F").mtz_object().write(mtz_file)
  args = [
    mtz_file,
    pdb_file,
    seq_file,
    "log=tst_xtriage_3.log",
  ]
  try :
    xtriage.run(args=args, out=null_out())
  except Sorry :
    pass
  else :
    raise Exception_expected
  # with only a handful of acentrics
  sel = flex.bool(acentrics.size(), False)
  for i in range(10):
    sel[i] = True
  f_obs_4 = centrics.concatenate(acentrics.select(sel))
  f_obs_4.as_mtz_dataset(column_root_label="F").mtz_object().write(mtz_file)
  try :
    xtriage.run(args=args, out=null_out())
  except Sorry :
    pass
  else :
    raise Exception_expected

# XXX code for debugging pickle size issues
def show_pickled_object_sizes(result):
  result_pkl = dumps(result)
  print("result", len(result_pkl))
  show_pickle_sizes(result, "  ")

# test consistency of output after pickling and unpickling
def test_pickle_consistency_and_size(result):
  all_out = StringIO()
  result.show(out=all_out)
  result_pkl_str = dumps(result)
  pkl_size = len(result_pkl_str)
  if (pkl_size >= 100000):
    print("Oversized pickle:", pkl_size)
    show_pickled_object_sizes(result)
    raise OverflowError()
  assert (pkl_size < 100000), pkl_size # limit pickle size
  result_pkl = loads(result_pkl_str)
  all_out_pkl = StringIO()
  result_pkl.show(out=all_out_pkl)
  assert not show_diff(all_out.getvalue(), all_out_pkl.getvalue())

def exercise_analyze_resolution_limits():
  for x in range(1, 231):
    sg = sgtbx.space_group_info(number=x)
    #sg = sgtbx.space_group_info("P222")
    uc = sg.any_compatible_unit_cell(80000)
    ms = miller.build_set(
      crystal_symmetry=crystal.symmetry(
        space_group_info=sg,
        unit_cell=uc),
      anomalous_flag=True,
      d_min=1.5)
    arl = ds.analyze_resolution_limits(ms)
    if (x > 2):
      assert (arl.max_d_min_delta() < 0.1)

# real data
def exercise_2():
  hkl_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/wizards/data/p9_se_w2.sca",
    test=os.path.isfile)
  if (hkl_file is None):
    warnings.warn("phenix_regression not available, skipping test")
    return
  hkl_in = file_reader.any_file(hkl_file).assert_file_type("hkl")
  i_obs_raw = hkl_in.file_object.as_miller_arrays(
    merge_equivalents=False,
    crystal_symmetry=crystal.symmetry(
      space_group_symbol="I4",
      unit_cell=(113.949,113.949,32.474,90,90,90)))[0]
  i_obs = i_obs_raw.merge_equivalents().array()
  # completeness and data strength
  cstats = ds.i_sigi_completeness_stats(i_obs)
  d_min_cut = cstats.resolution_cut
  assert approx_equal(d_min_cut, 2.150815)
  ws = ds.wilson_scaling(
    miller_array=i_obs,
    n_residues=120)
  # outliers - this shouldn't actually work, since it requires additional
  # processing steps on the input data
  try :
    outliers = ds.possible_outliers(i_obs)
  except AssertionError :
    pass
  else :
    raise Exception_expected
  ######################################################################
  # OVERALL ANALYSIS
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_examples/p9-build/p9.pdb",
    test=os.path.isfile)
  f_calc = None
  if (pdb_file is not None):
    pdb_in = iotbx.pdb.input(pdb_file)
    hierarchy = pdb_in.construct_hierarchy()
    xrs = pdb_in.xray_structure_simple(
      crystal_symmetry=i_obs)
    f_calc = xrs.structure_factors(d_min=i_obs.d_min()).f_calc()
    f_calc = abs(f_calc).generate_bijvoet_mates()
    f_calc = f_calc.set_observation_type_xray_amplitude()
    i_obs, f_calc = i_obs.common_sets(other=f_calc)
    with open("tmp_xtriage.pdb", "w") as f:
      f.write(hierarchy.as_pdb_string(crystal_symmetry=i_obs))
    pdb_file = "tmp_xtriage.pdb"
  params = xtriage.master_params.extract()
  params.scaling.input.asu_contents.n_residues = 141
  text_out = open("logfile3.log", "w")
  result = xtriage.xtriage_analyses(
    miller_obs=i_obs,
    miller_calc=f_calc,
    params=params,
    unmerged_obs=i_obs_raw,
    text_out=text_out)#sys.stdout)
  text_out.close()
  # XXX there appears to be some system-dependence here, hence sloppy limits
  assert (15.5 < result.aniso_b_min < 15.9)
  assert (10 < result.aniso_range_of_b < 11)
  # check relative Wilson
  if (pdb_file is not None):
    assert (result.relative_wilson is not None)
    # FIXME
    #assert (result.relative_wilson.n_outliers() == 34)
  #show_pickled_object_sizes(result)
  test_pickle_consistency_and_size(result)
  # XXX PDB validation server
  assert approx_equal(result.iso_b_wilson, 18.33, eps=0.1)
  assert approx_equal(result.aniso_b_ratio, 0.546, eps=0.1)
  assert (result.number_of_wilson_outliers == 0)
  assert approx_equal(result.l_test_mean_l, 0.493, eps=0.1)
  assert approx_equal(result.l_test_mean_l_squared, 0.326, eps=0.1)
  assert approx_equal(result.i_over_sigma_outer_shell, 3.25, eps=0.1)
  assert approx_equal(result.overall_i_sig_i,10.34,eps=0.1)
  assert approx_equal(result.anomalous_info.plan_sad_experiment_stats.get_overall(
      item="i_over_sigma_dict"),10.61,eps=0.1)
  assert approx_equal(result.anomalous_info.plan_sad_experiment_stats.get_overall(
      item="anom_signal_dict"),15.35,eps=0.1)
  assert ("No significant pseudotranslation is detected" in
          result.patterson_verdict)
  # test consistency of output after pickling and unpickling
  try :
    from phenix_dev.phenix_cloud import xtriage_json
  except ImportError :
    pass
  else :
    json_out = xtriage_json.json_output("p9.sca")
    result.show(out=json_out)
    with open("xtriage.json", "w") as f:
      f.write(json_out.export())
  # unmerged data
  assert result.merging_stats is not None
  out = StringIO()
  result.merging_stats.show(out=out)
  assert ("R-merge: 0.073" in out.getvalue())
  assert approx_equal(result.estimate_d_min(min_i_over_sigma=10), 1.9645,
    eps=0.001)
  # FIXME PDB doesn't actually have unit cell!
  # test detection of symmetry in reference file
  if (pdb_file is not None):
    args = [hkl_file, pdb_file]
    result = xtriage.run(args=args, out=null_out())

def exercise_loggraph():
  mtz_file = libtbx.env.find_in_repositories(
    relative_path="mmtbx/regression/mtz/tst_xtriage_fmodel.mtz",
    test=os.path.isfile)
  result = xtriage.run(
    args=[mtz_file, 'scaling.input.parameters.reporting.loggraphs=True'],
    out=null_out()
  )

if (__name__ == "__main__"):
  #exercise_2()
  exercise_1()
  exercise_analyze_resolution_limits()
  exercise_loggraph()
  print("OK")
