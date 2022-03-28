from __future__ import absolute_import, division, print_function
import time, sys, os
from six.moves import cStringIO as StringIO
import libtbx.load_env
import mmtbx.model
import iotbx.pdb
from mmtbx.monomer_library import pdb_interpretation
from cctbx.array_family import flex


def run():
  verbose = "--verbose" in sys.argv[1:]
  exercise_00(verbose=verbose)
  exercise_01(verbose=verbose)
  exercise_02(verbose=verbose)
  exercise_03(verbose=verbose)

def get_model(file_name, log, pdb_str=None, make_restraints=True):
  pdb_interpretation_params = iotbx.phil.parse(
          input_string=pdb_interpretation.grand_master_phil_str, process_includes=True).extract()
  pdb_interpretation_params.pdb_interpretation.sort_atoms=False
  if not pdb_str:
    pdb_inp = iotbx.pdb.input(file_name=file_name)
  else:
    pdb_inp = iotbx.pdb.input(lines=pdb_str, source_info=None)
  model = mmtbx.model.manager(
      model_input = pdb_inp,
      stop_for_unknowns = False,
      log=log)
  model.process(pdb_interpretation_params=pdb_interpretation_params,
                make_restraints = make_restraints)
  return model


def setup_riding(model):
  model.setup_riding_h_manager(ignore_h_with_dof = True)
  riding_h_manager = model.get_riding_h_manager()
  h_para = riding_h_manager.h_parameterization
  sel_h = model.get_hd_selection()
  sel_h_in_para = flex.bool(
    [bool(x) for x in riding_h_manager.h_parameterization])
  sel_h_not_in_para = sel_h_in_para.exclusive_or(sel_h)
  return(list(sel_h_not_in_para.iselection()))


def exercise_00(verbose):
  if (verbose): log = sys.stdout
  else: log = StringIO()
  model = get_model(file_name = None,
                    log       = log,
                    pdb_str   = pdb_str_1)
  rotatable_hd_sel = model.rotatable_hd_selection(use_shortcut = True)
  list_sel_h_not_in_para = setup_riding(model = model)
  s1 = set(list(rotatable_hd_sel))
  s2 = set(list_sel_h_not_in_para)
  assert(s1 == s2)
  ph = model.get_hierarchy()
  atoms = ph.atoms()


def exercise_01(verbose):
  if (verbose): log = sys.stdout
  else: log = StringIO()
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/ala_hd.pdb",
    test=os.path.isfile)
  model = get_model(pdb_file, log)
  # gives wrong result here: ignores D1,D2,D3
  #rotatable_hd_sel = model.rotatable_hd_selection(use_shortcut = True)
  list_sel_h_not_in_para = setup_riding(model = model)
  s1 = set([5, 6, 7, 9, 10, 11, 12, 13, 14])
  s2 = set(list_sel_h_not_in_para)
  assert(s1 == s2)


def exercise_02(verbose):
  if (verbose): log = sys.stdout
  else: log = StringIO()
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/ala_lys_arg_ser_tyr_neutron_hd.pdb",
    test=os.path.isfile)
  model = get_model(pdb_file, log)
  # here rotatable has 6,7,8 twice in the selection???
  rotatable_hd_sel = model.rotatable_hd_selection(use_shortcut = True)
  list_sel_h_not_in_para = setup_riding(model = model)
  s1 = set(list(rotatable_hd_sel))
  s2 = set(list_sel_h_not_in_para)
  assert(s1 == s2)


def exercise_03(verbose):
  """
  Check incomplete propellers in H/D exchanged sites
  """
  if (verbose): log = sys.stdout
  else: log = StringIO()
  model = get_model(file_name = None,
                    log       = log,
                    pdb_str   = pdb_str_2)
  rotatable_hd_sel = model.rotatable_hd_selection(use_shortcut = True)
  list_sel_h_not_in_para = setup_riding(model = model)
  #s1 = set(list(rotatable_hd_sel))
  s1 = set([12, 13, 14, 15, 18, 19, 20, 21, 31, 34, 41, 51])
  s2 = set(list_sel_h_not_in_para)
  #print(list(rotatable_hd_sel))
  #print(list_sel_h_not_in_para)
  assert(s1 == s2)
  ph = model.get_hierarchy()
  atoms = ph.atoms()

  #ph = model.get_hierarchy()
  #atoms = ph.atoms()
  #for bla in list_sel_h_not_in_para:
  #  print(atoms[bla].id_str())


# Ideal amino acids
pdb_str_1 = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
SCALE1      0.033333  0.000000  0.000000        0.00000
SCALE2      0.000000  0.033333  0.000000        0.00000
SCALE3      0.000000  0.000000  0.033333        0.00000
ATOM      1  N   ARG A   2       3.050  19.481  25.856  1.00  0.00           N
ATOM      2  CA  ARG A   2       3.759  20.265  24.852  1.00  0.00           C
ATOM      3  C   ARG A   2       5.253  20.304  25.155  1.00  0.00           C
ATOM      4  O   ARG A   2       5.659  20.516  26.298  1.00  0.00           O
ATOM      5  CB  ARG A   2       3.197  21.687  24.787  1.00  0.00           C
ATOM      6  CG  ARG A   2       3.800  22.542  23.684  1.00  0.00           C
ATOM      7  CD  ARG A   2       3.128  23.904  23.607  1.00  0.00           C
ATOM      8  NE  ARG A   2       3.721  24.751  22.576  1.00  0.00           N
ATOM      9  CZ  ARG A   2       3.311  25.983  22.290  1.00  0.00           C
ATOM     10  NH1 ARG A   2       2.301  26.525  22.957  1.00  0.00           N1+
ATOM     11  NH2 ARG A   2       3.914  26.676  21.334  1.00  0.00           N
ATOM     12  HA  ARG A   2       3.638  19.852  23.983  1.00  0.00           H
ATOM     13  HB2 ARG A   2       2.240  21.637  24.633  1.00  0.00           H
ATOM     14  HB3 ARG A   2       3.370  22.129  25.633  1.00  0.00           H
ATOM     15  HG2 ARG A   2       4.743  22.680  23.864  1.00  0.00           H
ATOM     16  HG3 ARG A   2       3.682  22.096  22.831  1.00  0.00           H
ATOM     17  HD2 ARG A   2       2.189  23.783  23.396  1.00  0.00           H
ATOM     18  HD3 ARG A   2       3.224  24.355  24.460  1.00  0.00           H
ATOM     19  HE  ARG A   2       4.380  24.432  22.124  1.00  0.00           H
ATOM     20 HH11 ARG A   2       1.905  26.082  23.579  1.00  0.00           H
ATOM     21 HH12 ARG A   2       2.041  27.322  22.767  1.00  0.00           H
ATOM     22 HH21 ARG A   2       4.569  26.329  20.898  1.00  0.00           H
ATOM     23 HH22 ARG A   2       3.649  27.473  21.149  1.00  0.00           H
ATOM     24  N   HIS A   3       6.067  20.098  24.123  1.00  0.00           N
ATOM     25  CA  HIS A   3       7.516  20.107  24.260  1.00  0.00           C
ATOM     26  C   HIS A   3       8.133  20.732  23.018  1.00  0.00           C
ATOM     27  O   HIS A   3       7.638  20.536  21.904  1.00  0.00           O
ATOM     28  CB  HIS A   3       8.069  18.691  24.464  1.00  0.00           C
ATOM     29  CG  HIS A   3       7.630  18.050  25.743  1.00  0.00           C
ATOM     30  ND1 HIS A   3       8.154  18.399  26.969  1.00  0.00           N
ATOM     31  CD2 HIS A   3       6.716  17.080  25.988  1.00  0.00           C
ATOM     32  CE1 HIS A   3       7.582  17.674  27.913  1.00  0.00           C
ATOM     33  NE2 HIS A   3       6.706  16.865  27.344  1.00  0.00           N
ATOM     34  H   HIS A   3       5.797  19.948  23.320  1.00  0.00           H
ATOM     35  HA  HIS A   3       7.764  20.646  25.027  1.00  0.00           H
ATOM     36  HB2 HIS A   3       7.768  18.129  23.732  1.00  0.00           H
ATOM     37  HB3 HIS A   3       9.038  18.731  24.470  1.00  0.00           H
ATOM     38  HD2 HIS A   3       6.194  16.642  25.355  1.00  0.00           H
ATOM     39  HE1 HIS A   3       7.765  17.723  28.824  1.00  0.00           H
ATOM     40  N   LYS A   4       9.213  21.483  23.217  1.00  0.00           N
ATOM     41  CA  LYS A   4       9.902  22.139  22.116  1.00  0.00           C
ATOM     42  C   LYS A   4      11.340  22.409  22.530  1.00  0.00           C
ATOM     43  O   LYS A   4      11.602  22.774  23.678  1.00  0.00           O
ATOM     44  CB  LYS A   4       9.205  23.447  21.722  1.00  0.00           C
ATOM     45  CG  LYS A   4       9.799  24.121  20.495  1.00  0.00           C
ATOM     46  CD  LYS A   4       9.008  25.357  20.102  1.00  0.00           C
ATOM     47  CE  LYS A   4       9.606  26.035  18.880  1.00  0.00           C
ATOM     48  NZ  LYS A   4       8.842  27.250  18.484  1.00  0.00           N1+
ATOM     49  H   LYS A   4       9.568  21.628  23.987  1.00  0.00           H
ATOM     50  HA  LYS A   4       9.909  21.551  21.344  1.00  0.00           H
ATOM     51  HB2 LYS A   4       8.272  23.258  21.533  1.00  0.00           H
ATOM     52  HB3 LYS A   4       9.270  24.070  22.462  1.00  0.00           H
ATOM     53  HG2 LYS A   4      10.710  24.392  20.688  1.00  0.00           H
ATOM     54  HG3 LYS A   4       9.784  23.500  19.750  1.00  0.00           H
ATOM     55  HD2 LYS A   4       8.096  25.101  19.892  1.00  0.00           H
ATOM     56  HD3 LYS A   4       9.018  25.991  20.837  1.00  0.00           H
ATOM     57  HE2 LYS A   4      10.517  26.302  19.078  1.00  0.00           H
ATOM     58  HE3 LYS A   4       9.594  25.415  18.135  1.00  0.00           H
ATOM     59  HZ1 LYS A   4       9.215  27.624  17.768  1.00  0.00           H
ATOM     60  HZ2 LYS A   4       8.001  27.031  18.292  1.00  0.00           H
ATOM     61  HZ3 LYS A   4       8.842  27.840  19.151  1.00  0.00           H
ATOM     62  N   ASP A   5      12.263  22.226  21.588  1.00  0.00           N
ATOM     63  CA  ASP A   5      13.678  22.447  21.845  1.00  0.00           C
ATOM     64  C   ASP A   5      14.374  22.763  20.531  1.00  0.00           C
ATOM     65  O   ASP A   5      14.035  22.197  19.489  1.00  0.00           O
ATOM     66  CB  ASP A   5      14.324  21.224  22.508  1.00  0.00           C
ATOM     67  CG  ASP A   5      15.729  21.504  23.005  1.00  0.00           C
ATOM     68  OD1 ASP A   5      16.213  22.640  22.819  1.00  0.00           O
ATOM     69  OD2 ASP A   5      16.350  20.587  23.582  1.00  0.00           O1-
ATOM     70  H   ASP A   5      12.091  21.972  20.785  1.00  0.00           H
ATOM     71  HA  ASP A   5      13.783  23.208  22.438  1.00  0.00           H
ATOM     72  HB2 ASP A   5      13.785  20.954  23.267  1.00  0.00           H
ATOM     73  HB3 ASP A   5      14.372  20.503  21.861  1.00  0.00           H
ATOM     74  N   GLU A   6      15.346  23.668  20.589  1.00  0.00           N
ATOM     75  CA  GLU A   6      16.096  24.064  19.403  1.00  0.00           C
ATOM     76  C   GLU A   6      17.433  24.684  19.794  1.00  0.00           C
ATOM     77  O   GLU A   6      17.829  24.647  20.959  1.00  0.00           O
ATOM     78  CB  GLU A   6      15.283  25.051  18.561  1.00  0.00           C
ATOM     79  CG  GLU A   6      15.942  25.433  17.245  1.00  0.00           C
ATOM     80  CD  GLU A   6      15.076  26.352  16.405  1.00  0.00           C
ATOM     81  OE1 GLU A   6      13.969  26.708  16.861  1.00  0.00           O
ATOM     82  OE2 GLU A   6      15.503  26.718  15.290  1.00  0.00           O1-
ATOM     83  H   GLU A   6      15.592  24.070  21.309  1.00  0.00           H
ATOM     84  HA  GLU A   6      16.273  23.279  18.862  1.00  0.00           H
ATOM     85  HB2 GLU A   6      14.424  24.651  18.356  1.00  0.00           H
ATOM     86  HB3 GLU A   6      15.153  25.864  19.073  1.00  0.00           H
ATOM     87  HG2 GLU A   6      16.776  25.893  17.428  1.00  0.00           H
ATOM     88  HG3 GLU A   6      16.112  24.628  16.730  1.00  0.00           H
TER
ATOM     89  N   SER B   2      10.846   4.462  24.534  1.00  0.00           N
ATOM     90  CA  SER B   2      11.518   5.747  24.391  1.00  0.00           C
ATOM     91  C   SER B   2      13.032   5.575  24.457  1.00  0.00           C
ATOM     92  O   SER B   2      13.532   4.595  25.009  1.00  0.00           O
ATOM     93  CB  SER B   2      11.052   6.719  25.477  1.00  0.00           C
ATOM     94  OG  SER B   2      11.455   6.279  26.762  1.00  0.00           O
ATOM     95  HA  SER B   2      11.295   6.129  23.528  1.00  0.00           H
ATOM     96  HB2 SER B   2      11.440   7.591  25.307  1.00  0.00           H
ATOM     97  HB3 SER B   2      10.084   6.777  25.453  1.00  0.00           H
ATOM     98  HG  SER B   2      11.125   5.523  26.924  1.00  0.00           H
ATOM     99  N   THR B   3      13.758   6.536  23.889  1.00  0.00           N
ATOM    100  CA  THR B   3      15.212   6.498  23.880  1.00  0.00           C
ATOM    101  C   THR B   3      15.741   7.920  23.778  1.00  0.00           C
ATOM    102  O   THR B   3      15.212   8.733  23.015  1.00  0.00           O
ATOM    103  CB  THR B   3      15.746   5.652  22.717  1.00  0.00           C
ATOM    104  OG1 THR B   3      15.230   4.319  22.813  1.00  0.00           O
ATOM    105  CG2 THR B   3      17.268   5.606  22.740  1.00  0.00           C
ATOM    106  H   THR B   3      13.425   7.226  23.499  1.00  0.00           H
ATOM    107  HA  THR B   3      15.530   6.111  24.711  1.00  0.00           H
ATOM    108  HB  THR B   3      15.465   6.046  21.876  1.00  0.00           H
ATOM    109  HG1 THR B   3      15.519   3.852  22.178  1.00  0.00           H
ATOM    110 HG21 THR B   3      17.587   4.933  22.118  1.00  0.00           H
ATOM    111 HG22 THR B   3      17.633   6.467  22.484  1.00  0.00           H
ATOM    112 HG23 THR B   3      17.580   5.381  23.631  1.00  0.00           H
ATOM    113  N   ASN B   4      16.785   8.212  24.552  1.00  0.00           N
ATOM    114  CA  ASN B   4      17.402   9.533  24.564  1.00  0.00           C
ATOM    115  C   ASN B   4      18.904   9.364  24.725  1.00  0.00           C
ATOM    116  O   ASN B   4      19.359   8.747  25.692  1.00  0.00           O
ATOM    117  CB  ASN B   4      16.835  10.400  25.693  1.00  0.00           C
ATOM    118  CG  ASN B   4      15.357  10.690  25.518  1.00  0.00           C
ATOM    119  OD1 ASN B   4      14.528  10.254  26.317  1.00  0.00           O
ATOM    120  ND2 ASN B   4      15.019  11.430  24.468  1.00  0.00           N
ATOM    121  H   ASN B   4      17.158   7.651  25.086  1.00  0.00           H
ATOM    122  HA  ASN B   4      17.232   9.979  23.719  1.00  0.00           H
ATOM    123  HB2 ASN B   4      16.953   9.937  26.537  1.00  0.00           H
ATOM    124  HB3 ASN B   4      17.307  11.247  25.709  1.00  0.00           H
ATOM    125 HD21 ASN B   4      14.193  11.622  24.326  1.00  0.00           H
ATOM    126 HD22 ASN B   4      15.626  11.717  23.931  1.00  0.00           H
ATOM    127  N   GLN B   5      19.666   9.908  23.781  1.00  0.00           N
ATOM    128  CA  GLN B   5      21.121   9.817  23.821  1.00  0.00           C
ATOM    129  C   GLN B   5      21.756  11.122  23.352  1.00  0.00           C
ATOM    130  O   GLN B   5      21.299  11.734  22.387  1.00  0.00           O
ATOM    131  CB  GLN B   5      21.610   8.654  22.957  1.00  0.00           C
ATOM    132  CG  GLN B   5      21.217   7.283  23.483  1.00  0.00           C
ATOM    133  CD  GLN B   5      21.709   6.155  22.598  1.00  0.00           C
ATOM    134  OE1 GLN B   5      22.296   6.391  21.542  1.00  0.00           O
ATOM    135  NE2 GLN B   5      21.472   4.920  23.025  1.00  0.00           N
ATOM    136  H   GLN B   5      19.362  10.339  23.102  1.00  0.00           H
ATOM    137  HA  GLN B   5      21.406   9.656  24.734  1.00  0.00           H
ATOM    138  HB2 GLN B   5      21.235   8.747  22.067  1.00  0.00           H
ATOM    139  HB3 GLN B   5      22.578   8.686  22.909  1.00  0.00           H
ATOM    140  HG2 GLN B   5      21.600   7.161  24.366  1.00  0.00           H
ATOM    141  HG3 GLN B   5      20.250   7.227  23.530  1.00  0.00           H
ATOM    142 HE21 GLN B   5      21.732   4.246  22.559  1.00  0.00           H
ATOM    143 HE22 GLN B   5      21.059   4.795  23.769  1.00  0.00           H
TER
ATOM    144  N   CYS C   2      14.789  22.141   7.021  1.00  0.00           N
ATOM    145  CA  CYS C   2      15.170  23.547   7.063  1.00  0.00           C
ATOM    146  C   CYS C   2      16.660  23.696   7.350  1.00  0.00           C
ATOM    147  O   CYS C   2      17.074  23.780   8.506  1.00  0.00           O
ATOM    148  CB  CYS C   2      14.354  24.292   8.122  1.00  0.00           C
ATOM    149  SG  CYS C   2      14.698  26.064   8.213  1.00  0.00           S
ATOM    150  HA  CYS C   2      14.987  23.952   6.201  1.00  0.00           H
ATOM    151  HB2 CYS C   2      13.411  24.185   7.920  1.00  0.00           H
ATOM    152  HB3 CYS C   2      14.548  23.908   8.991  1.00  0.00           H
ATOM    153  HG  CYS C   2      14.012  26.553   9.067  1.00  0.00           H
ATOM    154  N   GLY C   3      17.463  23.727   6.287  1.00  0.00           N
ATOM    155  CA  GLY C   3      18.892  23.864   6.428  1.00  0.00           C
ATOM    156  C   GLY C   3      19.301  25.288   6.737  1.00  0.00           C
ATOM    157  O   GLY C   3      18.476  26.206   6.840  1.00  0.00           O
ATOM    158  H   GLY C   3      17.194  23.670   5.472  1.00  0.00           H
ATOM    159  HA2 GLY C   3      19.202  23.291   7.146  1.00  0.00           H
ATOM    160  HA3 GLY C   3      19.325  23.590   5.604  1.00  0.00           H
ATOM    161  N   PRO C   4      20.618  25.490   6.891  1.00  0.00           N
ATOM    162  CA  PRO C   4      21.173  26.813   7.198  1.00  0.00           C
ATOM    163  C   PRO C   4      21.113  27.768   6.010  1.00  0.00           C
ATOM    164  O   PRO C   4      21.345  27.333   4.882  1.00  0.00           O
ATOM    165  CB  PRO C   4      22.622  26.503   7.578  1.00  0.00           C
ATOM    166  CG  PRO C   4      22.943  25.252   6.843  1.00  0.00           C
ATOM    167  CD  PRO C   4      21.667  24.460   6.790  1.00  0.00           C
ATOM    168  HA  PRO C   4      20.715  27.207   7.957  1.00  0.00           H
ATOM    169  HB2 PRO C   4      23.199  27.228   7.291  1.00  0.00           H
ATOM    170  HB3 PRO C   4      22.689  26.366   8.536  1.00  0.00           H
ATOM    171  HG2 PRO C   4      23.244  25.471   5.947  1.00  0.00           H
ATOM    172  HG3 PRO C   4      23.628  24.761   7.322  1.00  0.00           H
ATOM    173  HD2 PRO C   4      21.599  23.988   5.945  1.00  0.00           H
ATOM    174  HD3 PRO C   4      21.617  23.851   7.543  1.00  0.00           H
TER
ATOM    175  N   ALA D   2       2.619   6.211   7.571  1.00  0.00           N
ATOM    176  CA  ALA D   2       3.623   7.047   8.219  1.00  0.00           C
ATOM    177  C   ALA D   2       4.961   6.321   8.294  1.00  0.00           C
ATOM    178  O   ALA D   2       5.056   5.233   8.863  1.00  0.00           O
ATOM    179  CB  ALA D   2       3.159   7.450   9.610  1.00  0.00           C
ATOM    180  HA  ALA D   2       3.748   7.855   7.697  1.00  0.00           H
ATOM    181  HB1 ALA D   2       3.840   8.005  10.021  1.00  0.00           H
ATOM    182  HB2 ALA D   2       2.329   7.947   9.534  1.00  0.00           H
ATOM    183  HB3 ALA D   2       3.018   6.650  10.139  1.00  0.00           H
ATOM    184  N   LEU D   3       6.002   6.934   7.712  1.00  0.00           N
ATOM    185  CA  LEU D   3       7.348   6.371   7.695  1.00  0.00           C
ATOM    186  C   LEU D   3       8.344   7.510   7.935  1.00  0.00           C
ATOM    187  O   LEU D   3       9.150   7.868   7.078  1.00  0.00           O
ATOM    188  CB  LEU D   3       7.634   5.644   6.379  1.00  0.00           C
ATOM    189  CG  LEU D   3       6.792   4.396   6.105  1.00  0.00           C
ATOM    190  CD1 LEU D   3       7.036   3.888   4.692  1.00  0.00           C
ATOM    191  CD2 LEU D   3       7.088   3.307   7.125  1.00  0.00           C
ATOM    192  H   LEU D   3       5.945   7.694   7.313  1.00  0.00           H
ATOM    193  HA  LEU D   3       7.438   5.734   8.421  1.00  0.00           H
ATOM    194  HB2 LEU D   3       7.479   6.264   5.649  1.00  0.00           H
ATOM    195  HB3 LEU D   3       8.565   5.371   6.377  1.00  0.00           H
ATOM    196  HG  LEU D   3       5.853   4.628   6.180  1.00  0.00           H
ATOM    197 HD11 LEU D   3       6.493   3.099   4.543  1.00  0.00           H
ATOM    198 HD12 LEU D   3       6.791   4.583   4.061  1.00  0.00           H
ATOM    199 HD13 LEU D   3       7.975   3.668   4.594  1.00  0.00           H
ATOM    200 HD21 LEU D   3       6.619   2.498   6.867  1.00  0.00           H
ATOM    201 HD22 LEU D   3       8.044   3.144   7.143  1.00  0.00           H
ATOM    202 HD23 LEU D   3       6.785   3.600   7.998  1.00  0.00           H
ATOM    203  N   ILE D   4       8.286   8.091   9.130  1.00  0.00           N
ATOM    204  CA  ILE D   4       9.174   9.189   9.499  1.00  0.00           C
ATOM    205  C   ILE D   4      10.531   8.595   9.863  1.00  0.00           C
ATOM    206  O   ILE D   4      10.691   7.990  10.926  1.00  0.00           O
ATOM    207  CB  ILE D   4       8.599  10.013  10.657  1.00  0.00           C
ATOM    208  CG1 ILE D   4       7.284  10.674  10.236  1.00  0.00           C
ATOM    209  CG2 ILE D   4       9.604  11.067  11.108  1.00  0.00           C
ATOM    210  CD1 ILE D   4       6.474  11.224  11.392  1.00  0.00           C
ATOM    211  H   ILE D   4       7.736   7.865   9.751  1.00  0.00           H
ATOM    212  HA  ILE D   4       9.293   9.777   8.736  1.00  0.00           H
ATOM    213  HB  ILE D   4       8.421   9.417  11.401  1.00  0.00           H
ATOM    214 HG12 ILE D   4       7.482  11.410   9.637  1.00  0.00           H
ATOM    215 HG13 ILE D   4       6.734  10.018   9.778  1.00  0.00           H
ATOM    216 HG21 ILE D   4       9.158  11.702  11.690  1.00  0.00           H
ATOM    217 HG22 ILE D   4      10.326  10.633  11.590  1.00  0.00           H
ATOM    218 HG23 ILE D   4       9.953  11.525  10.328  1.00  0.00           H
ATOM    219 HD11 ILE D   4       5.603  11.498  11.065  1.00  0.00           H
ATOM    220 HD12 ILE D   4       6.373  10.531  12.063  1.00  0.00           H
ATOM    221 HD13 ILE D   4       6.939  11.986  11.771  1.00  0.00           H
ATOM    222  N   MET D   5      11.514   8.768   8.978  1.00  0.00           N
ATOM    223  CA  MET D   5      12.863   8.258   9.192  1.00  0.00           C
ATOM    224  C   MET D   5      13.872   9.277   8.657  1.00  0.00           C
ATOM    225  O   MET D   5      14.680   9.000   7.775  1.00  0.00           O
ATOM    226  CB  MET D   5      13.042   6.896   8.525  1.00  0.00           C
ATOM    227  CG  MET D   5      12.226   5.783   9.162  1.00  0.00           C
ATOM    228  SD  MET D   5      12.466   4.191   8.350  1.00  0.00           S
ATOM    229  CE  MET D   5      11.320   3.166   9.268  1.00  0.00           C
ATOM    230  H   MET D   5      11.418   9.186   8.232  1.00  0.00           H
ATOM    231  HA  MET D   5      13.009   8.144  10.144  1.00  0.00           H
ATOM    232  HB2 MET D   5      12.772   6.965   7.596  1.00  0.00           H
ATOM    233  HB3 MET D   5      13.977   6.644   8.577  1.00  0.00           H
ATOM    234  HG2 MET D   5      12.490   5.689  10.091  1.00  0.00           H
ATOM    235  HG3 MET D   5      11.284   6.009   9.107  1.00  0.00           H
ATOM    236  HE1 MET D   5      11.356   2.262   8.918  1.00  0.00           H
ATOM    237  HE2 MET D   5      11.574   3.169  10.204  1.00  0.00           H
ATOM    238  HE3 MET D   5      10.425   3.525   9.167  1.00  0.00           H
ATOM    239  N   PHE D   6      13.825  10.488   9.205  1.00  0.00           N
ATOM    240  CA  PHE D   6      14.731  11.547   8.787  1.00  0.00           C
ATOM    241  C   PHE D   6      16.111  11.331   9.393  1.00  0.00           C
ATOM    242  O   PHE D   6      16.244  10.831  10.514  1.00  0.00           O
ATOM    243  CB  PHE D   6      14.186  12.914   9.206  1.00  0.00           C
ATOM    244  CG  PHE D   6      12.873  13.276   8.562  1.00  0.00           C
ATOM    245  CD1 PHE D   6      12.416  12.605   7.437  1.00  0.00           C
ATOM    246  CD2 PHE D   6      12.093  14.294   9.087  1.00  0.00           C
ATOM    247  CE1 PHE D   6      11.210  12.943   6.852  1.00  0.00           C
ATOM    248  CE2 PHE D   6      10.887  14.635   8.506  1.00  0.00           C
ATOM    249  CZ  PHE D   6      10.445  13.959   7.388  1.00  0.00           C
ATOM    250  H   PHE D   6      13.274  10.720   9.823  1.00  0.00           H
ATOM    251  HA  PHE D   6      14.827  11.530   7.822  1.00  0.00           H
ATOM    252  HB2 PHE D   6      14.055  12.916  10.167  1.00  0.00           H
ATOM    253  HB3 PHE D   6      14.833  13.595   8.962  1.00  0.00           H
ATOM    254  HD1 PHE D   6      12.924  11.919   7.069  1.00  0.00           H
ATOM    255  HD2 PHE D   6      12.386  14.753   9.841  1.00  0.00           H
ATOM    256  HE1 PHE D   6      10.914  12.486   6.099  1.00  0.00           H
ATOM    257  HE2 PHE D   6      10.373  15.320   8.869  1.00  0.00           H
ATOM    258  HZ  PHE D   6       9.634  14.188   6.995  1.00  0.00           H
ATOM    259  N   TRP D   7      17.142  11.711   8.641  1.00  0.00           N
ATOM    260  CA  TRP D   7      18.523  11.570   9.084  1.00  0.00           C
ATOM    261  C   TRP D   7      19.334  12.742   8.558  1.00  0.00           C
ATOM    262  O   TRP D   7      19.223  13.098   7.381  1.00  0.00           O
ATOM    263  CB  TRP D   7      19.129  10.247   8.600  1.00  0.00           C
ATOM    264  CG  TRP D   7      18.518   9.038   9.242  1.00  0.00           C
ATOM    265  CD1 TRP D   7      17.534   8.246   8.727  1.00  0.00           C
ATOM    266  CD2 TRP D   7      18.852   8.483  10.520  1.00  0.00           C
ATOM    267  NE1 TRP D   7      17.234   7.232   9.605  1.00  0.00           N
ATOM    268  CE2 TRP D   7      18.030   7.356  10.713  1.00  0.00           C
ATOM    269  CE3 TRP D   7      19.767   8.831  11.519  1.00  0.00           C
ATOM    270  CZ2 TRP D   7      18.095   6.573  11.864  1.00  0.00           C
ATOM    271  CZ3 TRP D   7      19.829   8.053  12.661  1.00  0.00           C
ATOM    272  CH2 TRP D   7      18.998   6.937  12.824  1.00  0.00           C
ATOM    273  H   TRP D   7      17.064  12.059   7.859  1.00  0.00           H
ATOM    274  HA  TRP D   7      18.555  11.585  10.053  1.00  0.00           H
ATOM    275  HB2 TRP D   7      18.996  10.173   7.642  1.00  0.00           H
ATOM    276  HB3 TRP D   7      20.077  10.244   8.803  1.00  0.00           H
ATOM    277  HD1 TRP D   7      17.125   8.373   7.901  1.00  0.00           H
ATOM    278  HE1 TRP D   7      16.644   6.620   9.479  1.00  0.00           H
ATOM    279  HE3 TRP D   7      20.322   9.570  11.418  1.00  0.00           H
ATOM    280  HZ2 TRP D   7      17.544   5.832  11.976  1.00  0.00           H
ATOM    281  HZ3 TRP D   7      20.434   8.274  13.332  1.00  0.00           H
ATOM    282  HH2 TRP D   7      19.063   6.433  13.603  1.00  0.00           H
ATOM    283  N   TYR D   8      20.146  13.338   9.431  1.00  0.00           N
ATOM    284  CA  TYR D   8      20.988  14.478   9.066  1.00  0.00           C
ATOM    285  C   TYR D   8      22.308  14.337   9.818  1.00  0.00           C
ATOM    286  O   TYR D   8      22.402  14.703  10.993  1.00  0.00           O
ATOM    287  CB  TYR D   8      20.302  15.800   9.392  1.00  0.00           C
ATOM    288  CG  TYR D   8      19.076  16.075   8.551  1.00  0.00           C
ATOM    289  CD1 TYR D   8      19.192  16.579   7.263  1.00  0.00           C
ATOM    290  CD2 TYR D   8      17.801  15.832   9.046  1.00  0.00           C
ATOM    291  CE1 TYR D   8      18.075  16.832   6.490  1.00  0.00           C
ATOM    292  CE2 TYR D   8      16.678  16.082   8.282  1.00  0.00           C
ATOM    293  CZ  TYR D   8      16.820  16.582   7.005  1.00  0.00           C
ATOM    294  OH  TYR D   8      15.704  16.833   6.240  1.00  0.00           O
ATOM    295  H   TYR D   8      20.228  13.097  10.252  1.00  0.00           H
ATOM    296  HA  TYR D   8      21.173  14.454   8.115  1.00  0.00           H
ATOM    297  HB2 TYR D   8      20.027  15.788  10.322  1.00  0.00           H
ATOM    298  HB3 TYR D   8      20.931  16.523   9.244  1.00  0.00           H
ATOM    299  HD1 TYR D   8      20.037  16.749   6.913  1.00  0.00           H
ATOM    300  HD2 TYR D   8      17.703  15.494   9.907  1.00  0.00           H
ATOM    301  HE1 TYR D   8      18.168  17.170   5.629  1.00  0.00           H
ATOM    302  HE2 TYR D   8      15.831  15.914   8.626  1.00  0.00           H
ATOM    303  HH  TYR D   8      15.010  16.638   6.671  1.00  0.00           H
ATOM    304  N   VAL D   9      23.320  13.807   9.139  1.00  0.00           N
ATOM    305  CA  VAL D   9      24.632  13.618   9.743  1.00  0.00           C
ATOM    306  C   VAL D   9      25.377  14.947   9.786  1.00  0.00           C
ATOM    307  O   VAL D   9      25.756  15.491   8.749  1.00  0.00           O
ATOM    308  CB  VAL D   9      25.449  12.559   8.980  1.00  0.00           C
ATOM    309  CG1 VAL D   9      26.843  12.427   9.578  1.00  0.00           C
ATOM    310  CG2 VAL D   9      24.726  11.220   8.995  1.00  0.00           C
ATOM    311  H   VAL D   9      23.271  13.548   8.320  1.00  0.00           H
ATOM    312  HA  VAL D   9      24.518  13.309  10.655  1.00  0.00           H
ATOM    313  HB  VAL D   9      25.544  12.838   8.056  1.00  0.00           H
ATOM    314 HG11 VAL D   9      27.251  11.612   9.245  1.00  0.00           H
ATOM    315 HG12 VAL D   9      27.376  13.193   9.316  1.00  0.00           H
ATOM    316 HG13 VAL D   9      26.769  12.386  10.544  1.00  0.00           H
ATOM    317 HG21 VAL D   9      25.275  10.560   8.542  1.00  0.00           H
ATOM    318 HG22 VAL D   9      24.578  10.952   9.915  1.00  0.00           H
ATOM    319 HG23 VAL D   9      23.876  11.314   8.536  1.00  0.00           H
TER
END
"""

pdb_str_2 = """
CRYST1   43.705   48.284   96.292  90.00  90.00  90.00 P 1
SCALE1      0.022881  0.000000  0.000000        0.00000
SCALE2      0.000000  0.020711  0.000000        0.00000
SCALE3      0.000000  0.000000  0.010385        0.00000
ATOM     37  N   LEU A   3      23.374 -21.865 -28.222  1.00 15.90           N
ATOM     38  CA  LEU A   3      23.342 -21.224 -29.530  1.00 15.47           C
ATOM     39  C   LEU A   3      23.637 -19.727 -29.491  1.00 18.18           C
ATOM     40  O   LEU A   3      24.174 -19.175 -30.462  1.00 20.24           O
ATOM     41  CB  LEU A   3      22.021 -21.521 -30.238  1.00 15.76           C
ATOM     42  CG  LEU A   3      21.790 -23.005 -30.533  1.00 19.30           C
ATOM     43  CD1 LEU A   3      20.460 -23.197 -31.290  1.00 20.74           C
ATOM     44  CD2 LEU A   3      22.957 -23.627 -31.336  1.00 17.02           C
ATOM     45  DA  LEU A   3      24.132 -21.653 -30.101  1.00 12.76           D
ATOM     46  DB3 LEU A   3      21.233 -21.196 -29.610  1.00 17.93           D
ATOM     48  DB2ALEU A   3      21.964 -20.976 -31.155  0.67 16.51           D
ATOM     49  DG ALEU A   3      21.719 -23.530 -29.581  0.68 18.31           D
REMARK ATOM     50 DD11ALEU A   3      19.620 -22.925 -30.652  0.74 20.52           D
ATOM     51 DD12ALEU A   3      20.364 -24.239 -31.597  0.77 19.78           D
ATOM     52 DD13ALEU A   3      20.458 -22.568 -32.171  0.49 20.26           D
ATOM     53 DD21ALEU A   3      23.466 -22.855 -31.890  0.40 16.87           D
REMARK ATOM     54 DD22ALEU A   3      22.563 -24.357 -32.028  0.00 16.68           D
ATOM     55 DD23ALEU A   3      23.653 -24.109 -30.673  0.89 14.39           D
ATOM     57  HB2BLEU A   3      21.964 -20.976 -31.155  0.33 14.51           H
ATOM     58  HG BLEU A   3      21.719 -23.530 -29.581  0.32 18.31           H
ATOM     59 HD11BLEU A   3      19.620 -22.925 -30.652  0.26 20.52           H
ATOM     60 HD12BLEU A   3      20.364 -24.239 -31.597  0.23 19.78           H
REMARK ATOM     61 HD13BLEU A   3      20.458 -22.568 -32.171  0.51 20.26           H
REMARK ATOM     62 HD21BLEU A   3      23.466 -22.855 -31.890  0.60 16.87           H
ATOM     63 HD22BLEU A   3      22.563 -24.357 -32.028  0.58 16.68           H
ATOM     64 HD23BLEU A   3      23.653 -24.109 -30.673  0.11 14.39           H
ATOM     65  N   SER A   4       9.243 -32.509 -36.471  1.00 38.02           N
ATOM     66  CA  SER A   4      10.063 -33.604 -36.981  1.00 39.55           C
ATOM     67  C   SER A   4      11.459 -33.153 -37.367  1.00 37.45           C
ATOM     68  O   SER A   4      11.633 -32.304 -38.239  1.00 38.82           O
ATOM     69  CB  SER A   4       9.384 -34.279 -38.180  1.00 41.66           C
ATOM     70  OG  SER A   4       8.375 -35.181 -37.753  1.00 44.91           O
ATOM     71  DA  SER A   4      10.165 -34.324 -36.187  1.00 38.72           D
ATOM     73  DB2ASER A   4      10.116 -34.818 -38.774  0.63 41.36           D
ATOM     74  DB3ASER A   4       8.918 -33.518 -38.793  0.36 41.73           D
ATOM     75  DG ASER A   4       8.699 -36.084 -37.842  0.36 43.61           D
ATOM     77  HB2BSER A   4       8.918 -33.518 -38.793  0.64 41.73           H
ATOM     78  HB3BSER A   4      10.116 -34.818 -38.774  0.37 41.36           H
ATOM     79  HG BSER A   4       8.699 -36.084 -37.842  0.64 43.61           H
ATOM     80  N   SER A   5      26.478 -19.398  45.776  1.00 19.78           N
ATOM     81  C   SER A   5      27.940 -19.562  43.818  1.00 21.99           C
ATOM     82  O   SER A   5      29.005 -19.770  43.306  1.00 21.20           O
ATOM     83  CA ASER A   5      27.849 -19.441  45.318  0.64 18.99           C
ATOM     84  CB ASER A   5      28.611 -20.590  45.987  0.64 21.66           C
ATOM     85  OG ASER A   5      28.819 -20.337  47.363  0.64 23.21           O
ATOM     86  DG ASER A   5      29.639 -20.155  47.499  0.64 22.10           D
ATOM     88  HA ASER A   5      28.296 -18.610  45.561  0.64 20.89           H
ATOM     89  HB2ASER A   5      28.094 -21.408  45.903  0.00 22.17           H
ATOM     90  HB3ASER A   5      29.468 -20.686  45.544  0.64 22.38           H
ATOM     91  CA BSER A   5      27.851 -19.476  45.322  0.36 19.02           C
ATOM     92  CB BSER A   5      28.561 -20.680  45.932  0.36 21.98           C
ATOM     93  OG BSER A   5      28.005 -21.879  45.425  0.36 25.40           O
ATOM     95  HA BSER A   5      28.327 -18.675  45.599  0.36 21.03           H
ATOM     96  HB2BSER A   5      29.498 -20.638  45.691  0.36 22.55           H
ATOM     97  HB3BSER A   5      28.462 -20.665  46.898  0.36 23.83           H
ATOM     98  HG BSER A   5      27.486 -22.196  46.000  0.36 22.10           H
TER
"""

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
