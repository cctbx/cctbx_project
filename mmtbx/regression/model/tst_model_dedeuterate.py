from __future__ import absolute_import, division, print_function
import time, sys, os
from six.moves import cStringIO as StringIO
import libtbx.load_env
import mmtbx.model
import iotbx.pdb
from mmtbx.monomer_library import pdb_interpretation

def run():
  verbose = "--verbose" in sys.argv[1:]
  exercise_00(verbose=verbose)
  exercise_01(verbose=verbose)
  exercise_02(verbose=verbose)
  exercise_03(verbose=verbose)
  exercise_04(verbose=verbose)
  exercise_05(verbose=verbose)

def get_model(file_name, log, pdb_str=None):
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
  model.process(pdb_interpretation_params=pdb_interpretation_params)
  return model

def exercise_00(verbose):
  if (verbose): log = sys.stdout
  else: log = StringIO()
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/occ_mix1.pdb",
    test=os.path.isfile)
  model = get_model(pdb_file, log)
  model.de_deuterate()
  assert(len(model.hd_group_selections())==0)
  ph = model.get_hierarchy()
  atoms = ph.atoms()
  assert(atoms.extract_element().count('D') == 0)
  assert(model.get_hd_selection().count(True) == 14)
  atom_names = [x.strip() for x in atoms.extract_name()]
  assert ('DA' not in atom_names)

def exercise_01(verbose):
  if (verbose): log = sys.stdout
  else: log = StringIO()
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/ala_hd.pdb",
    test=os.path.isfile)
  model = get_model(pdb_file, log)
  model.de_deuterate()
  assert(len(model.hd_group_selections())==0)
  ph = model.get_hierarchy()
  atoms = ph.atoms()
  assert(atoms.extract_element().count('D') == 0)
  assert(model.get_hd_selection().count(True) == 7)

def exercise_02(verbose):
  if (verbose): log = sys.stdout
  else: log = StringIO()
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/ala_lys_arg_ser_tyr_neutron_hd.pdb",
    test=os.path.isfile)
  model = get_model(pdb_file, log)
  model.de_deuterate()
  assert(len(model.hd_group_selections())==0)
  ph = model.get_hierarchy()
  atoms = ph.atoms()
  assert(atoms.extract_element().count('D') == 0)
  assert(model.get_hd_selection().count(True) == 47)

def exercise_03(verbose):
  if (verbose): log = sys.stdout
  else: log = StringIO()
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/NAD_594_HD.pdb",
    test=os.path.isfile)
  model = get_model(pdb_file, log)
  model.de_deuterate()
  assert(len(model.hd_group_selections())==0)
  ph = model.get_hierarchy()
  atoms = ph.atoms()
  assert(atoms.extract_element().count('D') == 0)
  assert(model.get_hd_selection().count(True) == 43)

def exercise_04(verbose):
  if (verbose): log = sys.stdout
  else: log = StringIO()
  model = get_model(file_name=None, log = log, pdb_str = pdb_str_1)
  model.de_deuterate()
  assert(len(model.hd_group_selections())==0)
  ph = model.get_hierarchy()
  atoms = ph.atoms()
  assert(atoms.extract_element().count('D') == 0)
  assert(model.get_hd_selection().count(True) == 15)
  atom_names = [x.strip() for x in atoms.extract_name()]
  assert ('D' not in atom_names)
  assert ('DG' not in atom_names)

def exercise_05(verbose):
  if (verbose): log = sys.stdout
  else: log = StringIO()
  model = get_model(file_name=None, log = log, pdb_str = pdb_str_2)
  #model.get_hierarchy().write_pdb_file(file_name="haha.pdb")
  model.de_deuterate()
  assert(len(model.hd_group_selections())==0)
  ph = model.get_hierarchy()
  atoms = ph.atoms()
  assert(atoms.extract_element().count('D') == 0)
  #print(model.get_hd_selection().count(True))
  assert(model.get_hd_selection().count(True) == 44)
  assert(model.selection('resname DOD').count(True) == 0)
  #ph.write_pdb_file(file_name="test.pdb")

pdb_str_1 = '''
CRYST1   43.705   48.284   96.292  90.00  90.00  90.00 P 1
SCALE1      0.022881  0.000000  0.000000        0.00000
SCALE2      0.000000  0.020711  0.000000        0.00000
SCALE3      0.000000  0.000000  0.010385        0.00000
ATOM     65  N   SER A   4       9.243 -32.509 -36.471  1.00 38.02           N
ATOM     66  CA  SER A   4      10.063 -33.604 -36.981  1.00 39.55           C
ATOM     67  C   SER A   4      11.459 -33.153 -37.367  1.00 37.45           C
ATOM     68  O   SER A   4      11.633 -32.304 -38.239  1.00 38.82           O
ATOM     69  CB  SER A   4       9.384 -34.279 -38.180  1.00 41.66           C
ATOM     70  OG  SER A   4       8.375 -35.181 -37.753  1.00 44.91           O
ATOM     71  DA  SER A   4      10.165 -34.324 -36.187  1.00 38.72           D
ATOM     72  D  ASER A   4       8.730 -31.950 -37.103  0.62 38.48           D
ATOM     73  DB2ASER A   4      10.116 -34.818 -38.774  0.63 41.36           D
ATOM     74  DB3ASER A   4       8.918 -33.518 -38.793  0.36 41.73           D
ATOM     75  DG ASER A   4       8.699 -36.084 -37.842  0.36 43.61           D
ATOM     76  H  BSER A   4       8.730 -31.950 -37.103  0.28 38.48           H
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
ATOM     87  H  ASER A   5      26.064 -20.150  45.742  0.06 18.55           H
ATOM     88  HA ASER A   5      28.296 -18.610  45.561  0.64 20.89           H
ATOM     89  HB2ASER A   5      28.094 -21.408  45.903  0.00 22.17           H
ATOM     90  HB3ASER A   5      29.468 -20.686  45.544  0.64 22.38           H
ATOM     91  CA BSER A   5      27.851 -19.476  45.322  0.36 19.02           C
ATOM     92  CB BSER A   5      28.561 -20.680  45.932  0.36 21.98           C
ATOM     93  OG BSER A   5      28.005 -21.879  45.425  0.36 25.40           O
ATOM     94  D  BSER A   5      26.023 -20.128  45.688  0.94 18.55           D
ATOM     95  HA BSER A   5      28.327 -18.675  45.599  0.36 21.03           H
ATOM     96  HB2BSER A   5      29.498 -20.638  45.691  0.36 22.55           H
ATOM     97  HB3BSER A   5      28.462 -20.665  46.898  0.36 23.83           H
ATOM     98  HG BSER A   5      27.486 -22.196  46.000  0.36 22.10           H
TER
'''

pdb_str_2 = '''
CRYST1   43.705   48.284   96.292  90.00  90.00  90.00 P 1
SCALE1      0.022881  0.000000  0.000000        0.00000
SCALE2      0.000000  0.020711  0.000000        0.00000
SCALE3      0.000000  0.000000  0.010385        0.00000
ATOM      1  N   ASN A   1      25.174 -24.761 -21.940  1.00 16.59           N
ATOM      2  CA  ASN A   1      25.971 -25.578 -22.850  1.00 17.25           C
ATOM      3  C   ASN A   1      26.037 -24.941 -24.225  1.00 16.66           C
ATOM      4  O   ASN A   1      27.077 -24.979 -24.874  1.00 18.42           O
ATOM      5  CB  ASN A   1      25.498 -27.040 -22.913  1.00 18.30           C
ATOM      6  CG  ASN A   1      24.187 -27.208 -23.672  1.00 20.42           C
ATOM      7  OD1 ASN A   1      23.167 -26.590 -23.335  1.00 20.87           O
ATOM      8  ND2 ASN A   1      24.200 -28.069 -24.685  1.00 18.15           N
ATOM      9  DA  ASN A   1      26.975 -25.581 -22.461  1.00 15.99           D
ATOM     10 DD21 ASN A   1      25.030 -28.545 -24.882  1.00 16.07           D
ATOM     11  D  AASN A   1      24.379 -25.105 -21.517  0.82 16.71           D
ATOM     12  DB2AASN A   1      25.358 -27.400 -21.898  0.65 18.82           D
ATOM     13  DB3AASN A   1      26.255 -27.641 -23.402  0.56 18.30           D
ATOM     14 DD22AASN A   1      23.370 -28.163 -25.192  0.72 16.98           D
ATOM     15  H  BASN A   1      24.379 -25.105 -21.517  0.18 16.71           H
ATOM     16  HB2BASN A   1      26.255 -27.641 -23.402  0.44 18.30           H
ATOM     17  HB3BASN A   1      25.358 -27.400 -21.898  0.35 18.82           H
ATOM     18 HD22BASN A   1      23.370 -28.163 -25.192  0.28 16.98           H
ATOM     19  N   ASN A   2      23.615 -22.281 -25.492  1.00 16.59           N
ATOM     20  CA  ASN A   2      24.412 -23.098 -26.402  1.00 17.25           C
ATOM     21  C   ASN A   2      24.478 -22.461 -27.777  1.00 16.66           C
ATOM     22  O   ASN A   2      25.518 -22.499 -28.426  1.00 18.42           O
ATOM     23  CB  ASN A   2      23.939 -24.560 -26.465  1.00 18.30           C
ATOM     24  CG  ASN A   2      22.628 -24.728 -27.224  1.00 20.42           C
ATOM     25  OD1 ASN A   2      21.608 -24.110 -26.887  1.00 20.87           O
ATOM     26  ND2 ASN A   2      22.641 -25.589 -28.237  1.00 18.15           N
ATOM     27  DA  ASN A   2      25.416 -23.101 -26.013  1.00 15.99           D
ATOM     28 HD21 ASN A   2      23.471 -26.065 -28.434  1.00 16.07           H
ATOM     29  D  AASN A   2      22.820 -22.625 -25.069  0.82 16.71           D
ATOM     30  DB2AASN A   2      23.799 -24.920 -25.450  0.65 18.82           D
ATOM     31  DB3AASN A   2      24.696 -25.161 -26.954  0.56 18.30           D
ATOM     32 DD22AASN A   2      21.811 -25.683 -28.744  0.72 16.98           D
ATOM     33  H  BASN A   2      22.820 -22.625 -25.069  0.18 16.71           H
ATOM     34  HB2BASN A   2      24.696 -25.161 -26.954  0.44 18.30           H
ATOM     35  HB3BASN A   2      23.799 -24.920 -25.450  0.35 18.82           H
ATOM     36 HD22BASN A   2      21.811 -25.683 -28.744  0.28 16.98           H
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
ATOM     47  D  ALEU A   3      22.567 -21.865 -27.663  0.91 16.63           D
ATOM     48  DB2ALEU A   3      21.964 -20.976 -31.155  0.67 16.51           D
ATOM     49  DG ALEU A   3      21.719 -23.530 -29.581  0.68 18.31           D
ATOM     50 DD11ALEU A   3      19.620 -22.925 -30.652  0.74 20.52           D
ATOM     51 DD12ALEU A   3      20.364 -24.239 -31.597  0.77 19.78           D
ATOM     52 DD13ALEU A   3      20.458 -22.568 -32.171  0.49 20.26           D
ATOM     53 DD21ALEU A   3      23.466 -22.855 -31.890  0.40 16.87           D
ATOM     54 DD22ALEU A   3      22.563 -24.357 -32.028  0.00 16.68           D
ATOM     55 DD23ALEU A   3      23.653 -24.109 -30.673  0.89 14.39           D
ATOM     56  H  BLEU A   3      22.567 -21.865 -27.663  0.09 16.63           H
ATOM     57  HB2BLEU A   3      21.964 -20.976 -31.155  0.33 14.51           H
ATOM     58  HG BLEU A   3      21.719 -23.530 -29.581  0.32 18.31           H
ATOM     59 HD11BLEU A   3      19.620 -22.925 -30.652  0.26 20.52           H
ATOM     60 HD12BLEU A   3      20.364 -24.239 -31.597  0.23 19.78           H
ATOM     61 HD13BLEU A   3      20.458 -22.568 -32.171  0.51 20.26           H
ATOM     62 HD21BLEU A   3      23.466 -22.855 -31.890  0.60 16.87           H
ATOM     63 HD22BLEU A   3      22.563 -24.357 -32.028  0.58 16.68           H
ATOM     64 HD23BLEU A   3      23.653 -24.109 -30.673  0.11 14.39           H
TER
HETATM  100  O   DOD B   1      25.202 -21.329 -23.306  1.00 30.00           O
HETATM  101  D1  DOD B   1      25.963 -20.763 -23.255  1.00 30.00           D
HETATM  102  D2  DOD B   1      24.439 -20.764 -23.305  1.00 30.00           D
HETATM  103  O   HOH B   2      26.442 -21.454 -25.349  1.00 30.00           O
HETATM  104  D1  HOH B   2      27.249 -20.953 -25.339  1.00 30.00           D
HETATM  105  O   HOH B   3      24.970 -18.184 -24.460  1.00 30.00           O
HETATM  106  O   HOH B   4      21.573 -23.174 -24.023  1.00 30.00           O
HETATM  107  D2  HOH B   4      20.731 -22.736 -24.071  1.00 30.00           D
HETATM  108  D1 AHOH B   4      22.045 -22.784 -23.297  1.00 30.00           D
HETATM  109  D1 BHOH B   4      22.145 -21.784 -22.297  1.00 30.00           D
HETATM  110  D1  HOH B   5      21.772 -27.767 -28.880  1.00 30.00           D
HETATM  111  O   DOD B   6      17.493 -24.723 -25.229  1.00 20.00           O
HETATM  112  D1  DOD B   6      17.492 -25.660 -25.470  1.00 20.00           D
HETATM  113  D2  DOD B   6      18.304 -24.255 -25.470  1.00 20.00           D
TER     114      HOH B   6
HETATM    1  PG  ATP C   1       3.232  -0.948  -8.256  1.00 20.00      A    P
HETATM    2  O1G ATP C   1       4.490  -0.168  -7.985  1.00 20.00      A    O
HETATM    3  O2G ATP C   1       2.572  -0.420  -9.500  1.00 20.00      A    O
HETATM    4  O3G ATP C   1       3.573  -2.392  -8.447  1.00 20.00      A    O
HETATM    5  PB  ATP C   1       2.368   0.364  -5.856  1.00 20.00      A    P
HETATM    6  O1B ATP C   1       3.702   0.219  -5.162  1.00 20.00      A    O
HETATM    7  O2B ATP C   1       2.287   1.719  -6.506  1.00 20.00      A    O
HETATM    8  O3B ATP C   1       2.218  -0.793  -6.996  1.00 20.00      A    O
HETATM    9  PA  ATP C   1       0.926   1.250  -3.494  1.00 20.00      A    P
HETATM   10  O1A ATP C   1      -0.200   2.200  -3.825  1.00 20.00      A    O
HETATM   11  O2A ATP C   1       2.197   2.037  -3.253  1.00 20.00      A    O
HETATM   12  O3A ATP C   1       1.155   0.210  -4.753  1.00 20.00      A    O
HETATM   13  O5' ATP C   1       0.538   0.381  -2.115  1.00 20.00      A    O
HETATM   14  C5' ATP C   1      -0.378   0.948  -1.194  1.00 20.00      A    C
HETATM   15  C4' ATP C   1      -1.094  -0.162  -0.442  1.00 20.00      A    C
HETATM   16  O4' ATP C   1      -0.263  -0.783   0.344  1.00 20.00      A    O
HETATM   17  C3' ATP C   1      -2.258   0.478   0.549  1.00 20.00      A    C
HETATM   18  O3' ATP C   1      -3.508   0.701  -0.219  1.00 20.00      A    O
HETATM   19  C2' ATP C   1      -2.450  -0.373   1.436  1.00 20.00      A    C
HETATM   20  O2' ATP C   1      -3.564  -1.296   1.031  1.00 20.00      A    O
HETATM   21  C1' ATP C   1      -1.075  -1.183   1.553  1.00 20.00      A    C
HETATM   22  N9  ATP C   1      -0.436  -0.847   2.669  1.00 20.00      A    N
HETATM   23  C8  ATP C   1       0.891  -0.729   2.488  1.00 20.00      A    C
HETATM   24  N7  ATP C   1       1.445  -0.426   3.672  1.00 20.00      A    N
HETATM   25  C5  ATP C   1       0.455  -0.356   4.593  1.00 20.00      A    C
HETATM   26  C6  ATP C   1       0.464  -0.066   5.953  1.00 20.00      A    C
HETATM   27  N6  ATP C   1       1.576   0.468   6.779  1.00 20.00      A    N
HETATM   28  N1  ATP C   1      -0.678  -0.046   6.637  1.00 20.00      A    N
HETATM   29  C2  ATP C   1      -1.855  -0.311   6.004  1.00 20.00      A    C
HETATM   30  N3  ATP C   1      -1.865  -0.595   4.675  1.00 20.00      A    N
HETATM   31  C4  ATP C   1      -0.703  -0.615   3.975  1.00 20.00      A    C
HETATM   32 H5'1 ATP C   1       0.163   1.581  -0.485  1.00 20.00      A    H
HETATM   33 H5'2 ATP C   1      -1.105   1.545  -1.731  1.00 20.00      A    H
HETATM   34  H4' ATP C   1      -1.548  -0.870  -1.145  1.00 20.00      A    H
HETATM   35  H3' ATP C   1      -1.912   1.393   0.983  1.00 20.00      A    H
HETATM   36 HO3' ATP C   1      -3.892   1.548   0.031  1.00 20.00      A    H
HETATM   37  H2' ATP C   1      -2.680   0.114   2.384  1.00 20.00      A    H
HETATM   38 HO2' ATP C   1      -4.066  -1.553   1.810  1.00 20.00      A    H
HETATM   39  H1' ATP C   1      -1.273  -2.266   1.534  1.00 20.00      A    H
HETATM   40  H8  ATP C   1       1.431  -0.922   1.560  1.00 20.00      A    H
HETATM   41 HN61 ATP C   1       2.477   0.646   6.359  1.00 20.00      A    H
HETATM   42 HN62 ATP C   1       1.445   0.610   7.769  1.00 20.00      A    H
HETATM   43  H2  ATP C   1      -2.804  -0.294   6.573  1.00 20.00      A    H
TER
END
'''

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
