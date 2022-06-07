from __future__ import absolute_import, division, print_function
from mmtbx.validation.molprobity import mp_geo
import time
from libtbx.test_utils import approx_equal

def test_normal_val():
  """Regular valine. Chiral centers at CA and CB."""
  pdb_str = """\
ATOM     37  N   VAL A   5      28.260  33.943  11.096  1.00  4.44           N
ATOM     38  CA  VAL A   5      28.605  33.965  12.503  1.00  3.87           C
ATOM     39  C   VAL A   5      28.638  35.461  12.900  1.00  4.93           C
ATOM     40  O   VAL A   5      29.522  36.103  12.320  1.00  6.84           O
ATOM     41  CB  VAL A   5      29.963  33.317  12.814  1.00  2.99           C
ATOM     42  CG1 VAL A   5      30.211  33.394  14.304  1.00  5.28           C
ATOM     43  CG2 VAL A   5      29.957  31.838  12.352  1.00  9.13           C
TER
END
"""
  with open('chiral_volume_validation_cases_val_normal.pdb', 'w') as f:
    f.write(pdb_str)
  args = ['pdb=chiral_volume_validation_cases_val_normal.pdb',
          'out_file=chiral_volume_validation_cases_val_normal.out',
          'outliers_only=False',
          'bonds_and_angles=True']
  mp_geo.run(args)
  with open('chiral_volume_validation_cases_val_normal.out', 'r') as f:
    lines = [i.rstrip('\n\r') for i in f.readlines()]
  #chiral_volume_validation_cases_val_base.pdb: A:   5: : :VAL:CA:2.628:0.932:PROTEIN
  #chiral_volume_validation_cases_val_base.pdb: A:   5: : :VAL:CB:-2.762:0.664:PROTEIN
  ca_volume = None
  ca_sigma = None
  cb_volume = None
  cb_sigma = None
  for line in lines:
    #print(line)
    chiral_center = line.split(":")[6]
    if chiral_center == "CA":
      ca_volume = float(line.split(":")[7])
      ca_sigma = float(line.split(":")[8])
    elif chiral_center == "CB":
      cb_volume = float(line.split(":")[7])
      cb_sigma = float(line.split(":")[8])
  assert approx_equal(ca_volume, 2.628, eps = 0.01)
  assert approx_equal(ca_sigma, 0.932, eps = 0.01)
  assert approx_equal(cb_volume, -2.762, eps = 0.01)
  assert approx_equal(cb_sigma, 0.664, eps = 0.01)

def test_misnamed_val():
  """Atom positions for CG1 and CG2 have been swapped.
This should create a >20sigma chiral volume outlier at CB.
This is a 'pseudochiral naming error'.
Validation wants to be able to catch these naming errors."""
  pdb_str = """\
ATOM     37  N   VAL A   5      28.260  33.943  11.096  1.00  4.44           N
ATOM     38  CA  VAL A   5      28.605  33.965  12.503  1.00  3.87           C
ATOM     39  C   VAL A   5      28.638  35.461  12.900  1.00  4.93           C
ATOM     40  O   VAL A   5      29.522  36.103  12.320  1.00  6.84           O
ATOM     41  CB  VAL A   5      29.963  33.317  12.814  1.00  2.99           C
ATOM     42  CG1 VAL A   5      29.957  31.838  12.352  1.00  9.13           C
ATOM     43  CG2 VAL A   5      30.211  33.394  14.304  1.00  5.28           C
END
"""
  with open('chiral_volume_validation_cases_val_misnamed.pdb', 'w') as f:
    f.write(pdb_str)
  args = ['pdb=chiral_volume_validation_cases_val_misnamed.pdb',
          'out_file=chiral_volume_validation_cases_val_misnamed.out',
          'outliers_only=False',
          'bonds_and_angles=True']
  mp_geo.run(args)
  with open('chiral_volume_validation_cases_val_misnamed.out', 'r') as f:
    lines = [i.rstrip('\n\r') for i in f.readlines()]
  #chiral_volume_validation_cases_val_misnamed.pdb: A:   5: : :VAL:CA:2.628:0.932:PROTEIN
  #chiral_volume_validation_cases_val_misnamed.pdb: A:   5: : :VAL:CB:2.762:26.956:PROTEIN
  ca_volume = None
  ca_sigma = None
  cb_volume = None
  cb_sigma = None
  for line in lines:
    #print(line)
    chiral_center = line.split(":")[6]
    if chiral_center == "CA":
      ca_volume = float(line.split(":")[7])
      ca_sigma = float(line.split(":")[8])
    elif chiral_center == "CB":
      cb_volume = float(line.split(":")[7])
      cb_sigma = float(line.split(":")[8])
  assert approx_equal(ca_volume, 2.628, eps = 0.01) #same
  assert approx_equal(ca_sigma, 0.932, eps = 0.01) #same
  assert approx_equal(cb_volume, 2.762, eps = 0.01) #sign change
  assert approx_equal(cb_sigma, 26.956, eps = 0.01) #>20sigma outlier

################################################################################

def test_normal_leu():
  """Regular leucine. Chiral centers at CA and CG"""
  pdb_str = """\
ATOM     60  N   LEU A   8      30.132  40.069  18.642  1.00  9.84           N
ATOM     61  CA  LEU A   8      29.607  41.180  19.467  1.00 14.15           C
ATOM     62  C   LEU A   8      30.075  42.538  18.984  1.00 17.37           C
ATOM     63  O   LEU A   8      29.586  43.570  19.483  1.00 17.01           O
ATOM     64  CB  LEU A   8      29.919  40.890  20.938  1.00 16.63           C
ATOM     65  CG  LEU A   8      29.183  39.722  21.581  1.00 18.88           C
ATOM     66  CD1 LEU A   8      29.308  39.750  23.095  1.00 19.31           C
ATOM     67  CD2 LEU A   8      27.700  39.721  21.228  1.00 18.59           C
TER
END
"""
  with open('chiral_volume_validation_cases_leu_normal.pdb', 'w') as f:
    f.write(pdb_str)
  args = ['pdb=chiral_volume_validation_cases_leu_normal.pdb',
          'out_file=chiral_volume_validation_cases_leu_normal.out',
          'outliers_only=False',
          'bonds_and_angles=True']
  mp_geo.run(args)
  with open('chiral_volume_validation_cases_leu_normal.out', 'r') as f:
    lines = [i.rstrip('\n\r') for i in f.readlines()]
  #chiral_volume_validation_cases_leu_normal.pdb: A:   8: : :LEU:CA:2.369:0.707:PROTEIN
  #chiral_volume_validation_cases_leu_normal.pdb: A:   8: : :LEU:CG:-2.604:0.070:PROTEIN
  ca_volume = None
  ca_sigma = None
  cg_volume = None
  cg_sigma = None
  for line in lines:
    #print(line)
    chiral_center = line.split(":")[6]
    if chiral_center == "CA":
      ca_volume = float(line.split(":")[7])
      ca_sigma = float(line.split(":")[8])
    elif chiral_center == "CG":
      cg_volume = float(line.split(":")[7])
      cg_sigma = float(line.split(":")[8])
  assert approx_equal(ca_volume, 2.369, eps = 0.01)
  assert approx_equal(ca_sigma, 0.707, eps = 0.01)
  assert approx_equal(cg_volume, -2.604, eps = 0.01)
  assert approx_equal(cg_sigma, 0.070, eps = 0.01)

def test_misnamed_leu():
  """Atom positions for CD1 and CD2 have been swapped.
This should create a >20sigma chiral volume outlier at CG.
This is a 'pseudochiral naming error'.
Validation wants to be able to catch these naming errors."""
  pdb_str = """\
ATOM     60  N   LEU A   8      30.132  40.069  18.642  1.00  9.84           N
ATOM     61  CA  LEU A   8      29.607  41.180  19.467  1.00 14.15           C
ATOM     62  C   LEU A   8      30.075  42.538  18.984  1.00 17.37           C
ATOM     63  O   LEU A   8      29.586  43.570  19.483  1.00 17.01           O
ATOM     64  CB  LEU A   8      29.919  40.890  20.938  1.00 16.63           C
ATOM     65  CG  LEU A   8      29.183  39.722  21.581  1.00 18.88           C
ATOM     66  CD1 LEU A   8      27.700  39.721  21.228  1.00 18.59           C
ATOM     67  CD2 LEU A   8      29.308  39.750  23.095  1.00 19.31           C
END
"""
  with open('chiral_volume_validation_cases_leu_misnamed.pdb', 'w') as f:
    f.write(pdb_str)
  args = ['pdb=chiral_volume_validation_cases_leu_misnamed.pdb',
          'out_file=chiral_volume_validation_cases_leu_misnamed.out',
          'outliers_only=False',
          'bonds_and_angles=True']
  mp_geo.run(args)
  with open('chiral_volume_validation_cases_leu_misnamed.out', 'r') as f:
    lines = [i.rstrip('\n\r') for i in f.readlines()]
  #chiral_volume_validation_cases_leu_misnamed.pdb: A:   8: : :LEU:CA:2.369:0.707:PROTEIN
  #chiral_volume_validation_cases_leu_misnamed.pdb: A:   8: : :LEU:CG:2.604:25.967:PROTEIN
  ca_volume = None
  ca_sigma = None
  cg_volume = None
  cg_sigma = None
  for line in lines:
    #print(line)
    chiral_center = line.split(":")[6]
    if chiral_center == "CA":
      ca_volume = float(line.split(":")[7])
      ca_sigma = float(line.split(":")[8])
    elif chiral_center == "CG":
      cg_volume = float(line.split(":")[7])
      cg_sigma = float(line.split(":")[8])
  assert approx_equal(ca_volume, 2.369, eps = 0.01) #same
  assert approx_equal(ca_sigma, 0.707, eps = 0.01) #same
  assert approx_equal(cg_volume, 2.604, eps = 0.01) #sign change
  assert approx_equal(cg_sigma, 25.967, eps = 0.01) #>20sigma outlier

################################################################################

def test_normal_sf4():
  """Regular sf4 iron-sulfur cluster. Chiral centers at FEs"""
  pdb_str = """\
HETATM14724 FE1  SF4 E 501       0.575  47.227  37.154  1.00 44.16          FE3+
HETATM14725 FE2  SF4 E 501       1.672  45.200  38.813  1.00 42.81          FE3+
HETATM14726 FE3  SF4 E 501       0.846  44.650  36.268  1.00 41.68          FE3+
HETATM14727 FE4  SF4 E 501      -0.960  45.318  38.094  1.00 49.33          FE3+
HETATM14728  S1  SF4 E 501       0.351  43.431  38.164  1.00 42.06           S
HETATM14729  S2  SF4 E 501      -0.963  46.038  35.912  1.00 42.31           S
HETATM14730  S3  SF4 E 501       0.157  46.850  39.397  1.00 45.58           S
HETATM14731  S4  SF4 E 501       2.547  46.073  36.852  1.00 45.73           S
TER
END
"""
  with open('chiral_volume_validation_cases_sf4_normal.pdb', 'w') as f:
    f.write(pdb_str)
  args = ['pdb=chiral_volume_validation_cases_sf4_normal.pdb',
          'out_file=chiral_volume_validation_cases_sf4_normal.out',
          'outliers_only=False',
          'bonds_and_angles=True']
  mp_geo.run(args)
  with open('chiral_volume_validation_cases_sf4_normal.out', 'r') as f:
    lines = [i.rstrip('\n\r') for i in f.readlines()]
  #chiral_volume_validation_cases_sf4_normal.pdb: E: 501: : :SF4:FE1:-10.788:1.164:UNK
  #chiral_volume_validation_cases_sf4_normal.pdb: E: 501: : :SF4:FE2:11.095:2.699:UNK
  #chiral_volume_validation_cases_sf4_normal.pdb: E: 501: : :SF4:FE3:-10.559:0.020:UNK
  #chiral_volume_validation_cases_sf4_normal.pdb: E: 501: : :SF4:FE4:10.148:2.036:UNK
  fe1_volume = None
  fe1_sigma = None
  for line in lines:
    #print(line)
    chiral_center = line.split(":")[6]
    if chiral_center == "FE1":
      fe1_volume = float(line.split(":")[7])
      fe1_sigma = float(line.split(":")[8])
  assert approx_equal(fe1_volume, -10.788, eps = 0.01)
  assert approx_equal(fe1_sigma, 1.164, eps = 0.01)

def test_misnamed_sf4():
  """Atom positions for S2 and S3 have been swapped.
This should create a >>>20sigma chiral volume outlier at FE1.
Other FEs are affected but not tested here
This is a 'pseudochiral naming error'.
Validation wants to be able to catch these naming errors."""
  pdb_str = """\
HETATM14724 FE1  SF4 E 501       0.575  47.227  37.154  1.00 44.16          FE3+
HETATM14725 FE2  SF4 E 501       1.672  45.200  38.813  1.00 42.81          FE3+
HETATM14726 FE3  SF4 E 501       0.846  44.650  36.268  1.00 41.68          FE3+
HETATM14727 FE4  SF4 E 501      -0.960  45.318  38.094  1.00 49.33          FE3+
HETATM14728  S1  SF4 E 501       0.351  43.431  38.164  1.00 42.06           S
HETATM14729  S3  SF4 E 501      -0.963  46.038  35.912  1.00 42.31           S
HETATM14730  S2  SF4 E 501       0.157  46.850  39.397  1.00 45.58           S
HETATM14731  S4  SF4 E 501       2.547  46.073  36.852  1.00 45.73           S
END
"""
  with open('chiral_volume_validation_cases_sf4_misnamed.pdb', 'w') as f:
    f.write(pdb_str)
  args = ['pdb=chiral_volume_validation_cases_sf4_misnamed.pdb',
          'out_file=chiral_volume_validation_cases_sf4_misnamed.out',
          'outliers_only=False',
          'bonds_and_angles=True']
  mp_geo.run(args)
  with open('chiral_volume_validation_cases_sf4_misnamed.out', 'r') as f:
    lines = [i.rstrip('\n\r') for i in f.readlines()]
  #chiral_volume_validation_cases_sf4_misnamed.pdb: E: 501: : :SF4:FE1:10.788:106.712:UNK
  #chiral_volume_validation_cases_sf4_misnamed.pdb: E: 501: : :SF4:FE2:14.425:19.352:UNK
  #chiral_volume_validation_cases_sf4_misnamed.pdb: E: 501: : :SF4:FE3:-14.365:19.049:UNK
  #chiral_volume_validation_cases_sf4_misnamed.pdb: E: 501: : :SF4:FE4:-10.148:103.512:UNK
  fe1_volume = None
  fe1_sigma = None
  for line in lines:
    #print(line)
    chiral_center = line.split(":")[6]
    if chiral_center == "FE1":
      fe1_volume = float(line.split(":")[7])
      fe1_sigma = float(line.split(":")[8])
  assert approx_equal(fe1_volume, 10.788, eps = 0.01) #sign changes
  assert approx_equal(fe1_sigma, 106.712, eps = 0.01) #>>>20 sigma outlier

################################################################################

def exercise():
  test_normal_val()
  test_misnamed_val()
  test_normal_leu()
  test_misnamed_leu()
  test_normal_sf4()
  test_misnamed_sf4()

if (__name__ == "__main__"):
  t0 = time.time()
  exercise()
  print("OK. Time: %8.3f"%(time.time()-t0))
