
from __future__ import absolute_import, division, print_function
from mmtbx.validation.molprobity import mp_geo
from iotbx.data_manager import DataManager
from mmtbx.model import manager
from libtbx.utils import null_out
from iotbx import pdb
from libtbx.test_utils import approx_equal
import time
import json

pdb_str_1 = """
ATOM     60  N   GLU A   9       7.757   9.623  27.829  1.00 12.74           N
ATOM     61  CA  GLU A   9       7.695  10.769  26.927  1.00 12.72           C
ATOM     62  C   GLU A   9       7.351  12.056  27.684  1.00 12.59           C
ATOM     63  O   GLU A   9       6.321  12.134  28.348  1.00 14.32           O
ATOM     64  CB  GLU A   9       6.633  10.459  25.859  1.00 13.52           C
ATOM     65  CG  GLU A   9       6.323  11.587  24.899  1.00 15.93           C
ATOM     66  CD  GLU A   9       7.533  12.006  24.105  1.00 16.48           C
ATOM     67  OE1 GLU A   9       8.059  11.171  23.337  1.00 19.96           O
ATOM     68  OE2 GLU A   9       7.960  13.164  24.257  1.00 18.90           O
ATOM     69  N   ASP A  10       8.216  13.062  27.589  1.00 14.68           N
ATOM     70  CA  ASP A  10       7.992  14.328  28.283  1.00 15.47           C
ATOM     71  C   ASP A  10       6.987  15.255  27.598  1.00 15.82           C
ATOM     72  O   ASP A  10       6.315  16.049  28.258  1.00 15.43           O
ATOM     73  CB  ASP A  10       9.322  15.075  28.453  1.00 16.97           C
ATOM     74  CG AASP A  10       9.139  16.482  28.993  0.50 16.48           C
ATOM     75  CG BASP A  10       9.941  15.472  27.128  0.50 17.05           C
ATOM     76  OD1AASP A  10       9.189  17.442  28.197  0.50 17.34           O
ATOM     77  OD1BASP A  10      11.349  14.573  26.359  0.50 18.03           O
ATOM     78  OD2AASP A  10       8.934  16.630  30.215  0.50 17.88           O
ATOM     79  OD2BASP A  10      10.016  16.689  26.852  0.50 20.20           O
"""

pdb_str_2 = """
ATOM    516  P     C B 115      25.417  97.388  92.531  1.00161.67           P
ATOM    517  OP1   C B 115      25.620  96.786  93.873  1.00161.61           O
ATOM    518  OP2   C B 115      26.429  97.148  91.471  1.00157.52           O
ATOM    519  O5'   C B 115      23.989  96.929  91.981  1.00165.73           O
ATOM    520  C5'   C B 115      22.830  97.049  92.795  1.00166.81           C
ATOM    521  C4'   C B 115      21.568  96.608  92.088  1.00164.72           C
ATOM    522  O4'   C B 115      21.307  97.449  90.932  1.00165.07           O
ATOM    523  C3'   C B 115      21.553  95.203  91.511  1.00161.21           C
ATOM    524  O3'   C B 115      21.376  94.188  92.483  1.00162.42           O
ATOM    525  C2'   C B 115      20.419  95.291  90.500  1.00160.96           C
ATOM    526  O2'   C B 115      19.151  95.246  91.141  1.00163.14           O
ATOM    527  C1'   C B 115      20.634  96.695  89.939  1.00161.73           C
ATOM    528  N1    C B 115      21.449  96.670  88.703  1.00156.58           N
ATOM    529  C2    C B 115      20.847  96.216  87.521  1.00153.43           C
ATOM    530  O2    C B 115      19.660  95.854  87.547  1.00152.48           O
ATOM    531  N3    C B 115      21.572  96.181  86.379  1.00150.07           N
ATOM    532  C4    C B 115      22.845  96.576  86.385  1.00149.02           C
ATOM    533  N4    C B 115      23.520  96.526  85.235  1.00146.94           N
ATOM    534  C5    C B 115      23.484  97.040  87.572  1.00150.48           C
ATOM    535  C6    C B 115      22.758  97.069  88.698  1.00153.67           C
ATOM    547  P  A  A B 116      21.844  92.681  92.170  0.51179.08           P
ATOM    548  OP1A  A B 116      21.573  91.863  93.379  0.51179.56           O
ATOM    549  OP2A  A B 116      23.220  92.732  91.614  0.51177.92           O
ATOM    550  O5'A  A B 116      20.860  92.191  91.014  0.51167.56           O
ATOM    551  C5'A  A B 116      19.551  91.726  91.313  0.51164.82           C
ATOM    552  C4'A  A B 116      18.742  91.493  90.060  0.51155.92           C
ATOM    553  O4'A  A B 116      18.898  92.626  89.165  0.51149.57           O
ATOM    554  C3'A  A B 116      19.141  90.293  89.209  0.51151.37           C
ATOM    555  O3'A  A B 116      18.607  89.067  89.684  0.51153.31           O
ATOM    556  C2'A  A B 116      18.628  90.683  87.831  0.51140.91           C
ATOM    557  O2'A  A B 116      17.227  90.473  87.724  0.51147.22           O
ATOM    558  C1'A  A B 116      18.903  92.185  87.824  0.51133.46           C
ATOM    559  N9 A  A B 116      20.215  92.506  87.229  0.51115.10           N
ATOM    560  C8 A  A B 116      21.352  92.925  87.874  0.51103.22           C
ATOM    561  N7 A  A B 116      22.369  93.136  87.075  0.51 94.11           N
ATOM    562  C5 A  A B 116      21.867  92.833  85.817  0.51 93.20           C
ATOM    563  C6 A  A B 116      22.449  92.855  84.538  0.51 84.33           C
ATOM    564  N6 A  A B 116      23.716  93.210  84.308  0.51 83.97           N
ATOM    565  N1 A  A B 116      21.676  92.496  83.490  0.51 83.45           N
ATOM    566  C2 A  A B 116      20.406  92.141  83.722  0.51 89.16           C
ATOM    567  N3 A  A B 116      19.746  92.081  84.877  0.51100.01           N
ATOM    568  C4 A  A B 116      20.543  92.444  85.897  0.51106.03           C
ATOM    580  P  B  A B 116      21.673  92.654  92.107  0.49176.43           P
ATOM    581  OP1B  A B 116      22.052  91.946  93.355  0.49179.31           O
ATOM    582  OP2B  A B 116      22.589  92.619  90.939  0.49177.26           O
ATOM    583  O5'B  A B 116      20.256  92.106  91.631  0.49166.83           O
ATOM    584  C5'B  A B 116      20.155  90.953  90.811  0.49162.34           C
ATOM    585  C4'B  A B 116      18.753  90.785  90.285  0.49156.23           C
ATOM    586  O4'B  A B 116      18.380  91.954  89.503  0.49153.24           O
ATOM    587  C3'B  A B 116      18.530  89.625  89.328  0.49153.93           C
ATOM    588  O3'B  A B 116      18.398  88.370  89.971  0.49155.29           O
ATOM    589  C2'B  A B 116      17.290  90.074  88.574  0.49152.57           C
ATOM    590  O2'B  A B 116      16.122  89.934  89.370  0.49160.43           O
ATOM    591  C1'B  A B 116      17.592  91.558  88.399  0.49146.63           C
ATOM    592  N9 B  A B 116      18.368  91.795  87.168  0.49131.17           N
ATOM    593  C8 B  A B 116      19.613  92.360  87.046  0.49119.64           C
ATOM    594  N7 B  A B 116      20.046  92.418  85.809  0.49108.75           N
ATOM    595  C5 B  A B 116      19.018  91.847  85.070  0.49102.94           C
ATOM    596  C6 B  A B 116      18.857  91.610  83.694  0.49 95.23           C
ATOM    597  N6 B  A B 116      19.769  91.933  82.774  0.49 88.20           N
ATOM    598  N1 B  A B 116      17.711  91.021  83.287  0.49 91.58           N
ATOM    599  C2 B  A B 116      16.795  90.696  84.207  0.49 93.21           C
ATOM    600  N3 B  A B 116      16.833  90.868  85.525  0.49101.29           N
ATOM    601  C4 B  A B 116      17.982  91.457  85.895  0.49111.07           C
ATOM    613  C4'   G B 117      18.093  86.821  85.795  1.00149.42           C
ATOM    614  O4'   G B 117      18.411  88.142  85.277  1.00150.38           O
ATOM    615  C3'   G B 117      19.188  85.938  85.221  1.00152.20           C
ATOM    616  O3'   G B 117      18.860  84.561  85.207  1.00157.98           O
ATOM    617  C2'   G B 117      19.390  86.537  83.840  1.00148.15           C
ATOM    618  O2'   G B 117      18.348  86.156  82.953  1.00150.85           O
ATOM    619  C1'   G B 117      19.255  88.022  84.148  1.00142.90           C
ATOM    620  N9    G B 117      20.565  88.616  84.472  1.00129.09           N
ATOM    621  C8    G B 117      21.040  88.952  85.716  1.00125.15           C
ATOM    622  N7    G B 117      22.244  89.453  85.691  1.00119.89           N
ATOM    623  C5    G B 117      22.593  89.440  84.347  1.00115.13           C
ATOM    624  C6    G B 117      23.787  89.862  83.702  1.00106.78           C
ATOM    625  O6    G B 117      24.807  90.345  84.211  1.00103.03           O
ATOM    626  N1    G B 117      23.723  89.674  82.323  1.00103.58           N
ATOM    627  C2    G B 117      22.649  89.144  81.651  1.00105.31           C
ATOM    628  N2    G B 117      22.777  89.041  80.320  1.00103.13           N
ATOM    629  N3    G B 117      21.530  88.746  82.240  1.00112.17           N
ATOM    630  C4    G B 117      21.569  88.921  83.580  1.00119.07           C
ATOM    640  P  A  G B 117      19.377  87.681  89.414  0.57155.42           P
ATOM    641  OP1A  G B 117      18.636  86.602  90.115  0.57154.98           O
ATOM    642  OP2A  G B 117      20.818  87.886  89.703  0.57156.38           O
ATOM    643  O5'A  G B 117      19.216  87.435  87.847  0.57148.01           O
ATOM    644  C5'A  G B 117      18.022  86.888  87.300  0.57148.38           C
ATOM    647  P  B  G B 117      19.331  87.141  89.518  0.43155.13           P
ATOM    648  OP1B  G B 117      18.722  85.885  90.026  0.43152.30           O
ATOM    649  OP2B  G B 117      20.734  87.474  89.871  0.43155.64           O
ATOM    650  O5'B  G B 117      19.213  87.135  87.929  0.43147.38           O
ATOM    651  C5'B  G B 117      17.973  86.851  87.299  0.43148.35           C
ATOM    654  P     C B 118      20.023  83.460  85.359  1.00146.31           P
ATOM    655  OP1   C B 118      19.377  82.123  85.410  1.00144.77           O
ATOM    656  OP2   C B 118      20.921  83.881  86.465  1.00146.69           O
ATOM    657  O5'   C B 118      20.835  83.565  83.989  1.00123.79           O
ATOM    658  C5'   C B 118      20.203  83.271  82.751  1.00120.61           C
ATOM    659  C4'   C B 118      20.991  83.775  81.563  1.00114.71           C
ATOM    660  O4'   C B 118      21.385  85.159  81.765  1.00115.11           O
ATOM    661  C3'   C B 118      22.301  83.069  81.260  1.00107.01           C
ATOM    662  O3'   C B 118      22.124  81.826  80.602  1.00104.55           O
ATOM    663  C2'   C B 118      23.035  84.105  80.419  1.00107.93           C
ATOM    664  O2'   C B 118      22.544  84.120  79.087  1.00108.92           O
ATOM    665  C1'   C B 118      22.610  85.406  81.104  1.00111.02           C
ATOM    666  N1    C B 118      23.615  85.886  82.083  1.00107.14           N
ATOM    667  C2    C B 118      24.827  86.373  81.586  1.00106.46           C
ATOM    668  O2    C B 118      25.014  86.373  80.361  1.00109.38           O
ATOM    669  N3    C B 118      25.774  86.826  82.437  1.00103.02           N
ATOM    670  C4    C B 118      25.543  86.811  83.747  1.00103.40           C
ATOM    671  N4    C B 118      26.507  87.269  84.549  1.00103.30           N
ATOM    672  C5    C B 118      24.316  86.326  84.289  1.00103.70           C
ATOM    673  C6    C B 118      23.387  85.879  83.432  1.00104.01           C
"""

def exercise_mp_geo():
  with open('mp_geo.pdb', 'w') as f:
    f.write(pdb_str_1)
  args = ['pdb=mp_geo.pdb',
          'out_file=mp_geo.out',
          'outliers_only=False',
          'bonds_and_angles=True']
  mp_geo.run(args)
  f = open('mp_geo.out', 'r')
  # Strip out newline and carriage return chars to
  # prevent platform-specific errors.
  lines = [i.rstrip('\n\r') for i in f.readlines()]
  assert 'mp_geo.pdb: A:  10: :B:ASP:CG--OD1:1.839:31.054:PROTEIN' in lines
  assert \
    'mp_geo.pdb: A:  10: :B:ASP:OD1-CG-OD2:109.733:5.486:PROTEIN' in lines
  f.close()

  with open('mp_geo.pdb', 'w') as f:
    f.write(pdb_str_2)
  args = ['pdb=mp_geo.pdb',
          'out_file=mp_geo.out',
          'rna_backbone=True']
  mp_geo.run(args)
  f = open('mp_geo.out', 'r')
  lines = [i.rstrip('\n\r') for i in f.readlines()]
  f.close()
  assert lines[0] == \
    ' :1: B: 115: : :  C:__?__:178.072:55.525:76.414:-158.236:-67.172'
  assert lines[1] == \
    ' :1: B: 116: :A:  A:-80.906:172.347:71.412:81.732:-151.720:-70.053'
  assert lines[2] == \
    ' :1: B: 117: :A:  G:-80.508:177.077:61.989:80.676:-149.321:-69.591'
  assert lines[3] == \
    ' :1: B: 118: : :  C:-62.565:164.517:69.624:78.058:__?__:__?__'

def exercise_mp_geo_bonds_angles_json():
  dm = DataManager()
  #print(help(dm))
  dm.process_model_str("1",pdb_str_1)
  model = dm.get_model("1")
  model.set_stop_for_unknowns(False)
  hierarchy = model.get_hierarchy()
  p = manager.get_default_pdb_interpretation_params()
  ##print(dir(p.pdb_interpretation))
  p.pdb_interpretation.allow_polymer_cross_special_position=True
  p.pdb_interpretation.flip_symmetric_amino_acids=False
  p.pdb_interpretation.clash_guard.nonbonded_distance_threshold = None
  p.pdb_interpretation.restraints_library.cdl = True
  model.set_log(log = null_out())
  model.process(make_restraints=True, pdb_interpretation_params=p)
  geometry = model.get_restraints_manager().geometry
  atoms = hierarchy.atoms()
  rc = mp_geo.get_bond_and_angle_outliers(
           pdb_hierarchy=hierarchy,
           xray_structure=model.get_xray_structure(),
           geometry_restraints_manager=geometry,
           use_segids=False,
           outliers_only=False)
  bonds = rc.bonds
  angles = rc.angles
  import pprint
  bonds_json = json.loads(bonds.as_JSON())
  #pprint.pprint(bonds_json)
  assert len(bonds_json['flat_results'])==19, "tst_mp_validate_bonds total number of bonds changed, now: "+str(len(bonds_json['flat_results']))
  assert approx_equal(bonds_json['flat_results'][18]["sigma"], 0.019), "tst_mp_validate_bonds json output last sigma value changed, now: "+str(bonds_json['flat_results'][18]["sigma"])
  assert bonds_json['summary_results'][""]["num_outliers"] == 1, "tst_mp_validate_bonds json summary output total number of outliers changed, now: "+str(bonds_json['summary_results'][""]["num_outliers"])
  assert bonds_json['summary_results'][""]["num_total"]==19, "tst_mp_validate_bonds json summary output total number of bonds changed, now: "+str(bonds_json['summary_results'][""]["num_total"])
  angles_json = json.loads(angles.as_JSON())
  #pprint.pprint(angles_json)
  assert len(angles_json['flat_results'])==24, "tst_mp_validate_bonds total number of angles changed, now: "+str(len(angles_json['flat_results']))
  assert approx_equal(angles_json['flat_results'][23]["sigma"], 2.4), "tst_mp_validate_bonds json output last sigma value changed, now: "+str(angles_json['flat_results'][23]["sigma"])
  assert angles_json['summary_results'][""]["num_outliers"] == 1, "tst_mp_validate_bonds json summary output total number of outliers changed, now: "+str(angles_json['summary_results'][""]["num_outliers"])
  assert angles_json['summary_results'][""]["num_total"]==24, "tst_mp_validate_bonds json summary output total number of angles changed, now: "+str(angles_json['summary_results'][""]["num_total"])
  return bonds_json, angles_json

if (__name__ == "__main__"):
  t0 = time.time()
  exercise_mp_geo()
  exercise_mp_geo_bonds_angles_json()
  print("OK. Time: %8.3f"%(time.time()-t0))
