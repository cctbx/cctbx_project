from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from libtbx import group_args
from libtbx.utils import user_plus_sys_time
import random
import mmtbx.model
import iotbx.pdb
from mmtbx.refinement.real_space import individual_sites

pdb_str_1 = """\
CRYST1   19.547   20.035   19.435  90.00  90.00  90.00 P 1
ATOM      1  N   ASP A  18      16.300  13.354  13.167  1.00 21.41           N
ATOM      2  CA  ASP A  18      15.058  14.099  13.335  1.00 34.70           C
ATOM      3  C   ASP A  18      14.270  14.157  12.030  1.00 35.50           C
ATOM      4  O   ASP A  18      14.795  14.579  10.999  1.00 38.48           O
ATOM      5  CB  ASP A  18      15.347  15.514  13.839  1.00 38.15           C
ATOM      6  CG  ASP A  18      14.084  16.327  14.048  1.00 27.39           C
ATOM      7  OD1 ASP A  18      13.508  16.260  15.154  1.00 35.78           O
ATOM      8  OD2 ASP A  18      13.669  17.035  13.106  1.00 30.67           O
ATOM      9  N   ASN A  19      13.011  13.728  12.091  1.00 22.94           N
ATOM     10  CA  ASN A  19      12.119  13.717  10.933  1.00 22.82           C
ATOM     11  C   ASN A  19      12.680  12.928   9.751  1.00 20.25           C
ATOM     12  O   ASN A  19      13.252  13.500   8.824  1.00 33.24           O
ATOM     13  CB  ASN A  19      11.763  15.145  10.505  1.00 28.83           C
ATOM     14  CG  ASN A  19      10.717  15.181   9.407  1.00 26.44           C
ATOM     15  OD1 ASN A  19       9.934  14.245   9.247  1.00 22.64           O
ATOM     16  ND2 ASN A  19      10.700  16.267   8.643  1.00 28.91           N
ATOM     17  N   TYR A  20      12.510  11.610   9.792  1.00 23.81           N
ATOM     18  CA  TYR A  20      12.997  10.741   8.728  1.00 35.69           C
ATOM     19  C   TYR A  20      11.911   9.771   8.272  1.00 28.94           C
ATOM     20  O   TYR A  20      11.633   8.780   8.948  1.00 25.26           O
ATOM     21  CB  TYR A  20      14.236   9.970   9.191  1.00 33.37           C
ATOM     22  CG  TYR A  20      14.850   9.092   8.123  1.00 31.48           C
ATOM     23  CD1 TYR A  20      15.713   9.623   7.173  1.00 36.91           C
ATOM     24  CD2 TYR A  20      14.573   7.732   8.069  1.00 21.36           C
ATOM     25  CE1 TYR A  20      16.279   8.825   6.196  1.00 36.56           C
ATOM     26  CE2 TYR A  20      15.134   6.927   7.096  1.00 29.56           C
ATOM     27  CZ  TYR A  20      15.986   7.478   6.162  1.00 35.08           C
ATOM     28  OH  TYR A  20      16.547   6.680   5.192  1.00 38.73           O
ATOM     29  N   ARG A  21      11.308  10.071   7.124  1.00 38.95           N
ATOM     30  CA  ARG A  21      10.249   9.249   6.540  1.00 38.77           C
ATOM     31  C   ARG A  21       9.065   9.062   7.488  1.00 27.69           C
ATOM     32  O   ARG A  21       8.999   8.084   8.233  1.00 22.82           O
ATOM     33  CB  ARG A  21      10.797   7.891   6.085  1.00 20.00           C
ATOM     34  CG  ARG A  21       9.800   7.043   5.310  1.00 20.00           C
ATOM     35  CD  ARG A  21      10.411   5.715   4.894  1.00 20.00           C
ATOM     36  NE  ARG A  21       9.467   4.889   4.147  1.00 20.00           N
ATOM     37  CZ  ARG A  21       9.747   3.681   3.668  1.00 20.00           C
ATOM     38  NH1 ARG A  21      10.949   3.154   3.857  1.00 20.00           N
ATOM     39  NH2 ARG A  21       8.826   3.000   3.000  1.00 20.00           N
ATOM     40  N   GLY A  22       8.132  10.008   7.452  1.00 24.85           N
ATOM     41  CA  GLY A  22       6.954   9.950   8.297  1.00 29.53           C
ATOM     42  C   GLY A  22       7.262  10.260   9.749  1.00 33.22           C
ATOM     43  O   GLY A  22       7.136  11.402  10.190  1.00 30.06           O
ATOM     44  N   TYR A  23       7.668   9.236  10.494  1.00 27.25           N
ATOM     45  CA  TYR A  23       7.994   9.398  11.906  1.00 34.16           C
ATOM     46  C   TYR A  23       9.344  10.084  12.086  1.00 23.48           C
ATOM     47  O   TYR A  23      10.142  10.160  11.151  1.00 39.30           O
ATOM     48  CB  TYR A  23       7.994   8.042  12.615  1.00 29.65           C
ATOM     49  CG  TYR A  23       6.662   7.327  12.574  1.00 34.88           C
ATOM     50  CD1 TYR A  23       5.698   7.556  13.547  1.00 30.77           C
ATOM     51  CD2 TYR A  23       6.369   6.422  11.562  1.00 32.29           C
ATOM     52  CE1 TYR A  23       4.479   6.905  13.514  1.00 39.91           C
ATOM     53  CE2 TYR A  23       5.153   5.766  11.520  1.00 30.45           C
ATOM     54  CZ  TYR A  23       4.212   6.011  12.498  1.00 37.13           C
ATOM     55  OH  TYR A  23       3.000   5.360  12.460  1.00 38.59           O
ATOM     56  N   SER A  24       9.593  10.581  13.293  1.00 36.25           N
ATOM     57  CA  SER A  24      10.847  11.260  13.598  1.00 26.44           C
ATOM     58  C   SER A  24      11.747  10.391  14.470  1.00 20.40           C
ATOM     59  O   SER A  24      11.284   9.756  15.417  1.00 39.27           O
ATOM     60  CB  SER A  24      10.576  12.597  14.291  1.00 21.07           C
ATOM     61  OG  SER A  24       9.802  13.451  13.466  1.00 34.95           O
ATOM     62  N   LEU A  25      13.035  10.369  14.144  1.00 22.77           N
ATOM     63  CA  LEU A  25      14.003   9.579  14.896  1.00 36.60           C
ATOM     64  C   LEU A  25      14.997  10.476  15.627  1.00 28.89           C
ATOM     65  O   LEU A  25      14.608  11.319  16.435  1.00 38.59           O
ATOM     66  CB  LEU A  25      14.750   8.606  13.976  1.00 27.50           C
ATOM     67  CG  LEU A  25      14.011   7.358  13.478  1.00 38.78           C
ATOM     68  CD1 LEU A  25      13.007   7.690  12.381  1.00 22.10           C
ATOM     69  CD2 LEU A  25      15.004   6.308  13.000  1.00 26.92           C
TER
END
"""

def get_pdb_inputs(pdb_str):
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  model = mmtbx.model.manager(model_input=pdb_inp)
  model.process(make_restraints=True)
  return group_args(
    ph  = model.get_hierarchy(),
    grm = model.get_restraints_manager(),
    xrs = model.get_xray_structure())

def exercise(d_min = 3.5):
  pi = get_pdb_inputs(pdb_str=pdb_str_1)
  selection = flex.bool(pi.xrs.scatterers().size(), True)
  f_obs = abs(pi.xrs.structure_factors(d_min = d_min).f_calc())
  r_free_flags = f_obs.generate_r_free_flags()
  for d_min in [1, 2, 3]:
    print("d_min:", d_min)
    f_calc = pi.xrs.structure_factors(d_min = d_min).f_calc()
    fft_map = f_calc.fft_map(resolution_factor=0.25)
    fft_map.apply_sigma_scaling()
    target_map = fft_map.real_map_unpadded()
    rsr_simple_refiner = individual_sites.simple(
      target_map                  = target_map,
      selection                   = selection,
      real_space_gradients_delta  = d_min/4,
      max_iterations              = 150,
      geometry_restraints_manager = pi.grm.geometry)
    for shake_size in [1,]:
      print("  shake_size:", shake_size)
      for p in [(0.01, 1.0), (0.03, 3.0)]:
        print("    target:", p)
        w_opt = flex.double()
        for start_value in [0, 1000]:
          xrs_poor = pi.xrs.deep_copy_scatterers()
          random.seed(0)
          flex.set_random_seed(0)
          xrs_poor.shake_sites_in_place(mean_distance = shake_size)
          #
          refined = individual_sites.refinery(
            refiner                  = rsr_simple_refiner,
            xray_structure           = xrs_poor,
            start_trial_weight_value = start_value,
            rms_bonds_limit          = p[0],
            rms_angles_limit         = p[1])
          w_opt.append(refined.weight_final)
          dist = flex.mean(flex.sqrt((pi.xrs.sites_cart() -
            refined.sites_cart_result).dot()))
          print("      start_value:", start_value,refined.weight_final, \
            refined.rms_bonds_final, refined.rms_angles_final, dist)
          assert refined.rms_bonds_final  <= p[0]
          assert refined.rms_angles_final <= p[1]

if(__name__ == "__main__"):
  timer = user_plus_sys_time()
  exercise()
  print("Time: %6.2f" % timer.elapsed())
