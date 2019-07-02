from __future__ import absolute_import, division, print_function
import mmtbx.monomer_library.pdb_interpretation
import iotbx.mtz
from cctbx.array_family import flex
import time
from mmtbx import monomer_library
import iotbx.pdb
import sys
import mmtbx.refinement.real_space.fit_water

pdb_answer = """\
CRYST1   29.709   30.639   28.034  90.00  90.00  90.00 P 1
ATOM      1  CB  MET A   1      11.971  15.000  12.927  1.00 10.00      A    C
ATOM      2  CG  MET A   1      12.074  15.930  14.126  1.00 10.00      A    C
ATOM      3  SD  MET A   1      12.820  17.522  13.722  1.00 10.00      A    S
ATOM      4  CE  MET A   1      14.460  17.016  13.210  1.00 10.00      A    C
ATOM      5  C   MET A   1      11.191  12.831  11.962  1.00 10.00      A    C
ATOM      6  O   MET A   1      12.202  12.489  11.349  1.00 10.00      A    O
ATOM      7  CA  MET A   1      11.294  13.663  13.235  1.00 10.00      A    C
ATOM      8  N   MET A   1      12.057  12.926  14.282  1.00 10.00      A    N
ATOM      9  D1  MET A   1      11.949  13.331  15.067  1.00 10.00      A    D
ATOM     10  D2  MET A   1      11.756  12.090  14.335  1.00 10.00      A    D
ATOM     11  D3  MET A   1      12.921  12.918  14.068  1.00 10.00      A    D
ATOM     12  DA  MET A   1      10.401  13.848  13.565  1.00 10.00      A    D
ATOM     13  DB2 MET A   1      12.870  14.826  12.608  1.00 10.00      A    D
ATOM     14  DB3 MET A   1      11.461  15.457  12.240  1.00 10.00      A    D
ATOM     15  DG2 MET A   1      11.184  16.097  14.473  1.00 10.00      A    D
ATOM     16  DG3 MET A   1      12.621  15.508  14.806  1.00 10.00      A    D
ATOM     17  DE1 MET A   1      15.002  17.807  13.062  1.00 10.00      A    D
ATOM     18  DE2 MET A   1      14.854  16.470  13.909  1.00 10.00      A    D
ATOM     19  DE3 MET A   1      14.391  16.504  12.389  1.00 10.00      A    D
TER      20      MET A   1
ATOM     23  O   HOH S   2      16.264  14.181  13.359  1.00 10.00      A    O
ATOM     24  D1  HOH S   2      16.839  14.787  13.448  1.00 10.00      A    D
ATOM     25  D2  HOH S   2      15.607  14.525  12.964  1.00 10.00      A    D
ATOM     26  O   HOH S   3      12.726  17.501  10.452  1.00 10.00      A    O
ATOM     27  D1  HOH S   3      13.400  18.002  10.452  1.00 10.00      A    D
ATOM     28  D2  HOH S   3      12.051  18.001  10.452  1.00 10.00      A    D
ATOM     29  O   HOH S   4      15.741  14.515  10.999  1.00 10.00      A    O
ATOM     30  D1  HOH S   4      16.265  14.645  10.355  1.00 10.00      A    D
ATOM     31  D2  HOH S   4      15.049  14.965  10.844  1.00 10.00      A    D
ATOM     32  O   HOH S   5      14.953  14.255  15.525  1.00 10.00      A    O
ATOM     33  D1  HOH S   5      15.650  14.629  15.243  1.00 10.00      A    D
ATOM     34  D2  HOH S   5      14.518  14.838  15.946  1.00 10.00      A    D
ATOM     35  O   HOH S   6      16.821  13.421  16.856  1.00 10.00      A    O
ATOM     36  D1  HOH S   6      17.631  13.624  16.944  1.00 10.00      A    D
ATOM     37  D2  HOH S   6      16.377  14.028  17.230  1.00 10.00      A    D
ATOM     38  O   HOH S   7      14.155  17.061  16.899  1.00 10.00      A    O
ATOM     39  D1  HOH S   7      14.829  17.562  16.899  1.00 10.00      A    D
ATOM     40  D2  HOH S   7      13.480  17.561  16.899  1.00 10.00      A    D
ATOM     41  O   HOH S   8      13.608  20.138  14.103  1.00 10.00      A    O
ATOM     42  D1  HOH S   8      14.282  20.639  14.103  1.00 10.00      A    D
ATOM     43  D2  HOH S   8      12.933  20.638  14.103  1.00 10.00      A    D
ATOM     44  O   HOH S   9      19.039  15.775  13.706  1.00 10.00      A    O
ATOM     45  D1  HOH S   9      19.709  16.282  13.708  1.00 10.00      A    D
ATOM     46  D2  HOH S   9      18.362  16.267  13.776  1.00 10.00      A    D
ATOM     59  O   HOH S  14      14.398  12.604  10.286  1.00 10.00      A    O
ATOM     60  D1  HOH S  14      15.141  12.872  10.000  1.00 10.00      A    D
ATOM     61  D2  HOH S  14      13.834  13.195  10.088  1.00 10.00      A    D
ATOM     62  O   HOH S  15      14.467  12.610  12.855  1.00 10.00      A    O
ATOM     63  D1  HOH S  15      15.185  12.287  12.563  1.00 10.00      A    D
ATOM     64  D2  HOH S  15      14.064  12.942  12.196  1.00 10.00      A    D
ATOM     65  O   HOH S  16      14.507  11.895  15.986  1.00 10.00      A    O
ATOM     66  D1  HOH S  16      15.110  12.251  16.450  1.00 10.00      A    D
ATOM     67  D2  HOH S  16      13.801  12.325  16.133  1.00 10.00      A    D
ATOM     68  O   HOH S  17      13.717  10.000  14.614  1.00 10.00      A    O
ATOM     69  D1  HOH S  17      14.356  10.543  14.556  1.00 10.00      A    D
ATOM     70  D2  HOH S  17      13.019  10.457  14.711  1.00 10.00      A    D
ATOM     71  O   HOH S  18      13.592  14.007  17.533  1.00 10.00      A    O
ATOM     72  D1  HOH S  18      14.264  13.943  18.034  1.00 10.00      A    D
ATOM     73  D2  HOH S  18      12.988  14.378  17.983  1.00 10.00      A    D
ATOM     74  O   HOH S  19      16.563  12.061  14.539  1.00 10.00      A    O
ATOM     75  D1  HOH S  19      17.279  12.499  14.505  1.00 10.00      A    D
ATOM     76  D2  HOH S  19      15.937  12.621  14.553  1.00 10.00      A    D
TER      63      HOH S  19
END
"""

pdb_poor = """\
CRYST1   29.709   30.639   28.034  90.00  90.00  90.00 P 1
ATOM      1  CB  MET A   1      11.971  15.000  12.927  1.00 10.00      A    C
ATOM      2  CG  MET A   1      12.074  15.930  14.126  1.00 10.00      A    C
ATOM      3  SD  MET A   1      12.820  17.522  13.722  1.00 10.00      A    S
ATOM      4  CE  MET A   1      14.460  17.016  13.210  1.00 10.00      A    C
ATOM      5  C   MET A   1      11.191  12.831  11.962  1.00 10.00      A    C
ATOM      6  O   MET A   1      12.202  12.489  11.349  1.00 10.00      A    O
ATOM      7  CA  MET A   1      11.294  13.663  13.235  1.00 10.00      A    C
ATOM      8  N   MET A   1      12.057  12.926  14.282  1.00 10.00      A    N
ATOM      9  D1  MET A   1      11.949  13.331  15.067  1.00 10.00      A    D
ATOM     10  D2  MET A   1      11.756  12.090  14.335  1.00 10.00      A    D
ATOM     11  D3  MET A   1      12.921  12.918  14.068  1.00 10.00      A    D
ATOM     12  DA  MET A   1      10.401  13.848  13.565  1.00 10.00      A    D
ATOM     13  DB2 MET A   1      12.870  14.826  12.608  1.00 10.00      A    D
ATOM     14  DB3 MET A   1      11.461  15.457  12.240  1.00 10.00      A    D
ATOM     15  DG2 MET A   1      11.184  16.097  14.473  1.00 10.00      A    D
ATOM     16  DG3 MET A   1      12.621  15.508  14.806  1.00 10.00      A    D
ATOM     17  DE1 MET A   1      15.002  17.807  13.062  1.00 10.00      A    D
ATOM     18  DE2 MET A   1      14.854  16.470  13.909  1.00 10.00      A    D
ATOM     19  DE3 MET A   1      14.391  16.504  12.389  1.00 10.00      A    D
TER      20      MET A   1
ATOM     23  O   HOH S   2      16.970  14.660  13.630  1.00 10.00      A    O
ATOM     24  D1  HOH S   2      16.521  14.044  13.276  1.00 10.00      A    D
ATOM     25  D2  HOH S   2      16.766  15.370  13.231  1.00 10.00      A    D
ATOM     26  O   HOH S   3      12.348  17.911  10.504  1.00 10.00      A    O
ATOM     27  D1  HOH S   3      13.109  18.034  10.836  1.00 10.00      A    D
ATOM     28  D2  HOH S   3      11.846  17.664  11.131  1.00 10.00      A    D
ATOM     29  O   HOH S   4      16.209  14.723  11.119  1.00 10.00      A    O
ATOM     30  D1  HOH S   4      16.097  13.893  11.190  1.00 10.00      A    D
ATOM     31  D2  HOH S   4      15.524  15.028  10.742  1.00 10.00      A    D
ATOM     32  O   HOH S   5      15.249  14.287  16.163  1.00 10.00      A    O
ATOM     33  D1  HOH S   5      15.149  14.415  15.339  1.00 10.00      A    D
ATOM     34  D2  HOH S   5      15.157  15.034  16.537  1.00 10.00      A    D
ATOM     35  O   HOH S   6      17.607  13.216  16.419  1.00 10.00      A    O
ATOM     36  D1  HOH S   6      17.291  12.646  16.949  1.00 10.00      A    D
ATOM     37  D2  HOH S   6      17.622  13.948  16.831  1.00 10.00      A    D
ATOM     38  O   HOH S   7      13.968  17.831  16.964  1.00 10.00      A    O
ATOM     39  D1  HOH S   7      14.750  17.530  17.006  1.00 10.00      A    D
ATOM     40  D2  HOH S   7      13.867  18.338  17.627  1.00 10.00      A    D
ATOM     41  O   HOH S   8      13.431  20.300  13.702  1.00 10.00      A    O
ATOM     42  D1  HOH S   8      13.991  20.901  13.877  1.00 10.00      A    D
ATOM     43  D2  HOH S   8      12.704  20.536  14.050  1.00 10.00      A    D
ATOM     44  O   HOH S   9      19.173  15.952  13.349  1.00 10.00      A    O
ATOM     45  D1  HOH S   9      18.824  16.704  13.212  1.00 10.00      A    D
ATOM     46  D2  HOH S   9      18.730  15.586  13.961  1.00 10.00      A    D
ATOM     59  O   HOH S  14      14.246  13.051   9.708  1.00 10.00      A    O
ATOM     60  D1  HOH S  14      15.046  12.794   9.718  1.00 10.00      A    D
ATOM     61  D2  HOH S  14      13.805  12.489  10.150  1.00 10.00      A    D
ATOM     62  O   HOH S  15      14.459  12.209  12.550  1.00 10.00      A    O
ATOM     63  D1  HOH S  15      15.250  12.283  12.277  1.00 10.00      A    D
ATOM     64  D2  HOH S  15      14.028  12.850  12.221  1.00 10.00      A    D
ATOM     65  O   HOH S  16      13.958  12.275  16.789  1.00 10.00      A    O
ATOM     66  D1  HOH S  16      14.738  12.322  17.099  1.00 10.00      A    D
ATOM     67  D2  HOH S  16      13.977  12.633  16.030  1.00 10.00      A    D
ATOM     68  O   HOH S  17      13.943  10.162  15.404  1.00 10.00      A    O
ATOM     69  D1  HOH S  17      14.704  10.494  15.279  1.00 10.00      A    D
ATOM     70  D2  HOH S  17      13.509  10.274  14.694  1.00 10.00      A    D
ATOM     71  O   HOH S  18      13.082  14.321  17.799  1.00 10.00      A    O
ATOM     72  D1  HOH S  18      13.768  14.252  18.278  1.00 10.00      A    D
ATOM     73  D2  HOH S  18      12.594  13.660  17.976  1.00 10.00      A    D
ATOM     74  O   HOH S  19      16.233  12.348  14.233  1.00 10.00      A    O
ATOM     75  D1  HOH S  19      16.721  12.270  14.913  1.00 10.00      A    D
ATOM     76  D2  HOH S  19      15.438  12.353  14.504  1.00 10.00      A    D
TER      66      HOH S  19
END
"""

def exercise(d_min = 1.0, resolution_factor=0.2):
  # Fit HOH or DOD into density map
  #
  # answer model and map
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_answer)
  pdb_inp.write_pdb_file(file_name = "answer.pdb")
  xrs_answer = pdb_inp.xray_structure_simple()
  xrs_answer.switch_to_neutron_scattering_dictionary()
  f_calc_answer = xrs_answer.structure_factors(d_min = d_min).f_calc()
  fft_map = f_calc_answer.fft_map(resolution_factor=resolution_factor)
  fft_map.apply_volume_scaling()
  target_map = fft_map.real_map_unpadded()
  mtz_dataset = f_calc_answer.as_mtz_dataset(column_root_label = "FCmap")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = "answer.mtz")
  # poor model (DOD randomly shifted)
  mon_lib_srv = monomer_library.server.server()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv              = mon_lib_srv,
    ener_lib                 = monomer_library.server.ener_lib(),
    raw_records              = flex.std_string(pdb_poor.splitlines()),
    strict_conflict_handling = True,
    force_symmetry           = True,
    log                      = None)
  pdb_hierarchy_poor = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  xrs_poor = processed_pdb_file.xray_structure()
  xrs_poor.switch_to_neutron_scattering_dictionary()
  pdb_hierarchy_poor.write_pdb_file(file_name = "poor.pdb")
  # real-space refine DOD to fit target_map best
  mmtbx.refinement.real_space.fit_water.run(
    pdb_hierarchy = pdb_hierarchy_poor,
    target_map    = target_map,
    unit_cell     = xrs_poor.unit_cell(),
    real_space_gradients_delta = d_min*resolution_factor,
    log = sys.stdout)
  # best refined model
  pdb_hierarchy_poor.write_pdb_file(file_name = "refined.pdb")
  xrs_refined = pdb_hierarchy_poor.extract_xray_structure(
    crystal_symmetry = xrs_poor.crystal_symmetry())
  xrs_refined.switch_to_neutron_scattering_dictionary()
  # check results
  def r_factor(x,y):
    x = abs(x).data()
    y = abs(y).data()
    n = flex.sum(flex.abs(x-y))
    d = flex.sum(flex.abs(x+y))/2
    return n/d
  fc_poor = f_calc_answer.structure_factors_from_scatterers(
    xray_structure = xrs_poor).f_calc()
  fc_refined = f_calc_answer.structure_factors_from_scatterers(
    xray_structure = xrs_refined).f_calc()
  assert r_factor(f_calc_answer, fc_poor) > 0.45
  assert r_factor(f_calc_answer, fc_refined) < 0.08

if(__name__ == "__main__"):
  t0 = time.time()
  exercise()
  print("Time: %6.4f"%(time.time()-t0))
