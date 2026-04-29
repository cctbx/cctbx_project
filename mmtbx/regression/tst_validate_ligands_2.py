from __future__ import absolute_import, division, print_function
import time, traceback, os
import libtbx.load_env
from libtbx.utils import null_out
import mmtbx.model
import iotbx.pdb
from libtbx.test_utils import approx_equal
from iotbx.cli_parser import run_program
from mmtbx.programs import validate_ligands as val_lig
from mmtbx.regression.tst_validate_ligands import find_lr

from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL) # Only show critical errors

# ------------------------------------------------------------------------------

def run():
  run_test01()
  run_test_get_results_fallback()

# ------------------------------------------------------------------------------

def run_test01():
  '''
  Several tests for AQS (intercalator) in a DNA fragment (386D):
    - occupancy and ADP metrics for AQS
    - overlap metrics for AQS
    - DT is excluded because it is further than 3 A (even sym copies)
    - sites: some dna residues invole symmetry
    - sites: don't count HOH or EDO or sym copy of the ligand itself
  '''
  print('test01')
  model_fn = "tst_01_fragment.pdb"
  with open(model_fn, "w") as f:
    f.write(pdb_str_tst_01)
  args = [model_fn, 'run_reduce2=False']
  try:
    result = run_program(program_class=val_lig.Program,args=args,
     logger = null_out())
  except Exception as e:
    msg = traceback.format_exc()
  vl_manager = result.ligand_manager

  assert len(vl_manager) == 2

  # --- AQS A 7 ---
  lr = find_lr(vl_manager, 'chain A and resseq 7 and resname AQS')

  # get_occupancies(): all AQS atoms have occ=1.00
  occs = lr.get_occupancies()
  assert approx_equal(occs.occ_mean, 1.0, eps=0.01)
  assert occs.zero_count == 0
  assert occs.negative_count == 0

  # get_adps(): all AQS atoms have isotropic B factors
  adps = lr.get_adps()
  assert adps.n_aniso == 0
  assert adps.n_iso > 0
  assert approx_equal(adps.b_min, 13.12, eps=0.01)
  assert approx_equal(adps.b_max, 40.19, eps=0.01)
  assert approx_equal(adps.b_mean, 19.18, eps=0.05)
  assert adps.isel_within_noH.size() > 0
  assert approx_equal(adps.b_min_within, 6.90, eps=0.01)
  assert approx_equal(adps.b_max_within, 38.24, eps=0.01)
  assert approx_equal(adps.b_mean_within, 18.60, eps=0.02)

  # isel_within_noH: verify neighboring sites composition
  sites_model = lr.model.select(adps.isel_within_noH)
  sites_hierarchy = sites_model.get_hierarchy()
  resnames_in_sites = []
  for rg in sites_hierarchy.residue_groups():
    resnames_in_sites.extend(rg.unique_resnames())
  resnames_in_sites = [resname.strip() for resname in resnames_in_sites]
  # DT is further than 3A from AQS -> correctly excluded
  assert 'DT' not in resnames_in_sites
  # HOH and other ligands are excluded
  assert 'EDO' not in resnames_in_sites
  assert 'HOH' not in resnames_in_sites
  # Exact DNA neighbor set: DA, DC, DG are within 3.0 A of AQS; DT is not
  assert set(resnames_in_sites) == {'DA', 'DC', 'DG'}

  overlaps = lr.get_overlaps()
  assert overlaps is not None
  assert overlaps.n_clashes == 9
  assert overlaps.n_hbonds == 3
  assert approx_equal(overlaps.clashscore, 40.0, eps=0.5)

  # --- EDO A 43 ---
  lr_edo = find_lr(vl_manager, 'chain A and resseq 43 and resname EDO')
  assert lr_edo.resname == 'EDO'
  assert lr_edo.altloc == ''
  occs_edo = lr_edo.get_occupancies()
  assert approx_equal(occs_edo.occ_mean, 1.0, eps=0.01)
  assert occs_edo.zero_count == 0
  assert occs_edo.negative_count == 0

# ------------------------------------------------------------------------------

def run_test_get_results_fallback():
  '''
  When run_reduce2=False, get_results() must return the original model
  filename (not a _newH.cif path that was never written to disk).
  '''
  from libtbx import env
  pdb_fname = env.find_in_repositories(
    relative_path="phenix_regression/pdb/pdb1avd.ent.gz",
    test=os.path.isfile)
  if pdb_fname is None:
    print('phenix_regression not available, skipping')
    return
  result = run_program(
    program_class=val_lig.Program,
    args=[pdb_fname, 'run_reduce2=False'],
    logger=null_out())
  assert result.working_model_fn is not None, \
    "working_model_fn should not be None when run_reduce2=False"
  assert result.working_model_fn == pdb_fname, \
    "working_model_fn=%s expected=%s" % (result.working_model_fn, pdb_fname)
  assert os.path.isfile(result.working_model_fn), \
    "working_model_fn path does not exist on disk: %s" % result.working_model_fn
  print('OK: get_results fallback to original model when run_reduce2=False')

# ------------------------------------------------------------------------------

pdb_str_tst_01 = '''
REMARK from 386D, edited in Coot to add EDO and two HOH to test specific
REMARK scenarios
CRYST1   31.020   31.020   64.910  90.00  90.00 120.00 P 65 2 2
ATOM      1  O5'  DC A   1       7.131   3.724  25.592  1.00 32.17           O
ATOM      2  C5'  DC A   1       5.796   3.577  25.048  1.00 30.73           C
ATOM      3  C4'  DC A   1       5.320   4.758  24.148  1.00 26.68           C
ATOM      4  O4'  DC A   1       5.469   6.029  24.820  1.00 24.91           O
ATOM      5  C3'  DC A   1       6.099   4.837  22.827  1.00 25.57           C
ATOM      6  O3'  DC A   1       5.217   4.948  21.690  1.00 24.42           O
ATOM      7  C2'  DC A   1       6.989   6.044  23.005  1.00 21.10           C
ATOM      8  C1'  DC A   1       6.303   6.894  24.053  1.00 20.95           C
ATOM      9  N1   DC A   1       7.316   7.575  24.901  1.00 20.49           N
ATOM     10  C2   DC A   1       7.388   8.967  24.821  1.00 17.60           C
ATOM     11  O2   DC A   1       6.641   9.643  24.105  1.00 16.14           O
ATOM     12  N3   DC A   1       8.340   9.593  25.551  1.00 14.96           N
ATOM     13  C4   DC A   1       9.170   8.910  26.333  1.00 17.05           C
ATOM     14  N4   DC A   1      10.053   9.625  27.014  1.00 14.59           N
ATOM     15  C5   DC A   1       9.115   7.476  26.445  1.00 18.40           C
ATOM     16  C6   DC A   1       8.170   6.860  25.709  1.00 18.05           C
ATOM     17  H4*  DC A   1       4.370   4.680  23.970  1.00 26.68           H
ATOM     18  H3*  DC A   1       6.647   4.042  22.732  1.00 25.57           H
ATOM     19  H1*  DC A   1       5.704   7.549  23.662  1.00 20.95           H
ATOM     20  H41  DC A   1      10.066  10.481  26.935  1.00 14.59           H
ATOM     21  H42  DC A   1      10.614   9.233  27.535  1.00 14.59           H
ATOM     22  H5   DC A   1       9.698   7.002  26.994  1.00 18.40           H
ATOM     23 H5*1  DC A   1       5.175   3.491  25.788  1.00 30.73           H
ATOM     24 H5*2  DC A   1       5.772   2.764  24.519  1.00 30.73           H
ATOM     25  H6   DC A   1       8.093   5.934  25.748  1.00 18.05           H
ATOM     26 H2*1  DC A   1       7.872   5.777  23.306  1.00 21.10           H
ATOM     27 H2*2  DC A   1       7.074   6.532  22.171  1.00 21.10           H
ATOM     28  P    DG A   2       5.526   4.262  20.255  1.00 25.54           P
ATOM     29  OP1  DG A   2       4.247   4.029  19.558  1.00 26.58           O
ATOM     30  OP2  DG A   2       6.461   3.137  20.450  1.00 28.36           O
ATOM     31  O5'  DG A   2       6.321   5.391  19.492  1.00 25.76           O
ATOM     32  C5'  DG A   2       5.647   6.627  19.322  1.00 24.81           C
ATOM     33  C4'  DG A   2       6.542   7.672  18.709  1.00 23.30           C
ATOM     34  O4'  DG A   2       7.674   8.007  19.530  1.00 23.12           O
ATOM     35  C3'  DG A   2       6.982   7.195  17.329  1.00 20.87           C
ATOM     36  O3'  DG A   2       6.174   7.918  16.421  1.00 19.21           O
ATOM     37  C2'  DG A   2       8.413   7.588  17.325  1.00 21.26           C
ATOM     38  C1'  DG A   2       8.629   8.426  18.580  1.00 21.19           C
ATOM     39  N9   DG A   2       9.963   8.173  19.124  1.00 19.16           N
ATOM     40  C8   DG A   2      10.624   6.994  19.388  1.00 16.10           C
ATOM     41  N7   DG A   2      11.833   7.180  19.851  1.00 17.26           N
ATOM     42  C5   DG A   2      11.956   8.568  19.881  1.00 15.70           C
ATOM     43  C6   DG A   2      13.013   9.382  20.297  1.00 16.71           C
ATOM     44  O6   DG A   2      14.099   9.051  20.744  1.00 17.49           O
ATOM     45  N1   DG A   2      12.696  10.715  20.152  1.00 16.35           N
ATOM     46  C2   DG A   2      11.516  11.199  19.673  1.00 12.93           C
ATOM     47  N2   DG A   2      11.347  12.495  19.613  1.00 15.24           N
ATOM     48  N3   DG A   2      10.521  10.464  19.285  1.00 12.07           N
ATOM     49  C4   DG A   2      10.813   9.163  19.428  1.00 16.07           C
ATOM     50  H4*  DG A   2       6.119   8.544  18.657  1.00 23.30           H
ATOM     51  H3*  DG A   2       6.932   6.239  17.175  1.00 20.87           H
ATOM     52  H1*  DG A   2       8.480   9.369  18.406  1.00 21.19           H
ATOM     53  H8   DG A   2      10.250   6.154  19.251  1.00 16.10           H
ATOM     54 H5*1  DG A   2       5.341   6.941  20.187  1.00 24.81           H
ATOM     55 H5*2  DG A   2       4.879   6.491  18.745  1.00 24.81           H
ATOM     56  H1   DG A   2      13.293  11.289  20.383  1.00 16.35           H
ATOM     57  H21  DG A   2      11.979  13.023  19.860  1.00 15.24           H
ATOM     58  H22  DG A   2      10.604  12.819  19.327  1.00 15.24           H
ATOM     59 H2*1  DG A   2       8.980   6.801  17.342  1.00 21.26           H
ATOM     60 H2*2  DG A   2       8.619   8.108  16.532  1.00 21.26           H
ATOM     61  P    DT A   3       6.408   8.085  14.861  1.00 17.40           P
ATOM     62  OP1  DT A   3       5.074   8.142  14.234  1.00 19.07           O
ATOM     63  OP2  DT A   3       7.429   7.097  14.415  1.00 19.75           O
ATOM     64  O5'  DT A   3       7.033   9.548  14.753  1.00 17.91           O
ATOM     65  C5'  DT A   3       6.252  10.657  15.187  1.00 13.04           C
ATOM     66  C4'  DT A   3       7.052  11.928  15.089  1.00 14.03           C
ATOM     67  O4'  DT A   3       8.273  11.792  15.854  1.00 13.63           O
ATOM     68  C3'  DT A   3       7.436  12.240  13.641  1.00 15.29           C
ATOM     69  O3'  DT A   3       7.116  13.618  13.391  1.00 16.95           O
ATOM     70  C2'  DT A   3       8.936  11.936  13.620  1.00 13.01           C
ATOM     71  C1'  DT A   3       9.388  12.171  15.045  1.00 11.20           C
ATOM     72  N1   DT A   3      10.502  11.322  15.487  1.00  9.48           N
ATOM     73  C2   DT A   3      11.544  11.946  16.089  1.00 11.98           C
ATOM     74  O2   DT A   3      11.610  13.149  16.237  1.00 12.14           O
ATOM     75  N3   DT A   3      12.571  11.171  16.527  1.00 12.38           N
ATOM     76  C4   DT A   3      12.649   9.824  16.428  1.00 11.62           C
ATOM     77  O4   DT A   3      13.622   9.256  16.888  1.00 14.06           O
ATOM     78  C5   DT A   3      11.526   9.213  15.788  1.00 13.71           C
ATOM     79  C7   DT A   3      11.508   7.695  15.613  1.00 11.79           C
ATOM     80  C6   DT A   3      10.498   9.965  15.354  1.00 12.71           C
ATOM     81  H4*  DT A   3       6.580  12.683  15.473  1.00 14.03           H
ATOM     82  H3*  DT A   3       6.998  11.662  12.997  1.00 15.29           H
ATOM     83  H1*  DT A   3       9.625  13.106  15.144  1.00 11.20           H
ATOM     84 H5*1  DT A   3       5.979  10.518  16.107  1.00 13.04           H
ATOM     85 H5*2  DT A   3       5.461  10.731  14.631  1.00 13.04           H
ATOM     86  H3   DT A   3      13.232  11.574  16.902  1.00 12.38           H
ATOM     87 H5M1  DT A   3      12.420   7.376  15.526  1.00 11.79           H
ATOM     88 H5M2  DT A   3      11.004   7.475  14.814  1.00 11.79           H
ATOM     89 H5M3  DT A   3      11.089   7.292  16.390  1.00 11.79           H
ATOM     90  H6   DT A   3       9.768   9.550  14.954  1.00 12.71           H
ATOM     91 H2*1  DT A   3       9.098  11.018  13.353  1.00 13.01           H
ATOM     92 H2*2  DT A   3       9.400  12.529  13.008  1.00 13.01           H
ATOM     93  P    DA A   4       7.364  14.327  11.996  1.00 18.25           P
ATOM     94  OP1  DA A   4       6.446  15.485  11.860  1.00 20.86           O
ATOM     95  OP2  DA A   4       7.421  13.256  10.967  1.00 20.85           O
ATOM     96  O5'  DA A   4       8.837  14.890  12.208  1.00 14.98           O
ATOM     97  C5'  DA A   4       9.049  15.880  13.194  1.00 11.64           C
ATOM     98  C4'  DA A   4      10.477  16.382  13.167  1.00 15.41           C
ATOM     99  O4'  DA A   4      11.407  15.411  13.650  1.00 11.43           O
ATOM    100  C3'  DA A   4      10.867  16.747  11.748  1.00 14.39           C
ATOM    101  O3'  DA A   4      11.365  18.097  11.802  1.00 18.56           O
ATOM    102  C2'  DA A   4      11.876  15.657  11.436  1.00 14.26           C
ATOM    103  C1'  DA A   4      12.510  15.325  12.766  1.00 11.86           C
ATOM    104  N9   DA A   4      12.916  13.932  12.926  1.00 14.49           N
ATOM    105  C8   DA A   4      12.254  12.801  12.528  1.00 13.30           C
ATOM    106  N7   DA A   4      12.873  11.697  12.827  1.00 15.15           N
ATOM    107  C5   DA A   4      14.019  12.126  13.489  1.00 11.48           C
ATOM    108  C6   DA A   4      15.090  11.443  14.064  1.00  8.19           C
ATOM    109  N6   DA A   4      15.186  10.126  14.103  1.00  7.41           N
ATOM    110  N1   DA A   4      16.075  12.166  14.589  1.00  6.90           N
ATOM    111  C2   DA A   4      15.973  13.482  14.586  1.00 11.11           C
ATOM    112  N3   DA A   4      15.022  14.252  14.094  1.00 13.38           N
ATOM    113  C4   DA A   4      14.051  13.489  13.555  1.00 13.68           C
ATOM    114  H4*  DA A   4      10.602  17.125  13.778  1.00 15.41           H
ATOM    115  H3*  DA A   4      10.165  16.685  11.081  1.00 14.39           H
ATOM    116  H1*  DA A   4      13.253  15.927  12.926  1.00 11.86           H
ATOM    117  H8   DA A   4      11.437  12.821  12.085  1.00 13.30           H
ATOM    118 H5*1  DA A   4       8.860  15.504  14.068  1.00 11.64           H
ATOM    119 H5*2  DA A   4       8.447  16.623  13.033  1.00 11.64           H
ATOM    120  H61  DA A   4      15.874   9.755  14.462  1.00  7.41           H
ATOM    121  H62  DA A   4      14.560   9.639  13.769  1.00  7.41           H
ATOM    122  H2   DA A   4      16.677  13.934  14.992  1.00 11.11           H
ATOM    123 H2*1  DA A   4      11.436  14.878  11.061  1.00 14.26           H
ATOM    124 H2*2  DA A   4      12.543  15.976  10.808  1.00 14.26           H
ATOM    125  P    DC A   5      11.902  18.930  10.533  1.00 23.98           P
ATOM    126  OP1  DC A   5      12.006  20.359  10.906  1.00 22.46           O
ATOM    127  OP2  DC A   5      11.142  18.529   9.327  1.00 20.57           O
ATOM    128  O5'  DC A   5      13.391  18.392  10.353  1.00 23.11           O
ATOM    129  C5'  DC A   5      14.382  18.788  11.307  1.00 22.16           C
ATOM    130  C4'  DC A   5      15.659  18.015  11.079  1.00 22.14           C
ATOM    131  O4'  DC A   5      15.402  16.642  11.363  1.00 18.16           O
ATOM    132  C3'  DC A   5      16.081  18.121   9.601  1.00 25.05           C
ATOM    133  O3'  DC A   5      17.439  18.553   9.451  1.00 32.54           O
ATOM    134  C2'  DC A   5      15.926  16.697   9.172  1.00 22.88           C
ATOM    135  C1'  DC A   5      16.187  15.953  10.456  1.00 16.78           C
ATOM    136  N1   DC A   5      15.867  14.528  10.446  1.00 15.07           N
ATOM    137  C2   DC A   5      16.683  13.686  11.183  1.00 13.67           C
ATOM    138  O2   DC A   5      17.627  14.123  11.836  1.00 15.16           O
ATOM    139  N3   DC A   5      16.436  12.352  11.173  1.00 11.27           N
ATOM    140  C4   DC A   5      15.416  11.870  10.474  1.00 15.90           C
ATOM    141  N4   DC A   5      15.192  10.569  10.511  1.00 12.59           N
ATOM    142  C5   DC A   5      14.544  12.718   9.726  1.00 14.42           C
ATOM    143  C6   DC A   5      14.809  14.040   9.751  1.00 17.89           C
ATOM    144  H4*  DC A   5      16.388  18.287  11.659  1.00 22.14           H
ATOM    145  H3*  DC A   5      15.495  18.731   9.126  1.00 25.05           H
ATOM    146  H1*  DC A   5      17.124  15.992  10.705  1.00 16.78           H
ATOM    147  H41  DC A   5      15.701  10.058  10.980  1.00 12.59           H
ATOM    148  H42  DC A   5      14.537  10.233  10.067  1.00 12.59           H
ATOM    149  H5   DC A   5      13.827  12.375   9.244  1.00 14.42           H
ATOM    150 H5*1  DC A   5      14.055  18.614  12.203  1.00 22.16           H
ATOM    151 H5*2  DC A   5      14.558  19.737  11.214  1.00 22.16           H
ATOM    152  H6   DC A   5      14.258  14.626   9.285  1.00 17.89           H
ATOM    153 H2*1  DC A   5      15.033  16.523   8.837  1.00 22.88           H
ATOM    154 H2*2  DC A   5      16.572  16.461   8.488  1.00 22.88           H
ATOM    155  P    DG A   6      18.025  19.264   8.114  1.00 35.70           P
ATOM    156  OP1  DG A   6      17.884  20.721   8.311  1.00 38.24           O
ATOM    157  OP2  DG A   6      17.527  18.613   6.883  1.00 30.14           O
ATOM    158  O5'  DG A   6      19.552  18.878   8.241  1.00 33.46           O
ATOM    159  C5'  DG A   6      20.457  19.425   9.203  1.00 28.71           C
ATOM    160  C4'  DG A   6      21.892  19.172   8.715  1.00 26.64           C
ATOM    161  O4'  DG A   6      22.297  17.784   8.690  1.00 22.80           O
ATOM    162  C3'  DG A   6      21.902  19.642   7.277  1.00 23.84           C
ATOM    163  O3'  DG A   6      22.061  21.059   7.281  1.00 27.15           O
ATOM    164  C2'  DG A   6      22.999  18.816   6.667  1.00 19.30           C
ATOM    165  C1'  DG A   6      23.000  17.512   7.470  1.00 16.41           C
ATOM    166  N9   DG A   6      22.292  16.384   6.879  1.00 10.92           N
ATOM    167  C8   DG A   6      21.175  16.408   6.094  1.00 10.19           C
ATOM    168  N7   DG A   6      20.764  15.235   5.748  1.00  8.81           N
ATOM    169  C5   DG A   6      21.657  14.369   6.357  1.00 11.88           C
ATOM    170  C6   DG A   6      21.689  12.951   6.348  1.00 12.55           C
ATOM    171  O6   DG A   6      20.924  12.155   5.789  1.00 13.14           O
ATOM    172  N1   DG A   6      22.733  12.468   7.085  1.00 10.66           N
ATOM    173  C2   DG A   6      23.647  13.239   7.743  1.00 13.46           C
ATOM    174  N2   DG A   6      24.597  12.572   8.378  1.00 13.95           N
ATOM    175  N3   DG A   6      23.638  14.572   7.761  1.00 11.92           N
ATOM    176  C4   DG A   6      22.598  15.065   7.058  1.00  8.35           C
ATOM    177  H4*  DG A   6      22.556  19.554   9.310  1.00 26.64           H
ATOM    178  H3*  DG A   6      21.100  19.423   6.777  1.00 23.84           H
ATOM    179  H1*  DG A   6      23.931  17.281   7.611  1.00 16.41           H
ATOM    180  H8   DG A   6      20.751  17.194   5.834  1.00 10.19           H
ATOM    181 H5*1  DG A   6      20.321  18.998  10.063  1.00 28.71           H
ATOM    182 H5*2  DG A   6      20.303  20.378   9.296  1.00 28.71           H
ATOM    183  H1   DG A   6      22.821  11.614   7.138  1.00 10.66           H
ATOM    184  H21  DG A   6      24.601  11.712   8.365  1.00 13.95           H
ATOM    185  H22  DG A   6      25.211  12.998   8.804  1.00 13.95           H
ATOM    186 H2*1  DG A   6      22.818  18.642   5.730  1.00 19.30           H
ATOM    187 H2*2  DG A   6      23.854  19.268   6.745  1.00 19.30           H
TER
HETATM  188  C1  AQS A   7       9.719   8.942  22.465  1.00 13.67           C
HETATM  189  C16 AQS A   7      10.649  15.148  22.810  1.00 14.63           C
HETATM  190  C17 AQS A   7      10.002  13.963  22.422  1.00 16.51           C
HETATM  191  C18 AQS A   7      10.506  12.717  22.850  1.00 15.60           C
HETATM  192  C19 AQS A   7       9.806  11.448  22.429  1.00 16.51           C
HETATM  193  C1A AQS A   7      12.886   4.698  21.885  1.00 27.18           C
HETATM  194  C1B AQS A   7      14.462   6.333  22.815  1.00 20.26           C
HETATM  195  C2  AQS A   7      10.220   7.693  22.858  1.00 16.06           C
HETATM  196  C20 AQS A   7      10.374  10.121  22.863  1.00 15.37           C
HETATM  197  C2A AQS A   7      12.995   3.182  22.111  1.00 32.13           C
HETATM  198  C2B AQS A   7      15.789   5.790  23.395  1.00 19.86           C
HETATM  199  C3  AQS A   7      11.355   7.634  23.675  1.00 17.39           C
HETATM  200  C3A AQS A   7      11.918   2.223  21.589  1.00 36.72           C
HETATM  201  C3B AQS A   7      15.776   4.696  24.500  1.00 20.09           C
HETATM  202  C4  AQS A   7      12.015   8.806  24.087  1.00 13.12           C
HETATM  203  C5  AQS A   7      11.530  10.057  23.675  1.00 14.50           C
HETATM  204  C6  AQS A   7      12.219  11.319  24.116  1.00 14.61           C
HETATM  205  C7  AQS A   7      11.662  12.657  23.671  1.00 15.41           C
HETATM  206  C8  AQS A   7      12.290  13.849  24.083  1.00 15.49           C
HETATM  207  C9  AQS A   7      11.779  15.085  23.646  1.00 15.94           C
HETATM  208  N1  AQS A   7      13.190   5.602  23.005  1.00 22.91           N
HETATM  209  N3A AQS A   7      12.256   1.247  20.541  1.00 40.19           N
HETATM  210  N3B AQS A   7      17.002   4.021  24.966  1.00 14.84           N
HETATM  211  O19 AQS A   7       8.816  11.479  21.720  1.00 15.69           O
HETATM  212  O1S AQS A   7      12.567   6.182  25.473  1.00 18.76           O
HETATM  213  O2S AQS A   7      10.888   5.091  24.169  1.00 15.16           O
HETATM  214  O6  AQS A   7      13.206  11.264  24.826  1.00 17.92           O
HETATM  215  S1  AQS A   7      11.959   6.054  24.175  1.00 20.48           S
HETATM  216  H1  AQS A   7       8.952   8.986  21.940  1.00 13.67           H
HETATM  217  H2  AQS A   7       9.793   6.918  22.571  1.00 16.06           H
HETATM  218 HA11 AQS A   7      11.987   4.921  21.597  1.00 27.18           H
HETATM  219 HA12 AQS A   7      13.474   4.962  21.160  1.00 27.18           H
HETATM  220 HA21 AQS A   7      13.054   3.036  23.068  1.00 32.13           H
HETATM  221 HA22 AQS A   7      13.835   2.892  21.722  1.00 32.13           H
HETATM  222 HA31 AQS A   7      11.180   2.726  21.211  1.00 36.72           H
HETATM  223 HA32 AQS A   7      11.578   1.678  22.316  1.00 36.72           H
HETATM  224 HB11 AQS A   7      14.343   7.228  23.170  1.00 20.26           H
HETATM  225 HB12 AQS A   7      14.596   6.448  21.861  1.00 20.26           H
HETATM  226 HB21 AQS A   7      16.273   6.552  23.750  1.00 19.86           H
HETATM  227 HB22 AQS A   7      16.299   5.439  22.648  1.00 19.86           H
HETATM  228 HB31 AQS A   7      15.409   5.065  25.319  1.00 20.09           H
HETATM  229 HB32 AQS A   7      15.216   3.955  24.218  1.00 20.09           H
HETATM  230  H16 AQS A   7      10.331  15.970  22.515  1.00 14.63           H
HETATM  231  H17 AQS A   7       9.244  13.999  21.885  1.00 16.51           H
HETATM  232  H4  AQS A   7      12.769   8.765  24.629  1.00 13.12           H
HETATM  233  H8  AQS A   7      13.035  13.818  24.639  1.00 15.49           H
HETATM  234  H9  AQS A   7      12.196  15.871  23.915  1.00 15.94           H
HETATM  235 HNA1 AQS A   7      11.531   0.874  20.184  1.00 40.19           H
HETATM  236 HNA2 AQS A   7      12.759   0.582  20.853  1.00 40.19           H
HETATM  237 HNA3 AQS A   7      12.709   1.619  19.871  1.00 40.19           H
HETATM  238 HNB1 AQS A   7      16.831   3.421  25.601  1.00 14.84           H
HETATM  239 HNB2 AQS A   7      17.590   4.596  25.306  1.00 14.84           H
HETATM  240 HNB3 AQS A   7      17.415   3.589  24.307  1.00 14.84           H
HETATM  241  O   HOH A  21      16.983   9.805  27.044  0.50 31.01           O
HETATM  242  O   HOH A  22       5.351   9.268  21.635  0.50 25.40           O
HETATM  243  O   HOH A  23      22.145  -4.498  25.629  1.00 30.11           O
HETATM  244  O   HOH A  24      10.868  14.023   7.319  1.00 51.27           O
HETATM  245  O   HOH A  25      12.480  -1.766  21.445  1.00 27.04           O
HETATM  246  O   HOH A  26      12.713   4.345  15.922  1.00 27.20           O
HETATM  247  O   HOH A  27       9.823   2.873   8.807  1.00 38.36           O
HETATM  248  O   HOH A  28      22.345   0.729  26.121  1.00 39.98           O
HETATM  249  O   HOH A  29      10.046  12.524  10.249  1.00 29.29           O
HETATM  250  O   HOH A  30      12.432   9.664   9.319  1.00 45.05           O
HETATM  251  O   HOH A  31      12.053   9.121  12.292  1.00 29.10           O
HETATM  252  O   HOH A  32      14.386   6.495  17.140  1.00 24.20           O
HETATM  253  O   HOH A  33      10.180   3.995  17.932  1.00 41.63           O
HETATM  254  O   HOH A  34       7.156  -1.975   6.531  1.00 31.13           O
HETATM  255  O   HOH A  35      12.648  16.288   7.958  1.00 33.18           O
HETATM  256  O   HOH A  37       9.147   2.516  23.341  1.00 36.30           O
HETATM  257  O   HOH A  38      16.074   6.252  19.738  1.00 46.75           O
HETATM  258  O   HOH A  39      14.808  13.484  31.658  1.00 52.17           O
HETATM  259  O   HOH A  40      19.750  23.133   7.666  1.00 54.35           O
HETATM  260  O   HOH A  41       9.508   4.297  20.950  1.00 40.92           O
HETATM  261  O   HOH A  42       9.535   7.427  12.344  1.00 55.13           O
HETATM  262  C1  EDO A  43      15.183  13.071  27.804  1.00 20.00           C
HETATM  263  C2  EDO A  43      14.343  14.091  27.127  1.00 20.00           C
HETATM  264  O1  EDO A  43      16.506  13.024  27.304  1.00 20.00           O
HETATM  265  O2  EDO A  43      14.352  13.973  25.717  1.00 20.00           O
HETATM  266  H11 EDO A  43      15.184  13.269  28.754  1.00 20.00           H
HETATM  267  H12 EDO A  43      14.748  12.210  27.699  1.00 20.00           H
HETATM  268  H21 EDO A  43      14.671  14.957  27.415  1.00 20.00           H
HETATM  269  H22 EDO A  43      13.447  13.997  27.488  1.00 20.00           H
HETATM  270  HO1 EDO A  43      16.491  12.625  26.554  1.00 20.00           H
HETATM  271  HO2 EDO A  43      14.808  13.285  25.512  1.00 20.00           H
HETATM  272  O   HOH A  44      13.141   7.938  24.804  1.00 30.00           O
HETATM  273  O   HOH A  45      14.315   9.704  25.675  1.00 30.00           O
END
'''

# ------------------------------------------------------------------------------

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f" % (time.time()-t0))
