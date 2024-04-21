from __future__ import absolute_import, division, print_function

from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
import sys
try:
  from six.moves import cPickle as pickle
except ImportError:
  import pickle
from six.moves import cStringIO as StringIO
from libtbx.test_utils import show_diff
import iotbx.pdb
from libtbx.utils import null_out # import dependency
from time import time
from libtbx.test_utils import approx_equal


raw_records1 = """\
CRYST1   60.800   60.800   97.000  90.00  90.00 120.00 P 32 2 1      6
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.016447  0.009496  0.000000        0.00000
SCALE2      0.000000  0.018992  0.000000        0.00000
SCALE3      0.000000  0.000000  0.010309        0.00000
ATOM   1050  N   LYS A 135      31.992  14.930  -7.233  1.00  9.47           N
ATOM   1051  CA  LYS A 135      31.388  16.216  -7.637  1.00 12.89           C
ATOM   1052  C   LYS A 135      30.807  16.840  -6.406  1.00  6.47           C
ATOM   1053  O   LYS A 135      29.583  16.869  -6.191  1.00 15.74           O
ATOM   1054  CB  LYS A 135      30.263  16.059  -8.655  1.00 13.51           C
ATOM   1055  CG  LYS A 135      30.742  15.277  -9.843  1.00 16.23           C
ATOM   1056  CD  LYS A 135      29.612  15.131 -10.835  1.00 28.55           C
ATOM   1057  CE  LYS A 135      30.173  14.812 -12.216  1.00 34.52           C
ATOM   1058  NZ  LYS A 135      29.396  13.756 -12.899  1.00 46.18           N
TER    1294      LYS A 162
END
"""

raw_records2 = """\
CRYST1   24.627   42.717   46.906  90.00  90.00  90.00 P 21 21 21    8
ATOM    184  P    DG A   9       9.587  13.026  19.037  1.00  6.28           P
ATOM    185  OP1  DG A   9       9.944  14.347  19.602  1.00  8.07           O
ATOM    186  OP2  DG A   9      10.654  12.085  18.639  1.00  8.27           O
ATOM    187  O5'  DG A   9       8.717  12.191  20.048  1.00  5.88           O
ATOM    188  C5'  DG A   9       7.723  12.833  20.854  1.00  5.45           C
ATOM    189  C4'  DG A   9       7.145  11.818  21.807  1.00  5.40           C
ATOM    190  O4'  DG A   9       6.435  10.777  21.087  1.00  5.77           O
ATOM    191  C3'  DG A   9       8.142  11.036  22.648  1.00  5.10           C
ATOM    192  O3'  DG A   9       8.612  11.838  23.723  1.00  5.90           O
ATOM    193  C2'  DG A   9       7.300   9.857  23.068  1.00  5.97           C
ATOM    194  C1'  DG A   9       6.619   9.536  21.805  1.00  5.97           C
ATOM    195  N9   DG A   9       7.390   8.643  20.931  1.00  5.97           N
ATOM    196  C8   DG A   9       8.074   8.881  19.775  1.00  6.62           C
ATOM    197  N7   DG A   9       8.647   7.820  19.249  1.00  6.57           N
ATOM    198  C5   DG A   9       8.308   6.806  20.141  1.00  6.22           C
ATOM    199  C6   DG A   9       8.620   5.431  20.136  1.00  6.03           C
ATOM    200  O6   DG A   9       9.297   4.803  19.296  1.00  7.21           O
ATOM    201  N1   DG A   9       8.101   4.773  21.247  1.00  6.10           N
ATOM    202  C2   DG A   9       7.365   5.351  22.260  1.00  6.24           C
ATOM    203  N2   DG A   9       6.948   4.569  23.241  1.00  7.88           N
ATOM    204  N3   DG A   9       7.051   6.652  22.257  1.00  6.53           N
ATOM    205  C4   DG A   9       7.539   7.295  21.184  1.00  5.69           C
ATOM    206  P    DC A  10      10.081  11.538  24.300  1.00  5.91           P
ATOM    207  OP1  DC A  10      10.273  12.645  25.291  1.00  7.27           O
ATOM    208  OP2  DC A  10      11.063  11.363  23.228  1.00  6.84           O
ATOM    209  O5'  DC A  10       9.953  10.128  25.026  1.00  5.75           O
ATOM    210  C5'  DC A  10       9.077   9.959  26.149  1.00  5.87           C
ATOM    211  C4'  DC A  10       9.188   8.549  26.672  1.00  5.56           C
ATOM    212  O4'  DC A  10       8.708   7.612  25.667  1.00  5.70           O
ATOM    213  C3'  DC A  10      10.580   8.059  27.007  1.00  5.27           C
ATOM    214  O3'  DC A  10      11.010   8.447  28.315  1.00  5.83           O
ATOM    215  C2'  DC A  10      10.422   6.549  26.893  1.00  5.34           C
ATOM    216  C1'  DC A  10       9.436   6.405  25.754  1.00  5.23           C
ATOM    217  N1   DC A  10      10.113   6.168  24.448  1.00  5.30           N
ATOM    218  C2   DC A  10      10.514   4.871  24.152  1.00  5.28           C
ATOM    219  O2   DC A  10      10.283   3.972  25.000  1.00  5.75           O
ATOM    220  N3   DC A  10      11.131   4.627  22.965  1.00  5.65           N
ATOM    221  C4   DC A  10      11.395   5.628  22.138  1.00  5.80           C
ATOM    222  N4   DC A  10      12.034   5.327  21.005  1.00  6.75           N
ATOM    223  C5   DC A  10      11.029   6.970  22.449  1.00  5.99           C
ATOM    224  C6   DC A  10      10.394   7.203  23.612  1.00  5.56           C
ATOM    226  O5'  DG B  11      12.424  -4.393  18.427  1.00 22.70           O
ATOM    227  C5'  DG B  11      12.380  -5.516  19.282  1.00 14.75           C
ATOM    228  C4'  DG B  11      11.969  -5.112  20.676  1.00 10.42           C
ATOM    229  O4'  DG B  11      12.972  -4.192  21.210  1.00 10.51           O
ATOM    230  C3'  DG B  11      10.649  -4.394  20.782  1.00  8.57           C
ATOM    231  O3'  DG B  11       9.618  -5.363  20.846  1.00  8.69           O
ATOM    232  C2'  DG B  11      10.822  -3.597  22.051  1.00  8.63           C
ATOM    233  C1'  DG B  11      12.236  -3.233  21.980  1.00  9.81           C
ATOM    234  N9   DG B  11      12.509  -1.902  21.305  1.00  8.66           N
ATOM    235  C8   DG B  11      13.175  -1.667  20.135  1.00  9.57           C
ATOM    236  N7   DG B  11      13.255  -0.407  19.824  1.00  9.04           N
ATOM    237  C5   DG B  11      12.613   0.235  20.869  1.00  7.63           C
ATOM    238  C6   DG B  11      12.388   1.612  21.119  1.00  7.05           C
ATOM    239  O6   DG B  11      12.723   2.590  20.419  1.00  7.81           O
ATOM    240  N1   DG B  11      11.715   1.819  22.317  1.00  6.27           N
ATOM    241  C2   DG B  11      11.264   0.828  23.159  1.00  6.05           C
ATOM    242  N2   DG B  11      10.611   1.219  24.248  1.00  5.85           N
ATOM    243  N3   DG B  11      11.483  -0.457  22.942  1.00  6.55           N
ATOM    244  C4   DG B  11      12.150  -0.687  21.797  1.00  6.84           C
ATOM    245  P    DC B  12       8.134  -5.009  20.350  1.00  8.13           P
ATOM    246  OP1  DC B  12       7.367  -6.252  20.459  1.00 10.02           O
ATOM    247  OP2  DC B  12       8.172  -4.307  19.052  1.00  9.79           O
ATOM    248  O5'  DC B  12       7.564  -3.912  21.389  1.00  8.18           O
ATOM    249  C5'  DC B  12       7.275  -4.296  22.719  1.00  8.00           C
ATOM    250  C4'  DC B  12       6.856  -3.057  23.487  1.00  8.01           C
ATOM    251  O4'  DC B  12       8.006  -2.146  23.615  1.00  7.35           O
ATOM    252  C3'  DC B  12       5.763  -2.208  22.890  1.00  7.04           C
ATOM    253  O3'  DC B  12       4.456  -2.800  23.100  1.00  9.82           O
ATOM    254  C2'  DC B  12       6.019  -0.916  23.630  1.00  6.50           C
ATOM    255  C1'  DC B  12       7.467  -0.808  23.608  1.00  7.35           C
ATOM    256  N1   DC B  12       8.040  -0.143  22.396  1.00  6.64           N
ATOM    257  C2   DC B  12       8.017   1.257  22.382  1.00  5.68           C
ATOM    258  O2   DC B  12       7.524   1.832  23.357  1.00  6.32           O
ATOM    259  N3   DC B  12       8.543   1.930  21.312  1.00  6.18           N
ATOM    260  C4   DC B  12       9.009   1.236  20.266  1.00  6.48           C
ATOM    261  N4   DC B  12       9.518   1.926  19.243  1.00  7.43           N
ATOM    262  C5   DC B  12       9.012  -0.198  20.248  1.00  6.83           C
ATOM    263  C6   DC B  12       8.502  -0.825  21.311  1.00  6.80           C
"""

raw_records3 = """\
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1           1
ATOM      1  N   ALA A   2      -7.656   2.923   3.155  1.00 15.02           N
ATOM      2  CA  ALA A   2      -6.522   2.038   2.831  1.00 14.10           C
ATOM      3  C   ALA A   2      -5.241   2.537   3.427  1.00 13.13           C
ATOM      4  O   ALA A   2      -4.978   3.742   3.426  1.00 11.91           O
ATOM      5  CB  ALA A   2      -6.346   1.881   1.341  1.00 15.38           C
ATOM      7  N   ALA A   3      -4.438   1.590   3.905  1.00 12.26           N
ATOM      8  CA  ALA A   3      -3.193   1.904   4.589  1.00 11.74           C
ATOM      9  C   ALA A   3      -1.955   1.332   3.895  1.00 11.10           C
ATOM     10  O   ALA A   3      -1.872   0.119   3.648  1.00 10.42           O
ATOM     11  CB  ALA A   3      -3.259   1.378   6.042  1.00 12.15           C
ATOM     13  N   ALA A   4      -1.005   2.228   3.598  1.00 10.29           N
ATOM     14  CA  ALA A   4       0.384   1.888   3.199  1.00 10.53           C
ATOM     15  C   ALA A   4       1.435   2.606   4.088  1.00 10.24           C
ATOM     16  O   ALA A   4       1.547   3.843   4.115  1.00  8.86           O
ATOM     17  CB  ALA A   4       0.656   2.148   1.711  1.00  9.80           C
END
"""

def make_initial_grm(mon_lib_srv, ener_lib, records):
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv    = mon_lib_srv,
    ener_lib       = ener_lib,
    raw_records    = records,
    force_symmetry = True)

  geometry = processed_pdb_file.geometry_restraints_manager(
    show_energies      = True,
    plain_pairs_radius = 5.0)
  xrs = processed_pdb_file.xray_structure()
  return geometry, xrs

def make_geo_pickle_unpickle(geometry, xrs, prefix):
  init_out = StringIO()
  from_file_out = StringIO()
  geometry.show_sorted(
      sites_cart=xrs.sites_cart(),
      site_labels=xrs.scatterers().extract_labels(),
      f=init_out)
  energy_original = geometry.energies_sites(sites_cart=xrs.sites_cart())
  t0 = time()
  #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
  pklstr = pickle.dumps(geometry)
  t1 = time()
  grm_from_file = pickle.loads(pklstr)
  t2 = time()
  # Fails here:
  energy_from_pickle = grm_from_file.energies_sites(sites_cart=xrs.sites_cart())
  assert approx_equal(energy_original.target, energy_from_pickle.target)
  print("Time pickling/unpickling: %.4f, %.4f" % (t1-t0, t2-t1))
  grm_from_file.show_sorted(
      sites_cart=xrs.sites_cart(),
      site_labels=xrs.scatterers().extract_labels(),
      f=from_file_out)
  # print "INITIAL"
  init_v = init_out.getvalue()
  # print init_v
  # print "="*50
  # print "From disc"
  from_file_v = from_file_out.getvalue()
  # print from_file_v
  # STOP()
  assert not show_diff(init_v, from_file_v)
  return grm_from_file

def test_simple_protein(
    mon_lib_srv, ener_lib, prefix="tst_grm_pickling_simple_protein"):
  geometry, xrs = make_initial_grm(mon_lib_srv, ener_lib, raw_records1)
  make_geo_pickle_unpickle(geometry, xrs, prefix)

def test_nucleic_acid(mon_lib_srv, ener_lib, prefix="tst_grm_pickling_na"):
  with open("%s.pdb" % prefix, "w") as f:
    f.write(raw_records2)
  from mmtbx import monomer_library
  params = monomer_library.pdb_interpretation.master_params.extract()
  params.secondary_structure.enabled=True
  processed_pdb_file = monomer_library.pdb_interpretation.run(
    args=["%s.pdb" % prefix],
    params=params,
    strict_conflict_handling=False,
    log=sys.stdout)
  geo = processed_pdb_file.geometry_restraints_manager()
  assert geo.get_n_hbond_proxies() == 6
  assert geo.get_n_hangle_proxies() == 12
  assert geo.get_n_stacking_proxies() == 2
  assert geo.get_n_parallelity_bp_proxies() == 2
  make_geo_pickle_unpickle(geo, processed_pdb_file.xray_structure(), prefix)

def test_ramachandran(mon_lib_srv, ener_lib, prefix="tst_grm_pickling_rama"):
  with open("%s.pdb" % prefix, "w") as f:
    f.write(raw_records3)
  from mmtbx import monomer_library
  params = monomer_library.pdb_interpretation.master_params.extract()
  params.peptide_link.ramachandran_restraints=True
  processed_pdb_file = monomer_library.pdb_interpretation.run(
    args=["%s.pdb" % prefix],
    params=params,
    strict_conflict_handling=False,
    log=sys.stdout)
  geo = processed_pdb_file.geometry_restraints_manager()
  assert geo.get_n_ramachandran_proxies() == 1
  make_geo_pickle_unpickle(geo, processed_pdb_file.xray_structure(), prefix)

def test_cbeta(mon_lib_srv, ener_lib, prefix="tst_grm_pickling_cbeta"):
  with open("%s.pdb" % prefix, "w") as f:
    f.write(raw_records3)
  from mmtbx import monomer_library
  params = monomer_library.pdb_interpretation.master_params.extract()
  params.c_beta_restraints=True
  processed_pdb_file = monomer_library.pdb_interpretation.run(
    args=["%s.pdb" % prefix],
    params=params,
    strict_conflict_handling=False,
    log=sys.stdout)
  geo = processed_pdb_file.geometry_restraints_manager()
  assert geo.get_n_c_beta_torsion_proxies() == 6
  make_geo_pickle_unpickle(geo, processed_pdb_file.xray_structure(), prefix)

def test_reference_coordinate(mon_lib_srv, ener_lib, prefix="tst_grm_pickling_ref_coor"):
  from mmtbx.geometry_restraints import reference
  # for some strange reason without importing this the code doesn't work...
  from cctbx import adp_restraints # import dependency

  pdb_inp = iotbx.pdb.input(source_info=None, lines=raw_records3)
  params = monomer_library.pdb_interpretation.master_params.extract()
  params.reference_coordinate_restraints.enabled=False
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    params=params,
    strict_conflict_handling=False,
    pdb_inp=pdb_inp,
    log=sys.stdout)
  geo = processed_pdb_file.geometry_restraints_manager()
  pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  sites_cart = pdb_hierarchy.atoms().extract_xyz()
  rcp = reference.add_coordinate_restraints(sites_cart=sites_cart)
  geo.adopt_reference_coordinate_restraints_in_place(rcp)
  # print "number of rcr proxies:", geo.get_n_reference_coordinate_proxies()
  make_geo_pickle_unpickle(geo, processed_pdb_file.xray_structure(), prefix)

def test_secondary_structure(mon_lib_srv, ener_lib, prefix="tst_grm_pickling_ss"):
  pdb_str = """
HELIX    2   2 ASP A   37  GLY A   48  1                                  12
CRYST1  113.068  113.068   53.292  90.00  90.00  90.00 I 41          8
ATOM    266  N   ASP A  37       6.265  61.752  14.145  1.00 35.17           N
ATOM    267  CA  ASP A  37       5.251  62.335  15.056  1.00 37.08           C
ATOM    268  C   ASP A  37       5.433  61.900  16.511  1.00 37.79           C
ATOM    269  O   ASP A  37       6.443  61.316  16.858  1.00 37.54           O
ATOM    270  CB  ASP A  37       3.827  62.120  14.521  1.00 37.53           C
ATOM    271  CG  ASP A  37       3.427  60.683  14.400  1.00 38.76           C
ATOM    272  OD1 ASP A  37       4.001  59.819  15.070  1.00 38.84           O
ATOM    273  OD2 ASP A  37       2.506  60.327  13.624  1.00 41.78           O
ATOM    274  N   ASP A  38       4.467  62.205  17.382  1.00 38.31           N
ATOM    275  CA  ASP A  38       4.609  61.829  18.786  1.00 38.69           C
ATOM    276  C   ASP A  38       4.781  60.335  18.955  1.00 37.78           C
ATOM    277  O   ASP A  38       5.598  59.886  19.760  1.00 38.31           O
ATOM    278  CB  ASP A  38       3.376  62.258  19.608  1.00 39.51           C
ATOM    279  CG  ASP A  38       3.378  63.724  19.972  1.00 42.76           C
ATOM    280  OD1 ASP A  38       4.462  64.343  20.161  1.00 48.07           O
ATOM    281  OD2 ASP A  38       2.295  64.337  20.144  1.00 47.65           O
ATOM    282  N   ALA A  39       4.003  59.561  18.209  1.00 36.68           N
ATOM    283  CA  ALA A  39       4.065  58.107  18.287  1.00 36.58           C
ATOM    284  C   ALA A  39       5.433  57.607  17.773  1.00 35.91           C
ATOM    285  O   ALA A  39       6.014  56.661  18.319  1.00 35.28           O
ATOM    286  CB  ALA A  39       2.947  57.491  17.483  1.00 36.33           C
ATOM    287  N   GLY A  40       5.948  58.257  16.745  1.00 35.33           N
ATOM    288  CA  GLY A  40       7.296  57.938  16.267  1.00 35.20           C
ATOM    289  C   GLY A  40       8.386  58.218  17.295  1.00 34.81           C
ATOM    290  O   GLY A  40       9.320  57.432  17.456  1.00 34.92           O
ATOM    291  N   ARG A  41       8.297  59.351  17.981  1.00 35.65           N
ATOM    292  CA  ARG A  41       9.300  59.698  18.970  1.00 35.75           C
ATOM    293  C   ARG A  41       9.257  58.681  20.093  1.00 37.10           C
ATOM    294  O   ARG A  41      10.295  58.291  20.642  1.00 37.65           O
ATOM    295  CB  ARG A  41       9.090  61.118  19.494  1.00 36.15           C
ATOM    296  CG  ARG A  41       9.575  62.196  18.563  1.00 35.51           C
ATOM    297  CD  ARG A  41       9.383  63.592  19.134  1.00 38.98           C
ATOM    298  NE  ARG A  41       7.999  64.012  18.913  1.00 40.46           N
ATOM    299  CZ  ARG A  41       7.537  64.446  17.753  1.00 41.44           C
ATOM    300  NH1 ARG A  41       8.326  64.534  16.682  1.00 42.62           N
ATOM    301  NH2 ARG A  41       6.261  64.776  17.649  1.00 43.05           N
ATOM    302  N   ALA A  42       8.053  58.238  20.441  1.00 38.18           N
ATOM    303  CA  ALA A  42       7.878  57.270  21.524  1.00 38.42           C
ATOM    304  C   ALA A  42       8.398  55.909  21.116  1.00 38.32           C
ATOM    305  O   ALA A  42       8.952  55.181  21.927  1.00 37.15           O
ATOM    306  CB  ALA A  42       6.387  57.158  21.948  1.00 38.91           C
ATOM    307  N   THR A  43       8.209  55.567  19.842  1.00 37.57           N
ATOM    308  CA  THR A  43       8.756  54.324  19.328  1.00 37.21           C
ATOM    309  C   THR A  43      10.284  54.321  19.472  1.00 36.33           C
ATOM    310  O   THR A  43      10.842  53.315  19.824  1.00 36.44           O
ATOM    311  CB  THR A  43       8.316  54.130  17.873  1.00 37.54           C
ATOM    312  OG1 THR A  43       6.890  53.948  17.829  1.00 38.41           O
ATOM    313  CG2 THR A  43       8.897  52.837  17.280  1.00 36.05           C
ATOM    314  N   LEU A  44      10.948  55.436  19.192  1.00 36.74           N
ATOM    315  CA  LEU A  44      12.410  55.504  19.283  1.00 36.66           C
ATOM    316  C   LEU A  44      12.877  55.316  20.729  1.00 37.24           C
ATOM    317  O   LEU A  44      13.840  54.613  20.978  1.00 36.26           O
ATOM    318  CB  LEU A  44      12.957  56.819  18.725  1.00 36.22           C
ATOM    319  CG  LEU A  44      12.786  57.061  17.209  1.00 34.97           C
ATOM    320  CD1 LEU A  44      13.386  58.400  16.795  1.00 33.89           C
ATOM    321  CD2 LEU A  44      13.399  55.928  16.404  1.00 33.15           C
ATOM    322  N   ARG A  45      12.147  55.914  21.675  1.00 38.31           N
ATOM    323  CA  ARG A  45      12.485  55.801  23.095  1.00 39.49           C
ATOM    324  C   ARG A  45      12.296  54.381  23.589  1.00 39.97           C
ATOM    325  O   ARG A  45      13.113  53.864  24.338  1.00 40.63           O
ATOM    326  CB  ARG A  45      11.614  56.757  23.935  1.00 39.94           C
ATOM    327  N   ARG A  46      11.186  53.775  23.179  1.00 41.00           N
ATOM    328  CA  ARG A  46      10.849  52.397  23.503  1.00 41.33           C
ATOM    329  C   ARG A  46      11.912  51.412  23.025  1.00 40.34           C
ATOM    330  O   ARG A  46      12.278  50.485  23.731  1.00 39.81           O
ATOM    331  CB  ARG A  46       9.524  52.063  22.835  1.00 41.72           C
ATOM    332  CG  ARG A  46       8.773  50.911  23.395  1.00 46.36           C
ATOM    333  CD  ARG A  46       7.352  50.836  22.851  1.00 51.59           C
ATOM    334  NE  ARG A  46       7.345  50.162  21.548  1.00 57.79           N
ATOM    335  CZ  ARG A  46       6.851  50.659  20.399  1.00 61.01           C
ATOM    336  NH1 ARG A  46       6.282  51.872  20.344  1.00 62.67           N
ATOM    337  NH2 ARG A  46       6.918  49.916  19.290  1.00 61.73           N
ATOM    338  N   LEU A  47      12.402  51.620  21.809  1.00 39.47           N
ATOM    339  CA  LEU A  47      13.439  50.765  21.223  1.00 38.25           C
ATOM    340  C   LEU A  47      14.826  51.006  21.800  1.00 37.39           C
ATOM    341  O   LEU A  47      15.742  50.247  21.530  1.00 38.19           O
ATOM    342  CB  LEU A  47      13.502  51.010  19.712  1.00 38.57           C
ATOM    343  CG  LEU A  47      12.264  50.556  18.951  1.00 38.58           C
ATOM    344  CD1 LEU A  47      12.346  51.046  17.517  1.00 38.92           C
ATOM    345  CD2 LEU A  47      12.101  49.038  19.050  1.00 38.51           C
ATOM    346  N   GLY A  48      14.997  52.083  22.557  1.00 36.96           N
ATOM    347  CA  GLY A  48      16.262  52.383  23.191  1.00 35.72           C
ATOM    348  C   GLY A  48      17.323  52.969  22.286  1.00 34.43           C
ATOM    349  O   GLY A  48      18.512  52.912  22.607  1.00 34.93           O
  """

  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  params = monomer_library.pdb_interpretation.master_params.extract()
  params.secondary_structure.enabled=True
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    params=params,
    strict_conflict_handling=False,
    pdb_inp=pdb_inp,
    log=sys.stdout)
  geo = processed_pdb_file.geometry_restraints_manager()
  assert geo.get_n_hbond_proxies() == 8
  make_geo_pickle_unpickle(geo, processed_pdb_file.xray_structure(), prefix)

def test_secondary_structure_2(mon_lib_srv, ener_lib, prefix="tst_grm_pickling_ss2"):
  from iotbx.pdb.tst_secondary_structure import pdb_1ywf_sample_strings
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_1ywf_sample_strings)
  params = monomer_library.pdb_interpretation.master_params.extract()
  params.secondary_structure.enabled=True
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    params=params,
    strict_conflict_handling=False,
    pdb_inp=pdb_inp,
    log=sys.stdout)
  geo = processed_pdb_file.geometry_restraints_manager()
  assert geo.get_n_hbond_proxies() == 103, geo.get_n_hbond_proxies()
  make_geo_pickle_unpickle(geo, processed_pdb_file.xray_structure(), prefix)

def test_across_symmetry(mon_lib_srv, ener_lib, prefix="tst_across_symmetry"):
  raw_records1 = """\
CRYST1   60.800   60.800   97.000  90.00  90.00 120.00 P 32 2 1      6
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.016447  0.009496  0.000000        0.00000
SCALE2      0.000000  0.018992  0.000000        0.00000
SCALE3      0.000000  0.000000  0.010309        0.00000
ATOM   1050  N   LYS A 135      31.992  14.930  -7.233  1.00  9.47           N
ATOM   1051  CA  LYS A 135      31.388  16.216  -7.637  1.00 12.89           C
ATOM   1052  C   LYS A 135      30.807  16.840  -6.406  1.00  6.47           C
ATOM   1053  O   LYS A 135      29.583  16.869  -6.191  1.00 15.74           O
ATOM   1054  CB  LYS A 135      30.263  16.059  -8.655  1.00 13.51           C
ATOM   1055  CG  LYS A 135      30.742  15.277  -9.843  1.00 16.23           C
ATOM   1056  CD  LYS A 135      29.612  15.131 -10.835  1.00 28.55           C
ATOM   1057  CE  LYS A 135      30.173  14.812 -12.216  1.00 34.52           C
ATOM   1058  NZ  LYS A 135      29.396  13.756 -12.899  1.00 46.18           N
HETATM 1406  O   HOH A 282      32.366  19.942  24.727  1.00 38.09           O
END
"""
  geometry, xrs = make_initial_grm(mon_lib_srv, ener_lib, raw_records1)
  make_geo_pickle_unpickle(geometry, xrs, prefix)

def test_reference_model(mon_lib_srv, ener_lib, prefix="tst_reference_model"):
  from mmtbx.geometry_restraints.torsion_restraints.tst_reference_model import \
      model_raw_records, reference_raw_records
  from mmtbx.geometry_restraints.torsion_restraints.reference_model import \
    reference_model
  import mmtbx.model
  # mstream = StringIO()
  from libtbx.utils import multi_out
  mstream = multi_out()
  mstream.register("stdout", sys.stdout)
  mstream_file_name = "polder.log"
  mstreamfile = open(mstream_file_name, "w")
  mstream.register("logfile", mstreamfile)

  work_params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  work_params.reference_model.enabled = True
  work_params.reference_model.use_starting_model_as_reference = True
  work_params.reference_model.fix_outliers = False
  pdb_inp = iotbx.pdb.input(lines=model_raw_records.split('\n'), source_info=None)
  model = mmtbx.model.manager(model_input = pdb_inp)
  model.process(pdb_interpretation_params=work_params,
    make_restraints=True)
  reference_hierarchy_list = []
  tmp_hierarchy = iotbx.pdb.input(
    source_info=None,
    lines=reference_raw_records.split('\n')).construct_hierarchy()
  reference_hierarchy_list.append(tmp_hierarchy)
  rm = reference_model(
         model = model,
         reference_hierarchy_list=reference_hierarchy_list,
         params=work_params.reference_model,
         log=null_out())  # XXX changed from mstream which cannot be pickled
  assert rm.get_n_proxies() == 5, "Got %d, expected 5" % rm.get_n_proxies()
  geometry, xrs = make_initial_grm(mon_lib_srv, ener_lib, model_raw_records)
  geometry.adopt_reference_dihedral_manager(rm)
  make_geo_pickle_unpickle(geometry, xrs, prefix)
  mstream.close()

def exercise_all(args):
  mon_lib_srv = None
  ener_lib = None
  try:
    mon_lib_srv = monomer_library.server.server()
    ener_lib = monomer_library.server.ener_lib()
  except Exception:
    print("Can not initialize monomer_library, skipping test.")
    return 0
  import libtbx.load_env
  if libtbx.env.find_in_repositories(relative_path="chem_data") is None:
    print("Skipping exercise(): chem_data directory not available")
    return

  test_simple_protein(mon_lib_srv, ener_lib)
  test_nucleic_acid(mon_lib_srv, ener_lib)
  # test_ramachandran(mon_lib_srv, ener_lib)
  # test_cbeta(mon_lib_srv, ener_lib)
  # test_reference_coordinate(mon_lib_srv, ener_lib)
  # test_secondary_structure(mon_lib_srv, ener_lib)
  # test_secondary_structure_2(mon_lib_srv, ener_lib)
  # test_across_symmetry(mon_lib_srv, ener_lib)
  test_reference_model(mon_lib_srv, ener_lib)

if (__name__ == "__main__"):
  exercise_all(sys.argv[1:])
