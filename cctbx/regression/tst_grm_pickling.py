from __future__ import division
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
import sys
import pickle
import StringIO
from libtbx.test_utils import show_diff
from libtbx.utils import null_out
from mmtbx.monomer_library import pdb_interpretation
from time import time


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
  init_out = StringIO.StringIO()
  from_file_out = StringIO.StringIO()
  geometry.show_sorted(
      sites_cart=xrs.sites_cart(),
      site_labels=xrs.scatterers().extract_labels(),
      f=init_out)

  t0 = time()
  pklfile = open("%s.pkl" % prefix, 'wb')
  pickle.dump(geometry, pklfile)
  pklfile.close()
  t1 = time()
  pklfile = open("%s.pkl" % prefix, 'rb')
  grm_from_file = pickle.load(pklfile)
  pklfile.close()
  t2 = time()
  print "Time pickling/unpickling: %.4f, %.4f" % (t1-t0, t2-t1)
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
  assert not show_diff(init_v, from_file_v)

def test_simple_protein(
    mon_lib_srv, ener_lib, prefix="tst_grm_pickling_simple_protein"):
  geometry, xrs = make_initial_grm(mon_lib_srv, ener_lib, raw_records1)
  make_geo_pickle_unpickle(geometry, xrs, prefix)

def test_nucleic_acid(mon_lib_srv, ener_lib, prefix="tst_grm_pickling_na"):
  open("%s.pdb" % prefix, "w").\
      write(raw_records2)
  from mmtbx import monomer_library
  params = monomer_library.pdb_interpretation.master_params.extract()
  params.secondary_structure.enabled=True
  processed_pdb_file = pdb_interpretation.run(
    args=["%s.pdb" % prefix],
    params=params,
    strict_conflict_handling=False,
    log=null_out())
  geo = processed_pdb_file.geometry_restraints_manager()
  assert geo.get_n_hbond_proxies() == 6
  assert geo.get_n_hangle_proxies() == 12
  assert geo.get_n_stacking_proxies() == 2
  assert geo.get_n_parallelity_bp_proxies() == 2
  make_geo_pickle_unpickle(geo, processed_pdb_file.xray_structure(), prefix)

def exercise_all(args):
  mon_lib_srv = None
  ener_lib = None
  try:
    mon_lib_srv = monomer_library.server.server()
    ener_lib = monomer_library.server.ener_lib()
  except Exception:
    print "Can not initialize monomer_library, skipping test."
    return 0
  test_simple_protein(mon_lib_srv, ener_lib)
  # test_nucleic_acid(mon_lib_srv, ener_lib)

if (__name__ == "__main__"):
  exercise_all(sys.argv[1:])
