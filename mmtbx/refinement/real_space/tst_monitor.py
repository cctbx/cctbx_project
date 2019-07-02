from __future__ import absolute_import, division, print_function
import mmtbx.monomer_library.pdb_interpretation
import iotbx.mtz
from cctbx.array_family import flex
import time
from mmtbx import monomer_library
import iotbx.pdb
from libtbx import group_args
import mmtbx.restraints
import mmtbx.refinement.real_space

pdb_str_answer="""\
CRYST1   16.960   19.455   19.841  90.00  90.00  90.00 P 1
ATOM      1  N   ASP A  18      14.727  11.381  11.191  1.00 21.41           N
ATOM      2  CA  ASP A  18      14.022  12.655  11.131  1.00 34.70           C
ATOM      3  C   ASP A  18      13.805  13.085   9.684  1.00 35.50           C
ATOM      4  O   ASP A  18      14.760  13.408   8.976  1.00 38.48           O
ATOM      5  CB  ASP A  18      14.798  13.731  11.893  1.00 38.15           C
ATOM      6  CG  ASP A  18      14.068  15.059  11.937  1.00 27.39           C
ATOM      7  OD1 ASP A  18      13.282  15.275  12.883  1.00 35.78           O
ATOM      8  OD2 ASP A  18      14.282  15.889  11.028  1.00 30.67           O
ATOM      9  N   ASN A  19      12.544  13.082   9.257  1.00 22.94           N
ATOM     10  CA  ASN A  19      12.169  13.463   7.896  1.00 22.82           C
ATOM     11  C   ASN A  19      12.886  12.623   6.841  1.00 20.25           C
ATOM     12  O   ASN A  19      13.859  13.069   6.233  1.00 33.24           O
ATOM     13  CB  ASN A  19      12.419  14.959   7.663  1.00 28.83           C
ATOM     14  CG  ASN A  19      11.659  15.505   6.465  1.00 26.44           C
ATOM     15  OD1 ASN A  19      11.356  14.780   5.518  1.00 22.64           O
ATOM     16  ND2 ASN A  19      11.350  16.797   6.504  1.00 28.91           N
ATOM     17  N   TYR A  20      12.397  11.404   6.631  1.00 23.81           N
ATOM     18  CA  TYR A  20      12.992  10.492   5.662  1.00 35.69           C
ATOM     19  C   TYR A  20      11.909   9.888   4.773  1.00 28.94           C
ATOM     20  O   TYR A  20      11.281  10.589   3.980  1.00 25.26           O
ATOM     21  CB  TYR A  20      13.765   9.382   6.379  1.00 33.37           C
ATOM     22  CG  TYR A  20      13.470   9.266   7.861  1.00 31.48           C
ATOM     23  CD1 TYR A  20      12.437   8.463   8.327  1.00 36.91           C
ATOM     24  CD2 TYR A  20      14.228   9.963   8.794  1.00 21.36           C
ATOM     25  CE1 TYR A  20      12.169   8.355   9.678  1.00 36.56           C
ATOM     26  CE2 TYR A  20      13.967   9.862  10.147  1.00 29.56           C
ATOM     27  CZ  TYR A  20      12.936   9.057  10.583  1.00 35.08           C
ATOM     28  OH  TYR A  20      12.671   8.953  11.930  1.00 38.73           O
ATOM     29  N   ARG A  21      11.695   8.584   4.914  1.00 38.95           N
ATOM     30  CA  ARG A  21      10.695   7.881   4.120  1.00 38.77           C
ATOM     31  C   ARG A  21       9.825   6.998   5.010  1.00 27.69           C
ATOM     32  O   ARG A  21       8.635   6.819   4.750  1.00 22.82           O
ATOM     33  CB  ARG A  21      11.373   7.048   3.027  1.00 20.00           C
ATOM     34  CG  ARG A  21      10.419   6.410   2.025  1.00 20.00           C
ATOM     35  CD  ARG A  21      10.098   4.968   2.391  1.00 20.00           C
ATOM     36  NE  ARG A  21       9.189   4.346   1.434  1.00 20.00           N
ATOM     37  CZ  ARG A  21       8.730   3.103   1.542  1.00 20.00           C
ATOM     38  NH1 ARG A  21       9.094   2.347   2.569  1.00 20.00           N
ATOM     39  NH2 ARG A  21       7.906   2.617   0.625  1.00 20.00           N
TER
ATOM     41  N   ASP B   3      10.807  10.584  13.275  1.00 21.41           N
ATOM     42  CA  ASP B   3      10.009  11.773  13.002  1.00 34.70           C
ATOM     43  C   ASP B   3       9.681  11.887  11.519  1.00 35.50           C
ATOM     44  O   ASP B   3      10.577  11.850  10.675  1.00 38.48           O
ATOM     45  CB  ASP B   3      10.741  13.031  13.473  1.00 38.15           C
ATOM     46  CG  ASP B   3       9.956  14.298  13.202  1.00 27.39           C
ATOM     47  OD1 ASP B   3       9.139  14.688  14.063  1.00 35.78           O
ATOM     48  OD2 ASP B   3      10.155  14.906  12.130  1.00 30.67           O
ATOM     49  N   ASN B   4       8.392  12.025  11.215  1.00 22.94           N
ATOM     50  CA  ASN B   4       7.914  12.139   9.839  1.00 22.82           C
ATOM     51  C   ASN B   4       8.362  10.974   8.963  1.00 20.25           C
ATOM     52  O   ASN B   4       8.974  11.171   7.914  1.00 33.24           O
ATOM     53  CB  ASN B   4       8.347  13.472   9.221  1.00 28.83           C
ATOM     54  CG  ASN B   4       7.800  14.668   9.976  1.00 26.44           C
ATOM     55  OD1 ASN B   4       6.697  14.622  10.520  1.00 22.64           O
ATOM     56  ND2 ASN B   4       8.572  15.747  10.013  1.00 28.91           N
ATOM     57  N   ARG B   5       8.056   9.759   9.405  1.00 23.81           N
ATOM     58  CA  ARG B   5       8.452   8.557   8.682  1.00 35.69           C
ATOM     59  C   ARG B   5       7.240   7.898   8.028  1.00 28.94           C
ATOM     60  O   ARG B   5       7.119   6.673   8.018  1.00 25.26           O
ATOM     61  CB  ARG B   5       9.150   7.569   9.619  1.00 20.00           C
ATOM     62  CG  ARG B   5      10.498   8.049  10.133  1.00 20.00           C
ATOM     63  CD  ARG B   5      11.214   7.048  11.026  1.00 20.00           C
ATOM     64  NE  ARG B   5      12.468   7.581  11.550  1.00 20.00           N
ATOM     65  CZ  ARG B   5      13.227   6.959  12.441  1.00 20.00           C
ATOM     66  NH1 ARG B   5      12.864   5.773  12.913  1.00 20.00           N
ATOM     67  NH2 ARG B   5      14.354   7.519  12.862  1.00 20.00           N
ATOM     68  N   ARG B   6       6.348   8.727   7.490  1.00 38.95           N
ATOM     69  CA  ARG B   6       5.152   8.255   6.792  1.00 38.77           C
ATOM     70  C   ARG B   6       4.304   7.308   7.640  1.00 27.69           C
ATOM     71  O   ARG B   6       3.395   7.741   8.349  1.00 22.82           O
ATOM     72  CB  ARG B   6       5.531   7.599   5.460  1.00 20.00           C
ATOM     73  CG  ARG B   6       4.458   6.697   4.873  1.00 20.00           C
ATOM     74  CD  ARG B   6       4.877   5.236   4.948  1.00 20.00           C
ATOM     75  NE  ARG B   6       3.881   4.343   4.364  1.00 20.00           N
ATOM     76  CZ  ARG B   6       4.001   3.020   4.317  1.00 20.00           C
ATOM     77  NH1 ARG B   6       5.077   2.431   4.821  1.00 20.00           N
ATOM     78  NH2 ARG B   6       3.045   2.284   3.766  1.00 20.00           N
TER
END
"""

pdb_str_poor="""\n
CRYST1   16.960   19.455   19.841  90.00  90.00  90.00 P 1
ATOM      1  N   ASP A  18      14.727  11.381  11.191  1.00 21.41           N
ATOM      2  CA  ASP A  18      14.022  12.655  11.131  1.00 34.70           C
ATOM      3  C   ASP A  18      13.805  13.085   9.684  1.00 35.50           C
ATOM      4  O   ASP A  18      14.760  13.408   8.976  1.00 38.48           O
ATOM      5  CB  ASP A  18      14.798  13.731  11.893  1.00 38.15           C
ATOM      6  CG  ASP A  18      14.068  15.059  11.937  1.00 27.39           C
ATOM      7  OD1 ASP A  18      13.282  15.275  12.883  1.00 35.78           O
ATOM      8  OD2 ASP A  18      14.282  15.889  11.028  1.00 30.67           O
ATOM      9  N   ASN A  19      12.544  13.082   9.257  1.00 22.94           N
ATOM     10  CA  ASN A  19      12.169  13.463   7.896  1.00 22.82           C
ATOM     11  C   ASN A  19      12.886  12.623   6.841  1.00 20.25           C
ATOM     12  O   ASN A  19      13.859  13.069   6.233  1.00 33.24           O
ATOM     13  CB  ASN A  19      12.419  14.959   7.663  1.00 28.83           C
ATOM     14  CG  ASN A  19      11.659  15.505   6.465  1.00 26.44           C
ATOM     15  OD1 ASN A  19      11.356  14.780   5.518  1.00 22.64           O
ATOM     16  ND2 ASN A  19      11.350  16.797   6.504  1.00 28.91           N
ATOM     17  N   TYR A  20      12.397  11.404   6.631  1.00 23.81           N
ATOM     18  CA  TYR A  20      12.992  10.492   5.662  1.00 35.69           C
ATOM     19  C   TYR A  20      11.909   9.888   4.773  1.00 28.94           C
ATOM     20  O   TYR A  20      11.281  10.589   3.980  1.00 25.26           O
ATOM     21  CB  TYR A  20      13.765   9.382   6.379  1.00 33.37           C
ATOM     22  CG  TYR A  20      13.470   9.266   7.861  0.10  5.48           C
ATOM     23  CD1 TYR A  20      12.437   8.463   8.327  1.50  5.91           C
ATOM     24  CD2 TYR A  20      14.228   9.963   8.794  0.50 95.36           C
ATOM     25  CE1 TYR A  20      12.169   8.355   9.678  0.50  5.56           C
ATOM     26  CE2 TYR A  20      13.967   9.862  10.147  0.50  5.56           C
ATOM     27  CZ  TYR A  20      12.936   9.057  10.583  0.30  5.08           C
ATOM     28  OH  TYR A  20      12.671   8.953  11.930  1.00 38.73           O
ATOM     29  N   ARG A  21      11.695   8.584   4.914  1.00 38.95           N
ATOM     30  CA  ARG A  21      10.695   7.881   4.120  1.00 38.77           C
ATOM     31  C   ARG A  21       9.825   6.998   5.010  1.00 27.69           C
ATOM     32  O   ARG A  21       8.635   6.819   4.750  1.00 22.82           O
ATOM     33  CB  ARG A  21      11.373   7.048   3.027  1.00 20.00           C
ATOM     34  CG  ARG A  21      10.419   6.410   2.025  1.00 20.00           C
ATOM     35  CD  ARG A  21      10.098   4.968   2.391  1.00 20.00           C
ATOM     36  NE  ARG A  21       9.189   4.346   1.434  1.00 20.00           N
ATOM     37  CZ  ARG A  21       8.730   3.103   1.542  1.00 20.00           C
ATOM     38  NH1 ARG A  21       9.094   2.347   2.569  1.00 20.00           N
ATOM     39  NH2 ARG A  21       7.906   2.617   0.625  1.00 20.00           N
TER
ATOM     41  N   ASP B   3      10.807  10.584  13.275  1.00 21.41           N
ATOM     42  CA  ASP B   3      10.009  11.773  13.002  1.00 34.70           C
ATOM     43  C   ASP B   3       9.681  11.887  11.519  1.00 35.50           C
ATOM     44  O   ASP B   3      10.577  11.850  10.675  1.00 38.48           O
ATOM     45  CB  ASP B   3      10.741  13.031  13.473  1.00 38.15           C
ATOM     46  CG  ASP B   3       9.956  14.298  13.202  1.00 27.39           C
ATOM     47  OD1 ASP B   3       9.139  14.688  14.063  1.00 35.78           O
ATOM     48  OD2 ASP B   3      10.155  14.906  12.130  1.00 30.67           O
ATOM     49  N   ASN B   4       8.392  12.025  11.215  1.00 22.94           N
ATOM     50  CA  ASN B   4       7.914  12.139   9.839  1.00 22.82           C
ATOM     51  C   ASN B   4       8.362  10.974   8.963  1.00 20.25           C
ATOM     52  O   ASN B   4       8.974  11.171   7.914  1.00 33.24           O
ATOM     53  CB  ASN B   4       8.347  13.472   9.221  1.00 28.83           C
ATOM     54  CG  ASN B   4       7.800  14.668   9.976  1.00 26.44           C
ATOM     55  OD1 ASN B   4       6.697  14.622  10.520  1.00 22.64           O
ATOM     56  ND2 ASN B   4       8.572  15.747  10.013  1.00 28.91           N
ATOM     57  N   ARG B   5       8.056   9.759   9.405  1.00 23.81           N
ATOM     58  CA  ARG B   5       8.452   8.557   8.682  1.00 35.69           C
ATOM     59  C   ARG B   5       7.240   7.898   8.028  1.00 28.94           C
ATOM     60  O   ARG B   5       7.119   6.673   8.018  1.00 25.26           O
ATOM     61  CB  ARG B   5       9.150   7.569   9.619  1.00 20.00           C
ATOM     62  CG  ARG B   5      10.498   8.049  10.133  1.00 20.00           C
ATOM     63  CD  ARG B   5      11.214   7.048  11.026  1.00 20.00           C
ATOM     64  NE  ARG B   5      12.468   7.581  11.550  1.00 20.00           N
ATOM     65  CZ  ARG B   5      13.227   6.959  12.441  1.00 20.00           C
ATOM     66  NH1 ARG B   5      12.864   5.773  12.913  1.00 20.00           N
ATOM     67  NH2 ARG B   5      14.354   7.519  12.862  1.00 20.00           N
ATOM     68  N   ARG B   6       6.348   8.727   7.490  1.00 38.95           N
ATOM     69  CA  ARG B   6       5.152   8.255   6.792  1.00 38.77           C
ATOM     70  C   ARG B   6       4.304   7.308   7.640  1.00 27.69           C
ATOM     71  O   ARG B   6       3.395   7.741   8.349  1.00 22.82           O
ATOM     72  CB  ARG B   6       5.531   7.599   5.460  1.00 20.00           C
ATOM     73  CG  ARG B   6       4.458   6.697   4.873  1.00 20.00           C
ATOM     74  CD  ARG B   6       4.877   5.236   4.948  1.00 20.00           C
ATOM     75  NE  ARG B   6       3.881   4.343   4.364  1.00 20.00           N
ATOM     76  CZ  ARG B   6       4.001   3.020   4.317  1.00 20.00           C
ATOM     77  NH1 ARG B   6       5.077   2.431   4.821  1.00 20.00           N
ATOM     78  NH2 ARG B   6       3.045   2.284   3.766  1.00 20.00           N
TER
END
"""

def exercise(d_min=1.5, resolution_factor = 0.25):
  # answer
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str_answer)
  pdb_inp.write_pdb_file(file_name = "answer.pdb")
  xrs_answer = pdb_inp.xray_structure_simple()
  f_calc = xrs_answer.structure_factors(d_min = d_min).f_calc()
  fft_map = f_calc.fft_map(resolution_factor=resolution_factor)
  fft_map.apply_sigma_scaling()
  target_map = fft_map.real_map_unpadded()
  mtz_dataset = f_calc.as_mtz_dataset(column_root_label = "FCmap")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = "answer.mtz")
  # poor
  mon_lib_srv = monomer_library.server.server()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv              = mon_lib_srv,
    ener_lib                 = monomer_library.server.ener_lib(),
    raw_records              = flex.std_string(pdb_str_poor.splitlines()),
    strict_conflict_handling = True,
    force_symmetry           = True,
    log                      = None)
  pdb_hierarchy_poor = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  xrs_poor = processed_pdb_file.xray_structure()
  sites_cart_poor = xrs_poor.sites_cart()
  pdb_hierarchy_poor.write_pdb_file(file_name = "poor.pdb")
  ####
  target_map_object = group_args(
    map_data         = target_map,
    miller_array     = f_calc,
    crystal_gridding = fft_map,
    d_min            = d_min)
  grm = mmtbx.restraints.manager(
    geometry=processed_pdb_file.geometry_restraints_manager(show_energies=False),
    normalization = True)
  sm = mmtbx.refinement.real_space.structure_monitor(
    pdb_hierarchy               = pdb_hierarchy_poor,
    xray_structure              = xrs_poor,
    target_map_object           = target_map_object,
    geometry_restraints_manager = grm.geometry)
  sm.show(prefix="start")
  sm.show_residues(map_cc_all=2) # XXX assert printed output table:

if(__name__ == "__main__"):
  t0 = time.time()
  exercise()
  print("Time: %6.4f"%(time.time()-t0))
