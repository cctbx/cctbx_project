from __future__ import absolute_import, division, print_function
import iotbx.pdb
import mmtbx.f_model

pdb_str="""\
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
"""

def run():
  # get xray_structure object from input PDB file
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  xray_structure = pdb_inp.xray_structure_simple()
  # make up "F-obs" and R-free-flags
  f_obs = abs(xray_structure.structure_factors(d_min=1.5).f_calc())
  r_free_flags = f_obs.generate_r_free_flags()
  # instantiate fmodel object
  fmodel = mmtbx.f_model.manager(
    f_obs          = f_obs,
    r_free_flags   = r_free_flags,
    xray_structure = xray_structure)
  # now we can caclualte R-factors and show some statistics
  print("r_work:", fmodel.r_work())
  print("r_free:", fmodel.r_free())
  fmodel.show()
  # we can get get 3d array of map values
  fft_map = fmodel.electron_density_map().fft_map(map_type = "mFobs-DFmodel")
  fft_map.apply_sigma_scaling()
  map_data = fft_map.real_map_unpadded(in_place=False)
  # output map in X-plor formatted file
  fft_map.as_xplor_map(
    file_name="mfo-dfm.xplor",
    title_lines=["mFobs-DFmodel"],
    gridding_first=(0,0,0),
    gridding_last=fft_map.n_real())
  # we can compute another map
  map_coeffs = fmodel.electron_density_map(
    ).map_coefficients(map_type = "2mFobs-DFmodel")
  # output it as MTZ file
  mtz_dataset = map_coeffs.as_mtz_dataset(column_root_label="2mFoDFc")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = "amap.mtz")
  # and get actual map values
  fft_map = map_coeffs.fft_map(resolution_factor=0.25)
  fft_map.apply_volume_scaling()
  map_data = fft_map.real_map_unpadded()

if (__name__ == "__main__"):
  run()
  print("OK")
