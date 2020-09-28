from __future__ import absolute_import, division, print_function
import time
import iotbx.pdb
import mmtbx.utils
from scitbx.array_family import flex
from cctbx import maptbx
from libtbx import adopt_init_args
from mmtbx import monomer_library
import mmtbx.refinement.real_space.explode_and_refine

pdb_str_answer = """\
CRYST1   29.475   46.191   27.490  90.00  90.00  90.00 P 21 21 21
HELIX    1   1 ALA E    1  ALA E   16  1                                  16
ATOM      1  N   ALA E   1      14.622  35.477  13.274  1.00 40.00           N
ATOM      2  CA  ALA E   1      14.323  35.055  14.635  1.00 40.00           C
ATOM      3  C   ALA E   1      13.442  33.808  14.656  1.00 40.00           C
ATOM      4  O   ALA E   1      13.848  32.796  15.227  1.00 40.00           O
ATOM      5  CB  ALA E   1      13.675  36.191  15.414  1.00 40.00           C
ATOM      6  N   ALA E   2      12.256  33.870  14.012  1.00 40.00           N
ATOM      7  CA  ALA E   2      11.288  32.768  13.916  1.00 40.00           C
ATOM      8  C   ALA E   2      11.847  31.573  13.125  1.00 40.00           C
ATOM      9  O   ALA E   2      11.553  30.422  13.462  1.00 40.00           O
ATOM     10  CB  ALA E   2      10.000  33.267  13.276  1.00 40.00           C
ATOM     11  N   ALA E   3      12.664  31.856  12.084  1.00 40.00           N
ATOM     12  CA  ALA E   3      13.313  30.854  11.231  1.00 40.00           C
ATOM     13  C   ALA E   3      14.391  30.082  11.994  1.00 40.00           C
ATOM     14  O   ALA E   3      14.616  28.906  11.687  1.00 40.00           O
ATOM     15  CB  ALA E   3      13.915  31.517  10.000  1.00 40.00           C
ATOM     16  N   ALA E   4      15.063  30.744  12.976  1.00 40.00           N
ATOM     17  CA  ALA E   4      16.089  30.116  13.822  1.00 40.00           C
ATOM     18  C   ALA E   4      15.441  29.043  14.703  1.00 40.00           C
ATOM     19  O   ALA E   4      15.888  27.897  14.671  1.00 40.00           O
ATOM     20  CB  ALA E   4      16.785  31.160  14.683  1.00 40.00           C
ATOM     21  N   ALA E   5      14.338  29.398  15.423  1.00 40.00           N
ATOM     22  CA  ALA E   5      13.573  28.488  16.298  1.00 40.00           C
ATOM     23  C   ALA E   5      13.042  27.288  15.513  1.00 40.00           C
ATOM     24  O   ALA E   5      12.981  26.179  16.048  1.00 40.00           O
ATOM     25  CB  ALA E   5      12.433  29.233  16.975  1.00 40.00           C
ATOM     26  N   ALA E   6      12.715  27.518  14.221  1.00 40.00           N
ATOM     27  CA  ALA E   6      12.258  26.514  13.260  1.00 40.00           C
ATOM     28  C   ALA E   6      13.440  25.620  12.865  1.00 40.00           C
ATOM     29  O   ALA E   6      13.361  24.405  13.060  1.00 40.00           O
ATOM     30  CB  ALA E   6      11.664  27.193  12.032  1.00 40.00           C
ATOM     31  N   ALA E   7      14.556  26.234  12.374  1.00 40.00           N
ATOM     32  CA  ALA E   7      15.788  25.543  11.962  1.00 40.00           C
ATOM     33  C   ALA E   7      16.411  24.694  13.075  1.00 40.00           C
ATOM     34  O   ALA E   7      16.954  23.625  12.782  1.00 40.00           O
ATOM     35  CB  ALA E   7      16.806  26.540  11.433  1.00 40.00           C
ATOM     36  N   ALA E   8      16.327  25.168  14.343  1.00 40.00           N
ATOM     37  CA  ALA E   8      16.846  24.459  15.515  1.00 40.00           C
ATOM     38  C   ALA E   8      16.009  23.197  15.768  1.00 40.00           C
ATOM     39  O   ALA E   8      16.568  22.097  15.858  1.00 40.00           O
ATOM     40  CB  ALA E   8      16.833  25.368  16.742  1.00 40.00           C
ATOM     41  N   ALA E   9      14.664  23.357  15.803  1.00 40.00           N
ATOM     42  CA  ALA E   9      13.716  22.261  15.994  1.00 40.00           C
ATOM     43  C   ALA E   9      13.754  21.279  14.822  1.00 40.00           C
ATOM     44  O   ALA E   9      13.654  20.075  15.051  1.00 40.00           O
ATOM     45  CB  ALA E   9      12.314  22.805  16.193  1.00 40.00           C
ATOM     46  N   ALA E  10      13.942  21.791  13.578  1.00 40.00           N
ATOM     47  CA  ALA E  10      14.063  20.985  12.357  1.00 40.00           C
ATOM     48  C   ALA E  10      15.322  20.116  12.421  1.00 40.00           C
ATOM     49  O   ALA E  10      15.235  18.910  12.169  1.00 40.00           O
ATOM     50  CB  ALA E  10      14.105  21.881  11.130  1.00 40.00           C
ATOM     51  N   ALA E  11      16.480  20.722  12.803  1.00 40.00           N
ATOM     52  CA  ALA E  11      17.754  20.018  12.967  1.00 40.00           C
ATOM     53  C   ALA E  11      17.631  18.976  14.089  1.00 40.00           C
ATOM     54  O   ALA E  11      18.074  17.837  13.908  1.00 40.00           O
ATOM     55  CB  ALA E  11      18.870  21.004  13.279  1.00 40.00           C
ATOM     56  N   ALA E  12      16.972  19.356  15.217  1.00 40.00           N
ATOM     57  CA  ALA E  12      16.721  18.483  16.371  1.00 40.00           C
ATOM     58  C   ALA E  12      15.851  17.277  15.972  1.00 40.00           C
ATOM     59  O   ALA E  12      16.201  16.140  16.310  1.00 40.00           O
ATOM     60  CB  ALA E  12      16.052  19.269  17.490  1.00 40.00           C
ATOM     61  N   ALA E  13      14.748  17.535  15.215  1.00 40.00           N
ATOM     62  CA  ALA E  13      13.821  16.516  14.708  1.00 40.00           C
ATOM     63  C   ALA E  13      14.530  15.563  13.752  1.00 40.00           C
ATOM     64  O   ALA E  13      14.266  14.360  13.793  1.00 40.00           O
ATOM     65  CB  ALA E  13      12.636  17.173  14.008  1.00 40.00           C
ATOM     66  N   ALA E  14      15.437  16.112  12.903  1.00 40.00           N
ATOM     67  CA  ALA E  14      16.246  15.367  11.939  1.00 40.00           C
ATOM     68  C   ALA E  14      17.277  14.511  12.676  1.00 40.00           C
ATOM     69  O   ALA E  14      17.512  13.363  12.286  1.00 40.00           O
ATOM     70  CB  ALA E  14      16.942  16.329  10.986  1.00 40.00           C
ATOM     71  N   ALA E  15      17.863  15.069  13.767  1.00 40.00           N
ATOM     72  CA  ALA E  15      18.847  14.407  14.629  1.00 40.00           C
ATOM     73  C   ALA E  15      18.213  13.257  15.427  1.00 40.00           C
ATOM     74  O   ALA E  15      18.922  12.322  15.816  1.00 40.00           O
ATOM     75  CB  ALA E  15      19.475  15.418  15.573  1.00 40.00           C
ATOM     76  N   ALA E  16      16.879  13.325  15.652  1.00 40.00           N
ATOM     77  CA  ALA E  16      16.099  12.317  16.370  1.00 40.00           C
ATOM     78  C   ALA E  16      15.822  11.106  15.488  1.00 40.00           C
ATOM     79  O   ALA E  16      15.658  10.000  15.999  1.00 40.00           O
ATOM     80  CB  ALA E  16      14.790  12.918  16.852  1.00 40.00           C
TER
END
"""

pdb_str_poor = """\
CRYST1   29.475   46.191   27.490  90.00  90.00  90.00 P 21 21 21
ATOM      1  N   ALA E   1      17.517  34.071   6.538  1.00 40.00           N
ATOM      2  CA  ALA E   1      17.992  35.348   7.004  1.00 40.00           C
ATOM      3  C   ALA E   1      16.863  35.986   7.792  1.00 40.00           C
ATOM      4  O   ALA E   1      16.713  37.242   7.909  1.00 40.00           O
ATOM      5  CB  ALA E   1      18.477  36.116   5.812  1.00 40.00           C
ATOM      6  N   ALA E   2      16.080  35.113   8.348  1.00 40.00           N
ATOM      7  CA  ALA E   2      15.254  35.581   9.495  1.00 40.00           C
ATOM      8  C   ALA E   2      15.628  35.262  10.863  1.00 40.00           C
ATOM      9  O   ALA E   2      16.066  36.143  11.629  1.00 40.00           O
ATOM     10  CB  ALA E   2      13.821  35.144   9.297  1.00 40.00           C
ATOM     11  N   ALA E   3      15.670  34.035  11.221  1.00 40.00           N
ATOM     12  CA  ALA E   3      15.746  33.569  12.729  1.00 40.00           C
ATOM     13  C   ALA E   3      16.380  32.197  13.043  1.00 40.00           C
ATOM     14  O   ALA E   3      17.557  32.017  13.009  1.00 40.00           O
ATOM     15  CB  ALA E   3      14.390  33.574  13.297  1.00 40.00           C
ATOM     16  N   ALA E   4      15.437  31.334  13.516  1.00 40.00           N
ATOM     17  CA  ALA E   4      15.810  29.953  13.900  1.00 40.00           C
ATOM     18  C   ALA E   4      14.876  28.948  13.219  1.00 40.00           C
ATOM     19  O   ALA E   4      15.088  28.561  12.135  1.00 40.00           O
ATOM     20  CB  ALA E   4      15.807  29.747  15.562  1.00 40.00           C
ATOM     21  N   ALA E   5      13.827  28.582  13.930  1.00 40.00           N
ATOM     22  CA  ALA E   5      12.803  27.643  13.582  1.00 40.00           C
ATOM     23  C   ALA E   5      11.925  27.972  12.295  1.00 40.00           C
ATOM     24  O   ALA E   5      11.311  29.016  12.217  1.00 40.00           O
ATOM     25  CB  ALA E   5      11.884  27.496  14.800  1.00 40.00           C
ATOM     26  N   ALA E   6      11.824  27.002  11.326  1.00 40.00           N
ATOM     27  CA  ALA E   6      11.579  27.271   9.948  1.00 40.00           C
ATOM     28  C   ALA E   6      10.688  26.172   9.315  1.00 40.00           C
ATOM     29  O   ALA E   6       9.810  26.435   8.479  1.00 40.00           O
ATOM     30  CB  ALA E   6      12.892  27.434   9.163  1.00 40.00           C
ATOM     31  N   ALA E   7      10.985  24.936   9.568  1.00 40.00           N
ATOM     32  CA  ALA E   7      10.223  23.915   8.805  1.00 40.00           C
ATOM     33  C   ALA E   7       9.182  23.205   9.603  1.00 40.00           C
ATOM     34  O   ALA E   7       7.945  23.417   9.434  1.00 40.00           O
ATOM     35  CB  ALA E   7      11.183  22.804   8.200  1.00 40.00           C
ATOM     36  N   ALA E   8       9.745  22.435  10.542  1.00 40.00           N
ATOM     37  CA  ALA E   8       9.210  21.425  11.476  1.00 40.00           C
ATOM     38  C   ALA E   8      10.256  20.419  11.733  1.00 40.00           C
ATOM     39  O   ALA E   8       9.992  19.490  12.460  1.00 40.00           O
ATOM     40  CB  ALA E   8       7.965  20.817  10.947  1.00 40.00           C
ATOM     41  N   ALA E   9      11.466  20.737  11.268  1.00 40.00           N
ATOM     42  CA  ALA E   9      12.607  19.889  11.232  1.00 40.00           C
ATOM     43  C   ALA E   9      12.924  19.328  12.523  1.00 40.00           C
ATOM     44  O   ALA E   9      12.790  20.050  13.539  1.00 40.00           O
ATOM     45  CB  ALA E   9      13.739  20.662  10.656  1.00 40.00           C
ATOM     46  N   ALA E  10      13.495  18.150  12.482  1.00 40.00           N
ATOM     47  CA  ALA E  10      14.007  17.595  13.654  1.00 40.00           C
ATOM     48  C   ALA E  10      15.157  16.724  13.724  1.00 40.00           C
ATOM     49  O   ALA E  10      15.159  16.000  14.620  1.00 40.00           O
ATOM     50  CB  ALA E  10      12.870  16.878  14.366  1.00 40.00           C
ATOM     51  N   ALA E  11      16.092  16.782  12.766  1.00 40.00           N
ATOM     52  CA  ALA E  11      17.076  15.791  12.430  1.00 40.00           C
ATOM     53  C   ALA E  11      16.462  14.506  11.790  1.00 40.00           C
ATOM     54  O   ALA E  11      15.459  14.602  11.077  1.00 40.00           O
ATOM     55  CB  ALA E  11      17.943  15.352  13.622  1.00 40.00           C
ATOM     56  N   ALA E  12      17.191  13.458  11.773  1.00 40.00           N
ATOM     57  CA  ALA E  12      16.985  12.246  11.074  1.00 40.00           C
ATOM     58  C   ALA E  12      17.956  11.169  11.712  1.00 40.00           C
ATOM     59  O   ALA E  12      18.839  10.657  11.173  1.00 40.00           O
ATOM     60  CB  ALA E  12      17.278  12.440   9.655  1.00 40.00           C
ATOM     61  N   ALA E  13      17.693  10.781  12.936  1.00 40.00           N
ATOM     62  CA  ALA E  13      18.559   9.862  13.648  1.00 40.00           C
ATOM     63  C   ALA E  13      17.678   8.677  14.105  1.00 40.00           C
ATOM     64  O   ALA E  13      16.810   8.236  13.332  1.00 40.00           O
ATOM     65  CB  ALA E  13      19.221  10.551  14.826  1.00 40.00           C
ATOM     66  N   ALA E  14      17.647   8.453  15.439  1.00 40.00           N
ATOM     67  CA  ALA E  14      16.666   7.551  16.066  1.00 40.00           C
ATOM     68  C   ALA E  14      15.298   8.147  16.400  1.00 40.00           C
ATOM     69  O   ALA E  14      14.633   7.785  17.376  1.00 40.00           O
ATOM     70  CB  ALA E  14      17.224   7.034  17.308  1.00 40.00           C
ATOM     71  N   ALA E  15      14.855   9.000  15.485  1.00 40.00           N
ATOM     72  CA  ALA E  15      13.517   9.710  15.635  1.00 40.00           C
ATOM     73  C   ALA E  15      12.304   8.858  15.586  1.00 40.00           C
ATOM     74  O   ALA E  15      11.625   8.812  14.564  1.00 40.00           O
ATOM     75  CB  ALA E  15      13.302  10.842  14.526  1.00 40.00           C
ATOM     76  N   ALA E  16      12.014   8.140  16.649  1.00 40.00           N
ATOM     77  CA  ALA E  16      11.279   6.863  16.560  1.00 40.00           C
ATOM     78  C   ALA E  16      10.003   6.986  17.409  1.00 40.00           C
ATOM     79  O   ALA E  16       9.883   8.020  18.086  1.00 40.00           O
ATOM     80  CB  ALA E  16      12.188   5.653  16.999  1.00 40.00           C
TER
"""

def ccp4_map(crystal_symmetry, file_name, map_data):
  from iotbx import mrcfile
  mrcfile.write_ccp4_map(
      file_name=file_name,
      unit_cell=crystal_symmetry.unit_cell(),
      space_group=crystal_symmetry.space_group(),
      #gridding_first=(0,0,0),# This causes a bug (map gets shifted)
      #gridding_last=n_real,  # This causes a bug (map gets shifted)
      map_data=map_data,
      labels=flex.std_string([""]))

class scorer(object):
  def __init__(self, pdb_hierarchy, unit_cell, map_data):
    adopt_init_args(self, locals())
    self.sites_cart = self.pdb_hierarchy.atoms().extract_xyz()
    self.target = maptbx.real_space_target_simple(
      unit_cell   = self.unit_cell,
      density_map = self.map_data,
      sites_cart  = self.sites_cart)

  def update(self, sites_cart):
    target = maptbx.real_space_target_simple(
      unit_cell   = self.unit_cell,
      density_map = self.map_data,
      sites_cart  = sites_cart)
    if(target > self.target):
      self.target = target
      self.sites_cart = sites_cart
    print(self.target, target) # XXX for debugging

def run(prefix="tst_00"):
  # Good answer model
  pdb_file_name_answer = "%s_answer.pdb"%prefix
  pdb_inp = iotbx.pdb.input(source_info=None, lines = pdb_str_answer)
  pdb_inp.write_pdb_file(file_name="%s_answer.pdb"%prefix)
  ph_answer = pdb_inp.construct_hierarchy()
  ph_answer.atoms().reset_i_seq()
  xrs_answer = pdb_inp.xray_structure_simple()
  # Poor model that we want to refine so it matches the answer
  pdb_file_name_poor = "%s_poor.pdb"%prefix
  pdb_inp = iotbx.pdb.input(source_info=None, lines = pdb_str_poor)
  pdb_inp.write_pdb_file(file_name="%s_poor.pdb"%prefix)
  ph_poor = pdb_inp.construct_hierarchy()
  ph_poor.atoms().reset_i_seq()
  xrs_poor = pdb_inp.xray_structure_simple()
  # Initialize states accumulator
  states = mmtbx.utils.states(pdb_hierarchy=ph_answer)
  states.add(sites_cart = xrs_poor.sites_cart())
  # Compute target map
  fc = xrs_answer.structure_factors(d_min=3.5).f_calc()
  fft_map = fc.fft_map(resolution_factor = 0.25)
  fft_map.apply_sigma_scaling()
  target_map_data = fft_map.real_map_unpadded()
  ccp4_map(crystal_symmetry=fc.crystal_symmetry(), file_name="map.ccp4",
    map_data=target_map_data)
  # Output map coefficients
  mtz_dataset = fc.as_mtz_dataset(column_root_label="FC")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = "map.mtz")
  # Build geometry restraints
  params = monomer_library.pdb_interpretation.master_params.extract()
  #params.nonbonded_weight=200
  #params.peptide_link.ramachandran_restraints=True
  #params.peptide_link.rama_potential="emsley"
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv              = monomer_library.server.server(),
    ener_lib                 = monomer_library.server.ener_lib(),
    raw_records              = pdb_str_poor,
    params                   = params,
    strict_conflict_handling = True,
    force_symmetry           = True,
    log                      = None)

  geometry = processed_pdb_file.geometry_restraints_manager(
    show_energies                = False,
    plain_pairs_radius           = 5,
    assume_hydrogens_all_missing = True)
  restraints_manager = mmtbx.restraints.manager(
    geometry      = geometry,
    normalization = True)
  # Do real-space refinement
  t0=time.time()
  ear = mmtbx.refinement.real_space.explode_and_refine.run(
    xray_structure          = xrs_poor,
    pdb_hierarchy           = ph_poor,
    map_data                = target_map_data,
    restraints_manager      = restraints_manager,
    states                  = states,
    nproc=1)
  print("Time: %6.4f"%(time.time()-t0))
  ear.pdb_hierarchy.write_pdb_file(file_name="%s_refined.pdb"%prefix)
  states.write(file_name="%s_refined_all_states.pdb"%prefix)

if (__name__ == "__main__"):
  run()
  print("OK")
