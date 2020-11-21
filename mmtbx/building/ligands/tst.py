from __future__ import absolute_import, division, print_function
import time
import iotbx.pdb
import mmtbx.utils
from scitbx.array_family import flex
from cctbx import maptbx
from libtbx import adopt_init_args
from mmtbx import monomer_library
import mmtbx.refinement.real_space.explode_and_refine
import mmtbx.building.ligands
from libtbx.utils import null_out

pdb_str_answer = """\
CRYST1   25.172   35.024   28.879  90.00  90.00  90.00 P 1
HETATM    1  C1' ATP A   1      11.188  20.531  13.981  1.00 20.00      A    C
HETATM    2  C2  ATP A   1      11.407  24.221  16.551  1.00 20.00      A    C
HETATM    3  C2' ATP A   1      11.602  20.999  12.594  1.00 20.00      A    C
HETATM    4  C3' ATP A   1      11.507  19.778  11.698  1.00 20.00      A    C
HETATM    5  C4  ATP A   1      12.085  22.140  15.684  1.00 20.00      A    C
HETATM    6  C4' ATP A   1      11.047  18.640  12.596  1.00 20.00      A    C
HETATM    7  C5  ATP A   1      13.263  22.174  16.562  1.00 20.00      A    C
HETATM    8  C5' ATP A   1      11.885  17.382  12.400  1.00 20.00      A    C
HETATM    9  C6  ATP A   1      13.415  23.355  17.444  1.00 20.00      A    C
HETATM   10  C8  ATP A   1      13.279  20.341  15.375  1.00 20.00      A    C
HETATM   11  N1  ATP A   1      12.465  24.314  17.380  1.00 20.00      A    N
HETATM   12  N3  ATP A   1      11.213  23.175  15.727  1.00 20.00      A    N
HETATM   13  N6  ATP A   1      14.471  23.470  18.285  1.00 20.00      A    N
HETATM   14  N7  ATP A   1      13.952  21.044  16.323  1.00 20.00      A    N
HETATM   15  N9  ATP A   1      12.165  20.998  14.997  1.00 20.00      A    N
HETATM   16  O1A ATP A   1      11.768  14.190  12.067  1.00 20.00      A    O
HETATM   17  O1B ATP A   1      13.387  12.899  15.592  1.00 20.00      A    O
HETATM   18  O1G ATP A   1      12.562  11.112  17.675  1.00 20.00      A    O
HETATM   19  O2' ATP A   1      10.742  22.048  12.137  1.00 20.00      A    O
HETATM   20  O2A ATP A   1      13.631  15.127  13.643  1.00 20.00      A    O1-
HETATM   21  O2B ATP A   1      11.677  11.753  13.980  1.00 20.00      A    O1-
HETATM   22  O2G ATP A   1      10.080  11.175  18.205  1.00 20.00      A    O1-
HETATM   23  O3' ATP A   1      10.575  19.988  10.632  1.00 20.00      A    O
HETATM   24  O3A ATP A   1      11.472  14.201  14.578  1.00 20.00      A    O
HETATM   25  O3B ATP A   1      11.008  12.447  16.317  1.00 20.00      A    O
HETATM   26  O3G ATP A   1      10.901  10.000  16.110  1.00 20.00      A    O1-
HETATM   27  O4' ATP A   1      11.142  19.102  13.947  1.00 20.00      A    O
HETATM   28  O5' ATP A   1      11.455  16.370  13.308  1.00 20.00      A    O
HETATM   29  PA  ATP A   1      12.169  14.928  13.321  1.00 20.00      A    P
HETATM   30  PB  ATP A   1      11.978  12.754  15.070  1.00 20.00      A    P
HETATM   31  PG  ATP A   1      11.149  11.073  17.144  1.00 20.00      A    P
HETATM   32  H2  ATP A   1      10.680  25.024  16.548  1.00 20.00      A    H
HETATM   33  H1' ATP A   1      10.190  20.932  14.208  1.00 20.00      A    H
HETATM   34  H2' ATP A   1      12.645  21.345  12.626  1.00 20.00      A    H
HETATM   35  H3' ATP A   1      12.503  19.544  11.296  1.00 20.00      A    H
HETATM   36  H4' ATP A   1      10.000  18.405  12.357  1.00 20.00      A    H
HETATM   37  H8  ATP A   1      13.594  19.383  14.978  1.00 20.00      A    H
HETATM   38 H5'1 ATP A   1      12.938  17.610  12.575  1.00 20.00      A    H
HETATM   39 H5'2 ATP A   1      11.779  17.026  11.373  1.00 20.00      A    H
HETATM   40 HN61 ATP A   1      14.558  24.282  18.879  1.00 20.00      A    H
HETATM   41 HN62 ATP A   1      15.172  22.743  18.319  1.00 20.00      A    H
HETATM   42 HO2' ATP A   1      11.084  22.410  11.308  1.00 20.00      A    H
HETATM   43 HO3' ATP A   1      10.942  20.622  10.000  1.00 20.00      A    H
END
"""

pdb_str_poor = """\
CRYST1   25.172   35.024   28.879  90.00  90.00  90.00 P 1
HETATM    1  C1' ATP A   1      25.651  36.544 -15.828  1.00 20.00      A    C
HETATM    2  C2  ATP A   1      30.090  36.377 -15.092  1.00 20.00      A    C
HETATM    3  C2' ATP A   1      25.036  35.778 -14.666  1.00 20.00      A    C
HETATM    4  C3' ATP A   1      23.538  35.791 -14.905  1.00 20.00      A    C
HETATM    5  C4  ATP A   1      27.939  37.331 -15.163  1.00 20.00      A    C
HETATM    6  C4' ATP A   1      23.332  36.571 -16.196  1.00 20.00      A    C
HETATM    7  C5  ATP A   1      28.496  38.591 -14.652  1.00 20.00      A    C
HETATM    8  C5' ATP A   1      22.231  37.616 -16.062  1.00 20.00      A    C
HETATM    9  C6  ATP A   1      29.952  38.625 -14.382  1.00 20.00      A    C
HETATM   10  C8  ATP A   1      26.345  38.807 -14.964  1.00 20.00      A    C
HETATM   11  N1  ATP A   1      30.663  37.501 -14.623  1.00 20.00      A    N
HETATM   12  N3  ATP A   1      28.774  36.282 -15.358  1.00 20.00      A    N
HETATM   13  N6  ATP A   1      30.551  39.745 -13.909  1.00 20.00      A    N
HETATM   14  N7  ATP A   1      27.468  39.452 -14.553  1.00 20.00      A    N
HETATM   15  N9  ATP A   1      26.630  37.541 -15.327  1.00 20.00      A    N
HETATM   16  O1A ATP A   1      19.729  38.563 -14.285  1.00 20.00      A    O
HETATM   17  O1B ATP A   1      16.569  37.739 -16.745  1.00 20.00      A    O
HETATM   18  O1G ATP A   1      14.136  36.242 -16.499  1.00 20.00      A    O
HETATM   19  O2' ATP A   1      25.543  34.440 -14.617  1.00 20.00      A    O
HETATM   20  O2A ATP A   1      19.402  38.620 -16.879  1.00 20.00      A    O1-
HETATM   21  O2B ATP A   1      16.712  37.867 -14.138  1.00 20.00      A    O1-
HETATM   22  O2G ATP A   1      14.335  34.112 -15.132  1.00 20.00      A    O1-
HETATM   23  O3' ATP A   1      23.020  34.463 -15.030  1.00 20.00      A    O
HETATM   24  O3A ATP A   1      18.506  36.706 -15.490  1.00 20.00      A    O
HETATM   25  O3B ATP A   1      16.217  35.684 -15.314  1.00 20.00      A    O
HETATM   26  O3G ATP A   1      14.276  36.367 -13.967  1.00 20.00      A    O1-
HETATM   27  O4' ATP A   1      24.579  37.192 -16.518  1.00 20.00      A    O
HETATM   28  O5' ATP A   1      20.987  36.973 -15.796  1.00 20.00      A    O
HETATM   29  PA  ATP A   1      19.643  37.840 -15.608  1.00 20.00      A    P
HETATM   30  PB  ATP A   1      16.950  37.113 -15.424  1.00 20.00      A    P
HETATM   31  PG  ATP A   1      14.612  35.594 -15.221  1.00 20.00      A    P
HETATM   32  H2  ATP A   1      30.716  35.510 -15.264  1.00 20.00      A    H
HETATM   33  H1' ATP A   1      26.151  35.830 -16.498  1.00 20.00      A    H
HETATM   34  H2' ATP A   1      25.258  36.305 -13.727  1.00 20.00      A    H
HETATM   35  H3' ATP A   1      23.043  36.321 -14.079  1.00 20.00      A    H
HETATM   36  H4' ATP A   1      23.052  35.865 -16.991  1.00 20.00      A    H
HETATM   37  H8  ATP A   1      25.357  39.249 -14.996  1.00 20.00      A    H
HETATM   38 H5'1 ATP A   1      22.158  38.194 -16.986  1.00 20.00      A    H
HETATM   39 H5'2 ATP A   1      22.473  38.304 -15.249  1.00 20.00      A    H
HETATM   40 HN61 ATP A   1      31.545  39.750 -13.730  1.00 20.00      A    H
HETATM   41 HN62 ATP A   1      30.003  40.576 -13.735  1.00 20.00      A    H
HETATM   42 HO2' ATP A   1      25.232  34.004 -13.813  1.00 20.00      A    H
HETATM   43 HO3' ATP A   1      23.040  34.024 -14.168  1.00 20.00      A    H
END
"""

def ccp4_map(crystal_symmetry, file_name, map_data):
  from iotbx import mrcfile
  mrcfile.write_ccp4_map(
    file_name=file_name,
    unit_cell=crystal_symmetry.unit_cell(),
    space_group=crystal_symmetry.space_group(),
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
  model_answer = mmtbx.model.manager(model_input = pdb_inp)
  pdb_inp.write_pdb_file(file_name=pdb_file_name_answer)
  # Poor model that we want to refine so it matches the answer
  pdb_file_name_poor = "%s_poor.pdb"%prefix
  pdb_inp = iotbx.pdb.input(source_info=None, lines = pdb_str_poor)
  model_poor = mmtbx.model.manager(model_input = pdb_inp, build_grm = True,
    log = null_out())
  pdb_inp.write_pdb_file(file_name=pdb_file_name_poor)
  # Initialize states accumulator
  states = mmtbx.utils.states(
    pdb_hierarchy  = model_poor.get_hierarchy(),
    xray_structure = model_poor.get_xray_structure())
  states.add(sites_cart = model_poor.get_sites_cart())
  # Compute target map
  fc = model_answer.get_xray_structure().structure_factors(d_min=3.5).f_calc()
  fft_map = fc.fft_map(resolution_factor = 0.25)
  fft_map.apply_sigma_scaling()
  map_data = fft_map.real_map_unpadded()
  ccp4_map(crystal_symmetry=fc.crystal_symmetry(), file_name="map.ccp4",
    map_data=map_data)
  # Output map coefficients
  mtz_dataset = fc.as_mtz_dataset(column_root_label="FC")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = "map.mtz")
  # Do fitting
  result = mmtbx.building.ligands.run(model = model_poor, map_data = map_data)


#  # Build geometry restraints
#  params = monomer_library.pdb_interpretation.master_params.extract()
#  #params.nonbonded_weight=200
#  #params.peptide_link.ramachandran_restraints=True
#  #params.peptide_link.rama_potential="emsley"
#  processed_pdb_file = monomer_library.pdb_interpretation.process(
#    mon_lib_srv              = monomer_library.server.server(),
#    ener_lib                 = monomer_library.server.ener_lib(),
#    raw_records              = pdb_str_poor,
#    params                   = params,
#    strict_conflict_handling = True,
#    force_symmetry           = True,
#    log                      = None)
#
#  geometry = processed_pdb_file.geometry_restraints_manager(
#    show_energies                = False,
#    plain_pairs_radius           = 5,
#    assume_hydrogens_all_missing = True)
#  restraints_manager = mmtbx.restraints.manager(
#    geometry      = geometry,
#    normalization = True)
#  # Do real-space refinement
#  t0=time.time()
#  ear = mmtbx.refinement.real_space.explode_and_refine.run(
#    xray_structure          = xrs_poor,
#    pdb_hierarchy           = ph_poor,
#    map_data                = target_map_data,
#    restraints_manager      = restraints_manager,
#    states                  = states,
#    nproc=1)
#  print "Time: %6.4f"%(time.time()-t0)
#  ear.pdb_hierarchy.write_pdb_file(file_name="%s_refined.pdb"%prefix)
#  states.write(file_name="%s_refined_all_states.pdb"%prefix)

if (__name__ == "__main__"):
  run()
  print("OK")
