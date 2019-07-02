
from __future__ import absolute_import, division, print_function
from libtbx import easy_run
import libtbx.load_env
import os

def exercise_simple():
  mtz_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/reflection_files/ha_patterson.mtz",
    test=os.path.isfile)
  if (mtz_file is None):
    print("phenix_regression not available, skipping")
    return
  result = easy_run.fully_buffered("cctbx.patterson_map \"%s\"" % mtz_file
    ).raise_if_errors()
  assert (result.return_code == 0)

def exercise_isomorphous():
  import iotbx.pdb.hierarchy
  pdb_in = iotbx.pdb.hierarchy.input(pdb_string="""\
ATOM     47  N   TYR A   7       8.292   1.817   6.147  1.00 14.70           N
ATOM     48  CA  TYR A   7       9.159   2.144   7.299  1.00 15.18           C
ATOM     49  C   TYR A   7      10.603   2.331   6.885  1.00 15.91           C
ATOM     50  O   TYR A   7      11.041   1.811   5.855  1.00 15.76           O
ATOM     51  CB  TYR A   7       9.061   1.065   8.369  1.00 15.35           C
ATOM     52  CG  TYR A   7       7.665   0.929   8.902  1.00 14.45           C
ATOM     53  CD1 TYR A   7       6.771   0.021   8.327  1.00 15.68           C
ATOM     54  CD2 TYR A   7       7.210   1.756   9.920  1.00 14.80           C
ATOM     55  CE1 TYR A   7       5.480  -0.094   8.796  1.00 13.46           C
ATOM     56  CE2 TYR A   7       5.904   1.649  10.416  1.00 14.33           C
ATOM     57  CZ  TYR A   7       5.047   0.729   9.831  1.00 15.09           C
ATOM     58  OH  TYR A   7       3.766   0.589  10.291  1.00 14.39           O
ATOM     59  OXT TYR A   7      11.358   2.999   7.612  1.00 17.49           O
ATOM    100 CL   CL  B   1      12.000   5.000  10.000  1.00 10.00          CL
""")
  xrs_1 = pdb_in.input.xray_structure_simple()
  sel = pdb_in.hierarchy.atom_selection_cache().selection("not resname CL")
  assert (sel.count(False) == 1)
  pdb_in_2 = pdb_in.hierarchy.select(sel).as_pdb_input(
    crystal_symmetry=xrs_1)
  xrs_2 = pdb_in_2.xray_structure_simple()
  fc1 = abs(xrs_1.structure_factors(d_min=1.5).f_calc())
  fc2 = abs(xrs_2.structure_factors(d_min=1.5).f_calc())
  ic1 = fc1.f_as_f_sq().set_observation_type_xray_intensity()
  #ic1 = ic1.customized_copy(data=ic1.data()*5)
  ic1.as_mtz_dataset(
    column_root_label="I").mtz_object().write("tmp_patterson_1.mtz")
  fc2.as_mtz_dataset(
    column_root_label="F").mtz_object().write("tmp_patterson_2.mtz")
  args = ["tmp_patterson_1.mtz", "tmp_patterson_2.mtz"]
  result = easy_run.fully_buffered("cctbx.isomorphous_difference_patterson %s" %
    " ".join(args)).raise_if_errors()
  assert (result.return_code == 0)

if (__name__ == "__main__"):
  exercise_simple()
  exercise_isomorphous()
  print("OK")
