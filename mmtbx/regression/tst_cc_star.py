
# XXX This is a very minimal test with no file dependencies - a more thorough
# set of tests is located in phenix_regression.

from __future__ import absolute_import, division, print_function
from mmtbx.regression import model_1yjp
from mmtbx.command_line import cc_star
import iotbx.pdb
from cctbx import sgtbx
from cctbx import uctbx
from scitbx.array_family import flex
from libtbx.test_utils import Exception_expected
from libtbx.utils import Sorry, null_out
import random

def exercise_space_group_handling():
  flex.set_random_seed(123456)
  random.seed(123456)
  base = "tst_cc_star_space_group"
  pdb_in = iotbx.pdb.input(source_info=None, lines=model_1yjp)
  xrs = pdb_in.xray_structure_simple()
  xrs.set_inelastic_form_factors(
    photon=1.54,
    table="sasaki")
  fc = abs(xrs.structure_factors(d_min=1.5).f_calc()).average_bijvoet_mates()
  fc.set_observation_type_xray_amplitude()
  flags = fc.generate_r_free_flags()
  mtz = fc.as_mtz_dataset(column_root_label="F")
  mtz.add_miller_array(flags, column_root_label="FreeR_flag")
  mtz.mtz_object().write(base + ".mtz")
  xrs_p1 = xrs.expand_to_p1()
  xrs_p1.shake_sites_in_place(rms_difference=0.1)
  fc_p1 = xrs_p1.structure_factors(d_min=1.4).f_calc()
  fc_p1_extra = fc_p1.randomize_amplitude_and_phase(amplitude_error=1.0,
    phase_error_deg=0,
    random_seed=123456)
  fc_p1 = abs(fc_p1.concatenate(other=fc_p1_extra)).sort(
    by_value="packed_indices")
  fc_p1.set_observation_type_xray_amplitude()
  sg_p2 = sgtbx.space_group_info("P2")
  ic = fc_p1.f_as_f_sq().customized_copy(
    space_group_info=sg_p2,
    sigmas=flex.double(fc_p1.size(), 10.0))
  ic.export_as_scalepack_unmerged(file_name=base + ".sca")
  with open(base + ".pdb", "w") as f:
    f.write(model_1yjp)
  args = [
    base + ".mtz",
    base + ".pdb",
    "unmerged_data=%s.sca" % base,
  ]
  cc_star.run(args=args, out=null_out())
  # now with .sca in P1 (raises Sorry)
  ic2 = fc_p1.f_as_f_sq().customized_copy(
    sigmas=flex.double(fc_p1.size(), 10.0))
  ic2.export_as_scalepack_unmerged(file_name=base + "_p1.sca")
  args = [
    base + ".mtz",
    base + ".pdb",
    "unmerged_data=%s_p1.sca" % base,
  ]
  try :
    cc_star.run(args=args, out=null_out())
  except Sorry as s :
    assert (str(s) == "Incompatible space groups in merged and unmerged data:P 1 21 1 versus P 1"), s
  else :
    raise Exception_expected
  # now with CIF (complete symmetry)
  f = open(base + ".cif", "w")
  ic.as_cif_simple(array_type="meas",
    out=f)
  f.close()
  args = [
    base + ".mtz",
    base + ".pdb",
    "unmerged_data=%s.cif" % base,
  ]
  cc_star.run(args=args, out=null_out())
  # bad unit cell
  uc2 = uctbx.unit_cell((23,6.5,23.5,90,108,90))
  ic3 = ic.customized_copy(unit_cell=uc2)
  f = open(base + "_new_uc.cif", "w")
  ic3.as_cif_simple(array_type="meas",
    out=f)
  f.close()
  args = [
    base + ".mtz",
    base + ".pdb",
    "unmerged_data=%s_new_uc.cif" % base,
  ]
  try :
    cc_star.run(args=args, out=null_out())
  except Sorry as s :
    assert ("Incompatible symmetry definitions:" in str(s)), s
  else :
    raise Exception_expected

if (__name__ == "__main__"):
  exercise_space_group_handling()
  print("OK")
