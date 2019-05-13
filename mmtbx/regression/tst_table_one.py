from __future__ import division
from __future__ import print_function
from mmtbx.regression import model_1yjp
from mmtbx.command_line import table_one
import iotbx.pdb.hierarchy
from cctbx import sgtbx
from scitbx.array_family import flex
from libtbx.utils import null_out
import random

def exercise():
  flex.set_random_seed(123456)
  random.seed(123456)
  base = "tst_table_one"
  pdb_in = iotbx.pdb.hierarchy.input(pdb_string=model_1yjp)
  xrs = pdb_in.xray_structure_simple()
  xrs.set_inelastic_form_factors(
    photon=1.54,
    table="sasaki")
  fc = abs(xrs.structure_factors(d_min=1.5).f_calc()).average_bijvoet_mates()
  fc.set_observation_type_xray_amplitude()
  flags = fc.generate_r_free_flags()
  mtz = fc.as_mtz_dataset(column_root_label="F",
    wavelength=1.54)
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
  open(base + ".pdb", "w").write(model_1yjp)
  args = [
    base + ".mtz",
    base + ".pdb",
    "unmerged_data=%s.sca" % base,
    "prefix=tst_table_one_1",
  ]
# BROKEN
#  table_one.run(args=args, out=null_out(),
#    use_current_directory_if_not_specified=True)
#  # now with unmerged data in SHELX format
#  f = open(base + ".hkl", "w")
#  ic.export_as_shelx_hklf(file_object=f)
#  f.close()
#  args = [
#    base + ".mtz",
#    base + ".pdb",
#    "unmerged_data=%s.hkl=hklf4" % base,
#    "prefix=tst_table_one_2",
#  ]
#  table_one.run(args=args, out=null_out(),
#    use_current_directory_if_not_specified=True)
#  # now with phil file
#  f = open("tst_table_one_3.eff", "w")
#  f.write("""\
#table_one {
#  structure {
#    name = %(base)s
#    pdb_file = %(base)s.pdb
#    mtz_file = %(base)s.mtz
#    unmerged_data = %(base)s.hkl=hklf4
#  }
#  output {
#    directory = os.getcwd()
#    base_name = %(base)s_3
#  }
#}""" % {"base" : base })
#  args = [ "tst_table_one_3.eff" ]
#  table_one.run(args=args, out=null_out(),
#    use_current_directory_if_not_specified=True)

if (__name__ == "__main__"):
  exercise()
  print("OK")
