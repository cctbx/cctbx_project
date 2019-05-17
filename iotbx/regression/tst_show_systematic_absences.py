
from __future__ import absolute_import, division, print_function
from iotbx.command_line import show_systematic_absences
from cctbx.development import random_structure
from cctbx import sgtbx
from scitbx.array_family import flex
from libtbx.utils import Sorry
from libtbx.test_utils import Exception_expected
from six.moves import cStringIO as StringIO
import random

def exercise():
  flex.set_random_seed(12345)
  random.seed(12345)
  xrs = random_structure.xray_structure(
    sgtbx.space_group_info("P21212"),
    elements=["const"]*100)
  f_calc = xrs.structure_factors(d_min=2.5).f_calc()
  i_calc = abs(f_calc).f_as_f_sq().set_observation_type_xray_intensity()
  i_calc = i_calc.customized_copy(
    space_group_info=sgtbx.space_group_info("P222"),
    sigmas=flex.double(i_calc.size(), 1.0))
  complete_set = i_calc.complete_set()
  lone_set = complete_set.lone_set(other=i_calc)
  i_abs = lone_set.array(data=flex.double(lone_set.size(), 5),
    sigmas=flex.double(lone_set.size(), 10))
  i_calc = i_calc.concatenate(other=i_abs).set_observation_type_xray_intensity()
  i_calc.export_as_scalepack_unmerged(file_name="tst_sys_absent.sca")
  args = ["tst_sys_absent.sca"]
  out = StringIO()
  try :
    show_systematic_absences.run(args=args, out=out)
  except Sorry :
    pass
  else :
    raise Exception_expected
  args.append(",".join([ str(x) for x in xrs.unit_cell().parameters() ]))
  show_systematic_absences.run(args=args, out=out)
  assert (out.getvalue().count("  (   0,    3,    0): i/sigi =    0.5") == 4)
  i_calc = i_calc.customized_copy(sigmas=None)
  i_calc.as_mtz_dataset(column_root_label="I").mtz_object().write("tst_sys_absent.mtz")
  args = ["tst_sys_absent.mtz"]
  try :
    show_systematic_absences.run(args=args, out=out)
  except Sorry :
    pass
  else :
    raise Exception_expected

if (__name__ == "__main__"):
  exercise()
  print("OK")
