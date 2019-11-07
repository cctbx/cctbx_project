
from __future__ import absolute_import, division, print_function
from mmtbx.command_line import massage_data
from iotbx import file_reader
from cctbx.development import random_structure
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
import os.path as op
import random
from six.moves import zip

def exercise_twin_detwin():
  random.seed(12345)
  flex.set_random_seed(12345)
  xrs = random_structure.xray_structure(
    unit_cell=(12,5,12,90,90,90),
    space_group_symbol="P1",
    n_scatterers=12,
    elements="random")
  fc = abs(xrs.structure_factors(d_min=1.5).f_calc())
  fc = fc.set_observation_type_xray_amplitude()
  mtz_file = "tmp_massage_in.mtz"
  fc.as_mtz_dataset(column_root_label="F").mtz_object().write(mtz_file)
  massage_data.run(
    args=[
      mtz_file,
      "aniso.action=None",
      "outlier.action=None",
      "symmetry.action=twin",
      "twin_law='l,-k,h'",
      "fraction=0.3",
      "hklout=tmp_massage_twinned.mtz",
    ],
    out=null_out())
  assert op.isfile("tmp_massage_twinned.mtz")
  mtz_in = file_reader.any_file("tmp_massage_twinned.mtz")
  fc_twin = mtz_in.file_server.miller_arrays[0].f_sq_as_f()
  fc_twin, fc_tmp = fc_twin.common_sets(other=fc)
  for hkl, f1, f2 in zip(fc_tmp.indices(), fc_tmp.data(), fc_twin.data()):
    if (abs(hkl[0]) != abs(hkl[2])):
      assert not approx_equal(f1, f2, eps=0.01, out=null_out()), (hkl, f1, f2)
  massage_data.run(
    args=[
      mtz_file,
      "aniso.action=None",
      "outlier.action=None",
      "symmetry.action=twin",
      "twin_law='l,-k,h'",
      "fraction=0.3",
      "hklout=tmp_massage_twinned.sca",
    ],
    out=null_out())
  assert op.isfile("tmp_massage_twinned.sca")
  massage_data.run(
    args=[
      "tmp_massage_twinned.mtz",
      "aniso.action=None",
      "outlier.action=None",
      "symmetry.action=detwin",
      "twin_law='l,-k,h'",
      "fraction=0.3",
      "hklout=tmp_massage_detwinned.mtz",
    ],
    out=null_out())
  mtz_in = file_reader.any_file("tmp_massage_detwinned.mtz")
  fc_detwin = mtz_in.file_server.miller_arrays[0].f_sq_as_f()
  fc_detwin, fc_tmp = fc_detwin.common_sets(other=fc)
  # XXX we appear to lose some accuracy here, possibly due to the use of
  # MTZ format
  for hkl, f1, f2 in zip(fc_tmp.indices(), fc_tmp.data(), fc_detwin.data()):
    assert approx_equal(f1, f2, eps=0.01), hkl

if (__name__ == "__main__"):
  exercise_twin_detwin()
  print("OK")
