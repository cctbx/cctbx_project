from __future__ import absolute_import, division, print_function

import iotbx.symmetry
from cctbx import sgtbx, uctbx
from libtbx.test_utils import Exception_expected
from libtbx.utils import Sorry
from six.moves import cStringIO as StringIO

def exercise():
  m = iotbx.symmetry.manager(prefer_pdb_space_group=True)
  (uc_mismatch, sg_mismatch) = m.add_reflections_file(
    file_name="data.mtz",
    space_group=sgtbx.space_group_info("P222"),
    unit_cell=uctbx.unit_cell("50 60 70 90 90 90"))
  assert (m.get_current_as_strings() == ('P 2 2 2', '50 60 70 90 90 90'))
  (uc_mismatch, sg_mismatch) = m.add_pdb_file(
    file_name="model.pdb",
    space_group=sgtbx.space_group_info("P212121"),
    unit_cell=uctbx.unit_cell("50 60 70 90 90 90"))
  assert (not (uc_mismatch or sg_mismatch))
  (uc_mismatch, sg_mismatch) = m.add_pdb_file(
    file_name="reference_model.pdb",
    space_group=sgtbx.space_group_info("P63"),
    unit_cell=uctbx.unit_cell("40 40 75 90 90 120"))
  assert ((uc_mismatch, sg_mismatch) == (True, True))
  assert (m.get_current_as_strings() == ('P 21 21 21', '50 60 70 90 90 90'))
  (uc_mismatch, sg_mismatch) = m.add_reflections_file(
    file_name="data_neutron.mtz",
    space_group=sgtbx.space_group_info("P222"),
    unit_cell=uctbx.unit_cell("50.1 60 70.1 90 90 90"))
  assert (not (uc_mismatch or sg_mismatch))
  (uc_mismatch, sg_mismatch) = m.add_reflections_file(
    file_name="data_rfree.hkl",
    space_group=None,
    unit_cell=None)
  assert (not (uc_mismatch or sg_mismatch))
  assert (m.get_current_as_strings() == ('P 21 21 21', '50 60 70 90 90 90'))
  assert (m.check_cell_compatibility("phenix.refine"))
  symm_choices = m.get_symmetry_choices()
  assert (symm_choices.space_group_files == [('model.pdb', 'P 21 21 21'),
    ('reference_model.pdb', 'P 63'), ('data.mtz', 'P 2 2 2'),
    ('data_neutron.mtz', 'P 2 2 2')])
  assert (symm_choices.unit_cell_files == [
    ('model.pdb', '(50, 60, 70, 90, 90, 90)'),
    ('reference_model.pdb', '(40, 40, 75, 90, 90, 120)'),
    ('data.mtz', '(50, 60, 70, 90, 90, 90)'),
    ('data_neutron.mtz', '(50.1, 60, 70.1, 90, 90, 90)')])
  m.set_current_as_strings("P63", "50 60 70 90 90 90")
  try :
    m.check_cell_compatibility(
      program_name="phenix.refine",
      raise_error_if_incomplete=True)
  except Sorry :
    pass
  else :
    raise Exception_expected
  out = StringIO()
  m.show(out=out)
  assert (out.getvalue() == """\
model.pdb: (50, 60, 70, 90, 90, 90) P 21 21 21
reference_model.pdb: (40, 40, 75, 90, 90, 120) P 63
data.mtz: (50, 60, 70, 90, 90, 90) P 2 2 2
data_neutron.mtz: (50.1, 60, 70.1, 90, 90, 90) P 2 2 2
data_rfree.hkl: None None
""")

if (__name__ == "__main__"):
  exercise()
  print("OK")
