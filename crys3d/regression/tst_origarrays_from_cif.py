from __future__ import absolute_import, division, print_function
from crys3d.hklviewer.hklview_frame import origarrays_from_cif_block
from libtbx.test_utils import approx_equal
from cctbx.array_family import flex
from iotbx import cif

# Two wavelength ids, one crystal and one scale group; a '?' in the
# second-wavelength rows of F_meas_au and a '.' for every intensity there.
# The wavelength table of one block carries over to later blocks in
# iotbx.cif, so the single-id case below is a separate cif text.
multi_text = """\
data_multi
_cell.length_a 10
_cell.length_b 11
_cell.length_c 12
_cell.angle_alpha 90
_cell.angle_beta 90
_cell.angle_gamma 90
_symmetry.space_group_name_H-M 'P 1'
loop_
_diffrn_radiation_wavelength.id
_diffrn_radiation_wavelength.wavelength
1 0.9795
2 0.9793
loop_
_refln.crystal_id
_refln.wavelength_id
_refln.scale_group_code
_refln.index_h
_refln.index_k
_refln.index_l
_refln.status
_refln.F_meas_au
_refln.F_meas_sigma_au
_refln.intensity_meas
1 1 1 1 0 0 o 10.5 1.5 110.25
1 1 1 2 0 0 f 20.5 2.5 420.25
1 1 1 3 0 0 o 30.5 3.5 930.25
1 2 1 1 0 0 o 11.5 1.6 .
1 2 1 2 0 0 o ? 2.6 .
1 2 1 3 0 0 f 31.5 3.6 .
data_nodata
_cell.length_a 10
_cell.length_b 11
_cell.length_c 12
_cell.angle_alpha 90
_cell.angle_beta 90
_cell.angle_gamma 90
_symmetry.space_group_name_H-M 'P 1'
"""

single_text = """\
data_single
_cell.length_a 10
_cell.length_b 11
_cell.length_c 12
_cell.angle_alpha 90
_cell.angle_beta 90
_cell.angle_gamma 90
_symmetry.space_group_name_H-M 'P 1'
loop_
_refln.crystal_id
_refln.wavelength_id
_refln.scale_group_code
_refln.index_h
_refln.index_k
_refln.index_l
_refln.status
_refln.F_meas_au
_refln.F_meas_sigma_au
1 1 1 1 0 0 o 10.5 1.5
1 1 1 2 0 0 f 20.5 2.5
1 1 1 3 0 0 o . 3.5
"""

def check_type(table, key, flex_type):
  assert key in table, (key, sorted(table.keys()))
  assert type(table[key]) is flex_type, (key, type(table[key]).__name__)

def exercise_multi_id():
  block = cif.reader(input_string=multi_text).model()["multi"]
  table = origarrays_from_cif_block(block)
  expected_keys = ["HKLs"]
  for label in ["crystal_id", "wavelength_id", "scale_group_code", "status",
                "F_meas_au", "F_meas_sigma_au", "intensity_meas"]:
    expected_keys.append(label + ",wavelength_id=1")
    expected_keys.append(label + ",wavelength_id=2")
  assert sorted(table.keys()) == sorted(expected_keys), sorted(table.keys())
  check_type(table, "HKLs", flex.miller_index)
  assert list(table["HKLs"]) == [
    (1, 0, 0), (2, 0, 0), (3, 0, 0), (1, 0, 0), (2, 0, 0), (3, 0, 0)]
  # Numeric columns: the rows of that wavelength only.
  check_type(table, "F_meas_au,wavelength_id=1", flex.double)
  assert approx_equal(table["F_meas_au,wavelength_id=1"], [10.5, 20.5, 30.5])
  check_type(table, "F_meas_sigma_au,wavelength_id=2", flex.double)
  assert approx_equal(table["F_meas_sigma_au,wavelength_id=2"], [1.6, 2.6, 3.6])
  check_type(table, "intensity_meas,wavelength_id=1", flex.double)
  assert approx_equal(table["intensity_meas,wavelength_id=1"], [110.25, 420.25, 930.25])
  check_type(table, "crystal_id,wavelength_id=2", flex.int)
  assert list(table["crystal_id,wavelength_id=2"]) == [1, 1, 1]
  check_type(table, "wavelength_id,wavelength_id=2", flex.int)
  assert list(table["wavelength_id,wavelength_id=2"]) == [2, 2, 2]
  # A null anywhere in the selected rows keeps the whole column as strings.
  check_type(table, "F_meas_au,wavelength_id=2", flex.std_string)
  assert list(table["F_meas_au,wavelength_id=2"]) == [
    "10.5", "20.5", "30.5", "11.5", "?", "31.5"]
  check_type(table, "intensity_meas,wavelength_id=2", flex.std_string)
  assert list(table["intensity_meas,wavelength_id=2"]) == [
    "110.25", "420.25", "930.25", ".", ".", "."]
  check_type(table, "status,wavelength_id=1", flex.std_string)
  assert list(table["status,wavelength_id=1"]) == ["o", "f", "o", "o", "o", "f"]

def exercise_single_id():
  block = cif.reader(input_string=single_text).model()["single"]
  table = origarrays_from_cif_block(block)
  assert sorted(table.keys()) == sorted([
    "HKLs", "crystal_id", "wavelength_id", "scale_group_code", "status",
    "F_meas_au", "F_meas_sigma_au"]), sorted(table.keys())
  check_type(table, "HKLs", flex.miller_index)
  assert table["HKLs"].size() == 3
  check_type(table, "crystal_id", flex.int)
  assert list(table["crystal_id"]) == [1, 1, 1]
  check_type(table, "F_meas_sigma_au", flex.double)
  assert approx_equal(table["F_meas_sigma_au"], [1.5, 2.5, 3.5])
  check_type(table, "F_meas_au", flex.std_string)
  assert list(table["F_meas_au"]) == ["10.5", "20.5", "."]
  check_type(table, "status", flex.std_string)
  assert list(table["status"]) == ["o", "f", "o"]

def exercise_no_reflections():
  block = cif.reader(input_string=multi_text).model()["nodata"]
  table = origarrays_from_cif_block(block)
  assert len(table) == 0, sorted(table.keys())

def run():
  exercise_multi_id()
  exercise_single_id()
  exercise_no_reflections()
  print("OK")

if __name__ == "__main__":
  run()
