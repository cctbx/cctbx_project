from __future__ import absolute_import, division, print_function
from libtbx import easy_run
from libtbx.test_utils import approx_equal
from mmtbx.command_line import density_modification
import libtbx.load_env

import os
op = os.path

params = """\
density_modification {
  input {
    reflection_data.file_name = %s
    reflection_data.labels = f_nat,s_nat
    experimental_phases.file_name = %s
    unit_cell = 64.897,78.323,38.792,90,90,90
    space_group = P212121
  }
  protein_solvent_ratio = 1.31
  density_truncation {
    fraction_min = 0.35
    fraction_max = 1
  }
  solvent_modification {
    method = *flipping flattening
  }
  solvent_adjust = True
  solvent_mask {
    averaging_radius {
      initial = 3.5
      final = 2.5
    }
  }
  solvent_fraction = 0.45
  initial_steps = 10
  shrink_steps = 20
  final_steps = 10
  d_min = 2.5
  grid_resolution_factor = 1/4
  asu_contents {
    n_residues = 96
  }
  anisotropic_correction = True
}
"""

def exercise_density_modification():
  pdb_path = libtbx.env.find_in_repositories(
      relative_path="phenix_regression/density_modification/1RGE.pdb",
      test=os.path.isfile)

  reflection_path = libtbx.env.find_in_repositories(
      relative_path="phenix_regression/density_modification/scale.hkl",
      test=os.path.isfile)

  phases_path = libtbx.env.find_in_repositories(
      relative_path="phenix_regression/density_modification/ir_phase.hkl",
      test=os.path.isfile)

  if pdb_path is None or reflection_path is None or phases_path is None:
    print("Skipping exercise_density_modification(): phenix_regression is not available")
    return

  from libtbx.test_utils import open_tmp_file
  tmp_file = open_tmp_file(suffix=".params")
  tmp_file.write(params % (reflection_path, phases_path))
  tmp_file.close()
  #dm = density_modification.run(
    #args=[tmp_file.name, "%s/model/1RGE.pdb"%rnase_s_path])

  args = (
    tmp_file.name,
    pdb_path)
  for arg in args:
    assert arg.find('"') < 0
  cmd = 'mmtbx.density_modification "%s" "%s"' % args
  print(cmd)
  result = easy_run.fully_buffered(command=cmd).raise_if_errors()
  assert result.stdout_lines[-5].startswith('Starting dm/model correlation:')
  assert result.stdout_lines[-4].startswith('Final dm/model correlation:')
  assert approx_equal(float(result.stdout_lines[-5].split()[-1]), 0.59, 0.01)
  assert approx_equal(float(result.stdout_lines[-4].split()[-1]), 0.80, 0.01)

def run():
  exercise_density_modification()

if __name__ == '__main__':
  run()
  print("OK")
