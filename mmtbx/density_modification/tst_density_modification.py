from libtbx import easy_run
from mmtbx.command_line import density_modification

import os

params = """\
density_modification {
 input {
   reflection_data.file_name = %s/data/scale.hkl
   reflection_data.labels = f_nat,s_nat
   experimental_phases.file_name = %s/phasing/ir_phase.hkl
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
}
"""

def exercise_density_modification():
  cci_structure_lib = os.environ.get("CCI_STRUCTURE_LIB")
  if cci_structure_lib is None:
    print "Skipping exercise_density_modification(): $CCI_STRUCTURE_LIB is not set"
    return
  rnase_s_path = os.path.join(cci_structure_lib, "rnase-s")
  from libtbx.test_utils import open_tmp_file
  tmp_file = open_tmp_file(suffix=".params")
  tmp_file.write(params % (rnase_s_path, rnase_s_path))
  tmp_file.close()
  #density_modification.run(
    #args=[tmp_file.name, "%s/model/1RGE.pdb"%rnase_s_path])

  cmd = "mmtbx.density_modification %s %s/model/1RGE.pdb" %(
    tmp_file.name, rnase_s_path)
  print cmd
  result = easy_run.fully_buffered(command=cmd).raise_if_errors()
  assert result.stdout_lines[-4] == 'Starting dm/model correlation: 0.631999'
  assert result.stdout_lines[-3] == 'Final dm/model correlation:    0.777594'

def run():
  exercise_density_modification()

if __name__ == '__main__':
  run()
  print "OK"
