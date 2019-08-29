
"""
Tool for modifying experimental data, spun off from Xtriage.
"""

from __future__ import absolute_import, division, print_function
from libtbx.utils import Sorry
import os.path as op
import sys

master_phil_str = """
crystal_symmetry {
  space_group = None
    .type = space_group
  unit_cell = None
    .type = unit_cell
  symm_file = None
    .type = path
}
input {
  data = None
    .type = path
    .help = Data file
  labels = None
    .type = strings
}
options {
  include scope mmtbx.scaling.massage_twin_detwin_data.master_params
}
output {
  include scope mmtbx.scaling.massage_twin_detwin_data.output_params_str
}
"""

def run(args, out=sys.stdout):
  import iotbx.phil
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil_string=master_phil_str,
    pdb_file_def="crystal_symmetry.symm_file",
    reflection_file_def="input.data",
    space_group_def="crystal_symmetry.space_group",
    unit_cell_def="crystal_symmetry.unit_cell",
    usage_string="""\
mmtbx.massage_data data.mtz [labels=I,SIGI] [options]

Modification of experimental data: remove anisotropy, apply B-factor, filter
negative intensities, add or remove twinning.  For expert use (and extreme
cases) only.
""")
  params = cmdline.work.extract()
  if (params.input.data is None):
    raise Sorry("No data file supplied.")
  from mmtbx.scaling import massage_twin_detwin_data
  from iotbx import crystal_symmetry_from_any
  from iotbx import reflection_file_utils
  from cctbx import crystal
  crystal_symmetry = space_group = unit_cell = None
  if (params.crystal_symmetry.space_group is not None):
    space_group = params.crystal_symmetry.space_group
  if (params.crystal_symmetry.unit_cell is not None):
    unit_cell = params.crystal_symmetry.unit_cell
  crystal_symmetry = None
  if (params.crystal_symmetry.symm_file is not None):
    crystal_symmetry = crystal_symmetry_from_any.extract_from(
      file_name=params.crystal_symmetry.symm_file)
    if (crystal_symmetry is None):
      raise Sorry("No crystal symmetry defined in %s" %
        params.crystal_symmetry.symm_file)
  if (crystal_symmetry is None) and (not None in [space_group, unit_cell]):
    crystal_symmetry = crystal.symmetry(
      space_group_info=space_group,
      unit_cell=unit_cell)
  hkl_in = cmdline.get_file(params.input.data)
  hkl_server = reflection_file_utils.reflection_file_server(
    crystal_symmetry=crystal_symmetry,
    force_symmetry=True,
    reflection_files=[hkl_in.file_object],
    err=sys.stderr)
  data = hkl_server.get_xray_data(
    file_name=params.input.data,
    labels=params.input.labels,
    ignore_all_zeros=True,
    parameter_scope="input",
    prefer_anomalous=True,
    prefer_amplitudes=False)
  result = massage_twin_detwin_data.massage_data(
    miller_array=data,
    parameters=params.options,
    out=out)
  if (params.output.hklout is None):
    file_base = op.splitext(op.basename(params.input.data))[0]
    if (params.output.hklout_type in ["Auto", "mtz"]):
      params.output.hklout = file_base + ".mtz"
    else :
      params.output.hklout = file_base + ".sca"
  result.write_data(
    file_name=params.output.hklout,
    output_type=params.output.hklout_type,
    label_extension=params.output.label_extension)

if (__name__ == "__main__"):
  run(sys.argv[1:])
