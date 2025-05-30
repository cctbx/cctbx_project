"""Show the systematic absences in a reflection file"""
# LIBTBX_SET_DISPATCHER_NAME cctbx.show_systematic_absences

from __future__ import absolute_import, division, print_function
from iotbx.reflection_file_utils import reflection_file_server
import iotbx.phil
from cctbx import crystal
from libtbx.utils import Sorry
import sys

master_phil_str = """
data = None
  .type = path
labels = None
  .type = strings
space_group = None
  .type = space_group
unit_cell = None
  .type = unit_cell
symmetry = None
  .type = path
"""

def run(args, out=sys.stdout):
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil_string=master_phil_str,
    reflection_file_def="data",
    pdb_file_def="symmetry",
    space_group_def="space_group",
    unit_cell_def="unit_cell",
    usage_string="""\
iotbx.show_systematic_absences data.hkl [space_group] [unit_cell]""")
  params = cmdline.work.extract()
  if (params.data is None):
    raise Sorry("Data file not specified.")
  hkl_file = cmdline.get_file(params.data)
  hkl_server = reflection_file_server(
    crystal_symmetry=None,
    force_symmetry=True,
    reflection_files=[hkl_file.file_object],
    err=sys.stderr)
  data = hkl_server.get_xray_data(
    file_name=params.data,
    labels=params.labels,
    ignore_all_zeros=False,
    parameter_scope="",
    minimum_score=4,
    prefer_amplitudes=False)
  symm = data.crystal_symmetry()
  space_group = unit_cell = None
  if (params.symmetry is not None):
    from iotbx import crystal_symmetry_from_any
    symm = crystal_symmetry_from_any.extract_from(file_name=params.symmetry)
  if (symm is not None):
    space_group = symm.space_group_info()
    unit_cell = symm.unit_cell()
  if (space_group is None):
    if (params.space_group is not None):
      space_group = params.space_group
    else :
      raise Sorry("No space group defined.")
  if (unit_cell is None):
    if (params.unit_cell is not None):
      unit_cell = params.unit_cell
    else :
      raise Sorry("No unit cell defined.")
  symm = crystal.symmetry(
    space_group_info=space_group,
    unit_cell=unit_cell)
  data = data.customized_copy(crystal_symmetry=symm)
  if (data.sigmas() is None):
    raise Sorry("Input data are missing experimental sigmas.")
  data.show_all_possible_systematic_absences(out=out)

if (__name__ == "__main__"):
  run(sys.argv[1:])
