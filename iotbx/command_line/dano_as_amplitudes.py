"""Convert DANO values (F+ - F-) to pseudo-amplitudes"""
# TODO tests

from __future__ import absolute_import, division, print_function
from iotbx.reflection_file_utils import reflection_file_server
import iotbx.phil
from libtbx.utils import Sorry
import sys

master_phil_str = """
data = None
  .type = path
labels = None
  .type = strings
symmetry = None
  .type = path
mtz_out = dano_as_fobs.mtz
  .type = path
"""

def run(args):
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil_string=master_phil_str,
    reflection_file_def="data",
    pdb_file_def="symmetry")
  params = cmdline.work.extract()
  if (params.data is None):
    raise Sorry("Data file not defined")
  hkl_file = cmdline.get_file(params.data)
  crystal_symmetry = None
  if (params.symmetry is not None):
    from iotbx import crystal_symmetry_from_any
    crystal_symmetry = crystal_symmetry_from_any.extract_from(
      file_name=params.symmetry)
  hkl_server = reflection_file_server(
    crystal_symmetry=crystal_symmetry,
    force_symmetry=True,
    reflection_files=[hkl_file.file_object],
    err=sys.stderr)
  data = hkl_server.get_xray_data(
    file_name=params.data,
    labels=params.labels,
    ignore_all_zeros=True,
    parameter_scope="",
    minimum_score=4,
    prefer_anomalous=True)
  if (not data.anomalous_flag()):
    raise Sorry("Must provide anomalous data.")
  if data.is_xray_intensity_array():
    data = data.f_sq_as_f()
  dano = abs(data.anomalous_differences())
  dano.set_observation_type_xray_amplitude()
  dano.as_mtz_dataset(column_root_label="F").mtz_object().write(params.mtz_out)
  print("Wrote DANO to %s" % params.mtz_out)

if (__name__ == "__main__"):
  run(sys.argv[1:])
