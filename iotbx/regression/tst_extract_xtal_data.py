from __future__ import absolute_import, division, print_function
from iotbx.reflection_file_reader import any_reflection_file
import iotbx.reflection_file_reader
from iotbx import extract_xtal_data as ed
import sys, os
import cctbx.miller

def run():
  data_dir = os.path.dirname(os.path.abspath(__file__))
  data_mtz = os.path.join(data_dir, 'data', 'phaser_1.mtz')

  params = ed.data_and_flags_master_params().extract()
  params.file_name = data_mtz
  params.labels=['FP,SIGFP']

  phases_params = ed.experimental_phases_master_params().extract()
  phases_params.labels="HLA,HLB,HLC,HLD"

  reflection_files = [
    iotbx.reflection_file_reader.any_reflection_file(data_mtz)]

  server = iotbx.reflection_file_utils.reflection_file_server(
    crystal_symmetry = None,
    reflection_files = reflection_files,
    err              = sys.stdout)

  fo = ed.run(
    reflection_file_server     = server,
    parameters                 = params,
    experimental_phases_params = phases_params,
    extract_r_free_flags       = False).f_obs
  assert cctbx.miller.array(fo, cctbx.miller.array)

if(__name__ == "__main__"):
  run()
