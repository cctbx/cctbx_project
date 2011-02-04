# LIBTBX_SET_DISPATCHER_NAME phenix.french_wilson

from cctbx import french_wilson
import libtbx.phil
import iotbx.utils
from iotbx import reflection_file_utils
from iotbx import crystal_symmetry_from_any
import sys

master_phil = libtbx.phil.parse("""
french_wilson {
  file_name = None
    .type = path
    .help = '''input intensity data file (mtz)'''
  intensity_labels = None
    .type = str
  r_free_label = None
    .type = str
  output_file = None
    .type = path
    .optional = True
    .help = '''Enter a .mtz output name'''
  keep_r_free_flags = True
    .type = bool
    .help = '''Keep R-free flag data if present'''
  include scope cctbx.french_wilson.master_phil
}
""", process_includes=True)

def run(args):
  input_objects = iotbx.utils.process_command_line_inputs(
    args=args,
    master_phil=master_phil,
    input_types=("mtz",))
  work_phil = master_phil.fetch(sources=input_objects["phil"])
  work_params = work_phil.extract()
  if work_params.french_wilson.file_name == None:
    if len(input_objects["mtz"]) != 1:
      raise Usage("phenix.french_wilson data.mtz")
    file_obj = input_objects["mtz"][0]
    work_params.french_wilson.file_name = file_obj.file_name
  params = work_params.french_wilson
  crystal_symmetry = crystal_symmetry_from_any.extract_from(params.file_name)
  if(crystal_symmetry is None):
    raise Sorry("No crystal symmetry found.")
  xray_data_server =  reflection_file_utils.reflection_file_server(
      crystal_symmetry = crystal_symmetry,
      force_symmetry = True,
      reflection_files=[])
  i_obs = None
  i_obs = xray_data_server.get_xray_data(
      file_name = params.file_name,
      labels = None,
      ignore_all_zeros = True,
      parameter_scope = '',
      parameter_name = 'obs_labels'
  )
  r_free_flags = None
  r_free_flags_array = xray_data_server.get_r_free_flags(
      file_name = params.file_name,
      label = None,
      test_flag_value = None,
      disable_suitability_test = False,
      parameter_scope = '',
      return_all_valid_arrays=True
  )
  #TODO - improve handling of R-free arrays
  # currently only works for the simple case
  if len(r_free_flags_array) == 1:
    if r_free_flags_array[0][1] == 1:
      r_free_flags = r_free_flags_array[0][0]
  assert (i_obs is not None), "Couldn't find intensities!"
  f_obs = french_wilson.french_wilson_scale(miller_array=i_obs,
                                            params=params)
  if params.output_file == None:
    output_file = "french_wilson.mtz"
  else:
    output_file = params.output_file
  mtz_dataset = i_obs.as_mtz_dataset(
    column_root_label = "I")
  mtz_dataset.add_miller_array(
    miller_array      = f_obs,
    column_root_label = "F")
  if(r_free_flags is not None and params.keep_r_free_flags):
    mtz_dataset.add_miller_array(
      miller_array      = r_free_flags,
      column_root_label = "R-free-flags")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = output_file)

if(__name__ == "__main__"):
  run(sys.argv[1:])
