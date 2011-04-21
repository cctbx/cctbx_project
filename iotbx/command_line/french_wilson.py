# LIBTBX_SET_DISPATCHER_NAME phenix.french_wilson

import libtbx.phil
from libtbx import runtime_utils
from libtbx.utils import Sorry, Usage
import os
import sys

master_phil = libtbx.phil.parse("""
french_wilson {
  file_name = None
    .type = path
    .short_caption = Reflections
    .help = '''input intensity data file (mtz)'''
    .style = bold file_type:hkl process_hkl child:iobs:intensity_labels \
      child:rfree:r_free_flags.label
  intensity_labels = None
    .type = strings
    .style = bold renderer:draw_fobs_label_widget
  r_free_flags.label = None
    .type = str
    .short_caption = R-free label
    .style = bold renderer:draw_rfree_label_widget
  output_file = None
    .type = path
    .optional = True
    .help = '''Enter a .mtz output name'''
    .style = bold
  include scope libtbx.phil.interface.tracking_params
  keep_r_free_flags = True
    .type = bool
    .help = "Keep R-free flag data if present"
    .short_caption = Keep R-free flags if present in input file
  include scope cctbx.french_wilson.master_phil
}
""", process_includes=True)

def run (args, out=sys.stdout) :
  from cctbx import french_wilson
  from iotbx import file_reader
  hkl_file = None
  sources = []
  interpreter = master_phil.command_line_argument_interpreter()
  for arg in args :
    if os.path.isfile(arg) :
      input_file = file_reader.any_file(arg)
      if (input_file.file_type == "hkl") :
        hkl_file = input_file
        sources.append(interpreter.process(arg="file_name=\"%s\"" % arg))
      elif (input_file.file_type == "phil") :
        sources.append(input_file.file_object)
    else :
      arg_phil = interpreter.process(arg=arg)
      sources.append(arg_phil)
  work_phil = master_phil.fetch(sources=sources)
  work_params = work_phil.extract()
  if (work_params.french_wilson.file_name is None) :
    if (hkl_file is None) :
      raise Usage("phenix.french_wilson data.mtz [params.eff] [options ...]")
    else :
      work_params.french_wilson.file_name = hkl_file.file_name
  elif (hkl_file is None) :
    hkl_file = file_reader.any_file(work_params.french_wilson.file_name)
  params = work_params.french_wilson
  xray_data_server = hkl_file.file_server
  crystal_symmetry = xray_data_server.miller_arrays[0].crystal_symmetry()
  if (crystal_symmetry is None) :
    raise Sorry("No crystal symmetry found.  This program requires an input "+
      "format with complete symmetry information.")
  unit_cell = xray_data_server.miller_arrays[0].unit_cell()
  if (unit_cell is None) :
    raise Sorry("No unit cell found.  This program requires an input "+
      "format with complete unit cell information.")
  i_obs = None
  i_obs = xray_data_server.get_xray_data(
    file_name = params.file_name,
    labels = params.intensity_labels,
    ignore_all_zeros = True,
    parameter_scope = 'french_wilson',
    parameter_name = 'intensity_labels')
  import cStringIO
  xray_data_server.err = cStringIO.StringIO()
  try :
    r_free_flags, test_flag_value = xray_data_server.get_r_free_flags(
      file_name = params.file_name,
      label = params.r_free_flags.label,
      test_flag_value = None,
      disable_suitability_test = False,
      parameter_scope = "french_wilson.r_free_flags")
  except Sorry, e :
    r_free_flags = None
  if (i_obs is None) :
    raise Sorry("Couldn't find intensities!")
  f_obs = french_wilson.french_wilson_scale(miller_array=i_obs,
    params=params,
    log=out)
  if f_obs is None:
    raise Sorry("Not enough data to accurately apply the French-Wilson method."+\
                " Exiting.")
  if params.output_file == None:
    output_file = "french_wilson.mtz"
  else:
    output_file = params.output_file
  mtz_dataset = i_obs.as_mtz_dataset(
    column_root_label = "I")
  mtz_dataset.add_miller_array(
    miller_array      = f_obs,
    column_root_label = "F")
  if (r_free_flags is not None) and (params.keep_r_free_flags):
    mtz_dataset.add_miller_array(
      miller_array      = r_free_flags,
      column_root_label = "R-free-flags")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = output_file)
  return output_file

class launcher (runtime_utils.simple_target) :
  def __call__ (self) :
    return run(args=list(self.args), out=sys.stdout)

def validate_params (params) :
  if (params.french_wilson.file_name is None) :
    raise Sorry("Please specify a reflections file.")
  elif (params.french_wilson.intensity_labels is None) :
    raise Sorry("No intensity labels selected; are you sure the input file "+
      "contains appropriate data?")
  return True

def finish_job (result) :
  output_files = []
  if (result is not None) :
    output_files.append((result, "Corrected amplitudes"))
  return (output_files, [])

if(__name__ == "__main__"):
  run(sys.argv[1:])
