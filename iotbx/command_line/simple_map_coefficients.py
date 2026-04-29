"""Create map coefficients from a reflection file"""

from __future__ import absolute_import, division, print_function
from libtbx.utils import Sorry, Usage
from libtbx import Auto
import os
import sys

master_phil = """
file_name = None
  .type = path
  .short_caption = MTZ file
  .style = file_type:hkl
data_labels = None
  .type = str
  .short_caption = Data labels
  .style = None
phase_labels = None
  .type = str
weight_labels = None
  .type = str
map_type = *Fo anom
  .type = choice
use_weights = Auto
  .type = bool
output_file = None
  .type = path
"""

def run(args, out=sys.stdout):
  from iotbx import file_reader
  import iotbx.phil
  if (len(args) == 0):
    raise Usage("""\
iotbx.simple_map_coefficients data_phases.mtz [options]

Full parameters:
%s""" % iotbx.phil.parse(master_phil).as_str(attributes_level=1, prefix=" "))
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil_string=master_phil,
    reflection_file_def="file_name")
  params = cmdline.work.extract()
  if (params.file_name is None):
    raise Sorry("No reflection file specified.")
  hkl_in = file_reader.any_file(params.file_name).check_file_type("hkl")
  hkl_server = hkl_in.file_server
  data = hkl_server.get_xray_data(
    file_name=params.file_name,
    labels=params.data_labels,
    ignore_all_zeros=False,
    parameter_name="data_labels",
    parameter_scope="",
    prefer_anomalous=True,
    prefer_amplitudes=True)
  data_labels = data.info().label_string()
  if (data.is_xray_intensity_array()):
    from cctbx.french_wilson import french_wilson_scale
    data = french_wilson_scale(data, log=out)
  phases = hkl_server.get_phases_deg(
    file_name=params.file_name,
    labels=params.phase_labels,
    convert_to_phases_if_necessary=True,
    original_phase_units=None,
    parameter_scope="",
    parameter_name="phase_labels",
    minimum_score=2)
  assert (not phases.anomalous_flag())
  deg = True # FIXME
  weights = None
  if (params.use_weights in [Auto, True]):
    # FIXME centralize this in iotbx.reflection_file_utils
    for array in hkl_server.miller_arrays :
      if (array.is_real_array()):
        label_string = array.info().label_string()
        if ((label_string == params.weight_labels) or
            ((params.weight_labels is None) and ("FOM" in label_string))):
          weights = array
          break
  amplitudes = data
  if (params.map_type == "anom"):
    if (not data.anomalous_flag()):
      raise Sorry("Anomalous map requested, but selected data are merged.")
    amplitudes = data.anomalous_differences()
    print("Using anomalous differences in %s" % data_labels, file=out)
  else :
    print("Using amplitudes in %s" % data_labels, file=out)
    if (data.anomalous_flag()):
      amplitudes = data.average_bijvoet_mates()
  if (params.use_weights is Auto) and (weights is not None):
    if (params.map_type != "anom"):
      params.use_weights = True
  elif (params.use_weights == True) and (weights is None):
    raise Sorry("No weights (FOM, etc.) found in input file.")
  if (params.use_weights == True):
    assert (not weights.anomalous_flag())
    print("Applying weights in %s" % weights.info().label_string(), file=out)
    amplitudes, weights = amplitudes.common_sets(other=weights)
    amplitudes = amplitudes.customized_copy(
      data=amplitudes.data()*weights.data())
  amplitudes = amplitudes.customized_copy(sigmas=None)
  print("Applying phases in %s" % phases.info().label_string(), file=out)
  amplitudes, phases = amplitudes.common_sets(phases)
  coeffs = amplitudes.phase_transfer(phases,
    deg=deg).set_observation_type(None) # FIXME
  if (params.map_type == "anom") : # apply 90-degree phase shift
    coeffs = coeffs.customized_copy(data=coeffs.data()/(2j))
  assert (coeffs.is_complex_array())
  column_root_label = "F"
  decorator = None
  if (params.map_type == "anom"):
    column_root_label = "ANOM"
  elif (params.use_weights == True):
    column_root_label = "FWT"
    decorator = iotbx.mtz.ccp4_label_decorator()
  import iotbx.mtz
  mtz_dataset = coeffs.as_mtz_dataset(
    column_root_label=column_root_label,
    label_decorator=decorator)
  if (params.output_file is None):
    params.output_file = "map_coeffs.mtz"
  mtz_dataset.mtz_object().write(params.output_file)
  print("Wrote %s" % params.output_file, file=out)
  return os.path.abspath(params.output_file)

if (__name__ == "__main__"):
  run(sys.argv[1:])
