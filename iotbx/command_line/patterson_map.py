# LIBTBX_SET_DISPATCHER_NAME cctbx.patterson_map

from __future__ import division
from libtbx.utils import Sorry, Usage
import libtbx.phil
import math
import os
import sys

master_phil = libtbx.phil.parse("""
data = None
  .type = path
labels = None
  .type = str
map_type = *auto anom native
  .type = choice
scaling = *sigma volume
  .type = choice
resolution_factor = 1/4.
  .type = float
  .optional = False
high_resolution = None
  .type = float
sharpening = False
  .type = bool
min_sigma_ratio = 3.0
  .type = float
  .optional = False
diff_limit = None
  .type = float
remove_origin_peak = False
  .type = bool
map_file_name = None
  .type = path
""")

def run (args, out=None) :
  if (len(args) == 0) or ("--help" in args) :
    raise Usage("""\
cctbx.patterson_map data.mtz [options]

Calculates a simple or anomalous difference Patterson map.  Output is in CCP4
format.

Full options:

%s
""" % master_phil.as_str(prefix="  "))
  import iotbx.phil
  from iotbx import file_reader
  from cctbx.array_family import flex
  from cctbx import miller
  from cctbx import xray
  if (out is None) : out = sys.stdout
  hkl_in = None
  sources = []
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil=master_phil,
    reflection_file_def="data")
  params = cmdline.work.extract()
  if (params.data is None) :
    raise Usage("mmtbx.patterson_map data [options]\n\nFull parameters:\n%s" %
      master_phil.as_str(prefix="  "))
  hkl_in = file_reader.any_file(params.data,
    force_type="hkl",
    raise_sorry_if_errors=True)
  miller_arrays = hkl_in.file_server.miller_arrays
  f_obs = None
  all_fs = []
  all_fs_anom = []
  all_is = []
  all_is_anom = []
  for array in miller_arrays :
    labels = array.info().label_string()
    if (labels == params.labels) :
      f_obs = array
      break
    elif (params.labels is None) :
      if (array.is_xray_amplitude_array()) :
        if (array.anomalous_flag()) :
          all_fs_anom.append(array)
        else :
          all_fs.append(array)
      elif (array.is_xray_intensity_array()) :
        if (array.anomalous_flag()) :
          all_is_anom.append(array)
        else :
          all_is.append(array)
  if (f_obs is None) :
    for choices in [all_fs_anom,all_is_anom,all_fs,all_is] :
      if (len(choices) > 1) :
        raise Sorry(("Multiple equally good candidates for data:\n%s\n"+
          "Please specify labels=<your choice here>.") %
          "\n".join([ a.info().label_string() for a in choices ]))
      elif (len(choices) == 1) :
        f_obs = choices[0]
        print >> out, "Defaulting to data in %s" % f_obs.info().label_string()
        break
  if (f_obs is None) :
    raise Sorry("No suitable data found.")
  if (f_obs.is_xray_intensity_array()) :
    f_obs = f_obs.f_sq_as_f()
  f_obs = f_obs.map_to_asu()
  if (not f_obs.is_unique_set_under_symmetry()) :
    f_obs = f_obs.merge_equivalents().array()
  if (f_obs.sigmas() is not None) :
    f_obs = f_obs.select(f_obs.sigmas() > 0)
    f_obs = f_obs.select((f_obs.data() / f_obs.sigmas()) > params.min_sigma_ratio)
  final_array = f_obs # default
  if (f_obs.anomalous_flag()) :
    if (params.map_type in ["auto", "anom"]) :
      print >> out, "Output will be anomalous Patterson map"
      final_array = f_obs.anomalous_differences()
      final_array = abs(final_array)
      if (params.diff_limit is not None) :
        final_array = final_array.select(final_array.data() < params.diff_limit)
    else :
      final_array = f_obs.average_bijvoet_mates()
  final_array.setup_binner(auto_binning=True)
  e = final_array.quasi_normalize_structure_factors()
  # XXX not sure what this does - see cctbx/regression/tst_miller.py
  #u_base = xray.calc_u_base(e.d_min(), params.resolution_factor)
  #d_star_sq = e.unit_cell().d_star_sq(e.indices())
  #dw = flex.exp(d_star_sq*2*(math.pi**2)*u_base)
  #eb = miller.array(miller_set=e, data=e.data()/dw)
  map = e.patterson_map(
    resolution_factor=params.resolution_factor,
    d_min=params.high_resolution,
    sharpening=params.sharpening,
    origin_peak_removal=params.remove_origin_peak)
  if (params.scaling == "sigma") :
    map.apply_sigma_scaling()
  else :
    map.apply_volume_scaling()
  if (params.map_file_name is None) :
    base = os.path.splitext(os.path.basename(params.data))[0]
    params.map_file_name = base + "_patt.ccp4"
  map.as_ccp4_map(
    file_name=params.map_file_name)
  print >> out, "Wrote %s" % params.map_file_name

def validate_params (params) :
  if (params.data is None) :
    raise Sorry("Data file not specified.")
  elif (not os.path.isfile(params.data)) :
    raise Sorry("The path %s is not a valid file." % params.data)
  return True

if (__name__ == "__main__") :
  run(sys.argv[1:])
