from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.map_to_structure_factors

from libtbx.test_utils import approx_equal
import iotbx.ccp4_map
from cctbx import miller, crystal
from cctbx.array_family import flex
import mmtbx.utils
import sys
from libtbx.utils import Sorry

master_params_str = """
output_file_name = map_to_structure_factors.mtz
  .type=str
d_min = None
  .type=float
k_blur = 1
  .type = float
b_blur  = 100
  .type = float
"""

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=True)

def broadcast(m, log):
  print >> log, "-"*79
  print >> log, m
  print >> log, "*"*len(m)

def run(args, log=None):
  inputs = mmtbx.utils.process_command_line_args(args = args,
    master_params = master_params())
  got_map = False
  broadcast(m="Parameters:", log=log)
  inputs.params.show(prefix="  ")
  params = inputs.params.extract()
  if(inputs.ccp4_map is not None):
    broadcast(m="Processing input CCP4 map file: %s"%inputs.ccp4_map_file_name,
      log=log)
    ccp4_map = inputs.ccp4_map
    ccp4_map.show_summary(prefix="  ")
    got_map = True
  if(not got_map):
    raise Sorry("Map file is needed.")
  #
  m = ccp4_map
  broadcast(m="Input map information:", log=log)
  print "m.all()   :", m.data.all()
  print "m.focus() :", m.data.focus()
  print "m.origin():", m.data.origin()
  print "m.nd()    :", m.data.nd()
  print "m.size()  :", m.data.size()
  print "m.focus_size_1d():", m.data.focus_size_1d()
  print "m.is_0_based()   :", m.data.is_0_based()
  print "map: min/max/mean:", flex.min(m.data), flex.max(m.data), flex.mean(m.data)
  print "unit cell:", m.unit_cell_parameters
  #
  if(not m.data.is_0_based()):
    raise Sorry("Map must have origin at (0,0,0): recenter the map and try again.")
  # generate complete set of Miller indices up to given high resolution d_min
  cs = crystal.symmetry(m.unit_cell_parameters, 1)
  if(params.d_min is None):
    raise Sorry(
      "High resolution limit (d_min) for structure factors calculation must be given.")
  complete_set = miller.build_set(
    crystal_symmetry = cs,
    anomalous_flag   = False,
    d_min            = params.d_min)
  broadcast(m="Complete set information:", log=log)
  complete_set.show_comprehensive_summary(prefix="  ")
  try:
    f_obs_cmpl = complete_set.structure_factors_from_map(
      map            = m.data.as_double(),
      use_scale      = True,
      anomalous_flag = False,
      use_sg         = True)
  except Exception, e:
    if(str(e) == "cctbx Error: Miller index not in structure factor map."):
      msg = "Too high resolution requested. Try running with larger d_min."
      raise Sorry(msg)
  #
  mtz_dataset = f_obs_cmpl.as_mtz_dataset(column_root_label="F")
  mtz_dataset.add_miller_array(
    miller_array      = abs(f_obs_cmpl),
    column_root_label = "F_ampl")
  # convert phases into HL coefficeints
  f_model_phases = f_obs_cmpl.phases().data()
  sin_f_model_phases = flex.sin(f_model_phases)
  cos_f_model_phases = flex.cos(f_model_phases)
  ss = 1./flex.pow2(f_obs_cmpl.d_spacings().data()) / 4.
  t = 2*params.k_blur * flex.exp(-params.b_blur*ss)
  hl_a_model = t * cos_f_model_phases
  hl_b_model = t * sin_f_model_phases
  hl_data = flex.hendrickson_lattman(a = hl_a_model, b = hl_b_model)
  hl = f_obs_cmpl.customized_copy(data = hl_data)
  mtz_dataset.add_miller_array(
    miller_array      = hl,
    column_root_label = "HL")
  # write output MTZ file with all the data
  broadcast(m="Writing output MTZ file:", log=log)
  print >> log, "  file name:", params.output_file_name
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = params.output_file_name)
  # sanity check
  map_coeffs = abs(f_obs_cmpl).phase_transfer(phase_source = hl)
  assert approx_equal(map_coeffs.map_correlation(other=f_obs_cmpl), 1)

if(__name__ == "__main__"):
  run(sys.argv[1:])
  print "All done."
