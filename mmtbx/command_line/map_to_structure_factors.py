from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.map_to_structure_factors

import iotbx.ccp4_map
from cctbx import miller, crystal
from cctbx.array_family import flex
import mmtbx.utils
import sys
from libtbx.utils import Sorry
from cctbx import maptbx

master_params_str = """
output_file_name = map_to_structure_factors.mtz
  .type=str
d_min = None
  .type=float
k_blur = 1
  .type = float
b_blur  = 100
  .type = float
box = False
  .type = bool
"""

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=True)

def broadcast(m, log):
  print >> log, "-"*79
  print >> log, m
  print >> log, "*"*len(m)

def get_hl(f_obs_cmpl, k_blur, b_blur):
  f_model_phases = f_obs_cmpl.phases().data()
  sin_f_model_phases = flex.sin(f_model_phases)
  cos_f_model_phases = flex.cos(f_model_phases)
  ss = 1./flex.pow2(f_obs_cmpl.d_spacings().data()) / 4.
  t = 2*k_blur * flex.exp(-b_blur*ss)
  hl_a_model = t * cos_f_model_phases
  hl_b_model = t * sin_f_model_phases
  hl_data = flex.hendrickson_lattman(a = hl_a_model, b = hl_b_model)
  hl = f_obs_cmpl.customized_copy(data = hl_data)
  return hl

def get_cc(f, hl):
  map_coeffs = abs(f).phase_transfer(phase_source = hl)
  return map_coeffs.map_correlation(other=f)

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
  n_real = m.data.focus()
  cs = crystal.symmetry(m.unit_cell_parameters, 1)
  crystal_gridding = maptbx.crystal_gridding(
    unit_cell         = cs.unit_cell(),
    space_group_info  = cs.space_group_info(),
    #symmetry_flags     = maptbx.use_space_group_symmetry,
    pre_determined_n_real = n_real)
  #
  d_min = params.d_min
  if(d_min is None and not params.box):
    d_min = maptbx.d_min_from_map(
      map_data  = m.data,
      unit_cell = cs.unit_cell())
  if(d_min is None):
    # box of reflections in |h|<N1/2, |k|<N2/2, 0<=|l|<N3/2
    max_index = [(i-1)//2 for i in n_real]
    print "max_index:", max_index
    complete_set = miller.build_set(
      crystal_symmetry = cs,
      anomalous_flag   = False,
      max_index        = max_index)
    indices = complete_set.indices()
    indices.append((0,0,0))
    complete_set = complete_set.customized_copy(indices = indices)
    #if(not params.box):
    #  # XXX What is sphere resolution corresponding to given box?
    #  uc = complete_set.unit_cell()
    #  d1 = uc.d([0,0,max_index[2]])
    #  d2 = uc.d([0,max_index[1],0])
    #  d3 = uc.d([max_index[0],1,0])
    #  print d1,d2,d3
    #  complete_set_sp = miller.build_set(
    #    crystal_symmetry = cs,
    #    anomalous_flag   = False,
    #    d_min            = min(d1,d2,d3))
    #  complete_set = complete_set.common_set(complete_set_sp)
  else:
    complete_set = miller.build_set(
      crystal_symmetry = cs,
      anomalous_flag   = False,
      d_min            = d_min)
  broadcast(m="Complete set information:", log=log)
  complete_set.show_comprehensive_summary(prefix="  ")
  try:
    f_obs_cmpl = complete_set.structure_factors_from_map(
      map            = m.data.as_double(),
      use_scale      = True,
      anomalous_flag = False,
      use_sg         = False)
  except Exception, e:
    if(str(e) == "cctbx Error: Miller index not in structure factor map."):
      msg = "Too high resolution requested. Try running with larger d_min."
      raise Sorry(msg)
    else:
      raise Sorry(str(e))
  mtz_dataset = f_obs_cmpl.as_mtz_dataset(column_root_label="F")
  mtz_dataset.add_miller_array(
    miller_array      = abs(f_obs_cmpl),
    column_root_label = "F_ampl")
  mtz_dataset.add_miller_array(
    miller_array      = f_obs_cmpl.customized_copy(data=flex.double(f_obs_cmpl.data().size(),1)),
    column_root_label = "SIGF_ampl")
  # convert phases into HL coefficeints
  broadcast(m="Convert phases into HL coefficeints:", log=log)
  hl = get_hl(f_obs_cmpl=f_obs_cmpl, k_blur=params.k_blur, b_blur=params.b_blur)
  cc = get_cc(f = f_obs_cmpl, hl = hl)
  print "cc:", cc
  if(abs(1.-cc)>1.e-3):
    print "Supplied b_blur is not good. Attempting to find optimal b_blur."
    cc_best = 999.
    b_blur_best = params.b_blur
    for b_blur in range(1, 100):
      hl = get_hl(f_obs_cmpl=f_obs_cmpl, k_blur=params.k_blur, b_blur=b_blur)
      cc = get_cc(f = f_obs_cmpl, hl = hl)
      if(cc<cc_best):
        cc_best = cc
        b_blur_best = b_blur
      if(abs(1.-cc)<1.e-3):
        b_blur_best = b_blur
        break
    hl = get_hl(f_obs_cmpl=f_obs_cmpl, k_blur=params.k_blur, b_blur=b_blur_best)
    print "cc:", get_cc(f = f_obs_cmpl, hl = hl)
    print "b_blur_best:", b_blur_best
  mtz_dataset.add_miller_array(
    miller_array      = hl,
    column_root_label = "HL")
  # write output MTZ file with all the data
  broadcast(m="Writing output MTZ file:", log=log)
  print >> log, "  file name:", params.output_file_name
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = params.output_file_name)

if(__name__ == "__main__"):
  run(sys.argv[1:])
  print "All done."
