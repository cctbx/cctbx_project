from __future__ import division
from __future__ import print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.map_to_structure_factors

from builtins import str
from builtins import range
import iotbx.ccp4_map
from cctbx import miller, crystal
from cctbx.array_family import flex
import mmtbx.utils
import sys
from libtbx.utils import Sorry
from cctbx import maptbx
from cctbx import miller

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
  print("-"*79, file=log)
  print(m, file=log)
  print("*"*len(m), file=log)

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

def run(args, log=None, ccp4_map=None,
    return_as_miller_arrays=False, nohl=False,
    space_group_number=None,
    out=sys.stdout):
  if log is None: log=out
  inputs = mmtbx.utils.process_command_line_args(args = args,
    master_params = master_params())
  got_map = False
  if ccp4_map: got_map=True
  broadcast(m="Parameters:", log=log)
  inputs.params.show(prefix="  ",out=out)
  params = inputs.params.extract()
  if(ccp4_map is None and inputs.ccp4_map is not None):
    broadcast(m="Processing input CCP4 map file: %s"%inputs.ccp4_map_file_name,
      log=log)
    ccp4_map = inputs.ccp4_map
    ccp4_map.show_summary(prefix="  ",out=out)
    got_map = True
  if(not got_map):
    raise Sorry("Map file is needed.")
  #
  m = ccp4_map
  if(m.space_group_number > 1):
    raise Sorry("Input map space group: %d. Must be P1."%m.space_group_number)
  broadcast(m="Input map information:", log=log)
  print("m.all()   :", m.data.all(), file=out)
  print("m.focus() :", m.data.focus(), file=out)
  print("m.origin():", m.data.origin(), file=out)
  print("m.nd()    :", m.data.nd(), file=out)
  print("m.size()  :", m.data.size(), file=out)
  print("m.focus_size_1d():", m.data.focus_size_1d(), file=out)
  print("m.is_0_based()   :", m.data.is_0_based(), file=out)
  print("map: min/max/mean:", flex.min(m.data), flex.max(m.data), flex.mean(m.data), file=out)
  print("unit cell:", m.unit_cell_parameters, file=out)
  #
  map_data=m.data
  map_data = maptbx.shift_origin_if_needed(map_data = map_data).map_data
  # generate complete set of Miller indices up to given high resolution d_min
  n_real = map_data.focus()
  if not space_group_number:
    space_group_number=1
  if space_group_number <=1:
     symmetry_flags=None
  else:
    symmetry_flags = maptbx.use_space_group_symmetry,

  cs = crystal.symmetry(m.unit_cell_parameters, space_group_number)
  crystal_gridding = maptbx.crystal_gridding(
    unit_cell         = cs.unit_cell(),
    space_group_info  = cs.space_group_info(),
    symmetry_flags     = symmetry_flags,
    pre_determined_n_real = n_real)
  #
  d_min = params.d_min
  if(d_min is None and not params.box):
    d_min = maptbx.d_min_from_map(
      map_data  = map_data,
      unit_cell = cs.unit_cell())
  if(d_min is None):
    # box of reflections in |h|<N1/2, |k|<N2/2, 0<=|l|<N3/2
    f_obs_cmpl = miller.structure_factor_box_from_map(
      map              = map_data.as_double(),
      crystal_symmetry = cs,
      include_000      = True)
  else:
    complete_set = miller.build_set(
      crystal_symmetry = cs,
      anomalous_flag   = False,
      d_min            = d_min)
    try:
      f_obs_cmpl = complete_set.structure_factors_from_map(
        map            = map_data.as_double(),
        use_scale      = True,
        anomalous_flag = False,
        use_sg         = False)
    except Exception as e:
      if(str(e) == "cctbx Error: Miller index not in structure factor map."):
        msg = "Too high resolution requested. Try running with larger d_min."
        raise Sorry(msg)
      else:
        raise Sorry(str(e))

  if nohl and return_as_miller_arrays:
    return f_obs_cmpl

  mtz_dataset = f_obs_cmpl.as_mtz_dataset(column_root_label="F")
  f_obs = abs(f_obs_cmpl)
  f_obs.set_sigmas(sigmas=flex.double(f_obs_cmpl.data().size(),1))
  mtz_dataset.add_miller_array(
    miller_array      = f_obs,
    column_root_label = "F-obs")
  mtz_dataset.add_miller_array(
    miller_array      = f_obs.generate_r_free_flags(),
    column_root_label = "R-free-flags")
  if not nohl:
    # convert phases into HL coefficeints
    broadcast(m="Convert phases into HL coefficients:", log=log)
    hl = get_hl(f_obs_cmpl=f_obs_cmpl, k_blur=params.k_blur, b_blur=params.b_blur)
    cc = get_cc(f = f_obs_cmpl, hl = hl)
    print("cc:", cc, file=out)
    if(abs(1.-cc)>1.e-3):
      print("Supplied b_blur is not good. Attempting to find optimal b_blur.", file=out)
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
      print("cc:", get_cc(f = f_obs_cmpl, hl = hl), file=out)
      print("b_blur_best:", b_blur_best, file=out)
    mtz_dataset.add_miller_array(
      miller_array      = hl,
      column_root_label = "HL")
  else:
    hl=None
  if return_as_miller_arrays:
    return f_obs_cmpl,hl
  else:
    # write output MTZ file with all the data
    broadcast(m="Writing output MTZ file:", log=log)
    print("  file name:", params.output_file_name, file=log)
    mtz_object = mtz_dataset.mtz_object()
    mtz_object.write(file_name = params.output_file_name)

if(__name__ == "__main__"):
  run(sys.argv[1:])
  print("All done.")
