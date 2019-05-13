from __future__ import division
from __future__ import print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.map_to_structure_factors

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
  .help = Resolution of output structure factors. Default is based on the\
           gridding of the map and can lead to map coefficients that are \
           at much higher resolution than the map.
  .short_caption = Resolution
resolution_factor = 0.333
  .type = float
  .help = Scale factor to guess resolution of output structure factors.\
          A scale factor of 0.5 gives the highest-resolution data allowed by \
          the map.  Usual is 0.33 or 0.25.
  .short_caption = Resolution factor
scale_max = 99999
  .type = float
  .help = Maximum value of output map coefficients amplitudes.  If None \
           use volume scaling
  .short_caption = Scale max
k_blur = 1
  .type = float
  .help = Scale applied to HL coefficients.  The HL coefficients are arbitrary\
          as no error information is available. The HL coefficients will have\
          values of k_blur at low resolution, falling off with an effective\
          B-value of b_blur at higher resolution.
  .short_caption = Scale on HL coefficients
b_blur  = 100
  .type = float
  .help = Blurring applied to HL coefficients.  The HL coefficients are \
          arbitrary\
          as no error information is available. The HL coefficients will have\
          values of k_blur at low resolution, falling off with an effective\
          B-value of b_blur at higher resolution.
  .short_caption = Blurring of HL coefficients
box = False
  .type = bool
  .help = You can choose to generate a full box of map coefficients based on\
          the gridding of the map.  Default is to generate map coefficients to\
          a specific resolution
  .short_caption = Box of Fourier coefficients
keep_origin = True
  .type = bool
  .help = Default (keep_origin=True) is to set the origin for the output map\
           coefficients to be the same as the input map. A map calculated from\
           the output map coefficients will superimpose on the input map. If \
           keep_origin=False then the new origin will be at (0,0,0). \
           Note: keep_origin=True is only available if the map covers an \
           entire unit cell. It will be automatically set to False if \
           less than a full unit cell is available.
  .short_caption = Keep origin
output_origin_grid_units = None
  .type = ints
  .help = You can set the origin of the output map (in grid units)
  .short_caption = Output origin (grid units)

# for wx GUI
map_file = None
  .type = path
  .help = Input map file
include scope libtbx.phil.interface.tracking_params
gui
  .help = "GUI-specific parameter required for output directory"
{
  output_dir = None
  .type = path
  .style = output_dir
}
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

def get_shift_cart(map_data=None,crystal_symmetry=None,
     origin=None):
  # this is the shift applied to the map when origin moves to (0,0,0)
  N = map_data.all()
  if origin is None:
    O = map_data.origin()
  else:
    O = tuple(origin)
  if(not crystal_symmetry.space_group().type().number() in [0,1]):
      raise Sorry("Space groups other than P1 are not supported.")
  a,b,c = crystal_symmetry.unit_cell().parameters()[:3]
  sx,sy,sz = O[0]/N[0],O[1]/N[1], O[2]/N[2]
  shift_frac = [-sx,-sy,-sz]
  shift_cart = crystal_symmetry.unit_cell().orthogonalize(shift_frac)
  return shift_cart

def run(args, log=None, ccp4_map=None,
    return_as_miller_arrays=False, nohl=False, return_f_obs=False,
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
  if not space_group_number:
    space_group_number=1
  if space_group_number <=1:
     symmetry_flags=None
  else:
    symmetry_flags = maptbx.use_space_group_symmetry,

  #cs = crystal.symmetry(m.unit_cell_parameters, space_group_number)
  # this will not work if m.unit_cell_grid != m.data.all()

  # Instead use ccp4 map crystal_symmetry and classify according to the case
  cs = m.crystal_symmetry()

  if m.unit_cell_grid == m.data.all():
    print("\nOne unit cell of data is present in map", file=out)
  else:
    if params.keep_origin:
      print("\nNOTE: This map does not have exactly one unit cell of data, so \n"+\
        "keep_origin is not available\n", file=out)
      print("--> Setting keep_origin=False and creating a new cell and gridding.\n", file=out)
      params.keep_origin=False
    print("Moving origin of input map to (0,0,0)", file=out)
    print("New cell will be: (%.3f, %.3f, %.3f, %.1f, %.1f, %.1f) A " %(
       cs.unit_cell().parameters()), file=out)
    print("New unit cell grid will be: (%s, %s, %s) "%(
      m.data.all()), file=out)

  map_data=m.data

  # Get origin in grid units and new position of origin in grid units
  original_origin=map_data.origin()
  print("\nInput map has origin at grid point (%s,%s,%s)" %(
        tuple(original_origin)), file=out)

  if params.output_origin_grid_units is not None:
    params.keep_origin=False
    new_origin=tuple(params.output_origin_grid_units)
    print("User-specified origin at grid point (%s,%s,%s)" %(
        tuple(params.output_origin_grid_units)), file=out)
    if tuple(params.output_origin_grid_units)==tuple(original_origin):
      print("This is the same as the input origin. No origin shift.", file=out)
  elif params.keep_origin:
    new_origin=original_origin
    print("Keeping origin at grid point  (%s,%s,%s)" %(
        tuple(original_origin)), file=out)
  else:
    new_origin=(0,0,0,)
    print("New origin at grid point (%s,%s,%s)" %(
        tuple((0,0,0,))), file=out)

  # shift_cart is shift away from (0,0,0)
  if new_origin != (0,0,0,):
    shift_cart=get_shift_cart(map_data=map_data,crystal_symmetry=cs,
      origin=new_origin)
  else:
    shift_cart=(0,0,0,)

  map_data = maptbx.shift_origin_if_needed(map_data = map_data).map_data
  # generate complete set of Miller indices up to given high resolution d_min
  n_real = map_data.focus()
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
      unit_cell = cs.unit_cell(),
      resolution_factor = params.resolution_factor)
    print("\nResolution of map coefficients using "+\
       "resolution_factor of %.2f: %.1f A\n" %(params.resolution_factor,d_min), file=out)
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
    except Exception, e:
      if(str(e) == "cctbx Error: Miller index not in structure factor map."):
        msg = "Too high resolution requested. Try running with larger d_min."
        raise Sorry(msg)
      else:
        raise Sorry(str(e))


  if params.scale_max is not None:
    f_obs_cmpl = f_obs_cmpl.apply_scaling(target_max=params.scale_max)

  from scitbx.matrix import col
  if col(shift_cart) != col((0,0,0,)):
    print("Output origin is at: (%.3f, %.3f, %.3f) A "%(
      tuple(-col(shift_cart))), file=out)
    f_obs_cmpl=f_obs_cmpl.translational_shift(
        cs.unit_cell().fractionalize(-col(shift_cart)), deg=False)
  else:
    print("Output origin is at (0.000, 0.000, 0.000) A", file=out)

  if nohl and return_as_miller_arrays and not return_f_obs:
    return f_obs_cmpl

  mtz_dataset = f_obs_cmpl.as_mtz_dataset(column_root_label="F")
  f_obs = abs(f_obs_cmpl)
  f_obs.set_sigmas(sigmas=flex.double(f_obs_cmpl.data().size(),1))
  if nohl and return_as_miller_arrays and return_f_obs:
     return f_obs
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
    if return_f_obs:
      return f_obs,hl
    else:
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
