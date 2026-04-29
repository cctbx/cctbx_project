"""Calculate structure factors and HL coefficients from ccp4 map file and save in MTZ"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.map_to_structure_factors

import iotbx.ccp4_map
from cctbx.array_family import flex
import mmtbx.utils
import sys
from libtbx.utils import Sorry
from six.moves import range

master_params_str = """
output_file_name = map_to_structure_factors.mtz
  .type=str
d_min = None
  .type=float
  .help = Resolution of output structure factors. Default is based on the\
           gridding of the map and can lead to map coefficients that are \
           at much higher resolution than the map.
  .short_caption = Resolution
resolution_factor = 1./3
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

def get_cc(f, hl):
  map_coeffs = abs(f).phase_transfer(phase_source = hl)
  return map_coeffs.map_correlation(other=f)

def get_shift_cart(map_data, crystal_symmetry, origin=None):
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
  if(m.unit_cell_crystal_symmetry().space_group_number()> 1):
    raise Sorry("Input map space group: %d. Must be P1."%m.unit_cell_crystal_symmetry().space_group_number())
  broadcast(m="Input map information:", log=log)
  print("m.all()   :", m.map_data().all(), file=out)
  print("m.focus() :", m.map_data().focus(), file=out)
  print("m.origin():", m.map_data().origin(), file=out)
  print("m.nd()    :", m.map_data().nd(), file=out)
  print("m.size()  :", m.map_data().size(), file=out)
  print("m.focus_size_1d():", m.map_data().focus_size_1d(), file=out)
  print("m.is_0_based()   :", m.map_data().is_0_based(), file=out)
  print("map: min/max/mean:", flex.min(m.map_data()), flex.max(m.map_data()), flex.mean(m.map_data()), file=out)
  print("unit cell:", m.unit_cell().parameters(), file=out)
  #
  if m.unit_cell_grid == m.map_data().all():
    print("\nOne unit cell of data is present in map", file=out)
  else:
    if params.keep_origin:
      print("\nNOTE: This map does not have exactly one unit cell of data, so \n"+\
        "keep_origin is not available\n", file=out)
      print("--> Setting keep_origin=False and creating a new cell and gridding.\n", file=out)
      params.keep_origin=False
    print("Moving origin of input map to (0,0,0)", file=out)
    print("New cell will be: (%.3f, %.3f, %.3f, %.1f, %.1f, %.1f) A " %(
       m.crystal_symmetry().unit_cell().parameters()), file=out)
    print("New unit cell grid will be: (%s, %s, %s) "%(
      m.map_data().all()), file=out)
    m.unit_cell_grid = m.map_data().all()

  import iotbx.map_manager
  mm = iotbx.map_manager.map_manager(
    map_data                   = m.map_data().as_double(),
    unit_cell_grid             = m.unit_cell_grid,
    unit_cell_crystal_symmetry = m.crystal_symmetry(),
    wrapping                   = True)

  # Get origin in grid units and new position of origin in grid units
  original_origin=mm.map_data().origin()
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
    shift_cart=get_shift_cart(map_data=mm.map_data(), crystal_symmetry=mm.crystal_symmetry(),
      origin=new_origin)
  else:
    shift_cart=(0,0,0,)

  # Shift the map data if necessary
  mm.shift_origin()

  f_obs_cmpl = mm.map_as_fourier_coefficients(
     d_min             = params.d_min,
     box               = params.box,
     resolution_factor = params.resolution_factor)

  if params.scale_max is not None:
    f_obs_cmpl = f_obs_cmpl.apply_scaling(target_max=params.scale_max)

  from scitbx.matrix import col
  if col(shift_cart) != col((0,0,0,)):
    print("Output origin is at: (%.3f, %.3f, %.3f) A "%(
      tuple(-col(shift_cart))), file=out)
    f_obs_cmpl=f_obs_cmpl.translational_shift(
        mm.crystal_symmetry().unit_cell().fractionalize(-col(shift_cart)), deg=False)
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
  if not nohl and params.k_blur is not None and params.b_blur is not None:
    # convert phases into HL coefficeints
    broadcast(m="Convert phases into HL coefficients:", log=log)
    hl = f_obs_cmpl.make_up_hl_coeffs(k_blur=params.k_blur, b_blur=params.b_blur)
    cc = get_cc(f = f_obs_cmpl, hl = hl)
    print("cc:", cc, file=out)
    if(abs(1.-cc)>1.e-3):
      print("Supplied b_blur is not good. Attempting to find optimal b_blur.", file=out)
      cc_best = 999.
      b_blur_best = params.b_blur
      for b_blur in range(1, 100):
        hl = f_obs_cmpl.make_up_hl_coeffs(k_blur=params.k_blur, b_blur=b_blur)
        cc = get_cc(f = f_obs_cmpl, hl = hl)
        if(cc<cc_best):
          cc_best = cc
          b_blur_best = b_blur
        if(abs(1.-cc)<1.e-3):
          b_blur_best = b_blur
          break
      hl = f_obs_cmpl.make_up_hl_coeffs(k_blur=params.k_blur, b_blur=b_blur_best)
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

