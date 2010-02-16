from __future__ import division

from cctbx import crystal, miller, uctbx, xray
from cctbx.array_family import flex
from iotbx import reflection_file_utils, reflection_file_reader
from iotbx import shelx
from smtbx import masks
from cctbx.masks import flood_fill
from libtbx.test_utils import approx_equal
from libtbx.utils import time_log
from iotbx.option_parser import option_parser

from cStringIO import StringIO

def exercise_masks(xs, fo_sq,
                   solvent_radius,
                   shrink_truncation_radius,
                   resolution_factor,
                   use_space_group_symmetry,
                   debug=False,
                   timing=False):
  xs_ref = xs.deep_copy_scatterers()
  time_total = time_log("masks total").start()
  fo_sq_merged = fo_sq.merge_equivalents().array()
  mask = masks.mask(xs, fo_sq_merged)
  time_compute_mask = time_log("compute mask").start()
  mask.compute(solvent_radius=solvent_radius,
               shrink_truncation_radius=shrink_truncation_radius,
               resolution_factor=resolution_factor,
               atom_radii_table={'C':1.70, 'B':1.63, 'N':1.55, 'O':1.52},
               use_space_group_symmetry=use_space_group_symmetry)
  time_compute_mask.stop()
  time_structure_factors = time_log("structure factors").start()
  f_mask = mask.structure_factors()
  time_structure_factors.stop()
  mask.show_summary()
  print "F000 void: %.1f" %mask.f_000_s
  f_model = mask.f_model()
  # write modified structure factors as shelxl hkl
  out = StringIO()
  modified_fo_sq = mask.modified_structure_factors()
  modified_fo_sq.export_as_shelx_hklf(out)
  out_file = open('modified.hkl', 'w')
  out_file.write(out.getvalue())
  out_file.close()

  if timing:
    print
    print time_log.legend
    print time_compute_mask.report()
    print time_structure_factors.report()
    print time_total.log()

  if debug:
    f_obs = fo_sq_merged.as_amplitude_array()
    sf = xray.structure_factors.from_scatterers(
      miller_set=f_obs,
      cos_sin_table=True)
    f_calc = sf(xs, f_obs).f_calc()
    f_model = mask.f_model()
    scale_factor = f_obs.quick_scale_factor_approximation(
      f_model, cutoff_factor=0)
    # f_obs - f_calc
    k = f_obs.quick_scale_factor_approximation(f_calc, cutoff_factor=0)
    f_obs_minus_f_calc = f_obs.f_obs_minus_f_calc(1/k, f_calc)
    diff_map_calc = miller.fft_map(mask.crystal_gridding, f_obs_minus_f_calc)
    diff_map_calc.apply_volume_scaling()
    # f_mask
    mask_map = miller.fft_map(mask.crystal_gridding, f_mask)
    mask_map.apply_volume_scaling()
    # f_model
    model_map = miller.fft_map(mask.crystal_gridding, f_model)
    model_map.apply_volume_scaling()
    # f_obs - f_model
    f_obs_minus_f_model = f_obs.f_obs_minus_f_calc(1/scale_factor, f_model)
    diff_map_model = miller.fft_map(mask.crystal_gridding, f_obs_minus_f_model)
    diff_map_model.apply_volume_scaling()
    # modified f_obs
    modified_fo_sq_map = miller.fft_map(
      mask.crystal_gridding, modified_fo_sq.as_amplitude_array().phase_transfer(f_calc))
    modified_fo_sq_map.apply_volume_scaling()
    # view the maps
    from crys3d import wx_map_viewer
    wx_map_viewer.display(
      title="Mask",
      raw_map=mask.mask.data.as_double(),
      unit_cell=f_obs.unit_cell())
    wx_map_viewer.display(
      title="f_obs - f_calc",
      raw_map=diff_map_calc.real_map(),
      unit_cell=f_obs.unit_cell())
    wx_map_viewer.display(
      title="f_mask",
      raw_map=mask_map.real_map(),
      unit_cell=f_obs.unit_cell())
    wx_map_viewer.display(
      title="f_model",
      raw_map=model_map.real_map(),
      unit_cell=f_obs.unit_cell())
    wx_map_viewer.display(
      title="f_obs - f_model",
      raw_map=diff_map_model.real_map(),
      unit_cell=f_obs.unit_cell())
    wx_map_viewer.display(
      title="modified_fo_sq",
      raw_map=modified_fo_sq_map.real_map(),
      unit_cell=f_obs.unit_cell())
  return

def run(args):
  command_line = (option_parser().enable_symmetry_comprehensive()
                  .option(None, "--structure",
                          action="store")
                  .option(None, "--solvent_radius",
                          action="store",
                          type="float",
                          default=1.3)
                  .option(None, "--shrink_truncation_radius",
                          action="store",
                          type="float",
                          default=1.3)
                  .option(None, "--debug",
                          action="store_true")
                  .option(None, "--timing",
                          action="store_true")
                  .option(None, "--resolution_factor",
                          action="store",
                          type="float",
                          default=1/3)
                  .option(None, "--use_space_group_symmetry",
                          action="store_true")).process(args=args)
  xs = xray.structure.from_shelx(filename=command_line.options.structure)
  reflections_server = reflection_file_utils.reflection_file_server(
    crystal_symmetry = xs.crystal_symmetry(),
    reflection_files = [
      reflection_file_reader.any_reflection_file(command_line.args[0])
    ]
  )
  fo_sq = reflections_server.get_miller_arrays(None)[0]

  print "structure file: %s" %command_line.options.structure
  print "reflection file: %s" %command_line.args[0]
  if command_line.options.debug:
    print "debug: %s" %command_line.options.debug
  print

  xs.show_summary()
  print

  exercise_masks(
    xs, fo_sq,
    solvent_radius=command_line.options.solvent_radius,
    shrink_truncation_radius=command_line.options.shrink_truncation_radius,
    resolution_factor=command_line.options.resolution_factor,
    use_space_group_symmetry=command_line.options.use_space_group_symmetry,
    debug=command_line.options.debug,
    timing=command_line.options.timing)
  print "OK"

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
