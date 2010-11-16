# LIBTBX_SET_DISPATCHER_NAME smtbx.masks

from __future__ import division

from cctbx import miller, sgtbx, uctbx, xray
from iotbx import reflection_file_utils, reflection_file_reader
from iotbx import shelx
from smtbx import masks
from libtbx.utils import time_log
from iotbx.option_parser import option_parser

from cStringIO import StringIO
import os

def exercise_masks(xs, fo_sq,
                   solvent_radius,
                   shrink_truncation_radius,
                   resolution_factor=None,
                   grid_step=None,
                   resolution_cutoff=None,
                   atom_radii_table=None,
                   use_space_group_symmetry=False,
                   debug=False,
                   verbose=False):
  assert resolution_factor is None or grid_step is None
  xs_ref = xs.deep_copy_scatterers()
  time_total = time_log("masks total").start()
  fo_sq = fo_sq.customized_copy(anomalous_flag=True)
  fo_sq = fo_sq.eliminate_sys_absent()
  merging = fo_sq.merge_equivalents()
  fo_sq_merged = merging.array()
  if resolution_cutoff is not None:
    fo_sq_merged = fo_sq_merged.resolution_filter(d_min=resolution_cutoff)
  if verbose:
    print "Merging summary:"
    print "R-int, R-sigma: %.4f, %.4f" %(merging.r_int(), merging.r_sigma())
    merging.show_summary()
    print
    fo_sq_merged.show_comprehensive_summary()
    print
  mask = masks.mask(xs, fo_sq_merged)
  time_compute_mask = time_log("compute mask").start()
  mask.compute(solvent_radius=solvent_radius,
               shrink_truncation_radius=shrink_truncation_radius,
               resolution_factor=resolution_factor,
               grid_step=grid_step,
               atom_radii_table=atom_radii_table,
               use_space_group_symmetry=use_space_group_symmetry)
  time_compute_mask.stop()
  time_structure_factors = time_log("structure factors").start()
  f_mask = mask.structure_factors()
  time_structure_factors.stop()
  mask.show_summary()
  f_model = mask.f_model()
  # write modified intensities as shelxl hkl
  out = StringIO()
  modified_fo_sq = mask.modified_intensities()
  modified_fo_sq.export_as_shelx_hklf(out)
  out_file = open('modified.hkl', 'w')
  out_file.write(out.getvalue())
  out_file.close()

  if verbose:
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
    scale_factor = f_obs.scale_factor(f_model)
    # f_obs - f_calc
    k = f_obs.scale_factor(f_calc)
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
  return mask

def run(args):
  def vdw_radii_callback(option, opt_str, value, parser):
    # create a dict from space separated string of element types and radii
    radii = {}
    items = value.split()
    assert len(items) % 2 == 0
    for i in range(int(len(items) / 2)):
      radii.setdefault(items[i*2], float(items[i*2+1]))
    setattr(parser.values, option.dest, radii)
  command_line = (option_parser(
    usage="smtbx.masks structure reflections [options]")
                  .enable_symmetry_comprehensive()
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
                  .option(None, "--verbose",
                          action="store_true")
                  .option(None, "--resolution_factor",
                          action="store",
                          type="float",
                          default=1/4)
                  .option(None, "--grid_step",
                          action="store",
                          type="float")
                  .option(None, "--d_min",
                          action="store",
                          type="float")
                  .option(None, "--two_theta_max",
                          action="store",
                          type="float")
                  .option(None, "--cb_op",
                          action="store",
                          type="string")
                  .option(None, "--vdw_radii",
                          action="callback",
                          callback=vdw_radii_callback,
                          type="string",
                          nargs=1)
                  .option(None, "--use_space_group_symmetry",
                          action="store_true")).process(args=args)
  structure_file = command_line.args[0]
  ext = os.path.splitext(structure_file)[-1].lower()
  if ext in ('.res', '.ins'):
    xs = xray.structure.from_shelx(filename=structure_file)
  elif ext == '.cif':
    xs = xray.structure.from_cif(filename=structure_file)
  else:
    print "%s: unsupported structure file format {shelx|cif}" %ext
    return
  reflections_server = reflection_file_utils.reflection_file_server(
    crystal_symmetry = xs.crystal_symmetry(),
    reflection_files = [
      reflection_file_reader.any_reflection_file(command_line.args[1])
    ]
  )
  fo_sq = reflections_server.get_miller_arrays(None)[0]

  if command_line.options.cb_op is not None:
    cb_op = sgtbx.change_of_basis_op(sgtbx.rt_mx(command_line.options.cb_op))
    fo_sq = fo_sq.change_basis(cb_op).customized_copy(
      crystal_symmetry=xs)

  print "structure file: %s" %command_line.args[0]
  print "reflection file: %s" %command_line.args[1]
  if command_line.options.debug:
    print "debug: %s" %command_line.options.debug
  print

  xs.show_summary()
  print

  d_min = command_line.options.d_min
  two_theta_max = command_line.options.two_theta_max
  assert [d_min, two_theta_max].count(None) > 0
  if two_theta_max is not None:
    d_min = uctbx.two_theta_as_d(two_theta_max, wavelength=0.71073, deg=True)
  exercise_masks(
    xs, fo_sq,
    solvent_radius=command_line.options.solvent_radius,
    shrink_truncation_radius=command_line.options.shrink_truncation_radius,
    resolution_factor=command_line.options.resolution_factor,
    grid_step=command_line.options.grid_step,
    resolution_cutoff=d_min,
    atom_radii_table=command_line.options.vdw_radii,
    use_space_group_symmetry=command_line.options.use_space_group_symmetry,
    debug=command_line.options.debug,
    verbose=command_line.options.verbose)
  print "OK"

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
