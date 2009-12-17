from iotbx import shelx
from smtbx import masks
from cctbx import miller
from cctbx import xray
from cctbx.array_family import flex
from iotbx import reflection_file_utils
from iotbx import reflection_file_reader

from libtbx.test_utils import approx_equal

from cStringIO import StringIO

def exercise_masks(xs, fo_sq,
                   solvent_radius=1.3,
                   shrink_truncation_radius=1.3,
                   debug=False):
  fo_sq_merged = fo_sq.merge_equivalents().array()
  mask = masks.mask(xs, fo_sq_merged)
  print "solvent_radius: %f" %solvent_radius
  print "shrink_truncation_radius: %f" %shrink_truncation_radius
  mask.compute(solvent_radius=solvent_radius,
               shrink_truncation_radius=shrink_truncation_radius)
  f_mask = mask.structure_factors(max_cycles=10)
  f_model = mask.f_model()
  # write modified structure factors as shelxl hkl
  out = StringIO()
  f_mask
  modified_fo_sq = mask.modified_structure_factors()
  modified_fo_sq.export_as_shelx_hklf(out)
  out_file = open('modified.hkl', 'w')
  out_file.write(out.getvalue())
  out_file.close()
  #
  if debug:
    f_obs = fo_sq_merged.as_amplitude_array()
    sf = xray.structure_factors.from_scatterers(
      miller_set=f_obs,
      cos_sin_table=True)
    f_calc = sf(xs, f_obs).f_calc()
    f_model = mask.f_model()
    scale_factor = f_obs.quick_scale_factor_approximation(
      f_model, cutoff_factor=0.8)
    # f_obs - f_calc
    k = f_obs.quick_scale_factor_approximation(f_calc, cutoff_factor=0.8)
    f_obs_minus_f_calc = f_obs.f_obs_minus_f_calc(1./k, f_calc)
    diff_map_calc = miller.fft_map(mask.crystal_gridding, f_obs_minus_f_calc)
    diff_map_calc.apply_volume_scaling()
    # f_obs - f_model
    f_obs_minus_f_model = f_obs.f_obs_minus_f_calc(1./scale_factor, f_model)
    diff_map_model = miller.fft_map(mask.crystal_gridding, f_obs_minus_f_model)
    diff_map_model.apply_volume_scaling()
    # f_model
    model_map = miller.fft_map(mask.crystal_gridding, f_model)
    model_map.apply_volume_scaling()
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
      title="f_obs - f_model",
      raw_map=diff_map_model.real_map(),
      unit_cell=f_obs.unit_cell())
    wx_map_viewer.display(
      title="f_model",
      raw_map=model_map.real_map(),
      unit_cell=f_obs.unit_cell())
  return

def run(ins_name, hkl_name,
        solvent_radius=1.4,
        shrink_truncation_radius=1.4,
        debug=False):
  xs = xray.structure.from_shelx(filename=ins_name)
  reflections_server = reflection_file_utils.reflection_file_server(
    crystal_symmetry = xs.crystal_symmetry(),
    reflection_files = [
      reflection_file_reader.any_reflection_file('hklf4=%s' %hkl_name)
    ]
  )
  fo_sq = reflections_server.get_miller_arrays(None)[0]
  exercise_masks(xs, fo_sq,
                 solvent_radius=float(solvent_radius),
                 shrink_truncation_radius=float(shrink_truncation_radius),
                 debug=debug)
  print "OK"

if __name__ == '__main__':
  import sys
  run(*(sys.argv[1:]))
