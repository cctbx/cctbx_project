from __future__ import absolute_import, division, print_function
from cctbx import xray
from cctbx import adptbx
from cctbx import miller
from cctbx import crystal
from cctbx.development import debug_utils
from cctbx.array_family import flex
from libtbx.utils import user_plus_sys_time
import sys
from six.moves import range

def dummy_structure(space_group_info, volume, n_scatterers):
  structure = xray.structure(
    special_position_settings=crystal.special_position_settings(
      crystal_symmetry=crystal.symmetry(
        unit_cell=space_group_info.any_compatible_unit_cell(volume=volume),
        space_group_info=space_group_info)))
  b_iso = 20
  u_iso = adptbx.b_as_u(b_iso)
  u_star = adptbx.u_iso_as_u_star(structure.unit_cell(), u_iso)
  scatterer = xray.scatterer(label="C", site=(0.123,0.234,0.345), u=u_star)
  for i in range(n_scatterers):
    structure.add_scatterer(scatterer)
  return structure

def dummy_miller_set(crystal_symmetry, log_n_reflections):
  indices = flex.miller_index(((1,2,3),))
  for i in range(log_n_reflections):
    indices.append(indices)
  return miller.set(
    crystal_symmetry=crystal_symmetry,
    indices=indices,
    anomalous_flag=False)

def get_times(space_group_info, volume=100000, d_min=1.5):
  crystal_symmetry = dummy_structure(space_group_info, volume, 0)
  crystal_symmetry.show_summary()
  structure_factors_from_scatterers = xray.structure_factors.from_scatterers(
    crystal_symmetry=crystal_symmetry,
    d_min=d_min)
  for log_n_reflections in range(15,16):
    miller_set = dummy_miller_set(crystal_symmetry, log_n_reflections)
    for log_n_scatterers in range(1,10):
      structure = dummy_structure(space_group_info,volume,2**log_n_scatterers)
      timer = user_plus_sys_time()
      f_calc = structure_factors_from_scatterers(
        xray_structure=structure,
        miller_set=miller_set)
      print("%6d %6d %10.3f direct=%10.3f fft=%10.3f" % (
        2**log_n_scatterers, 2**log_n_reflections, timer.delta(),
        structure_factors_from_scatterers.estimate_time_direct(
          2**log_n_scatterers * 2**log_n_reflections),
        structure_factors_from_scatterers.estimate_time_fft(
          2**log_n_scatterers, 2**log_n_reflections)))
  e = structure_factors_from_scatterers.estimate_time_fft
  print("time_sampling:", e.time_sampling)
  print("time_fft:", e.time_fft)
  print("time_from_or_to_map:", e.time_from_or_to_map)
  print("time_apply_u_extra:", e.time_apply_u_extra)

def run_call_back(flags, space_group_info):
  get_times(space_group_info)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

if (__name__ == "__main__"):
  run()
