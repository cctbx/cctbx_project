from __future__ import absolute_import, division, print_function

from cctbx.development import random_structure
from cctbx.sgtbx import space_group_info
import boost_adaptbx.boost.python as bp
from six.moves import range
asu_map_ext = bp.import_ext("cctbx_asymmetric_map_ext")
from cctbx_asymmetric_map_ext import *
from cctbx.array_family import flex
from cctbx import maptbx
import random
import iotbx.map_manager

if (1):
  random.seed(0)
  flex.set_random_seed(0)

def run():
  result = flex.bool()
  for sgn in range(1,231):
    group = space_group_info(sgn)
    xrs = random_structure.xray_structure(
      space_group_info       = group,
      volume_per_atom        = 25.,
      general_positions_only = False,
      elements               = ('C', 'N', 'O', 'H')*30,
      min_distance           = 1.0)
    sgt = xrs.space_group().type()
    fc = xrs.structure_factors(d_min=2).f_calc()
    fft_map = fc.fft_map(symmetry_flags = maptbx.use_space_group_symmetry,
      resolution_factor=1./3)
    map_data = fft_map.real_map_unpadded()
    n_real = map_data.all()
    mm = iotbx.map_manager.map_manager(
      map_data                   = map_data,
      unit_cell_grid             = map_data.accessor().all(),
      unit_cell_crystal_symmetry = xrs.crystal_symmetry(),
      origin_shift_grid_units    = [0,0,0],
      wrapping                   = True)
    assert mm.is_consistent_with_wrapping() in [True,None]
    assert mm.map_data().all() == n_real
    new_mm=mm.as_full_size_map()
    assert new_mm.map_data().all() == n_real

    # Now cut off edges and should not work:
    from iotbx.map_model_manager import map_model_manager
    mmm=map_model_manager(map_manager=mm)
    upper_bounds=tuple([x-1 for x in mm.map_data().all()])
    mmm.box_all_maps_with_bounds_and_shift_origin(lower_bounds=(1,1,1), upper_bounds=upper_bounds)
    mm=mmm.map_manager()
    assert mm.is_consistent_with_wrapping() in [None,False]
    assert mm.map_data().all() != n_real
    new_mm=mm.as_full_size_map()
    assert new_mm.map_data().all() == n_real

if (__name__ == "__main__"):
  run()
  print ("OK")
