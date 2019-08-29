from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
import time, random
from mmtbx_rsr_ext import *
from libtbx.utils import user_plus_sys_time
from cctbx.development import random_structure
from cctbx import sgtbx, adptbx

if (1):
  random.seed(0)
  flex.set_random_seed(0)


def exercise_01(grid_step = 0.03, d_min = 1.0, wing_cutoff = 1.e-9):
  xrs = random_structure.xray_structure(
    space_group_info       = sgtbx.space_group_info("P 1"),
    elements               = ["O","N","C","P","S","U","AU"]*1,
    random_u_iso           = True,
    general_positions_only = False)
  # avoid excessive_range_error_limit crash
  bs = xrs.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1)
  sel = bs < 1
  bs = bs.set_selected(sel, 1)
  xrs.set_b_iso(values = bs)
  #
  p = xrs.unit_cell().parameters()
  timer = user_plus_sys_time()
  res = manager(nx = int(p[0]/grid_step),
                ny = int(p[1]/grid_step),
                nz = int(p[2]/grid_step),
                scattering_type_registry = xrs.scattering_type_registry(),
                unit_cell = xrs.unit_cell(),
                scatterers = xrs.scatterers(),
                wing_cutoff = wing_cutoff)
  print("time: %10.4f" % (timer.elapsed()))
  f_calc_dir = xrs.structure_factors(
    d_min     = d_min,
    algorithm = "direct").f_calc()
  #
  f_calc_den = f_calc_dir.structure_factors_from_map(map = res.density_array,
    use_scale = True)
  f1 = flex.abs(f_calc_dir.data())
  f2 = flex.abs(f_calc_den.data())
  r = flex.sum(flex.abs(f1-f2))/flex.sum(f2)
  print("r-factor:", r)
  assert r < 1.e-4, r


if (__name__=="__main__"):
  exercise_01()
