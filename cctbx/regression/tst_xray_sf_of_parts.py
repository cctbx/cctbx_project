import random
import time
from cctbx.development import random_structure
from cctbx import sgtbx
from cctbx import xray
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal

def exercise(space_group_symbol, anomalous_flag, b_shift):
  """
  With sampling fixed by the whole structure, FFT structure factors of its
  parts add up to those of the whole, and gradients calculated for a part
  are equal to those calculated for the whole. With large b_shift u_extra
  is negative.
  """
  xrs = random_structure.xray_structure(
    space_group_info = sgtbx.space_group_info(space_group_symbol),
    elements         = ["C","N","O","S"]*10,
    volume_per_atom  = 25,
    min_distance     = 1.2,
    use_u_iso        = True,
    use_u_aniso      = True,
    random_u_iso     = True,
    random_f_double_prime = anomalous_flag,
    random_occupancy = True)
  xrs.shift_us(b_shift = b_shift)
  n = xrs.scatterers().size()
  i_sel = flex.size_t(sorted(random.sample(range(n), 12)))
  b_sel = flex.bool(n, i_sel)
  xrs.scatterers().flags_set_grads(state = False)
  xrs.scatterers().flags_set_grad_occupancy(iselection = i_sel)
  xrs_sel = xrs.select(i_sel)
  xrs_rest = xrs.select(~b_sel)
  miller_set = xrs.build_miller_set(
    anomalous_flag = anomalous_flag, d_min = 1.5)
  d_target_d_f_calc = flex.polar(
    flex.random_double(miller_set.size()),
    flex.random_double(miller_set.size())*6)
  #
  def f_calc(manager, xray_structure):
    return manager(
      xray_structure = xray_structure,
      miller_set     = miller_set,
      algorithm      = "fft").f_calc().data()
  def gradients(manager, xray_structure):
    return manager(
      xray_structure         = xray_structure,
      u_iso_refinable_params = None,
      miller_set             = miller_set,
      d_target_d_f_calc      = d_target_d_f_calc,
      n_parameters           = 0,
      algorithm              = "fft").d_target_d_occupancy()
  # default sampling depends on scatterers passed
  from_scatterers = xray.structure_factors.from_scatterers(
    miller_set = miller_set)
  f_all = f_calc(from_scatterers, xrs)
  scale = flex.mean(flex.abs(f_all))
  f_parts = f_calc(from_scatterers, xrs_sel) \
          + f_calc(from_scatterers, xrs_rest)
  assert flex.max(flex.abs(f_all - f_parts))/scale > 1.e-4
  gradients_manager = xray.structure_factors.gradients(miller_set = miller_set)
  g_all = gradients(gradients_manager, xrs).select(i_sel)
  # sampling fixed by the whole structure
  from_scatterers.fix_sampling(xray_structure = xrs)
  assert approx_equal(f_calc(from_scatterers, xrs), f_all, 1.e-12*scale)
  f_parts = f_calc(from_scatterers, xrs_sel) \
          + f_calc(from_scatterers, xrs_rest)
  assert flex.max(flex.abs(f_all - f_parts))/scale < 1.e-10
  gradients_manager.fix_sampling(xray_structure = xrs)
  g_sel = gradients(gradients_manager, xrs_sel)
  assert flex.max(flex.abs(g_all)) > 1.e-3
  assert approx_equal(g_sel, g_all, 1.e-10*flex.max(flex.abs(g_all)))

def run():
  random.seed(0)
  flex.set_random_seed(0)
  for space_group_symbol in ["P1", "P21", "P212121"]:
    for anomalous_flag in [False, True]:
      for b_shift in [0, 50]:
        exercise(space_group_symbol, anomalous_flag, b_shift)

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("Time: %6.2f"%(time.time()-t0))
  print("OK")
