import random
import time
from cctbx.array_family import flex
from cctbx.development import random_structure
from cctbx import sgtbx

def finite_difference_slopes(d_min, algorithm, steps,
      exp_table_one_over_step_size=None):
  """
  Slopes of a least-squares target along a fixed direction of atom shifts,
  from central finite differences with given steps (A, norm of the shift of
  all atoms). Sharp atoms (small B) at high resolution.
  """
  random.seed(0)
  flex.set_random_seed(0)
  xrs = random_structure.xray_structure(
    space_group_info       = sgtbx.space_group_info("P212121"),
    elements               = ["C", "N", "O"]*10,
    volume_per_atom        = 30,
    min_distance           = 1.2,
    general_positions_only = True,
    random_u_iso           = True,
    u_iso                  = 0.03)
  f_obs = abs(xrs.structure_factors(
    d_min = d_min, algorithm = "direct").f_calc())
  # Start away from the minimum
  xrs.shake_sites_in_place(mean_distance = 0.05)
  sites_cart = xrs.sites_cart()
  direction = flex.vec3_double(flex.random_double(sites_cart.size()*3)-0.5)
  direction = direction * (1./direction.as_double().norm())
  def target(step):
    xrs_shifted = xrs.deep_copy_scatterers()
    xrs_shifted.set_sites_cart(sites_cart + direction*step)
    f_calc = abs(f_obs.structure_factors_from_scatterers(
      xray_structure               = xrs_shifted,
      algorithm                    = algorithm,
      exp_table_one_over_step_size = exp_table_one_over_step_size).f_calc())
    return flex.sum(flex.pow2(f_obs.data() - f_calc.data()))
  return [(target(step)-target(-step))/(2*step) for step in steps]

def exercise():
  """
  Structure factors computed with FFT are a continuous function of atomic
  coordinates: finite difference slopes of a target do not depend on the
  step. This is what a line search of a minimizer relies on. It takes
  interpolation in the table of exponents used to sample atomic density
  (the default, exp_table_one_over_step_size=-100) or no table (0).
  """
  steps = [1.e-4, 1.e-5, 1.e-6]
  for algorithm, exp_table_one_over_step_size in [
      ("direct", None), ("fft", None), ("fft", -100), ("fft", 0)]:
    slopes = finite_difference_slopes(
      d_min                        = 0.7,
      algorithm                    = algorithm,
      steps                        = steps,
      exp_table_one_over_step_size = exp_table_one_over_step_size)
    for slope in slopes:
      assert abs(slope-slopes[0]) < 0.01*abs(slopes[0]), (
        algorithm, exp_table_one_over_step_size, slopes)

if (__name__ == "__main__"):
  t0 = time.time()
  exercise()
  print("OK time =%8.3f"%(time.time() - t0))
