from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx.development import random_structure
from cctbx import sgtbx
import mmtbx.f_model
from libtbx import group_args
from libtbx.test_utils import approx_equal
import random

def build_fmodel(n_atoms=50, d_min=1.8, seed=0, space_group="P 21 21 21"):
  random.seed(seed)
  flex.set_random_seed(seed)
  x = random_structure.xray_structure(
    space_group_info       = sgtbx.space_group_info(space_group),
    elements                = (("O", "N", "C") * n_atoms),
    volume_per_atom         = 200,
    min_distance            = 1.5,
    general_positions_only  = True,
    random_u_iso            = True)
  fc = x.structure_factors(d_min=d_min, algorithm="direct").f_calc()
  f_obs = abs(fc)
  r_free_flags = f_obs.generate_r_free_flags(fraction=0.1)
  fmodel = mmtbx.f_model.manager(
    xray_structure = x,
    f_obs          = f_obs,
    r_free_flags   = r_free_flags)
  fmodel.update_all_scales()
  return fmodel

def synthetic_llgi_data(fmodel, seed=1, feff_scale=1.0):
  """ feff deliberately rescaled/perturbed relative to f_obs (not merely
  a copy of it) so a test can tell whether a result actually depends on
  feff, rather than being coincidentally identical to the f_obs-based
  answer. """
  f_obs = fmodel.f_obs()
  n = f_obs.size()
  rnd = random.Random(seed)
  dobs = f_obs.array(data=flex.double([0.5 + 0.4 * rnd.random()
    for i in range(n)]))
  feff_data = f_obs.data() * feff_scale * flex.double(
    [0.85 + 0.3 * rnd.random() for i in range(n)])
  feff = f_obs.array(data=feff_data)
  teps = f_obs.array(data=flex.double(n, 1.0))
  resn = f_obs.array(data=flex.double(n, 1.0))
  return group_args(dobs=dobs, feff=feff, teps=teps, resn=resn, info=None)

def exercise_llgi_target_active_gating():
  # llgi_target_active() requires BOTH target_name=="llgi" AND
  # llgi_data attached -- neither alone is sufficient.
  fmodel = build_fmodel(seed=9)
  assert not fmodel.llgi_target_active()  # neither

  llgi_data = synthetic_llgi_data(fmodel, seed=9)
  fmodel.set_llgi_data(llgi_data)
  assert not fmodel.llgi_target_active()  # llgi_data but target=ml

  fmodel2 = build_fmodel(seed=9)
  fmodel2._target_name = "llgi"  # bypass set_target_name's own Sorry
  assert fmodel2.llgi_data() is None
  assert not fmodel2.llgi_target_active()  # target=llgi but no data

  fmodel.set_target_name("llgi")
  assert fmodel.llgi_target_active()  # both

def exercise_info_in_llgi_mode():
  # In llgi mode info() reports the ordinary r_work()/r_free()/r_all()
  # and takes its likelihood statistics (FOM, phase error, coordinate
  # error, D and V in place of alpha and beta) from the LLGI fit; the ML
  # alpha/beta machinery must not run at all. info() builds a real
  # target_functor() internally, which requires sigmaa/scatfrac to be
  # attached, so run the real per-macrocycle estimator first.
  from six.moves import cStringIO as StringIO
  fmodel = build_fmodel(n_atoms=70, d_min=1.7, seed=18)
  llgi_data = synthetic_llgi_data(fmodel, seed=19, feff_scale=1.3)
  fmodel.set_llgi_data(llgi_data)
  fmodel.set_target_name("llgi")
  fmodel.update_llgi_sigmaa_scatfrac()
  def forbidden(*args, **kwargs):
    raise AssertionError("ML alpha/beta computed in LLGI mode")
  for name in ["alpha_beta", "alpha_beta_w", "alpha_beta_t",
               "figures_of_merit", "phase_errors", "model_error_ml"]:
    setattr(fmodel, name, forbidden)

  info = fmodel.info(n_bins=5)
  assert info._llgi
  assert len(info.bins) > 0
  assert approx_equal(info.r_work, fmodel.r_work(), eps=1.e-10)
  assert approx_equal(info.r_free, fmodel.r_free(), eps=1.e-10)
  assert approx_equal(info.r_all, fmodel.r_all(), eps=1.e-10)
  mch = fmodel.map_calculation_helper_llgi()
  assert approx_equal(info.ml_phase_error,
    flex.mean(fmodel.phase_errors_llgi(mch)), eps=1.e-10)
  assert 0 < info.ml_phase_error < 90
  assert info.ml_coordinate_error > 0
  assert approx_equal(info.alpha_work_mean,
    flex.mean(mch.d.select(fmodel.arrays.work_sel)))
  assert 0 < info.fom_work_mean <= 1
  out = StringIO()
  info.show_all(out=out)
  text = out.getvalue()
  assert "LLGI (E-scale) estimates" in text
  assert "Acta Cryst. (1995)" not in text

def exercise_info_ml_mode_unchanged():
  # Ordinary ml-target info() (no llgi_data at all) still reports the ML
  # alpha/beta/phase error, and the ordinary R-factors.
  from six.moves import cStringIO as StringIO
  fmodel = build_fmodel(n_atoms=50, d_min=2.1, seed=43)
  info = fmodel.info(n_bins=5)
  assert not info._llgi
  assert approx_equal(info.r_work, fmodel.r_work(), eps=1.e-10)
  assert approx_equal(info.r_free, fmodel.r_free(), eps=1.e-10)
  assert approx_equal(info.r_all, fmodel.r_all(), eps=1.e-10)
  alpha_w, beta_w = fmodel.alpha_beta_w()
  assert approx_equal(info.alpha_work_mean, flex.mean(alpha_w.data()))
  assert approx_equal(info.ml_phase_error, flex.mean(fmodel.phase_errors()))
  out = StringIO()
  info.show_all(out=out)
  text = out.getvalue()
  assert "Acta Cryst. (1995)" in text
  assert "LLGI (E-scale)" not in text

def exercise():
  exercise_llgi_target_active_gating()
  exercise_info_in_llgi_mode()
  exercise_info_ml_mode_unchanged()
  print("OK")

if (__name__ == "__main__"):
  exercise()
