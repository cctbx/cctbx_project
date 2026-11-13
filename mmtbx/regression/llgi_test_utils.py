from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx.development import random_structure
from cctbx import sgtbx
from libtbx import group_args
import mmtbx.f_model
import numpy as np
import random

""" Shared fixtures for the LLGI regression tests (tst_llgi_*.py). """

def build_fmodel(n_atoms, d_min, seed=0, space_group="P 21 21 21",
      update_scales=True):
  """ fmodel for a random structure of 3*n_atoms O/N/C atoms, with
  f_obs = |Fcalc| to d_min and 10% free flags; update_all_scales() is run
  unless update_scales is False. """
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
  if(update_scales):
    fmodel.update_all_scales()
  return fmodel

def synthetic_llgi_data(fmodel, seed=1, feff_scale=1.0):
  """ Synthetic nacelle-like DOBS/FEFF/TEPS/RESN on fmodel.f_obs()'s
  current index set. Feff is f_obs perturbed by 0.85-1.15 and scaled by
  feff_scale (so a test can tell whether a result depends on Feff rather
  than f_obs); RESN is sqrt(epsilon)*O(mean f_obs), so Eeff = Feff/RESN is
  of order 1, as in real nacelle output (otherwise the Bessel arguments
  are far too large). TEPS = 1. """
  f_obs = fmodel.f_obs()
  n = f_obs.size()
  epsilons = f_obs.epsilons().data().as_double()
  rnd = random.Random(seed)
  dobs = f_obs.array(data=flex.double([0.5 + 0.4 * rnd.random()
    for i in range(n)]))
  feff_data = f_obs.data() * feff_scale * flex.double(
    [0.85 + 0.3 * rnd.random() for i in range(n)])
  feff = f_obs.array(data=feff_data)
  teps = f_obs.array(data=flex.double(n, 1.0))
  resn = f_obs.array(data=flex.sqrt(epsilons) * flex.double(
    [rnd.uniform(2.0, 6.0) for i in range(n)]) * flex.mean(f_obs.data()))
  return group_args(dobs=dobs, feff=feff, teps=teps, resn=resn, info=None)

def add_hybrid_data(fmodel, llgi_data, rice_kappa, seed=7):
  """ Synthetic intensities consistent with the synthetic Feff/RESN:
  E_obs^2 = Eeff^2 + noise, sigma(E_obs^2) in [0.2, 2]. """
  import mmtbx.refinement.llgi_hybrid as llgi_hybrid
  rnd = random.Random(seed)
  f_obs = fmodel.f_obs()
  eeff = llgi_data.feff.data() / llgi_data.resn.data()
  n = eeff.size()
  sig = flex.double([0.2 + 1.8 * rnd.random() for i in range(n)])
  e2 = eeff * eeff + sig * flex.double([rnd.gauss(0, 1) for i in range(n)])
  params = llgi_hybrid.llgi_hybrid_params.extract()
  params.rice_kappa = rice_kappa
  flags = llgi_hybrid.force_exact_flags(e2, sig,
    f_obs.centric_flags().data())
  return llgi_hybrid.replace_llgi_data(llgi_data,
    e_obs_sq=f_obs.array(data=e2), sig_e_obs_sq=f_obs.array(data=sig),
    force_exact=f_obs.array(data=flags.force_exact), hybrid_params=params)

def build_llgi_fmodel(n_atoms, d_min, seed=0, feff_scale=1.0,
      rice_kappa=None):
  """ A scaled fmodel with target=llgi, synthetic llgi_data attached and
  sigmaa/scatfrac fitted by the real update_llgi_sigmaa_scatfrac(), i.e.
  the state phenix.refine has when it computes maps and statistics. With
  rice_kappa, synthetic intensities are added (add_hybrid_data) so the
  hybrid exact/Rice likelihood is active. """
  fmodel = build_fmodel(n_atoms, d_min, seed=seed)
  llgi_data = synthetic_llgi_data(fmodel, seed=seed + 100,
    feff_scale=feff_scale)
  if(rice_kappa is not None):
    llgi_data = add_hybrid_data(fmodel, llgi_data, rice_kappa)
  fmodel.set_llgi_data(llgi_data)
  fmodel.set_target_name("llgi")
  fmodel.update_llgi_sigmaa_scatfrac()
  return fmodel

def random_theta_and_grid(rnd, k):
  """ Random D_model parameters theta = [a_1..a_k, b, B_defect] and a
  sorted B_k ladder, in ranges that keep D_model well below 1. """
  a = rnd.uniform(0.05, 0.6, size=k)
  b_k_grid = np.sort(rnd.uniform(1.0, 100.0, size=k))
  b = rnd.uniform(0.02, 0.3)
  b_defect = rnd.uniform(5.0, 150.0)
  theta = np.empty(k + 2, dtype=float)
  theta[:k] = a
  theta[-2] = b
  theta[-1] = b_defect
  return theta, b_k_grid
