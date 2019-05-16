from __future__ import absolute_import, division, print_function
from cctbx import crystal
from cctbx import sgtbx
from cctbx import xray
from cctbx.array_family import flex
from cctbx import miller
from mmtbx import max_lik
from mmtbx.max_lik import maxlik
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times
import time
from cctbx.development import random_structure
from cctbx.eltbx.xray_scattering import wk1995
import random
import math
from six.moves import range

if (1): # fixed random seed to avoid rare failures
  random.seed(0)
  flex.set_random_seed(0)

def test_1():
  tstart = time.time()
  fraction_missing = 0.1
  d_min = 1.5
# create dummy model
  symmetry = crystal.symmetry(unit_cell=(15.67, 25.37, 35.68, 90, 90, 90),
                              space_group_symbol="P 21 21 21")
  structure = xray.structure(crystal_symmetry=symmetry)
  for k in range(1000):
    scatterer = xray.scatterer(
                 site = ((1.+k*abs(math.sin(k)))/1000.0,
                         (1.+k*abs(math.cos(k)))/1000.0,
                         (1.+ k)/1000.0),
                 u = abs(math.cos(k))*100./(8.*math.pi**2),
                 occupancy = 1.0,
                 scattering_type = "C")
    structure.add_scatterer(scatterer)

# partial model
  n_keep = int(round(structure.scatterers().size() * (1-fraction_missing)))
  partial_structure = xray.structure(special_position_settings=structure)
  partial_structure.add_scatterers(structure.scatterers()[:n_keep])

# fcalc (partial model), fobs (fcalc full model)
  f_calc = structure.structure_factors(d_min=d_min,
                               anomalous_flag=False).f_calc()
  f_calc_partial = partial_structure.structure_factors(d_min=d_min,
                                               anomalous_flag=False).f_calc()
  f_obs = abs(f_calc)
  f_calc= abs(f_calc_partial)
  d_star_sq = 1./flex.pow2(f_obs.d_spacings().data())

  assert approx_equal(flex.max(f_calc.data()),6810.19834824)
  assert approx_equal(flex.min(f_calc.data()),0.019589341727)
  assert approx_equal(flex.mean(f_calc.data()),76.651506629)
  assert approx_equal(flex.max(f_obs.data()),6962.58343229)
  assert approx_equal(flex.min(f_obs.data()),0.00111552904935)
  assert approx_equal(flex.mean(f_obs.data()),74.5148786464)
  assert f_obs.size() == f_calc.size()

# define test set reflections
  flags=flex.bool(f_calc_partial.indices().size(), False)
  k=0
  for i in range(f_calc_partial.indices().size()):
    k=k+1
    if (k !=10):
      flags[i]=False
    else:
      k=0
      flags[i]=True
  assert flags.count(True) == 250
  assert flags.count(False) == 2258
  assert flags.size() == 2508

# *********************************************************TEST = 1
  alpha, beta = maxlik.alpha_beta_est_manager(
                    f_obs  = f_obs,
                    f_calc = f_calc,
                    free_reflections_per_bin = 1000,
                    flags = flags,
                    interpolation = False,
                    epsilons = f_obs.epsilons().data().as_double()).alpha_beta()

  assert alpha.data().size() == beta.data().size()
  assert alpha.data().size() == f_obs.size()
  assert approx_equal(flex.min(alpha.data()),0.914152454693)
  assert approx_equal(flex.max(alpha.data()),0.914152454693)
  assert approx_equal(flex.min(beta.data()),818.503411782)
  assert approx_equal(flex.max(beta.data()),818.503411782)
# *********************************************************TEST = 2
  alpha, beta = maxlik.alpha_beta_est_manager(
    f_obs  = f_obs,
    f_calc = f_calc,
    free_reflections_per_bin = 50,
    flags = flags,
    interpolation = False,
    epsilons = f_obs.epsilons().data().as_double()).alpha_beta()

  assert alpha.data().size() == beta.data().size()
  assert alpha.data().size() == f_obs.size()
  assert approx_equal(flex.min(alpha.data()), 0.910350007113)
  assert approx_equal(flex.max(alpha.data()), 1.07104387776)
  assert approx_equal(flex.min(beta.data()), 21.7374310013)
  assert approx_equal(flex.max(beta.data()), 4222.81104745)
# *********************************************************TEST = 3
  alpha, beta = maxlik.alpha_beta_calc(
                    f   = f_obs,
                    n_atoms_absent  = 100,
                    n_atoms_included= 900,
                    bf_atoms_absent = 25.0,
                    final_error     = 0.0,
                    absent_atom_type = "C").alpha_beta()

  fom = max_lik.fom_and_phase_error(
    f_obs          = f_obs.data(),
    f_model        = flex.abs(f_calc.data()),
    alpha          = alpha.data(),
    beta           = beta.data(),
    epsilons       = f_obs.epsilons().data().as_double(),
    centric_flags  = f_obs.centric_flags().data()).fom()

  assert flex.max(fom) <= 1.0
  assert flex.min(fom) >= 0.0
  assert flex.min(alpha.data()) == flex.max(alpha.data()) == 1.0
  assert approx_equal(flex.min(beta.data()),7.964134920)
  assert approx_equal(flex.max(beta.data()),13695.1589364)
# *********************************************************TEST = 4

  xs = crystal.symmetry((3,4,5), "P 2 2 2")
  mi = flex.miller_index(((1,-2,3), (0,0,-4)))
  ms = miller.set(xs, mi)
  fc  = flex.double((1.,2.))
  fo  = flex.double((1.,2.))
  mso = miller.set(xs, mi)
  mac = miller.array(ms, fc)
  mao = miller.array(ms, fo)

  alp = flex.double(2,0.0)
  bet = flex.double(2,1e+9)
  fom = max_lik.fom_and_phase_error(
    f_obs          = mao.data(),
    f_model        = mac.data(),
    alpha          = alp,
    beta           = bet,
    epsilons       = mao.epsilons().data().as_double(),
    centric_flags  = mao.centric_flags().data()).fom()
  assert approx_equal(fom,[0.0, 0.0])
  alp = flex.double(2,1.0)
  bet = flex.double(2,1e+9)
  fom = max_lik.fom_and_phase_error(
    f_obs          = mao.data(),
    f_model        = mac.data(),
    alpha          = alp,
    beta           = bet,
    epsilons       = mao.epsilons().data().as_double(),
    centric_flags  = mao.centric_flags().data()).fom()

  assert approx_equal(fom,[0.0, 0.0])
  alp = flex.double(2,0.0)
  bet = flex.double(2,1e-9)
  fom = max_lik.fom_and_phase_error(
    f_obs          = mao.data(),
    f_model        = mac.data(),
    alpha          = alp,
    beta           = bet,
    epsilons       = mao.epsilons().data().as_double(),
    centric_flags  = mao.centric_flags().data()).fom()
  assert approx_equal(fom,[0.0, 0.0])
  alp = flex.double(2,1.0)
  bet = flex.double(2,1e-9)
  fom = max_lik.fom_and_phase_error(
    f_obs          = mao.data(),
    f_model        = mac.data(),
    alpha          = alp,
    beta           = bet,
    epsilons       = mao.epsilons().data().as_double(),
    centric_flags  = mao.centric_flags().data()).fom()
  assert approx_equal(fom,[1.0, 1.0])

def test_2():
  n_sites = 1000
  d_min   = 2.0
  volume_per_atom = 50
  fraction_missing = (0.0,)
  scale = 5.0

# create dummy model
  space_group_info = sgtbx.space_group_info("P212121")
  structure = random_structure.xray_structure(
                                 space_group_info  = space_group_info,
                                 elements          = ["N"]*(n_sites),
                                 volume_per_atom   = volume_per_atom,
                                 random_u_iso      = False)
  structure.scattering_type_registry(table="wk1995")
  f_calc = structure.structure_factors(d_min          = d_min,
                                       anomalous_flag = False,
                                       algorithm      = "direct").f_calc()
  f_obs = abs(f_calc)

  for fm in fraction_missing:
    # partial model
      n_keep = int(round(structure.scatterers().size() * (1-fm)))
      partial_structure = xray.structure(special_position_settings=structure)
      partial_structure.add_scatterers(structure.scatterers()[:n_keep])

    # fcalc (partial model), fobs (fcalc full model)
      f_calc_partial = partial_structure.structure_factors(d_min=d_min,
                                                   anomalous_flag=False,
                                                   algorithm = "direct").f_calc()
      f_calc= abs(f_calc_partial)

    # define test set reflections
      flags=flex.bool(f_calc_partial.indices().size(), False)
      k=0
      for i in range(f_calc_partial.indices().size()):
        k=k+1
        if (k !=10):
          flags[i]=False
        else:
          k=0
          flags[i]=True

    # *********************************************************TEST = 1
      alpha, beta = maxlik.alpha_beta_est_manager(
        f_obs           = f_obs,
        f_calc          = f_calc,
        free_reflections_per_bin = f_obs.data().size(),
        flags           = flags,
        interpolation   = False,
        epsilons        = f_obs.epsilons().data().as_double()).alpha_beta()

      assert alpha.size() == beta.size()
      assert alpha.size() == f_obs.size()
      assert approx_equal(flex.min(alpha.data()),1.0, 1.e-2)
      assert approx_equal(flex.max(alpha.data()),1.0, 1.e-2)
      assert approx_equal(flex.min(beta.data()) ,0.0, 1.e-2)
      assert approx_equal(flex.max(beta.data()) ,0.0, 1.e-2)

      alpha, beta = maxlik.alpha_beta_est_manager(
        f_obs           = f_obs,
        f_calc          = f_calc,
        free_reflections_per_bin = f_obs.data().size(),
        flags           = flags,
        interpolation   = True,
        epsilons        = f_obs.epsilons().data().as_double()).alpha_beta()

      assert alpha.size() == beta.size()
      assert alpha.size() == f_obs.size()
      assert approx_equal(flex.min(alpha.data()),1.0, 1.e-2)
      assert approx_equal(flex.max(alpha.data()),1.0, 1.e-2)
      assert approx_equal(flex.min(beta.data()) ,0.0, 1.e-2)
      assert approx_equal(flex.max(beta.data()) ,0.0, 1.e-2)
    # *********************************************************TEST = 2

      alpha, beta = maxlik.alpha_beta_est_manager(
        f_obs           = miller.array(miller_set = f_obs,
                                       data       = f_obs.data() * scale),
        f_calc          = f_calc,
        free_reflections_per_bin = 200,
        flags           = flags,
        interpolation   = False,
        epsilons        = f_obs.epsilons().data().as_double()).alpha_beta()

      assert alpha.size() == beta.size()
      assert alpha.size() == f_obs.size()
      assert approx_equal(flex.min(alpha.data()),1.0/scale, 1.e-2)
      assert approx_equal(flex.max(alpha.data()),1.0/scale, 1.e-2)
      assert approx_equal(flex.min(beta.data()) ,0.0, 1.e-2)
      assert approx_equal(flex.max(beta.data()) ,0.0, 1.e-2)

      alpha, beta = maxlik.alpha_beta_est_manager(
        f_obs           = miller.array(miller_set = f_obs,
                                       data       = f_obs.data() * scale),
        f_calc          = f_calc,
        free_reflections_per_bin = 200,
        flags           = flags,
        interpolation   = True,
        epsilons        = f_obs.epsilons().data().as_double()).alpha_beta()

      assert alpha.size() == beta.size()
      assert alpha.size() == f_obs.size()
      assert approx_equal(flex.min(alpha.data()),1.0/scale, 1.e-2)
      assert approx_equal(flex.max(alpha.data()),1.0/scale, 1.e-2)
      assert approx_equal(flex.min(beta.data()) ,0.0, 1.e-2)
      assert approx_equal(flex.max(beta.data()) ,0.0, 1.e-2)
    # *********************************************************TEST = 3

      alpha, beta = maxlik.alpha_beta_est_manager(
        f_obs           = miller.array(miller_set = f_obs,
                                       data       = f_obs.data() * scale),
        f_calc          = f_calc,
        free_reflections_per_bin = 200,
        flags           = flex.bool(f_obs.data().size(), True),
        interpolation   = False,
        epsilons        = f_obs.epsilons().data().as_double()).alpha_beta()

      assert alpha.size() == beta.size()
      assert alpha.size() == f_obs.size()
      assert approx_equal(flex.min(alpha.data()),1.0/scale, 1.e-2)
      assert approx_equal(flex.max(alpha.data()),1.0/scale, 1.e-2)
      assert approx_equal(flex.min(beta.data()) ,0.0, 1.e-2)
      assert approx_equal(flex.max(beta.data()) ,0.0, 1.e-2)

      alpha, beta = maxlik.alpha_beta_est_manager(
        f_obs           = miller.array(miller_set = f_obs,
                                       data       = f_obs.data() * scale),
        f_calc          = f_calc,
        free_reflections_per_bin = 200,
        flags           = flex.bool(f_obs.data().size(), True),
        interpolation   = True,
        epsilons        = f_obs.epsilons().data().as_double()).alpha_beta()

      assert alpha.size() == beta.size()
      assert alpha.size() == f_obs.size()
      assert approx_equal(flex.min(alpha.data()),1.0/scale, 1.e-2)
      assert approx_equal(flex.max(alpha.data()),1.0/scale, 1.e-2)
      assert approx_equal(flex.min(beta.data()) ,0.0, 1.e-2)
      assert approx_equal(flex.max(beta.data()) ,0.0, 1.e-2)
    # *********************************************************TEST = 4

      alpha, beta = maxlik.alpha_beta_est_manager(
        f_obs           = miller.array(miller_set = f_obs,
                                       data       = f_obs.data() * scale),
        f_calc          = f_calc,
        free_reflections_per_bin = 200,
        flags           = flex.bool(f_obs.data().size(), False),
        interpolation   = False,
        epsilons        = f_obs.epsilons().data().as_double()).alpha_beta()

      assert alpha.size() == beta.size()
      assert alpha.size() == f_obs.size()
      assert approx_equal(flex.min(alpha.data()),1.0/scale, 1.e-2)
      assert approx_equal(flex.max(alpha.data()),1.0/scale, 1.e-2)
      assert approx_equal(flex.min(beta.data()) ,0.0, 1.e-2)
      assert approx_equal(flex.max(beta.data()) ,0.0, 1.e-2)

      alpha, beta = maxlik.alpha_beta_est_manager(
        f_obs           = miller.array(miller_set = f_obs,
                                       data       = f_obs.data() * scale),
        f_calc          = f_calc,
        free_reflections_per_bin = 200,
        flags           = flex.bool(f_obs.data().size(), False),
        interpolation   = True,
        epsilons        = f_obs.epsilons().data().as_double()).alpha_beta()

      assert alpha.size() == beta.size()
      assert alpha.size() == f_obs.size()
      assert approx_equal(flex.min(alpha.data()),1.0/scale, 1.e-2)
      assert approx_equal(flex.max(alpha.data()),1.0/scale, 1.e-2)
      assert approx_equal(flex.min(beta.data()) ,0.0, 1.e-2)
      assert approx_equal(flex.max(beta.data()) ,0.0, 1.e-2)
    # *********************************************************TEST = 5

def test_3():

   symmetry = crystal.symmetry(unit_cell          = (15.67, 25.37, 35.68, 90, 90, 90),
                               space_group_symbol = "P 21 21 21")
   structure = xray.structure(crystal_symmetry = symmetry)
   mi = structure.structure_factors(d_min          = 1.5,
                                    anomalous_flag = False).f_calc().indices()
   # ================================================================= TEST-1
   alpha  = flex.double(mi.size())
   beta   = flex.double(mi.size())
   d_obs  = flex.double(mi.size())
   d_calc = flex.double(mi.size())

   for i in range(1,mi.size()+1):
     d_obs [i-1] = i*1.0
     d_calc[i-1] = i*1.0
     beta  [i-1] = i*500.0
     alpha [i-1] = float(i) / float(i + 1)

   obj = max_lik.f_star_w_star_mu_nu(f_obs          = d_obs,
                                     f_model        = d_calc,
                                     alpha          = alpha,
                                     beta           = beta,
                                     space_group    = symmetry.space_group(),
                                     miller_indices = mi)
   f_star = obj.f_star()
   w_star = obj.w_star()
   mu     = obj.mu()
   nu     = obj.nu()
   nzero  = obj.number_of_f_star_zero()

   assert approx_equal(flex.max(f_star) ,          2505.77677201 , 1.e-4)
   assert approx_equal(flex.min(f_star) ,          0.0           , 1.e-4)
   assert approx_equal(flex.mean(f_star),          1085.99060715 , 1.e-4)
   assert approx_equal(flex.max(w_star) ,          1.0           , 1.e-4)
   assert approx_equal(flex.min(w_star) ,          0.0           , 1.e-4)
   assert approx_equal(flex.mean(w_star),          0.01782658613 , 1.e-4)
   assert approx_equal(flex.max(mu)     ,          2.23810354633 , 1.e-4)
   assert approx_equal(flex.min(mu)     ,          0.0           , 1.e-4)
   assert approx_equal(flex.mean(mu)    ,          1.20159615933 , 1.e-4)
   assert approx_equal(flex.max(nu)     ,          0.999107484116, 1.e-4)
   assert approx_equal(flex.min(nu)     ,          0.0           , 1.e-4)
   assert approx_equal(flex.mean(nu)    ,          0.745699513719, 1.e-4)
   assert approx_equal(nzero            ,          501           )

def test_4():

   symmetry = crystal.symmetry(unit_cell          = (15.67, 25.37, 35.68, 90, 90, 90),
                               space_group_symbol = "P 21 21 21")
   structure = xray.structure(crystal_symmetry = symmetry)
   ma = structure.structure_factors(d_min          = 1.5,
                                    anomalous_flag = False).f_calc()
   mi = ma.indices()
   # ================================================================= TEST-1
   alpha  = flex.double(mi.size())
   beta   = flex.double(mi.size())
   d_obs  = flex.double(mi.size())
   d_calc = flex.double(mi.size())

   # define test set reflections
   flags=flex.int(beta.size(), 0)
   k=0
   for i in range(flags.size()):
     k=k+1
     if (k !=10):
       flags[i]=0
     else:
       k=0
       flags[i]=1

   for i in range(1,mi.size()+1):
     d_obs [i-1] = i*1.5
     d_calc[i-1] = i*1.0
     beta  [i-1] = i*500.0
     alpha [i-1] = float(i) / float(i + 1)

   obj = max_lik.fom_and_phase_error(
     f_obs          = d_obs,
     f_model        = d_calc,
     alpha          = alpha,
     beta           = beta,
     epsilons       = ma.epsilons().data().as_double(),
     centric_flags  = ma.centric_flags().data())
   per = obj.phase_error()
   fom = obj.fom()
   assert approx_equal(flex.max(per) ,  89.9325000127     , 1.e-4)
   assert approx_equal(flex.min(per) ,  5.37565067746e-05 , 1.e-4)
   assert approx_equal(flex.mean(per),  20.7942460698     , 1.e-4)
   assert approx_equal(flex.max(fom) ,  0.999999402705    , 1.e-4)
   assert approx_equal(flex.min(fom) ,  0.000749999859375 , 1.e-4)
   assert approx_equal(flex.mean(fom),  0.858269037582    , 1.e-4)

def run():
  test_1()
  test_2()
  test_3()
  test_4()
  print(format_cpu_times())

if (__name__ == "__main__"):
  run()
