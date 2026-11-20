from __future__ import absolute_import, division, print_function
import time
import random
from cctbx.array_family import flex
from cctbx import adptbx
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
import iotbx.pdb
import mmtbx.model
import mmtbx.f_model
from mmtbx.refinement import wrappers
from mmtbx.regression import tst_ordered_solvent_altlocs as t

"""
mmtbx.refinement.wrappers.unrestrained_qbr_fsr: occupancy and B refinement
of part of the model with structure factors of the rest kept fixed
(freeze_unselected=True) gives the same result as the regular calculation.
"""

def get_model(space_group="P 21 21 21", n_water=30):
  random.seed(1)
  flex.set_random_seed(1)
  lines = t.ser_lines()
  k = 0
  for ix in range(5):
    for iy in range(3):
      for iz in range(2):
        if(k >= n_water): break
        xyz = (2.+ix*3.3, 14.+iy*3.1, 15.+iz*3.2)
        lines.append(t.hoh(k+1, xyz, occ=0.4+0.6*random.random(),
          b=10+20*random.random()))
        k += 1
  cs = "CRYST1   22.000   24.000   26.000  90.00  90.00  90.00 %s"%space_group
  return mmtbx.model.manager(model_input = iotbx.pdb.input(
    source_info=None, lines=[cs]+lines), log=null_out())

def get_fmodel(model, d_min, anomalous=False, twin_law=None,
               target_name="ls_wunit_kunit"):
  """
  Least-squares target: its minimum is at the true parameters whatever the
  starting model (not so for maximum-likelihood with error model of a poor
  start, when refined atoms are half of the model as they are here).
  """
  model.setup_scattering_dictionaries(scattering_table="n_gaussian")
  xrs = model.get_xray_structure().deep_copy_scatterers()
  if(anomalous):
    # make data anomalous
    for sc in xrs.scatterers():
      if(sc.scattering_type.strip() == "O"): sc.fdp = 0.3
  f_obs = abs(xrs.structure_factors(d_min=d_min,
    anomalous_flag=anomalous).f_calc())
  r_free_flags = f_obs.generate_r_free_flags()
  fmodel = mmtbx.f_model.manager(
    f_obs          = f_obs,
    r_free_flags   = r_free_flags,
    twin_law       = twin_law,
    target_name    = target_name,
    xray_structure = model.get_xray_structure())
  fmodel.update_all_scales(remove_outliers=False)
  return fmodel

def perturb(model, fmodel, selection):
  """ Start: occupancy 1, B = 25 for selected atoms """
  xrs = fmodel.xray_structure
  occ = xrs.scatterers().extract_occupancies()
  occ.set_selected(selection, 1.0)
  xrs.set_occupancies(value = occ, selection = selection)
  u = xrs.extract_u_iso_or_u_equiv()
  u.set_selected(selection, adptbx.b_as_u(25.))
  xrs.set_u_iso(values = u, selection = selection)
  fmodel.update_xray_structure(update_f_calc = True)
  model.set_xray_structure(xrs)

def run_one(freeze, anomalous=False, twin_law=None, refine_q=True,
            refine_b=True, d_min=1.5):
  model = get_model()
  fmodel = get_fmodel(model, d_min=d_min, anomalous=anomalous,
    twin_law=twin_law)
  assert fmodel.r_work() < 0.02, fmodel.r_work()
  selection = model.solvent_selection()
  occ_true = fmodel.xray_structure.scatterers().extract_occupancies()
  b_true = fmodel.xray_structure.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1)
  perturb(model, fmodel, selection)
  r_start = fmodel.r_work()
  assert r_start > 0.1, r_start
  o = wrappers.unrestrained_qbr_fsr(
    fmodel            = fmodel,
    model             = model,
    selection         = selection,
    refine_xyz        = False,
    refine_q          = refine_q,
    refine_b          = refine_b,
    q_min             = 0.004,
    q_max             = 1.0,
    b_min             = 1,
    b_max             = 60,
    freeze_unselected = freeze)
  xrs = fmodel.xray_structure
  occ = xrs.scatterers().extract_occupancies()
  b = xrs.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1)
  # nothing but selected atoms changed
  assert approx_equal(occ.select(~selection), occ_true.select(~selection))
  assert approx_equal(b.select(~selection), b_true.select(~selection))
  # bounds
  assert flex.min(occ.select(selection)) >= 0.004
  assert flex.max(occ.select(selection)) <= 1.0+1.e-9
  assert flex.max(b.select(selection)) <= 60+1.e-6
  # fmodel is in sync with its structure
  r = fmodel.r_work()
  fmodel.update_xray_structure(update_f_calc = True)
  assert approx_equal(r, fmodel.r_work(), 1.e-6)
  return group(r_start=r_start, r=r, occ=occ.select(selection),
    b=b.select(selection), occ_true=occ_true.select(selection),
    b_true=b_true.select(selection), o=o)

class group(object):
  def __init__(self, **kw): self.__dict__.update(kw)

def exercise_frozen_equals_full():
  """
  Joint occupancy and B refinement with fixed contribution of the rest of
  the model recovers the answer, with fewer structure factor calculations
  than alternating occupancy and B refinement that gets close to it.
  """
  full   = run_one(freeze = False)
  frozen = run_one(freeze = True)
  assert approx_equal(full.r_start, frozen.r_start, 1.e-6)
  assert frozen.r < 0.01, frozen.r
  assert approx_equal(frozen.occ, frozen.occ_true, 0.02)
  assert approx_equal(frozen.b, frozen.b_true, 1.5)
  assert full.r < 0.05, full.r
  assert approx_equal(full.occ, full.occ_true, 0.2)
  assert frozen.r <= full.r
  assert frozen.o.n_calls < 150, frozen.o.n_calls

def exercise_one_of_two():
  """ occupancy only or B only: the same result both ways """
  for refine_q, refine_b in [(True, False), (False, True)]:
    full   = run_one(freeze=False, refine_q=refine_q, refine_b=refine_b)
    frozen = run_one(freeze=True,  refine_q=refine_q, refine_b=refine_b)
    assert approx_equal(full.r, frozen.r, 3.e-3)
    assert approx_equal(full.occ, frozen.occ, 0.02)
    # B above 40 are poorly defined by the data: the two ways stop up to
    # 1.8 apart there (0.5 before the table of exponents used to compute
    # structure factors with FFT got interpolated)
    assert approx_equal(full.b, frozen.b, 3.0)
    if(not refine_q): assert approx_equal(frozen.occ, [1.0]*frozen.occ.size())
    if(not refine_b): assert approx_equal(frozen.b, [25.]*frozen.b.size())

def exercise_anomalous():
  full   = run_one(freeze = False, anomalous = True)
  frozen = run_one(freeze = True,  anomalous = True)
  assert frozen.r < 0.02, frozen.r
  assert approx_equal(frozen.occ, frozen.occ_true, 0.02)
  assert frozen.r <= full.r

def exercise_twin():
  """ twinned data or coordinate refinement: regular calculation is used """
  called = []
  class recorder(wrappers.unrestrained_qbr_fsr):
    def _run_all(self):    called.append("all")
    def _run_frozen(self): called.append("frozen")
  class fmodel(object):
    twin_law = None
    sfg_params = group(algorithm = "fft")
  selection = flex.bool([True, False])
  recorder(fmodel=fmodel, model=None, selection=selection, refine_xyz=False,
    freeze_unselected=True)
  recorder(fmodel=fmodel, model=None, selection=selection, refine_xyz=True,
    freeze_unselected=True)
  recorder(fmodel=fmodel, model=None, selection=selection, refine_xyz=False,
    freeze_unselected=False)
  fmodel.twin_law = "k,h,-l"
  recorder(fmodel=fmodel, model=None, selection=selection, refine_xyz=False,
    freeze_unselected=True)
  assert called == ["frozen", "all", "all", "all"], called

def exercise_zero_occupancy_and_b():
  """
  Atoms that start with occupancy 0 and B 0 (this is what a coarse search
  can leave) are refined as any other: used to crash.
  """
  model = get_model()
  fmodel = get_fmodel(model, d_min=1.5)
  selection = model.solvent_selection()
  i_seqs = selection.iselection()
  xrs = fmodel.xray_structure
  occ_true = xrs.scatterers().extract_occupancies()
  perturb(model, fmodel, selection)
  for i in i_seqs[:2]:
    xrs.scatterers()[i].occupancy = 0
    xrs.scatterers()[i].u_iso = 0
  fmodel.update_xray_structure(update_f_calc = True)
  model.set_xray_structure(xrs)
  wrappers.unrestrained_qbr_fsr(
    fmodel            = fmodel,
    model             = model,
    selection         = selection,
    refine_xyz        = False,
    q_min             = 0.004,
    q_max             = 1.0,
    b_min             = 1,
    b_max             = 60,
    freeze_unselected = True)
  assert fmodel.r_work() < 0.02, fmodel.r_work()
  occ = fmodel.xray_structure.scatterers().extract_occupancies()
  assert approx_equal(occ.select(selection), occ_true.select(selection), 0.05)

def exercise_gradients():
  """
  Calculator for joint occupancy and B refinement: at start fmodel is what
  it was, gradients are equal to finite differences.
  """
  for algorithm in ["direct", "fft"]:
    model = get_model(n_water=6)
    fmodel = get_fmodel(model, d_min=1.5)
    selection = model.solvent_selection()
    perturb(model, fmodel, selection)
    fmodel.sfg_params.algorithm = algorithm
    fmodel.update_xray_structure(update_f_calc = True)
    r_work = fmodel.r_work()
    selection_u = selection.deep_copy()
    selection_u[selection.iselection()[2]] = False # as anisotropic atom
    c = wrappers.occupancy_and_u_iso_calculator(
      fmodel          = fmodel,
      selection       = selection,
      selection_u_iso = selection_u,
      q_min           = 0.004,
      q_max           = 1.0,
      u_min           = adptbx.b_as_u(1),
      u_max           = adptbx.b_as_u(60))
    assert approx_equal(fmodel.r_work(), r_work, 2.e-3)
    x = c.x.deep_copy()
    n = selection.count(True)
    assert x.size() == 2*n-1
    assert approx_equal(c.target(), 1)
    g = c.gradients().deep_copy()
    for i in range(x.size()):
      # FFT: target is not smooth in u_iso at small scale (sampling of an
      # atom changes in steps), finite differences need a larger step
      eps = 1.e-4
      if(algorithm == "fft" and i >= n): eps = 1.e-2
      ts = []
      for sign in [1, -1]:
        x_ = x.deep_copy()
        x_[i] += sign*eps
        c.update(x = x_)
        ts.append(c.target())
      fd = (ts[0]-ts[1])/(2*eps)
      tolerance = 1.e-4
      if(algorithm == "fft"): tolerance = max(5.e-3, abs(g[i])*0.05)
      assert approx_equal(fd, g[i], tolerance), (algorithm, i, fd, g[i])

def run():
  exercise_zero_occupancy_and_b()
  exercise_gradients()
  exercise_frozen_equals_full()
  exercise_one_of_two()
  exercise_anomalous()
  exercise_twin()

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %6.2f"%(time.time()-t0))
