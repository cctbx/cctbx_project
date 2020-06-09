from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from scitbx import matrix
import math
from libtbx import adopt_init_args
import scitbx.lbfgs
from mmtbx.bulk_solvent import kbu_refinery
from cctbx import maptbx
import mmtbx.masks
import boost.python
asu_map_ext = boost.python.import_ext("cctbx_asymmetric_map_ext")
from libtbx import group_args
from mmtbx import bulk_solvent
from mmtbx.ncs import tncs
from libtbx.test_utils import approx_equal
from collections import OrderedDict

# Utilities used by algorithm 2 ------------------------------------------------

class minimizer(object):
  def __init__(self, max_iterations, calculator):
    adopt_init_args(self, locals())
    self.x = self.calculator.x
    self.cntr=0
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=scitbx.lbfgs.termination_parameters(
        max_iterations=max_iterations))

  def compute_functional_and_gradients(self):
    self.cntr+=1
    self.calculator.update_target_and_grads(x=self.x)
    t = self.calculator.target()
    g = self.calculator.gradients()
    #print "step: %4d"%self.cntr, "target:", t, "params:", \
    #  " ".join(["%10.6f"%i for i in self.x]), math.log(t)
    return t,g

class minimizer2(object):

  def __init__(self, calculator, min_iterations=0, max_iterations=2000):
    adopt_init_args(self, locals())
    self.x = self.calculator.x
    self.n = self.x.size()
    self.cntr=0

  def run(self, use_curvatures=0):
    self.minimizer = kbu_refinery.lbfgs_run(
      target_evaluator=self,
      min_iterations=self.min_iterations,
      max_iterations=self.max_iterations,
      use_curvatures=use_curvatures)
    self(requests_f_and_g=True, requests_diag=False)
    return self

  def __call__(self, requests_f_and_g, requests_diag):
    self.cntr+=1
    self.calculator.update_target_and_grads(x=self.x)
    if (not requests_f_and_g and not requests_diag):
      requests_f_and_g = True
      requests_diag = True
    if (requests_f_and_g):
      self.f = self.calculator.target()
      self.g = self.calculator.gradients()
      self.d = None
    if (requests_diag):
      self.d = self.calculator.curvatures()
      #assert self.d.all_ne(0)
      if(self.d.all_eq(0)): self.d=None
      else:
        self.d = 1 / self.d
    #print "step: %4d"%self.cntr, "target:", self.f, "params:", \
    #  " ".join(["%10.6f"%i for i in self.x]) #, math.log(self.f)
    return self.x, self.f, self.g, self.d

class tg(object):
  def __init__(self, x, i_obs, F, use_curvatures):
    self.x = x
    self.i_obs = i_obs
    self.F = F
    self.t = None
    self.g = None
    self.d = None
    self.sum_i_obs = flex.sum(self.i_obs.data())
    self.use_curvatures=use_curvatures
    self.update_target_and_grads(x=x)

  def update(self, x):
    self.update_target_and_grads(x = x)

  def update_target_and_grads(self, x):
    self.x = x
    s = 1 #180/math.pi
    i_model = flex.double(self.i_obs.data().size(),0)
    for n, kn in enumerate(self.x):
      for m, km in enumerate(self.x):
        tmp = self.F[n].data()*flex.conj(self.F[m].data())
        i_model += kn*km*flex.real(tmp)
        #pn = self.F[n].phases().data()*s
        #pm = self.F[m].phases().data()*s
        #Fn = flex.abs(self.F[n].data())
        #Fm = flex.abs(self.F[m].data())
        #i_model += kn*km*Fn*Fm*flex.cos(pn-pm)
    diff = i_model - self.i_obs.data()
    t = flex.sum(diff*diff)/4
    #
    g = flex.double()
    for j in range(len(self.F)):
      tmp = flex.double(self.i_obs.data().size(),0)
      for m, km in enumerate(self.x):
        tmp += km * flex.real( self.F[j].data()*flex.conj(self.F[m].data()) )
        #pj = self.F[j].phases().data()*s
        #pm = self.F[m].phases().data()*s
        #Fj = flex.abs(self.F[j].data())
        #Fm = flex.abs(self.F[m].data())
        #tmp += km * Fj*Fm*flex.cos(pj-pm)
      g.append(flex.sum(diff*tmp))
    self.t = t
    self.g = g
    #
    if self.use_curvatures:
      d = flex.double()
      for j in range(len(self.F)):
        tmp1 = flex.double(self.i_obs.data().size(),0)
        tmp2 = flex.double(self.i_obs.data().size(),0)
        for m, km in enumerate(self.x):
          zz = flex.real( self.F[j].data()*flex.conj(self.F[m].data()) )
          tmp1 += km * zz
          tmp2 += zz
          #pj = self.F[j].phases().data()*s
          #pm = self.F[m].phases().data()*s
          #Fj = flex.abs(self.F[j].data())
          #Fm = flex.abs(self.F[m].data())
          #tmp += km * Fj*Fm*flex.cos(pj-pm)
        d.append(flex.sum(tmp1*tmp1 + tmp2))
      self.d=d

  def target(self): return self.t/self.sum_i_obs

  def gradients(self): return self.g/self.sum_i_obs

  def gradient(self): return self.gradients()

  def curvatures(self): return self.d/self.sum_i_obs
#-------------------------------------------------------------------------------

def mosaic_f_mask(miller_array, xray_structure, step, volume_cutoff):
  # compute mask in p1 (via ASU)
  cg = maptbx.crystal_gridding(
    unit_cell        = xray_structure.unit_cell(),
    space_group_info = xray_structure.space_group_info(),
    symmetry_flags   = maptbx.use_space_group_symmetry,
    step             = step)
  n_real = cg.n_real()
  mask_p1 = mmtbx.masks.mask_from_xray_structure(
    xray_structure        = xray_structure,
    p1                    = True,
    for_structure_factors = True,
    n_real                = n_real,
    in_asu                = False).mask_data
  maptbx.unpad_in_place(map=mask_p1)
  # conn analysis
  co = maptbx.connectivity(
    map_data                   = mask_p1,
    threshold                  = 0.01,
    preprocess_against_shallow = True,
    wrapping                   = False)
  conn = co.result().as_double()
  z = zip(co.regions(),range(0,co.regions().size()))
  sorted_by_volume = sorted(z, key=lambda x: x[0], reverse=True)
  f_mask_data = flex.complex_double(miller_array.data().size(), 0)
  #f_masks  = []
  FM = OrderedDict()
  for p in sorted_by_volume:
    v, i = p
    volume = v*step**3
    if(volume_cutoff is not None):
      if volume < volume_cutoff: continue
    if(i==0): continue
    mask_i = flex.double(flex.grid(n_real), 0)
    mask_i_asu = asu_map_ext.asymmetric_map(
      miller_array.space_group().type(), mask_i.set_selected(conn==i, 1)).data()
    volume_asu = (mask_i_asu>0).count(True)*step**3
    print("volume p1: %12.3f"%volume, "volume asu: %12.3f"%volume_asu,
          "region_id", i)
    if(volume_asu<1.e-6): continue
    f_mask_i = miller_array.structure_factors_from_asu_map(
      asu_map_data = mask_i_asu, n_real = n_real)
    #f_masks.append(f_mask_i)
    f_mask_data += f_mask_i.data()
    FM.setdefault(round(volume, 3), []).append(f_mask_i.data())
  # group asu pices corresponding to the same region in P1
  F_MASKS = []
  for k,v in zip(FM.keys(), FM.values()):
    tmp = flex.complex_double(miller_array.data().size(), 0)
    for v_ in v:
      tmp+=v_
    F_MASKS.append(miller_array.customized_copy(data = tmp))
  #
  f_mask = miller_array.customized_copy(data = f_mask_data)
  #
  return group_args(f_mask = f_mask, f_masks = F_MASKS, n_real=n_real)

def algorithm_0(f_obs, fc, f_masks):
  """
  Grid search
  """
  k_mask_trial_range=[]
  s = 0
  while s<1:
    k_mask_trial_range.append(s)
    s+=0.001
  r = []
  fc_data = fc.data()
  for i, f_mask in enumerate(f_masks):
    #print("mask ",i)
    assert f_obs.data().size() == fc.data().size()
    assert f_mask.data().size() == fc.data().size()
    #print (bulk_solvent.r_factor(f_obs.data(),fc_data))
    kmask_, k_ = \
      bulk_solvent.k_mask_and_k_overall_grid_search(
        f_obs.data(),
        fc_data,
        f_mask.data(),
        flex.double(k_mask_trial_range),
        flex.bool(fc.data().size(),True))
    r.append(kmask_)
    fc_data += fc_data*k_ + kmask_*f_mask.data()
    #print (bulk_solvent.r_factor(f_obs.data(),fc_data + kmask_*f_mask.data(),k_))
  r = [1,]+r
  return r

def algorithm_2(i_obs, fc, f_masks, x, use_curvatures=True, macro_cycles=10):
  """
  Unphased one-step search
  """
  F = [fc]+f_masks
  calculator = tg(i_obs = i_obs, F=F, x = x, use_curvatures=use_curvatures)
  for it in xrange(macro_cycles):
    if(use_curvatures):
      m = minimizer(max_iterations=100, calculator=calculator)
    else:
      upper = flex.double([10.] + [0.65]*(x.size()-1))
      lower = flex.double([0.1] + [0]*(x.size()-1))
      m = tncs.minimizer(
        potential       = calculator,
        use_bounds      = 2,
        lower_bound     = lower,
        upper_bound     = upper,
        initial_values  = x).run()
    calculator = tg(i_obs = i_obs, F=F, x = m.x, use_curvatures=use_curvatures)
  if(use_curvatures):
    for it in range(10):
      m = minimizer(max_iterations=100, calculator=calculator)
      calculator = tg(i_obs = i_obs, F=F, x = m.x, use_curvatures=use_curvatures)
      m = minimizer2(max_iterations=100, calculator=calculator).run(use_curvatures=True)
      calculator = tg(i_obs = i_obs, F=F, x = m.x, use_curvatures=use_curvatures)
  return m.x

def algorithm_3(i_obs, fc, f_masks):
  """
  Unphased two-step search
  """
  F = [fc]+f_masks
  Gnm = []
  cs = {}
  cntr=0
  nm=[]
  # Compute and store Gnm
  for n, Fn in enumerate(F):
    for m, Fm in enumerate(F):
      if m < n:
        continue
      Gnm.append( flex.real( Fn.data()*flex.conj(Fm.data()) ) )
      cs[(n,m)] = cntr
      cntr+=1
      nm.append((n,m))
  # Keep track of indices for "upper triangular matrix vs full"
  for k,v in zip(cs.keys(), cs.values()):
    i,j=k
    if i==j: continue
    else: cs[(j,i)]=v
  # Generate and solve system Ax=b, x = A_1*b
  A = []
  b = []
  for u, Gnm_u in enumerate(Gnm):
    for v, Gnm_v in enumerate(Gnm):
      scale = 2
      n,m=nm[v]
      if n==m: scale=1
      A.append( flex.sum(Gnm_u*Gnm_v)*scale )
    b.append( flex.sum(Gnm_u * i_obs.data()) )
  A = matrix.sqr(A)
  A_1 = A.inverse()
  b = matrix.col(b)
  x = A_1 * b
  # Expand Xmn from solution x
  Xmn = []
  for n, Fn in enumerate(F):
    rows = []
    for m, Fm in enumerate(F):
      x_ = x[cs[(n,m)]]
      rows.append(x_)
    Xmn.append(rows)
  # Do formula (19)
  lnK = []
  for j, Fj in enumerate(F):
    t1 = flex.sum( flex.log( flex.double(Xmn[j]) ) )
    t2 = 0
    for n, Fn in enumerate(F):
      for m, Fm in enumerate(F):
        t2 += math.log(Xmn[n][m])
    t2 = t2 / (2*len(F))
    lnK.append( 1/len(F)*(t1-t2) )
  return [math.exp(x) for x in lnK]

def algorithm_4(f_obs, fc, f_masks, max_cycles=100, auto_converge_eps=1.e-7):
  """
  Phased simultaneous search
  """
  fc = fc.deep_copy()
  F = [fc]+f_masks
  x_res = None
  cntr = 0
  x_prev = None
  while True:
    f_obs_cmpl = f_obs.phase_transfer(phase_source=fc)
    A = []
    b = []
    for j, Fj in enumerate(F):
      A_rows = []
      for n, Fn in enumerate(F):
        Gjn = flex.real( Fj.data()*flex.conj(Fn.data()) )
        A_rows.append( flex.sum(Gjn) )
      Hj = flex.real( Fj.data()*flex.conj(f_obs_cmpl.data()) )
      b.append(flex.sum(Hj))
      A.extend(A_rows)
    A = matrix.sqr(A)
    A_1 = A.inverse()
    b = matrix.col(b)
    x = A_1 * b
    if x_res is None: x_res  = flex.double(x)
    else:             x_res += flex.double(x)
    x_ = [x[0]] + list(x_res[1:])
    #print "iteration:", cntr, " ".join(["%10.6f"%i for i in x_])
    #
    fc_d = fc.data()
    for i, f in enumerate(F):
      if i == 0: continue
      fc_d += x[i]*f.data()
    fc = fc.customized_copy(data = fc_d)
    cntr+=1
    if(cntr>max_cycles): break
    if(x_prev is None): x_prev = x_[:]
    else:
      max_diff = flex.max(flex.abs(flex.double(x_prev)-flex.double(x_)))
      if(max_diff<=auto_converge_eps): break
      x_prev = x_[:]
  return x_
