from __future__ import division
import cctbx.array_family.flex # import dependency
import boost.python
ext = boost.python.import_ext("mmtbx_ncs_ext")
from scitbx.array_family import flex
from cctbx import sgtbx
from libtbx import adopt_init_args
from scitbx import lbfgsb

def lbfgs_run(target_evaluator, use_bounds):
  minimizer = lbfgsb.minimizer(
    n   = target_evaluator.n,
    l   = flex.double(target_evaluator.n, 0), # lower bound
    u   = flex.double(target_evaluator.n, 1), # upper bound
    nbd = flex.int(target_evaluator.n, use_bounds))    # apply both bounds
  minimizer.error = None
  try:
    icall = 0
    while 1:
      icall += 1
      x, f, g = target_evaluator()
      print "x,f:", ",".join(["%6.3f"%x_ for x_ in x]), f, icall
      have_request = minimizer.process(x, f, g)
      if(have_request):
        requests_f_and_g = minimizer.requests_f_and_g()
        continue
      assert not minimizer.requests_f_and_g()
      if(minimizer.is_terminated()): break
  except RuntimeError, e:
    minimizer.error = str(e)
  minimizer.n_calls = icall
  return minimizer

class minimizer:

  def __init__(self,
               potential,
               use_bounds):
    adopt_init_args(self, locals())
    self.x = self.potential.ncs_pairs[0].rho_mn
    self.n = self.x.size()

  def run(self, use_curvatures=0):
    self.minimizer = lbfgs_run(
      target_evaluator=self,
      use_bounds=self.use_bounds)
    self()
    return self

  def __call__(self):
    self.potential.update(x = self.x)
    self.f = self.potential.target()
    self.g = self.potential.gradient()
    return self.x, self.f, self.g

class potential(object):

  def __init__(self, f_obs, ncs_pairs, reflections_per_bin = 250):
    adopt_init_args(self, locals())
    # Create bins and SigmaN
    f_obs.setup_binner(reflections_per_bin = reflections_per_bin)
    binner = f_obs.binner()
    n_bins = binner.n_bins_used()
    self.SigmaN = flex.double(f_obs.data().size(), -1)
    for i_bin in binner.range_used():
      bin_sel = f_obs.binner().selection(i_bin)
      f_obs_bin = f_obs.select(bin_sel)
      f_obs_bin_data = f_obs_bin.data()
      eps_bin = f_obs_bin.epsilons().data().as_double()
      sn = flex.sum(f_obs_bin_data*f_obs_bin_data/eps_bin)/f_obs_bin_data.size()
      self.SigmaN = self.SigmaN.set_selected(bin_sel, sn)
      print "bin: %d n_refl.: %d" % (i_bin, f_obs_bin.data().size()), \
        "%6.3f-%-6.3f"%f_obs_bin.d_max_min(), "SigmaN: %6.4f"%sn
    assert self.SigmaN.all_gt(0)
    #
    self.rbin = flex.int(f_obs.data().size(), -1)
    for i_bin in binner.range_used():
      for i_seq in binner.array_indices(i_bin):
        self.rbin[i_seq] = i_bin-1 # i_bin starts with 1, not 0 !
    assert flex.min(self.rbin)==0
    assert flex.max(self.rbin)==n_bins-1
    # Extract symmetry matrices
    self.sym_matrices = []
    for m_as_string in f_obs.space_group().smx():
      o = sgtbx.rt_mx(symbol=str(m_as_string), t_den=f_obs.space_group().t_den())
      m_as_double = o.r().as_double()
      print m_as_string, m_as_double
      self.sym_matrices.append(m_as_double)

  def update(self, x):
    assert x.size() == self.ncs_pairs[0].rho_mn.size() #XXX
    self.target_and_grads = ext.tncs_eps_factor_refinery(
      tncs_pairs               = self.ncs_pairs,
      f_obs                    = self.f_obs.data(),
      sigma_f_obs              = self.f_obs.sigmas(),
      rbin                     = self.rbin,
      SigmaN                   = self.SigmaN,
      space_group              = self.f_obs.space_group(),
      miller_indices           = self.f_obs.indices(),
      fractionalization_matrix = self.f_obs.unit_cell().fractionalization_matrix(),
      sym_matrices             = self.sym_matrices)

  def target(self):
    return self.target_and_grads.target()

  def gradient(self):
    return self.target_and_grads.gradient()
