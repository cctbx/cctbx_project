from __future__ import division
import cctbx.array_family.flex # import dependency
import boost.python
ext = boost.python.import_ext("mmtbx_ncs_ext")
from scitbx.array_family import flex
from cctbx import sgtbx
from libtbx import adopt_init_args
from scitbx import lbfgsb
import iotbx.ncs
import math
import scitbx.math

class groups(object):

  def __init__(self,
               pdb_hierarchy,
               crystal_symmetry,
               n_bins, # XXX remove later
               angular_difference_threshold_deg=10,
               rad=None):
    #asc = pdb_hierarchy.atom_selection_cache()
    #sel = asc.selection("pepnames and (name CA or name N or name O or name C)")
    #pdb_hierarchy = pdb_hierarchy.select(sel)
    pdb_hierarchy = pdb_hierarchy.expand_to_p1(
      crystal_symmetry=crystal_symmetry)
    pdb_hierarchy.atoms().reset_i_seq()
    sites_cart = pdb_hierarchy.atoms().extract_xyz()
    ncs_inp = iotbx.ncs.input(hierarchy=pdb_hierarchy)
    ncs_groups = ncs_inp.get_ncs_restraints_group_list()
    self.rta = []
    for g in ncs_groups:
      for c in g.copies:
        angle = math.acos((c.r.trace()-1)/2)*180./math.pi
        r = c.r
        t = crystal_symmetry.unit_cell().fractionalize(c.t)
        if(angle < angular_difference_threshold_deg):
          t_new = []
          for t_ in t:
            t_new.append(math.modf(t_)[0]) # same as t_-int(t_)
          # radius
          if(rad is None):
            sites_cart_sel = sites_cart.select(c.iselection)
            radius = 0
            for sc in sites_cart_sel - sites_cart_sel.mean():
              x,y,z = sc
              radius += math.sqrt(x**2+y**2+z**2)
            radius = radius/sites_cart_sel.size()*4./3.
          else:
            radius = rad
          #
          self.rta.append([c.r, t_new, angle, radius])
    self.ncs_pairs = []
    for i, it in enumerate(self.rta):
      r,t,a, rad = it
      ncs_pair = ext.pair(
        r = r,
        t = t,
        radius=rad,
        fracscat=1./(2*len(self.rta)),
        rho_mn=flex.double(n_bins,0.98) )
      self.ncs_pairs.append(ncs_pair)

  def show_summary(self):
    for i, it in enumerate(self.rta):
      print "tNCS group: %d"%i
      r,t,a, rad = it
      t = ",".join([("%6.3f"%t_).strip() for t_ in t]).strip()
      r = ",".join([("%8.6f"%r_).strip() for r_ in r]).strip()
      print "  Translation (fractional): (%s)"%t
      print "  Rotation (deg): %-5.2f"%a
      print "  Rotation matrix: (%s)"%r
      print "  Radius: %-6.1f"%rad
      print "  fracscat:", self.ncs_pairs[i].fracscat

def lbfgs_run(target_evaluator, use_bounds, lower_bound, upper_bound):
  minimizer = lbfgsb.minimizer(
    n   = target_evaluator.n,
    l   = flex.double(target_evaluator.n, lower_bound), # lower bound
    u   = flex.double(target_evaluator.n, upper_bound), # upper bound
    nbd = flex.int(target_evaluator.n, use_bounds)) # flag to apply both bounds
  minimizer.error = None
  try:
    icall = 0
    while 1:
      icall += 1
      x, f, g = target_evaluator()
      #print "x,f:", ",".join(["%6.3f"%x_ for x_ in x]), f, icall
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

class minimizer(object):

  def __init__(self,
               potential,
               use_bounds,
               lower_bound,
               upper_bound,
               initial_values):
    adopt_init_args(self, locals())
    self.x = initial_values
    self.n = self.x.size()

  def run(self):
    self.minimizer = lbfgs_run(
      target_evaluator=self,
      use_bounds=self.use_bounds,
      lower_bound = self.lower_bound,
      upper_bound = self.upper_bound)
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
    self.binner = f_obs.binner()
    n_bins = self.binner.n_bins_used()
    self.SigmaN = None
    self.update_SigmaN()
    #
    self.rbin = flex.int(f_obs.data().size(), -1)
    for i_bin in self.binner.range_used():
      for i_seq in self.binner.array_indices(i_bin):
        self.rbin[i_seq] = i_bin-1 # i_bin starts with 1, not 0 !
    assert flex.min(self.rbin)==0
    assert flex.max(self.rbin)==n_bins-1
    # Extract symmetry matrices
    self.sym_matrices = []
    for m_as_string in f_obs.space_group().smx():
      o = sgtbx.rt_mx(symbol=str(m_as_string), t_den=f_obs.space_group().t_den())
      m_as_double = o.r().as_double()
      self.sym_matrices.append(m_as_double)
    self.gradient_evaluator = None
    self.update()

  def update(self, x=None):
    if(self.gradient_evaluator=="rhoMN"):
      self.ncs_pairs[0].set_rhoMN(x)
    elif(self.gradient_evaluator=="radius"):
      self.ncs_pairs[0].set_radius(x[0])
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

  def update_SigmaN(self):
    if(self.SigmaN is None):
      eps = self.f_obs.epsilons().data().as_double()
    else:
      eps = self.target_and_grads.tncs_epsfac()
    self.SigmaN = flex.double(self.f_obs.data().size(), -1)
    for i_bin in self.binner.range_used():
      bin_sel = self.f_obs.binner().selection(i_bin)
      f_obs_bin = self.f_obs.select(bin_sel)
      f_obs_bin_data = f_obs_bin.data()
      eps_bin = eps.select(bin_sel)
      sn = flex.sum(f_obs_bin_data*f_obs_bin_data/eps_bin)/f_obs_bin_data.size()
      self.SigmaN = self.SigmaN.set_selected(bin_sel, sn)
    assert self.SigmaN.all_gt(0)

  def set_refine_radius(self):
    self.gradient_evaluator = "radius"
    return self

  def set_refine_rhoMN(self):
    self.gradient_evaluator = "rhoMN"
    return self

  def target(self):
    return self.target_and_grads.target()

  def gradient(self):
    if(self.gradient_evaluator=="rhoMN"):
      return self.target_and_grads.gradient_rhoMN()
    elif(self.gradient_evaluator=="radius"):
      return self.target_and_grads.gradient_radius()
    else: assert 0
