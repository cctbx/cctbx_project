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
from scitbx.math import matrix
import sys

class groups(object):

  def __init__(self,
               pdb_hierarchy,
               crystal_symmetry,
               n_bins, # XXX remove later
               angular_difference_threshold_deg=10,
               rad=None):
    pdb_hierarchy = pdb_hierarchy.expand_to_p1(
      crystal_symmetry=crystal_symmetry)
    pdb_hierarchy.atoms().reset_i_seq()
    sites_cart = pdb_hierarchy.atoms().extract_xyz()
    ncs_inp = iotbx.ncs.input(hierarchy=pdb_hierarchy)
    ncs_groups = ncs_inp.get_ncs_restraints_group_list()
    self.rta = []
    for g in ncs_groups:
      for c in g.copies:
        angle = c.r.rotation_angle()
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
        radius_estimate=rad,
        fracscat=1./(2*len(self.rta)),
        rho_mn=flex.double(n_bins,0.98))
      self.ncs_pairs.append(ncs_pair)

def lbfgs_run(target_evaluator, use_bounds, lower_bound, upper_bound):
  minimizer = lbfgsb.minimizer(
    n   = target_evaluator.n,
    l   = lower_bound, # lower bound
    u   = upper_bound, # upper bound
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

  def __init__(self, f_obs, ncs_pairs, reflections_per_bin):
    adopt_init_args(self, locals())
    # Create bins
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
    ### Initialize rho_mn
    ### rhoMN = exp(-(2*pi^2/3)*(rms/d)^2, and rms=0.4-0.8 is probably a good.
    rho_mn_initial = flex.double(n_bins, 0)
    d_spacings = self.f_obs.d_spacings().data()
    cntr=0
    for i_bin in self.binner.range_used():
      sel_bin = self.binner.selection(i_bin)
      if(sel_bin.count(True)>0):
        arg = (2*math.pi**2/3)*(0.5/flex.mean(d_spacings.select(sel_bin)))**2
        rho_mn_initial[cntr] = math.exp(-1*arg)
      cntr+=1
    ncs_pairs[0].set_rhoMN(rho_mn_initial)
    ###
    self.update()

  def update(self, x=None):
    if(self.gradient_evaluator=="rhoMN"):
      size = len(self.ncs_pairs)
      for i, ncs_pair in enumerate(self.ncs_pairs):
        ncs_pair.set_rhoMN(x[i:i+x.size()//size])
    elif(self.gradient_evaluator=="radius"):
      for ncs_pair, x_ in zip(self.ncs_pairs, x):
        ncs_pair.set_radius(x_)
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
    self.SigmaN = flex.double(self.f_obs.data().size(), 0)
    for i_bin in self.binner.range_used():
      bin_sel = self.f_obs.binner().selection(i_bin)
      f_obs_bin = self.f_obs.select(bin_sel)
      f_obs_bin_data = f_obs_bin.data()
      f_obs_bin_data_size = f_obs_bin_data.size()
      if(f_obs_bin_data_size>0):
        eps_bin = eps.select(bin_sel)
        sn = flex.sum(f_obs_bin_data*f_obs_bin_data/eps_bin)/f_obs_bin_data_size
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

class compute_eps_factor(object):

  def __init__(self, f_obs, pdb_hierarchy, reflections_per_bin):
    f_obs = f_obs.deep_copy()
    if(not f_obs.sigmas_are_sensible()):
      f_obs.set_sigmas(sigmas = flex.double(f_obs.data().size(), 0.0))
    reflections_per_bin = min(f_obs.data().size(), reflections_per_bin)
    f_obs.setup_binner(reflections_per_bin = reflections_per_bin)
    binner = f_obs.binner()
    n_bins = binner.n_bins_used()
    self.unit_cell = f_obs.unit_cell()
    #
    self.ncs_pairs = groups(
      pdb_hierarchy    = pdb_hierarchy,
      crystal_symmetry = f_obs.crystal_symmetry(),
      n_bins           = n_bins).ncs_pairs
    self.epsfac = None
    if(len(self.ncs_pairs)>0):
      # Radii
      radii = flex.double()
      rad_lower_bound = flex.double()
      rad_upper_bound = flex.double()
      for ncs_pair in self.ncs_pairs:
        radii.append(ncs_pair.radius)
        rad_lower_bound.append(ncs_pair.radius/3)
        rad_upper_bound.append(ncs_pair.radius*3)
      # Target and gradients evaluator
      pot = potential(f_obs = f_obs, ncs_pairs = self.ncs_pairs,
        reflections_per_bin = reflections_per_bin)
      for it in xrange(10):
        # refine eps fac
        rho_mn = flex.double()
        for ncs_pair in self.ncs_pairs:
          rho_mn.extend(ncs_pair.rho_mn)
        m = minimizer(
          potential      = pot.set_refine_rhoMN(),
          use_bounds     = 2,
          lower_bound    = flex.double(rho_mn.size(), 0.),
          upper_bound    = flex.double(rho_mn.size(), 1.),
          initial_values = rho_mn).run()
        # refine radius
        radii = flex.double()
        for ncs_pair in self.ncs_pairs:
          radii.append(ncs_pair.radius)
        m = minimizer(
          potential      = pot.set_refine_radius(),
          use_bounds     = 2,
          lower_bound    = rad_lower_bound,
          upper_bound    = rad_upper_bound,
          initial_values = radii).run()
      self.epsfac = pot.target_and_grads.tncs_epsfac()

  def show_summary(self, log=None):
    if(self.epsfac is None): return None
    if(log is None): log = sys.stdout
    for i, ncs_pair in enumerate(self.ncs_pairs):
      print >> log, "tNCS group: %d"%i
      angle = matrix.sqr(ncs_pair.r).rotation_angle()
      t = ",".join([("%6.3f"%t_).strip() for t_ in ncs_pair.t]).strip()
      t_cart = ",".join([("%6.3f"%t_).strip()
        for t_ in self.unit_cell.orthogonalize(ncs_pair.t)]).strip()
      r = ",".join([("%8.6f"%r_).strip() for r_ in ncs_pair.r]).strip()
      print >> log, "  Translation (fractional): (%s)"%t
      print >> log, "  Translation (Cartesian):  (%s)"%t_cart
      print >> log, "  Rotation (deg): %-5.2f"%angle
      print >> log, "  Rotation matrix: (%s)"%r
      print >> log, "  Radius: %-6.1f"%ncs_pair.radius
      print >> log, "  fracscat:", ncs_pair.fracscat
    print >> log, "tNCS eps factor: min,max,mean: %6.4f %6.4f %6.4f"%\
      self.epsfac.min_max_mean().as_tuple()
