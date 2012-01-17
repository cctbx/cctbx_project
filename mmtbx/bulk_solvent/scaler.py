from scitbx.array_family import flex
import sys
from mmtbx import bulk_solvent
from cctbx import adptbx
from libtbx.math_utils import iround
from libtbx import adopt_init_args
import boost.python
ext = boost.python.import_ext("mmtbx_f_model_ext")
from cctbx import sgtbx
from mmtbx.bulk_solvent import bulk_solvent_and_scaling
import scitbx.math
import mmtbx.f_model

def sliding_window_average(x, offset):
  result = flex.double(x.size(), -1)
  for i in xrange(x.size()):
    s = 0
    cntr = 0
    for j in range(max(0,i-offset), min(x.size()-1, i+offset)):
      s+=x[j]
      cntr+=1
    result[i]=s/cntr
  return result

class core(object):
  def __init__(self,
               f_obs,
               f_calc,
               f_mask,
               scalar_scale,
               overall_scale,
               overall_anisotropic_scale,
               bulk_solvent_scale,
               ss):
    adopt_init_args(self, locals())
    assert f_obs.indices().all_eq(f_calc.indices())
    assert f_obs.indices().all_eq(f_mask.indices())
    assert f_obs.indices().size() == ss.size()
    assert overall_scale.size() == overall_anisotropic_scale.size()
    assert overall_scale.size() == bulk_solvent_scale.size()
    self.core_data = ext.core(
      f_calc                    = self.f_calc.data(),
      f_mask                    = self.f_mask.data(),
      scale                     = 1,
      overall_scale             = self.overall_scale,
      overall_anisotropic_scale = self.overall_anisotropic_scale,
      bulk_solvent_scale        = self.bulk_solvent_scale)

  def select(self, selection):
    assert self.f_obs.indices().size() == selection.size()
    return core(
      f_obs              = self.f_obs.select(selection),
      f_calc             = self.f_calc.select(selection),
      f_mask             = self.f_mask.select(selection),
      scalar_scale       = 1,
      overall_scale      = self.overall_scale.select(selection),
      bulk_solvent_scale = self.bulk_solvent_scale.select(selection),
      overall_anisotropic_scale = self.overall_anisotropic_scale.select(selection),
      ss                 = self.ss.select(selection))

  def f_model(self):
    return self.f_calc.customized_copy(data=self.core_data.f_model)

  def f_model_no_scale(self):
    return self.f_calc.customized_copy(data =
      self.core_data.f_model_no_aniso_scale)

  def r_factor(self):
    return bulk_solvent.r_factor(self.f_obs.data(), self.f_model().data())

  def update(self, overall_scale=None, bulk_solvent_scale=None,
             overall_anisotropic_scale=None):
    if(overall_scale is not None):
      assert overall_scale.size() == self.overall_scale.size()
      self.overall_scale = overall_scale
    if(bulk_solvent_scale is not None):
      assert bulk_solvent_scale.size() == self.bulk_solvent_scale.size()
      self.bulk_solvent_scale = bulk_solvent_scale
    if(overall_anisotropic_scale is not None):
      assert overall_anisotropic_scale.size() == \
        self.overall_anisotropic_scale.size()
      self.overall_anisotropic_scale = overall_anisotropic_scale
    self.core_data = ext.core(
      f_calc                    = self.f_calc.data(),
      f_mask                    = self.f_mask.data(),
      scale                     = 1,
      overall_scale             = self.overall_scale,
      overall_anisotropic_scale = self.overall_anisotropic_scale,
      bulk_solvent_scale        = self.bulk_solvent_scale)

  def try_bulk_solvent_scale(self, overall_scale=None, bulk_solvent_scale=None):
    core_data = ext.core(
      f_calc                    = self.f_calc.data(),
      f_mask                    = self.f_mask.data(),
      scale                     = 1,
      overall_scale             = overall_scale,
      overall_anisotropic_scale = self.overall_anisotropic_scale,
      bulk_solvent_scale        = bulk_solvent_scale)
    return bulk_solvent.r_factor(self.f_obs.data(), core_data.f_model)

  def try_overall_anisotropic_scale(self, scale_array):
    core_data = ext.core(
      f_calc                    = self.f_calc.data(),
      f_mask                    = self.f_mask.data(),
      scale                     = 1,
      overall_scale             = self.overall_scale,
      overall_anisotropic_scale = scale_array,
      bulk_solvent_scale        = self.bulk_solvent_scale)
    return bulk_solvent.r_factor(self.f_obs.data(), core_data.f_model)

def determine_n_bins(
  n_refl,
  n_refl_per_bin,
  max_n_bins,
  min_refl_per_bin,
  min_n_bins=1):
    assert n_refl_per_bin > 0
    n_refl_per_bin = min(n_refl, iround(n_refl_per_bin))
    result = max(1, iround(n_refl / max(1, n_refl_per_bin)))
    result = max(result, min(min_n_bins, iround(n_refl / min_refl_per_bin)))
    result = min(max_n_bins, result)
    return result

def split(mi, r1=5, r2=10, m1=500, m2=100, m3=10):
  d_spacings = mi.d_spacings().data()
  s1 = d_spacings <= r1
  s2 = (d_spacings > r1) & (d_spacings <= r2)
  s3 = d_spacings > r2
  c1 = s1.count(True)
  c2 = s2.count(True)
  c3 = s3.count(True)
  #
  if(c1>=m1 and c2>=m2 and c3>=m3):
    s1 = d_spacings <= r1
    s2 = (d_spacings > r1) & (d_spacings <= r2)
    s3 = d_spacings > r2
  elif(c1<=m1 and c2>=m2 and c3>=m3):
    s1 = None
    s2 = d_spacings <= r2
    s3 = d_spacings > r2
  elif(c1>=m1 and c2>=m2 and c3<=m3):
    s1 = d_spacings <= r1
    s2 = d_spacings > r1
    s3 = None
  elif(c1>=m1 and c2<=m2):
    s1 = flex.bool(mi.data().size(),True)
    s2 = None
    s3 = None
  elif(c1<=m1 and c2<=m2):
    s1 = None
    s2 = None
    s3 = flex.bool(mi.data().size(),True)
  else:
    assert 0 # should never get here (if it does - it's a bug)
  nb1,nb2,nb3 = None,None,None
  if(s1 is not None):
    nb1 = determine_n_bins(
      n_refl           = s1.count(True),
      n_refl_per_bin   = m1,
      max_n_bins       = 100,
      min_refl_per_bin = m1)
  if(s2 is not None):
    nb2 = determine_n_bins(
      n_refl           = s2.count(True),
      n_refl_per_bin   = m2,
      max_n_bins       = 50,
      min_refl_per_bin = m2)
  if(s3 is not None):
    nb3 = determine_n_bins(
      n_refl           = s3.count(True),
      n_refl_per_bin   = m3,
      max_n_bins       = 1000, #XXX
      min_refl_per_bin = m3)
  cntr = 0
  for s in [s1,s2,s3]:
    if(s is not None):
      cntr += s.count(True)
  assert cntr == mi.data().size()
  return ((s1,nb1), (s2,nb2), (s3,nb3))

class run(object):
  def __init__(self,
               f_obs,
               f_calc, # can be a sum: f_calc=f_hydrogens+f_calc+f_part
               f_mask, # only one shell is supported
               ss,
               number_of_cycles=10, # termination occures much earlier
               estimate_k_sol_and_b_sol=False, # not used, for information only
               ANISO_FIRST=False, # XXX TESTING ONLY
               log=None,
               try_poly = True,
               try_expanal = True,
               try_expmin = False,
               verbose=False):
    self.log = log
    self.verbose = verbose
    self.ss = ss
    self.try_poly    = try_poly
    self.try_expanal = try_expanal
    self.try_expmin  = try_expmin
    if(log is None): log = sys.stdout
    if(verbose):
      print >> log, \
        "Overall, iso- and anisotropic scaling and bulk-solvent modeling:"
    bins_and_selections = split(mi=f_obs, r1=5, r2=10, m1=500, m2=100, m3=10)
    point_group = sgtbx.space_group_info(
      symbol=f_obs.space_group().type().lookup_symbol()
      ).group().build_derived_point_group()
    self.adp_constraints = sgtbx.tensor_rank_2_constraints(
      space_group=point_group,
      reciprocal_space=True)
    scalar_scale = bulk_solvent.scale(f_obs.data(), f_calc.data())
    if(verbose):
      print >> log,\
        "  k_overall (k_bs=0,k_anisotropic=1,k_isotropic=1): %6.4f"%scalar_scale
    self.core = core(
      f_obs              = f_obs,
      f_calc             = f_calc,
      f_mask             = f_mask,
      scalar_scale       = 1,
      overall_scale      = flex.double(f_obs.size(), 1),
      overall_anisotropic_scale = flex.double(f_obs.size(), 1),
      bulk_solvent_scale = flex.double(f_obs.size(), 0),
      ss                 = ss)
    self.cores_and_selections = []
    if(verbose):
      print >> log, "  Binning:"
    for i_bas, bas in enumerate(bins_and_selections):
      if(bas[0] is not None):
        core_selected = self.core.select(selection=bas[0])
        f_obs_data = core_selected.f_obs.data()
        core_selected.f_obs.setup_binner_counting_sorted(
          reflections_per_bin = int(f_obs_data.size()/bas[1]))
        #core_selected.f_obs.setup_binner(reflections_per_bin=0, n_bins = bas[1])
        fmt="    resolution: %6.2f-%-6.2f bins:%3d reflections:%4d"
        dmm = core_selected.f_obs.d_max_min()
        if(verbose):
          print >> log, fmt%(dmm[1],dmm[0],bas[1],f_obs_data.size())
        bin_selections = []
        n_outl = 0
        for i_bin in core_selected.f_obs.binner().range_used():
          sel = core_selected.f_obs.binner().selection(i_bin)
          f_obs_data_i = f_obs_data.select(sel)
          fodim = flex.mean(f_obs_data_i)
          sel_ = (f_obs_data < fodim*6) & (f_obs_data > fodim/6)
          sel_use = sel & sel_
          bin_selections.append([sel.iselection(), sel_use])
        self.cores_and_selections.append([bas[0],core_selected,bin_selections])
    #
    if(verbose):
      print >> log, "  r_start: %6.4f"%self.core.r_factor()
    for cycle in xrange(number_of_cycles):
      r_start = self.core.r_factor()
      if(verbose):
        print >> log, "  cycle %d"%cycle
        print >> log, "    r: %6.4f"%(r_start)
      # anisotropic scaling
      if(ANISO_FIRST):
        if(verbose):
          print >> log, "    anisotropic scaling:"
        self.anisotropic_scaling(r_start = r_start) # order IS important
      # bulk-solvent and overall isotropic scale
      if(cycle==0): # XXX
        self.bulk_solvent_simple(r_start=r_start)
        if(verbose):
          print >> self.log, "bulk_solvent_simple: ", self.core.r_factor()
      else:
        ssi,x = self.bulk_solvent_scaling(r_start = r_start)
      if(not ANISO_FIRST):
        if(verbose):
          print >> log, "    anisotropic scaling:"
        self.anisotropic_scaling(r_start = r_start) # XXX order IS important
      self.r_final = self.core.r_factor()
      if((r_start<=self.r_final) or
         (r_start>self.r_final and abs(r_start-self.r_final)<1.e-4)):
        break
    if(verbose):
      print >> log, "r-factor (final): %6.4f"%(self.core.r_factor())
    if(estimate_k_sol_and_b_sol):
      for a,b,c,d,e,f in zip(ssi,x,
                             sliding_window_average(x=x,offset=1),
                             sliding_window_average(x=x,offset=2),
                             sliding_window_average(x=x,offset=3),
                             sliding_window_average(x=x,offset=4)):
        print >> log, "%8.6f %8.6f %8.6f %8.6f %8.6f %8.6f"%(a,b, c,d,e,f)
      add = False
      for i in [0]:#[0,1,2,3,4]:
        if i == 0:
          sel = x>0
          xi = x.select(sel)
          ssii = ssi.select(sel)
        else:
          xi = sliding_window_average(x=x,offset=i)
          sel = xi>0
          xi = xi.select(sel)
          ssii = ssi.select(sel)
        r = scitbx.math.gaussian_fit_1d_analytical(x = ssii, y = xi)
        print >> log, "k_sol, b_sol: %5.3f %7.2f"%(r.a, r.b)
    #
    #try:
    #  d = f_obs.d_spacings().data()
    #  sel = (d>8) & (self.core.bulk_solvent_scale>0)
    #  ssii = self.core.ss.select(sel)
    #  xi = self.core.bulk_solvent_scale.select(sel)
    #  r = scitbx.math.gaussian_fit_1d_analytical(x = flex.sqrt(ssii), y = xi)
    #  print r.a, r.b
    #  bs = r.a*flex.exp(-ssii*r.b)
    #  nbs = self.core.bulk_solvent_scale.set_selected(sel, bs)
    #
    #  core_data = ext.core(
    #    f_calc                    = self.core.f_calc.data(),
    #    f_mask                    = self.core.f_mask.data(),
    #    scale                     = 1,
    #    overall_scale             = self.core.overall_scale.set_selected(sel, 1),
    #    overall_anisotropic_scale = self.core.overall_anisotropic_scale,
    #    bulk_solvent_scale        = nbs)
    #  scale_k1 = bulk_solvent.scale(
    #    f_obs.data().select(sel), core_data.f_model.select(sel))
    #  nos = self.core.overall_scale.set_selected(sel, scale_k1)
    #
    #  r = self.core.try_bulk_solvent_scale(
    #    overall_scale      = nos,
    #    bulk_solvent_scale = nbs)
    #  print "LOOK2:", r
    #except:
    #  print "FAILED"

  def bulk_solvent_simple(self, r_start):
    scale_k1 = bulk_solvent.scale(self.core.f_obs.data(), self.core.f_model().data())
    overall_scale = flex.double(self.ss.size(), scale_k1)
    bulk_solvent_scale = bulk_solvent.ksol_bsol_grid_search(
      self.core.f_obs.data(),
      self.core.f_calc.data(),
      self.core.f_mask.data(),
      flex.double([0.1,0.3,0.5]),
      flex.double([30,60,90]),
      self.ss,
      1,
      overall_scale,
      self.core.overall_anisotropic_scale,
      r_start)
    self.core.update(
      bulk_solvent_scale = bulk_solvent_scale,
      overall_scale = overall_scale)

  def bulk_solvent_scaling(self, r_start):
    overall_scale = flex.double(self.core.f_obs.size(), -1)
    bulk_solvent_scale = flex.double(self.core.f_obs.size(), -1)
    ssi = flex.double()
    x = flex.double()
    for i_cas, cas in enumerate(self.cores_and_selections):
      sel = cas[0]
      scale = self.core.overall_anisotropic_scale.select(sel)
      a,b,c,d = self.scale(core=cas[1], scale=scale, bin_selections=cas[2])
      if(i_cas == 0):
        ssi,x = a,b
        overall_scale.set_selected(sel, c)
        bulk_solvent_scale.set_selected(sel, d)
      else:
        ssi.extend(a)
        x.extend(b)
        overall_scale.set_selected(sel, c)
        bulk_solvent_scale.set_selected(sel, d)
    assert (overall_scale < 0).count(True) == 0
    assert (bulk_solvent_scale < 0).count(True) == 0
    r = self.core.try_bulk_solvent_scale(
      overall_scale      = overall_scale,
      bulk_solvent_scale = bulk_solvent_scale)
    suffix = ""
    if(r>r_start): suffix = "(result rejected due to r-factor increase)"
    if(self.verbose):
      print >> self.log, "    bulk-solvent:"
      print >> self.log, "      r        : %6.4f %s"%(r, suffix)
    if(r<r_start):
      self.core.update(
        overall_scale      = overall_scale,
        bulk_solvent_scale = bulk_solvent_scale)
    return ssi,x

  def anisotropic_scaling(self, r_start):
    scale_k1 = bulk_solvent.scale(self.core.f_obs.data(),
      self.core.f_model_no_scale().data())
    r_expanal, r_poly, r_expmin = None,None,None
    overall_anisotropic_scale_expanal, overall_anisotropic_scale_poly, \
      overall_anisotropic_scale_expmin = None, None, None
    scale_matrix_expanal, scale_matrix_poly, scale_matrix_expmin= None,None,None
    # try exp_anal
    if(self.try_expanal):
      obj = bulk_solvent.aniso_u_scaler(
        f_model        = self.core.f_model_no_scale().data()*scale_k1,
        f_obs          = self.core.f_obs.data(),
        miller_indices = self.core.f_obs.indices(),
        adp_constraint_matrix = self.adp_constraints.gradient_sum_matrix())
      u_star = self.adp_constraints.all_params(tuple(obj.u_star_independent))
      scale_matrix_expanal = adptbx.u_as_b(adptbx.u_star_as_u_cart(
        self.core.f_obs.unit_cell(), u_star))
      overall_anisotropic_scale_expanal = ext.overall_anisotropic_scale(
        self.core.f_obs.indices(), u_star)
      r_expanal = self.core.try_overall_anisotropic_scale(
        scale_array = overall_anisotropic_scale_expanal)
      if(self.verbose):
        print >> self.log, "      r_expanal: %6.4f"%r_expanal
    # try poly
    if(self.try_poly):
      obj = bulk_solvent.aniso_u_scaler(
        f_model        = self.core.f_model_no_scale().data()*scale_k1,
        f_obs          = self.core.f_obs.data(),
        miller_indices = self.core.f_obs.indices(),
        unit_cell      = self.core.f_obs.unit_cell())
      scale_matrix_poly = obj.a
      overall_anisotropic_scale_poly = ext.overall_anisotropic_scale(
        self.core.f_obs.indices(), obj.a, self.core.f_obs.unit_cell())
      r_poly = self.core.try_overall_anisotropic_scale(
        scale_array = overall_anisotropic_scale_poly)
      if(self.verbose):
        print >> self.log, "      r_poly   : %6.4f"%r_poly
    # try expmin
    if(self.try_expmin):
      zero = self.core.f_calc.customized_copy(data =
        flex.complex_double(self.core.f_obs.data().size(), 0))
      fm = mmtbx.f_model.core(
        f_calc  = self.core.f_model_no_scale(),
        f_mask  = zero,
        k_sols  = [0.],
        b_sol   = 0.,
        f_part1 = None,
        f_part2 = None,
        u_star  = [0,0,0,0,0,0],
        fmodel  = None,
        ss      = None)
      obj = bulk_solvent_and_scaling.u_star_minimizer(
        fmodel_core_data = fm,
        f_obs            = self.core.f_obs,
        u_initial        = [0,0,0,0,0,0],
        refine_u         = True,
        min_iterations   = 500,
        max_iterations   = 500,
        symmetry_constraints_on_b_cart = True,
        u_min_max = 500.,
        u_min_min =-500.)
      u_star = obj.u_min
      scale_matrix_expmin = adptbx.u_as_b(adptbx.u_star_as_u_cart(
        self.core.f_obs.unit_cell(), u_star))
      overall_anisotropic_scale_expmin = ext.overall_anisotropic_scale(
        self.core.f_obs.indices(), u_star)
      r_expmin = self.core.try_overall_anisotropic_scale(
        scale_array = overall_anisotropic_scale_expmin)
    # select best
    r = [(r_expanal, overall_anisotropic_scale_expanal, scale_matrix_expanal),
         (r_poly,    overall_anisotropic_scale_poly,    scale_matrix_poly),
         (r_expmin,  overall_anisotropic_scale_expmin,  scale_matrix_expmin)]
    r_best = r_start
    overall_aniso_scale_best = None
    scale_matrix_best = None
    for result in r:
      r_factor, overall_aniso_scale, scale_matrix = result
      if(r_factor is not None and r_factor < r_best):
        r_best = r_factor
        overall_aniso_scale_best = overall_aniso_scale.deep_copy()
        scale_matrix_best = scale_matrix[:]
    if(scale_matrix_best is None):
      if(self.verbose):
        print >> self.log, "      result rejected due to r-factor increase"
    else:
      self.core.update(overall_anisotropic_scale = overall_aniso_scale_best)
      r_aniso = self.core.r_factor()
      if(self.verbose):
        print >> self.log, "      r_final  : %6.4f"%r_aniso
        if(len(scale_matrix_best)<=6):
          print >> self.log, "      b_cart(11,22,33,12,13,23):",\
            ",".join([str("%8.4f"%i).strip() for i in scale_matrix_best])
        else:
          print >> self.log, "      a:",\
            ",".join([str("%8.4f"%i).strip() for i in scale_matrix_best])

  def scale(self, core, scale, bin_selections):
    f_obs  = core.f_obs
    f_calc = core.f_calc.customized_copy(data = core.f_calc.data() * scale)
    f_mask = core.f_mask.customized_copy(data = core.f_mask.data() * scale)
    ss     = core.ss
    overall_scale = flex.double(f_obs.size(), -1)
    bulk_solvent_scale = flex.double(f_obs.size(), -1)
    ssi = flex.double()
    x = flex.double()
    sels = []
    for sel in bin_selections:
      sel, sel_use = sel
      ss_ = flex.mean(flex.sqrt(ss.select(sel)))
      obj = bulk_solvent.overall_and_bulk_solvent_scale_coefficients_analytical(
        f_obs     = f_obs.data(),
        f_calc    = f_calc.data(),
        f_mask    = f_mask.data(),
        selection = sel_use)
      ssi.append(ss_)
      x.append(obj.x)
      overall_scale=overall_scale.set_selected(sel, obj.t)
      bulk_solvent_scale=bulk_solvent_scale.set_selected(sel, obj.x)
    assert (overall_scale < 0).count(True) == 0
    assert (bulk_solvent_scale < 0).count(True) == 0
    return ssi, x, overall_scale, bulk_solvent_scale
