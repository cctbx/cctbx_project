from scitbx.array_family import flex
import sys, math, time
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
from scitbx.math import curve_fitting

def set_bin_selections(d_spacings):
  d_spacings = d_spacings.data()
  selections = []
  cntr = 0
  # this approximately corresponds to regular splitting of s between 0 and 1.
  limits = [(0   , 1),
            (1   , 1.25),
            (1.25, 1.5),
            (1.5 , 1.75),
            (1.75, 2),
            (2  , 2.5),
            (2.5, 3),
            (3  , 3.5),
            (3.5, 4),
            (4  , 5),
            (5  , 6),
            (6  , 7),
            (7  , 8),
            (8  , 9),
            (9  , 10),
            (10 , 11),
            (11 , 12),
            (12 , 13),
            (13 , 14),
            (14 , 15),
            (15 , 20),
            (20 , 25),
            (25 , 30),
            (30 , 35),
            (35 , 40),
            (40 , 45),
            (45 , 50),
            (50 , 999)]
  limits.reverse()
  for s in limits:
    sel  = d_spacings >= s[0]
    sel &= d_spacings <  s[1]
    #if(sel.count(True)>0):
    #  print s, sel.count(True), flex.min(d_spacings.select(sel)),flex.max(d_spacings.select(sel))
    cntr += sel.count(True)
    if(sel.count(True)>0): selections.append(sel)
  assert cntr == d_spacings.size()
  for i in xrange(len(selections)):
    count_true = selections[i].count(True)
    if(count_true>0 and count_true<100 and i+1<len(selections)):
      s = selections[i] | selections[i+1]
      selections[i+1]=s
      selections[i]=None
    elif(count_true<500 and i==len(selections)-1 and selections[i-1] is not None):
      #print selections[i] , selections[i-1]
      s = selections[i] | selections[i-1]
      selections[i-1]=s
      selections[i]=None
  #print
  sn = []
  cntr = 0
  for s in selections:
    if(s is not None):
      sn.append(s)
      cntr += s.count(True)
      #print s.count(True)
  assert cntr == d_spacings.size()
  return sn

class core(object):
  def __init__(self,
               f_obs,
               f_calc,
               f_mask,
               scalar_scale,
               k_isotropic,
               k_anisotropic,
               k_mask,
               ss):
    adopt_init_args(self, locals())
    assert f_obs.indices().all_eq(f_calc.indices())
    assert f_obs.indices().all_eq(f_mask.indices())
    assert f_obs.indices().size() == ss.size()
    assert k_isotropic.size() == k_anisotropic.size()
    assert k_isotropic.size() == k_mask.size()
    self.core_data = ext.core(
      f_calc        = self.f_calc.data(),
      f_mask        = self.f_mask.data(),
      scale         = 1,
      k_isotropic   = self.k_isotropic,
      k_anisotropic = self.k_anisotropic,
      k_mask        = self.k_mask)

  def select(self, selection):
    assert self.f_obs.indices().size() == selection.size()
    return core(
      f_obs         = self.f_obs.select(selection),
      f_calc        = self.f_calc.select(selection),
      f_mask        = self.f_mask.select(selection),
      scalar_scale  = 1,
      k_isotropic   = self.k_isotropic.select(selection),
      k_mask        = self.k_mask.select(selection),
      k_anisotropic = self.k_anisotropic.select(selection),
      ss            = self.ss.select(selection))

  def f_model(self):
    return self.f_calc.customized_copy(data=self.core_data.f_model)

  def f_model_no_scale(self):
    return self.f_calc.customized_copy(data =
      self.core_data.f_model_no_aniso_scale)

  def r_factor(self):
    return bulk_solvent.r_factor(self.f_obs.data(), self.f_model().data())

  def update(self, k_isotropic=None, k_mask=None,
             k_anisotropic=None):
    if(k_isotropic is not None):
      assert k_isotropic.size() == self.k_isotropic.size()
      self.k_isotropic = k_isotropic
    if(k_mask is not None):
      assert k_mask.size() == self.k_mask.size()
      self.k_mask = k_mask
    if(k_anisotropic is not None):
      assert k_anisotropic.size() == \
        self.k_anisotropic.size()
      self.k_anisotropic = k_anisotropic
    self.core_data = ext.core(
      f_calc        = self.f_calc.data(),
      f_mask        = self.f_mask.data(),
      scale         = 1,
      k_isotropic   = self.k_isotropic,
      k_anisotropic = self.k_anisotropic,
      k_mask        = self.k_mask)

  def try_scale(self, k_isotropic=None, k_mask=None, k_anisotropic=None):
    if(k_isotropic is None): k_isotropic = self.k_isotropic
    if(k_mask is None): k_mask = self.k_mask
    if(k_anisotropic is None):
      k_anisotropic = self.k_anisotropic
    core_data = ext.core(
      f_calc        = self.f_calc.data(),
      f_mask        = self.f_mask.data(),
      scale         = 1,
      k_isotropic   = k_isotropic,
      k_anisotropic = k_anisotropic,
      k_mask        = k_mask)
    return bulk_solvent.r_factor(self.f_obs.data(), core_data.f_model)

class run(object):
  def __init__(self,
               f_obs,
               f_calc, # can be a sum: f_calc=f_hydrogens+f_calc+f_part
               f_mask, # only one shell is supported
               ss,
               number_of_cycles=20, # termination occures much earlier
               use_polynomial_fit=True,
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
    self.d_spacings = f_obs.d_spacings()
    self.polynomial_fit_params = None
    if(log is None): log = sys.stdout
    if(verbose):
      print >> log, "-"*80
      print >> log, \
        "Overall, iso- and anisotropic scaling and bulk-solvent modeling:"
    point_group = sgtbx.space_group_info(
      symbol=f_obs.space_group().type().lookup_symbol()
      ).group().build_derived_point_group()
    self.adp_constraints = sgtbx.tensor_rank_2_constraints(
      space_group=point_group,
      reciprocal_space=True)
    self.core = core(
      f_obs         = f_obs,
      f_calc        = f_calc,
      f_mask        = f_mask,
      scalar_scale  = 1,
      k_isotropic   = flex.double(f_obs.size(), 1),
      k_anisotropic = flex.double(f_obs.size(), 1),
      k_mask        = flex.double(f_obs.size(), 0),
      ss            = ss)
    self.cores_and_selections = []
    self.bin_selections = set_bin_selections(d_spacings = self.d_spacings)
    if(verbose):
      print >> log, "  Using %d resolution bins"%len(self.bin_selections)
    self.ss_bin_average=flex.double()
    for i_sel, sel in enumerate(self.bin_selections):
      core_selected = self.core.select(selection=sel)
      f_obs_data = core_selected.f_obs.data()
      fodim = flex.mean(f_obs_data)
      sel_use  = (f_obs_data < fodim*6)
      sel_use &= (f_obs_data > fodim/6)
      self.cores_and_selections.append([sel, core_selected, sel_use])
      self.ss_bin_average.append(flex.mean(core_selected.ss)) # unsquare, take mena, then square back ?
    fit_accepted = False
    for cycle in xrange(number_of_cycles):
      r_start = self.core.r_factor()
      if(verbose):
        print >> log, "  cycle %d:"%cycle
        print >> log, "    r(start): %6.4f"%(r_start)
      # bulk-solvent and overall isotropic scale
      k_mask_bin,k_overall_bin,k_mask,k_isotropic = self.bulk_solvent_scaling()
      bs_accepted = self.accept_scale(k_mask=k_mask, k_isotropic=k_isotropic,
        r_ref=r_start, prefix="r(bulk solvent)")
      if(use_polynomial_fit):
        k_mask, k_isotropic, self.polynomial_fit_params = self.fit_poly_k_mask(
          y=k_mask_bin)
        fit_accepted = self.accept_scale(k_mask=k_mask, k_isotropic=k_isotropic,
          r_ref=self.core.r_factor(), prefix="r(bulk solvent, polynomial approx.)")
      # overall scale
      self.apply_overall_scale()
      # anisotropic scale
      if(verbose):
        print >> log, "    anisotropic scaling:"
      self.anisotropic_scaling(r_start = r_start) # XXX order IS important
      if(self.is_converged(r_start=r_start)): break
    self.apply_overall_scale()
    if(verbose):
      print >> log, "  r(final): %6.4f"%(self.core.r_factor())
    if(verbose):
      self.show(k_mask=k_mask_bin, fit_accepted=fit_accepted)
      print >> log, "-"*80

  def apply_overall_scale(self):
    self.core.update(k_isotropic = self.core.k_isotropic*bulk_solvent.scale(
      self.core.f_obs.data(), self.core.f_model().data()))

  def is_converged(self, r_start, tolerance=1.e-4):
    self.r_final = self.core.r_factor()
    result = False
    if((r_start<=self.r_final) or
       (r_start>self.r_final and abs(r_start-self.r_final)<tolerance)):
      result = True
    return result

  def show(self, k_mask, fit_accepted):
    if(fit_accepted):
      a = self.polynomial_fit_params
      a = []
      for i, p in enumerate(self.polynomial_fit_params):
        if(p>0):
          if(i==0): s=""
          else: s="+"
        else: s="-"
        p = str("%12.4f"%abs(p)).strip()
        p = "%s %s"%(s,p)
        a.append(p)
      print >> self.log, "  Polynomial approximation for k_mask(s):"
      print >> self.log, "    k_mask(ssq) = %s %s*ssq^1 %s*ssq^2 %s*ssq^3"%(a[0],
        a[1],a[2],a[3])
    else:
      print >> self.log, "  Polynomial approximation for k_mask was not used."
    print >> self.log, "  Statistics in resolution bins:"
    assert k_mask.size() == len(self.bin_selections)
    fmt="  %7.5f %6.2f -%6.2f %5.1f %5d %6.3f %6.3f %6.3f %6.3f %8.2f %6.4f"
    f_model = self.core.f_model().data()
    print >> self.log, "  s^2      Resolution    Compl Nrefl  k_mask        k_iso  k_ani <Fobs>   R"
    print >> self.log, "                (A)        (%)        raw    poly"
    for i_sel, sel in enumerate(self.bin_selections):
      ss_                  = self.ss_bin_average[i_sel]
      k_mask_              = k_mask[i_sel]
      k_maks_bin_averaged_ = flex.mean(self.core.k_mask.select(sel))
      d_                   = self.d_spacings.data().select(sel)
      d_min_               = flex.min(d_)
      d_max_               = flex.max(d_)
      n_ref_               = d_.size()
      f_obs_               = self.core.f_obs.select(sel)
      f_obs_mean_          = flex.mean(f_obs_.data())
      k_isotropic_         = flex.mean(self.core.k_isotropic.select(sel))
      k_anisotropic_       = flex.mean(self.core.k_anisotropic.select(sel))
      cmpl_                = f_obs_.completeness(d_max=d_max_)*100.
      r_          = bulk_solvent.r_factor(f_obs_.data(), f_model.select(sel), 1)
      print >> self.log, fmt%(ss_, d_max_, d_min_, cmpl_, n_ref_, k_mask_,
        k_maks_bin_averaged_, k_isotropic_, k_anisotropic_, f_obs_mean_, r_)

  def fit_poly_k_mask(self, y):
    X = self.ss_bin_average
    Y = y
    for i in xrange(len(Y)):
      if(i!=0 and i!=len(Y)-1):
        if(Y[i]<=0 and Y[i-1]>0 and Y[i+1]>0):
          Y[i]=(Y[i-1]>0 + Y[i+1])/2.
    Y_ = flex.double()
    X_ = flex.double()
    for tmpY, tmpX in zip(Y,X):
      if(tmpY <= 0):
        X_.append(tmpX)
        Y_.append(0)
        break
      X_.append(tmpX)
      Y_.append(tmpY)
    cutoff_x = tmpX
    cfo = curve_fitting.univariate_polynomial_fit(x_obs=X_, y_obs=Y_, degree=3,
      max_iterations=1000, min_iterations=1000)
    k_mask = bulk_solvent.set_k_mask_to_cubic_polynom(self.core.ss, cutoff_x,
      cfo.params)
    k_isotropic = flex.double(self.core.ss.size(), -1)
    core_data = ext.core(
      f_calc        = self.core.f_calc.data(),
      f_mask        = self.core.f_mask.data(),
      scale         = 1,
      k_isotropic   = flex.double(self.core.ss.size(), 1),
      k_anisotropic = self.core.k_anisotropic,
      k_mask        = k_mask)
    for sel in self.bin_selections:
      scale_k1 = bulk_solvent.scale(
        self.core.f_obs.data().select(sel), core_data.f_model.select(sel))
      k_isotropic = k_isotropic.set_selected(sel, scale_k1)
    assert k_isotropic.count(-1.) == 0
    assert (k_mask<0).count(True)==0
    return k_mask, k_isotropic, cfo.params

  def accept_scale(self, r_ref, k_isotropic=None, k_mask=None,
                   k_anisotropic=None, prefix=""):
    assert [k_isotropic, k_mask, k_anisotropic].count(None)>0
    r = self.core.try_scale(
      k_isotropic   = k_isotropic,
      k_mask        = k_mask,
      k_anisotropic = k_anisotropic)
    better = r<r_ref
    suffix = ""
    if(not better): suffix = "(result rejected)"
    if(self.verbose):
      print >> self.log, "    %s: %6.4f %s"%(prefix, r, suffix)
    if(better):
      self.core.update(
        k_isotropic   = k_isotropic,
        k_mask        = k_mask,
        k_anisotropic = k_anisotropic)
    return better

  def bulk_solvent_scaling(self):
    k_isotropic       = flex.double(self.core.f_obs.size(), -1)
    kmask             = flex.double(self.core.f_obs.size(), -1)
    k_isotropic_bin   = flex.double()
    kmask_bin         = flex.double()
    for i_cas, cas in enumerate(self.cores_and_selections):
      selection, core, selection_use = cas
      scale = self.core.k_anisotropic.select(selection)
      obj = bulk_solvent.overall_and_bulk_solvent_scale_coefficients_analytical(
        f_obs     = core.f_obs.data(),
        f_calc    = core.f_calc.data() * scale,
        f_mask    = core.f_mask.data() * scale,
        selection = selection_use)
      kmask_bin.append(obj.x)
      k_isotropic_bin.append(obj.t)
      k_isotropic.set_selected(selection, obj.t)
      kmask.set_selected(selection, obj.x)
    assert (k_isotropic < 0).count(True) == 0
    assert (kmask < 0).count(True) == 0
    return kmask_bin,k_isotropic_bin, kmask, k_isotropic

  def anisotropic_scaling(self, r_start):
    r_expanal, r_poly, r_expmin = None,None,None
    k_anisotropic_expanal, k_anisotropic_poly, \
      k_anisotropic_expmin = None, None, None
    scale_matrix_expanal, scale_matrix_poly, scale_matrix_expmin= None,None,None
    f_model_no_scale_data = self.core.f_model_no_scale().data()
    # try exp_anal
    if(self.try_expanal):
      obj = bulk_solvent.aniso_u_scaler(
        f_model        = f_model_no_scale_data,
        f_obs          = self.core.f_obs.data(),
        miller_indices = self.core.f_obs.indices(),
        adp_constraint_matrix = self.adp_constraints.gradient_sum_matrix())
      u_star = self.adp_constraints.all_params(tuple(obj.u_star_independent))
      scale_matrix_expanal = adptbx.u_as_b(adptbx.u_star_as_u_cart(
        self.core.f_obs.unit_cell(), u_star))
      k_anisotropic_expanal = ext.overall_anisotropic_scale(
        self.core.f_obs.indices(), u_star)
      r_expanal = self.core.try_scale(
        k_anisotropic = k_anisotropic_expanal)
      if(self.verbose):
        print >> self.log, "      r_expanal: %6.4f"%r_expanal
    # try poly
    if(self.try_poly):
      obj = bulk_solvent.aniso_u_scaler(
        f_model        = f_model_no_scale_data,
        f_obs          = self.core.f_obs.data(),
        miller_indices = self.core.f_obs.indices(),
        unit_cell      = self.core.f_obs.unit_cell())
      scale_matrix_poly = obj.a
      k_anisotropic_poly = ext.overall_anisotropic_scale(
        self.core.f_obs.indices(), obj.a, self.core.f_obs.unit_cell())
      r_poly = self.core.try_scale(
        k_anisotropic = k_anisotropic_poly)
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
        symmetry_constraints_on_b_cart = False,
        u_min_max = 500.,
        u_min_min =-500.)
      u_star = obj.u_min
      scale_matrix_expmin = adptbx.u_as_b(adptbx.u_star_as_u_cart(
        self.core.f_obs.unit_cell(), u_star))
      k_anisotropic_expmin = ext.overall_anisotropic_scale(
        self.core.f_obs.indices(), u_star)
      r_expmin = self.core.try_scale(
        k_anisotropic = k_anisotropic_expmin)
      if(self.verbose):
        print >> self.log, "    r_expmin   : %6.4f"%r_expmin
    # select best
    r = [(r_expanal, k_anisotropic_expanal, scale_matrix_expanal),
         (r_poly,    k_anisotropic_poly,    scale_matrix_poly),
         (r_expmin,  k_anisotropic_expmin,  scale_matrix_expmin)]
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
      self.core.update(k_anisotropic = overall_aniso_scale_best)
      r_aniso = self.core.r_factor()
      if(self.verbose):
        print >> self.log, "      r_final  : %6.4f"%r_aniso
        if(len(scale_matrix_best)<=6):
          print >> self.log, "      b_cart(11,22,33,12,13,23):",\
            ",".join([str("%8.4f"%i).strip() for i in scale_matrix_best])
        else:
          v0=[]
          v1=[]
          for i, a in enumerate(scale_matrix_best):
            if(i in [0,2,4,6,8,10]): v1.append(a)
            else: v0.append(a)
          print >> self.log, "      V0:",\
            ",".join([str("%8.4f"%i).strip() for i in v0])
          print >> self.log, "      V1:",\
            ",".join([str("%8.4f"%i).strip() for i in v1])
