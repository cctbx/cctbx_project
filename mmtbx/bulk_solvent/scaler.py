from scitbx.array_family import flex
import sys
from mmtbx import bulk_solvent
from cctbx import adptbx
from libtbx import adopt_init_args
import boost.python
ext = boost.python.import_ext("mmtbx_f_model_ext")
from cctbx import sgtbx
from mmtbx.bulk_solvent import bulk_solvent_and_scaling
import mmtbx.f_model
from scitbx.math import curve_fitting
import math, time
from libtbx import group_args
import scitbx.math
from libtbx.test_utils import approx_equal
from cctbx import miller

import mmtbx.arrays

def moving_average(x, anchor_left, anchor_right):
  x_ = [anchor_left] + list(x) + [anchor_right]
  result = flex.double()
  for i in xrange(len(x_)):
    if(i!=0 and i!=len(x_)-1 and i!=1):
      result.append( (x_[i-1]+x_[i]+x_[i+1])/3. )
    elif(i==1): result.append(x_[i])
  return result

def moving_average0(x):
  try:
    result = flex.double(x.size(), -1)
    for i in xrange(x.size()):
      if(i==0): result[i]=(x[i]+x[i+1])/2.
      elif(i==x.size()-1): result[i]=(x[i-1]+x[i])/2.
      else: result[i]=(x[i-1]+x[i]+x[i+1])/3.
    assert (result<0).count(True)==0
  except: return x
  return result



def moving_average2(x):
  x_ = x_ = [x[0]] + list(x) + [x[len(x)-1]]

  #result = x_[:]
  for cycle in xrange(5):
    result = x_[:]
    #
    selection = flex.bool(len(result), False)
    for i, s in enumerate(selection):
      if(i!=0 and i!=len(result)-1):
        if((result[i-1]<result[i] and result[i+1]<result[i]) or
           (result[i-1]>result[i] and result[i+1]>result[i])):
          selection[i]=True
    #
    for i in xrange(len(result)):
      if(i!=0 and i!=len(result)-1 and selection[i]):
        result[i] = (x_[i-1]+x_[i]+x_[i+1])/3.
    x_ = result[:]

  return result[1:len(result)-1]


def set_bin_selections(d_spacings, min_ref_low):
  d_spacings = d_spacings.data()
  selections = []
  cntr = 0
  # this approximately corresponds to regular splitting of s between 0 and 1.
  l1 = [0]
  l2 = [float("%6.2f"%(i/100.)) for i in range(100, 400,    25)]
  l3 = [float("%6.2f"%(i/100.)) for i in range(400, 1000,   50)]
  l4 = [float("%6.2f"%(i/100.)) for i in range(1000,1500,   100)]
  l5 = [float("%6.2f"%(i/100.)) for i in range(1500,5000,   500)]
  l6 = [float("%6.2f"%(i/100.)) for i in range(5000,10001, 1000)]
  l7 = [999]
  l = l1+l2+l3+l4+l5+l6+l7
  limits = []
  for i_l, l_ in enumerate(l):
    if(i_l != len(l)-1): limits.append( [l_,l[i_l+1]] )
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
    if(count_true>0 and count_true<min_ref_low and i+1<len(selections)):
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


def gs(fo, fc, fm):
  def target(fo, fc, fm, x, scale):
    fmodel = flex.abs(scale * (fc + x * fm))
    diff = (fmodel*fmodel-(fo*fo))
    return flex.sum(diff*diff)
  k_mask_trial_range = flex.double([i/1000. for i in range(0,1000,10)])
  t_best = target(fo=fo, fc=fc, fm=fm, x=0, scale=1)
  x_best = 0
  for x in k_mask_trial_range:
    t = target(fo=fo, fc=fc, fm=fm, x=x, scale=1)
    if(t<t_best):
      t_best = t
      x_best = x
  return x_best

class run(object):
  def __init__(self,
               f_obs,
               f_calc, # can be a sum: f_calc=f_hydrogens+f_calc+f_part
               f_mask, # only one shell is supported
               r_free_flags,
               ss,
               k_mask_approximator,
               scale_method,
               smooth_k_mask=0,
               number_of_cycles=20, # termination occures much earlier
               log=None,
               try_poly = True,
               try_expanal = True,
               try_expmin = False,
               tmp1=True,tmp2=True,
               verbose=False):
    self.log = log
    self.scale_method = scale_method
    self.verbose = verbose
    self.r_free_flags = r_free_flags
    self.ss = ss
    self.try_poly    = try_poly
    self.try_expanal = try_expanal
    self.try_expmin  = try_expmin
    self.d_spacings = f_obs.d_spacings()
    self.k_mask_approximator = k_mask_approximator
    self.r_low = None
    self.smooth_k_mask = smooth_k_mask
    self.poly_approx_cutoff = None
    def init_result():
      return group_args(
        k_mask_bin_orig   = None,
        k_mask_bin_smooth = None,
        k_mask            = None,
        k_isotropic       = None,
        k_mask_fit_params = None)
    self.bss_result = init_result()
    assert k_mask_approximator in ["poly3", None]
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
    self.core = mmtbx.arrays.init(
      f_obs         = f_obs,
      f_calc        = f_calc,
      f_mask        = f_mask,
      r_free_flags  = r_free_flags,
      k_isotropic   = flex.double(f_obs.size(), 1),
      k_anisotropic = flex.double(f_obs.size(), 1),
      k_mask        = flex.double(f_obs.size(), 0))
    self.cores_and_selections = []
    self.bin_selections = self._bin_selections()
    self.low_res_selection = self._low_resolution_selection()
    if(verbose):
      print >> log, "  Using %d resolution bins"%len(self.bin_selections)
    self.ss_bin_values=[]
    for i_sel, sel in enumerate(self.bin_selections):
      core_selected = self.core.select(selection=sel)
      f_obs_data = core_selected.f_obs.data()
      fodim = flex.mean(f_obs_data)
      sel_use  = (f_obs_data < fodim*6)
      sel_use &= (f_obs_data > fodim/6)
      self.cores_and_selections.append([sel, core_selected, sel_use])
      self.ss_bin_values.append([
        flex.min(core_selected.ss),
        flex.max(core_selected.ss),
        flex.mean(core_selected.ss)])
    fit_accepted = False
    for cycle in xrange(number_of_cycles):
      r_start = self.r_factor()
      r_start0 = r_start
      if(verbose):
        print >> log, "  cycle %d:"%cycle
        print >> log, "    r(start): %6.4f"%(r_start)
      # bulk-solvent and overall isotropic scale
      if(cycle==0):
        r_start = self.bulk_solvent_simple(r_start=r_start)
        if(verbose):
          print >> self.log, "    r(bulk_solvent_grid_search): %6.4f"%r_start
      else:
        r_start = self.bulk_solvent_scaling(r_start = r_start)
      # anisotropic scale
      if([try_poly, try_expanal, try_expmin].count(True)):
        if(verbose): print >> log, "    anisotropic scaling:"
        self.anisotropic_scaling(r_start = r_start)
      if(self.is_converged(r_start=r_start0)): break
    self.apply_overall_scale()
    if(verbose):
      print >> log, "  r(final): %6.4f"%(self.r_factor())
      self.show(fit_accepted=fit_accepted)
    #
    self.r_low = self._r_low()
    if(verbose):
      d = self.d_spacings.data().select(self.low_res_selection)
      d1 = ("%7.4f"%flex.min(d)).strip()
      d2 = ("%7.4f"%flex.max(d)).strip()
      n = d.size()
      print >> self.log, "r(low-resolution: %s-%s A; %d reflections): %6.4f"%(
        d2,d1,n,self.r_low)
      print >> log, "-"*80
    self.r_final = self.r_factor()

  def try_scale(self,
                k_isotropic=None,
                k_mask=None,
                k_anisotropic=None,
                selection=None):
    if(k_isotropic is None):   k_isotropic   = self.core.k_isotropic
    if(k_mask is None):        k_mask        = self.core.k_mask
    if(k_anisotropic is None): k_anisotropic = self.core.k_anisotropic
    c = mmtbx.arrays.init(
      f_obs        = self.core.f_obs,
      f_calc       = self.core.f_calc,
      f_mask       = self.core.f_mask,
      r_free_flags = self.core.r_free_flags,
      k_isotropic  = k_isotropic,
      k_anisotropic= k_anisotropic,
      k_mask       = k_mask)
    if(selection is not None):
      return bulk_solvent.r_factor(c.f_obs.data(), c.f_model.data(), selection)
    else:
      return bulk_solvent.r_factor(c.f_obs.data(), c.f_model.data())

  def r_factor(self):
    return bulk_solvent.r_factor(self.core.f_obs.data(),
      self.core.f_model.data())#, self.core.selection_work)

  def _low_resolution_selection(self):
    d = self.d_spacings.data()
    if((d>=8).count(True)>=500):
      return d>=8
    else:
      sel = self.bin_selections[0].deep_copy()
      if(sel.count(True)>=500): return sel
      else:
        for i in range(1, len(self.bin_selections)): #XXX
          sel |= self.bin_selections[i].deep_copy()
          if(sel.count(True)>=300): break
        return sel

  def _r_low(self):
    return bulk_solvent.r_factor(self.core.f_obs.data(),
      self.core.f_model.data(), self.low_res_selection) ### exclude test set

  def _bin_selections(self):
    r = set_bin_selections(d_spacings = self.d_spacings, min_ref_low = 300)
    d_sp_0 = self.d_spacings.data().select(r[0])
    if(flex.min(d_sp_0)<6 and d_sp_0.size() > 100):
      r = set_bin_selections(d_spacings = self.d_spacings, min_ref_low = 50)
    elif(flex.min(d_sp_0)<8 and d_sp_0.size() > 100):
      r = set_bin_selections(d_spacings = self.d_spacings, min_ref_low = 100)
    return r

  def populate_bin_to_individual_k_mask_linear_interpolation(self, k_mask_bin):
    assert len(k_mask_bin) == len(self.cores_and_selections)
    def linear_interpolation(x1,x2,y1,y2):
      k=0
      if(x1!=x2): k=(y2-y1)/(x2-x1)
      b = y1-k*x1
      return k,b
    result = flex.double(self.core.f_obs.size(), -1)
    for i, cas in enumerate(self.cores_and_selections):
      selection, zzz, zzz = cas
      x1,x2 = self.ss_bin_values[i][0], self.ss_bin_values[i][1]
      y1 = k_mask_bin[i]
      if(i==len(k_mask_bin)-1):
        y2 = k_mask_bin[i-1]
      else:
        y2 = k_mask_bin[i+1]
      k,b = linear_interpolation(x1=x1,x2=x2,y1=y1,y2=y2)
      bulk_solvent.set_to_liear_interpolated(self.core.ss,k,b,selection,result)
    assert (result < 0).count(True) == 0
    return result

  def bulk_solvent_simple(self, r_start):
    k_mask_trial_range = flex.double([i/1000. for i in range(0,650,50)])
    k_mask      = flex.double(self.core.f_obs.size(), -1)
    k_isotropic = flex.double(self.core.f_obs.size(), -1)
    k_mask_bin = flex.double()
    for i_cas, cas in enumerate(self.cores_and_selections):
      selection, core, selection_use = cas
      scale = self.core.k_anisotropic.select(selection)
      k_mask_bin_, k_isotropic_bin_ = \
          bulk_solvent.k_mask_and_k_overall_grid_search(
            core.f_obs.data()/scale,
            core.f_calc.data(),
            core.f_mask.data(),
            k_mask_trial_range)
      k_mask_bin.append(k_mask_bin_)
      k_mask.set_selected(selection, k_mask_bin_)
      fm = core.k_anisotropic*(core.f_calc.data()+k_mask_bin_*core.f_mask.data())
      scale_k1 = bulk_solvent.scale(core.f_obs.data(), fm)
      k_isotropic.set_selected(selection, scale_k1)
    k_mask_bin_smooth = self.smooth(k_mask_bin)
    k_mask = self.populate_bin_to_individual_k_mask_linear_interpolation(
      k_mask_bin = k_mask_bin_smooth)
    k_isotropic = self._k_isotropic_as_scale_k1(k_mask = k_mask)
    self.core = self.core.update(k_mask = k_mask, k_isotropic = k_isotropic)
    self.bss_result.k_mask_bin_orig = k_mask_bin
    self.bss_result.k_mask_bin_smooth = k_mask_bin_smooth
    self.bss_result.k_mask            = k_mask
    self.bss_result.k_isotropic       = k_isotropic
    return self.r_factor()

  def apply_overall_scale(self):
    self.core = self.core.update(k_isotropic = self.core.k_isotropic*bulk_solvent.scale(
      self.core.f_obs.data(), self.core.f_model.data()))

  def is_converged(self, r_start, tolerance=1.e-4):
    self.r_final = self.r_factor()
    result = False
    if((r_start<=self.r_final) or
       (r_start>self.r_final and abs(r_start-self.r_final)<tolerance)):
      result = True
    return result

  def poly3(self, params, x):
    a = params
    return a[0] + a[1] * x**1 + a[2] * x**2 + a[3] * x**3

  def show(self, fit_accepted):
    b = self.bss_result
    kmfp = None
    if(b is not None and b.k_mask_fit_params is not None):
      kmfp = b.k_mask_fit_params

    if(kmfp is not None):
      print >> self.log, "k_mask approximation:",self.k_mask_approximator,\
      [str("%10.5f"%i).strip() for i in kmfp]
    else:
      print >> self.log, "  Polynomial approximation for k_mask was not used."
    print >> self.log, "  Statistics in resolution bins:"
    #assert k_mask.size() == len(self.bin_selections)
    fmt="  %7.5f %6.2f -%6.2f %5.1f %5d %-6s %-6s %-6s %-6s %6.3f %6.3f %8.2f %6.4f"
    f_model = self.core.f_model.data()
    print >> self.log, "  s^2      Resolution    Compl Nrefl k_mask                       k_iso  k_ani <Fobs>   R"
    print >> self.log, "                (A)        (%)       orig   smooth ave    approx"
    k_mask_bin_orig_   = str(None)
    k_mask_bin_smooth_ = str(None)
    k_mask_bin_approx_ = str(None)
    for i_sel, sel in enumerate(self.bin_selections):
      ss_ = self.ss_bin_values[i_sel][2]
      if(b is not None and self.bss_result.k_mask_bin_orig is not None):
        k_mask_bin_orig_ = "%6.4f"%self.bss_result.k_mask_bin_orig[i_sel]
      if(b is not None and self.bss_result.k_mask_bin_smooth is not None):
        k_mask_bin_smooth_ = "%6.4f"%self.bss_result.k_mask_bin_smooth[i_sel]
      k_maks_bin_averaged_ = "%6.4f"%flex.mean(self.core.k_mask.select(sel))
      if(kmfp is not None):
        if(ss_< self.poly_approx_cutoff):
          k_mask_bin_approx_ = "%6.4f"%self.poly3(params=kmfp, x=ss_)
        else: k_mask_bin_approx_ = "%6.4f"%0
      d_             = self.d_spacings.data().select(sel)
      d_min_         = flex.min(d_)
      d_max_         = flex.max(d_)
      n_ref_         = d_.size()
      f_obs_         = self.core.f_obs.select(sel)
      f_obs_mean_    = flex.mean(f_obs_.data())
      k_isotropic_   = flex.mean(self.core.k_isotropic.select(sel))
      k_anisotropic_ = flex.mean(self.core.k_anisotropic.select(sel))
      cmpl_          = f_obs_.completeness(d_max=d_max_)*100.
      r_             = bulk_solvent.r_factor(f_obs_.data(), f_model.select(sel), 1)
      print >> self.log, fmt%(ss_, d_max_, d_min_, cmpl_, n_ref_,
        k_mask_bin_orig_, k_mask_bin_smooth_,k_maks_bin_averaged_,k_mask_bin_approx_,
        k_isotropic_, k_anisotropic_, f_obs_mean_, r_)

# DISABLE  def fit_poly_k_mask(self, y):
# DISABLE    X = []
# DISABLE    for x in self.ss_bin_values:
# DISABLE      X.append(x[2])
# DISABLE    Y = y
# DISABLE    for i in xrange(len(Y)):
# DISABLE      if(i!=0 and i!=len(Y)-1):
# DISABLE        if(Y[i]<=0 and Y[i-1]>0 and Y[i+1]>0):
# DISABLE          Y[i]=(Y[i-1]>0 + Y[i+1])/2.
# DISABLE    Y_ = flex.double()
# DISABLE    X_ = flex.double()
# DISABLE    for tmpY, tmpX in zip(Y,X):
# DISABLE      d = 1/(2.*math.sqrt(tmpX)) # XXX
# DISABLE      if(tmpY <= 0 or d<3.): # XXX
# DISABLE        X_.append(tmpX)
# DISABLE        Y_.append(0)
# DISABLE        break
# DISABLE      X_.append(tmpX)
# DISABLE      Y_.append(tmpY)
# DISABLE    self.poly_approx_cutoff = tmpX
# DISABLE    cfo = curve_fitting.univariate_polynomial_fit(x_obs=X_, y_obs=Y_, degree=3,
# DISABLE        max_iterations=30, min_iterations=25, number_of_cycles=25)
# DISABLE    k_mask = bulk_solvent.set_k_mask_to_cubic_polynom(self.core.ss,
# DISABLE      self.poly_approx_cutoff, cfo.params)
# DISABLE    k_isotropic = self._k_isotropic_as_scale_k1(k_mask = k_mask)
# DISABLE    assert (k_mask<0).count(True)==0
# DISABLE    return k_mask, k_isotropic, cfo.params

  def _k_isotropic_as_scale_k1(self, k_mask):
    result = flex.double(self.core.ss.size(), -1)
    core_data = ext.core(
      f_calc        = self.core.f_calc.data(),
      f_mask        = self.core.f_mask.data(),
      k_isotropic   = flex.double(self.core.ss.size(), 1),
      k_anisotropic = self.core.k_anisotropic,
      k_mask        = k_mask)
    for sel in self.bin_selections:
      scale_k1 =bulk_solvent.scale(self.core.f_obs.data(),core_data.f_model,sel)
      result = result.set_selected(sel, scale_k1)
    assert result.count(-1.) == 0
    return result

  def accept_scale(self, r_ref, k_isotropic=None, k_mask=None,
                   k_anisotropic=None, prefix=""):
    assert [k_isotropic, k_mask, k_anisotropic].count(None)>0
    r = self.try_scale(
      k_isotropic   = k_isotropic,
      k_mask        = k_mask,
      k_anisotropic = k_anisotropic)
    better = r<r_ref
    suffix = ""
    if(not better): suffix = "(result rejected)"
    if(self.verbose):
      print >> self.log, "    %s: %6.4f %s"%(prefix, r, suffix)
    if(better):
      self.core = self.core.update(
        k_isotropic   = k_isotropic,
        k_mask        = k_mask,
        k_anisotropic = k_anisotropic)
    return better

  def estimate_scale_k1(self, cutoff=4, width=1, min_reflections=500):
    cutoff = min(cutoff, self.core.f_obs.d_min()+width)
    sel_high = self.d_spacings.data()<cutoff
    scale_k1 = 1
    if(sel_high.count(True)>min_reflections):
      core = self.core.select(selection = sel_high)
      fm = core.k_anisotropic * (core.f_calc.data() +
        core.k_mask * core.f_mask.data())
      scale_k1 = bulk_solvent.scale(core.f_obs.data(), fm)
    return scale_k1

  def bulk_solvent_scaling(self, r_start):
    k_mask     = flex.double(self.core.f_obs.size(), -1)
    k_mask_bin = flex.double()
    k_mask_trial_range = flex.double([i/1000. for i in range(0,1000,10)])
    scale_k1 = self.estimate_scale_k1()
    def get_k_mask_trial_range(x, shift=0.05):
      result = flex.double()
      inc = max(0,x-shift)
      while inc<=x+shift+1.e-3:
        result.append(inc)
        inc+=0.01
      return result
    for i_cas, cas in enumerate(self.cores_and_selections):
      selection, core, selection_use = cas
      scale = self.core.k_anisotropic.select(selection)*scale_k1
      f_obs  = core.f_obs.data()
      f_calc = core.f_calc.data()*scale
      f_mask = core.f_mask.data()*scale
      if(self.scale_method == "k_iso_k_mask_anal"):
        obj = bulk_solvent.overall_and_bulk_solvent_scale_coefficients_analytical(
          f_obs     = f_obs,
          f_calc    = f_calc,
          f_mask    = f_mask,
          selection = flex.bool(selection_use.size(), True))
        k_mask_bin.append(obj.x_best)
        k_mask.set_selected(selection, obj.x_best)
      elif(self.scale_method == "k_mask_anal"):
        obj = bulk_solvent.bulk_solvent_scale_coefficients_analytical(
          f_obs     = f_obs,
          f_calc    = f_calc,
          f_mask    = f_mask,
          selection = flex.bool(selection_use.size(), True))
        k_mask_bin.append(obj.x_best)
        k_mask.set_selected(selection, obj.x_best)
      elif(self.scale_method == "k_mask_r_grid_search"):
        k_mask_bin_, k_isotropic_bin_ = \
          bulk_solvent.k_mask_and_k_overall_grid_search(
            f_obs,
            f_calc,
            f_mask,
            k_mask_trial_range)
        k_mask_bin.append(k_mask_bin_)
        k_mask.set_selected(selection, k_mask_bin_)
      elif(self.scale_method == "combo"):
        obj1 = bulk_solvent.overall_and_bulk_solvent_scale_coefficients_analytical(
          f_obs     = f_obs,
          f_calc    = f_calc,
          f_mask    = f_mask,
          selection = flex.bool(selection_use.size(), True)) # use selection_use: never a good idea
        k_mask.set_selected(selection, obj1.x_best)
        r1 = self.try_scale(k_mask = k_mask, selection=selection)
        #
        obj2 = bulk_solvent.bulk_solvent_scale_coefficients_analytical(
          f_obs     = f_obs,
          f_calc    = f_calc,
          f_mask    = f_mask,
          selection = flex.bool(selection_use.size(), True)) # use selection_use: never a good idea
        k_mask.set_selected(selection, obj2.x_best)
        r2 = self.try_scale(k_mask = k_mask, selection=selection)
        if(r1<r2): x = obj1.x_best
        else: x = obj2.x_best
        # fine-sample k_mask around minimum of LS to fall into minimum of R
        k_mask_bin_, k_isotropic_bin_ = \
          bulk_solvent.k_mask_and_k_overall_grid_search(
            f_obs,
            f_calc,
            f_mask,
            get_k_mask_trial_range(x = x))
        k_mask_bin.append(k_mask_bin_)
        k_mask.set_selected(selection, k_mask_bin_)
      else: assert 0
    #
    k_mask_bin_smooth = self.smooth(k_mask_bin)
    k_mask = self.populate_bin_to_individual_k_mask_linear_interpolation(
      k_mask_bin = k_mask_bin_smooth)
    k_isotropic = self._k_isotropic_as_scale_k1(k_mask = k_mask)
    #
    r = self.try_scale(k_mask = k_mask, k_isotropic = k_isotropic)
    if(r<=r_start or (r>r_start and abs(r-r_start)*100<0.5)): # may be 0.5?
      self.core = self.core.update(k_isotropic = k_isotropic, k_mask = k_mask)
      self.bss_result.k_mask_bin_orig   = k_mask_bin
      self.bss_result.k_mask_bin_smooth = k_mask_bin_smooth
      self.bss_result.k_mask            = k_mask
      self.bss_result.k_isotropic       = k_isotropic
      r_start = r
      if(self.verbose):
        print >> self.log, "    %s: %6.4f"%("bulk-solvent", r_start)
      #
# Note, the result of approximation was not actually used ! Fix it later.
#
# DISABLE       k_mask_fit_parameters = None
# DISABLE       k_mask, k_isotropic, k_mask_fit_parameters = self.fit_poly_k_mask(
# DISABLE         y = self.bss_result.k_mask_bin_smooth)
# DISABLE       r_poly = self.core.try_scale(
# DISABLE         k_isotropic = k_isotropic,
# DISABLE         k_mask      = k_mask) # XXX may use low-res selection
# DISABLE       if(r_poly < r_start):
# DISABLE         self.bss_result.k_mask            = k_mask
# DISABLE         self.bss_result.k_isotropic       = k_isotropic
# DISABLE         self.bss_result.k_mask_fit_params = k_mask_fit_parameters
# DISABLE         r_start                           = r_poly
# DISABLE         if(self.verbose):
# DISABLE           print >> self.log, "    %s: %6.4f (%s)"%("bulk-solvent:", r_start, "poly-approx.")
    else:
      if(self.verbose):
        print >> self.log, "    %s: %6.4f (%s)"%("bulk-solvent", r, "result rejected")
    return self.r_factor()

  def smooth(self, x):
    result = moving_average2(x = x)
    result_ = flex.double(len(result), 0)
    for i, r in enumerate(result):
      if(r==0): break
      result_[i]=r
    return result_

  def anisotropic_scaling(self, r_start):
    r_expanal, r_poly, r_expmin = None,None,None
    k_anisotropic_expanal, k_anisotropic_poly, \
      k_anisotropic_expmin = None, None, None
    scale_matrix_expanal, scale_matrix_poly, scale_matrix_expmin= None,None,None
    f_model_no_scale_data = self.core.f_model_no_aniso_scale.data()
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
      r_expanal = self.try_scale(
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
      r_poly = self.try_scale(
        k_anisotropic = k_anisotropic_poly)
      if(self.verbose):
        print >> self.log, "      r_poly   : %6.4f"%r_poly
    # pre-analyze
    force_to_use_expmin=False
    if(r_poly<r_expanal and (k_anisotropic_poly<=0).count(True)>0):
      force_to_use_expmin = True
      self.try_expmin = True
    #
    # try expmin
    if(self.try_expmin):
      zero = self.core.f_calc.customized_copy(data =
        flex.complex_double(self.core.f_obs.data().size(), 0))
      fm = mmtbx.f_model.core(
        f_calc  = self.core.f_model_no_aniso_scale,
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
      r_expmin = self.try_scale(
        k_anisotropic = k_anisotropic_expmin)
      if(self.verbose):
        print >> self.log, "    r_expmin   : %6.4f"%r_expmin
      if(force_to_use_expmin):
        self.core = self.core.update(k_anisotropic = k_anisotropic_expmin)
        if(self.verbose):
          print >> self.log, "      b_cart(11,22,33,12,13,23):",\
            ",".join([str("%8.4f"%i).strip() for i in scale_matrix_expmin])
        return
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
      self.core = self.core.update(k_anisotropic = overall_aniso_scale_best)
      r_aniso = self.r_factor()
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
