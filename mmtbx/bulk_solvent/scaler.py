from scitbx.array_family import flex
import sys
from mmtbx import bulk_solvent
from cctbx import adptbx
import boost.python
ext = boost.python.import_ext("mmtbx_f_model_ext")
from cctbx import sgtbx
from mmtbx.bulk_solvent import bulk_solvent_and_scaling
import mmtbx.f_model
import math
from libtbx import group_args
import scitbx.math
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
  except Exception: return x
  return result

def moving_average2(x):
  x_ = x_ = [x[0]] + list(x) + [x[len(x)-1]]
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
  selections = []
  cntr = 0
  # this approximately corresponds to regular splitting of s between 0 and 1.
  l1 = [0,0.5,0.6,0.7,0.8,0.9]
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
  min_ref_low = min(min_ref_low, d_spacings.size())
  fh = min(500,d_spacings.size())
  for t in [1,2,3]:
    for i in xrange(len(selections)):
      count_true = selections[i].count(True)
      if(count_true>0 and count_true<min_ref_low and i+1<len(selections)):
        s = selections[i] | selections[i+1]
        selections[i+1]=s
        selections[i]=None
      elif(count_true<fh and i==len(selections)-1 and selections[i-1] is not None):
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
    selections = sn
  return selections

def binning(unit_cell, miller_indices):
  d_spacings = unit_cell.d(miller_indices)
  r = set_bin_selections(d_spacings = d_spacings, min_ref_low = 300)
  d_sp_0 = d_spacings.select(r[0])
  if(flex.min(d_sp_0)<6 and d_sp_0.size() > 100):
    r = set_bin_selections(d_spacings = d_spacings, min_ref_low = 50)
  elif(flex.min(d_sp_0)<8 and d_sp_0.size() > 100):
    r = set_bin_selections(d_spacings = d_spacings, min_ref_low = 100)
  return r

#def gs(fo, fc, fm):
#  def target(fo, fc, fm, x, scale):
#    fmodel = flex.abs(scale * (fc + x * fm))
#    diff = (fmodel*fmodel-(fo*fo))
#    return flex.sum(diff*diff)
#  k_mask_trial_range = flex.double([i/1000. for i in range(0,1000,10)])
#  t_best = target(fo=fo, fc=fc, fm=fm, x=0, scale=1)
#  x_best = 0
#  for x in k_mask_trial_range:
#    t = target(fo=fo, fc=fc, fm=fm, x=x, scale=1)
#    if(t<t_best):
#      t_best = t
#      x_best = x
#  return x_best

class run(object):
  def __init__(self,
               f_obs,
               f_calc, # can be a sum: f_calc=f_hydrogens+f_calc+f_part
               f_mask, # only one shell is supported
               r_free_flags,
               ss,
               bin_selections=None,
               scale_method="combo",
               number_of_cycles=20, # termination occures much earlier
               auto_convergence_tolerance = 1.e-4,
               log=None,
               auto=True,
               auto_convergence=True,
               bulk_solvent = True,
               try_poly = True,
               try_expanal = True,
               try_expmin = False,
               verbose=False):
    self.log = log
    self.scale_method = scale_method
    self.verbose = verbose
    self.r_free_flags = r_free_flags
    self.ss = ss
    self.bulk_solvent = bulk_solvent
    self.try_poly    = try_poly
    self.try_expanal = try_expanal
    self.try_expmin  = try_expmin
    self.d_spacings = f_obs.d_spacings()
    self.r_low = None
    self.poly_approx_cutoff = None
    self.f_obs = f_obs
    self.r_free_flags = r_free_flags
    self.auto_convergence = auto_convergence
    self.auto_convergence_tolerance = auto_convergence_tolerance
    self.scale_matrices = None
    self.auto = auto
    self.bin_selections = bin_selections
    self.selection_work = miller.array(
      miller_set=self.f_obs,
      data      =~self.r_free_flags.data())
    self.ss = ss
    def init_result():
      return group_args(
        k_mask_bin_orig   = None,
        k_mask_bin_smooth = None,
        k_mask            = None,
        k_isotropic       = None,
        k_mask_fit_params = None)
    self.bss_result = init_result()
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
    self.core = mmtbx.arrays.init(f_calc = f_calc, f_masks = f_mask)
    if(abs(self.core.f_mask()).data().all_eq(0)): self.bulk_solvent=False
    self.cores_and_selections = []
    if(self.bin_selections is None):
      self.bin_selections = mmtbx.bulk_solvent.scaler.binning(
        unit_cell      = self.f_obs.unit_cell(),
        miller_indices = self.f_obs.indices())
    self.low_resolution_selection = self._low_resolution_selection()
    self.high_resolution_selection = self._high_resolution_selection()
    if(verbose):
      print >> log, "  Using %d resolution bins"%len(self.bin_selections)
    self.ss_bin_values=[]
    sel_positive = self.f_obs.data()>0
    self.selection_work = self.selection_work.customized_copy(
      data = self.selection_work.data() & sel_positive)
    for i_sel, sel in enumerate(self.bin_selections):
      core_selected = self.core.select(selection=sel)
      sel_use = self.selection_work.data().select(sel)
      #
      #f_obs_data = self.core.f_obs.data()
      #fodim = flex.mean(core_selected.f_obs.data())
      #sel_use_  = (f_obs_data < fodim*2)
      #sel_use_ &= (f_obs_data > fodim/2)
      #
      #sel_use_ = self.f_obs.select(selection=sel).data()>0
      #sel_use &= sel_use_
      sel_work = sel & self.selection_work.data() #& sel_positive <<<
      #print sel_use.count(True), sel_work.count(True)
      #sel_work &= sel_use_
      #print sel_work.count(True)
      #print
      self.cores_and_selections.append([sel, core_selected, sel_use, sel_work])
      ss = self.ss.select(sel)
      self.ss_bin_values.append([
        flex.min(ss),
        flex.max(ss),
        flex.mean(ss)])
    for cycle in xrange(number_of_cycles):
      r_start = self.r_factor()
      r_start0 = r_start
      if(verbose):
        print >> log, "  cycle %d:"%cycle
        print >> log, "    r(start): %6.4f"%(r_start)
      # bulk-solvent and overall isotropic scale
      if(self.bulk_solvent):
        if(cycle==0):
          r_start = self.k_mask_grid_search(r_start=r_start)
          if(verbose):
            print >> self.log, "    r(bulk_solvent_grid_search): %6.4f"%r_start
        else:
          r_start = self.bulk_solvent_scaling(r_start = r_start)
      # anisotropic scale
      if([try_poly, try_expanal, try_expmin].count(True)):
        if(verbose): print >> log, "    anisotropic scaling:"
        self.anisotropic_scaling(r_start = r_start)
      if(self.auto_convergence and self.is_converged(r_start=r_start0,
         tolerance=self.auto_convergence_tolerance)):
        break
    self.apply_overall_scale()
    if(verbose):
      print >> log, "  r(final): %6.4f"%(self.r_factor())
      self.show()
    #
    self.r_low = self._r_low()
    self.r_high = self._r_high()
    if(verbose):
      d = self.d_spacings.data().select(self.low_resolution_selection)
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
    if(k_mask is None):        k_mask        = self.core.k_mask()
    if(k_anisotropic is None): k_anisotropic = self.core.k_anisotropic
    c = mmtbx.arrays.init(
      f_calc       = self.core.f_calc,
      f_masks      = self.core.f_mask(),
      k_isotropic  = k_isotropic,
      k_anisotropic= k_anisotropic,
      k_masks      = k_mask)
    sel = self.selection_work.data()
    if(selection is not None): sel = selection & sel
    return bulk_solvent.r_factor(self.f_obs.data(), c.f_model.data(), sel)

  def r_factor(self):
    return bulk_solvent.r_factor(self.f_obs.data(),
      self.core.f_model.data(), self.selection_work.data(), 1)

  def _r_high(self):
    return bulk_solvent.r_factor(self.f_obs.data(),
      self.core.f_model.data(), self.high_resolution_selection, 1)

  def _r_low(self):
    return bulk_solvent.r_factor(self.f_obs.data(),
      self.core.f_model.data(), self.low_resolution_selection, 1)

  def _high_resolution_selection(self):
    return self.bin_selections[len(self.bin_selections)-1] & \
      self.selection_work.data()

  def _low_resolution_selection(self):
    d = self.d_spacings.data()
    result = None
    if((d>=8).count(True)>=500):
      result = d>=8
    else:
      sel = self.bin_selections[0].deep_copy()
      if(sel.count(True)>=500): result = sel
      else:
        for i in range(1, len(self.bin_selections)): #XXX
          sel |= self.bin_selections[i].deep_copy()
          if(sel.count(True)>=300): break
        result = sel
    return result & self.selection_work.data()

  def populate_bin_to_individual_k_mask_linear_interpolation(self, k_mask_bin):
    assert len(k_mask_bin) == len(self.cores_and_selections)
    def linear_interpolation(x1,x2,y1,y2):
      k=0
      if(x1!=x2): k=(y2-y1)/(x2-x1)
      b = y1-k*x1
      return k,b
    result1 = flex.double(self.f_obs.size(), -1)
    result2 = flex.double(self.f_obs.size(), -1)
    result  = flex.double(self.f_obs.size(), -1)
    for i, cas in enumerate(self.cores_and_selections):
      selection, zzz, zzz, zzz = cas
      x1,x2 = self.ss_bin_values[i][0], self.ss_bin_values[i][1]
      y1 = k_mask_bin[i]
      if(i==len(k_mask_bin)-1):
        y2 = k_mask_bin[i-1]
      else:
        y2 = k_mask_bin[i+1]
      k,b = linear_interpolation(x1=x1,x2=x2,y1=y1,y2=y2)
      bulk_solvent.set_to_linear_interpolated(self.ss,k,b,selection,result1)
      result2.set_selected(selection, y1)
      r1 = self.try_scale(k_mask = result1, selection=selection) # XXX inefficient
      r2 = self.try_scale(k_mask = result2, selection=selection) # XXX inefficient
      if(r1<r2):
        bulk_solvent.set_to_linear_interpolated(self.ss,k,b,selection,result)
      else:
        result.set_selected(selection, y1)
    assert (result < 0).count(True) == 0
    return result

  def k_mask_grid_search(self, r_start):
    k_mask_trial_range = flex.double([i/1000. for i in range(0,650,50)])
    k_mask      = flex.double(self.f_obs.size(), -1)
    k_isotropic = flex.double(self.f_obs.size(), -1)
    k_mask_bin = flex.double()
    for i_cas, cas in enumerate(self.cores_and_selections):
      selection, core, selection_use, sel_work = cas
      scale = self.core.k_anisotropic.select(selection)
      f_obs = self.f_obs.select(selection)
      k_mask_bin_, k_isotropic_bin_ = \
          bulk_solvent.k_mask_and_k_overall_grid_search(
            f_obs.data()/scale,
            core.f_calc.data(),
            core.f_mask().data(),
            k_mask_trial_range,
            selection_use)
      k_mask_bin.append(k_mask_bin_)
      k_mask.set_selected(selection, k_mask_bin_)
      fm =core.k_anisotropic*(core.f_calc.data()+k_mask_bin_*core.f_mask().data())
      scale_k1 = bulk_solvent.scale(f_obs.data(), fm)
      k_isotropic.set_selected(selection, scale_k1)
    k_mask_bin_smooth = self.smooth(k_mask_bin)
    k_mask = self.populate_bin_to_individual_k_mask_linear_interpolation(
      k_mask_bin = k_mask_bin_smooth)
    k_isotropic = self._k_isotropic_as_scale_k1(k_mask = k_mask)
    self.core = self.core.update(k_masks = k_mask, k_isotropic = k_isotropic)
    self.bss_result.k_mask_bin_orig   = k_mask_bin
    self.bss_result.k_mask_bin_smooth = k_mask_bin_smooth
    self.bss_result.k_mask            = k_mask
    self.bss_result.k_isotropic       = k_isotropic
    ####
    if(len(self.cores_and_selections)>2):
      x=flex.double()
      y=flex.double()
      # XXX use much much finer bins !!!
      for i_sel, cas in enumerate(self.cores_and_selections):
        selection, core, selection_use, sel_work = cas
        sel = sel_work
        ss_ = self.ss_bin_values[i_sel][2]
        k_isotropic_   = flex.mean(self.core.k_isotropic.select(sel))
        x.append(ss_)
        y.append(k_isotropic_)
      import scitbx.math
      r = scitbx.math.gaussian_fit_1d_analytical(x = flex.sqrt(x), y = y)
      r_start = self.r_factor()
      k_isotropic = r.a*flex.exp(-self.ss*r.b)
      r = self.try_scale(k_isotropic = k_isotropic)
      if(r<r_start):
        self.core = self.core.update(k_isotropic = k_isotropic)
    return self.r_factor()

  def apply_overall_scale(self):
    scale_k1 = bulk_solvent.scale(self.f_obs.data(),
      self.core.f_model.data(), self.selection_work.data())
    self.core = self.core.update(k_isotropic = self.core.k_isotropic*scale_k1)

  def is_converged(self, r_start, tolerance=1.e-4):
    self.r_final = self.r_factor()
    result = False
    if((r_start<=self.r_final) or
       (r_start>self.r_final and abs(r_start-self.r_final)<tolerance)):
      result = True
    return result

  def show(self):
    b = self.bss_result
    print >> self.log, "  Statistics in resolution bins:"
    #assert k_mask.size() == len(self.bin_selections)
    fmt="  %7.5f %6.2f -%6.2f %5.1f %5d %-6s %-6s %-6s  %6.3f %6.3f %8.2f %6.4f"
    f_model = self.core.f_model.data()
    print >> self.log, "  s^2      Resolution    Compl Nrefl k_mask                 k_iso  k_ani <Fobs>   R"
    print >> self.log, "                (A)        (%)       orig   smooth average"
    k_mask_bin_orig_   = str(None)
    k_mask_bin_smooth_ = str(None)
    k_mask_bin_approx_ = str(None)
    for i_sel, cas in enumerate(self.cores_and_selections):
      selection, core, selection_use, sel_work = cas
      sel = sel_work
      ss_ = self.ss_bin_values[i_sel][2]
      if(b is not None and self.bss_result.k_mask_bin_orig is not None):
        k_mask_bin_orig_ = "%6.4f"%self.bss_result.k_mask_bin_orig[i_sel]
      if(b is not None and self.bss_result.k_mask_bin_smooth is not None):
        k_mask_bin_smooth_ = "%6.4f"%self.bss_result.k_mask_bin_smooth[i_sel]
      k_mask_bin_averaged_ = "%6.4f"%flex.mean(self.core.k_mask().select(sel))
      d_             = self.d_spacings.data().select(sel)
      d_min_         = flex.min(d_)
      d_max_         = flex.max(d_)
      n_ref_         = d_.size()
      f_obs_         = self.f_obs.select(sel)
      f_obs_mean_    = flex.mean(f_obs_.data())
      k_isotropic_   = flex.mean(self.core.k_isotropic.select(sel))
      k_anisotropic_ = flex.mean(self.core.k_anisotropic.select(sel))
      cmpl_          = f_obs_.completeness(d_max=d_max_)*100.
      r_             = bulk_solvent.r_factor(f_obs_.data(), f_model.select(sel), 1)
      print >> self.log, fmt%(ss_, d_max_, d_min_, cmpl_, n_ref_,
        k_mask_bin_orig_, k_mask_bin_smooth_,k_mask_bin_averaged_,
        k_isotropic_, k_anisotropic_, f_obs_mean_, r_)

  def _k_isotropic_as_scale_k1(self, k_mask):
    result = flex.double(self.ss.size(), -1)
    core_data = mmtbx.arrays.init(
      f_calc        = self.core.f_calc,
      f_masks       = self.core.f_mask(),
      k_anisotropic = self.core.k_anisotropic,
      k_masks       = k_mask).data
    for i_cas, cas in enumerate(self.cores_and_selections):
      selection, core, selection_use, sel_work = cas
      scale_k1 = bulk_solvent.scale(self.f_obs.data(),
        core_data.f_model, sel_work)
      result = result.set_selected(selection, scale_k1)
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
    cutoff = min(cutoff, self.f_obs.d_min()+width)
    sel_high = self.d_spacings.data()<cutoff
    sel_high = sel_high & self.selection_work.data()
    scale_k1 = 1
    if(sel_high.count(True)>min_reflections):
      core = self.core.select(selection = sel_high)
      f_obs = self.f_obs.select(sel_high)
      fm = core.k_anisotropic * (core.f_calc.data() +
        core.k_mask() * core.f_mask().data())
      scale_k1 = bulk_solvent.scale(f_obs.data(), fm)
    return scale_k1

  def bulk_solvent_scaling(self, r_start):
    k_mask     = flex.double(self.f_obs.size(), -1)
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
      selection, core, selection_use, sel_work = cas
      scale = self.core.k_anisotropic.select(selection)*scale_k1
      f_obs  = self.f_obs.select(selection).data()
      f_calc = core.f_calc.data()*scale
      f_mask = core.f_mask().data()*scale
      if(self.scale_method == "k_iso_k_mask_anal"):
        obj = bulk_solvent.overall_and_bulk_solvent_scale_coefficients_analytical(
          f_obs     = f_obs,
          f_calc    = f_calc,
          f_mask    = f_mask,
          selection = selection_use)
        k_mask_bin.append(obj.x_best)
        k_mask.set_selected(selection, obj.x_best)
      elif(self.scale_method == "k_mask_anal"):
        obj = bulk_solvent.bulk_solvent_scale_coefficients_analytical(
          f_obs     = f_obs,
          f_calc    = f_calc,
          f_mask    = f_mask,
          selection = selection_use)
        k_mask_bin.append(obj.x_best)
        k_mask.set_selected(selection, obj.x_best)
      elif(self.scale_method == "k_mask_r_grid_search"):
        k_mask_bin_, k_isotropic_bin_ = \
          bulk_solvent.k_mask_and_k_overall_grid_search(
            f_obs,
            f_calc,
            f_mask,
            k_mask_trial_range,
            selection_use)
        k_mask_bin.append(k_mask_bin_)
        k_mask.set_selected(selection, k_mask_bin_)
      elif(self.scale_method == "combo"):
        r = flex.double()
        k = flex.double()
        x0 = self.bss_result.k_mask_bin_orig[i_cas]
        k_mask.set_selected(selection, x0)
        r0 = self.try_scale(k_mask = k_mask, selection=selection)
        r.append(r0)
        k.append(x0)
        #
        obj1 = bulk_solvent.overall_and_bulk_solvent_scale_coefficients_analytical(
          f_obs     = f_obs,
          f_calc    = f_calc,
          f_mask    = f_mask,
          selection = selection_use)
        k_mask.set_selected(selection, obj1.x_best)
        r.append(self.try_scale(k_mask = k_mask, selection=selection))
        k.append(obj1.x_best)
        #
        obj2 = bulk_solvent.bulk_solvent_scale_coefficients_analytical(
          f_obs     = f_obs,
          f_calc    = f_calc,
          f_mask    = f_mask,
          selection = selection_use)
        k_mask.set_selected(selection, obj2.x_best)
        r.append(self.try_scale(k_mask = k_mask, selection=selection))
        k.append(obj2.x_best)
        s = flex.sort_permutation(r)
        x = k.select(s)[0]
        # fine-sample k_mask around minimum of LS to fall into minimum of R
        k_mask_bin_, k_isotropic_bin_ = \
          bulk_solvent.k_mask_and_k_overall_grid_search(
            f_obs,
            f_calc,
            f_mask,
            get_k_mask_trial_range(x = x),
            selection_use)
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
      self.core = self.core.update(k_isotropic = k_isotropic, k_masks = k_mask)
      self.bss_result.k_mask_bin_orig   = k_mask_bin
      self.bss_result.k_mask_bin_smooth = k_mask_bin_smooth
      self.bss_result.k_mask            = k_mask
      self.bss_result.k_isotropic       = k_isotropic
      r_start = r
      if(self.verbose):
        print >> self.log, "    %s: %6.4f"%("bulk-solvent", r_start)
    else:
      if(self.verbose):
        print >> self.log, "    %s: %6.4f (%s)"%("bulk-solvent", r, "result rejected")
    return self.r_factor()

  def smooth(self, x):
    result = moving_average2(x = x)
    result_ = flex.double(len(result), 0)
    for i, r in enumerate(result):
      d = 1/math.sqrt(self.ss_bin_values[i][1])/2
      if(r==0 and d<3): break
      result_[i]=r
    return result_

  def anisotropic_scaling(self, r_start):
    r_expanal, r_poly, r_expmin = None,None,None
    k_anisotropic_expanal, k_anisotropic_poly, \
      k_anisotropic_expmin = None, None, None
    scale_matrix_expanal, scale_matrix_poly, scale_matrix_expmin= None,None,None
    sel     = self.selection_work.data()
    f_model = self.core.f_model_no_aniso_scale.data().select(sel)
    f_obs   = self.f_obs.data().select(sel)
    mi      = self.f_obs.indices().select(sel)
    uc      = self.f_obs.unit_cell()
    mi_all  = self.f_obs.indices()
    # try exp_anal
    if(self.try_expanal):
      obj = bulk_solvent.aniso_u_scaler(
        f_model        = f_model,
        f_obs          = f_obs,
        miller_indices = mi,
        adp_constraint_matrix = self.adp_constraints.gradient_sum_matrix())
      u_star = self.adp_constraints.all_params(tuple(obj.u_star_independent))
      scale_matrix_expanal = adptbx.u_as_b(adptbx.u_star_as_u_cart(uc, u_star))
      k_anisotropic_expanal = ext.k_anisotropic(mi_all, u_star)
      r_expanal = self.try_scale(k_anisotropic = k_anisotropic_expanal)
      if(self.verbose):
        print >> self.log, "      r_expanal: %6.4f"%r_expanal
    # try poly
    if(self.try_poly):
      obj = bulk_solvent.aniso_u_scaler(
        f_model        = f_model,
        f_obs          = f_obs,
        miller_indices = mi,
        unit_cell      = uc)
      scale_matrix_poly = obj.a
      k_anisotropic_poly = ext.k_anisotropic(mi_all, obj.a, uc)
      r_poly = self.try_scale(k_anisotropic = k_anisotropic_poly)
      if(self.verbose):
        print >> self.log, "      r_poly   : %6.4f"%r_poly
    # pre-analyze
    force_to_use_expmin=False
    if(self.auto and r_poly<r_expanal and (k_anisotropic_poly<=0).count(True)>0):
      force_to_use_expmin = True
      self.try_expmin = True
    # try expmin
    if(self.try_expmin):
      zero = self.f_obs.select(sel).customized_copy(data =
        flex.complex_double(f_obs.size(), 0))
      fm = mmtbx.f_model.manager_kbu(
        f_obs         = self.f_obs.select(sel),
        f_calc        = self.core.f_model_no_aniso_scale.select(sel),
        f_masks       = [zero],
        f_part1       = zero,
        f_part2       = zero,
        ss            = self.ss)
      obj = bulk_solvent_and_scaling.u_star_minimizer(
        fmodel_core_data = fm,
        f_obs            = self.f_obs.select(sel),
        u_initial        = [0,0,0,0,0,0],
        refine_u         = True,
        min_iterations   = 500,
        max_iterations   = 500,
        symmetry_constraints_on_b_cart = True,
        u_min_max = 500.,
        u_min_min =-500.)
      u_star = obj.u_min
      scale_matrix_expmin = adptbx.u_as_b(adptbx.u_star_as_u_cart(uc, u_star))
      k_anisotropic_expmin = ext.k_anisotropic(mi_all, u_star)
      r_expmin = self.try_scale(k_anisotropic = k_anisotropic_expmin)
      if(self.verbose): print >> self.log, "    r_expmin   : %6.4f"%r_expmin
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
    k_anisotropic_best = None
    scale_matrix_best = None
    for result in r:
      r_factor, k_anisotropic, scale_matrix = result
      if(r_factor is not None and r_factor < r_best):
        r_best = r_factor
        k_anisotropic_best = k_anisotropic.deep_copy()
        scale_matrix_best = scale_matrix[:]
    if(scale_matrix_best is None):
      if(self.verbose):
        print >> self.log, "      result rejected due to r-factor increase"
    else:
      self.scale_matrices = scale_matrix_best
      self.core = self.core.update(k_anisotropic = k_anisotropic_best)
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
