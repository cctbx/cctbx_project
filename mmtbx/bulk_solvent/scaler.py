from scitbx.array_family import flex
import sys, math
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

#def moving_average(x, offset, d_spacings):
#  result = flex.double(x.size(), -1)
#  for i in xrange(x.size()):
#    d = d_spacings[i]
#    if(d>20): offset=2
#    if(d<=20 and d>10): offset = 3
#    if(d<=10 and d>5):  offset = 5
#    if(d<=5 and d>2):   offset = 10
#    if(d<=2):   offset = 20
#    s = 0
#    cntr = 0
#    for j in range(max(0,i-offset), min(x.size()-1, i+offset)):
#      s+=x[j]
#      cntr+=1
#    if(cntr!=0): result[i]=s/cntr
#    else: result[i]=0
#  return result

def moving_average(x, offset):
  result = flex.double(x.size(), -1)
  for i in xrange(x.size()):
    s = 0
    cntr = 0
    for j in range(max(0,i-offset), min(x.size()-1, i+offset)):
      s+=x[j]
      cntr+=1
    if(cntr!=0): result[i]=s/cntr
    else: result[i]=0
  return result

def set_bin_selections(d_spacings):
  selections = []
  cntr = 0
  #ranges = [(0   , 1),
  #          (1   , 1.25),
  #          (1.25, 1.5),
  #          (1.5 , 1.75),
  #          (1.75, 2),
  #
  #          (2  , 2.5),
  #          (2.5, 3),
  #          (3  , 3.5),
  #          (3.5, 4),
  #
  #          (4  , 5),
  #          (5  , 6),
  #          (6  , 7),
  #          (7  , 8),
  #          (8  , 9),
  #          (9  , 10),
  #
  #          (10 , 11),
  #          (11 , 12),
  #          (12 , 13),
  #          (13 , 14),
  #          (14 , 15),
  #
  #          (15 , 20),
  #          (20 , 25),
  #          (25 , 30),
  #          (30 , 35),
  #          (35 , 40),
  #          (40 , 45),
  #          (45 , 50),
  #          (50 , 999)
  #          ]
  ranges = [(0   , 1),
            (1   , 1.25),
            (1.25, 1.5),
            (1.5 , 1.75),
            (1.75, 2),

            #(2  , 2.5),
            #(2.5, 3),
            #(3  , 3.5),
            #(3.5, 4),
            (2  , 3),
            (3  , 4),

            (4  , 5),
            (5  , 6),
            #(6  , 7),
            #(7  , 8),
            #(8  , 9),

            (6  , 8),
            (8  , 10),

            #(10,20),
            #(20,999)

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
            (50 , 999)
            ]
  ranges.reverse()
  for s in ranges:
    sel  = d_spacings >= s[0]
    sel &= d_spacings <  s[1]
    #if(sel.count(True)>0):
    #  print s, sel.count(True), flex.min(d_spacings.select(sel)),flex.max(d_spacings.select(sel))
    cntr += sel.count(True)
    if(sel.count(True)>0): selections.append(sel)
  assert cntr == d_spacings.size()
  for i in xrange(len(selections)):
    if(selections[i].count(True)<100 and i+1<len(selections)):
      s = selections[i] | selections[i+1]
      selections[i+1]=s
      selections[i]=None
    elif(selections[i].count(True)<500 and i==len(selections)-1 and selections[i-1] is not None):
      #print selections[i] , selections[i-1]
      s = selections[i] | selections[i-1]
      selections[i-1]=s
      selections[i]=None
  #print
  sn = []
  for s in selections:
    if(s is not None): sn.append(s)
  cntr = 0
  for s in sn:
    #print s.count(True)
    cntr += s.count(True)
  assert cntr == d_spacings.size()
  return sn


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
    d_spacings = f_obs.d_spacings().data()
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
    self.bin_selections = set_bin_selections(d_spacings = self.core.f_obs.d_spacings().data())
    #print "Number of zones:", len(sels)
    for i_sel, sel in enumerate(self.bin_selections):
      core_selected = self.core.select(selection=sel)
      f_obs_data = core_selected.f_obs.data()
      fodim = flex.mean(f_obs_data)
      sel_  = (f_obs_data < fodim*3)
      sel_ &= (f_obs_data > fodim/3)
      bin_selections = [[flex.bool(f_obs_data.size(), True).iselection(), sel_]]
      self.cores_and_selections.append([sel,core_selected,bin_selections])
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
      if(cycle==0): #
        self.bulk_solvent_simple(r_start=r_start)
        if(verbose):
          print >> self.log, "r(bulk_solvent_simple): %6.4f"%self.core.r_factor()
      else:

        #self.bulk_solvent_dev(r_start = r_start)

        ssi,x,k = self.bulk_solvent_scaling(r_start = r_start)
        self.fit_poly_k_mask(ssi=ssi, x=x, k=k, r_start = self.core.r_factor(),cycle=cycle)

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
      assert x.size()==k.size() and x.size()==ssi.size()
      #
      #self.fit_poly_k_mask(ssi=ssi, x=x, k=k)
      #


      #for a,b,c,d,e,f in zip(ssi,x,
      #                       moving_average(x=x,offset=1),
      #                       moving_average(x=x,offset=2),
      #                       moving_average(x=x,offset=3),
      #                       moving_average(x=x,offset=4)):
      #  print >> log, "%8.6f %8.6f %8.6f %8.6f %8.6f %8.6f"%(a,b, c,d,e,f)
      #add = False
      #for i in [0]:#[0,1,2,3,4]:
      #  if i == 0:
      #    sel = x>0
      #    xi = x.select(sel)
      #    ssii = ssi.select(sel)
      #  else:
      #    xi = moving_average(x=x,offset=i)
      #    sel = xi>0
      #    xi = xi.select(sel)
      #    ssii = ssi.select(sel)
      #  r = scitbx.math.gaussian_fit_1d_analytical(x = ssii, y = xi)
      #  print >> log, "k_sol, b_sol: %5.3f %7.2f"%(r.a, r.b)
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

  def fit_poly_k_mask(self, ssi, x, k, r_start, cycle):
    X = ssi*ssi
    Y = x
    #Y = moving_average(x=x,offset=2)
    #
    for i in xrange(len(Y)):
      if(i!=0 and i!=len(Y)-1):
        if(Y[i]<=0 and Y[i-1]>0 and Y[i+1]>0):
          Y[i]=(Y[i-1]>0 + Y[i+1])/2.
    #
    sX = flex.sort_permutation(X)
    X = X.select(sX)
    Y = Y.select(sX)
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
    def foo(a, x):
      return a[0] + a[1] * x**1 + a[2] * x**2 + a[3] * x**3
    kmask = flex.double(self.core.ss.size(), 0)
    for i, ss_ in enumerate(self.core.ss):
      v = foo(a=cfo.params,x=ss_)
      if(ss_>=cutoff_x): break
      else: kmask[i]=v
    #for i, ss_ in enumerate(self.core.ss):
    #  print "%10.6f %10.6f"%(ss_,kmask[i])
    #
    if 1:
      overall_scale = flex.double(self.core.ss.size(), -1)
      core_data = ext.core(
        f_calc                    = self.core.f_calc.data(),
        f_mask                    = self.core.f_mask.data(),
        scale                     = 1,
        overall_scale             = flex.double(self.core.ss.size(), 1),
        overall_anisotropic_scale = self.core.overall_anisotropic_scale,
        bulk_solvent_scale        = kmask)
      for sel in self.bin_selections:
        scale_k1 = bulk_solvent.scale(
          self.core.f_obs.data().select(sel), core_data.f_model.select(sel))
        overall_scale = overall_scale.set_selected(sel, scale_k1)
      assert overall_scale.count(-1.) == 0
    else:
      core_data = ext.core(
        f_calc                    = self.core.f_calc.data(),
        f_mask                    = self.core.f_mask.data(),
        scale                     = 1,
        overall_scale             = self.core.overall_scale,
        overall_anisotropic_scale = self.core.overall_anisotropic_scale,
        bulk_solvent_scale        = kmask)
      scale_k1 = bulk_solvent.scale(self.core.f_obs.data(), core_data.f_model)
      overall_scale = self.core.overall_scale*scale_k1

    r = self.core.try_bulk_solvent_scale(
      overall_scale      = overall_scale,
      bulk_solvent_scale = kmask)
    suffix = ""
    if(r>r_start): suffix = "(result rejected due to r-factor increase)"
    if(self.verbose):
      print >> self.log, "    bulk-solvent(polynomial approximation):"
      print >> self.log, "      r        : %6.4f %s"%(r, suffix)
    if(r<r_start):
      self.core.update(
        overall_scale      = overall_scale,
        bulk_solvent_scale = kmask)
    #for ssi_,k_,x_ in zip(ssi,k,x):
    #  v = foo(a=cfo.params,x=ssi_*ssi_)
    #  if(ssi_*ssi_>=cutoff_x): v = 0
    #  print >> self.log, "%8.6f %8.6f %8.6f %8.6f"%(ssi_*ssi_,k_,x_, v)


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
    k = flex.double()
    for i_cas, cas in enumerate(self.cores_and_selections):
      sel = cas[0]
      scale = self.core.overall_anisotropic_scale.select(sel)
      a,b,c,d,e = self.scale(core=cas[1], scale=scale, bin_selections=cas[2])
      if(i_cas == 0):
        ssi,x,k = a,b,e
        overall_scale.set_selected(sel, c)
        bulk_solvent_scale.set_selected(sel, d)
      else:
        ssi.extend(a)
        x.extend(b)
        k.extend(e)
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
    return ssi,x,k

  #####
#  def bulk_solvent_dev(self, r_start):
#    def rr(min1, max1, step):
#      br = flex.double()
#      b = min1
#      while b <= max1:
#        br.append(b)
#        b+=step
#      return br
#    def sc(fo,fc,fm):
#      x_range = rr(min1=-1.5, max1=1.5, step=0.001)
#      r_best = 1.e+9
#      x_best = None
#      for x in x_range:
#        f_model = fc+x*fm
#        r = bulk_solvent.r_factor(fo, f_model, 1)
#        if(r<r_best):
#          r_best = r
#          x_best = x
#      return x_best
#    #
#    scale_k1 = bulk_solvent.scale(self.core.f_obs.data(), self.core.f_model().data())
#    scale_k1a = flex.double(self.core.f_obs.data().size(), scale_k1)
#    oas = self.core.overall_anisotropic_scale
#    overall_scale = scale_k1*oas
#    #
#    f_obs = self.core.f_obs
#    fo = f_obs.data()
#    fc = self.core.f_calc.data()*overall_scale
#    fm = self.core.f_mask.data()*overall_scale
#    f_obs.setup_binner(reflections_per_bin = 250, n_bins = 0)
#    d_spacings = f_obs.d_spacings().data()
#    a = flex.double()
#    b = flex.double()
#    c = flex.double()
#    bulk_solvent_scale = flex.double(f_obs.data().size(), 0)
#    #selections = []
#
#
#    #for i_bin in f_obs.binner().range_used():
#    for sel in xxx(d_spacings = d_spacings):
#      #sel = f_obs.binner().selection(i_bin)
#      fo_ = fo.select(sel)
#      fc_ = fc.select(sel)
#      fm_ = fm.select(sel)
#      d_  = flex.mean(d_spacings.select(sel))
#      ss_  = flex.mean(self.core.ss.select(sel))
#      x = sc(fo=fo_,fc=fc_,fm=fm_)
#      #print "%8.6f %8.4f %8.4f"%(ss_, d_, x)
#      a.append(d_)
#      c.append(ss_)
#      b.append(x)
#      bulk_solvent_scale = bulk_solvent_scale.set_selected(sel, x)
#      #selections.append(sel)
#    #
#    bb = moving_average(x=b, offset=5, d_spacings=d_spacings)
#    #for i, sel in enumerate(selections):
#    #  bulk_solvent_scale = bulk_solvent_scale.set_selected(sel, bb[i])
#    r = self.core.try_bulk_solvent_scale(
#      overall_scale      = scale_k1a,
#      bulk_solvent_scale = bulk_solvent_scale)
#    suffix = ""
#    if(r>r_start): suffix = "(result rejected due to r-factor increase)"
#    if(self.verbose):
#      print >> self.log, "    bulk-solvent:"
#      print >> self.log, "      r        : %6.4f %s"%(r, suffix)
#    if(r<r_start):
#      core_data = ext.core(
#        f_calc                    = self.core.f_calc.data(),
#        f_mask                    = self.core.f_mask.data(),
#        scale                     = 1,
#        overall_scale             = flex.double(oas.size(),1),
#        overall_anisotropic_scale = oas,
#        bulk_solvent_scale        = bulk_solvent_scale)
#      scale_k1a = flex.double(oas.size(), bulk_solvent.scale(f_obs.data(),
#        core_data.f_model))
#      self.core.update(
#        overall_scale      = scale_k1a,
#        bulk_solvent_scale = bulk_solvent_scale)
#
#
#    for a_, c_, b_, b1_ in zip(a, c, b, bb):
#      print "%9.6f %9.7f %9.6f %9.6f"%(a_, c_, b_, b1_)

  def anisotropic_scaling(self, r_start):
    sc=bulk_solvent.scale(self.core.f_obs.data(), self.core.f_model().data())
    scale_k1a = self.core.overall_scale*sc
    self.core.update(overall_scale = scale_k1a)
    #
    r_expanal, r_poly, r_expmin = None,None,None
    overall_anisotropic_scale_expanal, overall_anisotropic_scale_poly, \
      overall_anisotropic_scale_expmin = None, None, None
    scale_matrix_expanal, scale_matrix_poly, scale_matrix_expmin= None,None,None
    # try exp_anal
    if(self.try_expanal):
      obj = bulk_solvent.aniso_u_scaler(
        f_model        = self.core.f_model_no_scale().data(),
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
        f_model        = self.core.f_model_no_scale().data(),
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
        symmetry_constraints_on_b_cart = False,
        u_min_max = 500.,
        u_min_min =-500.)
      u_star = obj.u_min
      scale_matrix_expmin = adptbx.u_as_b(adptbx.u_star_as_u_cart(
        self.core.f_obs.unit_cell(), u_star))
      overall_anisotropic_scale_expmin = ext.overall_anisotropic_scale(
        self.core.f_obs.indices(), u_star)
      r_expmin = self.core.try_overall_anisotropic_scale(
        scale_array = overall_anisotropic_scale_expmin)
      if(self.verbose):
        print >> self.log, "    r_expmin   : %6.4f"%r_expmin
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
    k = flex.double()
    sels = []
    for sel in bin_selections:
      sel, sel_use = sel
      ss_ = flex.mean(flex.sqrt(ss.select(sel)))
      #print "-"*100
      #for i, f in enumerate(f_obs.data().select(sel_use)):
      #  print f, flex.mean(f_obs.data().select(sel_use))
      #print
      obj = bulk_solvent.overall_and_bulk_solvent_scale_coefficients_analytical(
        f_obs     = f_obs.data(),
        f_calc    = f_calc.data(),
        f_mask    = f_mask.data(),
        selection = sel_use)
      ssi.append(ss_)
      x.append(obj.x)
      k.append(obj.t)
      overall_scale=overall_scale.set_selected(sel, obj.t)
      bulk_solvent_scale=bulk_solvent_scale.set_selected(sel, obj.x)
    #assert (overall_scale < 0).count(True) == 0
    #assert (bulk_solvent_scale < 0).count(True) == 0
    return ssi, x, overall_scale, bulk_solvent_scale, k
