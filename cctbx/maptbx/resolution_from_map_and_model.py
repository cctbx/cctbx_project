from __future__ import division
from cctbx.array_family import flex
from libtbx import adopt_init_args
from libtbx.test_utils import approx_equal
from cctbx import maptbx
from cctbx import miller
import scitbx.math.curve_fitting
import iotbx.pdb
from libtbx import group_args

def canal(x, y, simple=False, assert_concave_up=False,
          use_longest_slope_criteria=False, show=False, smooth=False,
          xround=False, of=None):
#  xmin, xmax = flex.min(x), flex.max(x)
#  fit = scitbx.math.curve_fitting.univariate_polynomial_fit(x_obs=x, y_obs=y,
#    degree=2, min_iterations=50, number_of_cycles=10)
#  a0,a1,a2 = fit.params
#  xpmax=-a1/(2*a2)
#  inside=False
#  if(xpmax>xmin and xpmax<xmax): inside=True
#  # a0+a1*x+a2*x**2, 2*a2*x + a1=0, x = -a1/(2*a2)
#  #if(a2>=0. or abs(a2)<1.e-6): return [] # concave down or flat line
#  #if(a2>=0.): return [] # concave down or flat line
#  #
  # Remove duplicates
  d = {}
  for i,j in zip(x,y):
    d.setdefault(i,[]).append(j)
  x_,y_ = flex.double(), flex.double()
  for i,j in zip(d.keys(), d.values()):
    x_.append(i)
    y_.append(max(j))
  x,y = x_[:], y_[:]
  # Special cases of few values:
  if(x.size()==0): return []
  if(x.size()==1): return []
  if(x.size()==2): return []
  # Sort
  s = flex.sort_permutation(x)
  x = x.select(s)
  y = y.select(s)
  # Smooth and round
  #if(not simple):
  if(smooth): y = moving_average(y)
  if(xround): y = flex.double([round(i,4) for i in y])
  # Show
  if(show):
    if(of is not None):
      for i,j in zip(x,y):
        print >> of,  "%5.2f %9.6f"%(i,j)
  # Find maxima
  maxima = []
  if(simple):
    s = flex.sort_permutation(y, reverse=True)
    maxima.append([x.select(s)[0], y.select(s)[0], None])
    return maxima
  else:
    for i in xrange(y.size()):
      yi = y[i]
      cntr=0
      yp,ym=[],[]
      offsets = [1,]#[1,2,3]
      for j in offsets:
        if(i-j>=0 and y[i-j]<yi and not y[i-j] in ym):
          cntr+=1
          ym.append(y[i-j])
        if(i+j<=y.size()-1 and y[i+j]<yi and not y[i+j] in yp):
          cntr+=1
          yp.append(y[i+j])
      if(cntr>=2 and len(ym)==len(yp)): maxima.append([x[i],yi,i])
    #
    for i in xrange(y.size()):
      yi = y[i]
      if(i-1>=0 and y[i-1]<yi):
        offsets = [1,2,3,4,5,6,7,8,9,10]
        xi = [x[i]]
        yij=None
        for j in offsets:
          if(i+j<=y.size()-1):
            yij=y[i+j]
            if(abs(yij-yi)<1.e-6):
              xi.append(x[i+j])
            else:
              break
        if(yij is not None and not yij<yi): xi=[]
        if(len(xi)<6 and len(xi)>1):
          maxima.append([flex.mean(flex.double(xi)),yi])
    #
    if(len(maxima)>1):
      tmp = []
      for i,mi in enumerate(maxima):
        for j,mj in enumerate(maxima):
          if(i==j): continue
          if(mi[1]>mj[1] and abs(mi[1]-mj[1])>0.01):
            tmp.append(mi)
      maxima = tmp[:]
    # Final final
    if(len(maxima)==0 and x.size()>3):
      ymax = flex.max(y)
      if(y[0]==ymax and y[1]<ymax and y[2]<ymax and y[3]<ymax and (ymax-y[3])>0.01):
        maxima = [[x[0],y[0], 0]]
    return maxima

def moving_average(y):
    y_=[]
    for i in xrange(y.size()):
      if(i>=1 and i<=y.size()-2):
        y_ave = (y[i-1]+y[i]+y[i+1])/3.
      elif(i==0):
        y_ave = (y[0]+y[1])/2.
      elif(y==y.size()-1):
        y_ave = (y[y.size()-2]+y[y.size()-1])/2.
      else:
        assert 0
      #if(i>=2 and i<=y.size()-3):
      #  y_ave = (y[i-2]+y[i-1]+y[i]+y[i+1]+y[i+2])/5.
      #elif(i>=1 and i<=y.size()-2):
      #  y_ave = (y[i-1]+y[i]+y[i+1])/3.
      #else:
      #  y_ave = y[i]
      y_.append(y_ave)
    return flex.double(y_)


class run_loop(object):
  def __init__(self, fc, f_obs, b_iso, map, d_mins, d_spacings, ss):
    adopt_init_args(self, locals())
    result = maptbx.cc_complex_complex(
      f_1        = self.f_obs.data(),
      f_2        = self.fc.data(),
      d_spacings = self.d_spacings,
      ss         = self.ss,
      d_mins     = self.d_mins,
      b_iso      = self.b_iso)
    self.d_min, self.cc = result[0], result[1]

def get_fo_fc_adjust_d_min_start(map, fc, xrs, d_mins):
  d_min_start = flex.min(d_mins)
  while True:
    f_obs = None
    try:
      f_obs = fc.structure_factors_from_map(
        map            = map,
        use_scale      = True,
        anomalous_flag = False,
        use_sg         = False)
    except RuntimeError, e:
      pass
    if(f_obs is not None): break
    d_min_start += 0.1
    fc = fc.resolution_filter(d_min=d_min_start)
  d_min_start = fc.d_min()
  s = d_mins>=d_min_start
  d_mins = d_mins.select(s)
  return f_obs, fc, d_mins

def _resolution_from_map_and_model_helper_2(
      map,
      f_obs, fc,
      xray_structure,
      b_range,
      nproc,
      simple,
      d_mins,
      ofn=None):
  ccs    = flex.double()
  bs     = flex.double()
  d_min_opt = flex.double()
  cntr=0
  of=None
  if(ofn is not None): of = open(ofn,"w")
  d_spacings = fc.d_spacings().data()
  ss = 1./flex.pow2(d_spacings) / 4.
  for b in b_range:
    o = run_loop(fc=fc, f_obs=f_obs, b_iso=b, map=map,
      d_mins=d_mins, d_spacings=d_spacings, ss=ss)
    ccs.append(o.cc)
    d_min_opt.append(o.d_min)
    bs.append(b)
  maxima = canal(x=d_min_opt, y=ccs, simple=simple, smooth=True, xround=True,
                 show=True, of=of)
  if(ofn is not None): of.close()
  d_min_result = None
  cc_result = None
  b_result = None
  if(len(maxima)>0):
    d_min_result = maxima[0][0]
    cc_result = maxima[0][1]
    dist=1.e+9
    for i, cc in enumerate(ccs):
      dist_= abs(cc-cc_result)
      if(dist_<dist):
        dist=dist_
        b_result=bs[i]
  return d_min_result, b_result, cc_result

def get_trial_resolutions(map_data, unit_cell, d_min_min):
  d_min_end = round(maptbx.d_min_from_map(
    map_data=map_data, unit_cell=unit_cell, resolution_factor=0.1),1)
  d_min_start = round(maptbx.d_min_from_map(
    map_data=map_data, unit_cell=unit_cell, resolution_factor=0.5),1)
  if(d_min_min is not None and d_min_start<d_min_min):
    d_min_start=d_min_min
  d_min_end = max(8,d_min_end)
  l1 = [i/10. for i in range(int(d_min_min*10),30)]
  #l2 = [i/100. for i in range(300,625,25)]
  l2 = [i/100. for i in range(300,610,10)]
  l3 = [i/10. for i in range(65,105,1)]
  l4 = [i*1. for i in range(11,100,1)]
  l = flex.double(l1+l2+l3+l4)
  s  = l>=d_min_start
  s &= l<=d_min_end
  return l.select(s)

def run_at_b(b, map_data, xray_structure, d_min_min=None):
  m1 = flex.mean(map_data)
  xrs = xray_structure.deep_copy_scatterers().set_b_iso(value=0)
  d_mins = get_trial_resolutions(
    map_data  = map_data,
    unit_cell = xray_structure.unit_cell(),
    d_min_min = d_min_min)
  d_min_start = flex.min(d_mins)
  fc = xrs.structure_factors(d_min=d_min_start).f_calc()
  f_obs, fc, d_mins = get_fo_fc_adjust_d_min_start(
    map    = map_data,
    xrs    = xrs,
    fc     = fc,
    d_mins = d_mins)
  d_spacings = fc.d_spacings().data()
  ss = 1./flex.pow2(d_spacings) / 4.
  m2 = flex.mean(map_data)
  assert approx_equal(m1, m2)
  return run_loop(fc=fc, f_obs=f_obs, b_iso=b, map=map, d_mins=d_mins,
    d_spacings=d_spacings, ss=ss)

def run_fast(map_data, f_obs, d_fsc_model, xray_structure):
  #fo = f_obs.resolution_filter(d_min = d_fsc_model)
  xrs = xray_structure.deep_copy_scatterers()
  xrs = xrs.set_b_iso(value=0)
  fc = f_obs.structure_factors_from_scatterers(xray_structure = xrs).f_calc()
  fo = f_obs.resolution_filter(d_min = d_fsc_model)
  fo, fc, = fo.common_sets(fc)
  del xrs
  cc = -999
  b=None
  ss = 1./flex.pow2(fc.d_spacings().data()) / 4.
  data = fc.data()
  for b_ in range(-500,500,1):
    sc = flex.exp(-b_*ss)
    fc_ = fc.customized_copy(data = data*sc)
    cc_ = fo.map_correlation(other = fc_)
    if(cc_>cc):
      cc = cc_
      b = b_
  o = run_at_b(
    b                = b,
    map_data         = map_data,
    xray_structure   = xray_structure,
    d_min_min        = 1.7)
  return group_args(d_min = o.d_min, b_iso = b)

class run(object):
  def __init__(self, map_data, xray_structure, pdb_hierarchy, d_min_min=None,
                     nproc=1, ofn=None):
    """
    Given map and model estimate resolution by maximizing map CC(map, model-map).
    As a by-product, also provides CC and optimal overall B-factor.
    """
    # XXX xray_structure, pdb_hierarchy are modified!
    xray_structure = xray_structure.deep_copy_scatterers()
    pdb_hierarchy  = pdb_hierarchy.deep_copy()
    #
    unit_cell = xray_structure.unit_cell()
    d_mins = get_trial_resolutions(map_data=map_data,
      unit_cell=xray_structure.unit_cell(), d_min_min=d_min_min)
    ###
    xrs = xray_structure.deep_copy_scatterers().set_b_iso(value=0)
    assert approx_equal(flex.mean(xrs.extract_u_iso_or_u_equiv()),0.)
    d_min_start = flex.min(d_mins)
    fc = xrs.structure_factors(d_min=d_min_start).f_calc()
    f_obs, fc, d_mins = get_fo_fc_adjust_d_min_start(map=map_data, xrs=xrs,
      fc=fc, d_mins=d_mins)
    b_range=list(range(-500,510,10))
    self.d_min, self.b_iso, self.cc = _resolution_from_map_and_model_helper_2(
      map            = map_data,
      f_obs=f_obs, fc=fc,
      xray_structure = xray_structure,
      b_range        = b_range,
      d_mins         = d_mins,
      nproc          = nproc,
      simple=True,
      ofn=None)
