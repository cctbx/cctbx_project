from __future__ import division
from cctbx.array_family import flex
from libtbx import adopt_init_args
from libtbx.test_utils import approx_equal
from cctbx import maptbx
import scitbx.math.curve_fitting

def parabola_is_good(x, y, xmin, xmax, xstep):
  fit = scitbx.math.curve_fitting.univariate_polynomial_fit(x_obs=x, y_obs=y,
    degree=2, min_iterations=50, number_of_cycles=10)
  a0,a1,a2 = fit.params
  if(a2>=0.): return [] # concave down
  maxima = []
  for i in xrange(y.size()):
    yi = y[i]
    cntr=0
    yp,ym=[],[]
    offsets = [1,2,3]
    for j in offsets:
        if(i-j>=0 and y[i-j]<=yi and not y[i-j] in ym):
          cntr+=1
          ym.append(y[i-j])
        if(i+j<=y.size()-1 and y[i+j]<=yi and not y[i+j] in yp):
          cntr+=1
          yp.append(y[i+j])
    if(cntr==len(offsets)+len(offsets)): maxima.append([x[i],yi,i])
  #print "       ",len(maxima), maxima
  # Remove plateau
  tmp=[]
  tmp2 = []
  for m in maxima:
    if(not m[1] in tmp):
        tmp.append(m[1])
        tmp2.append(m)
  maxima = tmp2
  # Choose between several peaks
  if(len(maxima)==2):
    x_max = -a1/(2*a2)
    if(x_max<xmin or x_max>xmax): return []
    d=1.e+9
    maximum=None
    for m in maxima:
        d_=abs(x_max-m[0])
        if(d_<d):
          d=d_
          maximum=m[:]
    return [maximum]
  return maxima

def _resolution_from_map_and_model_helper(
      map,
      xray_structure,
      b_range,
      nproc,
      d_min_start=1.5,
      d_min_end=10.,
      radius=5.0,
      d_min_step=0.1,
      approximate=False):
  from cctbx import miller
  import mmtbx.bulk_solvent
  xrs = xray_structure.deep_copy_scatterers().set_b_iso(value=0)
  cg = maptbx.crystal_gridding(
    unit_cell             = xray_structure.unit_cell(),
    space_group_info      = xray_structure.space_group_info(),
    pre_determined_n_real = map.accessor().all())
  ### Embedded class to get radii (parallel enabled)
  class run_get_selection(object):
    def __init__(self, xrs, map, nproc):
      adopt_init_args(self, locals())
      self.uc = self.xrs.unit_cell()
      self.sites_cart = self.xrs.sites_cart()
      self.size = self.xrs.scatterers().size()
      self.selections = []
      self.radii = [2,2.25,2.5,2.75,3,3.5,4,4.5,5]
      if(self.nproc>1):
        from libtbx import easy_mp
        stdout_and_results = easy_mp.pool_map(
          processes    = self.nproc,
          fixed_func   = self.run,
          args         = self.radii,
          func_wrapper = "buffer_stdout_stderr")
        for so, sel in stdout_and_results:
          self.selections.append(sel)
      else:
        for radius in self.radii:
          self.selections.append(self.run(radius=radius))
    def run(self, radius):
      return maptbx.grid_indices_around_sites(
        unit_cell  = self.uc,
        fft_n_real = self.map.focus(),
        fft_m_real = self.map.all(),
        sites_cart = self.sites_cart,
        site_radii = flex.double(self.size, radius))
  ###
  so = run_get_selection(xrs=xrs, map=map, nproc=nproc)
  selections = so.selections
  assert approx_equal(flex.mean(xrs.extract_u_iso_or_u_equiv()),0.)
  fc = xrs.structure_factors(d_min=d_min_start).f_calc()
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
  #
  def get_cc(m1, m2, selections, radii):
    result = -9999.
    radius = None
    for sel, rad in zip(selections, radii):
      m1_ = m1.select(sel)
      m2_ = m2.select(sel)
      cc = flex.linear_correlation(
        x=m1_.as_1d(),
        y=m2_.as_1d()).coefficient()
      if(cc>result):
        result=cc
        radius = rad
    return result, radius
  ### Embedded class to run over resolutions (parallel enabled)
  class run_loop_body(object):
    def __init__(self, cg, fc, f_obs, b_range, map, selections, d_min_start,
                 d_min_end, d_min_step, nproc):
      adopt_init_args(self, locals())
      self.d_spacings = self.fc.d_spacings().data()
      self.ss = 1./flex.pow2(self.d_spacings) / 4.
      self.x=flex.double()
      self.y=flex.double()
      self.b=flex.double()
      self.radii=flex.double()
      d_mins = []
      d_min = self.d_min_start
      while d_min<self.d_min_end:
        d_mins.append(d_min)
        d_min+=self.d_min_step
      if(self.nproc>1):
        from libtbx import easy_mp
        stdout_and_results = easy_mp.pool_map(
          processes    = self.nproc,
          fixed_func   = self.run,
          args         = d_mins,
          func_wrapper = "buffer_stdout_stderr")
        for so, res in stdout_and_results:
          self.x.append(res[0])
          self.y.append(res[1])
          self.b.append(res[2])
          self.radii.append(res[3])
      else:
        for d_min in d_mins:
          res = self.run(d_min=d_min)
          #print d_min, res[1], res[2], res[3]
          self.x.append(res[0])
          self.y.append(res[1])
          self.b.append(res[2])
          self.radii.append(res[3])
    def run(self, d_min):
      sel = self.d_spacings > d_min
      fc_ = self.fc.select(sel)
      ss_ = self.ss.select(sel)
      o = mmtbx.bulk_solvent.complex_f_kb_scaled(
        self.f_obs.data().select(sel), fc_.data(),
        flex.double(self.b_range), ss_)
      fc__ = fc_.array(data = o.scaled())
      fft_map = miller.fft_map(
        crystal_gridding     = self.cg,
        fourier_coefficients = fc__)
      map_calc = fft_map.real_map_unpadded()
      cc, radius = get_cc(m1=map_calc, m2=self.map, selections=self.selections,
        radii=so.radii)
      return d_min, cc, o.b(), radius
  ###
  o = run_loop_body(cg=cg, fc=fc, f_obs=f_obs, b_range=b_range, map=map,
    selections=selections, d_min_start=d_min_start, d_min_end=d_min_end,
    d_min_step=d_min_step, nproc=nproc)
  x,y,b, radii = o.x, o.y, o.b, o.radii
  ###
  y = flex.double([round(i,4) for i in y])
  if(approximate):
    x_=[]
    y_=[]
    b_=[]
    r_=[]
    for i in xrange(y.size()):
      x_.append(x[i])
      b_.append(b[i])
      r_.append(radii[i])
      if(i>=2 and i<=y.size()-3):
        y_ave = (y[i-2]+y[i-1]+y[i]+y[i+1]+y[i+2])/5.
      else:
        y_ave = y[i]
      y_.append(y_ave)
      #print x[i],y_ave, b[i], radii[i]
    x,y,b, radii = x_,y_,b_, r_
    maxima = parabola_is_good(x=flex.double(x), y=flex.double(y),
      xmin=d_min_start, xmax=d_min_end, xstep=d_min_step)
    if(len(maxima)==0):
      return None,None,None,None
  if(approximate): cc_max = maxima[0][1]
  else: cc_max = max(y)
  r=1.e+9
  d_best, b_best, cc_best, rad_best = None,None,None,None
  for di,bi,cci,ri in zip(x,b,y,radii):
    r_ = abs(cci-cc_max)
    if(r_<r):
      r=r_
      d_best, b_best, cc_best,rad_best = di, bi, cci,ri
  return d_best, b_best, cc_best, rad_best

class run(object):
  def __init__(self, map_data, xray_structure, d_min_min=None, nproc=1):
    """
    Given map and model estimate resolution by maximizing map CC(map, model-map).
    As a by-product, also provides CC and optimal overall B-factor.
    """
    unit_cell = xray_structure.unit_cell()
    d_min_end = round(maptbx.d_min_from_map(
      map_data=map_data, unit_cell=unit_cell, resolution_factor=0.1),1)
    d_min_start = round(maptbx.d_min_from_map(
      map_data=map_data, unit_cell=unit_cell, resolution_factor=0.5),1)
    if(d_min_min is not None and d_min_start<d_min_min):
      d_min_start=d_min_min
    step = (d_min_end-d_min_start)/10.
    b_range=[0,10,20,30,40,50,60,70,80,90,100,150,200]
    if(d_min_end>10):
      b_range = b_range + [300,400,500]
    result, b_iso, cc, radius = _resolution_from_map_and_model_helper(
      map            = map_data,
      xray_structure = xray_structure,
      b_range        = b_range,
      d_min_start    = d_min_start,
      d_min_end      = d_min_end,
      d_min_step     = step,
      approximate    = False,
      nproc          = nproc)
    if(cc<0.5):
      self.d_min, self.b_iso, self.cc, self.radius = None,None,None,None
    else:
      if(d_min_min is not None and result<d_min_min):
        result=d_min_min
      scale=1
      if(result<3): scale=2
      b_range=[i for i in range(0,506,5)]
      d_min_start = round(result-result/2, 1)
      d_min_end   = round(result+result/2*scale, 1)
      if(d_min_min is not None and d_min_start<d_min_min):
        d_min_start=d_min_min
      self.d_min, self.b_iso, self.cc, self.radius = _resolution_from_map_and_model_helper(
        map              = map_data,
        xray_structure   = xray_structure,
        b_range          = b_range,
        d_min_start      = d_min_start,
        d_min_end        = d_min_end,
        d_min_step       = 0.1,
        approximate      = True,
        nproc            = nproc)
