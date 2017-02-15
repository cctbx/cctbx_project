from __future__ import division
from cctbx.array_family import flex
from libtbx import adopt_init_args
from libtbx.test_utils import approx_equal
from cctbx import maptbx
from cctbx import miller
import scitbx.math.curve_fitting

def parabola_is_good(x, y, assert_concave_up, use_longest_slope_criteria=False):
  if(assert_concave_up):
    fit = scitbx.math.curve_fitting.univariate_polynomial_fit(x_obs=x, y_obs=y,
      degree=2, min_iterations=50, number_of_cycles=10)
    a0,a1,a2 = fit.params
    if(a2>=0.): return [] # concave down
  maxima = []
  for i in xrange(y.size()):
    yi = y[i]
    cntr=0
    yp,ym=[],[]
    offsets = [1,] #XXX[1,2,3]
    for j in offsets:
      if(i-j>=0 and y[i-j]<=yi and not y[i-j] in ym):
        cntr+=1
        ym.append(y[i-j])
      if(i+j<=y.size()-1 and y[i+j]<=yi and not y[i+j] in yp):
        cntr+=1
        yp.append(y[i+j])
    #XXX if(cntr==len(offsets)+len(offsets)): maxima.append([x[i],yi,i])
    if(cntr==2 and len(offsets)==len(offsets)): maxima.append([x[i],yi,i])
  # Remove plateau
  tmp=[]
  tmp2 = []
  for m in maxima:
    if(not m[1] in tmp):
      tmp.append(m[1])
      tmp2.append(m)
  maxima = tmp2
  # Choose peaks with longest slope
  if(use_longest_slope_criteria):
    if(len(maxima)>1):
      maxima_plus = []
      for m in maxima:
        i = m[2]
        cntr=0
        for j in range(1,100):
          if(i-j>=0 and y[i-j]<=y[i-j+1]):
            cntr+=1
          else: break
          if(i+j<=y.size()-1 and y[i+j]<=y[i+j-1]):
            cntr+=1
          else: break
        maxima_plus.append([m,cntr])
      cntr_max = -1
      maximum = None
      for mp in maxima_plus:
        if(mp[1]>cntr_max):
          cntr_max = mp[1]
          maximum = mp[0]
      maxima = [maximum]
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
      thorough=False):
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
  d_min_start = fc.d_min()
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
          self.x.append(res[0])
          self.y.append(res[1])
          self.b.append(res[2])
          self.radii.append(res[3])
    def y_smooth(self):
      y = self.y
      y_=[]
      for i in xrange(y.size()):
        if(i>=1 and i<=y.size()-2):
          y_ave = (y[i-1]+y[i]+y[i+1])/3.
        else:
          y_ave = y[i]
        #if(i>=2 and i<=y.size()-3):
        #  y_ave = (y[i-2]+y[i-1]+y[i]+y[i+1]+y[i+2])/5.
        #elif(i>=1 and i<=y.size()-2):
        #  y_ave = (y[i-1]+y[i]+y[i+1])/3.
        #else:
        #  y_ave = y[i]
        y_.append(y_ave)
      return flex.double(y_)
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
  if(thorough):
    o = run_loop_body(cg=cg, fc=fc, f_obs=f_obs, b_range=b_range, map=map,
      selections=selections, d_min_start=d_min_start,
      d_min_end=f_obs.d_min()+1.e-6,
      d_min_step=d_min_step, nproc=1)
    junk,junk,b, junk = o.x, o.y, o.b, o.radii
    assert o.b.size()==1
    o1 = run_loop_body(cg=cg, fc=fc, f_obs=f_obs, b_range=o.b, map=map,
      selections=selections, d_min_start=d_min_start, d_min_end=d_min_end,
      d_min_step=d_min_step, nproc=nproc)
    maxima1 = parabola_is_good(x=o1.x, y=o1.y_smooth(), assert_concave_up=True,
      use_longest_slope_criteria=True)
    #print "maxima1",maxima1
    if(len(maxima1)!=1): return None,None,None,None # Exactly one peak expected
    o2 = run_loop_body(cg=cg, fc=fc, f_obs=f_obs, b_range=b_range, map=map,
      selections=selections, d_min_start=d_min_start,
      d_min_end=max(maxima1[0][0]+2., d_min_end),
      d_min_step=d_min_step, nproc=nproc)
    maxima2 = parabola_is_good(x=o2.x, y=o2.y_smooth(), assert_concave_up=False)
    #print "maxima2",maxima2
    if(len(maxima2)==0): return None,None,None,None # At least one peak expected
    # Match maxima2 against maxima1 to find the closest peak to maxima1
    d1 = maxima1[0][0]
    dist=1.e+9
    m_best = None
    for m in maxima2:
      dist_ = abs(d1-m[0])
      if(dist_<dist):
        dist=dist_
        m_best = m[:]
    #assert approx_equal(o1.x, o2.x)
    #for d_min, cc1,b1,r1 in zip(o1.x, o1.y_smooth(),o1.b,o1.radii):
    #  print "%4.1f %8.6f %3.0f %5.2f"%(d_min, cc1,b1,r1)
    #print
    #for d_min, cc1,b1,r1 in zip(o2.x, o2.y_smooth(),o2.b,o2.radii):
    #  print "%4.1f %8.6f %3.0f %5.2f"%(d_min, cc1,b1,r1)
    #
    return m_best[0], o2.b[m_best[2]], o2.y[m_best[2]], o2.radii[m_best[2]]
  else:
    o2 = run_loop_body(cg=cg, fc=fc, f_obs=f_obs, b_range=b_range, map=map,
      selections=selections, d_min_start=d_min_start, d_min_end=d_min_end,
      d_min_step=d_min_step, nproc=nproc)
    cc_max = flex.max(o2.y)
  # Finalize
  r=1.e+9
  d_best, b_best, cc_best, rad_best = None,None,None,None
  for di,bi,cci,ri in zip(o2.x, o2.b, o2.y, o2.radii):
    r_ = abs(cci-cc_max)
    if(r_<r):
      r=r_
      d_best, b_best, cc_best,rad_best = di, bi, cci,ri
  return d_best, b_best, cc_best, rad_best

def strip_model(map_data, xray_structure, pdb_hierarchy, d_min, radius, b_iso):
  from mmtbx.maps import correlation
  size_start = xray_structure.scatterers().size()
  xray_structure = xray_structure.set_b_iso(value=b_iso)
  crystal_gridding = maptbx.crystal_gridding(
    unit_cell             = xray_structure.unit_cell(),
    space_group_info      = xray_structure.space_group_info(),
    pre_determined_n_real = map_data.accessor().all(),
    symmetry_flags        = maptbx.use_space_group_symmetry)
  f_calc = xray_structure.structure_factors(d_min=d_min).f_calc()
  fft_map = miller.fft_map(
    crystal_gridding     = crystal_gridding,
    fourier_coefficients = f_calc)
  map_calc = fft_map.real_map_unpadded()
  selection = flex.size_t()
  for rg in pdb_hierarchy.residue_groups():
    sel = rg.atoms().extract_i_seq()
    xyz = rg.atoms().extract_xyz()
    cc = correlation.from_map_map_atoms(
      map_1=map_data, map_2=map_calc, sites_cart=xyz,
      unit_cell=f_calc.unit_cell(), radius=radius)
    if(cc>0.5):
      selection.extend(sel)
  xray_structure = xray_structure.select(selection)
  size_final = xray_structure.scatterers().size()
  if(size_final*100./size_start < 50):
    return None
  else:
    return xray_structure

class run(object):
  def __init__(self, map_data, xray_structure, pdb_hierarchy, d_min_min=None,
                     nproc=1):
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
      thorough       = False,
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
      self.d_min, self.b_iso, self.cc, self.radius = \
        _resolution_from_map_and_model_helper(
          map            = map_data,
          xray_structure = xray_structure,
          b_range        = b_range,
          d_min_start    = d_min_start,
          d_min_end      = d_min_end,
          d_min_step     = 0.1,
          thorough       = True,
          nproc          = nproc)
      xray_structure = strip_model(
        map_data       = map_data,
        xray_structure = xray_structure,
        pdb_hierarchy  = pdb_hierarchy,
        d_min          = self.d_min,
        radius         = self.radius,
        b_iso          = self.b_iso)
      if(xray_structure is not None):
        d_min, b_iso, cc, radius = \
          _resolution_from_map_and_model_helper(
            map            = map_data,
            xray_structure = xray_structure,
            b_range        = b_range,
            d_min_start    = d_min_start,
            d_min_end      = d_min_end,
            d_min_step     = 0.1,
            thorough       = True,
            nproc          = nproc)
        if(d_min is not None):
          self.d_min, self.b_iso, self.cc, self.radius = \
            d_min, b_iso, cc, radius
