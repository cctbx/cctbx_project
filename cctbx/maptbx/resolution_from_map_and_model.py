from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx import maptbx
from libtbx import group_args
from cctbx import miller
from six.moves import range

class run_loop(object):
  def __init__(self, fc, f_obs, b_iso, d_mins, d_spacings, ss):
    result = maptbx.cc_complex_complex(
      f_1        = f_obs.data(),
      f_2        = fc.data(),
      d_spacings = d_spacings,
      ss         = ss,
      d_mins     = d_mins,
      b_iso      = b_iso)
    self.d_min, self.cc = result[0], result[1]

def get_trial_resolutions(d_start):
  if(  d_start<5):  d_end = 10.
  elif(d_start<10): d_end = 20.
  elif(d_start<20): d_end = 40.
  else:             d_end = 100.
  l1 = [i/10. for i in range(int(d_start*10),30)]
  l2 = [i/100. for i in range(300,610,10)]
  l3 = [i/10. for i in range(65,105,1)]
  l4 = [i*1. for i in range(11,100,1)]
  l = flex.double(l1+l2+l3+l4)
  s  = l>=d_start
  s &= l<=d_end
  return l.select(s)

def run_at_b(b, f_map, f_calc=None, xray_structure=None):
  d_mins = get_trial_resolutions(d_start = f_map.d_min())
  if(f_calc is None):
    assert xray_structure is not None
    xrs = xray_structure.deep_copy_scatterers().set_b_iso(value=0)
    f_calc = f_map.structure_factors_from_scatterers(
      xray_structure = xrs).f_calc()
  else:
    assert xray_structure is None
  d_spacings = f_calc.d_spacings().data()
  ss = 1./flex.pow2(d_spacings) / 4.
  result = run_loop(fc=f_calc, f_obs=f_map, b_iso=b, d_mins=d_mins,
    d_spacings=d_spacings, ss=ss)
  return result

def tmp(f):
  import math
  f.setup_binner(reflections_per_bin = 100)
  data = abs(f).data()
  ds = f.d_spacings().data()
  for i_bin in f.binner().range_used():
    sel = f.binner().selection(i_bin)
    m = flex.mean(data.select(sel))
    dm = flex.mean(ds.select(sel))
    print("%5d %10.6f %10.6f"%(i_bin, dm,  math.log(m)))

def find_b(fo, fc):
  # TODO: Need to use linear
  #tmp(f=fo)
  #fo=fo.resolution_filter(d_min=5)
  #fo, fc, = fo.common_sets(fc)
  cc = -999
  b=None
  ss = 1./flex.pow2(fc.d_spacings().data()) / 4.
  data = fc.data()
  for b_ in range(-500,500,5):
    sc = flex.exp(-b_*ss)
    fc_ = fc.customized_copy(data = data*sc)
    cc_ = fo.map_correlation(other = fc_)
    if(cc_>cc):
      cc = cc_
      b = b_
  return b

def run(xray_structure, f_map=None, map_data=None, d_fsc_model=None):
  assert [f_map, map_data].count(None) == 1
  xrs = xray_structure.deep_copy_scatterers().set_b_iso(value=0)
  if(f_map is None):
    f_map = miller.structure_factor_box_from_map(
      map              = map_data,
      crystal_symmetry = xray_structure.crystal_symmetry())
  fc = f_map.structure_factors_from_scatterers(xray_structure = xrs).f_calc()
  d_model_b0 = run_at_b(b=0, f_map=f_map, f_calc=fc).d_min
  del xrs
  if(d_fsc_model is None):
    d_fsc_model = fc.d_min_from_fsc(other=f_map, fsc_cutoff=0).d_min
  fo = f_map.resolution_filter(d_min = d_fsc_model)
  fo, fc, = fo.common_sets(fc)
  b = find_b(fo=fo.deep_copy(), fc=fc.deep_copy())
  o = run_at_b(b = b, f_map = fo, f_calc = fc)
  return group_args(d_min = o.d_min, b_iso = b, d_model_b0 = d_model_b0,
    d_fsc_model=d_fsc_model)
