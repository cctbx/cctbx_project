from __future__ import absolute_import, division, print_function
from cctbx.development import random_structure
from cctbx.sgtbx import space_group_info
import boost_adaptbx.boost.python as bp
from six.moves import range
ext = bp.import_ext("cctbx_asymmetric_map_ext")
from cctbx_asymmetric_map_ext import *
from cctbx.array_family import flex
from cctbx import miller
from cctbx import maptbx  # import dependency
from libtbx.test_utils import approx_equal
import time, random

if(1):
  random.seed(0)
  flex.set_random_seed(0)

def get_map(xrs,ms):
  fc = ms.structure_factors_from_scatterers(
    xray_structure=xrs, algorithm="direct").f_calc()
  fft_map = fc.fft_map(resolution_factor=1./10)
  fft_map.apply_sigma_scaling()
  return fft_map.real_map_unpadded()

def get_err(x,y):
  return abs(x-y)/abs(x+y)*2*100

def run_group(symbol, err_l, err_q, err_t, err_lt1, err_lt2, err_lt12, err_qt1):
  group = space_group_info(symbol)
  elements = ('C', )*10
  xrs = random_structure.xray_structure(
    space_group_info       = group,
    volume_per_atom        = 25.,
    general_positions_only = False,
    elements               = elements,
    min_distance           = 1.0)
  ms = miller.set(xrs.crystal_symmetry(),
    flex.miller_index(((1,2,3),)), False).complete_set(d_min=2)
  m1 = get_map(xrs=xrs, ms=ms)
  #
  sites_frac = xrs.sites_frac()
  for i, shift in enumerate([[-2,-3,-4],[2,3,4]]):
    sf_sh = sites_frac+shift
    xrs.set_sites_frac(sites_frac=sf_sh)
    #
    m2 = get_map(xrs=xrs, ms=ms)
    #
    if(i==0): assert min(sf_sh.min())<0.
    else:     assert min(sf_sh.min())>1.
    if(i==0): assert max(sf_sh.min())<0.
    else:     assert max(sf_sh.min())>1.
    #
    for sf in sf_sh:
      l1 = abs(m1.eight_point_interpolation(sf))
      l2 = abs(m2.eight_point_interpolation(sf))
      err_l.append(get_err(l1,l2))
      #
      q1 = abs(m1.quadratic_interpolation_with_gradients(sf,[1,1,1])[0])
      q2 = abs(m2.quadratic_interpolation_with_gradients(sf,[1,1,1])[0])
      err_q.append(get_err(q1,q2))
      #
      t1 = abs(m1.tricubic_interpolation_with_gradients(sf,[1,1,1])[0])
      t2 = abs(m2.tricubic_interpolation_with_gradients(sf,[1,1,1])[0])
      err_t.append(get_err(t1,t2))
      #
      err_lt1.append( get_err(l1,t1))
      err_lt2.append( get_err(l2,t2))
      err_lt12.append(get_err(l1,t2))
      #
      err_qt1.append(get_err(l1,q1))

def run():
  err_l    = flex.double()
  err_q    = flex.double()
  err_t    = flex.double()
  err_lt1  = flex.double()
  err_lt2  = flex.double()
  err_lt12 = flex.double()
  err_qt1  = flex.double()
  for i in range(1,231):
    run_group(i, err_l, err_q, err_t, err_lt1, err_lt2, err_lt12, err_qt1)
  assert approx_equal([flex.mean(err_l), flex.mean(err_q), flex.mean(err_t)],
    [0,0,0], 1.e-3)
  assert flex.mean(err_lt1)  < 2.
  assert flex.mean(err_lt2)  < 2.
  assert flex.mean(err_lt12) < 2.
  assert flex.mean(err_qt1)  < 2.

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("Time: %6.2f"%(time.time()-t0))
  print("OK")
