from __future__ import absolute_import, division, print_function
from cctbx import maptbx
from cctbx import miller
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from libtbx import group_args
import time

def get_map_data(xrs, d_min):
  cg = maptbx.crystal_gridding(
    unit_cell         = xrs.unit_cell(),
    space_group_info  = xrs.space_group_info(),
    step              = 0.1)
  fc = xrs.structure_factors(d_min = d_min, algorithm = "direct").f_calc()
  fft_map = miller.fft_map(
    crystal_gridding     = cg,
    fourier_coefficients = fc,
    f_000                = xrs.f_000())
  fft_map.apply_volume_scaling()
  return fft_map.real_map_unpadded()

def exercise_00():
  o = maptbx.atom_curves(scattering_type="C")
  result = o.image(d_min = 3.0, b_iso = 50, radius_step=0.01)
  assert approx_equal(result.radius, 2.59, 0.01)
  #
  xrs = o.get_xray_structure(box=10, b=50)
  map_data = get_map_data(xrs=xrs, d_min=3.0)
  image=flex.double()
  dim=xrs.unit_cell().parameters()[0]
  for i, r in enumerate(result.radii):
    mv = map_data.eight_point_interpolation([r/dim,0,0])
    image.append(mv)
  cc = flex.linear_correlation(x=image, y=result.image_values).coefficient()
  assert cc > 0.99
  num = flex.sum(flex.abs(flex.abs(image)-flex.abs(result.image_values)))
  den = flex.sum(flex.abs(flex.abs(image)+flex.abs(result.image_values)))
  r = 100.*2.*num/den
  assert r < 6.

def _exercise_01(d_min=2.5, b=20, box=10, step=0.025, cut=True):
  o = maptbx.atom_curves(scattering_type="C", scattering_table="wk1995")
  cs = o.get_xray_structure(box=box,b=b).crystal_symmetry()
  #
  n = int(box/step)
  complete_set = miller.structure_factor_box_from_map(
    crystal_symmetry = cs,
    n_real           = [n,n,n],
    anomalous_flag   = False,
    include_000      = True)
  if(cut): complete_set = complete_set.resolution_filter(d_min=d_min)
  #
  r1 = o.image_from_miller_indices(miller_indices=complete_set.indices(),
    b_iso=b, uc = cs.unit_cell(), radius_max=box/box, radius_step=step/box)
  r2 = o.image_from_3d(box=box, b=b, step=step, unit_cell=complete_set.unit_cell(),
    space_group_info = complete_set.space_group_info(),
    miller_array = complete_set)
  r3 = o.image(d_min = d_min, b_iso = b, radius_max=box, radius_step=step)
  r4 = o.exact_density(b_iso=b, radius_max=box, radius_step=step)
  return group_args(r1=r1, r2=r2, r3=r3, r4=r4)

def exercise_01():
  def rf(x,y):
    n = flex.sum(flex.abs(flex.abs(x)-flex.abs(y)))
    d = flex.sum(flex.abs(flex.abs(x)+flex.abs(y))/2.)
    return n/d*100.
  r = _exercise_01(d_min=2.5, b=20, box=10, step=0.025, cut=True)
  i_half = int(r.r1.image_values.size()/3)
  r1 = r.r1.image_values[:i_half]
  r2 = r.r2.image_values[:i_half]
  r3 = r.r3.image_values[:i_half]
  assert rf(r1,r2) < 2.5
  assert rf(r1,r3) < 2
  assert flex.linear_correlation(r1,r2).coefficient() > 0.999
  assert flex.linear_correlation(r1,r3).coefficient() > 0.999

def exercise_02():
  def rf(x,y):
    n = flex.sum(flex.abs(flex.abs(x)-flex.abs(y)))
    d = flex.sum(flex.abs(flex.abs(x)+flex.abs(y))/2.)
    return n/d*100.
  r = _exercise_01(d_min=2.5, b=10, box=7, step=0.1, cut=False)
  i_half = int(r.r1.image_values.size()/3)
  r1 = r.r1.image_values[:i_half]
  r2 = r.r2.image_values[:i_half]
  r3 = r.r3.image_values[:i_half]
  r4 = r.r4.density[:i_half]
  assert rf(r1,r4) < 2
  assert flex.linear_correlation(r1,r4).coefficient() > 0.999
  assert rf(r1,r2) < 2
  assert flex.linear_correlation(r1,r2).coefficient() > 0.999

if (__name__ == "__main__"):
  t0 = time.time()
  exercise_00()
  exercise_01()
  exercise_02()
  print("Time: %6.3f"%(time.time()-t0))
  print("OK")
