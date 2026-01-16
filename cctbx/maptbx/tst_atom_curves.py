from __future__ import absolute_import, division, print_function
from cctbx import maptbx
from cctbx import miller
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from libtbx import group_args
import time, math


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
  assert approx_equal(result.radius, 2.59, 0.001)
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
  assert r < 6.5, r

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

def get_rel_err(m1, m2):
  mean_max = (m1[0]+m2[0])/2
  diff = flex.abs(m1-m2)
  rel_err = diff/mean_max*100.
  return flex.mean(rel_err), flex.max(rel_err)

def exercise_03(d_min = 1.0, n_grid = 2000, dist_max = 5.0):
  #
  # Image of S (fast), (8-9)
  #
  ff_S = [6.372157,  5.154568, 1.473732,  1.635073,  1.209372, 0.154722,
          1.514347, 22.092527, 0.061373, 55.445175,  0.646925, 0.000000]
  image1, radii = maptbx.atom_image_fast(
     ff_packed = ff_S,
     d_min     = d_min,
     n_grid    = n_grid,
     dist_max  = dist_max)
  assert approx_equal(image1.size(), n_grid)
  assert approx_equal(image1.size(), radii.size())
  assert approx_equal(radii[0], 0.0)
  assert approx_equal(radii[-1], 5.0)
  # control, (8-9) (same call, using different interface)
  o = maptbx.atom_curves(scattering_type="S", scattering_table="wk1995")
  im = o.image(
    d_min = d_min,
    b_iso = 0,
    radii = radii,
    fast  = True)
  assert approx_equal(image1, im.image_values)
  assert approx_equal(radii, im.radii)
  emean, emax = get_rel_err(m1=image1, m2=im.image_values)
  assert approx_equal(emean, 0)
  assert approx_equal(emax, 0)
  #
  # Image via integration (slow), (8-9)
  #
  im = o.image(
      d_min = d_min,
      b_iso = 0,
      radii = radii,
      fast  = False)
  emean, emax = get_rel_err(m1=image1, m2=im.image_values)
  assert emean < 1.e-5, emean
  assert emax  < 1.e-5, emax
  assert approx_equal(image1, im.image_values)
  assert approx_equal(radii, im.radii)
  #
  # BCR approximation, (1-2)
  #
  from cctbx.maptbx.bcr import qmap
  from cctbx.maptbx.bcr import bcr
  t = qmap.load_table(element="S", table="wk1995")
  d = t["1.0"]
  B = d["B"]
  C = d["C"]
  R = d["R"]
  approx = bcr.curve(B=B, C=C, R=R, radii=radii, b_iso=0)
  emean, emax = get_rel_err(m1=image1, m2=approx)
  print("emean, emax:", emean, emax)
  #
  with open("tst_image_and_approx_S_dmin1A_rmax5A_ngrid2000.log","w") as fo:
    print("  Distance Image(8-9)    Approx (1-2)", file=fo)
    for r, im, ap in zip(radii, image1, approx):
      print("%8.4f %13.8f %13.8f"%(r, im, ap), file=fo)

if (__name__ == "__main__"):
  t0 = time.time()
  exercise_00()
  exercise_01()
  exercise_02()
  exercise_03()
  print("Time: %6.3f"%(time.time()-t0))
  print("OK")
