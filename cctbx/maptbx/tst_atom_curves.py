from __future__ import division
from __future__ import print_function
from cctbx import maptbx
from cctbx import miller
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal

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

def run():
  o = maptbx.atom_curves(scattering_type="C")
  result = o.image(d_min = 3.0, b_iso = 50, radius_step=0.01)
  assert approx_equal(result.radius, 2.59, 0.01)
  #
  o.xray_structure.set_b_iso(value=50)
  map_data = get_map_data(xrs=o.xray_structure, d_min=3.0)
  image=flex.double()
  dim=o.xray_structure.unit_cell().parameters()[0]
  for i, r in enumerate(result.radii):
    mv = map_data.eight_point_interpolation([r/dim,0,0])
    image.append(mv)
  cc = flex.linear_correlation(x=image, y=result.image_values).coefficient()
  assert cc > 0.99
  num = flex.sum(flex.abs(flex.abs(image)-flex.abs(result.image_values)))
  den = flex.sum(flex.abs(flex.abs(image)+flex.abs(result.image_values)))
  r = 100.*2.*num/den
  assert r < 6.

if (__name__ == "__main__"):
  run()
  print("OK")
