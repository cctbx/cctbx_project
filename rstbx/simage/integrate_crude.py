from __future__ import absolute_import, division, print_function
from six.moves import range
def predict_spot_positions(
      work_params,
      miller_indices,
      unit_cell,
      crystal_rotation):
  from rstbx.simage import image_simple
  image = image_simple(
    store_spots=True,
    store_miller_index_i_seqs=True).compute(
      unit_cell=unit_cell,
      miller_indices=miller_indices,
      spot_intensity_factors=None,
      crystal_rotation_matrix=crystal_rotation,
      ewald_radius=1/work_params.wavelength,
      ewald_proximity=work_params.ewald_proximity,
      signal_max=work_params.signal_max,
      detector_distance=work_params.detector.distance,
      detector_size=work_params.detector.size,
      detector_pixels=work_params.detector.pixels,
      point_spread=work_params.point_spread,
      gaussian_falloff_scale=work_params.gaussian_falloff_scale)
  return image.spots, image.miller_index_i_seqs

def sum_pixels(pixels, point_spread_inner, point_spread_outer, center):
  circle_radius_sq_inner = point_spread_inner**2 / 4
  circle_radius_sq_outer = point_spread_outer**2 / 4
  dpx, dpy = pixels.focus()
  pxf, pyf, _ = center
  pxi = int(pxf)
  pyi = int(pyf)
  pxb = pxi - point_spread_outer // 2
  pyb = pyi - point_spread_outer // 2
  if (point_spread_outer % 2 == 0):
    if (pxf - pxi > 0.5): pxb += 1
    if (pyf - pyi > 0.5): pyb += 1
  n_outer = 0
  n_inner = 0
  sum_outer = 0
  sum_inner = 0
  for i in range(0, point_spread_outer+1):
    pi = pxb + i
    if (pi < 0 or pi >= dpx): return 0
    for j in range(0, point_spread_outer+1):
      pj = pyb + j
      if (pj < 0 or pj >= dpy): return 0
      pcx = (pi + 0.5) - pxf
      pcy = (pj + 0.5) - pyf
      pc_sq = pcx*pcx + pcy*pcy
      if (pc_sq > circle_radius_sq_outer): continue
      c = pixels[pi,pj]
      if (pc_sq > circle_radius_sq_inner):
        n_outer += 1
        sum_outer += c
      else:
        n_inner += 1
        sum_inner += c
  assert n_outer != 0
  return max(0, sum_inner - sum_outer * n_inner / n_outer)

def collect_spot_intensities(
      pixels,
      spot_positions,
      point_spread_inner,
      point_spread_outer):
  from scitbx.array_family import flex
  raw_sums = flex.double()
  raw_sums.reserve(spot_positions.size())
  for position in spot_positions:
    raw_sums.append(sum_pixels(
      pixels=pixels,
      point_spread_inner=point_spread_inner,
      point_spread_outer=point_spread_outer,
      center=position))
  return raw_sums
