from __future__ import division, print_function
from dxtbx.model import DetectorFactory, Crystal
from iotbx.phil import parse
from dials.command_line.dials_import import  phil_scope as detector_phil_scope
from scitbx.matrix import sqr, col
from dials.array_family import flex
from matplotlib import pyplot as plt
import sys, math

"""
This script uses several approaches to simulate a detector and crystal and predict
spot separation and resolution extent.

Invocation:
libtbx.python show_spot_separation.py detector.phil crystal.phil [parameters]

Where detector.phil looks like (for example):
geometry {
  detector {
    panel
    {
      id = 0
      name = Jungfrau Monolithic
      type = PAD
      gain = 1.47
      pixel_size = (0.075, 0.075)
      image_size = (4096, 4096)
      trusted_range = -1, 65535
      thickness = None
      material = Si
      fast_axis = 1,0,0
      slow_axis = 0,-1,0
      origin = -153.6, 153.6, -250
    }
  }
}

And crystal looks like (for example):
unit_cell = 77, 77, 37, 90, 90, 90
space_group = P43212

detector.phil describes a single panel according to the dxtbx specifications while
crystal.phil describes a unit cell and space group.

"""

crystal_scope = parse("""
space_group = None
  .type = space_group
  .help = "Target space group."
unit_cell = None
  .type = unit_cell
  .help = "Target unit cell."
""")

phil_scope = parse("""
energy = 9500.0
  .type = float
d_min = 2.0
  .type = float
reference_reflection = 20
  .type = int
bandpass = None
  .type = float
  .help = Full width
show_plots=True
  .type = bool
""")

def run(args):
  # read in phil files (detector and crystal)
  d_params = detector_phil_scope.fetch(parse(file_name = args[0])).extract()
  detector = DetectorFactory.from_phil(d_params.geometry)
  print(detector)
  assert len(detector) == 1; panel = detector[0]

  c_params = crystal_scope.fetch(parse(file_name = args[1])).extract()
  unit_cell = c_params.unit_cell
  sg_info = c_params.space_group
  a = sqr(unit_cell.orthogonalization_matrix()) * col((1,0,0))
  b = sqr(unit_cell.orthogonalization_matrix()) * col((0,1,0))
  c = sqr(unit_cell.orthogonalization_matrix()) * col((0,0,1))
  crystal = Crystal(a,b,c,sg_info.group())
  print(crystal)

  # load additional parameters
  user_phil = []
  for arg in args[2:]:
    user_phil.append(parse(arg))
  params = phil_scope.fetch(sources=user_phil).extract()

  energy = float(params.energy)
  wavelength = 12398.4/energy
  s0 = col((0,0,-1/wavelength))
  if params.bandpass is not None:
    wavelength1 = 12398.4/(energy-(params.bandpass/2))
    wavelength2 = 12398.4/(energy+(params.bandpass/2))
  vals = []
  print("Reference reflections 1 and 2, resolutions, two theta (deg) 1 and 2:")
  for axis in range(3):
    m1 = [0,0,0]; m1[axis] += params.reference_reflection
    m2 = [0,0,0]; m2[axis] += params.reference_reflection+1

    # n Lambda = 2dsin(theta)
    d = unit_cell.d(flex.miller_index([m1, m2]))
    try:
      if params.bandpass:
        tt_1 = math.asin(wavelength1/(2*d[0])) * 2
        tt_2 = math.asin(wavelength2/(2*d[1])) * 2
      else:
        tt_1 = math.asin(wavelength/(2*d[0])) * 2
        tt_2 = math.asin(wavelength/(2*d[1])) * 2
    except ValueError: # domain error if resolution is too high
      continue

    # Compute two s1 vectors
    s1_1 = s0.rotate(col((0,1,0)), -tt_1)
    s1_2 = s0.rotate(col((0,1,0)), -tt_2)

    print(m1, m2, list(d), tt_1*180/math.pi, tt_2*180/math.pi)

    # Get panel intersections and compute spacing
    v1 = col(panel.get_ray_intersection_px(s1_1))
    v2 = col(panel.get_ray_intersection_px(s1_2))
    vals.append((v1-v2).length())

  print("Spot separations:", vals)
  print("Smallest spot separation: %7.1f px"%(min(vals)))

  # Hack for quick tests
  assert len(detector)==1
  panel = detector[0]
  fast, slow = panel.get_image_size()
  f = fast//2; s = slow//2
  print("Inscribed resolution, assuming single panel centered detector %.3f:"% \
    min([panel.get_resolution_at_pixel(s0, p) for p in [(f,0),(fast,s),(f,slow),(0,s)]]))

  print("Computing pixel resolutions...")
  resolutions = []
  for panel in detector:
    fast, slow = panel.get_image_size()
    resolutions.append(flex.double(flex.grid(slow, fast)))

    for s in range(slow):
      for f in range(fast):
        resolutions[-1][s,f] = panel.get_resolution_at_pixel(s0, (f, s))

  print("Done")

  d_max = params.d_min * 1.1
  in_range = 0; total = 0
  for r in resolutions:
    in_range += len(r.as_1d().select((r.as_1d()>=params.d_min) & (r.as_1d() <= d_max)))
    total += len(r)

  print("%d of %d pixels between %.2f and %.2f angstroms (%.1f%%)"%(in_range, total, params.d_min, d_max, 100*in_range/total))
  two_theta_d_min = math.asin(wavelength/(2*params.d_min))*2
  d_min_radius_mm = math.tan(two_theta_d_min)*panel.get_distance()
  d_min_radius_px = d_min_radius_mm / panel.get_pixel_size()[0]
  possible_coverage_d_min = math.pi*d_min_radius_px**2
  two_theta_d_max = math.asin(wavelength/(2*d_max))*2
  d_max_radius_mm = math.tan(two_theta_d_max)*panel.get_distance()
  d_max_radius_px = d_max_radius_mm / panel.get_pixel_size()[0]
  possible_coverage_d_max = math.pi*d_max_radius_px**2
  possible_coverage = possible_coverage_d_min - possible_coverage_d_max
  print("Ideal detector would include %d pixels between %.2f-%.2f angstroms"%(possible_coverage, params.d_min, d_max))
  print("Coverage: %d/%d = %.1f%%"%(in_range, possible_coverage, 100*in_range/possible_coverage))

  two_theta_values = flex.double()
  step = (two_theta_d_max - two_theta_d_min)/10
  for i in range(11):
    two_theta_values.append(two_theta_d_max + (step*i))
  s0 = flex.vec3_double(len(two_theta_values), (0,0,-1))
  v = s0.rotate_around_origin((0,1,0), two_theta_values)
  all_v = flex.vec3_double()
  for i in range(720):
    i = i/2
    all_v.extend(v.rotate_around_origin((0,0,-1), i*math.pi/180))

  intersecting_rays = flex.bool()

  for i in range(len(all_v)):
    try:
      panel, mm = detector.get_ray_intersection(all_v[i])
    except RuntimeError:
      intersecting_rays.append(False)
    else:
      intersecting_rays.append(panel >=0 and panel < len(detector))

  print("%d rays out of %d projected between %f and %f intersected the detector (%.1f%%)"% \
    (intersecting_rays.count(True), len(intersecting_rays), params.d_min, d_max, intersecting_rays.count(True)*100/len(intersecting_rays)))

  resolutions[0].set_selected(resolutions[0] > 50, 50)
  if params.show_plots:
    plt.imshow(resolutions[0].as_numpy_array(), cmap='gray')
    plt.colorbar()

    plt.figure()
    r = resolutions[0]
    sel = (r.as_1d()>=params.d_min) & (r.as_1d() <= d_max)
    r.as_1d().set_selected(~sel, 0)
    plt.imshow(r.as_numpy_array(), cmap='gray')
    plt.colorbar()

    plt.show()

if __name__ == "__main__":
  run(sys.argv[1:])
