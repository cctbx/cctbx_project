from __future__ import division, print_function
from dxtbx.model import DetectorFactory, Crystal
from iotbx.phil import parse
from dials.command_line.dials_import import  phil_scope as detector_phil_scope
from scitbx.matrix import sqr, col
from dials.array_family import flex
from matplotlib import pyplot as plt
import sys, math, copy

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
separation_plot {
  distance_min = None
    .type = float
    .help = Minimum detector distance (mm) for the distance vs. energy plot. \
            Default is half the distance in the detector phil.
  distance_max = None
    .type = float
    .help = Maximum detector distance (mm) for the distance vs. energy plot. \
            Default is twice the distance in the detector phil.
  energy_min = None
    .type = float
    .help = Minimum energy (eV) for the distance vs. energy plot. Default is \
            half the energy parameter.
  energy_max = None
    .type = float
    .help = Maximum energy (eV) for the distance vs. energy plot. Default is \
            1.5 times the energy parameter.
  n_steps = 50
    .type = int(value_min=2)
    .help = Number of grid points along each of the distance and energy axes
  target_separation = 5.5
    .type = float
    .help = Draw a dashed line at this spot separation (px), the current \
            estimate of a good spot separation. Set to None to disable.
}
""")

def spot_separations(panel, unit_cell, energy, reference_reflection, bandpass=None, verbose=False):
  """
  Compute the separation in pixels between the reference reflection and its
  neighbor (e.g. 20,0,0 and 21,0,0) along each of the three reciprocal axes.
  Axes for which the reflections cannot be observed are skipped.
  @return list of separations in pixels
  """
  wavelength = 12398.4/energy
  s0 = col((0,0,-1/wavelength))
  if bandpass:
    wavelength1 = 12398.4/(energy-(bandpass/2))
    wavelength2 = 12398.4/(energy+(bandpass/2))
  else:
    wavelength1 = wavelength2 = wavelength
  vals = []
  for axis in range(3):
    m1 = [0,0,0]; m1[axis] += reference_reflection
    m2 = [0,0,0]; m2[axis] += reference_reflection+1

    # n Lambda = 2dsin(theta)
    d = unit_cell.d(flex.miller_index([m1, m2]))
    try:
      tt_1 = math.asin(wavelength1/(2*d[0])) * 2
      tt_2 = math.asin(wavelength2/(2*d[1])) * 2
    except ValueError: # domain error if resolution is too high
      continue

    # Compute two s1 vectors
    s1_1 = s0.rotate(col((0,1,0)), -tt_1)
    s1_2 = s0.rotate(col((0,1,0)), -tt_2)

    if verbose:
      print(m1, m2, list(d), tt_1*180/math.pi, tt_2*180/math.pi)

    # Get panel intersections and compute spacing
    try:
      v1 = col(panel.get_ray_intersection_px(s1_1))
      v2 = col(panel.get_ray_intersection_px(s1_2))
    except RuntimeError: # ray does not intersect the panel plane
      continue
    vals.append((v1-v2).length())
  return vals

def inscribed_resolution(panel, energy):
  """
  Resolution (angstroms) of the largest circle inscribed on the panel, taken
  as the lowest resolution among the midpoints of the four panel edges.
  Assumes a single panel, centered detector.
  """
  s0 = col((0,0,-energy/12398.4))
  fast, slow = panel.get_image_size()
  f = fast//2; s = slow//2
  return max([panel.get_resolution_at_pixel(s0, p) for p in [(f,0),(fast,s),(f,slow),(0,s)]])

def panel_at_distance(panel, distance):
  """
  Return a copy of panel translated along the beam (z) axis such that its
  distance (along the panel normal) is the requested value in mm.
  """
  beam_dir = col((0,0,-1))
  normal = col(panel.get_normal())
  origin = col(panel.get_origin())
  # panel distance is origin.dot(normal), so solve for the shift along the beam
  shift = float(distance - panel.get_distance()) / beam_dir.dot(normal)
  new_panel = copy.deepcopy(panel)
  new_panel.set_frame(panel.get_fast_axis(), panel.get_slow_axis(), origin + beam_dir*shift)
  return new_panel

def distance_for_separation(panel, unit_cell, energy, reference_reflection, target, bandpass=None):
  """
  Find the detector distance (mm) at which the smallest spot separation equals
  target (px), by bisection. The separation increases with distance.
  @return the distance, or None if it could not be found
  """
  def separation(distance):
    vals = spot_separations(panel_at_distance(panel, distance), unit_cell, energy, reference_reflection, bandpass)
    return min(vals) if vals else None
  lo = 1.0; hi = max(2*panel.get_distance(), 2*lo)
  for i in range(30): # expand the upper bound until it exceeds the target
    val = separation(hi)
    if val is None: return None
    if val >= target: break
    hi *= 2
  else:
    return None
  val = separation(lo)
  if val is None or val > target: return None
  while hi - lo > 1e-3:
    mid = (lo+hi)/2
    if separation(mid) < target: lo = mid
    else: hi = mid
  return (lo+hi)/2

def plot_separation_vs_distance_and_energy(panel, unit_cell, params):
  """
  Plot detector distance vs. energy, colored by the smallest spot separation
  in pixels over the three reciprocal axes.
  """
  import numpy as np
  sp = params.separation_plot
  current_distance = panel.get_distance()
  d_lo = sp.distance_min if sp.distance_min is not None else current_distance * 0.5
  d_hi = sp.distance_max if sp.distance_max is not None else current_distance * 2.0
  e_lo = sp.energy_min if sp.energy_min is not None else params.energy * 0.5
  e_hi = sp.energy_max if sp.energy_max is not None else params.energy * 1.5
  assert d_lo < d_hi and e_lo < e_hi
  distances = np.linspace(d_lo, d_hi, sp.n_steps)
  energies = np.linspace(e_lo, e_hi, sp.n_steps)

  separations = np.full((len(distances), len(energies)), np.nan)
  resolutions = np.full((len(distances), len(energies)), np.nan)
  for i, distance in enumerate(distances):
    test_panel = panel_at_distance(panel, distance)
    for j, energy in enumerate(energies):
      vals = spot_separations(test_panel, unit_cell, float(energy), params.reference_reflection, params.bandpass)
      if vals:
        separations[i,j] = min(vals)
      resolutions[i,j] = inscribed_resolution(test_panel, float(energy))

  fig, ax = plt.subplots()
  mesh = ax.pcolormesh(energies, distances, separations, shading='auto', cmap='viridis')
  fig.colorbar(mesh, ax=ax, label="Smallest spot separation (px)")
  if not np.all(np.isnan(separations)):
    contours = ax.contour(energies, distances, separations, colors='white', linewidths=0.75)
    ax.clabel(contours, fmt="%.1f", fontsize=8)
    ax.plot([], [], color='white', linewidth=0.75, label="Spot separation (px)")
    if sp.target_separation is not None and \
       np.nanmin(separations) < sp.target_separation < np.nanmax(separations):
      ax.contour(energies, distances, separations, levels=[sp.target_separation],
                 colors='red', linewidths=2, linestyles='dashed')
      ax.plot([], [], color='red', linewidth=2, linestyle='dashed',
              label="Target separation: %g px"%sp.target_separation)
  res_contours = ax.contour(energies, distances, resolutions, colors='cyan',
                            linewidths=1, linestyles='dashdot')
  ax.clabel(res_contours, fmt=u"%.1f \u00c5", fontsize=8)
  ax.plot([], [], color='cyan', linewidth=1, linestyle='dashdot',
          label=u"Inscribed resolution (\u00c5)")
  # Dashed lines for the requested d_min and energy. Where these cross the
  # target separation line is the desired detector distance.
  if np.nanmin(resolutions) < params.d_min < np.nanmax(resolutions):
    ax.contour(energies, distances, resolutions, levels=[params.d_min],
               colors='cyan', linewidths=2, linestyles='dashed')
    ax.plot([], [], color='cyan', linewidth=2, linestyle='dashed',
            label=u"d_min: %g \u00c5 inscribed"%params.d_min)
  if e_lo <= params.energy <= e_hi:
    ax.axvline(params.energy, color='orange', linewidth=2, linestyle='dashed',
               label="Energy: %g eV"%params.energy)
  if e_lo <= params.energy <= e_hi and d_lo <= current_distance <= d_hi:
    ax.plot([params.energy], [current_distance], marker='*', color='red', markersize=14,
            markeredgecolor='white', linestyle='none',
            label="Current: %.0f eV, %.1f mm"%(params.energy, current_distance))
  if ax.get_legend_handles_labels()[0]:
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, fontsize=8, facecolor='gray')
  ax.set_xlabel("Energy (eV)")
  ax.set_ylabel("Detector distance (mm)")
  title = "Spot separation, reflections %d and %d"%(params.reference_reflection, params.reference_reflection+1)
  if params.bandpass:
    title += ", %g eV bandpass"%params.bandpass
  ax.set_title(title)

  # Mouseover text: energy, distance, and the resolution and spot separation
  # computed exactly at the cursor position (not interpolated from the grid)
  def format_coord(x, y):
    text = "%.1f eV, %.1f mm"%(x, y)
    if x <= 0 or y <= 0:
      return text
    cursor_panel = panel_at_distance(panel, y)
    text += u", %.2f \u00c5"%inscribed_resolution(cursor_panel, x)
    vals = spot_separations(cursor_panel, unit_cell, x, params.reference_reflection, params.bandpass)
    if vals:
      text += ", %.1f px"%min(vals)
    return text
  ax.format_coord = format_coord
  mesh.get_cursor_data = lambda event: None # avoid appending the mesh value a second time
  fig.tight_layout()
  return separations

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
  print("Reference reflections 1 and 2, resolutions, two theta (deg) 1 and 2:")
  vals = spot_separations(panel, unit_cell, energy, params.reference_reflection, params.bandpass, verbose=True)

  print("Spot separations:", vals)
  print("Smallest spot separation: %7.1f px"%(min(vals)))

  target = params.separation_plot.target_separation
  if target is not None:
    target_distance = distance_for_separation(panel, unit_cell, energy, params.reference_reflection, target, params.bandpass)
    if target_distance is None:
      print("Could not find a detector distance giving a spot separation of %g px at %g eV"%(target, energy))
    else:
      target_resolution = inscribed_resolution(panel_at_distance(panel, target_distance), energy)
      print("Detector distance for a spot separation of %g px at %g eV: %.1f mm"%(target, energy, target_distance))
      print("Inscribed resolution at that distance: %.2f angstroms (%s than the d_min of %g)"%(
        target_resolution, "better" if target_resolution <= params.d_min else "worse", params.d_min))

  # Hack for quick tests
  assert len(detector)==1
  panel = detector[0]
  fast, slow = panel.get_image_size()
  print("Inscribed resolution, assuming single panel centered detector %.3f:"% \
    inscribed_resolution(panel, energy))

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
    plot_separation_vs_distance_and_energy(detector[0], unit_cell, params)

    plt.figure()
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
