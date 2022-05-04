# -*- coding: utf-8 -*-
#!/usr/bin/env python
#
# detector_congruence.py
#
#  Copyright (C) 2016 Lawrence Berkeley National Laboratory (LBNL)
#
#  Author: Aaron Brewster
#
#  This code is distributed under the X license, a copy of which is
#  included in the root directory of this package.
#
# LIBTBX_SET_DISPATCHER_NAME dev.cctbx.xfel.detector_residuals
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.detector_residuals
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1
#
from __future__ import absolute_import, division, print_function
from six.moves import range
from dials.array_family import flex
from dials.util import show_mail_on_error
from scitbx.matrix import col
from matplotlib import pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.colors import Normalize
from matplotlib import cm
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from libtbx.phil import parse
import math
from six.moves import zip

help_message = '''

This program is used to calculate statisical measurements of consistency
between observed and predicted reflections

Example:

  dev.cctbx.xfel.detector_residuals experiment.expt reflections.refl
'''

# Create the phil parameters
phil_scope = parse('''
dot_size = 10
  .type = int
  .help = Size of dots in detector plots
panel_numbers = True
  .type = bool
  .help = Whether to show panel numbers on each panel
verbose = True
  .type = bool
  .help = Whether to print statistics
repredict_input_reflections = True
  .type = bool
  .help = Whether to use the input models to repredict reflection positions \
          prior to making plots
residuals {
  plot_max=0.3
    .type = float
    .help = Maximum residual value to be shown in the detector plot
  histogram_max=None
    .type = float
    .help = Maximum x value to be used computing the histogram
  histogram_xmax=None
    .type = float
    .help = Maximum x value to be shown in the histogram
  histogram_ymax=None
    .type = float
    .help = Maximum y value to be shown in the histogram
  exclude_outliers=None
    .type = bool
    .deprecated = True
    .help = Deprecated. Use residuals.exclude_outliers_from_refinement instead
  exclude_outliers_from_refinement=True
    .type = bool
    .help = Whether to exclude outliers while computing statistics
  i_sigi_cutoff=None
    .type = float
    .help = Minimum I/sigI filter for RMSD plots
  recompute_outliers = False
    .type = bool
    .help = If True, use sauter_poon to recompute outliers and remove them
  mcd_filter {
    enable = False
      .type = bool
      .help = on the DeltaPsi plot, apply a 3-feature MCD filter
    mahalanobis_distance = 7
      .type = float(value_min=0.)
      .help = cutoff level expressed as a multivariate Gaussian std dev
    keep = *inliers outliers
      .type = choice
      .help = this should be obvious.  Either keep the good ones or the bad ones.
  }
  print_correlations = True
    .type = bool
    .help = For each panel group, print the correlation between radial offset\
            and delta_psi, and transverse offset and delta_psi.
}

repredict
  .expert_level = 2 {
  enable = False
    .type = bool
    .help = If True, repredict reflection positions using bandpass and mosaic  \
            approximations
  mode = tophat_mosaicity_and_bandpass *gaussian_mosaicity_and_bandpass
    .type = choice
    .help = Use tophat functions or gaussians to approximate the mosaic        \
            parameters of the crystal and spectral properties of the beam
  refine_mode = *per_experiment all None
    .type = choice
    .help = Refine mosaicity for each experiment individually (per_experiment) \
            or refine a single mosaicity of all experiments (all).  If None,   \
            do not refine mosaicity.
  refine_bandpass = True
    .type = bool
    .help = If True and if also refining mosaicity, then refine the bandpass.
  initial_mosaic_parameters = None
    .type = floats(size=2)
    .help = If set, provide two numbers: domain size (angstroms) and half      \
            mosaic angle (degrees). If unset, use the crystal model's          \
            mosaic parameters.
}

plots {
  all_plots = False
    .type = bool
    .help = Override the rest of the flags in the plots phil scope to show    \
            all plots available
  reflection_energies = False
    .type = bool
    .help = Plot the energy each reflection would need to be exactly on the   \
            Ewald sphere
  pos_vs_neg_delta_psi = False
    .type = bool
    .help = For each image, determine how many reflections have a delta psi   \
            > 0 and how many reflections have a delta psi < 0. It is expected \
            that for most images these numbers will be about the same. Then   \
            create a 2D histogram, where each pixel is how many images have   \
            Y reflections with delta psi < 0  vs. X reflections with delta    \
            psi > 0. It is expected most of the data lies on an x=y line.
  deltaXY_histogram = False
    .type = bool
    .help = Histogram of the distance between the observed and predicted      \
            reflections. Mean and mode are plotted, but it is expected that   \
            this forms a Rayleigh distribution, so the Rayleigh mean is also  \
            shown. Further, RMSD and Rayleigh RMSD obs - pred are also shown. \
            This plot is also made for each panel.
  per_image_RMSDs_histogram = True
    .type = bool
    .help = For each image, compute the RMSD of the observered - predicted    \
            reflection locations. Then histogram these RMSDs over all the     \
            images. This plot is also made for each panel.
  per_image_RMSDs_boxplot = False
    .type = bool
    .help = Box plot of per-image RMSDs
  positional_displacements = True
    .type = bool
    .help = Each reflection is plotted as observed on the detector. The color \
            of the reflection is the difference between the reflection's      \
            observed and predicted locations.
  include_radial_and_transverse = True
    .type = bool
    .help = Also plot radial and transverse displacements: the displacement   \
            vectors between observed and predicted spot locations are split   \
            into radial and transverse components, where radial is the        \
            component of the vector along the line from the refleciton to the \
            beam center, and the transvrse component is the component of the  \
            diplacement vector along the line orthogonal to the radial        \
            component.
  deltaXY_by_deltapsi = True
    .type = bool
    .help = For each reflection, compute the displacement vector deltaXY      \
            between the observed and predicted spot locations. Using the      \
            center of the panel as an origin, plot the reflection displaced   \
            from the center of its panel along its deltaXY vector. The        \
            reflections are colored by their delta psi values.
  deltaXY_by_reflection_energy = False
    .type = bool
    .help = As deltaXY_by_deltapsi, but reflections are colored mean pixel   \
            energy. Every pixel in a reflection could be exactly on the Ewald \
            sphere given a specific wavelength. The mean pixel energy for a   \
            reflection is the mean energy of the pixels in the reflection if  \
            each pixel was individually on the Ewald sphere.
  repredict_from_reflection_energies = False
    .type = bool
    .help = As deltaXY_by_deltapsi, but repredict each reflection using the  \
            wavelength needed to place that reflection on the Ewald sphere.
  deltaXY_by_deltaXY = False
    .type = bool
    .help = As deltaXY_by_deltapsi, but color reflections by the magnitude   \
            of the deltaXY vector.
  manual_cdf = False
    .type = bool
  radial_vs_deltaPsi_vs_deltaXY = False
    .type = bool
    .help = For each panel, each reflection is plotted given an origin in the \
            center of the panel. Y: radial displacement (obs-pred along the   \
            radial direction), X: delta psi. Colored by the magnitude of      \
            deltaXY
  radial_difference_histograms = False
    .type = bool
    .help = For each panel, the radial displacements of each reflection are   \
            histogrammed.  The center is marked with a red line. Asymmetry    \
            indicates a panel not properly aligned in the radial direciton.
  intensity_vs_radials_2dhist = False
    .type = bool
    .help = 2D histogram of reflection intensities vs radial displacements.   \
            Each pixel is the number of reflections with a give intensity and \
            a given radial displacement.
  delta2theta_vs_deltapsi_2dhist = False
    .type = bool
    .help = For each reflection, compute the difference in measured vs pred-   \
            icted two theta and the delta psi. 2D histogram is plotted where  \
            each pixel is the number of reflections with a given delta two    \
            theta and a given delta psi.  10 plots are shown, one for each of \
            10 resolution bins.
  delta2theta_vs_2theta_2dhist = False
    .type = bool
    .help = For each reflection, compute the two theta and the difference in  \
            measured vs predicted two theta. 2D histogram is plotted where    \
            each pixel is the number of reflections with a given delta two    \
            theta and a given two theta.
  deltaPsi_vs_2theta_2dhist = False
    .type = bool
    .help = For each reflection, compute the delta psi and two theta angles.  \
            2D histogram is plotted where each pixel is the number of         \
            reflections with a given delta psi and two theta. Result is       \
            similar to a trumpet plot but for the whole dataset.
  grouped_stats = False
    .type = bool
    .help = 5 plots are shown with different stats. For each panel group, the \
            stat is shown and the panel group is colored by that stat. Stats  \
            are: number of reflections, overall, radial and transverse RMSDs, \
            and the CC between delta 2 theta and delta psi among reflections  \
            in that panel group.
  unit_cell_histograms = True
    .type = bool
    .help = Unit cell histograms are shown for each of the a, b, and c axes.
  stats_by_2theta = False
    .type = bool
    .help = Reflections are binned by 2 theta, then various stats are computed\
            per bin. RMSD, R RMSD, T RMSD, RMSD delta 2theta: observed -      \
            predicted RMSDs, including overall, radial and transverse as well \
            as RMSD of delta 2theta. R/T RMSD: ratio of radial and transverse \
            RMSDs.
  stats_by_panelgroup = False
    .type = bool
    .help = As stats_by_2theta, but reflections are binned by panelgroup.
  trumpet_plot = False
    .type = bool
    .help = Show trumpet plot from Sauter 2014 for the first experiment
  ewald_offset_plot = False
    .type = bool
    .help = Show trumpet plot equivalent but using the Ewald offset for Y \
            (I.E. the reflection distance from the Ewald sphere)
  include_offset_dots = False
    .type = bool
    .help = Plot mean offset for spot predictions from observations in \
            deltaXY_by_deltapsi plot
  include_scale_bar_in_pixels = 0
    .type = float (value_min=0)
    .help = For each panel in the deltaXY_by_deltapsi plot, show xy scale \
            bars. Length of the scale bar is specified here in pixels.
}

save_pdf = False
  .type = bool
  .help = Whether to show the plots or save as a multi-page pdf
save_png = False
  .type = bool
  .help = Whether to show the plots or save as a series of pngs
include scope xfel.command_line.cspad_detector_congruence.phil_scope
''', process_includes=True)

def setup_stats(experiments, reflections, two_theta_only = False):
  # Compute a set of radial and transverse displacements for each reflection
  print("Setting up stats...")
  tmp = flex.reflection_table()
  # Need to construct a variety of vectors
  for expt_id, expt in enumerate(experiments):
    expt_refls = reflections.select(reflections['id'] == expt_id)
    if len(expt_refls) == 0: continue
    for panel_id, panel in enumerate(expt.detector):
      refls = expt_refls.select(expt_refls['panel'] == panel_id)
      if len(refls) == 0: continue
      # Compute the beam center in lab space (a vector pointing from the origin to where the beam would intersect
      # the panel, if it did intersect the panel)
      beam = expt.beam
      s0 = beam.get_s0()
      beam_centre = panel.get_beam_centre_lab(s0)
      bcl = flex.vec3_double(len(refls), beam_centre)
      cal_x, cal_y, _ = refls['xyzcal.px'].parts()
      ttc = flex.double([panel.get_two_theta_at_pixel(s0, (cal_x[i], cal_y[i])) for i in range(len(refls))])
      if 'xyzobs.px.value' in refls:
        obs_x, obs_y, _ = refls['xyzobs.px.value'].parts()
        tto = flex.double([panel.get_two_theta_at_pixel(s0, (obs_x[i], obs_y[i])) for i in range(len(refls))])
      refls['beam_centre_lab'] = bcl
      refls['two_theta_cal'] = ttc * (180/math.pi) #+ (0.5*panel_refls['delpsical.rad']*panel_refls['two_theta_obs'])
      if 'xyzobs.px.value' in refls:
        refls['two_theta_obs'] = tto * (180/math.pi)
      if not two_theta_only:
        # Compute obs in lab space
        x, y, _ = refls['xyzobs.mm.value'].parts()
        c = flex.vec2_double(x, y)
        refls['obs_lab_coords'] = panel.get_lab_coord(c)
        # Compute deltaXY in panel space. This vector is relative to the panel origin
        x, y, _ = (refls['xyzcal.mm'] - refls['xyzobs.mm.value']).parts()
        # Convert deltaXY to lab space, subtracting off of the panel origin
        refls['delta_lab_coords'] = panel.get_lab_coord(flex.vec2_double(x,y)) - panel.get_origin()
      tmp.extend(refls)
  reflections = tmp
  return reflections

def get_unweighted_rmsd(reflections, verbose=True):
  n = len(reflections)
  if n == 0:
    return 0
  #weights = 1/reflections['intensity.sum.variance']
  reflections = reflections.select(reflections['xyzobs.mm.variance'].norms() > 0)
  weights = 1/reflections['xyzobs.mm.variance'].norms()

  un_rmsd = math.sqrt( flex.sum(reflections['difference_vector_norms']**2)/n)
  w_rmsd = math.sqrt( flex.sum( weights*(reflections['difference_vector_norms']**2) )/flex.sum(weights))

  if verbose:
    print("%20s%7.3f"%("Unweighted RMSD (μm)", un_rmsd*1000))
    print("%20s%7.3f"%("Weighted RMSD (μm)", w_rmsd*1000))

  return un_rmsd


def reflection_wavelength_from_pixels(experiments, reflections):
  if 'shoebox' not in reflections:
    return reflections

  print("Computing per-reflection wavelengths from shoeboxes")

  from dials.algorithms.shoebox import MaskCode
  valid_code = MaskCode.Valid | MaskCode.Foreground

  table = flex.reflection_table()
  wavelengths = flex.double()
  for expt_id, expt in enumerate(experiments):
    refls = reflections.select(reflections['id'] == expt_id)
    table.extend(refls)

    beam = expt.beam
    unit_cell = expt.crystal.get_unit_cell()
    d = unit_cell.d(refls['miller_index'])

    for i in range(len(refls)):
      sb = refls['shoebox'][i]
      # find the coordinates with signal
      mask = flex.bool([(m & valid_code) != 0 for m in sb.mask])
      coords = sb.coords().select(mask)
      panel = expt.detector[refls['panel'][i]]

      # compute two theta angle for each pixel
      s1 = panel.get_lab_coord(panel.pixel_to_millimeter(flex.vec2_double(coords.parts()[0], coords.parts()[1])))
      two_theta = s1.angle(beam.get_s0())

      # nLambda = 2dST
      wavelengths.append(flex.mean(2*d[i]*flex.sin(two_theta/2)))

  table['reflection_wavelength_from_pixels'] = wavelengths
  return table

def trumpet_plot(experiment, reflections, axis = None):
  half_mosaicity_deg = experiment.crystal.get_half_mosaicity_deg()
  domain_size_ang = experiment.crystal.get_domain_size_ang()
  if not axis:
    fig = plt.figure()
    axis = plt.gca()
  two_thetas = reflections['two_theta_cal']
  delpsi = reflections['delpsical.rad']*180/math.pi
  axis.scatter(two_thetas, delpsi)

  LR = flex.linear_regression(two_thetas, delpsi)
  model_y = LR.slope()*two_thetas + LR.y_intercept()
  axis.plot(two_thetas, model_y, "k-")

  tan_phi_deg = (experiment.crystal.get_unit_cell().d(reflections['miller_index']) / domain_size_ang)*180/math.pi
  tan_outer_deg = tan_phi_deg + (half_mosaicity_deg/2)

  axis.set_title("Mosaicity FW=%4.2f deg, Dsize=%5.0fA on %d spots"%(2*half_mosaicity_deg, domain_size_ang, len(two_thetas)))
  axis.plot(two_thetas, tan_phi_deg, "r.")
  axis.plot(two_thetas, -tan_phi_deg, "r.")
  axis.plot(two_thetas, tan_outer_deg, "g.")
  axis.plot(two_thetas, -tan_outer_deg, "g.")

from xfel.command_line.cspad_detector_congruence import iterate_detector_at_level, iterate_panels, id_from_name, get_center, detector_plot_dict
from xfel.command_line.cspad_detector_congruence import Script as DCScript
class Script(DCScript):
  ''' Class to parse the command line options. '''

  def __init__(self):
    ''' Set the expected options. '''
    from dials.util.options import ArgumentParser
    import libtbx.load_env

    # Create the option parser
    usage = "usage: %s [options] /path/to/refined/json/file" % libtbx.env.dispatcher_name
    self.parser = ArgumentParser(
      usage=usage,
      sort_options=True,
      phil=phil_scope,
      read_experiments=True,
      read_reflections=True,
      check_format=False,
      epilog=help_message)

  def run(self):
    ''' Parse the options. '''
    from dials.util.options import flatten_experiments, flatten_reflections
    # Parse the command line arguments
    params, options = self.parser.parse_args(show_diff_phil=True)
    self.params = params
    if params.plots.all_plots:
      for attr in dir(params.plots):
        if attr.startswith('__'): continue
        setattr(params.plots, attr, True)

    experiments = flatten_experiments(params.input.experiments)

    # Find all detector objects
    detectors = experiments.detectors()

    # Verify inputs
    if len(params.input.reflections) == len(detectors) and len(detectors) > 1:
      # case for passing in multiple images on the command line
      assert len(params.input.reflections) == len(detectors)
      reflections = flex.reflection_table()
      for expt_id in range(len(detectors)):
        subset = params.input.reflections[expt_id].data
        subset['id'] = flex.int(len(subset), expt_id)
        reflections.extend(subset)
    else:
      # case for passing in combined experiments and reflections
      reflections = flatten_reflections(params.input.reflections)[0]

    ResidualsPlotter(params, experiments, reflections).plot_all()

class ResidualsPlotter(object):
  def __init__(self, params, experiments, reflections):
    self.params = params
    self.experiments = experiments
    self.reflections = reflections

  def get_normalized_colors(self, data, vmin=None, vmax=None):
    if vmax is None:
      vmax = self.params.residuals.plot_max
    if vmax is None:
      vmax = flex.max(data)
    if vmin is None:
      vmin = min(flex.min(data), 0)
    if len(data) > 1:
      assert vmax > vmin, "vmax: %f, vmin: %f"%(vmax, vmin)

    # initialize the color map
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.cm.get_cmap(self.params.colormap)
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    color_vals = np.linspace(vmin, vmax, 11)
    sm.set_array(color_vals) # needed for colorbar

    return norm, cmap, color_vals, sm

  def plot_deltas(self, reflections, panel = None, ax = None, bounds = None):
    assert panel is not None and ax is not None and bounds is not None

    data = (reflections['xyzcal.mm']-reflections['xyzobs.mm.value']).norms()
    norm, cmap, color_vals, sm = self.get_normalized_colors(data)
    deltas = (reflections['xyzcal.mm']-reflections['xyzobs.mm.value'])*self.delta_scalar

    x, y = panel.get_image_size_mm()
    offset = col((x, y, 0))/2
    deltas += offset
    mm_panel_coords = flex.vec2_double(deltas.parts()[0], deltas.parts()[1])

    lab_coords = panel.get_lab_coord(mm_panel_coords)

    ax.scatter(lab_coords.parts()[0], lab_coords.parts()[1], c = data, norm=norm, cmap = cmap, linewidths=0, s=self.params.dot_size)

    return sm, color_vals

  def plot_obs_colored_by_radial_deltas(self, reflections, panel = None, ax = None, bounds = None):
    return self.plot_obs_colored_by_data(flex.abs(reflections['radial_displacements']), reflections, panel, ax, bounds)

  def plot_obs_colored_by_transverse_deltas(self, reflections, panel = None, ax = None, bounds = None):
    return self.plot_obs_colored_by_data(flex.abs(reflections['transverse_displacements']), reflections, panel, ax, bounds)

  def plot_obs_colored_by_deltas(self, reflections, panel = None, ax = None, bounds = None):
    data = (reflections['xyzcal.mm']-reflections['xyzobs.mm.value']).norms()
    return self.plot_obs_colored_by_data(data, reflections, panel, ax, bounds)

  def plot_obs_colored_by_data(self, data, reflections, panel = None, ax = None, bounds = None):
    assert panel is not None and ax is not None and bounds is not None
    norm, cmap, color_vals, sm = self.get_normalized_colors(data)
    mm_panel_coords = flex.vec2_double(reflections['xyzobs.mm.value'].parts()[0], reflections['xyzobs.mm.value'].parts()[1])
    lab_coords = panel.get_lab_coord(mm_panel_coords)

    ax.scatter(lab_coords.parts()[0], lab_coords.parts()[1], c = data, norm=norm, cmap = cmap, linewidths=0, s=self.params.dot_size)

    return sm, color_vals

  def plot_obs_colored_by_deltapsi(self, reflections, panel = None, ax = None, bounds = None):
    if 'delpsical.rad' not in reflections:
      return
    assert panel is not None and ax is not None and bounds is not None
    data = reflections['delpsical.rad'] * (180/math.pi)
    norm, cmap, color_vals, sm = self.get_normalized_colors(data, vmin=-0.1, vmax=0.1)
    deltas = (reflections['xyzcal.mm']-reflections['xyzobs.mm.value'])*self.delta_scalar

    x, y = panel.get_image_size_mm()
    offset = col((x, y, 0))/2
    deltas += offset
    mm_panel_coords = flex.vec2_double(deltas.parts()[0], deltas.parts()[1])

    lab_coords = panel.get_lab_coord(mm_panel_coords)

    lab_coords_x, lab_coords_y, _ = lab_coords.parts()
    if self.params.residuals.mcd_filter.enable and len(reflections)>5:
      from xfel.metrology.panel_fitting import Panel_MCD_Filter
      MCD = Panel_MCD_Filter(lab_coords_x, lab_coords_y, data, i_panel = reflections["panel"][0],
                      delta_scalar = self.delta_scalar, params = self.params)
      sX,sY,sPsi = MCD.scatter_coords()
      MCD.plot_contours(ax,show=False) # run this to pre-compute the center position
      ax.scatter(sX, sY, c = sPsi, norm=norm, cmap = cmap, linewidths=0, s=self.params.dot_size)
    else:
      ax.scatter(lab_coords_x, lab_coords_y, c = data, norm=norm, cmap = cmap, linewidths=0, s=self.params.dot_size)

    if self.params.plots.include_offset_dots:
      panel_center = panel.get_lab_coord(offset[0:2])
      ax.scatter(panel_center[0], panel_center[1], c='k', s=self.params.dot_size)
      if self.params.plots.include_scale_bar_in_pixels > 0:
        pxlsz0,pxlsz1 = panel.get_pixel_size()
        ax.plot([panel_center[0],panel_center[0]],
                [panel_center[1],panel_center[1]+self.params.plots.include_scale_bar_in_pixels*pxlsz1*self.delta_scalar], 'k-')
        ax.plot([panel_center[0],panel_center[0]+self.params.plots.include_scale_bar_in_pixels*pxlsz0*self.delta_scalar],
                [panel_center[1],panel_center[1]], 'k-')
      ax.scatter(flex.mean(lab_coords_x), flex.mean(lab_coords_y), c='b', s=self.params.dot_size)
      if self.params.residuals.mcd_filter.enable:
        ax.scatter(MCD.robust_model_XY.location_[0],
                   MCD.robust_model_XY.location_[1], c='r', s=self.params.dot_size)
      print(panel.get_name(), (flex.mean(lab_coords_x) - panel_center[0])/self.delta_scalar, (flex.mean(lab_coords_y) - panel_center[0])/self.delta_scalar)

    return sm, color_vals

  def plot_obs_colored_by_deltapsi_pxlambda(self, reflections, panel = None, ax = None, bounds = None):
    if 'delpsical.rad' not in reflections \
      or 'delpsical.rad.pxlambda' not in reflections \
      or 'xyzcal.mm.pxlambda' not in reflections:
      return
    assert panel is not None and ax is not None and bounds is not None
    data = reflections['delpsical.rad.pxlambda'] * (180/math.pi)
    norm, cmap, color_vals, sm = self.get_normalized_colors(data, vmin=-0.1, vmax=0.1)
    deltas = (reflections['xyzcal.mm.pxlambda']-reflections['xyzobs.mm.value'])*self.delta_scalar

    x, y = panel.get_image_size_mm()
    offset = col((x, y, 0))/2
    deltas += offset
    mm_panel_coords = flex.vec2_double(deltas.parts()[0], deltas.parts()[1])

    lab_coords = panel.get_lab_coord(mm_panel_coords)

    ax.scatter(lab_coords.parts()[0], lab_coords.parts()[1], c = data, norm=norm, cmap = cmap, linewidths=0, s=self.params.dot_size)

    return sm, color_vals

  def plot_obs_colored_by_mean_pixel_wavelength(self, reflections, panel = None, ax = None, bounds = None):
    if 'reflection_wavelength_from_pixels' not in reflections:
      return
    assert panel is not None and ax is not None and bounds is not None
    data = 12398.4/reflections['reflection_wavelength_from_pixels']
    norm, cmap, color_vals, sm = self.get_normalized_colors(data, vmin=self.min_energy, vmax=self.max_energy)
    deltas = (reflections['xyzcal.mm']-reflections['xyzobs.mm.value'])*self.delta_scalar

    x, y = panel.get_image_size_mm()
    offset = col((x, y, 0))/2
    deltas += offset
    mm_panel_coords = flex.vec2_double(deltas.parts()[0], deltas.parts()[1])

    lab_coords = panel.get_lab_coord(mm_panel_coords)

    ax.scatter(lab_coords.parts()[0], lab_coords.parts()[1], c = data, norm=norm, cmap = cmap, linewidths=0, s=self.params.dot_size)

    return sm, color_vals

  def plot_radial_displacements_vs_deltapsi(self, reflections, panel = None, ax = None, bounds = None):
    if 'delpsical.rad' not in reflections:
      return
    assert panel is not None and ax is not None and bounds is not None
    data = reflections['difference_vector_norms']
    norm, cmap, color_vals, sm = self.get_normalized_colors(data)

    a = reflections['delpsical.rad']*180/math.pi
    b = reflections['radial_displacements']

    fake_coords = flex.vec2_double(a, b) * self.delta_scalar

    x, y = panel.get_image_size_mm()
    offset = col((x, y))/2

    lab_coords = fake_coords + panel.get_lab_coord(offset)[0:2]

    ax.scatter(lab_coords.parts()[0], lab_coords.parts()[1], c = data, norm=norm, cmap = cmap, linewidths=0, s=self.params.dot_size)

    return sm, color_vals

  def plot_unitcells(self, experiments):
    if len(experiments) == 1:
      return
    all_a = flex.double()
    all_b = flex.double()
    all_c = flex.double()
    for crystal in experiments.crystals():
      a, b, c = crystal.get_unit_cell().parameters()[0:3]
      all_a.append(a); all_b.append(b); all_c.append(c)

    fig, axes = plt.subplots(nrows=3, ncols=1)
    for ax, axis, data in zip(axes, ['A', 'B', 'C'], [all_a, all_b, all_c]):
      stats = flex.mean_and_variance(data)
      cutoff = 4*stats.unweighted_sample_standard_deviation()
      if cutoff < 0.5:
        cutoff = 0.5
      limits = stats.mean()-cutoff, stats.mean()+cutoff
      sel = (data >= limits[0]) & (data <= limits[1])
      subset = data.select(sel)
      h = flex.histogram(subset,n_slots=50)
      ax.plot(h.slot_centers().as_numpy_array(),h.slots().as_numpy_array(),'-')
      ax.set_title("%s axis histogram (showing %d of %d xtals). Mean: %7.2f Stddev: %7.2f"%(
        axis, len(subset), len(data), stats.mean(),
        stats.unweighted_sample_standard_deviation()))
      ax.set_ylabel("N lattices")
      ax.set_xlabel(r"$\AA$")
      ax.set_xlim(limits)
    plt.tight_layout()

  def plot_histograms(self, reflections, panel = None, ax = None, bounds = None):
    data = reflections['difference_vector_norms']
    colors = ['b-', 'g-', 'g--', 'r-', 'b-', 'b--']
    n_slots = 20
    if self.params.residuals.histogram_max is None:
      h = flex.histogram(data, n_slots=n_slots)
    else:
      h = flex.histogram(data.select(data <= self.params.residuals.histogram_max), n_slots=n_slots)

    n = len(reflections)
    rmsd_obs = math.sqrt((reflections['xyzcal.mm']-reflections['xyzobs.mm.value']).sum_sq()/n)
    sigma = mode = h.slot_centers()[list(h.slots()).index(flex.max(h.slots()))]
    mean_obs = flex.mean(data)
    median = flex.median(data)
    mean_rayleigh = math.sqrt(math.pi/2)*sigma
    rmsd_rayleigh = math.sqrt(2)*sigma

    data = flex.vec2_double([(i,j) for i, j in zip(h.slot_centers(), h.slots())])
    n = len(data)
    for i in [mean_obs, mean_rayleigh, mode, rmsd_obs, rmsd_rayleigh]:
      data.extend(flex.vec2_double([(i, 0), (i, flex.max(h.slots()))]))
    data = self.get_bounded_data(data, bounds)
    tmp = [data[:n]]
    for i in range(len(colors)):
      tmp.append(data[n+(i*2):n+((i+1)*2)])
    data = tmp

    for d, c in zip(data, colors):
      ax.plot(d.parts()[0], d.parts()[1], c)

    if ax.get_legend() is None:
      ax.legend([r"$\Delta$XY", "MeanObs", "MeanRayl", "Mode", "RMSDObs", "RMSDRayl"])

  def plot_cdf_manually(self, reflections, panel = None, ax = None, bounds = None):
    colors = ['blue', 'green']
    r = (reflections['xyzcal.mm']-reflections['xyzobs.mm.value']).norms()
    h = flex.histogram(r)
    sigma = h.slot_centers()[list(h.slots()).index(flex.max(h.slots()))] # mode

    x_extent = max(r)
    y_extent = len(r)
    xobs = [i/x_extent for i in sorted(r)]
    yobs = [i/y_extent for i in range(y_extent)]
    obs = [(x, y) for x, y in zip(xobs, yobs)]

    ncalc = 100
    xcalc = [i/ncalc for i in range(ncalc)]
    ycalc = [1-math.exp((-i**2)/(2*(sigma**2))) for i in xcalc]
    calc = [(x, y) for x, y in zip(xcalc, ycalc)]

    data = [flex.vec2_double(obs),
            flex.vec2_double(calc)]
    if bounds is None:
      ax.set_xlim((-1,1))
      ax.set_ylim((-1,1))
      ax.set_title("%s Outlier SP Manually"%self.params.tag)
    if bounds is not None:
      data = [self.get_bounded_data(d, bounds) for d in data]

    if ax is None:
      fig = plt.figure()
      ax = fig.add_subplot(111)

    for subset,c in zip(data, colors):
        ax.plot(subset.parts()[0], subset.parts()[1], '-', c=c)

  def plot_radial_difference_histograms(self, reflections, panel = None, ax = None, bounds = None):
    r = reflections['radial_displacements']*1000
    h = flex.histogram(r, n_slots=10, data_min=-300, data_max=300)

    x = h.slot_centers()
    y = h.slots().as_double()
    x.append(0); y.append(0)
    x.append(0); y.append(flex.max(y))

    if bounds is None:
      ax.set_title("%s Radial differences"%self.params.tag)
    else:
      d = flex.vec2_double(x, y)
      data = self.get_bounded_data(d, bounds)
      x, y = data.parts()

    linex = [x[-2],x[-1]]; liney = [y[-2],y[-1]]
    x = x[:-2]; y = y[:-2]

    if ax is None:
      fig = plt.figure()
      ax = fig.add_subplot(111)
    ax.plot(x.as_numpy_array(), y.as_numpy_array(), '-', c='blue')
    ax.plot(linex, liney, '-', c='red')

  def plot_difference_vector_norms_histograms(self, reflections, panel = None, ax = None, bounds = None):
    r = reflections['difference_vector_norms']*1000
    h = flex.histogram(r, n_slots=10, data_min=0, data_max=100)

    x = h.slot_centers()
    y = h.slots().as_double()

    if bounds is None:
      ax.set_title("%s Residual norms histogram"%self.params.tag)
    else:
      d = flex.vec2_double(x, y)
      data = self.get_bounded_data(d, bounds)
      x, y = data.parts()

    if ax is None:
      fig = plt.figure()
      ax = fig.add_subplot(111)
    ax.plot(x.as_numpy_array(), y.as_numpy_array(), '-', c='blue')

  def get_bounded_data(self, data, bounds):
    assert len(bounds) == 4
    x = [b[0] for b in bounds]
    y = [b[1] for b in bounds]
    left = sorted(x)[1]
    right = sorted(x)[2]
    top = sorted(y)[2]
    bottom = sorted(y)[1]
    origin = col((left, bottom))
    scale_x = right-left
    scale_y = top-bottom
    scale = min(scale_x, scale_y)

    data_max_x = flex.max(data.parts()[0])
    data_min_x = flex.min(data.parts()[0])
    data_max_y = flex.max(data.parts()[1])
    data_min_y = flex.min(data.parts()[1])
    data_scale_x = data_max_x - data_min_x
    data_scale_y = data_max_y - data_min_y

    if data_scale_x == 0 or data_scale_y == 0:
      print("WARNING bad scale")
      return data

    xscale = scale/abs(data_scale_x)
    yscale = scale/abs(data_scale_y)

    return flex.vec2_double(((data.parts()[0] * xscale) - (data_min_x * xscale)),
                            ((data.parts()[1] * yscale) - (data_min_y * yscale))) + origin

  def plot_all(self):
    params = self.params
    experiments = self.experiments
    reflections = self.reflections

    detector = experiments.detectors()[0]

    if params.repredict_input_reflections:
      from dials.algorithms.refinement.prediction.managed_predictors import ExperimentsPredictorFactory
      ref_predictor = ExperimentsPredictorFactory.from_experiments(experiments, force_stills=experiments.all_stills())
      reflections = ref_predictor(reflections)

    if params.verbose: print("N reflections total:", len(reflections))
    if params.residuals.exclude_outliers_from_refinement:
      reflections = reflections.select(reflections.get_flags(reflections.flags.used_in_refinement))
      if params.verbose: print("N reflections used in refinement:", len(reflections))
      if params.verbose: print("Reporting only on those reflections used in refinement")

    if params.residuals.recompute_outliers:
      print("Performing outlier rejection on %d reflections"%len(reflections))
      from dials.algorithms.refinement.outlier_detection.sauter_poon import SauterPoon
      outlier = SauterPoon(px_sz=experiments[0].detector[0].get_pixel_size(), separate_panels=False)
      rejection_occured = outlier(reflections)
      if rejection_occured:
        reflections = reflections.select(~reflections.get_flags(reflections.flags.centroid_outlier))
        if params.verbose: print("N reflections after outlier rejection:", len(reflections))
      else:
        if params.verbose: print("No rejections found")

    if self.params.residuals.i_sigi_cutoff is not None:
      sel = (reflections['intensity.sum.value']/flex.sqrt(reflections['intensity.sum.variance'])) >= self.params.residuals.i_sigi_cutoff
      reflections = reflections.select(sel)
      if params.verbose: print("After filtering by I/sigi cutoff of %f, there are %d reflections left"%(self.params.residuals.i_sigi_cutoff,len(reflections)))

    if 'shoebox' in reflections and (params.repredict.enable or (params.show_plots and params.plots.reflection_energies)):
      reflections = reflection_wavelength_from_pixels(experiments, reflections)
      stats = flex.mean_and_variance(12398.4/reflections['reflection_wavelength_from_pixels'])
      if params.verbose: print("Mean energy: %.1f +/- %.1f"%(stats.mean(), stats.unweighted_sample_standard_deviation()))
      self.min_energy = stats.mean() - stats.unweighted_sample_standard_deviation()
      self.max_energy = stats.mean() + stats.unweighted_sample_standard_deviation()

      try:
        from dials_scratch.asb.predictions_from_reflection_wavelengths import predictions_from_per_reflection_energies, tophat_vector_wavelengths, refine_wavelengths, wavelengths_from_gaussians
      except ImportError:
        if params.repredict.enable:
          from libtbx.utils import Sorry
          raise Sorry("dials_scratch not configured so cannot do reprediction")
      else:
        reflections = predictions_from_per_reflection_energies(experiments, reflections, 'reflection_wavelength_from_pixels', 'pxlambda')

        if params.show_plots and params.plots.reflection_energies:
          fig = plt.figure()
          stats = flex.mean_and_variance(12398.4/reflections['reflection_wavelength_from_pixels'])
          plt.title("Energies derived from indexed pixels, mean: %.1f +/- %.1f"%(stats.mean(), stats.unweighted_sample_standard_deviation()))
          plt.hist(12398.4/reflections['reflection_wavelength_from_pixels'], bins=100)
          plt.xlabel("Energy (eV)")
          plt.ylabel("Count")

      if params.repredict.enable:
        init_mp = params.repredict.initial_mosaic_parameters

        if params.repredict.mode == 'tophat_mosaicity_and_bandpass':
          tag = 'reflection_wavelength_from_mosaicity_and_bandpass'
          dest = 'mosbandp'
          func = tophat_vector_wavelengths
          gaussians = False
        elif params.repredict.mode == 'gaussian_mosaicity_and_bandpass':
          tag = 'reflection_wavelength_from_gaussian_mosaicity_and_bandpass'
          dest = 'gmosbandp'
          func = wavelengths_from_gaussians
          gaussians = True

        if params.repredict.refine_mode == 'per_experiment':
          refined_reflections = flex.reflection_table()
          for expt_id in range(len(experiments)):
            if params.verbose: print("*"*80, "EXPERIMENT", expt_id)
            refls = reflections.select(reflections['id']==expt_id)
            refls['id'] = flex.int(len(refls), 0)
            refls = refine_wavelengths(experiments[expt_id:expt_id+1], refls, init_mp, tag, dest,
              refine_bandpass=params.repredict.refine_bandpass, gaussians=gaussians)
            refls['id'] = flex.int(len(refls), expt_id)
            refined_reflections.extend(refls)
          reflections = refined_reflections
        elif params.repredict.refine_mode == 'all':
          reflections = refine_wavelengths(experiments, reflections, init_mp, tag, dest,
            refine_bandpass=params.repredict.refine_bandpass, gaussians=gaussians)
        elif params.repredict.refine_mode is None or params.repredict.refine_mode == 'None':
          reflections = func(experiments, reflections, init_mp)
        reflections = predictions_from_per_reflection_energies(experiments, reflections, tag, dest)
        stats = flex.mean_and_variance(12398.4/reflections[tag])
        if params.verbose: print("Mean energy: %.1f +/- %.1f"%(stats.mean(), stats.unweighted_sample_standard_deviation()))
        reflections['delpsical.rad'] = reflections['delpsical.rad.%s'%dest]
        reflections['xyzcal.mm'] = reflections['xyzcal.mm.%s'%dest]
        reflections['xyzcal.px'] = reflections['xyzcal.px.%s'%dest]

    if 'xyzobs.mm.value' not in reflections:
      reflections.centroid_px_to_mm(experiments)
    reflections['difference_vector_norms'] = (reflections['xyzcal.mm']-reflections['xyzobs.mm.value']).norms()

    n = len(reflections)
    rmsd = get_unweighted_rmsd(reflections, params.verbose)
    print("%20s%7.3f"%("Dataset RMSD (μm)", rmsd * 1000))

    if params.tag is None:
      tag = ''
    else:
      tag = '%s '%params.tag

    if 'delpsical.rad' in reflections and params.show_plots and params.plots.pos_vs_neg_delta_psi:
      # set up delta-psi ratio heatmap
      fig = plt.figure()
      p = flex.int() # positive
      n = flex.int() # negative
      for i in set(reflections['id']):
        exprefls = reflections.select(reflections['id']==i)
        p.append(len(exprefls.select(exprefls['delpsical.rad']>0)))
        n.append(len(exprefls.select(exprefls['delpsical.rad']<0)))
      plt.hist2d(p, n, bins=30)
      cb = plt.colorbar()
      cb.set_label("N images")
      plt.title(r"%s2D histogram of pos vs. neg $\Delta\Psi$ per image"%tag)
      plt.xlabel(r"N reflections with $\Delta\Psi$ > 0")
      plt.ylabel(r"N reflections with $\Delta\Psi$ < 0")

    self.delta_scalar = 50

    # Iterate through the detectors, computing detector statistics at the per-panel level (IE one statistic per panel)
    # Per panel dictionaries
    rmsds = {}
    refl_counts = {}
    transverse_rmsds = {}
    radial_rmsds = {}
    ttdpcorr = {}
    pg_bc_dists = {}
    mean_delta_two_theta = {}
    # per panelgroup flex arrays
    pg_rmsds = flex.double()
    pg_r_rmsds = flex.double()
    pg_t_rmsds = flex.double()
    pg_refls_count = flex.int()
    pg_refls_count_d = {}
    table_header = ["PG id", "RMSD","Radial", "Transverse", "N_refls"]
    table_header2 = ["","(μm)","RMSD(μm)","RMSD(μm)",""]
    if params.residuals.print_correlations:
      table_header += ["Correl", "Correl"]
      table_header2 += ["ΔR,ΔΨ","ΔT,ΔΨ"]
    table_data = []
    table_data.append(table_header)
    table_data.append(table_header2)

    reflections = setup_stats(experiments, reflections)

    # The radial vector points from the center of the reflection to the beam center
    radial_vectors = (reflections['obs_lab_coords'] - reflections['beam_centre_lab']).each_normalize()
    # The transverse vector is orthogonal to the radial vector and the beam vector
    transverse_vectors = radial_vectors.cross(reflections['beam_centre_lab']).each_normalize()
    # Compute the raidal and transverse components of each deltaXY
    reflections['radial_displacements']     = reflections['delta_lab_coords'].dot(radial_vectors)
    reflections['transverse_displacements'] = reflections['delta_lab_coords'].dot(transverse_vectors)

    # Iterate through the detector at the specified hierarchy level
    if hasattr(detector, 'hierarchy'):
      iterable = enumerate(iterate_detector_at_level(detector.hierarchy(), 0, params.hierarchy_level))
    else:
      iterable = enumerate(detector)

    if params.residuals.print_correlations:
      from xfel.metrology.panel_fitting import three_feature_fit
      pg_r_all = flex.double(); pg_t_all = flex.double(); pg_delpsi_all = flex.double()

    s0 = experiments[0].beam.get_s0()

    for pg_id, pg in iterable:
      pg_msd_sum = 0
      pg_r_msd_sum = 0
      pg_t_msd_sum = 0
      pg_refls = 0
      if params.residuals.print_correlations:
        pg_r = flex.double(); pg_t = flex.double()
      pg_delpsi = flex.double()
      pg_deltwotheta = flex.double()
      for p in iterate_panels(pg):
        panel_id = id_from_name(detector, p.get_name())
        panel_refls = reflections.select(reflections['panel'] == panel_id)
        n = len(panel_refls)
        pg_refls += n

        delta_x = panel_refls['xyzcal.mm'].parts()[0] - panel_refls['xyzobs.mm.value'].parts()[0]
        delta_y = panel_refls['xyzcal.mm'].parts()[1] - panel_refls['xyzobs.mm.value'].parts()[1]

        tmp = flex.sum((delta_x**2)+(delta_y**2))
        pg_msd_sum += tmp

        r = panel_refls['radial_displacements']
        t = panel_refls['transverse_displacements']
        if params.residuals.print_correlations:
          pg_r.extend(r); pg_t.extend(t)
        pg_r_msd_sum += flex.sum_sq(r)
        pg_t_msd_sum += flex.sum_sq(t)

        pg_delpsi.extend(panel_refls['delpsical.rad']*180/math.pi)
        pg_deltwotheta.extend(panel_refls['two_theta_obs'] - panel_refls['two_theta_cal'])

      bc = col(pg.get_beam_centre_lab(s0))
      ori = get_center(pg)
      pg_bc_dists[pg.get_name()] = (ori-bc).length()
      if len(pg_deltwotheta) > 0:
        mean_delta_two_theta[pg.get_name()] = flex.mean(pg_deltwotheta)
      else:
        mean_delta_two_theta[pg.get_name()] = 0

      if pg_refls == 0:
        pg_rmsd = pg_r_rmsd = pg_t_rmsd = 0
      else:
        pg_rmsd = math.sqrt(pg_msd_sum/pg_refls) * 1000
        pg_r_rmsd = math.sqrt(pg_r_msd_sum/pg_refls) * 1000
        pg_t_rmsd = math.sqrt(pg_t_msd_sum/pg_refls) * 1000
      pg_rmsds.append(pg_rmsd)
      pg_r_rmsds.append(pg_r_rmsd)
      pg_t_rmsds.append(pg_t_rmsd)
      pg_refls_count.append(pg_refls)
      pg_refls_count_d[pg.get_name()] = pg_refls
      table_data.append(["%d"%pg_id, "%.1f"%pg_rmsd, "%.1f"%pg_r_rmsd, "%.1f"%pg_t_rmsd, "%6d"%pg_refls])
      if params.residuals.print_correlations:
        pg_r_all.extend(pg_r); pg_t_all.extend(pg_t); pg_delpsi_all.extend(pg_delpsi)
        if len(pg_r)>2:
          TF = three_feature_fit(delta_radial = pg_r, delta_transverse = pg_t, delta_psi = pg_delpsi, i_panel=pg_id, verbose=False)
          pg_cc_Rpsi = 100.*TF.cross_correl[2]
          pg_cc_Tpsi = 100.*TF.cross_correl[1]
        else:
          pg_cc_Rpsi = 0.; pg_cc_Tpsi = 0.
        table_data[-1].extend(["%3.0f%%"%pg_cc_Rpsi, "%3.0f%%"%pg_cc_Tpsi])

      refl_counts[pg.get_name()] = pg_refls
      if pg_refls == 0:
        rmsds[p.get_name()] = -1
        radial_rmsds[p.get_name()] = -1
        transverse_rmsds[p.get_name()] = -1
        ttdpcorr[pg.get_name()] = -1
      else:
        rmsds[pg.get_name()] = pg_rmsd
        radial_rmsds[pg.get_name()]     = pg_r_rmsd
        transverse_rmsds[pg.get_name()] = pg_t_rmsd

        lc = flex.linear_correlation(pg_delpsi, pg_deltwotheta)
        ttdpcorr[pg.get_name()] = lc.coefficient()


    r1 = ["Weighted PG mean"]
    r2 = ["Weighted PG stddev"]
    if len(pg_rmsds) > 1:
      stats = flex.mean_and_variance(pg_rmsds, pg_refls_count.as_double())
      r1.append("%.1f"%stats.mean())
      r2.append("%.1f"%stats.gsl_stats_wsd())
      stats = flex.mean_and_variance(pg_r_rmsds, pg_refls_count.as_double())
      r1.append("%.1f"%stats.mean())
      r2.append("%.1f"%stats.gsl_stats_wsd())
      stats = flex.mean_and_variance(pg_t_rmsds, pg_refls_count.as_double())
      r1.append("%.1f"%stats.mean())
      r2.append("%.1f"%stats.gsl_stats_wsd())
    else:
      r1.extend([""]*3)
      r2.extend([""]*3)
    r1.append("")
    r2.append("")
    table_data.append(r1)
    table_data.append(r2)
    table_data.append(["PG Mean", "", "", "", "%8.1f"%flex.mean(pg_refls_count.as_double())])

    if params.residuals.print_correlations:
      TFA = three_feature_fit(delta_radial = pg_r_all, delta_transverse = pg_t_all, delta_psi = pg_delpsi_all,
         i_panel=pg_id, verbose=False)
      table_data.append(["Refls Mean", "", "", "", "",
                       "%3.0f%%"%(100.*TFA.cross_correl[2]), "%3.0f%%"%(100.*TFA.cross_correl[1]) ])

    from libtbx import table_utils
    if params.verbose: print("Detector statistics by panel group (PG)")
    if params.verbose: print(table_utils.format(table_data,has_header=2,justify='center',delim=" "))

    self.histogram(reflections, r"%s$\Delta$XY histogram (mm)"%tag, plots = params.show_plots and params.plots.deltaXY_histogram, verbose = params.verbose)

    if params.show_plots:
      if self.params.tag is None:
        t = ""
      else:
        t = "%s "%self.params.tag
      if params.plots.per_image_RMSDs_histogram: self.image_rmsd_histogram(reflections, tag, boxplot = params.plots.per_image_RMSDs_boxplot)

      # Plots! these are plots with callbacks to draw on individual panels
      if params.plots.positional_displacements:
        self.detector_plot_refls(detector, reflections, '%sOverall positional displacements (mm)'%tag,
                                 show=False, plot_callback=self.plot_obs_colored_by_deltas)
        if params.plots.include_radial_and_transverse:
          self.detector_plot_refls(detector, reflections, '%sRadial positional displacements (mm)'%tag,
                                   show=False, plot_callback=self.plot_obs_colored_by_radial_deltas)
          self.detector_plot_refls(detector, reflections, '%sTransverse positional displacements (mm)'%tag,
                                   show=False, plot_callback=self.plot_obs_colored_by_transverse_deltas)

      if params.plots.deltaXY_by_deltapsi:                self.detector_plot_refls(detector, reflections, r'%s$\Delta\Psi$'%tag,
                                                                                   show=False, plot_callback=self.plot_obs_colored_by_deltapsi, colorbar_units=r"$\circ$")
      if params.plots.deltaXY_by_reflection_energy:       self.detector_plot_refls(detector, reflections, '%sMean pixel energy'%tag,
                                                                                   show=False, plot_callback=self.plot_obs_colored_by_mean_pixel_wavelength, colorbar_units="eV")
      if params.plots.repredict_from_reflection_energies: self.detector_plot_refls(detector, reflections, r'%s$\Delta\Psi$ from mean pixel energies'%tag,
                                                                                   show=False, plot_callback=self.plot_obs_colored_by_deltapsi_pxlambda, colorbar_units=r"$\circ$")
      if params.plots.deltaXY_by_deltaXY:                 self.detector_plot_refls(detector, reflections, r'%s$\Delta$XY*%s'%(tag,
                                                                                   self.delta_scalar), show=False, plot_callback=self.plot_deltas)
      if params.plots.manual_cdf:                         self.detector_plot_refls(detector, reflections, '%sSP Manual CDF'%tag,
                                                                                   show=False, plot_callback=self.plot_cdf_manually)
      if params.plots.deltaXY_histogram:                  self.detector_plot_refls(detector, reflections, r'%s$\Delta$XY Histograms'%tag,
                                                                                   show=False, plot_callback=self.plot_histograms)
      if params.plots.radial_vs_deltaPsi_vs_deltaXY:      self.detector_plot_refls(detector, reflections, r'%sRadial displacements vs. $\Delta\Psi$, colored by $\Delta$XY'%tag,
                                                                                   show=False, plot_callback=self.plot_radial_displacements_vs_deltapsi)
      if params.plots.per_image_RMSDs_histogram:          self.detector_plot_refls(detector, reflections, r'%sPer image RMSD histograms'%tag,
                                                                                   show=False, plot_callback=self.plot_difference_vector_norms_histograms)
      if params.plots.radial_difference_histograms:       self.detector_plot_refls(detector, reflections, r'%sRadial differences'%tag,
                                                                                   show=False, plot_callback=self.plot_radial_difference_histograms)

      if params.plots.intensity_vs_radials_2dhist and 'intensity.sum.value' in reflections:
        # Plot intensity vs. radial_displacement
        fig = plt.figure()
        panel_id = 15
        panel_refls = reflections.select(reflections['panel'] == panel_id)
        a = panel_refls['radial_displacements']
        b = panel_refls['intensity.sum.value']
        sel = (a > -0.2) & (a < 0.2) & (b < 50000)
        plt.hist2d(a.select(sel), b.select(sel), bins=100)
        plt.title("%s2D histogram of intensity vs. radial displacement for panel %d"%(tag, panel_id))
        plt.xlabel("Radial displacement (mm)")
        plt.ylabel("Intensity")
        ax = plt.colorbar()
        ax.set_label("Counts")

      if params.plots.delta2theta_vs_deltapsi_2dhist:
        # Plot delta 2theta vs. deltapsi
        n_bins = 10
        bin_size = len(reflections)//n_bins
        bin_low = []
        bin_high = []
        data = flex.sorted(reflections['two_theta_obs'])
        for i in range(n_bins):
          bin_low = data[i*bin_size]
          if (i+1)*bin_size >= len(reflections):
            bin_high = data[-1]
          else:
            bin_high = data[(i+1)*bin_size]
          refls = reflections.select((reflections['two_theta_obs'] >= bin_low) &
                                     (reflections['two_theta_obs'] <= bin_high))
          a = refls['delpsical.rad']*180/math.pi
          b = refls['two_theta_obs'] - refls['two_theta_cal']
          fig = plt.figure()
          sel = (a > -0.2) & (a < 0.2) & (b > -0.05) & (b < 0.05)
          plt.hist2d(a.select(sel), b.select(sel), bins=50, range = [[-0.2, 0.2], [-0.05, 0.05]])
          cb = plt.colorbar()
          cb.set_label("N reflections")
          plt.title(r'%sBin %d (%.02f, %.02f 2$\Theta$) $\Delta2\Theta$ vs. $\Delta\Psi$. Showing %d of %d refls'%(tag,i,bin_low,bin_high,len(a.select(sel)),len(a)))
          plt.xlabel(r'$\Delta\Psi \circ$')
          plt.ylabel(r'$\Delta2\Theta \circ$')

      if params.plots.delta2theta_vs_2theta_2dhist:
        # Plot delta 2theta vs. 2theta
        a = reflections['two_theta_obs']#[:71610]
        b = reflections['two_theta_obs'] - reflections['two_theta_cal']
        fig = plt.figure()
        limits = -0.10, 0.10
        sel = (b > limits[0]) & (b < limits[1])
        plt.hist2d(a.select(sel), b.select(sel), bins=100, range=((0,45), limits))
        plt.clim((0,400))
        cb = plt.colorbar()
        cb.set_label("N reflections")
        plt.title(r'%s$\Delta2\Theta$ vs. 2$\Theta$. Showing %d of %d refls'%(tag,len(a.select(sel)),len(a)))
        plt.xlabel(r'2$\Theta \circ$')
        plt.ylabel(r'$\Delta2\Theta \circ$')

        # calc the trendline
        z = np.polyfit(a.select(sel), b.select(sel), 1)
        if params.verbose: print('y=%.7fx+(%.7f)'%(z[0],z[1]))

      if params.plots.deltaPsi_vs_2theta_2dhist:
        # Plot delta psi vs. 2theta
        from matplotlib.colors import LogNorm
        x = reflections['two_theta_obs'].as_numpy_array()
        y = (reflections['delpsical.rad']*180/math.pi).as_numpy_array()
        fig = plt.figure()
        plt.hist2d(x, y, bins=100, range=((0,45), (-1,1)), norm=LogNorm())
        cb = plt.colorbar()
        cb.set_label("N reflections")
        plt.title(r'%s$\Delta\Psi$ vs. 2$\Theta$. %d refls'%(tag,len(x)))
        plt.xlabel(r'2$\Theta \circ$')
        plt.ylabel(r'$\Delta\Psi \circ$')

      if params.plots.grouped_stats:
        # Plots with single values per panel
        detector_plot_dict(self.params, detector, refl_counts, u"%s N reflections"%t, u"%6d", show=False)
        detector_plot_dict(self.params, detector, rmsds, "%s Positional RMSDs (microns)"%t, u"%4.1f", show=False)
        detector_plot_dict(self.params, detector, radial_rmsds, "%s Radial RMSDs (microns)"%t, u"%4.1f", show=False)
        detector_plot_dict(self.params, detector, transverse_rmsds, "%s Transverse RMSDs (microns)"%t, u"%4.1f", show=False)
        detector_plot_dict(self.params, detector, ttdpcorr, r"%s $\Delta2\Theta$ vs. $\Delta\Psi$ CC"%t, u"%5.3f", show=False)

      if params.plots.unit_cell_histograms: self.plot_unitcells(experiments)
      if params.plots.stats_by_2theta:      self.plot_data_by_two_theta(reflections, tag)

      if params.plots.stats_by_panelgroup:
        # Plot data by panel group
        sorted_values = sorted(pg_bc_dists.values())
        vdict = {}
        for k in pg_bc_dists:
          vdict[pg_bc_dists[k]] = k
        sorted_keys = [vdict[v] for v in sorted_values if vdict[v] in rmsds]
        pg_bc_dists_keylist = list(pg_bc_dists.keys())
        x = [sorted_values[i] for i in range(len(sorted_values)) if pg_bc_dists_keylist[i] in rmsds]

        self.plot_multi_data(x,
                             [[pg_refls_count_d[k] for k in sorted_keys],
                              ([rmsds[k] for k in sorted_keys],
                               [radial_rmsds[k] for k in sorted_keys],
                               [transverse_rmsds[k] for k in sorted_keys]),
                              [radial_rmsds[k]/transverse_rmsds[k] for k in sorted_keys],
                              [mean_delta_two_theta[k] for k in sorted_keys]],
                             "Panel group distance from beam center (mm)",
                             ["N reflections",
                              ("Overall RMSD",
                               "Radial RMSD",
                               "Transverse RMSD"),
                              "R/T RMSD ratio",
                              "Delta two theta"],
                             ["N reflections",
                              "RMSD (microns)",
                              "R/T RMSD ratio",
                              "Delta two theta (degrees)"],
                             "%sData by panelgroup"%tag)

      # Trumpet plot
      if params.plots.trumpet_plot:
        expt_id = min(set(reflections['id']))
        refls = reflections.select(reflections['id'] == expt_id)
        trumpet_plot(experiments[expt_id], refls)

      if params.plots.ewald_offset_plot:
        n_bins = 10
        all_offsets = flex.double()
        all_twothetas = flex.double()

        binned_offsets = [flex.double() for _ in range(n_bins)]
        binned_isigi = [flex.double() for _ in range(n_bins)]
        for expt_id in range(len(experiments)):
          refls = reflections.select(reflections['id'] == expt_id)
          s1 = refls['s1']
          s0 = flex.vec3_double(len(refls), experiments[expt_id].beam.get_s0())
          q = experiments[expt_id].crystal.get_A() * refls['miller_index'].as_vec3_double()
          wavelength = experiments[expt_id].beam.get_wavelength()
          offset = (q+s0).norms() - (1/wavelength)
          two_thetas = refls['two_theta_cal']
          if expt_id == 0:
            fig = plt.figure()
            plt.scatter(two_thetas, offset)
            plt.title(u"%d: Ewald offset ($\AA^{-1}$) vs $2\\theta$ on %d spots"%(expt_id, len(two_thetas)))
            plt.xlabel(u"$2\\theta (\circ)$")
            plt.ylabel(u"Ewald offset ($\AA^{-1}$)")

          array = refls.as_miller_array(experiments[expt_id])
          binner = array.setup_binner(d_min=2.0, n_bins=n_bins)
          isigi = refls['intensity.sum.value']/flex.sqrt(refls['intensity.sum.variance'])
          for bin_id, bin_number in enumerate(binner.range_used()):
            sel = binner.selection(bin_number)
            offsetsel = offset.select(sel)
            all_offsets.extend(offsetsel)
            all_twothetas.extend(two_thetas.select(sel))
            binned_offsets[bin_id].extend(offsetsel)
            binned_isigi[bin_id].extend(isigi.select(sel))

        fig = plt.figure()
        legend = []
        h = flex.histogram(all_offsets, n_slots=20)
        for bin_id, bin_number in enumerate(binner.range_used()):
          legend.append("%5.2f-%5.2f"%binner.bin_d_range(bin_number))
          x, y = [], []
          for slot_id, slot in enumerate(h.slot_infos()):
            sel = (binned_offsets[bin_id] >= slot.low_cutoff) & (binned_offsets[bin_id] <= slot.high_cutoff)
            x.append(slot.center())
            y.append(flex.median(binned_isigi[bin_id].select(sel)) if sel.count(True) else 0)
            print (bin_id, binner.bin_d_range(bin_number), sel.count(True), x[-1], y[-1])
          plt.plot(x, y, '-')
        plt.legend(legend)
        plt.title(u"Binned $I/\sigma_I$ vs. Ewald offset")
        plt.xlabel(u"Ewald offset ($\AA^{-1}$)")
        plt.ylabel(u"Median $I/\sigma_I$")

        plt.figure()
        plt.title('Ewald offsets vs. two theta')
        plt.hist2d(all_twothetas.as_numpy_array(), all_offsets.as_numpy_array(), bins=100)

      if self.params.save_pdf:
        pp = PdfPages('residuals_%s.pdf'%(tag.strip()))
        for i in plt.get_fignums():
          pp.savefig(plt.figure(i))
        pp.close()
      elif self.params.save_png:
        if len(tag) == 0:
          prefix = ""
        else:
          prefix = "%s_"%tag.strip()
        for i in plt.get_fignums():
          print("Saving figure", i)
          plt.figure(i).savefig("%sfig%02d.png"%(prefix, i))
      else:
        plt.show()

  def plot_multi_data(self, x, ygroups, xlabel, legends, ylabels, title):
    fig = plt.figure()
    from mpl_toolkits.axes_grid1 import host_subplot
    import mpl_toolkits.axisartist as AA

    host = host_subplot(111, axes_class=AA.Axes)
    plt.subplots_adjust(right=0.75)

    y1 = ygroups[0]
    ygroups = ygroups[1:]
    current_offset = 0
    offset = 60

    host.set_title(title)
    host.set_xlabel(xlabel)
    host.set_ylabel(ylabels[0])
    p1, = host.plot(x, y1, label=legends[0])
    host.axis["left"].label.set_color(p1.get_color())

    for y, legend, ylabel in zip(ygroups, legends[1:], ylabels[1:]):
      par = host.twinx()
      new_fixed_axis = par.get_grid_helper().new_fixed_axis
      par.axis["right"] = new_fixed_axis(loc="right",
                                         axes=par,
                                         offset=(current_offset, 0))
      par.axis["right"].toggle(all=True)
      current_offset += offset

      if isinstance(y, tuple):
        for data, l in zip(y, legend):
          p, = par.plot(x, data, label=l)
      else:
        p, = par.plot(x, y, label=legend)
        par.axis["right"].label.set_color(p.get_color())
      par.set_ylabel(ylabel)

    host.legend(loc='upper center', ncol=2, fancybox=True, shadow=True, fontsize=8)
    plt.draw()

  def plot_data_by_two_theta(self, reflections, tag):
    n_bins = 30
    arbitrary_padding = 1
    sorted_two_theta = flex.sorted(reflections['two_theta_obs'])
    bin_low = [sorted_two_theta[int((len(sorted_two_theta)/n_bins) * i)] for i in range(n_bins)]
    bin_high = [bin_low[i+1] for i in range(n_bins-1)]
    bin_high.append(sorted_two_theta[-1]+arbitrary_padding)

    title = "%sBinned data by two theta (n reflections per bin: %.1f)"%(tag, len(sorted_two_theta)/n_bins)

    x = flex.double()
    x_centers = flex.double()
    n_refls = flex.double()
    rmsds = flex.double()
    radial_rmsds = flex.double()
    transverse_rmsds = flex.double()
    rt_ratio = flex.double()
    #delta_two_theta = flex.double()
    rmsd_delta_two_theta = flex.double()

    for i in range(n_bins):
      x_centers.append(((bin_high[i]-bin_low[i])/2) + bin_low[i])
      refls = reflections.select((reflections['two_theta_obs'] >= bin_low[i]) & (reflections['two_theta_obs'] < bin_high[i]))
      n = len(refls)
      n_refls.append(n)
      if n == 0:
        rmsds.append(0)
        radial_rmsds.append(0)
        transverse_rmsds.append(0)
        rt_ratio.append(0)
        rmsd_delta_two_theta.append(0)
      else:
        rmsds.append(1000*math.sqrt(flex.sum_sq(refls['difference_vector_norms'])/n))
        radial_rmsds.append(1000*math.sqrt(flex.sum_sq(refls['radial_displacements'])/n))
        transverse_rmsds.append(1000*math.sqrt(flex.sum_sq(refls['transverse_displacements'])/n))
        rt_ratio.append(radial_rmsds[-1]/transverse_rmsds[-1])
        rmsd_delta_two_theta.append(math.sqrt(flex.sum_sq(refls['two_theta_obs']-refls['two_theta_cal'])/n))
      #delta_two_theta.append(flex.mean(refls['two_theta_obs']-refls['two_theta_cal']))
    assert len(reflections) == flex.sum(n_refls)

    self.plot_multi_data(x_centers,
                         [rt_ratio, (rmsds, radial_rmsds, transverse_rmsds), rmsd_delta_two_theta],
                         "Two theta (degrees)",
                         ["R/T RMSD",
                          ("RMSD","R RMSD","T RMSD"),
                          r"RMSD $\Delta2\Theta$"],
                         ["R/T RMSD ratio",
                          "Overall, radial, transverse RMSD (microns)",
                          r'$\Delta2\Theta RMSD (\circ)$'],
                         title)

  def histogram(self, reflections, title, plots = True, verbose = True):
    data = reflections['difference_vector_norms']
    n_slots = 100
    if self.params.residuals.histogram_max is None:
      h = flex.histogram(data, n_slots=n_slots)
    else:
      h = flex.histogram(data.select(data <= self.params.residuals.histogram_max), n_slots=n_slots)

    n = len(reflections)
    rmsd = math.sqrt((reflections['xyzcal.mm']-reflections['xyzobs.mm.value']).sum_sq()/n)
    sigma = mode = h.slot_centers()[list(h.slots()).index(flex.max(h.slots()))]
    mean = flex.mean(data)
    median = flex.median(data)
    if verbose: print("RMSD (microns)", rmsd * 1000)
    if verbose: print("Histogram mode (microns):", mode * 1000)
    if verbose: print("Overall mean (microns):", mean * 1000)
    if verbose: print("Overall median (microns):", median * 1000)
    mean2 = math.sqrt(math.pi/2)*sigma
    rmsd2 = math.sqrt(2)*sigma
    if verbose: print("Rayleigh Mean (microns)", mean2 * 1000)
    if verbose: print("Rayleigh RMSD (microns)", rmsd2 * 1000)

    r = reflections['radial_displacements']
    t = reflections['transverse_displacements']
    if verbose: print("Overall radial RMSD (microns)", math.sqrt(flex.sum_sq(r)/len(r)) * 1000)
    if verbose: print("Overall transverse RMSD (microns)", math.sqrt(flex.sum_sq(t)/len(t)) * 1000)

    if not plots: return

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(h.slot_centers().as_numpy_array(), h.slots().as_numpy_array(), '-')

    vmax = self.params.residuals.plot_max
    if self.params.residuals.histogram_xmax is not None:
      ax.set_xlim((0,self.params.residuals.histogram_xmax))
    if self.params.residuals.histogram_ymax is not None:
      ax.set_ylim((0,self.params.residuals.histogram_ymax))
    plt.title(title)


    ax.plot((mean, mean), (0, flex.max(h.slots())), 'g-')
    ax.plot((mean2, mean2), (0, flex.max(h.slots())), 'g--')
    ax.plot((mode, mode), (0, flex.max(h.slots())), 'r-')
    ax.plot((rmsd, rmsd), (0, flex.max(h.slots())), 'b-')
    ax.plot((rmsd2, rmsd2), (0, flex.max(h.slots())), 'b--')

    ax.legend([r"$\Delta$XY", "MeanObs", "MeanRayl", "Mode", "RMSDObs", "RMSDRayl"])
    ax.set_xlabel("(mm)")
    ax.set_ylabel("Count")

  def image_rmsd_histogram(self, reflections, tag, boxplot = True):
    data = flex.double()
    for i in set(reflections['id']):
      refls = reflections.select(reflections['id']==i)
      if len(refls) == 0:
        continue
      rmsd = math.sqrt(flex.sum_sq(refls['difference_vector_norms'])/len(refls))
      data.append(rmsd)
    data *= 1000
    h = flex.histogram(data, n_slots=40)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(h.slot_centers().as_numpy_array(), h.slots().as_numpy_array(), '-')
    plt.title("%sHistogram of image RMSDs"%tag)
    ax.set_xlabel("RMSD (microns)")
    ax.set_ylabel("Count")

    if not boxplot: return

    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.boxplot(data, vert=False)
    plt.title("%sBoxplot of image RMSDs"%tag)
    ax.set_xlabel("RMSD (microns)")

  def detector_plot_refls(self, detector, reflections, title, show=True, plot_callback=None, colorbar_units=None, new_fig = True):
    """
    Use matplotlib to plot a detector, color coding panels according to callback
    @param detector detector reference detector object
    @param title title string for plot
    @param units_str string with a formatting statment for units on each panel
    """
    if new_fig:
      fig = plt.figure()
      ax = fig.add_subplot(111, aspect='equal')
    else:
      fig = plt.gcf()
      ax = plt.gca()
    max_dim = 0
    for panel_id, panel in enumerate(detector):
      # get panel coordinates
      size = panel.get_image_size()
      p0 = col(panel.get_pixel_lab_coord((0,0)))
      p1 = col(panel.get_pixel_lab_coord((size[0]-1,0)))
      p2 = col(panel.get_pixel_lab_coord((size[0]-1,size[1]-1)))
      p3 = col(panel.get_pixel_lab_coord((0,size[1]-1)))
      bounds = (p0[0:2],p1[0:2],p2[0:2],p3[0:2])

      v1 = p1-p0
      v2 = p3-p0
      vcen = ((v2/2) + (v1/2)) + p0

     # add the panel to the plot
      ax.add_patch(Polygon(bounds, closed=True, fill=False))
      if self.params.panel_numbers:
        ax.annotate("%d"%(panel_id), vcen[0:2], ha='center')

      # find the plot maximum dimensions
      for p in [p0, p1, p2, p3]:
        for c in p[0:2]:
          if abs(c) > max_dim:
            max_dim = abs(c)

      panel_refls = reflections.select(reflections['panel'] == panel_id)
      if len(panel_refls) == 0:
        sm = color_vals = None
        continue

      if plot_callback is None:
        sm = color_vals = None
      else:
        result = plot_callback(panel_refls, panel, ax, bounds)
        if result is not None:
          sm, color_vals = result
        else:
          sm = color_vals = None

    # plot the results
    ax.set_xlim((-max_dim,max_dim))
    ax.set_ylim((-max_dim,max_dim))
    ax.set_xlabel("mm")
    ax.set_ylabel("mm")
    if sm is not None and color_vals is not None:
      if colorbar_units is None:
        colorbar_units = "mm"
      cb = ax.figure.colorbar(sm, ticks=color_vals)
      cb.ax.set_yticklabels(["%3.2f %s"%(i,colorbar_units) for i in color_vals])

    plt.title(title)

    if show:
      plt.show()

if __name__ == '__main__':
  with show_mail_on_error():
    script = Script()
    script.run()
