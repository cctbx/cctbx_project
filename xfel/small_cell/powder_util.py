from __future__ import division
from dxtbx.model.experiment_list import DetectorComparison
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import numpy as np
import copy
from dials.array_family import flex
from cctbx import uctbx
from scitbx.math import five_number_summary
from cctbx.crystal import symmetry
import cctbx.miller
import itertools

def angle(v1, v2):
    """
    Compute the angle between two cartesian vectors.

    Parameters:
    v1 (numpy.ndarray): The first vector.
    v2 (numpy.ndarray): The second vector.

    Returns:
    float: The angle between the vectors in degrees.
    """
    dot_product = np.dot(v1, v2)
    magnitude_v1 = np.linalg.norm(v1)
    magnitude_v2 = np.linalg.norm(v2)
    cos_theta = dot_product / (magnitude_v1 * magnitude_v2)
    return np.degrees(np.arccos(cos_theta))

def create_pairwise_plots(points, labels):
    # Create a figure with 3 subplots
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))

    # Define custom color map
    color_map = {
        -1: 'lightgray',  # outliers in light gray
        0: 'red',
        1: 'blue',
        2: 'orange',
        3: 'green'
    }

    # Convert labels to colors
    colors = [color_map[label] if label in color_map else 'gray' for label in labels]

    # Define scatter plot properties
    scatter_kwargs = {'c': labels, 'cmap': 'viridis', 's': 1, 'alpha': 0.6}
    scatter_kwargs = {'c': colors,  's': 1, 'alpha': 0.6}

    # Plot ab
    ax1.scatter(points[:, 0], points[:, 1], **scatter_kwargs)
    ax1.set_xlabel('a')
    ax1.set_ylabel('b')
    ax1.set_title('a vs b')

    # Plot ac
    ax2.scatter(points[:, 0], points[:, 2], **scatter_kwargs)
    ax2.set_xlabel('a')
    ax2.set_ylabel('c')
    ax2.set_title('a vs c')

    # Plot bc
    ax3.scatter(points[:, 1], points[:, 2], **scatter_kwargs)
    ax3.set_xlabel('b')
    ax3.set_ylabel('c')
    ax3.set_title('b vs c')

    plt.tight_layout()
    plt.show()

class Spotfinder_radial_average:

  def __init__(self, experiments, reflections, params):
    self.reflections = reflections
    self.experiments = experiments
    self.params = params
    self.n_panels = len(experiments[0].detector)
    self.panelsums = [np.zeros(params.n_bins) for _ in range(self.n_panels)]
    self.filtered_panelsums = [
        np.zeros(params.n_bins) for _ in range(self.n_panels)
    ]
    self.antifiltered_panelsums = [
        np.zeros(params.n_bins) for _ in range(self.n_panels)
    ]
    self.expt_count = 0
    self.filtered_expt_count = 0
    self.antifiltered_expt_count = 0

  def _process_pixel(self, i_panel, s0, panel, xy, value):
    value -= self.params.downweight_weak
    d_max_inv = 1/self.params.d_max
    d_min_inv = 1/self.params.d_min
    res_inv = 1 / panel.get_resolution_at_pixel(s0, xy)
    res = 1/res_inv
    self.dvals.append(res)
    if self.params.filter.enable:
      for i, (dmax, dmin) in enumerate(zip(self.d_max_vals, self.d_min_vals)):
        if dmax > res > dmin:
          self.filter_counts[i] += 1
    n_bins = self.params.n_bins
    i_bin = int(
        n_bins * (res_inv - d_max_inv ) / (d_min_inv - d_max_inv)
        )
    if 0 <= i_bin < n_bins:
      self.current_panelsums[i_panel][i_bin] += value
    return res

  def _nearest_peak(self, x, xvalues, yvalues):
    i = np.searchsorted(xvalues, x, side="left")

    #Exclude (before) first and (after) last points
    if i < 1 or i >= len(xvalues):
      print ("Not a valid peak.")
      return None
    #Exclude troughs
    if yvalues[i-1] >= yvalues[i] and yvalues[i+1] >= yvalues[i]:
      print ("Not a valid peak.")
      return None
    #find the highest nearby yvalue
    direction = 1 if yvalues[i+1] > yvalues[i] else -1
    while yvalues[i+direction] > yvalues[i]:
      i += direction
    return 1/xvalues[i]

  def calculate(self):
    params = self.params

    # setup limits and bins
    n_bins = params.n_bins
    d_max, d_min = params.d_max, params.d_min
    d_inv_low, d_inv_high = 1/d_max, 1/d_min
    unit_wt = (params.peak_weighting == "unit")
    refls = self.reflections
    expts = self.experiments
    self.dvals = []

    #apply beam center correction to expts
    detector = expts[0].detector
    if not np.allclose(params.xyz_offset, [0,0,0]):
      assert params.single_detector
      ref_detector = copy.deepcopy(detector)
      hierarchy = detector.hierarchy()
      fast = hierarchy.get_local_fast_axis()
      slow = hierarchy.get_local_slow_axis()
      origin = hierarchy.get_local_origin()
      corrected_origin = (
              origin[0] + params.xyz_offset[0],
              origin[1] + params.xyz_offset[1],
              origin[2] + params.xyz_offset[2]
              )
      hierarchy.set_local_frame(fast, slow, corrected_origin)
    else:
      ref_detector = detector
    if params.single_detector:
      compare_detector = DetectorComparison()
      for expt in expts:
        if expt.detector is detector: continue
        assert compare_detector(ref_detector, expt.detector)
        expt.detector = detector

    if params.angle_histogram.enable:
      if 's1' not in refls.keys():
        refls.centroid_px_to_mm(expts)
        refls.map_centroids_to_reciprocal_space(expts)
      angles_12 = []
      angles_13 = []
      angles_23 = []
      detplot_counter = 0

    for i, expt in enumerate(expts):
      self.current_panelsums = [
          np.zeros(params.n_bins) for _ in range(self.n_panels)
      ]
      if self.params.filter.enable:
        assert self.params.filter.d_vals and len(self.params.filter.d_vals)%2==0
        self.filter_counts = [0 for _ in range(len(self.params.filter.d_vals)//2)]
        self.d_max_vals = self.params.filter.d_vals[::2]
        self.d_min_vals = self.params.filter.d_vals[1::2]

      if self.params.filter.enable:
        self.use_current_expt = False
      else:
        self.use_current_expt = True
      if i % 1000 == 0: print("experiment ", i)
      if self.params.n_max is not None and i>self.params.n_max:
        break
      s0 = expt.beam.get_s0()
      sel = refls['id'] == i
      refls_sel = refls.select(sel)
      xyzobses = refls_sel['xyzobs.px.value']
      intensities = refls_sel['intensity.sum.value']
      panels = refls_sel['panel']
      shoeboxes = refls_sel['shoebox']

      if params.angle_histogram.enable: 
        r1max, r1min = params.angle_histogram.range1
        r2max, r2min = params.angle_histogram.range2
        r3max, r3min = params.angle_histogram.range3

        i_r1, i_r2 ,i_r3 = [],[],[]
      for i_refl in range(len(refls_sel)):
        self.expt_count += 1
        i_panel = panels[i_refl]
        panel = expt.detector[i_panel]

        if params.peak_position=="xyzobs":
          xy = xyzobses[i_refl][0:2]
          if params.peak_weighting == "intensity":
            value = intensities[i_refl]
          else:
            value = 1
          res = self._process_pixel(i_panel, s0, panel, xy, value)
          if params.angle_histogram.enable and r1max > res > r1min:
            i_r1.append(i_refl)
          if params.angle_histogram.enable and r2max > res > r2min:
            i_r2.append(i_refl)
          if params.angle_histogram.enable and r3max > res > r3min:
            i_r3.append(i_refl)

      if params.angle_histogram.enable and i_r1 and i_r2 and i_r3:
        i_all = i_r1 + i_r2 + i_r3
        a12, a13, a23 = [],[],[]
        subsel_mask = flex.bool([n in i_all for n in range(len(refls_sel))])
        subsel_mask_inv = flex.bool([not x for x in subsel_mask])

#        subsel = refls_sel.select(subsel_mask)
#        subsel_inv = refls_sel.select(subsel_mask_inv)
#        xyz1 = np.array(subsel['xyzobs.mm.value'])
#        xyz2 = np.array(subsel_inv['xyzobs.mm.value'])
#        plt.scatter(xyz1[:,0], xyz1[:,1], c='red', s=2)
#        plt.scatter(xyz2[:,0], xyz2[:,1], c='blue', s=2)
#        bc = expt.detector[0].get_beam_centre(expt.beam.get_s0())
#        plt.scatter(*bc, c='k', s=5)
#        plt.xlim((100,240))
#        plt.ylim((100,240))
          
        s0 = np.array(expt.beam.get_s0())
        for i1, i2, i3 in itertools.product(i_r1, i_r2, i_r3):
          v1 = np.array(refls_sel[i1]['s1']) - s0
          v2 = np.array(refls_sel[i2]['s1']) - s0
          v3 = np.array(refls_sel[i3]['s1']) - s0
          a12.append(angle(v1, v2))
          a13.append(angle(v1, v3))
          a23.append(angle(v2, v3))
#        print('\n----------------------')
#        print('Pairwise angles:')
#        headers = '1,2: 59, 121', '1,3: 30, 150', '2,3: 25, 155'
#        for header, vals in zip(headers, (a12, a13, a23)):
#          print()
#          print(header)
#          for v in vals: print(round(v, 2))
#        plt.show()
        angles_12.extend(a12)
        angles_13.extend(a13)
        angles_23.extend(a23)
#        for i1, i2 in itertools.product(i_r1, i_r2):
#          v1 = np.array(refls_sel[i1]['s1']) - s0
#          v2 = np.array(refls_sel[i2]['s1']) - s0
#          angles_12.append(angle(v1, v2))
#        for i1, i2 in itertools.product(i_r1, i_r3):
#          v1 = np.array(refls_sel[i1]['s1']) - s0
#          v2 = np.array(refls_sel[i2]['s1']) - s0
#          angles_13.append(angle(v1, v2))
#        for i1, i2 in itertools.product(i_r2, i_r3):
#          v1 = np.array(refls_sel[i1]['s1']) - s0
#          v2 = np.array(refls_sel[i2]['s1']) - s0
#          angles_23.append(angle(v1, v2))

      for i in range(len(self.panelsums)):
        self.panelsums[i] = self.panelsums[i] + self.current_panelsums[i]
      if self.params.filter.enable and self.params.filter.select_mode=='any':
        use_current_expt = any(self.filter_counts)
      elif self.params.filter.enable and self.params.filter.select_mode=='all':
        use_current_expt = all(self.filter_counts)
      else:
        use_current_expt = True
      if use_current_expt:
        self.filtered_expt_count += 1
        for i in range(len(self.panelsums)):
          self.filtered_panelsums[i] = \
              self.filtered_panelsums[i] + self.current_panelsums[i]
      else:
        self.antifiltered_expt_count += 1
        for i in range(len(self.panelsums)):
          self.antifiltered_panelsums[i] = \
              self.antifiltered_panelsums[i] + self.current_panelsums[i]
    if params.angle_histogram.enable:
#      from sklearn.cluster import DBSCAN
#      data=np.vstack((angles_12, angles_13, angles_23)).transpose()
#      dbscan = DBSCAN(eps=6, min_samples=15)
#      labels = dbscan.fit_predict(data)
#      create_pairwise_plots(data, labels)
#      fig, (ax1, ax2, ax3) = plt.subplots(1,3)
#      ax1.scatter(angles_12, angles_13, s=.5)
#      ax2.scatter(angles_13, angles_23, s=.5)
#      ax3.scatter(angles_12, angles_23, s=.5)
      fig, (ax1, ax2, ax3) = plt.subplots(3,1)
      ax1.hist(angles_12, bins=180)
      ax2.hist(angles_13, bins=180)
      ax3.hist(angles_23, bins=180)

      plt.show()
    self.dvals = np.array(self.dvals)


  def plot(self):
    params = self.params
    d_max_inv = 1/params.d_max
    d_min_inv = 1/params.d_min
    xvalues = np.linspace(d_max_inv, d_min_inv, params.n_bins)
    fig, ax = plt.subplots()

    ps_maxes = [max(ps) for ps in self.panelsums]
    ps_max = max(ps_maxes)

    if params.split_panels:
      offset = 0.5*ps_max
      for i_sums, sums in enumerate(self.panelsums):
        yvalues = np.array(sums)
        plt.plot(xvalues, yvalues+0.5*i_sums*offset)
    elif params.filter.enable and params.filter.plot_mode=="ratio":
      print('filtered: ', self.filtered_expt_count)
      print('antifiltered: ', self.antifiltered_expt_count)
      for x in self.filtered_panelsums:
        x /= self.filtered_expt_count
      for x in self.antifiltered_panelsums:
        x /= self.antifiltered_expt_count
      yvalues = sum(self.filtered_panelsums) - sum(self.antifiltered_panelsums)
      plt.plot(xvalues, yvalues)
    elif params.filter.enable and params.filter.plot_mode=="simple":
      # same as below, but keeping this separate for flexibility
      yvalues = sum(self.filtered_panelsums)
      plt.plot(xvalues, yvalues)
    else:
      yvalues = sum(self.filtered_panelsums)
      plt.plot(xvalues, yvalues)

    if params.augment:
      plt.plot(*augment(self.experiments, self.reflections, params.d_min, params.d_max))
    ax.set_xlim(d_max_inv, d_min_inv)
#    ax.get_xaxis().set_major_formatter(tick.FuncFormatter(
#      lambda x, _: "{:.3f}".format(1/x)))

    if params.output.xy_file:
      with open(params.output.xy_file, 'w') as f:
        for x,y in zip(xvalues, yvalues):
          f.write("{:.6f}\t{}\n".format(1/x, y))

    # Now plot the predicted peak positions if requested
    if params.unit_cell or params.space_group:
      assert params.unit_cell and params.space_group
      sym = symmetry(
          unit_cell=params.unit_cell, space_group=params.space_group.group()
      )
      hkl_list = cctbx.miller.build_set(sym, True, d_min=params.d_min)
      dspacings = params.unit_cell.d(hkl_list.indices())
#      for hkl, d in sorted(zip(hkl_list.indices(), dspacings), key=lambda x:x[1]):
#        print('{:.3f}: {}'.format(d, hkl))
      dspacings_inv = 1/dspacings
      pplot_min = -.05*max(yvalues)
      for d in dspacings_inv:
        plt.plot((d,d),(pplot_min,0), 'r-', linewidth=1)

    if params.output.plot_file:
      plt.savefig(params.output.plot_file)

    if params.plot.interactive and params.output.peak_file:
      backend_list = ["TkAgg","QtAgg"]
      assert (plt.get_backend() in backend_list), """Matplotlib backend not compatible with interactive peak picking.
You can set the MPLBACKEND environment varibale to change this.
Currently supported options: %s""" %backend_list
      #If a peak list output file is specified, do interactive peak picking:
      with open(params.output.peak_file, 'w') as f:
        vertical_line = ax.axvline(color='r', lw=0.8, ls='--', x=xvalues[1])
        vertical_line.set_visible(False)
        def onmove(event):
          if fig.canvas.toolbar.mode: return
          x = event.xdata
          vertical_line.set_xdata(x)
          if plt.getp(vertical_line, 'visible'):
            ax.figure.canvas.draw()
        def onclick(event):
          if fig.canvas.toolbar.mode: return
          self.d1 = 1/event.xdata
          vertical_line.set_visible(True)
        def onrelease(event):
          if fig.canvas.toolbar.mode: return
          self.d2 = 1/event.xdata
          vertical_line.set_visible(False)
          left = max(self.d1, self.d2)
          right = min(self.d1, self.d2)
          if left==right:
            peak = left
          else:
            matching_dvals = self.dvals[
                np.logical_and(self.dvals<left, self.dvals>right)
            ]
            peak = np.median(matching_dvals)
          ax.figure.canvas.draw()
          print('Median=%f, writing to %s.' % (peak, params.output.peak_file))
          f.write(str(peak)+"\n")

        mmv = fig.canvas.mpl_connect('motion_notify_event', onmove)
        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        ciu = fig.canvas.mpl_connect('button_release_event', onrelease)

        plt.show()

    elif params.plot.interactive:
      import os
      plt.title(os.getcwd())
      plt.show()


class Center_scan:
  def __init__(self, experiments, reflections, params):
    self.params = params

    assert params.single_detector
    self.reflections = reflections
    self.experiments = experiments
    assert len(self.experiments.detectors())==1

    self.net_origin_shift = np.array([0.,0.,0.])
    self.centroid_px_mm_done = False
    self.px_size = self.experiments.detectors()[0][0].get_pixel_size()
    self.target_refl_count = 0
    self.prune_refls()


  def prune_refls(self, margin=0.02):
    self.refls_pruned = self.reflections
    self.update_dvals()
    d_min_inv = 1/self.params.center_scan.d_min + margin
    d_max_inv = max(1/self.params.center_scan.d_max - margin, 0.001)
    d_min = 1/d_min_inv
    d_max = 1/d_max_inv
    sel_lt = self.reflections['d'] < d_max
    sel_gt = self.reflections['d'] > d_min
    sel = sel_lt & sel_gt
    self.refls_pruned = self.reflections.select(sel)

  def update_dvals(self):
    refls = self.refls_pruned
    if not self.centroid_px_mm_done:
      refls.centroid_px_to_mm(self.experiments)
      self.centroid_px_mm_done = True
    refls.map_centroids_to_reciprocal_space(self.experiments)
    d_star_sq = flex.pow2(refls['rlp'].norms())
    refls['d'] = uctbx.d_star_sq_as_d(d_star_sq)
    sel_lt = refls['d'] < self.params.center_scan.d_max
    sel_gt = refls['d'] > self.params.center_scan.d_min
    sel = sel_lt & sel_gt
    count = sel.count(True)
    if count > self.target_refl_count:
      self.target_refl_count = count
    self.dvals = refls.select(sel)['d']

  def width(self):
    if len(self.dvals) < 3: return 999
    _, q1, _, q3, _ = five_number_summary(self.dvals)
    iqr = q3 - q1
    sel_lt = self.dvals < q3 + 1.5*iqr
    sel_gt = self.dvals > q1 - 1.5*iqr
    sel = sel_lt & sel_gt
    if sel.count(True) < 0.8*self.target_refl_count:
      return 999
    else:
      result = self.dvals.select(sel).sample_standard_deviation()
      print(f'width {result:.5f} from {sel.count(True)} dvals')
      return result

  def search_step(self, step_px, nsteps=3, update=True):
    step_size_mm = np.array(self.px_size + (0.,)) * step_px
    assert nsteps%2 == 1, "nsteps should be odd"
    step_min = -1 * (nsteps // 2)
    step_max = nsteps//2 + 1e-6 # make the range inclusive
    step_arange = np.arange(step_min, step_max)
    steps = np.array([(x, y, 0) for x in step_arange for y in step_arange])
    steps_mm = steps * step_size_mm

    detector = self.experiments.detectors()[0]
    hierarchy = detector.hierarchy()
    fast = hierarchy.get_local_fast_axis()
    slow = hierarchy.get_local_slow_axis()
    origin = hierarchy.get_local_origin()

    results = []
    self.update_dvals()
    width_start = self.width()
    print(f'start: {width_start:.5f}')

    for step_mm in steps_mm:
      new_origin = step_mm + origin
      hierarchy.set_local_frame(fast, slow, new_origin)
      self.update_dvals()
      result = self.width()
      results.append(result)

    width_end = min(results)
    i_best = results.index(width_end)
    origin_shift = steps_mm[i_best]
    if update:
      new_origin = origin + origin_shift
      self.net_origin_shift += origin_shift
    else:
      new_origin = origin
    hierarchy.set_local_frame(fast, slow, new_origin)
    print(f'step: {origin_shift}')
    print(f'end: {width_end:.5f}')
    print(f'net shift: {self.net_origin_shift}')
    return width_start, width_end, self.net_origin_shift

def augment(expts, refls, d_min, d_max):
  """ Add pairwise 3D spot distances to the d-spacing histogram """
  lab = flex.vec3_double()
  det = expts[0].detector
  filtered = flex.reflection_table()
  for panel_id, panel in enumerate(det):
    print("Processing panel", panel_id)
    subset = refls.select(refls['panel'] == panel_id)
    mm = panel.pixel_to_millimeter(flex.vec2_double(*subset['xyzobs.px.value'].parts()[0:2]))
    lab.extend(panel.get_lab_coord(mm))
    filtered.extend(subset)
  refls = filtered

  s0 = flex.vec3_double(len(refls))
  for expt_id, expt in enumerate(expts):
    if expt_id % 500 == 0:
      print("Processing experiment", expt_id)
    sel = refls['id'] == expt_id
    s0.set_selected(sel, expt.beam.get_s0())

  s1 = lab / lab.norms() * s0.norms()
  rlp = s1 - s0

  refls['rlp'] = rlp
  refls['d'] = 1/rlp.norms()

  sel = (refls['d'] >= d_min) & (refls['d'] <= d_max)
  refls = refls.select(sel)
  print("After filtering by resolution,", len(refls), "reflections remain")

  _, x_sf = np.histogram((1/refls['d']).as_numpy_array(), bins=2000, range = (1/d_max,1/d_min))
  y = None

  for expt_id, expt in enumerate(expts):
    if expt_id % 500 == 0:
      print("Processing experiment", expt_id)
    subset = refls.select(refls['id'] == expt_id)
    if len(subset) <= 1: continue
    for i in range(len(subset)):
      diffs = 1/(subset['rlp']-flex.vec3_double(len(subset), subset['rlp'][i])).norms()
      if y is None:
        y = np.histogram((1/diffs).as_numpy_array(), bins=2000, range = (1/d_max,1/d_min))[0]
      else:
        y += np.histogram((1/diffs).as_numpy_array(), bins=2000, range = (1/d_max,1/d_min))[0]

  return x_sf[:-1], y

