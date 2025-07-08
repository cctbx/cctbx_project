
from __future__ import division
import logging
from iotbx.phil import parse
from dials.util.options import ArgumentParser
from dials.util import show_mail_on_error
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
from matplotlib.colors import LogNorm, Normalize
from matplotlib.widgets import RectangleSelector
from dials.array_family import flex

help_message = """
Just gonna see what happens
"""

master_phil = parse(
"""
  plot {
    d_max = 20
      .type = float
      .help = x-axis left limit
    d_min = 2
      .type = float
      .help = x-axis right limit
    d_bins = 1000
      .type = int
      .help = x-axis bins
    th_min = 0
      .type = float
      .help = y-axis lower limit
    th_max = 180
      .type = float
      .help = y-axis upper limit
    th_bins = 180
      .type = int
      .help = y-axis bins
    zscale = *log linear
      .type = choice
  }
  filter {
    ranges = None
      .type = floats
      .help = Resolution ranges (pairs dmin, dmax) to construct angle heatmaps
  }
  n_max = None
    .type = int
    .help = Stop processing after this many experiments
  save_processed_reflections = True
    .type = bool
    .help = Modify the refl file in-place with extra columns.
"""
)
phil_scope = master_phil

def vector_angle(s1_a, s1_b):
  """
  Compute angle in degrees between two scattering vectors
  """
  cos_angle = np.dot(s1_a, s1_b) / (np.linalg.norm(s1_a) * np.linalg.norm(s1_b))
  # Handle numerical instability at bounds
  cos_angle = min(1.0, max(-1.0, cos_angle))
  return np.rad2deg(np.arccos(cos_angle))

class Angle_heatmap:
  def __init__(self, experiments, reflections, params):
    """
    Args:
      experiments: dxtbx ExperimentList
      reflections: dials reflection table
      params: phil scope containing plot and filter parameters
    """
    self.experiments = experiments
    self.reflections = reflections
    self.params = params

    # Validate and parse resolution ranges
    if not params.filter.ranges or len(params.filter.ranges) % 2 != 0:
      raise ValueError("filter.ranges must contain pairs of d_max, d_min values")

    # Convert flat list of floats into list of (d_max, d_min) tuples
    self.resolution_ranges = list(zip(
      params.filter.ranges[::2],  # d_max values
      params.filter.ranges[1::2]  # d_min values
    ))

    # Store the inverse d-limits
    self.plot_s_min = 1/self.params.plot.d_max
    self.plot_s_max = 1/self.params.plot.d_min

    # Validate each range
    for d_max, d_min in self.resolution_ranges:
      if d_max <= d_min:
        raise ValueError(f"Invalid resolution range: d_max ({d_max}) must be greater than d_min ({d_min})")

    # Initialize results dictionary
    self.results = {
      (d_max, d_min): {
        'd_spacings_A': [],
        'd_spacings_B': [],
        'angles': [],
        'i_expt': []
      }
      for d_max, d_min in self.resolution_ranges
    }
    self.powder_data = [[] for _ in range(len(self.experiments))]

  def on_select(self, eclick, erelease, d_range, current_ax):
    """Callback for box selection"""
    x1, y1 = eclick.xdata, eclick.ydata
    x2, y2 = erelease.xdata, erelease.ydata

    # Initialize selection box storage
    if not hasattr(self, '_selection_boxes'):
      self._selection_boxes = []

    # Clear any existing selection boxes
    for box, ax in self._selection_boxes:
      box.remove()
      ax.figure.canvas.draw()
    self._selection_boxes = []

    # Check for single click (reset condition)
    if abs(x1-x2) < 1e-6 and abs(y1-y2) < 1e-6:
      selected_experiments = None  # None means use all data
      title = 'Powder Pattern (all data)'
    else:
      # Get d-spacing and angle ranges from box
      d_min_b, d_max_b = sorted([1/x1, 1/x2])
      th_min, th_max = sorted([y1, y2])

      # Draw selection box on the current heatmap
      box = plt.Rectangle(
        (x1, y1), x2-x1, y2-y1,
        fill=False, color='red', linewidth=1
      )
      current_ax.add_patch(box)
      self._selection_boxes.append((box, current_ax))
      current_ax.figure.canvas.draw()

      # Get experiments that contributed in this range
      d_max_range, d_min_range = d_range  # resolution range for this heatmap
      data = self.results[(d_max_range, d_min_range)]

      mask = ((data['d_spacings_B'] >= d_min_b) &
              (data['d_spacings_B'] <= d_max_b) &
              (data['angles'] >= th_min) &
              (data['angles'] <= th_max))

      selected_experiments = set(np.array(data['i_expt'])[mask])

      positions = self.compute_spot_scatter(
        selected_experiments,
        range_a=(d_max_range, d_min_range),
        range_b=(d_max_b, d_min_b),
        range_theta=(th_min, th_max)
      )
      self.plot_spot_scatter(positions)

      title = f'Powder Pattern (selection: d={d_min_b:.2f}-{d_max_b:.2f}Å, θ={th_min:.1f}-{th_max:.1f}°)'

    # Update powder pattern
    full_hist, edges = self.make_powder_pattern(None)  # Get full pattern
    ymax = full_hist.max()
    selected_hist, _ = self.make_powder_pattern(selected_experiments, rescale_ymax=ymax)  # Get selected pattern

    self.powder_fig.clear()
    ax = self.powder_fig.add_axes([.1, .15, .8, .75]) # [left, bottom, width, height]

    # Plot full pattern in light blue
    ax.plot(edges[:-1], full_hist, '-', color='lightblue', label='All data')
    # Plot selected pattern in black
    ax.plot(edges[:-1], selected_hist, 'k-', label='Selection')

    # Add d-spacing limit markers if this is a box selection
    if selected_experiments is not None:
      ymin, ymax = ax.get_ylim()
      marker_height = ymax * 0.05  # 5% of plot height
      s_min, s_max = 1/d_max_b, 1/d_min_b  # Convert d-spacing limits to s
      ax.vlines([s_min, s_max], ymin, ymin + marker_height, color='red', linewidth=1)

    # Set x-limits
    ax.set_xlim(self.plot_s_min, self.plot_s_max)



    ax.get_xaxis().set_major_formatter(tick.FuncFormatter(
      lambda x, _: "{:.3f}".format(1/x)))
    ax.set_xlabel('d-spacing (Å)')
    ax.set_ylabel('Counts')
    ax.set_title(title)
    ax.legend()
    self.powder_fig.canvas.draw()

  def make_powder_pattern(self, selected_experiments=None, rescale_ymax=None):
    """
    Create histogram of d-spacings from selected experiments

    Args:
      selected_experiments: List of experiment indices to include. If None, use all.
    """
    if selected_experiments is None:
      selected_experiments = range(len(self.experiments))

    # Collect d-spacings from selected experiments
    all_d_spacings = []
    for i_expt in selected_experiments:
      all_d_spacings.extend(self.powder_data[i_expt])

    # Convert to s-magnitudes for binning
    s_mags = 1/np.array(all_d_spacings)

    # Create histogram with equal-width bins in s
    hist, edges = np.histogram(
      s_mags,
      bins=self.params.plot.d_bins,
      range=(self.plot_s_min, self.plot_s_max)
    )

    if rescale_ymax is not None and hist.max() > 0:
      hist = hist * (rescale_ymax / hist.max())

    return hist, edges

  def compute_spot_scatter(self, selected_experiments, range_a, range_b, range_theta):
    if selected_experiments is None:
      return None

    valid_positions = []
    ref_spots_A = []
    ref_spots_B = []

    for i_expt in selected_experiments:
      i_expt = int(i_expt)
      # Get all spots for this experiment
      exp_sel = flex.bool(self.reflections['id'] == i_expt)
      exp_refls = self.reflections.select(exp_sel)

      # Get beam center
      detector = self.experiments[i_expt].detector[0]  # first panel
      beam = self.experiments[i_expt].beam
      beam_x, beam_y = detector.get_beam_centre(beam.get_s0())

      # Get d-spacings and find spots A and B
      s0 = beam.get_s0()
      s1_vectors = exp_refls['s1'].as_numpy_array()
      d_spacings = np.array([1/np.linalg.norm(s1 - s0) for s1 in s1_vectors])

      # Find unique spot A
      d_max_A, d_min_A = range_a
      spots_A = (d_spacings >= d_min_A) & (d_spacings <= d_max_A)
      if np.sum(spots_A) != 1:
        continue

      # Find unique spot B using angles
      spot_A = s1_vectors[spots_A][0]
      angles = []
      for s1 in s1_vectors:
        vec_a = spot_A - s0
        vec_b = s1 - s0
        cos_angle = np.dot(vec_a, vec_b) / (np.linalg.norm(vec_a) * np.linalg.norm(vec_b))
        angles.append(np.rad2deg(np.arccos(min(1.0, max(-1.0, cos_angle)))))
      angles = np.array(angles)

      d_max_B, d_min_B = range_b
      th_min, th_max = range_theta
      spots_B = (d_spacings >= d_min_B) & (d_spacings <= d_max_B) & \
                (angles >= th_min) & (angles <= th_max)
      if np.sum(spots_B) != 1:
        continue

      # Get coordinates of all spots
      xyz = exp_refls['xyzobs.mm.value']
      x = xyz.parts()[0]
      y = xyz.parts()[1]

      # Translate to make beam center the origin
      x = x - beam_x
      y = y - beam_y

      # Get coordinates of spot A (after translation)
      spot_A_x = x.select(flex.bool(spots_A))[0]
      spot_A_y = y.select(flex.bool(spots_A))[0]

      # Calculate rotation angle to put spot A at 12 o'clock
      angle = np.arctan2(spot_A_x, spot_A_y)

      # Apply rotation to all spots
      cos_ang = np.cos(angle)
      sin_ang = np.sin(angle)
      x_rot = x * cos_ang - y * sin_ang
      y_rot = x * sin_ang + y * cos_ang

      # Check if spot B needs reflection
      spot_B_x = x_rot.select(flex.bool(spots_B))[0]
      if spot_B_x < 0:
        x_rot = -x_rot  # Reflect across y-axis

      valid_positions.append(np.column_stack((x_rot, y_rot)))
      ref_spots_A.append(np.array([x_rot.select(flex.bool(spots_A))[0], y_rot.select(flex.bool(spots_A))[0]]))
      ref_spots_B.append(np.array([x_rot.select(flex.bool(spots_B))[0], y_rot.select(flex.bool(spots_B))[0]]))

    if not valid_positions:
      return None

    return {
      'x': np.concatenate([p[:,0] for p in valid_positions]),
      'y': np.concatenate([p[:,1] for p in valid_positions]),
      'ref_A': np.array(ref_spots_A),
      'ref_B': np.array(ref_spots_B)
    }

  def plot_spot_scatter(self, positions):
    """Plot transformed spot positions"""
    if not hasattr(self, 'scatter_fig'):
      self.scatter_fig = plt.figure(figsize=(8,8))
      self.scatter_fig.canvas.manager.window.wm_geometry("+500+50")
      ax = self.scatter_fig.add_axes([0.1, 0.1, 0.85, 0.85])
      self.scatter_fig.show()
    else:
      self.scatter_fig.clear()
      ax = self.scatter_fig.add_axes([0.1, 0.1, 0.85, 0.85])

    if positions is not None:
      ax.scatter(positions['x'], positions['y'], alpha=0.2, s=3)
      ax.scatter(positions['ref_A'][:,0], positions['ref_A'][:,1],
                color='red', s=5, label='Spot A')
      ax.scatter(positions['ref_B'][:,0], positions['ref_B'][:,1],
                color='green', s=5, label='Spot B')

    ax.set_xlim((-75,75))
    ax.set_ylim((-75,75))
    ax.set_aspect('equal')
    ax.set_xlabel('x (mm)')
    ax.set_ylabel('y (mm)')
    ax.set_title('Transformed detector positions')
    self.scatter_fig.canvas.draw()
    self.scatter_fig.canvas.flush_events()


  def calculate(self):
    """Process all images and accumulate angle/d-spacing data for each resolution range"""

    if 's1' not in self.reflections.keys():
      self.reflections.centroid_px_to_mm(self.experiments)
      self.reflections.map_centroids_to_reciprocal_space(self.experiments)
      if self.params.save_processed_reflections:
        assert len(self.params.input.reflections)==1
        self.reflections.as_file(self.params.input.reflections[0].filename)

    for i_expt, expt in enumerate(self.experiments):
      if i_expt%1000==0:
        print(i_expt)
      if self.params.n_max is not None and i_expt > self.params.n_max:
        break
      # Get spots for this image
      img_sel = self.reflections['id'] == i_expt
      img_refls = self.reflections.select(img_sel)

      # Get s0 and compute d-spacings
      s0 = expt.beam.get_s0()
      s1_vectors = img_refls['s1'].as_numpy_array()
      d_spacings = np.array([1/np.linalg.norm(s1 - s0) for s1 in s1_vectors])

      # Store d-spacings
      self.powder_data[i_expt].extend(d_spacings)

      # For each resolution range, find spots and compute angles
      for d_max, d_min in self.resolution_ranges:
        # Find spots in this resolution range
        range_sel = (d_spacings >= d_min) & (d_spacings <= d_max)
        range_spots_idx = np.where(range_sel)[0]

        # For each spot in range
        for idx_a in range_spots_idx:
          s1_a = s1_vectors[idx_a]
          d_a = d_spacings[idx_a]

          # Compare to all other spots in image
          for idx_b in range(len(s1_vectors)):
            if idx_b == idx_a:
              continue

            s1_b = s1_vectors[idx_b]
            d_b = d_spacings[idx_b]

            angle = vector_angle(s1_a-s0, s1_b-s0)

            # Store results for this pair
            self.results[(d_max, d_min)]['d_spacings_A'].append(d_a)
            self.results[(d_max, d_min)]['d_spacings_B'].append(d_b)
            self.results[(d_max, d_min)]['angles'].append(angle)
            self.results[(d_max, d_min)]['i_expt'].append(i_expt)

    # Convert accumulated lists to numpy arrays
    for range_key in self.results:
      for key in ['d_spacings_A', 'd_spacings_B', 'angles']:
        self.results[range_key][key] = np.array(self.results[range_key][key])

  def plot(self):
    """Create heatmaps for each resolution range and initial powder pattern"""

    # Get the plot limits as scattering magnitudes
    s_min = self.plot_s_min
    s_max = self.plot_s_max

    # Create powder pattern figure first
    powder_fig = plt.figure(figsize=(16,3))
    powder_fig.canvas.manager.window.geometry("+50+400")
    ax = powder_fig.add_axes([.1, .15, .8, .75]) # [left, bottom, width, height]
    hist, edges = self.make_powder_pattern()
    plt.plot(edges[:-1], hist, 'k-')  # Plot against s-magnitude edges directly
    ax.set_xlim(s_min, s_max)
    ax.get_xaxis().set_major_formatter(tick.FuncFormatter(
      lambda x, _: "{:.3f}".format(1/x)))
    plt.xlabel('d-spacing (Å)')
    plt.ylabel('Counts')
    plt.title('Powder Pattern (all data)')

    # Store reference to powder figure for later updates
    self.powder_fig = powder_fig

    # Original heatmap plotting
    for (d_max, d_min), data in self.results.items():
      # Convert d-spacings to s-magnitudes for binning
      s_mags = 1/data['d_spacings_B']

      # Create full histogram
      hist, xedges, yedges = np.histogram2d(
        s_mags,
        data['angles'],
        bins=[self.params.plot.d_bins, self.params.plot.th_bins],
        range=[[s_min, s_max],
               [self.params.plot.th_min, self.params.plot.th_max]]
      )

      # Create masked version for color scaling
      s_centers = (xedges[:-1] + xedges[1:]) / 2
      d_centers = 1/s_centers
      mask = (d_centers >= d_min) & (d_centers <= d_max)
      mask_2d = mask[:, np.newaxis] * np.ones((1, hist.shape[1]))
      hist_masked = np.ma.array(hist, mask=mask_2d)

      fig = plt.figure(figsize=(16,3))
      fig.canvas.manager.window.geometry("+50+50")
      ax = fig.add_axes([.1, .15, .8, .75]) # [left, bottom, width, height]

      cmap = plt.cm.viridis.copy()
      ax.set_facecolor(cmap(0.0))

      if self.params.plot.zscale == 'log':
        norm = LogNorm(vmin=0.1, vmax=hist_masked.max())
      else:  # linear
        norm = Normalize(vmin=0, vmax=hist_masked.max())

      plt.imshow(hist.T,
                norm=norm,
                cmap=cmap,
                origin='lower',
                aspect='auto',
                extent=[s_min, s_max,
                       self.params.plot.th_min, self.params.plot.th_max])

      ax.get_xaxis().set_major_formatter(tick.FuncFormatter(
        lambda x, _: "{:.3f}".format(1/x)))

      # After creating the heatmap...
      selector = RectangleSelector(
        ax,
        lambda eclick, erelease, ax=ax: self.on_select(eclick, erelease, (d_max, d_min), ax),
        useblit=True,
        props=dict(facecolor='white', edgecolor='black', alpha=0.5, fill=True)
      )
      # Store selector to prevent garbage collection
      self._selectors = getattr(self, '_selectors', [])
      self._selectors.append(selector)

      #plt.colorbar(label='Counts')
      plt.xlabel('d-spacing (Å)')
      plt.ylabel('Angle (degrees)')
      plt.title(f'Angle vs reference peak ({d_max}-{d_min}Å)')


    plt.show()

		

class Script(object):
  def __init__(self):
    usage = "experimental"
    self.parser = ArgumentParser(
        usage=usage,
        phil=phil_scope,
        epilog=help_message,
        check_format=False,
        read_reflections=True,
        read_experiments=True,
    )


  def run(self):
    params, options = self.parser.parse_args(show_diff_phil=False)
    assert len(params.input.experiments) == len(params.input.reflections) == 1
    experiments = params.input.experiments[0].data
    reflections = params.input.reflections[0].data

    ahm = Angle_heatmap(experiments, reflections, params)
    ahm.calculate()
    ahm.plot()


if __name__ == "__main__":
  with show_mail_on_error():
    script = Script()
    script.run()
