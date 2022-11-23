from __future__ import absolute_import, print_function, division
import glob
import sys
import os
from dxtbx.model.experiment_list import ExperimentList
from libtbx.phil import parse
from scitbx.array_family import flex
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from .weather import params_from_phil

message = ''' This script aims to investigate the spatial drift of a detector
              as a function of experimental progress. It requires the directory
              structure to follow that resulting from indexing and ensemble
              refinement performed by cctbx.xfel. End result is a multiplot
              with detector origin position (vertical position) as a function
              of run number (horizontal position), colored according to tag
              name. Error bars are derived from the uncertainty of individual
              reflections' position in laboratory reference system.
              Example usage: `libtbx.python drifter.py input_glob=batch*TDER/`
'''
phil_scope = parse('''
  input_glob = batch*TDER/
    .type = str
    .help = glob which matches all directories after TDER to be investigated.
  pickle_plot = False
    .type = bool
    .help = If True, will pickle matplotlib session so that it can be opened later for analysis/viewing \
            https://stackoverflow.com/questions/29160177/matplotlib-save-file-to-be-reedited-later-in-ipython
  pickle_filename = fig_object.pickle
    .type = str
    .help = Default name of pickled matplotlib plot saved to disk
  show_plot = True
    .type = bool
    .help = If True, resulting plot will be displayed on screen
  uncertainties = True
    .type = bool
    .help = If True, origin uncertainties will be estimated using differences
            between predicted and observed positions of all reflections.
''')


class DriftRegistry(object):
  """Object responsible for storing information about all detector positions"""
  REQUIRED_KEYS = ['tag', 'run', 'rungroup', 'trial', 'task', 'x', 'y', 'z']
  EXTRA_KEYS = ['delta_x', 'delta_y', 'delta_z']

  def __init__(self):
    self.data = {k: [] for k in self.REQUIRED_KEYS + self.EXTRA_KEYS}

  def __len__(self):
    self.assert_all_active_keys_have_same_length()
    return len(self.data[self.REQUIRED_KEYS[0]])

  def __str__(self):
    lines = ['\t'.join(k.upper() for k in self.active_keys)]
    for row in self.rows:
      lines.append('\t'.join(row))
    return '\n'.join(lines)

  def add(self, **kwargs):
    for k, v in kwargs.items():
      if k not in self.data.keys():
        raise KeyError('Attempt to add unknown key to DriftRegistry: ' + k)
      self.data[k].append(v)

  def sort(self, by_key=REQUIRED_KEYS[0]):
    reference = self.data[by_key]
    for key in self.active_keys:
      sortee = self.data[key]
      self.data[key] = [x for _, x in sorted(zip(reference, sortee),
                                             key=lambda pair: pair[0])]

  def assert_all_active_keys_have_same_length(self):
    assert len(list({len(self.data[k]) for k in self.active_keys})) == 1

  @property
  def active_keys(self):
    return self.REQUIRED_KEYS + [k for k in self.REQUIRED_KEYS if self.data[k]]

  @property
  def rows(self):
    columns = [self.data.get(k) for k in self.active_keys]
    return zip(*[column for column in columns])


class DriftArtist(object):
  """Object responsible for plotting an instance of `DriftRegistry`."""
  def __init__(self, registry=DriftRegistry()):
    self.colormap = plt.get_cmap("tab10")
    self.colormap_period = 10
    self.color_by = 'tag'
    self.order_by = 'run'
    self.registry = registry
    self._init_figure()
    self._setup_figure()

  def _init_figure(self):
    self.fig = plt.figure(tight_layout=True)
    gs = GridSpec(3, 5)
    self.axx = self.fig.add_subplot(gs[0, :4])
    self.axy = self.fig.add_subplot(gs[1, :4], sharex=self.axx)
    self.axz = self.fig.add_subplot(gs[2, :4], sharex=self.axx)
    self.axl = self.fig.add_subplot(gs[:, 4])

  def _setup_figure(self):
    self.axx.tick_params(axis='x', which='both', labelbottom=False)
    self.axy.tick_params(axis='x', which='both', labelbottom=False)
    self.axz.tick_params(axis='x', which='both', rotation=90)
    self.axx.set_ylabel('Detector X')
    self.axy.set_ylabel('Detector Y')
    self.axz.set_ylabel('Detector Z')
    self.axz.set_xlabel(self.order_by.title())

  @property
  def color_array(self):
    """Registry-length color list with colors corresponding to self.color_by"""
    color_i = [0, ] * len(self.registry)
    color_by = self.registry.data[self.color_by]
    for i, cat in enumerate(color_by[1:], 1):
      color_i[i] = color_i[i-1] if cat in color_by[:i] else color_i[i-1] + 1
    return [self.colormap(i % self.colormap_period) for i in color_i]

  def _plot_axes(self, axes, values_key, deltas_key=None):
    x = self.registry.data[self.order_by]
    y = self.registry.data[values_key]
    y_err = self.registry.data.get(deltas_key, [])
    try:
      axes.errorbar(x, y, yerr=y_err, ecolor='gray')
    except IndexError:
      pass
    axes.scatter(x, y, c=self.color_array)

  def _plot_legend(self):
    unique_keys = []
    for key in self.registry.data[self.color_by]:
      if key not in unique_keys:
        unique_keys.append(key)
    handles = [Line2D([], [], c=self.colormap(i % 10), ls='', marker='.')
               for i in range(len(unique_keys))]
    self.axl.legend(handles, unique_keys, loc=7)
    self.axl.axis('off')

  def plot(self):
    self._plot_axes(self.axx, 'x', 'delta_x')
    self._plot_axes(self.axy, 'y', 'delta_y')
    self._plot_axes(self.axz, 'z', 'delta_z')
    self._plot_legend()


def extract_origin_from_expt_file(expt_path):
  """Read origin from lines 14-16 after 'hierarchy' text in expt file"""
  origin = []
  line_counter = 0
  with open(expt_path, 'r') as expt_file:
    for line in expt_file:
      if 'hierarchy' in line:
        line_counter = 1
      elif line_counter:
        line_counter += 1
      if 15 <= line_counter <= 17:
        origin.append(float(line.replace(',', '').strip()))
      elif line_counter > 17:
        break
  return origin


def extract_deltas_from_expt_and_refl_files(expt_path, refl_path):
  """Estimate origin deltas from individual reflection position deviations"""
  reflections = flex.reflection_table.from_file(refl_path)
  experiments = ExperimentList.from_file(expt_path, check_format=False)
  detector = experiments[0].detector
  for panel_index, panel in enumerate(detector):
    pr = reflections.select(reflections['panel'] == 0)  # panel reflections
    pr_obs = panel.get_lab_coord(flex.vec2_double(*pr['xyzobs.mm.value'].parts()[0:2]))
    pr_cal = panel.get_lab_coord(flex.vec2_double(*pr['xyzcal.mm'].parts()[0:2]))
  # get deltas for all 128 panels; calculate std dev for a list of all deltas
  # code below will fail; only rough idea
  # pr_deltas2 = abs(pr_obs**2 - pr_cal**2)
  # stddev = (sum(pr_deltas2) / (len(pr_deltas2) - 1)) ** 0.5
  return 0.1, 0.2, 0.3


def run(params):
  dr = DriftRegistry()
  da = DriftArtist(registry=dr)
  tag_pattern = 'batch*_TDER/'
  tag_list = glob.glob(tag_pattern)
  for tag in tag_list:
    version_list = glob.glob(os.path.join(tag, 'v*/'))
    version_last = sorted(version_list)[-1]
    phil_list = glob.glob(os.path.join(version_last, '*_params.phil'))
    input_list = []
    for phil_path in phil_list:
      with open(phil_path, "r") as phil_file:
        for line in phil_file:
          if line[:11] == 'input.path=':
            input_list.append(line[11:].strip())
    input_list = list(set(input_list))
    for input_path in input_list:
      expt_list = glob.glob(os.path.join(input_path, 'scaling_*.expt'))
      try:
        expt_first = sorted(expt_list)[0]
      except IndexError:
        continue
      task_dir = os.path.join(*os.path.normpath(expt_first).split(os.sep)[:-2])
      task = int(os.path.split(task_dir)[-1][-3:])
      trial_dir = os.path.join(*os.path.normpath(task_dir).split(os.sep)[:-1])
      trial = int(trial_dir[-9:-6])
      rungroup = int(trial_dir[-3:])
      run_name = os.path.normpath(trial_dir).split(os.sep)[-2]
      origin = extract_origin_from_expt_file(expt_first)
      dr.add(tag=tag, run=run_name, rungroup=rungroup, trial=trial, task=task,
             x=origin[0], y=origin[1], z=origin[2])

  dr.sort(by_key='run')
  print(dr)
  da.plot()
  if params.pickle_plot:
    from libtbx.easy_pickle import dump
    dump('%s' % params.pickle_filename, da.fig)
  if params.show_plot:
    plt.show()


if __name__ == '__main__':
  if '--help' in sys.argv[1:] or '-h' in sys.argv[1:]:
    print(message)
    exit()
  params = params_from_phil(sys.argv[1:])
run(params)
