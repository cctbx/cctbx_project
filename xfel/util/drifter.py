from __future__ import absolute_import, print_function, division
import glob
import sys
import os
from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentList
from libtbx.phil import parse
from libtbx.utils import Sorry
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

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
    .help = If True, matplotlib session will be pickled so that it can be \
            opened later for viewing: https://www.stackoverflow.com/q/29160177/
  pickle_filename = fig_object.pickle
    .type = str
    .help = Default name of pickled matplotlib plot saved to disk
  show_plot = True
    .type = bool
    .help = If True, resulting plot will be displayed on screen
  sizes = True
    .type = bool
    .help = If True, the size of points will depend on number of experiments
  uncertainties = True
    .type = bool
    .help = If True, origin uncertainties will be estimated using differences \
            between predicted and observed positions of all reflections.
''')


def params_from_phil(args):
  user_phil = []
  for arg in args:
    if os.path.isfile(arg):
      user_phil.append(parse(file_name=arg))
    else:
      try:
        user_phil.append(parse(arg))
      except Exception:
        raise Sorry("Unrecognized argument: %s" % arg)
  params_ = phil_scope.fetch(sources=user_phil).extract()
  return params_


def get_input_paths_from_phils(phil_paths):
  """Return reminder of lines which start with 'input.path=' in all phils"""
  input_list = []
  for phil_path in phil_paths:
    with open(phil_path, "r") as phil_file:
      for line in phil_file:
        if line.startswith('input.path='):
          input_list.append(line[11:].strip())
  return sorted(list(set(input_list)))


def get_origin_from_expt(expt_path):
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


def get_deltas_from_expts_and_refls(expt_paths, refl_paths):
  """Estimate origin deltas from individual reflection position deviations"""
  deltas = flex.vec3_double()
  for expt_path, refl_path in zip(expt_paths, refl_paths):
    reflections = flex.reflection_table.from_file(refl_path)
    experiments = ExperimentList.from_file(expt_path, check_format=False)
    detector = experiments[0].detector                 # pr:  panel reflections
    for panel in detector:                             # observed / calculated
      pr = reflections.select(reflections['panel'] == panel.index())
      pr_obs_det = pr['xyzobs.mm.value'].parts()[0:2]  # det: in detector space
      pr_cal_det = pr['xyzcal.mm'].parts()[0:2]        # lab: in lab space
      pr_obs_lab = panel.get_lab_coord(flex.vec2_double(*pr_obs_det))
      pr_cal_lab = panel.get_lab_coord(flex.vec2_double(*pr_cal_det))
      deltas.extend(pr_obs_lab - pr_cal_lab)
  return [deltas.parts()[i].sample_standard_deviation() for i in (0, 1, 2)]


def get_size_from_refls(refl_paths):
  """Retrieve a total number of indexed reflections from list of refl files"""
  return sum(len(flex.reflection_table.from_file(r)) for r in refl_paths)


def path_join(*path_elements):
  """Join path from elements, resolving all redundant or relative calls"""
  path_elements = [os.pardir if p == '..' else p for p in path_elements]
  return os.path.normpath(os.path.join(*path_elements))


def path_split(path):
  """Split path into directories and file basename using os separator"""
  return os.path.normpath(path).split(os.sep)


class DetectorDriftRegistry(object):
  """Object responsible for storing information about all detector positions"""
  REQUIRED_KEYS = ['tag', 'run', 'rungroup', 'trial', 'task', 'x', 'y', 'z']
  EXTRA_KEYS = ['delta_x', 'delta_y', 'delta_z', 'size']

  def __init__(self, parameters):
    self.data = {k: [] for k in self.REQUIRED_KEYS + self.EXTRA_KEYS}
    self.parameters = parameters

  def __len__(self):
    self.assert_all_active_keys_have_same_length()
    return len(self.data[self.REQUIRED_KEYS[0]])

  def __str__(self):
    lines = [' '.join('{!s:9.9}'.format(k.upper()) for k in self.active_keys)]
    for row in self.rows:
      lines.append(' '.join('{!s:9.9}'.format(cell) for cell in row))
    return '\n'.join(lines)

  def add(self, **kwargs):
    for k, v in kwargs.items():
      if k not in self.data.keys():
        raise KeyError('Attempt to add unknown key to DriftRegistry: ' + k)
      self.data[k].append(v)

  def sort(self, by_key=REQUIRED_KEYS[0]):
    self.assert_all_active_keys_have_same_length()
    reference = self.data[by_key]
    for key in self.active_keys:
      sortee = self.data[key]
      self.data[key] = [x for _, x in sorted(zip(reference, sortee),
                                             key=lambda pair: pair[0])]

  def assert_all_active_keys_have_same_length(self):
    if len({len(self.data[k]) for k in self.active_keys}) > 1:
      raise IndexError('Registry items have different lengths: ' +
                       str({k, len(self.data[k])} for k in self.active_keys))

  @property
  def active_keys(self):
    return self.REQUIRED_KEYS + [k for k in self.EXTRA_KEYS if self.data[k]]

  @property
  def rows(self):
    self.assert_all_active_keys_have_same_length()
    columns = [self.data[k] for k in self.active_keys]
    return zip(*[column for column in columns])

  def populate_from_merging_job(self, tag_pattern='*/'):
    tag_list = glob.glob(tag_pattern)
    for tag in tag_list:
      version_last = sorted(glob.glob(path_join(tag, 'v*/')))[-1]
      merging_phils = glob.glob(path_join(version_last, '*_params.phil'))
      merging_input_list = get_input_paths_from_phils(merging_phils)
      for input_path in merging_input_list:
        scaling_expts = glob.glob(os.path.join(input_path, 'scaling_*.expt'))
        try:
          first_scaling_expt_path = sorted(scaling_expts)[0]
        except IndexError:
          continue
        task_dir = path_join(first_scaling_expt_path, '..', '..')
        trial_dir = path_join(task_dir, '..')
        task = int(path_split(task_dir)[-1][-3:])
        trial = int(trial_dir[-9:-6])
        rungroup = int(trial_dir[-3:])
        run_name = path_split(trial_dir)[-2]
        print('Processing run {} in tag {}'.format(run_name, tag))
        origin = get_origin_from_expt(first_scaling_expt_path)
        self.add(tag=tag, run=run_name, rungroup=rungroup, trial=trial,
                 task=task, x=origin[0], y=origin[1], z=origin[2])
        if self.parameters.uncertainties or self.parameters.sizes:
          scaling_phils = glob.glob(path_join(task_dir, 'params_1.phil'))
          scaling_input_list = get_input_paths_from_phils(scaling_phils[-1:])
          tder_expts, tder_refls = [], []
          for input_path2 in scaling_input_list:
            tder_expts.extend(sorted(glob.glob(input_path2 + '.expt')))
            tder_refls.extend(sorted(glob.glob(input_path2 + '.refl')))
          if self.parameters.sizes:
            self.add(size=get_size_from_refls(tder_refls))
          if self.parameters.uncertainties:
            deltas = get_deltas_from_expts_and_refls(tder_expts, tder_refls)
            self.add(delta_x=deltas[0], delta_y=deltas[1], delta_z=deltas[2])


class DetectorDriftArtist(object):
  """Object responsible for plotting an instance of `DriftRegistry`."""
  def __init__(self, registry, parameters):
    self.colormap = plt.get_cmap("tab10")
    self.colormap_period = 10
    self.color_by = 'tag'
    self.order_by = 'run'
    self.registry = registry
    self.parameters = parameters
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

  @property
  def size_array(self):
    """Registry-length size list with sizes corr. to number of experiments"""
    return [0.0001 * s for s in self.registry.data['size']] \
      if self.parameters.sizes else [36] * len(self.registry)

  def _get_all_handles_and_labels(self):
    handles1, labels1 = self._get_standard_handles_and_labels()
    handles2, labels2 = self._get_size_handles_and_labels()
    return handles1 + handles2, labels1 + labels2

  def _get_standard_handles_and_labels(self):
    unique_keys = []
    for key in self.registry.data[self.color_by]:
      if key not in unique_keys:
        unique_keys.append(key)
    handles = [Line2D([], [], c=self.colormap(i % 10), ls='', marker='.', ms=36)
               for i in range(len(unique_keys))]
    return handles, unique_keys

  def _get_size_handles_and_labels(self):
    handles, labels = [], []
    if self.parameters.sizes:
      handles.append(Line2D([0], [0], alpha=0))
      labels.append('')
      size_min, size_max = min(self.size_array), max(self.size_array)
      for n in range(10):
        if size_min / 3.2 <= 10**n <= size_max * 3.2:
          handles.append(Line2D([], [], c='black', ls='', ms=(10**n)**0.5, marker='.'))
          labels.append('{} reflections'.format(10**n))
    return handles, labels

  def _plot_axes(self, axes, values_key, deltas_key=None):
    x = self.registry.data[self.order_by]
    y = self.registry.data[values_key]
    y_err = self.registry.data.get(deltas_key, [])
    axes.scatter(x, y, c=self.color_array, s=self.size_array)
    if self.parameters.uncertainties:
      axes.errorbar(x, y, yerr=y_err, ecolor='black', ls='')

  def _plot_legend(self):
    handles, labels = self._get_all_handles_and_labels()
    self.axl.legend(handles, labels, loc=7)
    self.axl.axis('off')

  def plot(self):
    self._plot_axes(self.axx, 'x', 'delta_x')
    self._plot_axes(self.axy, 'y', 'delta_y')
    self._plot_axes(self.axz, 'z', 'delta_z')
    self._plot_legend()


def run(params_):
  ddr = DetectorDriftRegistry(parameters=params_)
  ddr.populate_from_merging_job(tag_pattern=params_.input_glob)
  ddr.sort(by_key='run')
  print(ddr)
  dda = DetectorDriftArtist(registry=ddr, parameters=params_)
  dda.plot()
  if params_.pickle_plot:
    from libtbx.easy_pickle import dump
    dump(str(params_.pickle_filename), dda.fig)
  if params_.show_plot:
    plt.show()


if __name__ == '__main__':
  if '--help' in sys.argv[1:] or '-h' in sys.argv[1:]:
    print(message)
    exit()
  params = params_from_phil(sys.argv[1:])
run(params)
