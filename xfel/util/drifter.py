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
  uncertainties = False
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


DEFAULT_INPUT_SCOPE = parse("""
  input {
    path = None
      .multiple = True
      .type = str
    experiments_suffix = .expt
      .type = str
    reflections_suffix = .refl
      .type = str
    experiments = None
      .multiple = True
      .type = str
    reflections = None
      .multiple = True
      .type = str
}
""")


def get_scaling_directories(merging_phil_paths):
  """Return paths to all directories specified as input.path in phil file"""
  merging_phils = [parse(file_name=mpp) for mpp in merging_phil_paths]
  phil = DEFAULT_INPUT_SCOPE.fetch(sources=merging_phils).extract()
  return sorted(set(phil.input.path))


def get_tder_refined_expts_and_refls(scaling_phil_paths):
  """Return paths to all expt and refl files mentioned in phil files"""
  expt_paths = []
  refl_paths = []
  scaling_phils = [parse(file_name=spp) for spp in scaling_phil_paths]
  phil = DEFAULT_INPUT_SCOPE.fetch(sources=scaling_phils).extract()
  refined_tder_input_paths = [ip.replace('*reintegrated*', '*refined*')
                              for ip in phil.input.path]
  for rip in refined_tder_input_paths:
    expt_glob = os.path.join(rip + phil.input.experiments_suffix)
    refl_glob = os.path.join(rip + phil.input.reflections_suffix)
    expt_paths.extend(glob.glob(expt_glob))
    refl_paths.extend(glob.glob(refl_glob))
  return sorted(set(expt_paths)), sorted(set(refl_paths))


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


def load_experiments(*expt_paths):
  """Create an instance of ExperimentList from *expt_paths"""
  expts = ExperimentList()
  for expt_path in expt_paths:
    expts.extend(ExperimentList.from_file(expt_path, check_format=False))
  return expts


def load_reflections(*refl_paths):
  """Create an instance of flex.reflection_table from *refl_paths"""
  refls_list = [flex.reflection_table.from_file(rp) for rp in refl_paths]
  return flex.reflection_table.concat(refls_list)


def get_extra_info_dict_from_refinement_expts_and_refls(
        tder_expt_paths, tder_refl_paths, uncertainties=True):
  """Retrieve total number of spot-finding reflections as well as uncertainties
  of origin positions from refined TDER reflection position deviations"""
  tder_expts = load_experiments(*tder_expt_paths)
  tder_refls = load_reflections(*tder_refl_paths)
  return_dict = {'expts': len(tder_expts), 'refls': len(tder_refls)}
  deltas_flex = flex.vec3_double()
  if uncertainties:
    for panel in tder_expts[0].detector:               # pr:  panel reflections
      pr = tder_refls.select(tder_refls['panel'] == panel.index())
      pr_obs_det = pr['xyzobs.mm.value'].parts()[0:2]  # det: in detector space
      pr_cal_det = pr['xyzcal.mm'].parts()[0:2]        # lab: in lab space
      pr_obs_lab = panel.get_lab_coord(flex.vec2_double(*pr_obs_det))
      pr_cal_lab = panel.get_lab_coord(flex.vec2_double(*pr_cal_det))
      deltas_flex.extend(pr_obs_lab - pr_cal_lab)
    d = [flex.mean(flex.abs(deltas_flex.parts()[i])) for i in range(3)]
    return_dict.update({'delta_x': d[0], 'delta_y': d[1], 'delta_z': d[2]})
  return return_dict


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
  EXTRA_KEYS = ['delta_x', 'delta_y', 'delta_z', 'expts', 'refls']

  def __init__(self, parameters):
    self.data = {k: [] for k in self.REQUIRED_KEYS + self.EXTRA_KEYS}
    self.parameters = parameters

  def __len__(self):
    self.assert_all_active_keys_have_same_length()
    return len(self.data[self.REQUIRED_KEYS[0]])

  def __str__(self):
    lines = [' '.join('{!s:9.9}'.format(k.upper()) for k in self.active_keys)]
    for row in self.rows:
      cells = ['{:.20f}'.format(c) if isinstance(c, float) else c for c in row]
      lines.append(' '.join('{!s:9.9}'.format(cell) for cell in cells))
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
      for scaling_dir in get_scaling_directories(merging_phils):
        scaling_expts = glob.glob(os.path.join(scaling_dir, 'scaling_*.expt'))
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
        if self.parameters.uncertainties:
          scaling_phils = glob.glob(path_join(task_dir, 'params_1.phil'))
          tder_expts, tder_refls = get_tder_refined_expts_and_refls(scaling_phils)
          extra_info_dict = get_extra_info_dict_from_refinement_expts_and_refls(
            tder_expts, tder_refls, uncertainties=self.parameters.uncertainties)
          self.add(**extra_info_dict)


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
    gs = GridSpec(4, 2, hspace=0, wspace=0, height_ratios=[1, 2, 2, 2],
                  width_ratios=[4, 1])
    self.axh = self.fig.add_subplot(gs[0, 0])
    self.axx = self.fig.add_subplot(gs[1, 0], sharex=self.axh)
    self.axy = self.fig.add_subplot(gs[2, 0], sharex=self.axh)
    self.axz = self.fig.add_subplot(gs[3, 0], sharex=self.axh)
    self.axw = self.fig.add_subplot(gs[0, 1], sharey=self.axh)
    self.axl = self.fig.add_subplot(gs[1:, 1])

  def _setup_figure(self):
    self.axh.spines['top'].set_visible(False)
    self.axh.spines['right'].set_visible(False)
    self.axh.spines['bottom'].set_visible(False)
    self.axh.set_zorder(self.axx.get_zorder() - 1.0)
    self.axw.axis('off')
    common = {'direction': 'inout', 'top': True, 'bottom': True, 'length': 6}
    self.axx.tick_params(axis='x', labelbottom=False, **common)
    self.axy.tick_params(axis='x', labelbottom=False, **common)
    self.axz.tick_params(axis='x', rotation=90, **common)
    self.axx.ticklabel_format(useOffset=False)
    self.axy.ticklabel_format(useOffset=False)
    self.axz.ticklabel_format(useOffset=False)
    self.axh.set_ylabel('# expts')
    self.axx.set_ylabel('Detector X')
    self.axy.set_ylabel('Detector Y')
    self.axz.set_ylabel('Detector Z')
    self.axz.set_xlabel(self.order_by.title())

  @property
  def _refl_to_expt_ratios(self):
    expts, refls = self.registry.data['expts'], self.registry.data['refls']
    return [r / e for e, r in zip(expts, refls)]

  @property
  def color_array(self):
    """Registry-length color list with colors corresponding to self.color_by"""
    color_i = [0, ] * len(self.registry)
    color_by = self.registry.data[self.color_by]
    for i, cat in enumerate(color_by[1:], 1):
      color_i[i] = color_i[i-1] if cat in color_by[:i] else color_i[i-1] + 1
    return [self.colormap(i % self.colormap_period) for i in color_i]

  def _get_handles_and_labels(self):
    unique_keys = []
    for key in self.registry.data[self.color_by]:
      if key not in unique_keys:
        unique_keys.append(key)
    handles = [Line2D([], [], c=self.colormap(i % 10), ls='', ms=12, marker='.')
               for i in range(len(unique_keys))]
    return handles, unique_keys

  def _plot_drift(self, axes, values_key, deltas_key=None, top=False):
    x = self.registry.data[self.order_by]
    y = self.registry.data[values_key]
    y_err = self.registry.data.get(deltas_key, [])
    axes.scatter(x, y, c=self.color_array)
    if top:
      ax_t = axes.secondary_xaxis('top')
      ax_t.tick_params(rotation=90)
      ax_t.set_xticks(axes.get_xticks())
      ax_t.set_xticklabels(self.registry.data['expts'])
    if self.parameters.uncertainties:
      axes.errorbar(x, y, yerr=y_err, ecolor='black', ls='')

  def _plot_bars(self):
    x = self.registry.data[self.order_by]
    y = self.registry.data['expts']
    ratios = self._refl_to_expt_ratios
    self.axh.bar(x, y, width=[ratio / max(ratios) for ratio in ratios],
                 color=self.color_array, alpha=0.5)

  def _plot_legend(self):
    handles, labels = self._get_handles_and_labels()
    self.axl.legend(handles, labels, loc=7)
    self.axl.axis('off')

  def _plot_width_info(self):
    ex = [min(self.registry.data['expts']), max(self.registry.data['expts']),
          min(self._refl_to_expt_ratios), max(self._refl_to_expt_ratios),
          min(self.registry.data['refls']), max(self.registry.data['refls'])]
    s = "bar heights: # of expts ({} to {})\n" \
        "bar widths: # refl / expt ({:.2f} to {:.2f})\n" \
        "bar areas: # of refls ({} to {})".format(*ex)
    self.axw.annotate(text=s, xy=(0.5, 0.5), xycoords='axes fraction',
                      ha='center', ma='center', va='center')

  def plot(self):
    self._plot_bars()
    self._plot_width_info()
    self._plot_drift(self.axx, 'x', 'delta_x', top=True)
    self._plot_drift(self.axy, 'y', 'delta_y')
    self._plot_drift(self.axz, 'z', 'delta_z')
    self._plot_legend()
    self.fig.align_labels()


def run(params_):
  ddr = DetectorDriftRegistry(parameters=params_)
  ddr.populate_from_merging_job(tag_pattern=params_.input_glob)
  ddr.sort(by_key='run')
  print(ddr)
  dda = DetectorDriftArtist(registry=ddr, parameters=params_)
  dda.plot()
  if params_.pickle_plot:
    try:
      from libtbx.easy_pickle import dump
      dump(str(params_.pickle_filename), dda.fig)
    except TypeError as e:
      print(e)
  if params_.show_plot:
    plt.show()


if __name__ == '__main__':
  if '--help' in sys.argv[1:] or '-h' in sys.argv[1:]:
    print(message)
    exit()
  params = params_from_phil(sys.argv[1:])
run(params)
