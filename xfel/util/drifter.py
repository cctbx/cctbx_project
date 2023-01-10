from __future__ import absolute_import, print_function, division
import glob
import os
import sys
import tempfile
from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentList
from libtbx.phil import parse
from libtbx.utils import Sorry
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import PercentFormatter


message = ''' This script aims to investigate the spatial drift of a detector
              as a function of experimental progress. It requires the directory
              structure to follow that resulting from indexing and ensemble
              refinement performed by cctbx.xfel. End result is a multiplot
              with detector origin position (vertical position) as a function
              of run number (horizontal position), colored according to tag
              name. Error bars are derived from the uncertainty of individual
              reflections' position in laboratory reference system.
              Example usage: `libtbx.python `libtbx.find_in_repositories
              xfel`/util/drifter.py input.glob=batch*TDER/`, where `batch*TDER`
              describes folders (and thus dataset names) with merging results.
'''
phil_scope = parse('''
  input {
    glob = None
      .type = str
      .multiple = True
      .help = glob which matches directories after TDER to be investigated.
  }
  plot {
    show = True
      .type = bool
      .help = If False, do not display resulting plot interactively
    path = ""
      .type = str
      .help = if given, plot will be saved with this path and name (eg.: fig.png)
    height = 8.0
      .type = float
     .help = Height of saved plot in inches
    width = 10.0
      .type = float
      .help = Width of saved plot in inches
  }
  uncertainties = True
    .type = bool
    .help = If True, uncertainties will be estimated using differences \
          between predicted and observed refl positions and cell distribution
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


def average(sequence, weights=None):
  weights = [1] * len(sequence) if weights is None else weights
  return sum(s * w for s, w in zip(sequence, weights)) / sum(weights)


def correlation(xs, ys, weights=None):
  weights = [1] * len(xs) if weights is None else weights
  x_variance = variance(xs, xs, weights=weights)
  y_variance = variance(ys, ys, weights=weights)
  covariance = variance(xs, ys, weights=weights)
  return covariance / (x_variance * y_variance) ** 0.5


def variance(xs, ys, weights=None):
  weights = [1] * len(xs) if weights is None else weights
  x_avg = average(xs, weights=weights)
  y_avg = average(ys, weights=weights)
  x_dev = [x - x_avg for x in xs]
  y_dev = [y - y_avg for y in ys]
  return sum(w * x * y for w, x, y in zip(weights, x_dev, y_dev)) / sum(weights)


def normalise(sequence):
  min_, max_ = min(sequence), max(sequence)
  return sequence if min_ == max_ else \
    [(s - min_) / (max_ - min_) for s in sequence]


class CorrelationMatrix(object):
  def __init__(self, variables, weights=None):
    """Calculate corr. matrix between every pair in `variables` dict-like"""
    self.keys = variables.keys()
    self.corr = {k: {} for k in self.keys}
    for i1, k1 in enumerate(self.keys):
      for i2, k2 in enumerate(self.keys):
        if i1 == i2:
          self.corr[k1][k2] = self.corr[k2][k1] = 1.0
        elif i2 > i1:
          corr = correlation(variables[k1], variables[k2], weights=weights)
          self.corr[k1][k2] = self.corr[k2][k1] = corr

  def __str__(self):
    s = 'Correlation ' + ' '.join('{:>11}'.format(k) for k in self.keys)
    for k1 in self.keys:
      s += '\n + {:>11} '.format(k1)
      s += ' '.join('    {:+6.4f}'.format(self.corr[k1][k2] for k2 in self.keys))
    return s


class DriftScraper(object):
  """Class for scraping cctbx.xfel output into instance of `DriftTable`"""
  def __init__(self, table, parameters):
    self.table = table
    self.parameters = parameters

  @staticmethod
  def load_experiments(*expt_paths):
    """Create an instance of ExperimentList from *expt_paths"""
    expts = ExperimentList()
    for expt_path in expt_paths:
      expts.extend(ExperimentList.from_file(expt_path, check_format=False))
    return expts

  @staticmethod
  def load_reflections(*refl_paths):
    """Create an instance of flex.reflection_table from *refl_paths"""
    refls_list = [flex.reflection_table.from_file(rp) for rp in refl_paths]
    return flex.reflection_table.concat(refls_list)

  @staticmethod
  def path_join(*path_elements):
    """Join path from elements, resolving all redundant or relative calls"""
    path_elements = [os.pardir if p == '..' else p for p in path_elements]
    return os.path.normpath(os.path.join(*path_elements))

  def path_lookup(self, *path_elements):
    """Join path elements and return a list of all matching files/dirs"""
    return glob.glob(self.path_join(*path_elements))

  @staticmethod
  def path_split(path):
    """Split path into directories and file basename using os separator"""
    return os.path.normpath(path).split(os.sep)

  @staticmethod
  def extract_origin(expt_path):
    """Read origin from lines 14-16 after 'hierarchy' text in expt file"""
    with open(expt_path, 'r') as expt_file:
      expt_lines = expt_file.read().splitlines()
      hierarchy_word_pos = [i for i, li in enumerate(expt_lines)
                            if 'hierarchy' in li][0]
      x = float(expt_lines[hierarchy_word_pos + 14].replace(',', '').strip())
      y = float(expt_lines[hierarchy_word_pos + 15].replace(',', '').strip())
      z = float(expt_lines[hierarchy_word_pos + 16].replace(',', '').strip())
    return {'x': x, 'y': y, 'z': z}

  def extract_size_and_origin_deltas(self, expt_paths, refl_paths):
    """Get number of experiments and reflections, as well as uncertainties
    of origin positions from refined TDER reflection position deviations"""
    tder_expts = self.load_experiments(*expt_paths)
    tder_refls = self.load_reflections(*refl_paths)
    return_dict = {'expts': len(tder_expts), 'refls': len(tder_refls)}
    deltas_flex = flex.vec3_double()
    if self.parameters.uncertainties:
      for panel in tder_expts[0].detector:               # pr:  panel refls
        pr = tder_refls.select(tder_refls['panel'] == panel.index())
        pr_obs_det = pr['xyzobs.mm.value'].parts()[0:2]  # det: in det space
        pr_cal_det = pr['xyzcal.mm'].parts()[0:2]        # lab: in lab space
        pr_obs_lab = panel.get_lab_coord(flex.vec2_double(*pr_obs_det))
        pr_cal_lab = panel.get_lab_coord(flex.vec2_double(*pr_cal_det))
        deltas_flex.extend(pr_obs_lab - pr_cal_lab)
      d = [flex.mean(flex.abs(deltas_flex.parts()[i])) for i in range(3)]
      return_dict.update({'delta_x': d[0], 'delta_y': d[1], 'delta_z': d[2]})
    return return_dict

  def extract_unit_cell_distribution(self, scaling_expt_paths):
    """Retrieve average a, b, c and their deltas using expt paths"""
    af, bf, cf = flex.double(), flex.double(), flex.double()
    with tempfile.NamedTemporaryFile() as tdata_file:
      self._write_tdata(scaling_expt_paths, tdata_file.name)
      with open(tdata_file.name, 'r') as tdata:
        for line in tdata.read().splitlines():
          a, b, c = line.strip().split(' ')[:3]
          af.append(float(a))
          bf.append(float(b))
          cf.append(float(c))
    return {'a': flex.mean(af), 'b': flex.mean(bf), 'c': flex.mean(cf),
            'delta_a': af.standard_deviation_of_the_sample(),
            'delta_b': bf.standard_deviation_of_the_sample(),
            'delta_c': cf.standard_deviation_of_the_sample()}

  def locate_input_tags(self):
    input_paths = []
    for ig in self.parameters.input.glob:
        input_paths.extend(glob.glob(ig))
    return input_paths

  @staticmethod
  def locate_scaling_directories(merging_phil_paths):
    """Return paths to all directories specified as input.path in phil file"""
    merging_phils = [parse(file_name=mpp) for mpp in merging_phil_paths]
    phil = DEFAULT_INPUT_SCOPE.fetch(sources=merging_phils).extract()
    return sorted(set(phil.input.path))

  @staticmethod
  def locate_tder_refined_expts_and_refls(scaling_phil_paths):
    """Return paths to refined expt and refl files mentioned in scaling phil"""
    expt_paths, refl_paths = [], []
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

  def _write_tdata(self, expt_paths, tdata_path):
    """Read all expt_paths and write a tdata file with unit cells in lines"""
    s = '{:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {}'
    tdata_lines = []
    for expt in self.load_experiments(*expt_paths):
      uc_params = expt.crystal.get_unit_cell().parameters()
      sg = expt.crystal.get_space_group().type().universal_hermann_mauguin_symbol()
      tdata_lines.append(s.format(*uc_params, sg.replace(' ', '')))
    with open(tdata_path, 'w') as tdata_file:
      tdata_file.write('\n'.join(tdata_lines))

  def scrap(self):
    for tag in self.locate_input_tags():
      last_v_dir = sorted(self.path_lookup(tag, 'v*/'))[-1]
      merging_phil_files = self.path_lookup(last_v_dir, '*_params.phil')
      for scaling_dir in self.locate_scaling_directories(merging_phil_files):
        scaling_expts = self.path_lookup(scaling_dir, 'scaling_*.expt')
        try:
          first_scaling_expt_path = sorted(scaling_expts)[0]
        except IndexError:
          continue
        task_dir = self.path_join(first_scaling_expt_path, '..', '..')
        trial_dir = self.path_join(task_dir, '..')
        scrap_dict = {'tag': tag, 'trial': int(trial_dir[-9:-6]),
                      'task': int(self.path_split(task_dir)[-1].lstrip('task')),
                      'rungroup': int(trial_dir[-3:]),
                      'run': self.path_split(trial_dir)[-2]}
        print('Processing run {} in tag {}'.format(scrap_dict['run'], tag))
        scrap_dict.update(self.extract_origin(first_scaling_expt_path))
        if self.parameters.uncertainties:
          scaling_phils = self.path_lookup(task_dir, 'params_1.phil')
          expt_p, refl_p = self.locate_tder_refined_expts_and_refls(scaling_phils)
          scrap_dict.update(self.extract_size_and_origin_deltas(expt_p, refl_p))
          scrap_dict.update(self.extract_unit_cell_distribution(scaling_expts))
        self.table.add(**scrap_dict)


class DriftTable(object):
  """Class responsible for storing info about all DriftDataclass instances"""
  KEYS = ['tag', 'run', 'rungroup', 'trial', 'task', 'x', 'y', 'z',
          'delta_x', 'delta_y', 'delta_z', 'expts', 'refls',
          'a', 'b', 'c', 'delta_a', 'delta_b', 'delta_c']

  def __init__(self, parameters):
    self.data = []
    self.parameters = parameters

  def __getitem__(self, key):
    auxiliary_key = key in self.auxiliary.keys()
    return self.auxiliary[key] if auxiliary_key else [d[key] for d in self.data]

  def __str__(self):
    lines = [' '.join('{!s:9.9}'.format(k.upper()) for k in self.available_keys)]
    for i, d in enumerate(self.data):
      cells = ['{:.20f}'.format(d[k]) if isinstance(d[k], float) else d[k]
               for k in self.active_keys]
      for ak in self.auxiliary.keys():
        cells.append('{:.20f}'.format(self[ak][i]))
      lines.append(' '.join('{!s:9.9}'.format(cell) for cell in cells))
    return '\n'.join(lines)

  def add(self, **kwargs):
    self.data.append(kwargs)

  def get(self, key, default=None):
    return self[key] if key in self.available_keys else default

  def sort(self, by_key=KEYS[0]):
    self.data = sorted(self.data, key=lambda d: d[by_key])

  @property
  def active_keys(self):
    return [key for key in self.KEYS if not all(v is None for v in self[key])]

  @property
  def available_keys(self):
    return self.active_keys + list(self.auxiliary.keys())

  @property
  def auxiliary(self):
    density = [d['refls'] / d['expts'] if d['expts'] else 0 for d in self.data]
    return {'density': density}


class DriftArtist(object):
  """Object responsible for plotting an instance of `DriftTable`."""
  def __init__(self, table, parameters):
    self.colormap = plt.get_cmap('tab10')
    self.colormap_period = 10
    self.color_by = 'tag'
    self.order_by = 'run'
    self.cov_colormap = plt.get_cmap('seismic')
    self.table = table
    self.parameters = parameters
    self._init_figure()
    self._setup_figure()

  def _init_figure(self):
    self.fig = plt.figure(tight_layout=True)
    gs = GridSpec(7, 2, hspace=0, wspace=0, width_ratios=[4, 1],
                  height_ratios=[2, 3, 3, 3, 3, 3, 3])
    self.axh = self.fig.add_subplot(gs[0, 0])
    self.axx = self.fig.add_subplot(gs[1, 0], sharex=self.axh)
    self.axy = self.fig.add_subplot(gs[2, 0], sharex=self.axh)
    self.axz = self.fig.add_subplot(gs[3, 0], sharex=self.axh)
    self.axa = self.fig.add_subplot(gs[4, 0], sharex=self.axh)
    self.axb = self.fig.add_subplot(gs[5, 0], sharex=self.axh)
    self.axc = self.fig.add_subplot(gs[6, 0], sharex=self.axh)
    self.axw = self.fig.add_subplot(gs[0, 1])
    self.axl = self.fig.add_subplot(gs[1:, 1])

  def _setup_figure(self):
    self.axl.axis('off')
    self.axw.axis('off')
    self.axh.spines['top'].set_visible(False)
    self.axh.spines['right'].set_visible(False)
    self.axh.spines['bottom'].set_visible(False)
    self.axh.set_zorder(self.axx.get_zorder() - 1.0)
    common = {'direction': 'inout', 'top': True, 'bottom': True, 'length': 6}
    main_axes = self.axx, self.axy, self.axz, self.axa, self.axb, self.axc
    for ax, label in zip(main_axes, ['X', 'Y', 'Z', 'a', 'b', 'c']):
      ax.set_ylabel(label)
      ax.tick_params(axis='x', labelbottom=False, **common)
      ax.ticklabel_format(useOffset=False)
    self.axc.tick_params(axis='x', labelbottom=True, rotation=90)
    self.axc.set_xlabel(self.order_by.title())
    self.axh.set_ylabel('# expts')

  @property
  def _refl_to_expt_ratios(self):
    return [r / e for e, r in zip(self.table['expts'], self.table['refls'])]

  @property
  def color_array(self):
    """Registry-length color list with colors corresponding to self.color_by"""
    color_i = [0, ] * len(self.table.data)
    color_by = self.table[self.color_by]
    for i, cat in enumerate(color_by[1:], 1):
      color_i[i] = color_i[i-1] if cat in color_by[:i] else color_i[i-1] + 1
    return [self.colormap(i % self.colormap_period) for i in color_i]

  @property
  def x(self):
    return self.table[self.order_by]

  def _get_handles_and_labels(self):
    handles, unique_keys = [], []
    for key in self.table[self.color_by]:
      if key not in unique_keys:
        handles.append(Line2D([], [], c=self.colormap(len(unique_keys) % 10),
                              ls='', ms=12, marker='.'))
        unique_keys.append(key)
    return handles, unique_keys

  def _plot_drift(self, axes, values_key, deltas_key=None, top=False):
    y = self.table[values_key]
    y_err = self.table.get(deltas_key, [])
    axes.scatter(self.x, y, c=self.color_array)
    if top:
      ax_top = self.axx.secondary_xaxis('top')
      ax_top.tick_params(rotation=90)
      ax_top.set_xticks(self.axx.get_xticks())
      ax_top.set_xticklabels(self.table['expts'])
    if self.parameters.uncertainties:
      axes.errorbar(self.x, y, yerr=y_err, ecolor='black', ls='')
    axes2 = axes.twinx()
    avg_y = average(y, weights=self.table['refls'])
    axes2.set_ylim([lim / avg_y - 1 for lim in axes.get_ylim()])
    axes2.yaxis.set_major_formatter(PercentFormatter(xmax=1))

  def _plot_bars(self):
    y = self.table['expts']
    w = normalise([0, *self.table['density']])[1:]
    self.axh.bar(self.x, y, width=w, color=self.color_array, alpha=0.5)

  def _plot_correlations(self):
    keys = ['x', 'y', 'z', 'a', 'b', 'c']
    cm = CorrelationMatrix({k: self.table[k] for k in keys},
                           weights=self.table['refls'])
    print(cm)
    self.axw.set_xlim([0, len(keys)])
    self.axw.set_ylim([0, len(keys)])
    for ix, kx in enumerate(keys):
      for iy, ky in enumerate(keys):
        if ix == iy:
          self.axw.text(x=ix+0.5, y=len(keys)-iy-0.5, s=kx,
                        ha='center', va='center')
        if ix > iy:
          print('Calculating correlation between {} and {}:'.format(kx, ky))
          print('{} values: {}'.format(kx, list(self.table[kx])))
          print('{} values: {}'.format(ky, list(self.table[ky])))
          corr = cm.corr[kx][ky]
          color = self.cov_colormap(normalise([corr, -1, 1])[0])
          r = Rectangle(xy=(ix, len(keys) - iy), width=1, height=-1, fill=True,
                        ec='white', fc=color, linewidth=2)
          print('Calculated corr: {}; color: {}; pos: {}'.format(corr, color,
                                                        (ix, len(keys) - iy)))
          print('-' * 80 + '\n')
          self.axw.add_patch(r)

  def _plot_legend(self):
    handles, labels = self._get_handles_and_labels()
    self.axl.legend(handles, labels, loc=7)

  def _plot_width_info(self):
    extrema = [min(self.table['expts']), max(self.table['expts']),
               min(self.table['refls']), max(self.table['refls'])]
    s = "#expts/run: {} - {}\n#refls/run: {} - {}".format(*extrema)
    self.axl.text(x=0.5, y=0.0, s=s, clip_on=False, ha='center',
                  ma='center', va='top', transform=self.axl.transAxes)

  def publish(self):
    if self.table.data:
      self._plot_bars()
      self._plot_correlations()
      self._plot_width_info()
      self._plot_drift(self.axx, 'x', 'delta_x', top=True)
      self._plot_drift(self.axy, 'y', 'delta_y')
      self._plot_drift(self.axz, 'z', 'delta_z')
      self._plot_drift(self.axa, 'a', 'delta_a')
      self._plot_drift(self.axb, 'b', 'delta_b')
      self._plot_drift(self.axc, 'c', 'delta_c')
      self._plot_legend()
    self.fig.align_labels()
    if self.parameters.plot.path:
      self.fig.set_size_inches(self.parameters.plot.width,
                               self.parameters.plot.height)
      self.fig.savefig(self.parameters.plot.path)
    if self.parameters.plot.show:
      plt.show()


def run(params_):
  dt = DriftTable(parameters=params_)
  ds = DriftScraper(table=dt, parameters=params_)
  da = DriftArtist(table=dt, parameters=params_)
  ds.scrap()
  dt.sort(by_key='run')
  print(dt)
  da.publish()


if __name__ == '__main__':
  if '--help' in sys.argv[1:] or '-h' in sys.argv[1:]:
    print(message)
    exit()
  params = params_from_phil(sys.argv[1:])
run(params)
