from __future__ import absolute_import, print_function, division
import abc
import itertools
from collections import OrderedDict
from json.decoder import JSONDecodeError
import glob
import os
import six
import sys
import tempfile

from dials.array_family import flex  # noqa
from dxtbx.model.experiment_list import ExperimentList  # noqa
from libtbx.phil import parse
from libtbx.utils import Sorry

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import PercentFormatter
import pandas as pd


message = '''
This script aims to investigate the spatial drift of a detector as a function
of experimental progress. It requires the directory structure to follow that
resulting from indexing and ensemble refinement performed by cctbx.xfel.
End result is a multiplot with detector origin position (vertical position) as
a function of run number (horizontal position), colored according to tag name.
Error bars can be derived from the uncertainty of individual reflections'
position in laboratory reference system.

Example usage 1:
"libtbx.python `libtbx.find_in_repositories xfel`/util/drift.py
input.glob=batch*TDER/ input.kind=merging_directory"
where `batch*TDER` points to folders with merging results.

Example usage 2:
"libtbx.python `libtbx.find_in_repositories xfel`/util/drift.py
input.glob=r0*/039_rg084/task209 input.kind=tder_task_directory"
where `r0*/039_rg084/task209` points to folders with TDER results.
'''


################################ PHIL HANDLING ################################


phil_scope = parse('''
  input {
    glob = None
      .type = str
      .multiple = True
      .help = glob which matches files after TDER to be investigated.
    exclude = None
      .type = str
      .multiple = True
      .help = glob which matches files to exclude from input.glob
    kind = tder_task_directory *merging_directory
      .type = choice
      .help = The type of files located by input.glob
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
  origin = *first average distribution
    .type = choice
    .help = Use origin of first experiment only or average/distribution of all?
  uncertainties = True
    .type = bool
    .help = If True, uncertainties will be estimated using differences \
          between predicted and observed refl positions and cell distribution
  unit_cell = *average distribution
    .type = choice
    .help = Use average unit cell or distribution of all?
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


############################### MATHS FUNCTIONS ###############################


def average(sequence, weights=None):
  """Calculate weighted arithmetic mean of an iterable"""
  weights = [1] * len(sequence) if weights is None else weights
  return sum(s * w for s, w in zip(sequence, weights)) / sum(weights)


def correlation(xs, ys, weights=None):
  """Calculate weighted Pearson's correlation coefficient between iterables"""
  weights = [1] * len(xs) if weights is None else weights
  x_variance = variance(xs, xs, weights=weights)
  y_variance = variance(ys, ys, weights=weights)
  covariance = variance(xs, ys, weights=weights)
  return covariance / (x_variance * y_variance) ** 0.5


def variance(xs, ys, weights=None):
  """Calculate weighted variance between iterables"""
  weights = [1] * len(xs) if weights is None else weights
  x_avg = average(xs, weights=weights)
  y_avg = average(ys, weights=weights)
  x_dev = [x - x_avg for x in xs]
  y_dev = [y - y_avg for y in ys]
  return sum(w * x * y for w, x, y in zip(weights, x_dev, y_dev)) / sum(weights)


def normalize(sequence, floor=0, ceiling=1):
  """Normalize `sequence`'s values to lie between `floor` and `ceiling`"""
  min_, max_ = min(sequence), max(sequence)
  old_span = max_ - min_
  new_span = ceiling - floor
  normalized = [new_span * (s - min_) / old_span + floor for s in sequence]
  return sequence if min_ == max_ else normalized


class CorrelationMatrix(object):
  def __init__(self, variables, weights=None):
    """Calculate corr. matrix between every pair in `variables` dict-like"""
    self.keys = variables.keys()
    self.corr = {k: {} for k in self.keys}
    for i1, k1 in enumerate(self.keys):
      for i2, k2 in enumerate(self.keys):
        if i1 == i2:
          self.corr[k1][k2] = self.corr[k2][k1] = 1.0
        elif len(variables[k1]) < 2:
          self.corr[k1][k2] = self.corr[k2][k1] = 0.0
        elif i2 > i1:
          try:
            corr = correlation(variables[k1], variables[k2], weights=weights)
          except ZeroDivisionError:
            corr = 0
          self.corr[k1][k2] = self.corr[k2][k1] = corr

  def __str__(self):
    s = 'Correl. ' + ' '.join('{:>7}'.format(k) for k in self.keys)
    for k1 in self.keys:
      s += '\n{:7}'.format(k1)
      for k2 in self.keys:
        s += ' {:+6.4f}'.format(self.corr[k1][k2])
    return s


############################## UTILITY FUNCTIONS ##############################


def is_iterable(value):
  """Return `True` if object is iterable, else `False`"""
  try:
    iter(value)
  except TypeError:
    return False
  else:
    return True


def path_join(*path_elements):
  """Join path from elements, resolving all redundant or relative calls"""
  path_elements = [os.pardir if p == '..' else p for p in path_elements]
  return os.path.normpath(os.path.join(*path_elements))


def path_lookup(*path_elements):
  """Join path elements and return a list of all matching files/dirs"""
  return glob.glob(path_join(*path_elements), recursive=True)


def path_split(path):
  """Split path into directories and file basename using os separator"""
  return os.path.normpath(path).split(os.sep)


def read_experiments(*expt_paths):
  """Create an instance of ExperimentList from one or more `*expt_paths`"""
  expts = ExperimentList()
  for expt_path in unique_elements(expt_paths):
    expts.extend(ExperimentList.from_file(expt_path, check_format=False))
  return expts


def read_reflections(*refl_paths):
  """Create an instance of flex.reflection_table from one/more `*refl_paths`"""
  r = [flex.reflection_table.from_file(p) for p in unique_elements(refl_paths)]
  return flex.reflection_table.concat(r)


def represent_range_as_str(sorted_iterable):
  """Return str in only one in iterable, range e.g. "r00[81-94]" otherwise"""
  fs, ls = str(sorted_iterable[0]), str(sorted_iterable[-1])
  d = min([i for i, (fl, ll) in enumerate(zip(fs, ls)) if fl != ll] or [None])
  return fs if not d else fs[:d] + '[' + fs[d:] + '-' + ls[d:] + ']'


def unique_elements(sequence):
  """Return list of unique elements in sequence while preserving its order"""
  return list(OrderedDict.fromkeys(sequence))


############################### DRIFT SCRAPPING ###############################


class DriftScraperRegistrar(abc.ABCMeta):
  """Metaclass for `DriftScraper`s, auto-registers them by `input_kind`."""
  REGISTRY = {}
  def __new__(mcs, name, bases, attrs):
    new_cls = super().__new__(mcs, name, bases, attrs)
    if hasattr(new_cls, 'input_kind') and new_cls.input_kind:
      mcs.REGISTRY[new_cls.input_kind] = new_cls
    return new_cls


@six.add_metaclass(DriftScraperRegistrar)
class BaseDriftScraper(object):
  """Base class for scraping cctbx.xfel output into instance of `DriftTable`,
  with automatic registration into the `DriftScraperRegistrar`."""

  def __init__(self, table, parameters):
    self.table = table
    self.parameters = parameters

  @staticmethod
  def calc_expt_refl_len(expts, refls):
    expts_len = len(expts)
    refls_len = []
    for expt_id, expt in enumerate(expts):
      refls_len.append((expt_id == refls['id']).count(True))
    return expts_len, refls_len

  @staticmethod
  def extract_db_metadata(combine_phil_path):
    """Get trial, task, rungroup, chunk, run info based on combining phil"""
    parsed_combine_phil = parse(file_name=combine_phil_path)
    phil = DEFAULT_INPUT_SCOPE.fetch(sources=[parsed_combine_phil]).extract()
    index_dirs = [path_join(pie, '..') for pie in phil.input.experiments]
    rungroups = sorted(set(index_dir[-7:-4] for index_dir in index_dirs))
    trials = sorted(set(index_dir[-13:-10] for index_dir in index_dirs))
    runs = sorted(set(index_dir[-19:-14] for index_dir in index_dirs))
    return {'chunk': path_split(combine_phil_path)[-1][16:19],
            'run': represent_range_as_str(runs),
            'rungroup': represent_range_as_str(rungroups),
            'task': path_split(combine_phil_path)[-4],
            'trial': represent_range_as_str(trials)}

  def locate_input_paths(self):
    """Return all paths (either common files or directories, relative to
    working directory) of specified kind to be processed"""
    input_paths, exclude_paths = [], []
    for ig in self.parameters.input.glob:
      input_paths.extend(glob.glob(ig))
    for ie in self.parameters.input.exclude:
      exclude_paths.extend(glob.glob(ie))
    return [it for it in input_paths if it not in exclude_paths]

  @staticmethod
  def locate_combining_phil_paths(scaling_phil_paths):
    """Return paths to all phil files used to combine later-scaled expts"""
    parsed_scaling_phil = [parse(file_name=spp) for spp in scaling_phil_paths]
    phil = DEFAULT_INPUT_SCOPE.fetch(sources=parsed_scaling_phil).extract()
    combine_dirs = [path_join(ip, '..') for ip in phil.input.path]
    combine_phil_paths = []
    for cd in combine_dirs:
      combine_phil_paths.extend(path_lookup(cd, '*chunk*_combine_*.phil'))
    return sorted(set(combine_phil_paths))

  @staticmethod
  def locate_scaling_directories(merging_phil_paths):
    """Return paths to all directories specified as input.path in phil file"""
    merging_phils = [parse(file_name=mpp) for mpp in merging_phil_paths]
    phil = DEFAULT_INPUT_SCOPE.fetch(sources=merging_phils).extract()
    return sorted(set(phil.input.path))

  @staticmethod
  def locate_refined_expts_refls(combine_phil_path):
    """Return all refined expts and refls down-stream from combine_phil_path"""
    path_stem = combine_phil_path.replace('_combine_experiments.phil', '')
    expts_paths = path_lookup(path_stem + '_refined.expt')
    refls_paths = path_lookup(path_stem + '_refined.refl')
    expts = read_experiments(*expts_paths)
    refls = read_reflections(*refls_paths)
    return expts, refls

  @abc.abstractmethod
  def scrap(self):
    """Fill `self.table` based on `self.parameters` provided."""
    pass


class TderTaskDirectoryDriftScraper(BaseDriftScraper):
  input_kind = 'tder_task_directory'

  def scrap(self):
    combining_phil_paths = []
    for tder_task_directory in self.locate_input_paths():
      cpp = path_lookup(tder_task_directory, 'combine_experiments_t*',
                        'intermediates', '*chunk*_combine_*.phil')
      combining_phil_paths.extend(cpp)
    for cpp in unique_elements(combining_phil_paths):  # combine.phil paths
      try:
        scrap_dict = {'tag': "TDER"}
        scrap_dict.update(self.extract_db_metadata(cpp))
        print('Processing run {}'.format(scrap_dict['run']))
        refined_expts, refined_refls = self.locate_refined_expts_refls(cpp)
        elen, rlen = self.calc_expt_refl_len(refined_expts, refined_refls)
        print(f'Found {elen} expts and {sum(rlen)} refls')
        scrap_dict.update({'expts': elen, 'refls': rlen})
        scrap_dict.update(self.get_origin(refined_expts))
        scrap_dict.update(self.get_unit_cell(refined_expts))
        scrap_dict.update(self.get_origin_deltas(refined_expts, refined_refls))
      except (KeyError, IndexError, JSONDecodeError) as e:
        print(e)
      else:
        self.table.add(**scrap_dict)


class MergingDirectoryDriftScraper(BaseDriftScraper):
  input_kind = 'merging_directory'

  def scrap(self):
    for tag in self.locate_input_paths():
      merging_phil_paths = path_lookup(tag, '**', '*.phil')
      merging_phil_paths.sort(key=os.path.getmtime)
      for scaling_dir in self.locate_scaling_directories(merging_phil_paths):
        scaled_expt_paths = path_lookup(scaling_dir, 'scaling_*.expt')
        scaled_expts = read_experiments(*scaled_expt_paths)
        scaled_identifiers = list(scaled_expts.identifiers())
        scaling_phil_paths = []
        for sep in scaled_expt_paths:
          scaling_phil_paths.extend(path_lookup(sep, '..', '..', '*.phil'))
        comb_phil_paths = self.locate_combining_phil_paths(scaling_phil_paths)
        for cpp in unique_elements(comb_phil_paths):
          try:
            scrap_dict = {'tag': tag}
            scrap_dict.update(self.extract_db_metadata(cpp))
            print('Processing run {} in tag {}'.format(scrap_dict['run'], tag))
            refined_expts, refined_refls = self.locate_refined_expts_refls(cpp)
            elen, rlen = self.calc_expt_refl_len(refined_expts, refined_refls)
            print(f'Found {elen} expts and {sum(rlen)} refls')
            refined_expts.select_on_experiment_identifiers(scaled_identifiers)
            refined_refls = refined_refls.select(refined_expts)
            elen, rlen = self.calc_expt_refl_len(refined_expts, refined_refls)
            print(f'Accepted {elen} expts and {sum(rlen)} refls')
            scrap_dict.update({'expts': elen, 'refls': rlen})
            scrap_dict.update(self.get_origin(refined_expts))
            scrap_dict.update(self.get_unit_cell(refined_expts))
            scrap_dict.update(self.get_origin_deltas(refined_expts, refined_refls))
          except (KeyError, IndexError, JSONDecodeError) as e:
            print(e)
          else:
            self.table.add(**scrap_dict)


class FirstOriginMixin(object):
  @staticmethod
  def get_origin(expts):
    """Read detector origin (x, y, z) from the first expt file"""
    x, y, z = expts[0].detector.hierarchy().get_origin()
    return {'x': x, 'y': y, 'z': z}


class AverageOriginMixin(object):
  @staticmethod
  def get_origin(expts):
    """Read detector origin (x, y, z) from the first expt file"""
    xs, ys, zs = flex.double(), flex.double(), flex.double()
    for expt in expts:
      x, y, z = expt.detector.hierarchy().get_origin()
      xs.append(x)
      ys.append(y)
      zs.append(z)
    return {'x': flex.mean(xs), 'y': flex.mean(ys), 'z': flex.mean(zs)}


class DistributionOriginMixin(object):
  @staticmethod
  def get_origin(expts):
    """Read detector origin (x, y, z) from the first expt file"""
    xs, ys, zs = flex.double(), flex.double(), flex.double()
    for expt in expts:
      x, y, z = expt.detector.hierarchy().get_origin()
      xs.append(x)
      ys.append(y)
      zs.append(z)
    return {'x': xs, 'y': ys, 'z': zs}


class FalseUncertaintiesMixin(object):
  @staticmethod
  def get_origin_deltas(expts, refls):
    """If uncertainties=False, return origin uncertainties as zeros"""
    return {'delta_x': 0., 'delta_y': 0., 'delta_z': 0.}


class TrueUncertaintiesMixin(object):
  @staticmethod
  def get_origin_deltas(expts, refls):
    """Get uncertainties of origin positions from refl. position deviations"""
    deltas_flex = flex.vec3_double()
    for panel in expts[0].detector:
      pr = refls.select(refls['panel'] == panel.index())      # pr: panel refls
      pr_obs_det = pr['xyzobs.mm.value'].parts()[0:2]  # det: in detector space
      pr_cal_det = pr['xyzcal.mm'].parts()[0:2]        # lab: in labor. space
      pr_obs_lab = panel.get_lab_coord(flex.vec2_double(*pr_obs_det))
      pr_cal_lab = panel.get_lab_coord(flex.vec2_double(*pr_cal_det))
      deltas_flex.extend(pr_obs_lab - pr_cal_lab)
    d = [flex.mean(flex.abs(deltas_flex.parts()[i])) for i in range(3)]
    return {'delta_x': d[0], 'delta_y': d[1], 'delta_z': d[2]}


class BaseUnitCellMixin(object):
  @staticmethod
  def _write_tdata(expts, tdata_path):
    """Read all expt_paths and write a tdata file with unit cells in lines"""
    s = '{:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {}'
    tdata_lines = []
    for expt in expts:
      uc_params = expt.crystal.get_unit_cell().parameters()
      sg = expt.crystal.get_space_group().type().universal_hermann_mauguin_symbol()
      tdata_lines.append(s.format(*uc_params, sg.replace(' ', '')))
    with open(tdata_path, 'w') as tdata_file:
      tdata_file.write('\n'.join(tdata_lines))


class AverageUnitCellMixin(BaseUnitCellMixin):
  def get_unit_cell(self, expts):
    """Retrieve average a, b, c and their deltas using expt paths"""
    af, bf, cf = flex.double(), flex.double(), flex.double()
    with tempfile.NamedTemporaryFile() as tdata_file:
      self._write_tdata(expts, tdata_file.name)
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


class DistributionUnitCellMixin(BaseUnitCellMixin):
  def get_unit_cell(self, expts):
    """Retrieve average a, b, c and their deltas using expt paths"""
    af, bf, cf = flex.double(), flex.double(), flex.double()
    with tempfile.NamedTemporaryFile() as tdata_file:
      self._write_tdata(expts, tdata_file.name)
      with open(tdata_file.name, 'r') as tdata:
        for line in tdata.read().splitlines():
          a, b, c = line.strip().split(' ')[:3]
          af.append(float(a))
          bf.append(float(b))
          cf.append(float(c))
    return {'a': af, 'b': bf, 'c': cf,
            'delta_a': af.standard_deviation_of_the_sample(),
            'delta_b': bf.standard_deviation_of_the_sample(),
            'delta_c': cf.standard_deviation_of_the_sample()}


class DriftScraperFactory(object):
  """Produces appropriate DriftScraper class based on phil `parameters`."""
  ORIGIN_MIXINS = {
    'first': FirstOriginMixin,
    'average': AverageOriginMixin,
    'distribution': DistributionOriginMixin,
  }
  UNCERTAINTIES_MIXINS = {
    True: TrueUncertaintiesMixin,
    False: FalseUncertaintiesMixin,
  }
  UNIT_CELL_MIXINS = {
    'average': AverageUnitCellMixin,
    'distribution': DistributionUnitCellMixin,
  }

  @classmethod
  def get_drift_scraper(cls, table, parameters):
    base = DriftScraperRegistrar.REGISTRY[parameters.input.kind]
    mixins = [
      cls.ORIGIN_MIXINS[parameters.origin],
      cls.UNCERTAINTIES_MIXINS[parameters.uncertainties],
      cls.UNIT_CELL_MIXINS[parameters.unit_cell],
    ]
    class DriftScraper(base, *mixins):
      pass
    return DriftScraper(table=table, parameters=parameters)


################################ DRIFT STORAGE ################################


class DriftTable(object):
  """Class responsible for storing all info collected by `DriftScraper`"""
  STATIC_KEYS = ['tag', 'run', 'rungroup', 'trial', 'chunk', 'task',
                 'x', 'y', 'z', 'delta_x', 'delta_y', 'delta_z', 'expts',
                 'refls', 'a', 'b', 'c', 'delta_a', 'delta_b', 'delta_c']
  DYNAMIC_KEYS = ['density']

  def __init__(self):
    self.data = pd.DataFrame()

  def __getitem__(self, key):
    if key in self.DYNAMIC_KEYS:
      self.recalculate_dynamic_column(key)
    return self.data[key]

  def __len__(self):
    return len(self.data.index)

  def __str__(self):
    for key in self.DYNAMIC_KEYS:
      self.recalculate_dynamic_column(key)
    return str(self.data)

  def add(self, *dicts):
    dicts2 = []  # if only 'refls' column is iterable, immediately sum it
    for d in dicts:
      if any([is_iterable(d[k]) for k in d.keys() if k is not 'refls']):
        dicts2.append(d)
      else:
        dicts2.append({k: sum(v) if k is 'refls' else v for k, v in d.items()})
    new_rows = pd.DataFrame(*dicts2)
    self.data = pd.concat([self.data, new_rows], ignore_index=True)

  def get(self, key, default=None):
    return self[key] if key in self.data.columns else default

  def sort(self, by='index'):
    self.data.sort_values(by=by, ignore_index=True)

  @property
  def column_is_flat(self):
    return {k: not is_iterable(self[k][0]) for k in self.data.columns}

  @property
  def flat(self):
    """Drift table with all iterable fields expanded over rows"""
    c = self.column_is_flat
    col_names = self.data.columns
    flatteners = [itertools.repeat if c[k] else lambda x: x for k in col_names]
    flat_table = DriftTable()
    for row in self.data.itertuples(index=False):
      lens = [len(el) for el in row if is_iterable(el)]
      if len(unique_elements(lens)) > 1:
        raise ValueError('All row elements must be scalars or same-length')
      elif len(unique_elements(lens)) == 1:
        flat_columns = zip(*[f(cell) for cell, f in zip(row, flatteners)])
      else:
        flat_columns = [[cell] for cell in row]
      flat_table.add({k: v for k, v in zip(col_names, flat_columns)})
    flat_table.data['expts'] = 1
    return flat_table

  @property
  def is_flat(self):
    return all(self.column_is_flat[k] for k in self.data.columns)

  def recalculate_dynamic_column(self, key):
    if key == 'density':
      refls = self.data['refls'] if self.column_is_flat['refls'] \
        else pd.Series([sum(refl) for refl in self.data['refls']])
      self.data['density'] = self.data['refls'] / self.data['expts']
    else:
      raise KeyError(f'Unknown dynamic column key: {key}')


############################## DRIFT VISUALIZING ##############################


class DriftArtist(object):
  """Object responsible for plotting an instance of `DriftTable`."""
  def __init__(self, table: DriftTable, parameters):
    self.colormap = plt.get_cmap('tab10')
    self.colormap_period = 10
    self.cov_colormap = plt.get_cmap('seismic')
    self.order_by = ['run']
    self.table = table
    self.table.data.sort_values(by=self.order_by, ignore_index=True)
    self.table_flat = self.table.flat
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
    self.axc.set_xlabel(', '.join(ob for ob in self.order_by))
    self.axh.set_ylabel('# expts')

  @property
  def color_array(self):
    """Registry-length color list with colors corresponding to self.color_by"""
    color_i = [0, ] * len(self.table)
    color_by = self.table[self.color_by]
    for i, cat in enumerate(color_by[1:], 1):
      color_i[i] = color_i[i-1] if cat in color_by[:i] else color_i[i-1] + 1
    return [self.colormap(i % self.colormap_period) for i in color_i]

  @property
  def x(self):
    return self.table['index']

  @property
  def color_by(self):
    return 'tag' if len(set(self.table['tag'])) > 1 else 'task'

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
    if y and not self.table.column_is_flat[values_key]:
      self._plot_drift_distribution(axes, y, values_key)
    else:
      self._plot_drift_point(axes, y, deltas_key)
    if top:
      ax_top = self.axx.secondary_xaxis('top')
      ax_top.tick_params(rotation=90)
      ax_top.set_xticks(self.axx.get_xticks())
      ax_top.set_xticklabels(self.table['expts'])
    axes.set_xticklabels(self.table[self.order_by[0]])
    flattened_y = self.table_flat[values_key]
    flattened_weights = self.table_flat['refls']
    avg_y = average(flattened_y, weights=flattened_weights)
    if avg_y != 0:
      axes2 = axes.twinx()
      axes2.set_ylim([lim / avg_y - 1 for lim in axes.get_ylim()])
      axes2.yaxis.set_major_formatter(PercentFormatter(xmax=1))

  def _plot_drift_point(self, axes, y, deltas_key):
    axes.scatter(self.x, y, c=self.color_array)
    y_err = self.table.get(deltas_key, [])
    if self.parameters.uncertainties:
      axes.errorbar(self.x, y, yerr=y_err, ecolor='black', ls='')

  def _plot_drift_distribution(self, axes, y, values_key):
    x_flat = self.table_flat['index']
    y_flat = self.table_flat[values_key]
    b = (len(self.x), 100)
    r = [[-0.5, len(self.x) - 0.5], [min(y_flat), max(y_flat)]]
    axes.hist2d(x_flat, y_flat, bins=b, range=r, cmap=plt.cm.magma_r, cmin=0.5)
    axes.scatter(self.x, [average(y) for y in self.table[values_key]],
                 c=self.color_array, edgecolors='white')

  def _plot_bars(self):
    y = self.table['expts']
    w = normalize([0, *self.table['density']])[1:]
    self.axh.bar(self.x, y, width=w, color=self.color_array, alpha=0.5)

  def _plot_correlations(self):
    keys = ['x', 'y', 'z', 'a', 'b', 'c']
    flat_cols = self.table_flat[keys]
    weights = self.table_flat['refls']
    correlated = {k: f for k, f in zip(keys, flat_cols)}
    cm = CorrelationMatrix(correlated, weights=weights)
    print(cm)
    self.axw.set_xlim([0, len(keys)])
    self.axw.set_ylim([0, len(keys)])
    for ix, kx in enumerate(keys):
      for iy, ky in enumerate(keys):
        if ix == iy:
          self.axw.text(x=ix+0.5, y=len(keys)-iy-0.5, s=kx,
                        ha='center', va='center')
        if ix > iy:
          color = self.cov_colormap(normalize([cm.corr[kx][ky], -1, 1])[0])
          r = Rectangle(xy=(ix, len(keys) - iy), width=1, height=-1, fill=True,
                        ec='white', fc=color, linewidth=2)
          self.axw.add_patch(r)

  def _plot_legend(self):
    handles, labels = self._get_handles_and_labels()
    self.axl.legend(handles, labels, loc=7)

  def _plot_width_info(self):
    expt_lens = self.table['expts']
    refl_lens = self.table['refls'] if self.table.column_is_flat['refls'] \
      else [sum(refl) for refl in self.table['refls']]
    s = f"#expts/run: {min(expt_lens)} - {max(expt_lens)}\n" \
        f"#refls/run: {min(refl_lens)} - {max(refl_lens)}"
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


################################ ENTRY POINTS #################################


def run(params_):
  dt = DriftTable()
  ds = DriftScraperFactory.get_drift_scraper(table=dt, parameters=params_)
  da = DriftArtist(table=dt, parameters=params_)
  ds.scrap()
  dt.sort(by='run')
  print(dt)
  da.publish()


params = []
if __name__ == '__main__':
  if '--help' in sys.argv[1:] or '-h' in sys.argv[1:]:
    print(message)
    exit()
  params = params_from_phil(sys.argv[1:])
run(params)
