from __future__ import absolute_import, division, print_function

from collections import namedtuple
import glob
import os
from six.moves import UserDict, UserList

import numpy as np
from orderedset import OrderedSet


class StemLocator(UserDict):
  """Subclass of dict which raises an error when overwriting existing value"""

  class StemExistsError(FileExistsError):
    MSG = 'Multiple files with stem {key} detected. Matching expts/refls ' \
          'must EITHER be pairwise located in the same directory' \
          'OR have unique names across all input directories.'

    @classmethod
    def from_stem(cls, stem, *args, **kwargs):
      return cls(cls.MSG.format(stem), *args, **kwargs)  # noqa

  def __setitem__(self, key, value):
    if key in self.data:
      if value != self.data.__getitem__(key):
        raise self.StemExistsError.from_stem(key)
    else:
      self.data.__setitem__(key, value)


def path_pair_class_factory(params):
  """Factory function, dynamically makes `PathPair` class based on params"""
  es = params.input.experiments_suffix
  rs = params.input.reflections_suffix

  class PathPair(namedtuple('PathPair', ['expt_path', 'refl_path'])):
    """Named 2-el. tuple of expt & refl file paths with convenience methods"""
    __slots__ = ()

    @classmethod
    def from_dir_expt_name(cls, base_path, expt_filename):
      expt_path = os.path.join(base_path, expt_filename)
      refl_path = expt_path.split(es, 1)[0] + rs
      refl_path = refl_path if os.path.exists(refl_path) else None
      return cls(expt_path, refl_path)

    @classmethod
    def from_dir_refl_name(cls, base_path, refl_filename):
      refl_path = os.path.join(base_path, refl_filename)
      expt_path = refl_path.split(rs, 1)[0] + es
      expt_path = expt_path if os.path.exists(expt_path) else None
      return cls(expt_path, refl_path)

    @property
    def expt_stem(self):
      return os.path.basename(self.expt_path).split(es, 1)[0]

    @property
    def refl_stem(self):
      return os.path.basename(self.refl_path).split(rs, 1)[0]

  return PathPair


class PathPairList(UserList):
  """List of `PathPair` instances with some convenience creation methods"""
  def __init__(self, params):
    """On creation, fill self with expt/refl `PathPair`s based on `params`"""
    super(PathPairList, self).__init__()
    self.params = params
    self.PathPair = path_pair_class_factory(params)
    self.reject_alisters = self.params.input.alist.op == "reject"
    self.keep_alisters = not self.reject_alisters
    self.alist = set()  # contains image tags or absolute experiment paths
    self._load_alist()
    self._load_paths()
    self._match_singlets()

  def _load_alist(self):
    alist_files = self.params.input.alist.file
    if alist_files is not None:
      for f in alist_files:
        lines = np.loadtxt(f, str, ndmin=1)
        self.alist = self.alist.union(lines)
    if self.params.input.alist.type == "files":
      assert all(os.path.isabs(f) for f in self.alist)

  def _is_on_alist(self, exp_path):
    if self.params.input.alist.type == "tags":
      return any(tag in exp_path for tag in self.alist)
    else:  # alist.type == "files"
      return os.path.abspath(exp_path) in self.alist

  @property
  def num_alist(self):
    return len(self.alist)

  def _load_path_if_expt_or_refl(self, path, filename):
    if filename.endswith(self.params.input.experiments_suffix):
      self.data.append(self.PathPair.from_dir_expt_name(path, filename))
    if filename.endswith(self.params.input.reflections_suffix):
      self.data.append(self.PathPair.from_dir_refl_name(path, filename))

  def _load_paths(self):
    for pathstring in self.params.input.path:
      for path in glob.glob(pathstring):
        if os.path.isdir(path):
          for filename in os.listdir(path):
            self._load_path_if_expt_or_refl(path, filename)
        else:
          filename, dir_name = os.path.split(path)
          self._load_path_if_expt_or_refl(dir_name, filename)
    accepted_pairs = []
    for pp in OrderedSet(self.data):
      if pp.expt_path is None or self._is_accepted_expt(pp.expt_path):
        accepted_pairs.append(pp)
    self.data = accepted_pairs

  def _is_accepted_expt(self, experiment_path):
    """returns True if the experiment is deemed worthy"""
    if self.num_alist > 0:
      is_on_alist = self._is_on_alist(experiment_path)
      if is_on_alist and self.reject_alisters:
        return False
      if self.keep_alisters and not is_on_alist:
        return False
    return True

  def _match_singlets(self):
    """Merge every matching pair of PathPair(expt, None) + PathPair(None, refl)
    in self.data with common stub name into a single PathPair(expt, refl)"""
    doublets, expt_singlets, refl_singlets = [], StemLocator(), StemLocator()
    for path_pair in self.data:
      if path_pair.expt_path and path_pair.refl_path:
        doublets.append(path_pair)
      elif path_pair.expt_path and not path_pair.refl_path:
        expt_singlets[path_pair.expt_stem] = path_pair.expt_path
      elif path_pair.refl_path and not path_pair.expt_path:
        refl_singlets[path_pair.refl_stem] = path_pair.refl_path
    common = OrderedSet(expt_singlets).intersection(OrderedSet(refl_singlets))
    new = [self.PathPair(expt_singlets[c], refl_singlets[c]) for c in common]
    self.data = doublets + new
