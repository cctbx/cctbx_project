from __future__ import absolute_import, division, print_function
import glob, os
from six.moves import UserDict
import numpy as np


class ImmutableValuesDict(UserDict):
  def __setitem__(self, key, value):
    try:
      old_value = self.__getitem__(key)
      if value != old_value:
        raise ValueError("Attempting to change value of existing key: ", key)
    except KeyError:
      super(ImmutableValuesDict, self).__setitem__(key, value)


class file_lister(object):
  """ Class for matching jsons to reflection table pickles """
  def __init__(self, params):
    self.params = params
    self.reject_alisters = self.params.input.alist.op == "reject"
    self.keep_alisters = not self.reject_alisters
    self.alist = set()  # contains image tags or absolute experiment paths
    self._load_alist()

  def _load_alist(self):
    alist_files = self.params.input.alist.file
    if alist_files is not None:
      for f in alist_files:
        lines = np.loadtxt(f, str, ndmin=1)
        self.alist = self.alist.union(lines)

    if self.params.input.alist.type=="files":
      assert all(os.path.isabs(f) for f in self.alist)

  def _is_on_alist(self, exp_path):
    if self.params.input.alist.type=="tags":
      return any(tag in exp_path for tag in self.alist)
    else:  # alist.type=="files"
      return os.path.abspath(exp_path) in self.alist

  @property
  def num_alist(self):
    return len(self.alist)

  @property
  def filepair_generator(self):
    """ Used by rank 0, client/server to yield a list of json/reflection table pairs """
    if self.params.input.match_directories:
      return self._matching_path_name_filepair_generator
    else:
      return self._matching_file_name_filepath_generator

  def _matching_path_name_filepair_generator(self):
    """File pair generator used when input expt/refl must come from same dir"""
    for pathstring in self.params.input.path:
      for path in glob.glob(pathstring):
        if os.path.isdir(path):
          for filename in os.listdir(path):
            if filename.endswith(self.params.input.reflections_suffix):
              experiments_path = os.path.join(path, filename.split(self.params.input.reflections_suffix)[0] +
                 self.params.input.experiments_suffix)
              if self._is_accepted_expt(experiments_path):
                yield experiments_path, os.path.join(path, filename)
        else:
          dirname = os.path.dirname(path)
          filename = os.path.basename(path)
          if filename.endswith(self.params.input.reflections_suffix):
            experiments_path = os.path.join(dirname, filename.split(self.params.input.reflections_suffix)[0] +
               self.params.input.experiments_suffix)
            if self._is_accepted_expt(experiments_path):
              yield experiments_path, path

  def _matching_file_name_filepath_generator(self):
    """File pair generator used when input expt/refl can be in different dir"""
    expt_paths = ImmutableValuesDict()
    refl_paths = ImmutableValuesDict()
    expt_suffix_len = len(self.params.input.experiments_suffix)
    refl_suffix_len = len(self.params.input.reflections_suffix)
    for pathstring in self.params.input.path:
      for path in glob.glob(pathstring):
        if os.path.isdir(path):
          for filename in os.listdir(path):
            if filename.endswith(self.params.input.reflections_suffix):
              file_stem = filename[:-refl_suffix_len]
              refl_paths[file_stem] = os.path.join(path, filename)
            if filename.endswith(self.params.input.experiments_suffix):
              if self._is_accepted_expt(filename):
                file_stem = filename[:-expt_suffix_len]
                expt_paths[file_stem] = os.path.join(path, filename)
        else:
          filename = os.path.basename(path)
          if filename.endswith(self.params.input.reflections_suffix):
            file_stem = filename[:-refl_suffix_len]
            refl_paths[file_stem] = path
          if filename.endswith(self.params.input.experiments_suffix):
            if self._is_accepted_expt(filename):
              file_stem = filename[:-expt_suffix_len]
              expt_paths[file_stem] = path
    for common_file_stem in set(expt_paths).intersection(set(refl_paths)):
      yield expt_paths[common_file_stem], refl_paths[common_file_stem]


  def _is_accepted_expt(self, experiment_path):
    """returns True if the experiment is deemed worthy"""
    if not os.path.exists(experiment_path):
      return False
    if self.num_alist > 0:
      is_on_alist = self._is_on_alist(experiment_path)
      if is_on_alist and self.reject_alisters:
        return False
      if self.keep_alisters and not is_on_alist:
        return False

    return True

  def filename_lister(self):
    """ Used by rank 0, striping """
    return list(self.filepair_generator())

