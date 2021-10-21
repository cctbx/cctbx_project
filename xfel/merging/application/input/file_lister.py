from __future__ import absolute_import, division, print_function
import glob, os
import numpy as np

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
        lines = np.loadtxt(f, str)
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

  def filepair_generator(self):
    """ Used by rank 0, client/server to yield a list of json/reflection table pairs """
    for pathstring in self.params.input.path:
      for path in glob.glob(pathstring):
        if os.path.isdir(path):
          for filename in os.listdir(path):
            if filename.endswith(self.params.input.reflections_suffix):
              experiments_path = os.path.join(path, filename.split(self.params.input.reflections_suffix)[0] +
                 self.params.input.experiments_suffix)
              if self._is_accepted(experiments_path):
                yield experiments_path, os.path.join(path, filename)
        else:
          dirname = os.path.dirname(path)
          filename = os.path.basename(path)
          if filename.endswith(self.params.input.reflections_suffix):
            experiments_path = os.path.join(dirname, filename.split(self.params.input.reflections_suffix)[0] +
               self.params.input.experiments_suffix)
            if self._is_accepted(experiments_path):
              yield experiments_path, path

  def _is_accepted(self, experiment_path):
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

