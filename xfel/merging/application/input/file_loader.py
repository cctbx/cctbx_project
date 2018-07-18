from __future__ import division
import glob, os
from dxtbx.model.experiment_list import ExperimentListFactory
from libtbx import easy_pickle

"""
Utility functions used for reading input data
"""

class file_loader(object):
  def __init__(self, params):
    self.params = params

  # used by rank 0, client/server
  def filepair_generator(self):
    for pathstring in self.params.input.path:
      for path in glob.glob(pathstring):
        if os.path.isdir(path):
          for filename in os.listdir(path):
            if filename.endswith(self.params.input.reflections_suffix):
              experiments_path = os.path.join(path, filename.split(self.params.input.reflections_suffix)[0] +
                 self.params.input.experiments_suffix)
              if not os.path.exists(experiments_path): continue
              yield experiments_path, os.path.join(path, filename)
        else:
          dirname = os.path.dirname(path)
          filename = os.path.basename(path)
          if filename.endswith(self.params.input.reflections_suffix):
            experiments_path = os.path.join(dirname, filename.split(self.params.input.reflections_suffix)[0] +
               self.params.input.experiments_suffix)
            if not os.path.exists(experiments_path): continue
            yield experiments_path, path
          elif filename.endswith(self.params.input.experiments_suffix):
            pickle_path = os.path.join(dirname, filename.split(self.params.input.experiments_suffix)[0] +
               self.params.input.reflections_suffix)
            if not os.path.exists(pickle_path): continue
            yield path, pickle_path

  # used by rank 0, striping
  def filename_lister(self):
    return list(self.filepair_generator())

  # used by worker rank
  def load_data(self, experiments_filename, reflections_filename):
    return ExperimentListFactory.from_json_file(experiments_filename, check_format = False), \
           easy_pickle.load(reflections_filename)
