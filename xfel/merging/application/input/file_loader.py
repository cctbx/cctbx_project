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

  # used by rank 0, striping
  def filename_lister(self):
    return list(self.filepair_generator())

  # used by worker rank
  def load_data(self, experiments_filename, reflections_filename):

    experiments = ExperimentListFactory.from_json_file(experiments_filename, check_format = False)
    reflections = easy_pickle.load(reflections_filename)

    return experiments, reflections


from xfel.merging.application.phil.phil import Script as Script_Base

class Script(Script_Base):
  '''A class for running the script.'''

  def load_data(self):
    from dxtbx.model.experiment_list import ExperimentList
    from dials.array_family import flex
    all_experiments = ExperimentList()
    all_reflections = flex.reflection_table()

    # example showing what reading all the data into a single experiment list/
    # reflection table would look like
    loader = file_loader(self.params)
    for experiments_filename, reflections_filename in loader.filepair_generator():
      experiments, reflections = loader.load_data(experiments_filename, reflections_filename)
      sel = (reflections['id'] < 0) | (reflections['id'] >= len(experiments))
      assert sel.count(True) == 0, "Unindexed or invalid reflections found"
      reflections['id'] += len(all_experiments)
      all_reflections.extend(reflections)
      all_experiments.extend(experiments)

    return all_experiments, all_reflections

  def run(self):
    print('''Mock run, merge some data.''')
    self.initialize()

    self.validate()

    experiments, reflections = self.load_data()
    print ('Read %d experiments consisting of %d reflections'%(len(experiments), len(reflections)))

    experiments, reflections = self.modify(experiments, reflections)

    return

if __name__ == '__main__':
  script = Script()
  result = script.run()
  print ("OK")
