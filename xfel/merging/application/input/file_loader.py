from __future__ import division
import glob, os
from dxtbx.model.experiment_list import ExperimentListFactory
from libtbx import easy_pickle

"""
Utility functions used for reading input data
"""

class file_lister(object):
  """ Class for matching jsons to reflection table pickles """
  def __init__(self, params):
    self.params = params

  def filepair_generator(self):
    """ Used by rank 0, client/server to yield a list of json/reflection table pairs """
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

  def filename_lister(self):
    """ Used by rank 0, striping """
    return list(self.filepair_generator())

from xfel.merging.application.worker import worker
class simple_file_loader(worker):
  '''A class for running the script.'''

  def get_list(self):
    """ Read the list of json/reflection table pickle pairs """
    lister = file_lister(self.params)

    file_list = list(lister.filepair_generator())

    return file_list

  def run(self, all_experiments, all_reflections):
    """ Load all the data using MPI """
    from dxtbx.model.experiment_list import ExperimentList
    from dials.array_family import flex
    from libtbx.mpi4py import MPI
    comm = MPI.COMM_WORLD

    rank = comm.Get_rank()
    size = comm.Get_size()

    # Both must be none or not none
    test = [all_experiments is None, all_reflections is None].count(True)
    assert test in [0,2]
    if test == 2:
      all_experiments = ExperimentList()
      all_reflections = flex.reflection_table()
      starting_expts_count = starting_refls_count = 0
    else:
      starting_expts_count = len(all_experiments)
      starting_refls_count = len(all_reflections)

    # Send the list of file paths to each worker
    if rank == 0:
      file_list = self.get_list()
      self.params.input.path = None # the input is already parsed
      print ('Transmitting file list of length %d'%(len(file_list)))
      transmitted = file_list
    else:
      transmitted = None

    transmitted = comm.bcast(transmitted, root = 0)

    # Read the data
    new_file_list = transmitted[rank::size]
    for experiments_filename, reflections_filename in new_file_list:
      experiments = ExperimentListFactory.from_json_file(experiments_filename, check_format = False)
      reflections = easy_pickle.load(reflections_filename)
      for experiment_id, experiment in enumerate(experiments):
        refls = reflections.select(reflections['id'] == experiment_id)
        all_experiments.append(experiment)
        refls['id'] = flex.int(len(refls), len(all_experiments)-1)
        all_reflections.extend(refls)

    print ('Rank %d has read %d experiments consisting of %d reflections'%(rank,
      len(all_experiments)-starting_expts_count, len(all_reflections)-starting_refls_count))

    return all_experiments, all_reflections
