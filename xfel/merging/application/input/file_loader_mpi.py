from __future__ import division, print_function

from xfel.merging.application.input.file_loader import Script as Script_Base
from xfel.merging.application.input.file_loader import file_loader

class Script(Script_Base):
  '''A class for running the script.'''

  def get_list(self):

    # example showing what reading all the data into a single experiment list/
    # reflection table would look like

    loader = file_loader(self.params)

    file_list = list(loader.filepair_generator())

    return file_list

  def run(self, comm = None):
    print('''Mock run, merge some data.''')

    if comm is None:
      from mpi4py import MPI
      comm = MPI.COMM_WORLD

    rank = comm.Get_rank()
    size = comm.Get_size()

    if(rank==0):
      self.initialize()
      self.validate()
      file_list = self.get_list()
      self.params.input.path = None # the input is already parsed
      print ('Transmitting file list of length %d'%(len(file_list)))

      transmitted = dict(params = self.params, options = self.options, file_list = file_list)

    else:
      transmitted = None

    transmitted = comm.bcast(transmitted, root = 0)

    self.params = transmitted['params']
    self.options = transmitted['options']

    new_file_list = transmitted['file_list'][rank::size]

    #print ('\nRank: %d, of size %d, %d' %(rank, size, transmitted['params'].filter.n_obs.min))
    print (new_file_list[0])

    experiments, reflections = self.load_data(new_file_list)

    print ('Rank %d has read %d experiments consisting of %d reflections'%(rank, len(experiments), len(reflections)))

    experiments, reflections = self.modify(experiments, reflections)

    return experiments, reflections

  def load_data(self, file_list):
    from dxtbx.model.experiment_list import ExperimentList
    from dials.array_family import flex
    all_experiments = ExperimentList()
    all_reflections = flex.reflection_table()

    # example showing what reading all the data into a single experiment list/
    # reflection table would look like

    loader = file_loader(self.params)

    for experiments_filename, reflections_filename in file_list:
      experiments, reflections = loader.load_data(experiments_filename, reflections_filename)
      for experiment_id, experiment in enumerate(experiments):
        all_experiments.append(experiment)
        refls = reflections.select(reflections['id'] == experiment_id)
        refls['id'] = flex.int(len(refls), len(all_experiments)-1)
        all_reflections.extend(refls)

    return all_experiments, all_reflections

if __name__ == '__main__':

  from mpi4py import MPI

  comm = MPI.COMM_WORLD


  script = Script()
  result = script.run(comm=comm)
  print ("OK")
