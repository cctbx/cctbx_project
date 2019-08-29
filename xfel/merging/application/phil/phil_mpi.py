from __future__ import absolute_import, division, print_function

from xfel.merging.application.phil.phil import Script as Script_Base

class Script(Script_Base):
  '''A class for running the script.'''

  def run(self,comm):

    rank = comm.Get_rank()
    size = comm.Get_size()

    if(rank==0):
      self.initialize()
      self.validate()

      transmitted = dict(params = self.params, options = self.options)

    else:
      transmitted = None

    transmitted = comm.bcast(transmitted, root = 0)

    print ('\nRank: %d, of size %d, %d' %(rank, size, transmitted['params'].filter.n_obs.min))

    # do other stuff
    return

if __name__ == '__main__':

  from mpi4py import MPI

  comm = MPI.COMM_WORLD

  rank = comm.Get_rank()
  size = comm.Get_size()

  script = Script()
  result = script.run(comm=comm)

  print ("OK")
