from __future__ import absolute_import, division, print_function
import os

class mpi_logger(object):
  """A class to facilitate each rank writing to its own log file and (optionally) to a special log file for timing"""

  def __init__(self, params=None):
    from xfel.merging.application.mpi_helper import mpi_helper
    self.mpi_helper = mpi_helper()

    if params:
      self.set_log_file_paths(params)
    else:
      self.main_log_file_path = None
      self.rank_log_file_path = None
      self.timing_file_path = None

    self.log_buffer = self.main_log_buffer = '' # a buffer for handling a run-time case, when the file path isn't known yet, so we need to store the output somewhere

    self.timing_table = dict()

  def set_log_file_paths(self, params):

    if params.output.title:
      title = params.output.title
    else:
      title = 'merge'

    self.main_log_file_path = os.path.join(params.output.output_dir, title + ".log")
    self.rank_log_file_path = os.path.join(params.output.output_dir, 'rank_%06d_%06d.out'%(self.mpi_helper.size, self.mpi_helper.rank))

    if params.output.do_timing:
      self.timing_file_path = os.path.join(params.output.output_dir, 'timing_%06d_%06d.out'%(self.mpi_helper.size, self.mpi_helper.rank))
    else:
      self.timing_file_path = None

  def log(self, message, rank_prepend=True):
    '''Log a rank message'''

    # Prepend the message with the rank index - helpful for grep-ping
    if rank_prepend:
      rank_message = "Rank %d: %s\n"%(self.mpi_helper.rank, message)
    else:
      rank_message = message + "\n"
    self.log_buffer += rank_message

    if self.rank_log_file_path == None:
      return # the file path hasn't been parsed yet, so just keep the output in the buffer

    log_file_handle = open(self.rank_log_file_path, 'a')
    log_file_handle.write(self.log_buffer)
    self.log_buffer = ''
    log_file_handle.close()

  def main_log(self, message):
    '''Log a message to the main log file'''

    self.main_log_buffer += (message + '\n')

    if self.main_log_file_path == None:
      return

    log_file_handle = open(self.main_log_file_path, 'a')
    log_file_handle.write(self.main_log_buffer)
    self.main_log_buffer = ''
    log_file_handle.close()

  def log_step_time(self, step, step_finished=False):
    '''Log elapsed time for an execution step'''

    if not step_finished: # a step has started - cache its start time and return
      if not step in self.timing_table:
        self.timing_table[step] = dict({'single_step':dict({'start':self.mpi_helper.time(), 'elapsed':0.0}),
                                        'cumulative': 0.0
                                       }
                                      )
      else: # if a step is executed repeatedly - re-use the existent step key
        self.timing_table[step]['single_step']['start'] = self.mpi_helper.time()
      return

    # a step has finished - calculate its elapsed and cumulative time (the latter is needed when the step is executed repeatedly)
    #self.mpi_helper.comm.barrier()

    if not step in self.timing_table:
      assert False, "A step has finished, but doesn't have its starting time entry"
      return

    self.timing_table[step]['single_step']['elapsed'] = self.mpi_helper.time() - self.timing_table[step]['single_step']['start']
    self.timing_table[step]['cumulative'] += self.timing_table[step]['single_step']['elapsed']

    # log the elapsed time - single step and cumulative
    if self.timing_file_path == None:
      return

    log_file = open(self.timing_file_path,'a')
    log_file.write("RANK %d %s: %f s %f s\n"%(self.mpi_helper.rank, step, self.timing_table[step]['single_step']['elapsed'], self.timing_table[step]['cumulative']))
    log_file.close()
