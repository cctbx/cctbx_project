from __future__ import absolute_import, division, print_function
from six.moves import range
import sys
import os

b2GB = 1024 * 1024 * 1024 # byte to GB
memory_per_node_GB = 90 # node memory limit in GB; the maximum memory per cori knl node is 96GB, but the calculation here is approximate, so 6 GB is used as a cushion
ranks_per_node = 68 # ranks per cori knl node
pickle_to_memory = 3.5 # an empirical coefficient to convert pickle file size to anticipated run-time process memory on Cori KNL

debug_file_load_calculator = False
debug_log_path = None

def debug_log_write(string):

  if not debug_file_load_calculator:
    return

  debug_log_file_handle = open(debug_log_path, 'a')
  debug_log_file_handle.write(string)
  debug_log_file_handle.close()

from xfel.merging.application.mpi_logger import mpi_logger

class file_load_calculator(object):
  def __init__(self, params, file_list):
    self.params = params
    self.file_list = file_list
    self.mpi_logger = mpi_logger(params)

    global debug_log_path
    if debug_file_load_calculator:
      debug_log_path = os.path.join(self.params.output.output_dir, 'calculate_file_load.out')

  def calculate_file_load(self, available_rank_count=0):
    '''Calculate a load and build a dictionary {rank:file_list} for the input number of ranks'''
    self.mpi_logger.log_step_time("CALCULATE_FILE_LOAD")
    rank_files = {}
    if self.params.input.parallel_file_load.method == "uniform":
      rank_files = self.calculate_file_load_simple(available_rank_count)
    elif self.params.input.parallel_file_load.method == "node_memory":
      rank_files = self.calculate_file_load_node_memory_based(available_rank_count)

    if debug_file_load_calculator:
      for rank in range(len(rank_files)):
        debug_log_write("\nRank %d"%rank)
        for file_pair in rank_files[rank]:
          debug_log_write("\n%s"%str(file_pair))

    total_file_pairs = 0
    for key, value in rank_files.items():
      total_file_pairs += len(value)
    self.mpi_logger.log("Generated a list of %d file items for %d ranks"%(total_file_pairs, len(rank_files)))
    self.mpi_logger.log_step_time("CALCULATE_FILE_LOAD", True)

    return rank_files

  def calculate_file_load_simple(self, available_rank_count):
    '''Uniformly distribute json/pickle file pairs over the input number of ranks. Return a dictionary {rank:filepair_list}'''
    assert available_rank_count > 0, "Available rank count has to be greater than zero."
    rank_files = {} #{rank:[file_pair1, file_pair2, ...]}
    for rank in range(0, available_rank_count):
      rank_files[rank] = self.file_list[rank::available_rank_count]

    return rank_files

  def calculate_file_load_node_memory_based(self, available_rank_count):
    '''Assign json/pickle file pairs to nodes taking into account node memory limit. Then distribute node-assigned file pairs over the ranks within each node. Return a dictionary {rank:file_list}'''
    # get sizes of all files
    file_sizes = {} # {file_pair:file_pair_size_GB}
    for index in range(len(self.file_list)):
      file_sizes[self.file_list[index]] = os.stat(self.file_list[index][1]).st_size / b2GB # [1] means: use only the pickle file size for now

    # assign files to the anticipated nodes - based on the file sizes and the node memory limit
    node_files = {} # {node:[file_pair1, file_pair2,...]}
    node = 0
    node_files[node] = []
    anticipated_memory_usage_GB = 0
    for file_pair in file_sizes:
      anticipated_memory_usage_GB += (file_sizes[file_pair] * pickle_to_memory)
      if anticipated_memory_usage_GB < memory_per_node_GB: # keep appending the files as long as the total anticipated memory doesn't exceed the node allowance
        node_files[node].append(file_pair)
      else:
        node += 1
        node_files[node] = []
        node_files[node].append(file_pair)
        anticipated_memory_usage_GB = (file_sizes[file_pair] * pickle_to_memory)

    # now we know how many nodes are required
    required_number_of_nodes = len(node_files)
    print("\nRequired number of nodes: %d\n"%required_number_of_nodes)

    # for each node evenly distribute the files over the ranks
    rank_files = {} #{rank:[file_pair1, file_pair2, ...]}
    rank_base = 0
    required_number_of_ranks = 0
    for node in range(required_number_of_nodes):
      rank_base = (node * ranks_per_node)
      for rank in range(rank_base, rank_base + ranks_per_node):
        rank_files[rank] = node_files[node][rank - rank_base::ranks_per_node]
        if len(rank_files[rank]) > 0:
          required_number_of_ranks += 1

    if available_rank_count > 0: # if the caller has provided the available rank count, assert that we have enough ranks
      assert required_number_of_ranks <= available_rank_count, "Not enough ranks to load the pickle files: available %d rank(s), required %d rank(s)"%(available_rank_count, required_number_of_ranks)

    # optionally print out the anticipated memory load per node
    if debug_file_load_calculator:
      for node in range(required_number_of_nodes):
        anticipated_memory_usage_GB = 0
        for file_pair in node_files[node]:
          anticipated_memory_usage_GB += os.stat(file_pair[1]).st_size / b2GB * pickle_to_memory
        debug_log_write ("\nNode %d: anticipated memory usage %f GB"%(node, anticipated_memory_usage_GB))

    return rank_files

from xfel.merging.application.phil.phil import Script as Script_Base

class Script(Script_Base):
  '''A class for running the script.'''

  def get_file_list(self):

    # get json/pickle file list
    loader = file_loader(self.params)
    file_list = loader.filename_lister()

    print("Discovered %d input json/pickle file pairs"%len(file_list))

    return file_list

  def run(self, available_rank_count=0):

    # read and parse phil
    self.initialize()
    self.validate()

    load_calculator = file_load_calculator(self.params, self.get_file_list())
    rank_files = load_calculator.calculate_file_load(available_rank_count)

    return

use_rank_count = 0 # change this to the available rank count, when testing this script with the "uniform" file loading method
if __name__ == '__main__':

  script = Script()

  result = script.run(use_rank_count)

  if result is None:
    sys.exit(1)

  debug_log_write("\nOK")
