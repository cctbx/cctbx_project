from __future__ import absolute_import, division, print_function
from six.moves import range
import sys
import os

b2GB = 1024 * 1024 * 1024 # byte to GB
memory_per_node_GB = 90 # node memory limit; it's really 96, but let's have a 6 GB cushion for now
ranks_per_node = 68 # ranks per node
pickle_to_memory = 4 # a coefficient to convert pickle size to anticipated run-time memory

debug_file_load_calculator = False
log_path = None
def log_write(self, string):
  log_file_handle = open(log_path, 'a')
  log_file_handle.write(string)
  log_file_handle.close()

class file_load_calculator(object):
  def __init__(self, params, file_list):
    self.params = params
    self.file_list = file_list

  # Create a dictionary {rank:file_list} for the input number of ranks
  def calculate_file_load(self):
    if self.params.file_load == simple:
      return self.calculate_file_load_simple()
    else:
      return self.calculate_file_load_smart()

  # distribute files over input number of ranks. Return a dictionary {rank:file_list}
  def calculate_file_load_simple(self):

    for rank in range(0, self.params.rank_count):
        rank_files[rank] = self.file_list[rank::self.rank_count]

    return rank_files

  # assign files to nodes taking into account node memory limit. Distribute assigned files over the ranks within each node. Return a dictionary {rank:file_list}
  def calculate_file_load_smart(self):

    # get sizes of all files
    file_sizes = {} # {file_pair:file_pair_size_GB}
    for index in range(len(file_list)):
      file_sizes[file_list[index]] = os.stat(file_list[index][1]).st_size / b2GB # [1] means: use only the pickle file size for now

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
        anticipated_memory_usage_GB = 0
        node_files[node] = []
        node_files[node].append(file_pair)

    # now we know how many nodes are required
    required_number_of_nodes = len(node_files)
    if debug_file_load_calculator:
      print("\nRequired number of nodes: %d"%required_number_of_nodes)

    # for each node evenly distribute the files over the ranks
    rank_files = {} #{rank:[file_pair1, file_pair2, ...]}
    rank_base = 0
    for node in range(required_number_of_nodes):
      rank_base = (node * ranks_per_node)
      for rank in range(rank_base, rank_base + ranks_per_node):
        rank_files[rank] = node_files[node][rank - rank_base::ranks_per_node]

    # optionally print out the anticipated memory load per node
    if debug_file_load_calculator:
      for node in range(required_number_of_nodes):
        anticipated_memory_usage_GB = 0
        for file_pair in node_files[node]:
          anticipated_memory_usage_GB += os.stat(file_pair[1]).st_size / b2GB * pickle_to_memory
        self.log_write ("\nNode %d: anticipated memory usage %f GB"%(node, anticipated_memory_usage_GB))

    return rank_files

from xfel.merging.application.phil.phil import Script as Script_Base
from xfel.merging.application.input.file_loader import file_loader

class Script(Script_Base):
  '''A class for running the script.'''

  def get_file_list(self):

    # get json/pickle file list
    loader = file_loader(self.params)
    file_list = loader.filename_lister()

    if debug_file_load_calculator:
      self.log_write("\nDiscovered %d input json/pickle file pairs"%len(file_list))

    print("\nDiscovered %d input json/pickle file pairs"%len(file_list))

    return file_list


  def run(self):

    debug_file_load_calculator = True
    log_path = os.path(self.params.output.output_dir, 'calculate_load.out')

    # read and parse phil
    self.initialize()
    self.validate()

    load_calculator = file_load_calculator(self.params)

    rank_files = load_calculator.calculate_file_load()

    if debug_file_load_calculator:
      for rank in range(len(rank_files)):
        self.log_write("\nRank %d"%rank)
        for file_pair in rank_files[rank]:
          self.log_write("\n%s"%str(file_pair))

    return

if __name__ == '__main__':

  script = Script()

  result = script.run()

  if result is None:
    sys.exit(1)

  self.log_write ("\nOK")
