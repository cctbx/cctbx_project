from __future__ import division
import iotbx.phil
#from libtbx.utils import Usage, Sorry
import sys, os, shutil
import logging

master_phil = iotbx.phil.parse("""
description = Integration optimization and transfer app (IOTA) input file
  .type = str
  .help = Run description (optional).
  .multiple = False
  .optional = True
input = None
	.type = path
	.multiple = False
	.help = Path to folder with raw data in pickle format. Can be a tree w/ subfolders
	.optional = False
output = iota_output
	.type = str
	.help = Base name for output folder
	.optional = False
target = target.phil
	.type = str
	.multiple = False
	.help = Target (.phil) file with integration parameters
	.optional = False
flag_random = False
  .type = bool
  .help = Activate grid search / selection of random sample & analysis of parameters
grid_search
  .help = "Parameters for the grid search."
{
  flag_on = True
      .type = bool
      .help = Set to False to turn post-refinement in this section off.
  a_min = 1
    .type = int
    .help = Minimum spot area.
  a_max = 10
    .type = int
    .help = Maximum spot area.
  h_min = 1
    .type = int
    .help = Minimum spot height.
  h_max = 10
    .type = int
    .help = Maximum spot height.  
}
flag_prefilter = True
  .type = bool
  .help = Activate space group / unit cell pre-filter.
target_unit_cell = 79.4, 79.4, 38.1, 90.0, 90.0, 90.0
  .type = unit_cell
  .help = Target unit-cell parameters are used to discard outlier cells.
target_space_group = P 43 21 2
  .type = str
  .help = Target space group.
target_pointgroup = P 4
  .type = str
  .help = Target point group.
target_uc_tolerance = 0.05
  .type = float
  .help = Maximum allowed unit cell deviation from target
select_by = strong
	.type = str
  .help = Pickle selection method
min_sigma = 5.0
  .type = float
  .help = Minimum sigma for "strong" reflections
n_processors = 32
  .type = int
  .help = No. of processing units
""")

def process_input(input_file_list):

	user_phil = []
	for input_file in input_file_list:
		user_phil.append(iotbx.phil.parse(open(input_file).read()))

	working_phil = master_phil.fetch(sources=user_phil)
	params = working_phil.extract()

	#capture input read out by phil
	from cStringIO import StringIO
	class Capturing(list):
		def __enter__(self):
			self._stdout = sys.stdout
			sys.stdout = self._stringio = StringIO()
			return self
		def __exit__(self, *args):
			self.extend(self._stringio.getvalue().splitlines())
			sys.stdout = self._stdout

	with Capturing() as output:
		working_phil.show()

	txt_out = '{:-^100}\n\n'.format('IOTA Dry Run')
	for one_output in output:
		txt_out += one_output + '\n'

	return params, txt_out

# Read input directory tree (if any) and make lists of input folder, output folder and
# input files for use in everything 
def make_lists (input_dir, output_dir):
	input_list = []
	input_dir_list = []
	output_dir_list = []

	abs_inp_path = os.path.abspath(input_dir)
	abs_out_path = os.path.abspath(output_dir)

	# search for *.pickle files within the tree and record in a list w/
	# full absolute path and filanames
	for root, dirs, files in os.walk(abs_inp_path):
		for filename in files:
			if filename.endswith("pickle"): 
				pickle_file = root + '/' + filename
				input_list.append(pickle_file)
			
	# make lists of input and output directories and files
	for input_entry in input_list:
		path = os.path.dirname(input_entry)

		if os.path.relpath(path, abs_inp_path) == '.':  # in case of all input in one dir
			input_dir = abs_inp_path
			if input_dir not in input_dir_list: input_dir_list.append(input_dir)
			output_dir = abs_out_path
		else:											# in case of input in tree
			input_dir = abs_inp_path + '/' + os.path.relpath(path, abs_inp_path)
			output_dir = abs_out_path + '/' + os.path.relpath(path, abs_inp_path)
			
		if input_dir not in input_dir_list: input_dir_list.append(input_dir)
		if output_dir not in output_dir_list: output_dir_list.append(output_dir)

	return input_list, input_dir_list, output_dir_list, abs_out_path + '/logs'
	
# Setup for log files & display
def setup_logger(logger_name, log_file, level=logging.INFO):
	l = logging.getLogger(logger_name)
	formatter = logging.Formatter('%(message)s')
	fileHandler = logging.FileHandler(log_file, mode='w')
	fileHandler.setFormatter(formatter)
	streamHandler = logging.StreamHandler()
	streamHandler.setFormatter(formatter)

	l.setLevel(level)
	l.addHandler(fileHandler)
	l.addHandler(streamHandler)

# Make output directories preserving the tree structure
def make_dirs (output_dir_list, log_dir):

	# Make output directory structure
	for output_dir in output_dir_list:
		if os.path.exists(output_dir):     
			shutil.rmtree(output_dir)
			os.makedirs(output_dir) 
		else:
			os.makedirs(output_dir)

		# make output directory for selected pickles
		best_output_dir = output_dir + '/selected'
		if os.path.exists(best_output_dir):	
			shutil.rmtree(best_output_dir)
			os.makedirs(best_output_dir) 
		else:
			os.makedirs(best_output_dir)
			
		# make log folder for cxi.index logs (one per integration attempt)
		index_log_dir = output_dir + '/logs'
		if os.path.exists(index_log_dir):
			shutil.rmtree(index_log_dir)
			os.makedirs(index_log_dir) 
		else:
			os.makedirs(index_log_dir)

			
	# make log folder (under main output folder regardless of tree structure)
	if os.path.exists(log_dir):
		shutil.rmtree(log_dir)
		os.makedirs(log_dir) 
	else:
		os.makedirs(log_dir)