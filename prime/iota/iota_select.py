import os
import shutil
import os,cPickle as pickle,math
import numpy as np
import logging
from xfel.clustering.cluster import Cluster
from iotbx import reflection_file_reader
from cctbx.array_family import flex
from cctbx import miller

ps_logger = logging.getLogger('ps_log')

# Read integrated pickle file and return res limits, space group, unit cell, number of
# total reflections and number of reflections above specified I / sigI
def read_pickle(gs_params, temp_pickle_file):
	acceptable_pickle = False
	pickle_content = pickle.load(open(temp_pickle_file,"rb"))
	observations = pickle_content["observations"][0]

	# Count "strong" reflections, i.e. w/ I/sigma(I) > high minimum, e.g. 5
	I_over_sigI = observations.data()/ observations.sigmas()
	num_strong_obs = len([val for val in I_over_sigI if val >= gs_params.min_sigma])

	#Get resolution, space group and unit cell parameters
	sg = observations.space_group_info()
	uc_a, uc_b, uc_c, uc_alpha, uc_beta, uc_gamma = observations.unit_cell().parameters()

	return observations.d_max_min(), sg, observations.unit_cell().parameters(), len(observations.data()), num_strong_obs

# Selects only integrated pickles that fit sg / uc parameters specified in the .phil
# file. Also checks that the low-res cutoff is beyond 10A.
def prefilter(gs_params, total_tmp_pickles):

	acceptable_pickles = []
	for tmp_pickle in total_tmp_pickles:
		uc_tol = gs_params.target_uc_tolerance
		user_uc = []
		user_uc.append(gs_params.target_pointgroup)
		for prm in gs_params.target_unit_cell.parameters(): user_uc.append(prm)

		# read pickle info and determine uc differences
		pickle_res, pickle_sg, pickle_uc, pickle_reflections, pickle_strong_reflections = read_pickle(gs_params, tmp_pickle)
		uc_a, uc_b, uc_c, uc_alpha, uc_beta, uc_gamma = pickle_uc
		delta_a = abs(uc_a - user_uc[1])
		delta_b = abs(uc_b - user_uc[2])
		delta_c = abs(uc_c - user_uc[3])

		# Determine if pickle satisfies sg / uc parameters within given tolerance and low 
		# resolution cutoff
		if str(pickle_sg) == user_uc[0]:
			if  delta_a <= user_uc[1] * uc_tol and delta_b <= user_uc[2] * uc_tol and delta_c <= user_uc[3] * uc_tol:
				if pickle_res[0] >= 10:
					acceptable_pickles.append(tmp_pickle)

	return acceptable_pickles

# Select integrated pickle with the most reflections with I / sigI over a specified limit
def best_by_strong(gs_params, acceptable_pickles, dest_dir):
	sel_pickle_list = []
	for pickle in acceptable_pickles:
		sel_pickle_list.append(read_pickle(gs_params, pickle)[4])

	best_file_strong = acceptable_pickles[sel_pickle_list.index(max(sel_pickle_list))]
	destination_file = '{0}/{1}'.format(dest_dir, os.path.split(best_file_strong)[1])
	
	return best_file_strong, destination_file

# Main selection module. Looks through integrated pickles in a specified folder and
# copies the best ones to a file. Outputs a list to log file and marks the selected
# pickle file.
def best_file_selection(gs_params, output_dir, log_dir):

	# make a list of folders with integrated pickles
	tmp_dirs = [tmp_dir for tmp_dir in os.listdir(output_dir) if "tmp" in tmp_dir]
	
	for tmp_dir in tmp_dirs:
		abs_tmp_dir = os.path.join(output_dir, tmp_dir)
		total_tmp_pickles = [os.path.join(abs_tmp_dir, tmp_pickle) for tmp_pickle in os.listdir(abs_tmp_dir) if ".pickle" in tmp_pickle]
		
		# apply prefilter if specified and make a list of acceptable pickles
		if gs_params.flag_prefilter == True:
			acceptable_pickles = prefilter(gs_params, total_tmp_pickles)
		else:
			acceptable_pickles = total_tmp_pickles
		
		# Selection and copying of pickles, output of stats to log file
		if len(acceptable_pickles) == 0:
			ps_logger.info("Discarded all {0} integrated pickles in {1}:\n".format(len(acceptable_pickles), tmp_dir))
		else:
			ps_logger.info("Selecting from {0} out of {1} integrated pickles in {2}:\n".format(len(acceptable_pickles), len(total_tmp_pickles), tmp_dir))
			filename = str(os.path.split(acceptable_pickles[0])[1])
			categories = ' {:^{pwidth}}{:^16}{:^15}{:^45}{:^12}{:^10}'.format('Filename', 'resolution', 's.g.', 'unit cell', 'total', 'strong', pwidth=len(filename)+5)
			line = ' {:-^{pwidth}}{:-^16}{:-^15}{:-^45}{:-^12}{:-^12}'.format('', '', '', '', '', '', pwidth=len(filename)+5)
			ps_logger.info(categories)
			ps_logger.info(line)
		
			# Select best pickle (by number of strong reflections) and copy to new folder
			sel_pickle, dest_pickle = best_by_strong(gs_params, acceptable_pickles, output_dir + '/selected')
			shutil.copyfile(sel_pickle, dest_pickle)		
			
			# Report pickle stats. Mark selected pickle with asterisk for posterity
			for pickle in acceptable_pickles:
				pickle_name = os.path.split(pickle)[1]
				selection_tag = ' '
				if os.path.split(pickle)[1] == os.path.split(sel_pickle)[1]: selection_tag = "*"
				res, sg, uc, ref, sref = read_pickle(gs_params, pickle) 
				info_line = '{} {:<{pwidth}}{:>7.2f} - {:<5.2f}{:^15}{:>6.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}{:^12}{:^10}{}'.format(selection_tag, pickle_name, res[0],res[1], sg, uc[0], uc[1], uc[2], uc[3], uc[4], uc[5], ref, sref, selection_tag, pwidth=len(filename)+5)
				ps_logger.info(info_line)
			
		ps_logger.info('\n')