import os
import sys
import shutil
from libtbx.easy_mp import pool_map
import prime.iota.iota_gridsearch as gs
from prime.iota.iota_select import read_pickle, best_file_selection
from prime.iota.iota_input import process_input, setup_logger, make_lists, make_dirs
import logging
from datetime import datetime, time
import random

# Multiprocessor wrapper for grid search module
def index_mproc_wrapper(current_img):
        return gs.cxi_index(current_img, log_dir, gs_params)


# ------------------------------------------------------------------------------------------------------------- #

if __name__ == "__main__":

        gs_version = '0.1'
        ps_version = '0.1'

        print '\n{}'.format(datetime.now())
        print 'Starting IOTA ... \n\n'

        gs_params, txt_out = process_input(sys.argv[1:])


        if gs_params.input == None:
                print txt_out
        else:
                img_list, input_dir_list, output_dir_list, log_dir = make_lists(gs_params.input, gs_params.output)

                # ------------------ Grid Search ------------------

                # Check for grid search toggle, only do it turned on
                if gs_params.grid_search.flag_on == True:

                        make_dirs(output_dir_list, log_dir)
                        setup_logger('gs_log', log_dir + '/grid_search.log')
                        gs_logger = logging.getLogger('gs_log')

                        # Starting info
                        f = open(gs_params.target, 'r')
                        phil_file_contents = f.read()
                        f.closed
                        gs_logger.info('{:-^100} \n'.format(' GRID SEARCH AND PICKLE SELECTION '))
                        gs_logger.info("\nTarget file {0} contents:\n".format(gs_params.target))
                        gs_logger.info(phil_file_contents)

#                       if gs_params.flag_random == True:
#                               input_list = random.sample(img_list, 5)
#                               gs_logger.info('Selected the following random images for analysis:')
#                               for item in input_list: gs_logger.info(item)
#                               gs_logger.info('\nSpot-finding parameter grid search: {0} input files, spot height: {1} - {2}, spot area: {3} - {4} \n'.format(len(input_list), int(gs_params.grid_search.h_min/2), int(gs_params.grid_search.h_max/2), int(gs_params.grid_search.a_min*2), int(gs_params.grid_search.a_max*2)))
#
#                       else:
#                               input_list = img_list
#                               for folder in input_dir_list: gs_logger.info(str(os.path.abspath(folder)))
#                               gs_logger.info('Found image files in the following folder(s):')
#                               gs_logger.info('\nSpot-finding parameter grid search: {0} input files, spot height: {1} - {2}, spot area: {3} - {4} \n'.format(len(input_list), gs_params.grid_search.h_min, gs_params.grid_search.h_max, gs_params.grid_search.a_min, gs_params.grid_search.a_max))

                        input_list = img_list
                        gs_logger.info('Found image files in the following folder(s):')
                        for folder in input_dir_list: gs_logger.info(str(os.path.abspath(folder)))
                        gs_logger.info('\nSpot-finding parameter grid search: {0} input files, spot height: {1} - {2}, spot area: {3} - {4} \n'.format(len(input_list), gs_params.grid_search.h_min, gs_params.grid_search.h_max, gs_params.grid_search.a_min, gs_params.grid_search.a_max))
                        gs_logger.info('{:-^100} \n\n'.format(' STARTING GRID SEARCH '))

                        # run grid search on multiple processes
                        if len(input_list) <= gs_params.n_processors:
                                num_procs = len(input_list)
                        else:
                                num_procs = gs_params.n_processors
                        pool_map(args=input_list,
                                                func=index_mproc_wrapper,
                                                processes=num_procs)

                        gs_logger.info('\n\nIOTA grid search version {0}'.format(gs_version))


                # ------------------ Pickle Selection ------------------

                setup_logger('ps_log', log_dir + '/pickle_select.log')
                ps_logger = logging.getLogger('ps_log')

                ps_logger.info('\n\n{:-^100} \n'.format(' PICKLE SELECTION '))

                for output_dir in output_dir_list:
                        ps_logger.info('Found integrated pickles under {0}'.format(os.path.abspath(output_dir)))
                ps_logger.info('Selected by most reflections with I / sigI > {0}'.format(gs_params.min_sigma))

                if gs_params.flag_prefilter == True:
                        prefilter = "ON"
                else: prefilter = "OFF"

                ps_logger.info('Space group / unit cell prefilter turned {0} \n\n'.format(prefilter))
                ps_logger.info('{:-^100} \n'.format(' STARTING SELECTION '))

                for output_dir in output_dir_list:
                        ps_logger.info('Processing integrated pickles in {0}\n'.format(output_dir))
                        best_file_selection(gs_params, output_dir, log_dir)

                gs_logger.info('\n\nIOTA pickle select version {0}'.format(ps_version))
