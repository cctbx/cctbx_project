# Indexing and integration of raw images over a grid
import os
import logging
from iota_select import read_pickle
import shutil
from subprocess import check_output
logger = logging.getLogger('gs_log')

# Indexing and integration w/ grid search of spotfinding parameters
def cxi_index(current_img, log_dir, gs_params):

        # generate filenames, etc.
        path = os.path.dirname(current_img)
        img_filename = os.path.split(current_img)[1]

        if os.path.relpath(path, os.path.abspath(gs_params.input)) == '.':  # in case of all input in one dir
                input_dir = os.path.abspath(gs_params.input)
                output_dir = os.path.abspath(gs_params.output)
        else:                                                                                                                                                                                                                                                           # in case of input in tree
                input_dir = os.path.abspath(gs_params.input) + '/' + os.path.relpath(path, os.path.abspath(gs_params.input))
                output_dir = os.path.abspath(gs_params.output) + '/' + os.path.relpath(path, os.path.abspath(gs_params.input))

        current_output_dir = "{0}/tmp_{1}".format(output_dir, img_filename.split('.')[0])
        index_log_dir = "{0}/tmp_{1}".format(log_dir, img_filename.split('.')[0])

        # Make directories for output / log file for the image being integrated
        if not os.path.exists(current_output_dir): os.makedirs(current_output_dir)
        if not os.path.exists(index_log_dir): os.makedirs(index_log_dir)

        # Indexing / integration against the grid of spotfinding parameters
        for sig_height in range(gs_params.grid_search.h_min, gs_params.grid_search.h_max + 1):
                for spot_area in range (gs_params.grid_search.a_min, gs_params.grid_search.a_max + 1):
                        cmds = ["cxi.index",
                                                        "-d",
                                                        "-o" + current_output_dir,
                                                        "-b" + "int_h{0}_a{1}_".format(sig_height, spot_area),
                                                        'target={0}'.format(gs_params.target),
                                                        "distl.minimum_signal_height={0}".format(str(sig_height)),
                                                        "distl.minimum_spot_height={0}".format(str(sig_height)),
                                                        "distl.minimum_spot_area={0}".format(str(spot_area)),
                                                        current_img]

                        output = check_output(cmds)

                        # write integration logfile
                        current_log_file = "{0}/int_h{1}_a{2}_{3}.log".format(index_log_dir, sig_height, spot_area, img_filename.split('.')[0])
                        f = open(current_log_file, 'w')
                        f.write(output)
                        f.close()

                        # generate, display and log grid search output
                        current_file = "{0}/int_h{1}_a{2}_{3}".format(current_output_dir, sig_height, spot_area, img_filename)
                        if os.path.isfile(current_file):
                                pickle_res, sg, unit_cell, num_obs, num_strong_obs = read_pickle(gs_params, current_file)
                                int_status = "integrated with res = {:>6.2f} - {:<5.2f}, s.g.: {:^{wsg}}, u.c.: {:>6.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}".format(pickle_res[0], pickle_res[1], sg, unit_cell[0], unit_cell[1], unit_cell[2], unit_cell[3], unit_cell[4], unit_cell[5], wsg = len(str(sg)))
                        else:
                                int_status = "not integrated"
                        grid_search_output = "{:^{width}}: h = {:<3}, a = {:<3} ---> {}".format(current_img, sig_height, spot_area, int_status, width = len(current_img) + 2)
                        logger.info(grid_search_output)
