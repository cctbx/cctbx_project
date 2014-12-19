from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/10/2014
Last Changed: 12/15/2014
Description : Runs cxi.index with signal/spot height and area grid search
'''

import os
import logging

from xfel.clustering.singleframe import SingleFrame
from subprocess import check_output


gs_logger = logging.getLogger('gs_log')

# Indexing and integration w/ grid search of spotfinding parameters
def index_integrate(current_img, log_dir, gs_params):

  # generate filenames, etc.
  path = os.path.dirname(current_img)
  img_filename = os.path.basename(current_img)

  if os.path.relpath(path, os.path.abspath(gs_params.input)) == '.':
    input_dir = os.path.abspath(gs_params.input)
    output_dir = os.path.abspath(gs_params.output)
  else:
    input_dir = '{0}/{1}'.format(os.path.abspath(gs_params.input),
                                  os.path.relpath(path,
                                  os.path.abspath(gs_params.input)))
    output_dir = '{0}/{1}'.format(os.path.abspath(gs_params.output),
                                  os.path.relpath(path,
                                  os.path.abspath(gs_params.input)))

  current_output_dir = "{0}/tmp_{1}".format(output_dir,
                                            img_filename.split('.')[0])
  index_log_dir = "{0}/tmp_{1}".format(log_dir, img_filename.split('.')[0])

  # Make directories for output / log file for the image being integrated
  if not os.path.exists(current_output_dir): os.makedirs(current_output_dir)
  if not os.path.exists(index_log_dir): os.makedirs(index_log_dir)

  # Indexing / integration against the grid of spotfinding parameters
  for sig_height in range(gs_params.grid_search.h_min,
                          gs_params.grid_search.h_max + 1):
    for spot_area in range (gs_params.grid_search.a_min,
                            gs_params.grid_search.a_max + 1):
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
      current_log_file = "{0}/int_h{1}_a{2}_{3}.log".format(index_log_dir,
                        sig_height, spot_area, img_filename.split('.')[0])
      f = open(current_log_file, 'w')
      f.write(output)
      f.close()

      # generate, display and log grid search output
      current_file = "{0}/int_h{1}_a{2}_{3}".format(current_output_dir,
                                                    sig_height, spot_area,
                                                    img_filename)
      if os.path.isfile(current_file):
        observations = SingleFrame(current_file,
                                  os.path.split(current_file)[1]).miller_array
        pickle_res = observations.d_max_min()
        pg = observations.space_group_info()
        unit_cell = observations.unit_cell().parameters()
        int_status = 'integrated with res = {:>6.2f} - {:<5.2f}, ' \
                      's.g.: {:^{wsg}}, u.c.: {:>6.2f}, {:>6.2f}, {:>6.2f}, ' \
                      '{:>6.2f}, {:>6.2f}, {:>6.2f}'.format(pickle_res[0],
                      pickle_res[1], pg, unit_cell[0], unit_cell[1],
                      unit_cell[2], unit_cell[3], unit_cell[4], unit_cell[5],
                      wsg = len(str(pg)))
      else:
        int_status = "not integrated"
      grid_search_output = '{:^{width}}: h = {:<3}, ' \
                              'a = {:<3} ---> {}'.format(current_img,
                              sig_height, spot_area, int_status,
                              width = len(current_img) + 2)
      gs_logger.info(grid_search_output)
