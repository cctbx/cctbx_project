from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/12/2014
Last Changed: 01/05/2015
Description : IOTA command-line module. Version 0.7
'''

import os
import sys
import logging
from datetime import datetime

from libtbx.easy_mp import pool_map

import prime.iota.iota_input as inp
import prime.iota.iota_gridsearch as gs
from prime.iota.iota_select import best_file_selection


# Multiprocessor wrapper for grid search module
def index_mproc_wrapper(current_img):
        return gs.index_integrate(current_img, log_dir, gs_params)


# ---------------------------------------------------------------------------- #

if __name__ == "__main__":

  gs_version = '0.7'
  ps_version = '0.7'

  print '\n{}'.format(datetime.now())
  print 'Starting IOTA ... \n\n'

  gs_params, txt_out = inp.process_input(sys.argv[1:])

  if gs_params.input == None:
    print '{:-^100}\n\n'.format('IOTA Dry Run')
    print txt_out
    with open('iota_default.phil', 'w') as default_settings_file:
      default_settings_file.write(txt_out)
  else:
    input_list, input_dir_list, output_dir_list, log_dir = inp.make_lists(gs_params.input, gs_params.output)

    # Check for list of files and extract a) list of files, b) input paths, 
    # and c) output paths that will override earlier lists
    if gs_params.input_list != None:
      input_dir_list = []
      output_dir_list = []
      input_list = []
      with open(gs_params.input_list, 'r') as listfile:
        listfile_contents = listfile.read()
        input_list = listfile_contents.splitlines()

        for source_file in input_list:
          if os.path.relpath(os.path.dirname(source_file), 
                             os.path.abspath(gs_params.input)) == '.':
            input_dir = os.path.abspath(gs_params.input)
            output_dir = os.path.abspath(gs_params.output)
          else:
            input_dir = os.path.dirname(source_file)
            output_dir = os.path.join(os.path.abspath(gs_params.output), 
                                      os.path.relpath(input_dir, 
                                      os.path.abspath(gs_params.input)))

          if input_dir not in input_dir_list: 
            input_dir_list.append(input_dir)
          if output_dir not in output_dir_list: 
            output_dir_list.append(output_dir)


    # ------------------ Grid Search ------------------

    # Check for grid search toggle, only do it turned on
    if gs_params.grid_search.flag_on == True:

      inp.make_dirs(output_dir_list, log_dir)

      # Setup grid search logger
      gs_logger = logging.getLogger('gs_log')
      gs_formatter = logging.Formatter('%(message)s')
      gs_fileHandler = logging.FileHandler(log_dir + '/grid_search.log', mode='w')
      gs_fileHandler.setFormatter(gs_formatter)
      gs_streamHandler = logging.StreamHandler()
      gs_streamHandler.setFormatter(gs_formatter)

      gs_logger.setLevel(logging.INFO)
      gs_logger.addHandler(gs_fileHandler)
      gs_logger.addHandler(gs_streamHandler)

      # Starting info
      f = open(gs_params.target, 'r')
      phil_file_contents = f.read()
      f.closed
      gs_logger.info('{:-^100} \n'.format(' GRID SEARCH AND PICKLE SELECTION '))
      gs_logger.info("\nTarget file {0} contents:\n".format(gs_params.target))
      gs_logger.info(phil_file_contents)

      gs_logger.info('Found image files in the following folder(s):')
      for folder in input_dir_list:
        gs_logger.info(str(os.path.abspath(folder)))
      gs_logger.info('\nSpot-finding parameter grid search: ' \
                    '{0} input files, spot height: {1} - {2}, '\
                    'spot area: {3} - {4} \n'.format(len(input_list),
                    gs_params.grid_search.h_min, gs_params.grid_search.h_max,
                    gs_params.grid_search.a_min, gs_params.grid_search.a_max))
      gs_logger.info('{:-^100} \n\n'.format(' STARTING GRID SEARCH '))

      # run grid search on multiple processes
      if len(input_list) <= gs_params.n_processors:
        num_procs = len(input_list)
      else:
        num_procs = gs_params.n_processors

      pool_map(args=input_list, func=index_mproc_wrapper, processes=num_procs)

      gs_logger.info('\n\nIOTA grid search version {0}'.format(gs_version))


    # ------------------ Pickle Selection ------------------

    # Setup pickle selection logger
    ps_logger = logging.getLogger('ps_log')
    ps_formatter = logging.Formatter('%(message)s')
    ps_fileHandler = logging.FileHandler(log_dir + '/pickle_select.log', mode='w')
    ps_fileHandler.setFormatter(ps_formatter)
    ps_streamHandler = logging.StreamHandler()
    ps_streamHandler.setFormatter(ps_formatter)

    ps_logger.setLevel(logging.INFO)
    ps_logger.addHandler(ps_fileHandler)
    ps_logger.addHandler(ps_streamHandler)


    ps_logger.info('\n\n{:-^100} \n'.format(' PICKLE SELECTION '))
    
    print output_dir_list
    

    for output_dir in output_dir_list:
      ps_logger.info('Found integrated pickles ' \
                      'under {0}'.format(os.path.abspath(output_dir)))

    if gs_params.flag_prefilter == True:
      prefilter = "ON"
    else:
      prefilter = "OFF"

    ps_logger.info('Space group / unit cell prefilter '\
                   'turned {0} \n\n'.format(prefilter))
    ps_logger.info('{:-^100} \n'.format(' STARTING SELECTION '))

    for output_dir in output_dir_list:
         ps_logger.info('Processing integrated pickles in {0}\n'.format(output_dir))
         best_file_selection(gs_params, output_dir, log_dir)


    # This section checks for output and summarizes file integration and
    # selection results
    ps_logger.info('raw images processed:         {}'.format(len(input_list)))
    
    if os.path.isfile('{0}/not_integrated.lst'.format(os.path.abspath(gs_params.output))):
      with open('{0}/not_integrated.lst'.format(os.path.abspath(gs_params.output)), 
                'r') as int_fail_list:
        int_fail_list_contents = int_fail_list.read()
        int_fail_count = len(int_fail_list_contents.splitlines())
        ps_logger.info('raw images not integrated:    {}'.format(int_fail_count))    
    
    if os.path.isfile('{0}/prefilter_fail.lst'.format(os.path.abspath(gs_params.output))):
      with open('{0}/prefilter_fail.lst'.format(os.path.abspath(gs_params.output)), 
              'r') as bad_int_list:
        bad_int_list_contents = bad_int_list.read()
        bad_int_count = len(bad_int_list_contents.splitlines())
        ps_logger.info('images failed prefilter:      {}'.format(bad_int_count))

    if os.path.isfile('{0}/best_by_strong.lst'.format(os.path.abspath(gs_params.output))):
      with open('{0}/best_by_strong.lst'.format(os.path.abspath(gs_params.output)), 
                'r') as output_list:
        output_list_contents = output_list.read()
        final_count = len(output_list_contents.splitlines())    
    else:
      final_count = 0
    ps_logger.info('pickles in final selection:   {}'.format(final_count))
      
    ps_logger.info('\n\nIOTA pickle select version {0}'.format(ps_version))
