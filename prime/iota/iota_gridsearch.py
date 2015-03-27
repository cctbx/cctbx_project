from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/10/2014
Last Changed: 03/06/2015
Description : Runs cxi.index with signal/spot height and area grid search
'''

import os
from prime.iota.iota_input import main_log

import dials.util.command_line as cmd
import iota_vis_integration as viz

import sys
from cStringIO import StringIO


class Capturing(list):
  """ Class used to capture stdout from cctbx.xfel objects. Saves output in
  appendable list for potential logging.
  """
  def __enter__(self):
    self._stdout = sys.stdout
    sys.stdout = self._stringio = StringIO()
    return self
  def __exit__(self, *args):
    self.extend(self._stringio.getvalue().splitlines())
    sys.stdout = self._stdout


# =============================== INTEGRATION ================================ #

def integrate_one_image(mp_entry, n_int, log_dir, save_pickle, gs_params):
  """ Indexes / integrates using spotfinding parameters determined with the
      previous grid search. Outputs a single pickle file and returns a list of
      selected integration results with filename.


  Output:
  1. integrated pickle in common output folder
  2. integration logfile (from captured cctbx.xfel output, one per image)
  3. list entry for integrated, not integrated or incorrectly integrated file

  Parameters:
  1. mp_entry: raw integration image in *.pickle format
  2. current_output_dir: output directory for this image
  3. n_int: number of integrations (for progress bar)
  2. log_dir: main folder for logs
  3. gs_params: PHIL object containing grid-search and other parameters (from
                user input file)
  """

  from xfel.phil_preferences import load_cxi_phil
  from xfel.cxi.display_spots import run_one_index_core

  logfile = '{}/iota.log'.format(log_dir)

  current_img = mp_entry[0]
  sig_height = mp_entry[1]
  spot_height = mp_entry[2]
  spot_area = mp_entry[3]
  current_output_dir = mp_entry[4]
  target = gs_params.target
  int_results = None

# generate filenames, etc.
  path = os.path.dirname(current_img)
  img_filename = os.path.basename(current_img)
  img_no_ext = img_filename.split('.')[0]

  current_file = os.path.normpath("{0}/int_h{1}_a{2}_{3}"\
            "".format(current_output_dir, sig_height, spot_area, img_filename))
  result_file = os.path.normpath("{0}/int_{1}.lst"\
                                "".format(current_output_dir, img_no_ext))
  current_log_file = os.path.normpath("{0}/{1}.log".format(log_dir, img_no_ext))

# Add signal/spot height and spot area arguments, w/ target file
  if save_pickle or gs_params.advanced.save_tmp_pickles:
    arguments = ["target={0}".format(target),
                 "distl.minimum_signal_height={0}".format(str(sig_height)),
                 "distl.minimum_spot_height={0}".format(str(spot_height)),
                 "distl.minimum_spot_area={0}".format(str(spot_area)),
                 "indexing.completeness_pickle={0}".format(current_file)
                 ]
  else:
    arguments = ["target={0}".format(target),
                 "distl.minimum_signal_height={0}".format(str(sig_height)),
                 "distl.minimum_spot_height={0}".format(str(spot_height)),
                 "distl.minimum_spot_area={0}".format(str(spot_area))
                 ]

#Actually run integration (from run_one_index_core)
  with Capturing() as index_log:
    try:
      args = ["indexing.data={}".format(current_img),
            "beam_search_scope=0.5",
            "lepage_max_delta = 3.0",
            "spots_pickle = None",
            "subgroups_pickle = None",
            "refinements_pickle = None",
            "rmsd_tolerance = 5.0",
            "mosflm_rmsd_tolerance = 5.0",
            "difflimit_sigma_cutoff=2.0",
            "indexing.verbose_cv=True",
            "indexing.open_wx_viewer=False"
            ] + list(arguments[1:])

      horizons_phil = load_cxi_phil(target, args)
      info = run_one_index_core(horizons_phil)
      int_final = info.organizer.info['best_integration']['integration']
      int_results = info.organizer.info['best_integration']['integration']['AD14_parameters']

    except Exception, e:
      if hasattr(e, "classname"):
        print e.classname, "for %s:"%file,
      else:
        print "Indexing error for %s:"%file,
      print e

  # Output results of integration (from the "info" object returned by
  # run_one_index_core)
  if int_results != None:
    # **** FOR TRYING THINGS ONLY!! ****
    if gs_params.advanced.debug:
      for item in int_results: print item
      print "\n*****\n"
      for item in int_final.keys(): print item
      return

    init_res = int_final['resolution']
    uc = int_final['cell'].split()
    sg = int_final['spacegroup']
    res = int_final['resolution']
    
    Is = int_final['I_Observations'].data()
    sigmas = int_final['I_Observations'].sigmas()
    I_over_sigI = Is / sigmas
    spots = len(Is)
    strong_spots = len([i for i in I_over_sigI if i >= gs_params.min_sigma])


    num_corr = int_results['N_correction_vectors']
    dom_size = int_results['domain_sz_ang']
    mosaicity = int_results['fw_mos_deg']
    mos_quality = int_results['mosaic_model_area_under_green_curve_sampled']
    pos_rmsd = int_results['rmsd_px']

    # Write results file
    with open(result_file, 'a') as res_file:
      res_file.write('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},'\
                     '{13}\n'.format(img_filename, spot_height, spot_area,
                                     sg, uc[0], uc[1], uc[2], uc[3], uc[4],
                                     uc[5], strong_spots, res,
                                     mosaicity, mos_quality))

    cell = "{:>6.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}"\
           "".format(float(uc[0]), float(uc[1]), float(uc[2]),
                     float(uc[3]), float(uc[4]), float(uc[5]))
    int_status = 'RES: {:>4.2f} SG: {:^{wsg}}  CELL: {}'\
                 ''.format(init_res, sg, cell, wsg = len(str(sg)))
  else:
    int_status = "not integrated"

  grid_search_output = '{:^{width}}: H = {:<3}, ' \
                       'A = {:<3} ---> {}'.format(img_filename,
                        sig_height, spot_area, int_status,
                        width = len(img_filename) + 2)

# write integration logfile
  with open(current_log_file, 'a') as index_logfile:
    index_logfile.write("{:-^100}\n".format(" INTEGRATION: "\
                        "H={:>2}, A={:>2} ".format(spot_height, spot_area)))
    for item in index_log:
      index_logfile.write("{}\n".format(item))

    index_logfile.write("\n{}\n\n".format(grid_search_output))

  main_log(logfile, grid_search_output)

# Progress bar for integration
  with (open ('{0}/logs/progress.log'.format(gs_params.output), 'a')) as prog_log:
    prog_log.write("{0}\n".format(grid_search_output))

  with (open ('{0}/logs/progress.log'.format(gs_params.output), 'r')) as prog_log:
    prog_content = prog_log.read()
    prog_count = len(prog_content.splitlines())

  gs_prog = cmd.ProgressBar(title='GRID SEARCH', estimate_time=False, spinner=False)
  if prog_count >= n_int:
    gs_prog.finished()
  else:
    prog_step = 100 / n_int
    gs_prog.update(prog_count * prog_step)



# ============================================================================ #
# ============================ FINAL INTEGRATION ============================= #

def final_integrate_one(mp_entry, log_dir, n_int, gs_params):
  """ Performs grid search of spotfinding parameters. Indexes / integrates w/o
  outputting a pickle file, but a list of spotfinding parameters with relevant
  integration results.

  Output:
  1. integrated pickle in correct output folder
  2. integration logfile (from captured cctbx.xfel output, one per image)
  3. list entry for integrated, not integrated or incorrectly integrated file

  Parameters:
  1. mp_entry: raw integration image in *.pickle format
  2. log_dir: main folder for logs
  3. gs_params: PHIL object containing grid-search and other parameters (from
                user input file)
  """


  logfile = '{}/iota.log'.format(log_dir)

  current_img = mp_entry[0]
  sig_height = mp_entry[1]
  spot_height = mp_entry[2]
  spot_area = mp_entry[3]
  target = gs_params.target

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

  current_output_dir = "{0}/tmp_{1}".format(output_dir, img_filename.split('.')[0])
  current_file = os.path.normpath("{0}/int_h{1}_a{2}_{3}".format(current_output_dir,
                               sig_height, spot_area, img_filename))
  current_log_file = os.path.normpath("{0}/{1}.log".format(log_dir,
                                      img_filename.split('.')[0]))
  current_error_file = os.path.normpath("{0}/err_{1}.log".format(log_dir,
                                      img_filename.split('.')[0]))


  from xfel.phil_preferences import load_cxi_phil
  from xfel.cxi.display_spots import run_one_index_core

  # Generate additional parameters for chart output
  # for all charts
  if gs_params.advanced.charts:
    advanced_args = ["integration.enable_residual_map=True",
                     "integration.enable_residual_map_deltapsi=True",
                     "integration.enable_residual_scatter=True",
                     "integration.mosaic.enable_AD14F7B=True",
                     "integration.graphics_backend=pdf",
                     "integration.pdf_output_dir={0}/pdf_h{1}_a{2}".format(current_output_dir, spot_height, spot_area)
                     ]
    if not os.path.exists("{0}/pdf_h{1}_a{2}".format(current_output_dir, spot_height, spot_area)):
      os.makedirs("{0}/pdf_h{1}_a{2}".format(current_output_dir, spot_height, spot_area))

  # for mosaicity chart only
  elif gs_params.advanced.mosaicity_plot:
    advanced_args = ["integration.enable_residual_map=False",
                     "integration.enable_residual_map_deltapsi=True",
                     "integration.enable_residual_scatter=False",
                     "integration.mosaic.enable_AD14F7B=True",
                     "integration.graphics_backend=pdf",
                     "integration.pdf_output_dir={0}/pdf_h{1}_a{2}".format(current_output_dir, spot_height, spot_area)
                     ]
    if not os.path.exists("{0}/pdf_h{1}_a{2}".format(current_output_dir, spot_height, spot_area)):
      os.makedirs("{0}/pdf_h{1}_a{2}".format(current_output_dir, spot_height, spot_area))
  else:
    advanced_args = []

  # Add signal/spot height and spot area arguments, w/ target file
  # Assemble full set of arguments from everything
  arguments = ["target={0}".format(target),
           "distl.minimum_signal_height={0}".format(str(sig_height)),
           "distl.minimum_spot_height={0}".format(str(spot_height)),
           "distl.minimum_spot_area={0}".format(str(spot_area)),
           "indexing.completeness_pickle={0}".format(current_file)
           ] + list(advanced_args[1:])

  int_results = None

  #Actually run integration (from run_one_index_core)
  with Capturing() as index_log:
    try:
      display = False
      args = ["indexing.data={}".format(current_img),
            "beam_search_scope=0.5",
            "lepage_max_delta = 3.0",
            "spots_pickle = None",
            "subgroups_pickle = None",
            "refinements_pickle = None",
            "rmsd_tolerance = 5.0",
            "mosflm_rmsd_tolerance = 5.0",
            "difflimit_sigma_cutoff=2.0",
            "indexing.verbose_cv=True",
            "indexing.open_wx_viewer=%s"%display
            ] + list(arguments[1:])

      horizons_phil = load_cxi_phil(target, args)
      info = run_one_index_core(horizons_phil)
      int_final = info.organizer.info['best_integration']['integration']
      int_results = info.organizer.info['best_integration']['integration']['AD14_parameters']

      with open('{}/integrated.lst'.format(os.path.abspath(gs_params.output)), \
                                                            'a') as final_int:
        final_int.write('{}\n'.format(current_file))

    except Exception, e:
      if hasattr(e, "classname"):
        print e.classname, "for %s:"%file,
      else:
        print "Indexing error for %s:"%file,
      print e

  # Output results of integration (from the "info" object returned by
  # run_one_index_core)
  if int_results != None:
    num_corr = int_results['N_correction_vectors']
    dom_size = int_results['domain_sz_ang']
    mosaicity = int_results['fw_mos_deg']
    mos_quality = int_results['mosaic_model_area_under_green_curve_sampled']
    pos_rmsd = int_results['rmsd_px']
    resol = int_final['resolution']

    cell = int_final['cell'].split()
    uc = [float(cell[0]), float(cell[1]), float(cell[2]),
          float(cell[3]), float(cell[4]), float(cell[5])]

    Is = int_final['I_Observations'].data()
    sigmas = int_final['I_Observations'].sigmas()
    I_over_sigI = Is / sigmas
    spots = len(Is)
    strong_spots = len([i for i in I_over_sigI if i >= gs_params.min_sigma])

    results = [img_filename, spot_height, spot_area, resol, strong_spots, 
              dom_size, mosaicity, pos_rmsd, uc]

    int_status = 'RES: {:>4.2f} SPOTS: {:>5} DOM: {:>8.2f} '\
                 'MOS: {:>6.4f} MQ:{:>6.2f} RMSD: {:>4.2f}'\
                 ''.format(resol, strong_spots, dom_size, mosaicity, 
                           mos_quality, pos_rmsd)

    if gs_params.advanced.pred_img.flag == True:
      viz.make_png(current_img, current_file)
      if gs_params.advanced.pred_img.cv_vectors == True:
        viz.cv_png(current_img, current_file)

  else:
    int_status = "not integrated"
    results = []

  grid_search_output = '{:^{width}}: h = {:<3}, ' \
                       'a = {:<3} ---> {}'.format(img_filename,
                        sig_height, spot_area, int_status,
                        width = len(img_filename) + 2)

# write integration logfile
  with open(current_log_file, 'a') as index_logfile:
    index_logfile.write("\n\n*****\n{:-^100}\n".format(" FINAL INTEGRATION: "\
                        "H={:>2}, A={:>2} ".format(spot_height, spot_area)))
    for item in index_log:
      index_logfile.write("{}\n".format(item))

    index_logfile.write("\n{}\n\n".format(grid_search_output))

  main_log(logfile, grid_search_output)

  # Progress bar for integration
  with (open ('{0}/logs/progress.log'.format(gs_params.output), 'a')) as prog_log:
    prog_log.write("{0}\n".format(grid_search_output))

  with (open ('{0}/logs/progress.log'.format(gs_params.output), 'r')) as prog_log:
    prog_content = prog_log.read()
    prog_count = len(prog_content.splitlines())

  gs_prog = cmd.ProgressBar(title='INTEGRATING', estimate_time=False, spinner=False)
  if prog_count >= n_int:
    gs_prog.finished()
  else:
    prog_step = 100 / n_int
    gs_prog.update(prog_count * prog_step)

  return results
  
  
def int_single_image(mp_entry, log_dir, n_int, gs_params):
  """ Performs grid search of spotfinding parameters. Set aside for single-image
      mode ONLY. Need to clean this up later such that have only one integration
      module.

  Output:
  1. integrated pickle output folder
  2. integration logfile (from captured cctbx.xfel output, one per image)
  3. list entry for integrated, not integrated or incorrectly integrated file

  Parameters:
  1. mp_entry: raw integration image in *.pickle format
  2. log_dir: main folder for logs
  3. n_int: number of integration operations for progress bar
  4. gs_params: PHIL object containing grid-search and other parameters (from
                user input file)
  """


  logfile = '{}/iota.log'.format(log_dir)

  current_img = mp_entry[0]
  sig_height = mp_entry[1]
  spot_height = mp_entry[2]
  spot_area = mp_entry[3]
  target = gs_params.target

# generate filenames, etc.
  path = os.path.dirname(current_img)
  img_filename = os.path.basename(current_img)
  output_dir = os.path.abspath(gs_params.output)
  current_output_dir = "{0}/tmp_{1}".format(output_dir, img_filename.split('.')[0])
  current_file = os.path.normpath("{0}/int_h{1}_a{2}_{3}".format(current_output_dir,
                               sig_height, spot_area, img_filename))
  current_log_file = os.path.normpath("{0}/{1}.log".format(log_dir,
                                      img_filename.split('.')[0]))
  current_error_file = os.path.normpath("{0}/err_{1}.log".format(log_dir,
                                      img_filename.split('.')[0]))


  from xfel.phil_preferences import load_cxi_phil
  from xfel.cxi.display_spots import run_one_index_core

  # Generate additional parameters for chart output
  # for all charts
  if gs_params.advanced.charts:
    advanced_args = ["integration.enable_residual_map=True",
                     "integration.enable_residual_map_deltapsi=True",
                     "integration.enable_residual_scatter=True",
                     "integration.mosaic.enable_AD14F7B=True",
                     "integration.graphics_backend=pdf",
                     "integration.pdf_output_dir={0}/pdf_h{1}_a{2}".format(current_output_dir, spot_height, spot_area)
                     ]
    if not os.path.exists("{0}/pdf_h{1}_a{2}".format(current_output_dir, spot_height, spot_area)):
      os.makedirs("{0}/pdf_h{1}_a{2}".format(current_output_dir, spot_height, spot_area))

  # for mosaicity chart only
  elif gs_params.advanced.mosaicity_plot:
    advanced_args = ["integration.enable_residual_map=False",
                     "integration.enable_residual_map_deltapsi=True",
                     "integration.enable_residual_scatter=False",
                     "integration.mosaic.enable_AD14F7B=True",
                     "integration.graphics_backend=pdf",
                     "integration.pdf_output_dir={0}/pdf_h{1}_a{2}".format(current_output_dir, spot_height, spot_area)
                     ]
    if not os.path.exists("{0}/pdf_h{1}_a{2}".format(current_output_dir, spot_height, spot_area)):
      os.makedirs("{0}/pdf_h{1}_a{2}".format(current_output_dir, spot_height, spot_area))
  else:
    advanced_args = []

  # Add signal/spot height and spot area arguments, w/ target file
  # Assemble full set of arguments from everything
  arguments = ["target={0}".format(target),
           "distl.minimum_signal_height={0}".format(str(sig_height)),
           "distl.minimum_spot_height={0}".format(str(spot_height)),
           "distl.minimum_spot_area={0}".format(str(spot_area)),
           "indexing.completeness_pickle={0}".format(current_file)
           ] + list(advanced_args[1:])

  int_results = None

  #Actually run integration (from run_one_index_core)
  with Capturing() as index_log:
    try:
      display = False
      args = ["indexing.data={}".format(current_img),
            "beam_search_scope=0.5",
            "lepage_max_delta = 3.0",
            "spots_pickle = None",
            "subgroups_pickle = None",
            "refinements_pickle = None",
            "rmsd_tolerance = 5.0",
            "mosflm_rmsd_tolerance = 5.0",
            "difflimit_sigma_cutoff=2.0",
            "indexing.verbose_cv=True",
            "indexing.open_wx_viewer=%s"%display
            ] + list(arguments[1:])

      horizons_phil = load_cxi_phil(target, args)
      info = run_one_index_core(horizons_phil)
      int_final = info.organizer.info['best_integration']['integration']
      int_results = info.organizer.info['best_integration']['integration']['AD14_parameters']

      with open('{}/integrated.lst'.format(os.path.abspath(gs_params.output)), \
                                                            'a') as final_int:
        final_int.write('{}\n'.format(current_file))

    except Exception, e:
      if hasattr(e, "classname"):
        print e.classname, "for %s:"%file,
      else:
        print "Indexing error for %s:"%file,
      print e

  # Output results of integration (from the "info" object returned by
  # run_one_index_core)
  if int_results != None:
    init_res = int_final['resolution']
    uc = int_final['cell'].split()
    sg = int_final['spacegroup']
    res = int_final['resolution']
    
    Is = int_final['I_Observations'].data()
    sigmas = int_final['I_Observations'].sigmas()
    I_over_sigI = Is / sigmas
    spots = len(Is)
    strong_spots = len([i for i in I_over_sigI if i >= gs_params.min_sigma])

    num_corr = int_results['N_correction_vectors']
    dom_size = int_results['domain_sz_ang']
    mosaicity = int_results['fw_mos_deg']
    mos_quality = int_results['mosaic_model_area_under_green_curve_sampled']
    pos_rmsd = int_results['rmsd_px']

    results = [img_filename, spot_height, spot_area, res, num_corr, dom_size,\
         mosaicity, pos_rmsd, uc]

    cell = "{:>6.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}"\
           "".format(float(uc[0]), float(uc[1]), float(uc[2]),
                     float(uc[3]), float(uc[4]), float(uc[5]))
    int_status = 'RES: {:>4.2f} SG: {:^{wsg}}  CELL: {}'\
                 ''.format(init_res, sg, cell, wsg = len(str(sg)))
# 
# 
#     int_status = 'RES: {:>4.2f} SPOTS: {:>5} DOM: {:>8.2f} '\
#                  'MOS: {:>6.4f} MQ:{:>6.2f} RMSD: {:>4.2f}'\
#                  ''.format(resol, num_spots, dom_size, mosaicity, mos_quality,
#                            pos_rmsd)

    if gs_params.advanced.pred_img.flag == True:
      viz.make_png(current_img, current_file)
      if gs_params.advanced.pred_img.cv_vectors == True:
        viz.cv_png(current_img, current_file)

  else:
    int_status = "not integrated"
    results = []

  grid_search_output = '{:^{width}}: h = {:<3}, ' \
                       'a = {:<3} ---> {}'.format(img_filename,
                        sig_height, spot_area, int_status,
                        width = len(img_filename) + 2)

# write integration logfile
  with open(current_log_file, 'a') as index_logfile:
    index_logfile.write("\n\n*****\n{:-^100}\n".format(" FINAL INTEGRATION: "\
                        "H={:>2}, A={:>2} ".format(spot_height, spot_area)))
    for item in index_log:
      index_logfile.write("{}\n".format(item))

    index_logfile.write("\n{}\n\n".format(grid_search_output))

  print grid_search_output

  return results
  
# ============================================================================ #
# =========================== EXPERIMENTAL SECTION =========================== #

def exp_integrate_one(mp_entry, log_dir, n_int, gs_params):
  """ Performs grid search of spotfinding parameters. Indexes / integrates w/o
  outputting a pickle file, but a list of spotfinding parameters with relevant
  integration results.

  Output:
  1. integrated pickle in correct output folder
  2. integration logfile (from captured cctbx.xfel output, one per image)
  3. list entry for integrated, not integrated or incorrectly integrated file

  Parameters:
  1. mp_entry: raw integration image in *.pickle format
  2. log_dir: main folder for logs
  3. gs_params: PHIL object containing grid-search and other parameters (from
                user input file)
  """


  logfile = '{}/iota.log'.format(log_dir)

  current_img = mp_entry[0]
  sig_height = mp_entry[1]
  spot_height = mp_entry[2]
  spot_area = mp_entry[3]
  target = gs_params.target

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

  current_output_dir = "{0}/tmp_{1}".format(output_dir, img_filename.split('.')[0])
  current_file = os.path.normpath("{0}/int_h{1}_a{2}_{3}".format(current_output_dir,
                               sig_height, spot_area, img_filename))
  current_log_file = os.path.normpath("{0}/{1}.log".format(log_dir,
                                      img_filename.split('.')[0]))
  current_error_file = os.path.normpath("{0}/err_{1}.log".format(log_dir,
                                      img_filename.split('.')[0]))


  from xfel.phil_preferences import load_cxi_phil
  from xfel.cxi.display_spots import run_one_index_core

  # Generate additional parameters for chart output
  # for all charts
  if gs_params.advanced.charts:
    advanced_args = ["integration.enable_residual_map=True",
                     "integration.enable_residual_map_deltapsi=True",
                     "integration.enable_residual_scatter=True",
                     "integration.mosaic.enable_AD14F7B=True",
                     "integration.graphics_backend=pdf",
                     "integration.pdf_output_dir={0}/pdf_h{1}_a{2}".format(current_output_dir, spot_height, spot_area)
                     ]
    if not os.path.exists("{0}/pdf_h{1}_a{2}".format(current_output_dir, spot_height, spot_area)):
      os.makedirs("{0}/pdf_h{1}_a{2}".format(current_output_dir, spot_height, spot_area))

  # for mosaicity chart only
  elif gs_params.advanced.mosaicity_plot:
    advanced_args = ["integration.enable_residual_map=False",
                     "integration.enable_residual_map_deltapsi=True",
                     "integration.enable_residual_scatter=False",
                     "integration.mosaic.enable_AD14F7B=True",
                     "integration.graphics_backend=pdf",
                     "integration.pdf_output_dir={0}/pdf_h{1}_a{2}".format(current_output_dir, spot_height, spot_area)
                     ]
    if not os.path.exists("{0}/pdf_h{1}_a{2}".format(current_output_dir, spot_height, spot_area)):
      os.makedirs("{0}/pdf_h{1}_a{2}".format(current_output_dir, spot_height, spot_area))
  else:
    advanced_args = []

  # Add signal/spot height and spot area arguments, w/ target file
  # Assemble full set of arguments from everything
  arguments = ["target={0}".format(target),
           "distl.minimum_signal_height={0}".format(str(sig_height)),
           "distl.minimum_spot_height={0}".format(str(spot_height)),
           "distl.minimum_spot_area={0}".format(str(spot_area)),
           "indexing.completeness_pickle={0}".format(current_file)
           ] + list(advanced_args[1:])

  int_results = None

  #Actually run integration (from run_one_index_core)
  with Capturing() as index_log:
    try:
      display = False
      args = ["indexing.data={}".format(current_img),
            "beam_search_scope=0.5",
            "lepage_max_delta = 3.0",
            "spots_pickle = None",
            "subgroups_pickle = None",
            "refinements_pickle = None",
            "rmsd_tolerance = 5.0",
            "mosflm_rmsd_tolerance = 5.0",
            "difflimit_sigma_cutoff=2.0",
            "indexing.verbose_cv=True",
            "indexing.open_wx_viewer=%s"%display
            ] + list(arguments[1:])

      horizons_phil = load_cxi_phil(target, args)
      info = run_one_index_core(horizons_phil)
      int_final = info.organizer.info['best_integration']['integration']
      int_results = info.organizer.info['best_integration']['integration']['AD14_parameters']

      
      with open('{}/integrated.lst'.format(os.path.abspath(gs_params.output)), \
                                                            'a') as final_int:
        final_int.write('{}\n'.format(current_file))

    except Exception, e:
      if hasattr(e, "classname"):
        print e.classname, "for %s:"%file,
      else:
        print "Indexing error for %s:"%file,
      print e

  # Output results of integration (from the "info" object returned by
  # run_one_index_core)
  if int_results != None:
    num_corr = int_results['N_correction_vectors']
    dom_size = int_results['domain_sz_ang']
    mosaicity = int_results['fw_mos_deg']
    mos_quality = int_results['mosaic_model_area_under_green_curve_sampled']
    pos_rmsd = int_results['rmsd_px']
    resol = int_final['resolution']
    num_spots = len(int_final['all_good_spots'])
    cell = int_final['cell'].split()
    uc = [float(cell[0]), float(cell[1]), float(cell[2]),
          float(cell[3]), float(cell[4]), float(cell[5])]

    for i in sorted(int_final.keys()): print i
    #observations = [a.get_obs(int_final['spacegroup']) for a in int_final['results']]
    #print len(observations)

    results_one_image =  info.organizer.info['best_integration']['integration']['results'][0]
    spacegroup =  info.organizer.info['best_integration']['integration']['spacegroup']
    observations = results_one_image.get_obs(spacegroup)
    number_of_reflections = len(observations)
    
    print "results_one_image = {}".format(results_one_image)
    print "spacegroup = {}".format(spacegroup)
    print "observations = {}".format(observations)
    print "number_of_reflections = {}".format(number_of_reflections)

    results = [img_filename, spot_height, spot_area, resol, num_corr, dom_size,\
         mosaicity, pos_rmsd, uc]

    int_status = 'RES: {:>4.2f} SPOTS: {:>5} DOM: {:>8.2f} '\
                 'MOS: {:>6.4f} MQ:{:>6.2f} RMSD: {:>4.2f}'\
                 ''.format(resol, num_spots, dom_size, mosaicity, mos_quality,
                           pos_rmsd)

    if gs_params.advanced.pred_img.flag == True:
      viz.make_png(current_img, current_file)
      if gs_params.advanced.pred_img.cv_vectors == True:
        viz.cv_png(current_img, current_file)

  else:
    int_status = "not integrated"
    results = []

  grid_search_output = '{:^{width}}: h = {:<3}, ' \
                       'a = {:<3} ---> {}'.format(img_filename,
                        sig_height, spot_area, int_status,
                        width = len(img_filename) + 2)

# write integration logfile
  with open(current_log_file, 'a') as index_logfile:
    index_logfile.write("\n\n*****\n{:-^100}\n".format(" FINAL INTEGRATION: "\
                        "H={:>2}, A={:>2} ".format(spot_height, spot_area)))
    for item in index_log:
      index_logfile.write("{}\n".format(item))

    index_logfile.write("\n{}\n\n".format(grid_search_output))

  main_log(logfile, grid_search_output)

  return results
