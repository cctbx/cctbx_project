from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/10/2014
Last Changed: 05/19/2015
Description : Runs cctbx.xfel integration module either in grid-search or final
              integration mode. Has options to output diagnostic visualizations
'''

import os
from prime.iota.iota_input import main_log

import iota_vis_integration as viz
import csv

import sys
from cStringIO import StringIO


class Capturing(list):
  """ Class used to capture stdout from cctbx.xfel objects. Saves output in
  appendable list for potential logging.
  """
  def __enter__(self):
    self._stdout = sys.stdout
    self._stderr = sys.stderr
    sys.stdout = self._stringio_stdout = StringIO()
    sys.stderr = self._stringio_stderr = StringIO()
    return self
  def __exit__(self, *args):
    self.extend(self._stringio_stdout.getvalue().splitlines())
    sys.stdout = self._stdout
    self.extend(self._stringio_stderr.getvalue().splitlines())
    sys.stderr = self._stderr


def organize_parameters(int_type, mp_entry, gs_params):
  """ Generates a list of additional arguments for integration per user input.

      input: int_type - grid search, mosaicity scan or final
             mp_entry - list of parameters for integration
             gs_params - user-defined parameters in PHIL format


      output: advanced_args - list of additional arguments for integration
  """

  current_img = mp_entry[0]
  current_output_dir = mp_entry[1]
  img_filename = os.path.basename(current_img)

  sig_height = mp_entry[2]
  spot_height = mp_entry[3]
  spot_area = mp_entry[4]

  current_file = os.path.normpath("{}/int_s{}_h{}_a{}_{}"\
                                  "".format(current_output_dir, sig_height,
                                     spot_height, spot_area, img_filename))

  # Generate additional parameters for chart output
  # for all charts
  if gs_params.advanced.charts:
    advanced_args = ["integration.enable_residual_map=True",
                     "integration.enable_residual_map_deltapsi=True",
                     "integration.enable_residual_scatter=True",
                     "integration.mosaic.enable_AD14F7B=True",
                     "integration.graphics_backend=pdf",
                     "integration.pdf_output_dir={}/pdf_s{}_h{}_a{}"\
                     "".format(current_output_dir, sig_height, spot_height,
                               spot_area)
                     ]
    if not os.path.exists("{}/pdf_s{}_h{}_a{}"\
                        "".format(current_output_dir, sig_height,
                                  spot_height, spot_area)):
      os.makedirs("{}/pdf_s{}_h{}_a{}"
                  "".format(current_output_dir, sig_height,
                            spot_height, spot_area))

  # for mosaicity chart only
  elif gs_params.advanced.mosaicity_plot:
    advanced_args = ["integration.enable_residual_map=False",
                     "integration.enable_residual_map_deltapsi=True",
                     "integration.enable_residual_scatter=False",
                     "integration.mosaic.enable_AD14F7B=True",
                     "integration.graphics_backend=pdf",
                     "integration.pdf_output_dir={}/pdf_s{}_h{}_a{}"\
                     "".format(current_output_dir, sig_height, spot_height,
                               spot_area)
                     ]
    if not os.path.exists("{}/pdf_s{}_h{}_a{}"\
                        "".format(current_output_dir, sig_height,
                                  spot_height, spot_area)):
      os.makedirs("{}/pdf_s{}_h{}_a{}"
                  "".format(current_output_dir, sig_height,
                            spot_height, spot_area))
  else:
    advanced_args = []

  if int_type == "grid":
    if gs_params.advanced.save_tmp_pickles:
      arguments = ["target={0}".format(gs_params.target),
                   "distl.minimum_signal_height={0}".format(str(sig_height)),
                   "distl.minimum_spot_height={0}".format(str(spot_height)),
                   "distl.minimum_spot_area={0}".format(str(spot_area)),
                   "indexing.completeness_pickle={0}".format(current_file)
                   ] + list(advanced_args[1:])
    else:
      arguments = ["target={0}".format(gs_params.target),
                   "distl.minimum_signal_height={0}".format(str(sig_height)),
                   "distl.minimum_spot_height={0}".format(str(spot_height)),
                   "distl.minimum_spot_area={0}".format(str(spot_area))
                   ] + list(advanced_args[1:])
  elif int_type == "final":
    arguments = ["target={0}".format(gs_params.target),
                 "distl.minimum_signal_height={0}".format(str(sig_height)),
                 "distl.minimum_spot_height={0}".format(str(spot_height)),
                 "distl.minimum_spot_area={0}".format(str(spot_area)),
                 "indexing.completeness_pickle={0}".format(current_file),
                 ] + list(advanced_args[1:])


  return arguments


def integrate_image(mp_entry, current_log_file, arguments, ptitle, gs_params):
  """ Runs the integration module in cctbx.xfel; used by either grid-search or
      final integration function.

      input: mp_entry - list of parameters for integration
             current_log_file - verbose integration log output by cctbx.xfel
             arguments - list of additional arguments for integration
             ptitle - title of the integration run, for progress bar
             gs_params - general parameters in PHIL format

      output:  int_results - dictionary of integration results
               int_status - status of the integration attempt
  """

  from xfel.phil_preferences import load_cxi_phil
  from xfel.cxi.display_spots import run_one_index_core

  current_img = mp_entry[0]
  sig_height = mp_entry[2]
  spot_height = mp_entry[3]
  spot_area = mp_entry[4]
  int_final = None
  int_results = {}

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

      horizons_phil = load_cxi_phil(gs_params.target, args)
      info = run_one_index_core(horizons_phil)
      int_final = info.organizer.info['best_integration']['integration']
      int_AD14 = int_final['AD14_parameters']

    except Exception, e:
      if hasattr(e, "classname"):
        print e.classname, "for %s:"%file,
      else:
        print "Indexing error for %s:"%file,
      print e


  # Output results of integration (from the "info" object returned by
  # run_one_index_core)
  if int_final == None:
    int_status = 'not integrated'
    int_results = {}
  elif int_final['I_Observations'] == None:
    int_status = 'no data recorded'
    int_results = {}
  else:
    try:
      # Unit cell / resolution:
      uc = int_final['cell'].split()
      cell = (float(uc[0]), float(uc[1]), float(uc[2]),
              float(uc[3]), float(uc[4]), float(uc[5]))
      sg = int_final['spacegroup']
      res = round(int_final['I_Observations'].d_min(), 4)

      # Calculate number of spots w/ high I / sigmaI
      Is = int_final['I_Observations'].data()
      sigmas = int_final['I_Observations'].sigmas()
      I_over_sigI = Is / sigmas
      spots = len(Is)
      strong_spots = len([i for i in I_over_sigI if i >= gs_params.min_sigma])

      # Mosaicity parameters
      dom_size = int_AD14['domain_sz_ang']
      mosaicity = round(int_AD14['fw_mos_deg'], 6)
      mos_quality = round(int_AD14['mosaic_model_area_under_green_curve_sampled'], 6)

      # Assemble output for log file and/or integration result file
      int_results = {'img':current_img, 'sih':sig_height, 'sph':spot_height,
                     'spa':spot_area, 'sg':sg, 'a':cell[0], 'b':cell[1],
                     'c':cell[2], 'alpha':cell[3], 'beta':cell[4],
                     'gamma':cell[5], 'strong':strong_spots, 'res':res,
                     'mos':mosaicity, 'mq':mos_quality}

      p_cell = "{:>6.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}"\
             "".format(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5])

      int_status = 'RES: {:>4.2f} SG: {:^{wsg}}  CELL: {}'\
                   ''.format(res, sg, p_cell, wsg = len(str(sg)))
    except ValueError:
      import traceback
      print
      print current_img
      raise Exception("".join(traceback.format_exception(*sys.exc_info())))
      sys.exit()


  # write integration logfile
  with open(current_log_file, 'a') as index_logfile:
    if ptitle == 'INTEGRATING':
        index_logfile.write("{:-^100}\n{:-^100}\n{:-^100}\n"\
                            "".format("", " FINAL INTEGRATION: ", ""\
                            "S = {:>2}, H ={:>2}, A ={:>2} "\
                            "".format(sig_height, spot_height, spot_area)))
    else:
      index_logfile.write("{:-^100}\n".format(" INTEGRATION: "\
                          "S = {:>2}, H ={:>2}, A ={:>2} "\
                          "".format(sig_height, spot_height, spot_area)))
    for item in index_log:
      index_logfile.write("{}\n".format(item))

    index_logfile.write("\n[ {:^100} ]\n\n".format(int_status))

  return int_results, int_status

def integration(int_type, mp_entry, log_dir, gs_params):
  """ Integration unit. Calls on integrate_image(). For grid search and
      mosaicity scan integration, saves the results in a CSV-formatted file for
      the selection step. For final integration, outputs the integration result
      in pickle format

      input: int_type - grid search, mosaicity scan, or final
             mp_entry - list of parameters for integration
               a.   current_img - raw image to integrate
               b-d. signal height (b), spot height (c) and spot area (d)
               e.   current_output_dir - folder to output results
             log_dir - general log folder
             gs-params - general parameters in PHIL format

      output: results file in CSV format and (optional) integration pickle
  """

  logfile = '{}/iota.log'.format(log_dir)

  current_img = mp_entry[0]
  current_output_dir = mp_entry[1]
  img_filename = os.path.basename(current_img)
  img_no_ext = img_filename.split('.')[0]

  sig_height = mp_entry[2]
  spot_height = mp_entry[3]
  spot_area = mp_entry[4]
  target = gs_params.target

  current_file = os.path.normpath("{}/int_s{}_h{}_a{}_{}"\
                                  "".format(current_output_dir, sig_height,
                                     spot_height, spot_area, img_filename))
  current_log_file = os.path.normpath("{0}/{1}.log".format(log_dir, img_no_ext))

  if int_type == "grid":
    prog_label = "GRID SEARCH"
    result_file = os.path.normpath("{0}/int_gs_{1}.lst"\
                                  "".format(current_output_dir, img_no_ext))
  elif int_type == "final":
    prog_label = "INTEGRATING"

  arguments = organize_parameters(int_type, mp_entry, gs_params)

  if int_type == "grid" and gs_params.advanced.debug:
    debug_file = '{}/s{}_h{}_a{}_{}.debug'.format(gs_params.output, sig_height,
                                                  spot_height, spot_area,
                                                  img_no_ext)
    with open(debug_file, 'w') as f:
      f.write('')

  # run integration
  results, int_status = integrate_image(mp_entry, current_log_file, arguments,
                                        prog_label, gs_params)

  if int_type == "grid" and gs_params.advanced.debug:
    debug_file = '{}/s{}_h{}_a{}_{}.debug'.format(gs_params.output, sig_height,
                                                  spot_height, spot_area,
                                                  img_no_ext)
    with open(debug_file, 'a') as f:
      f.write('PROCESSED: {}, S = {}, H = {}, A = {}'\
              ''.format(current_img, sig_height, spot_height, spot_area))

  # output results to log file
  if int_type == "grid":
    if results != {}:
      if os.path.isfile(result_file):
        with open(result_file, 'ab') as res_file:
          fn = ['img', 'sih', 'sph', 'spa', 'sg', 'a', 'b', 'c', 'alpha', 'beta',
                'gamma', 'strong', 'res', 'mos', 'mq']
          writer = csv.DictWriter(res_file, fieldnames=fn)
          writer.writerow(results)
      else:
        with open(result_file, 'wb') as res_file:
          fn = ['img', 'sih', 'sph', 'spa', 'sg', 'a', 'b', 'c', 'alpha', 'beta',
                'gamma', 'strong', 'res', 'mos', 'mq']
          writer = csv.DictWriter(res_file, fieldnames=fn)
          writer.writeheader()
          writer.writerow(results)


  grid_search_output = '{:^{width}}: S = {:<3}, H = {:<3}, ' \
                       'A = {:<3} ---> {}'.format(img_filename,
                        sig_height, spot_height, spot_area, int_status,
                        width = len(img_filename) + 2)

  main_log(logfile, grid_search_output)

  if int_type == "final":
    if int_status != 'not integrated':
      with open('{}/integrated.lst'.format(os.path.abspath(gs_params.output)),\
                                                              'a') as f_int:
          f_int.write('{}\n'.format(current_file))

    if gs_params.advanced.pred_img.flag == True:
      viz.make_png(current_img, current_file)
      if gs_params.advanced.pred_img.cv_vectors == True:
        viz.cv_png(current_img, current_file)

    return results
