from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 04/07/2015
Last Changed: 05/19/2015
Description : Analyzes integration results and outputs them in an accessible
              format. Includes warnings in case of detected non-isomorphism
              and other anomalies that require a more careful processing
'''
import os
import numpy as np
from collections import Counter

import prime.iota.iota_input as inp

def sort_uc(clean_results):
  """Sorts to entries to compile a table of point groups with consensus unit
     cells encountered in the selection. Should work regardless of whether the
     unit cell prefilter was activated.

     input: sel_clean - list of integrated results
     output: uc_table - entries with point group and unit cell info
  """

  uc_table = []
  sgs = [results['sg'] for results in clean_results]

  ranked_sgs = Counter(sgs).most_common()
  cons_sg = ranked_sgs[0][0]
  reduced_sg_list = [i[0] for i in ranked_sgs]

  for sg in reduced_sg_list:
    batch = [e for e in clean_results if e['sg'] == sg]

    avg_a = np.mean([results['a'] for results in batch])
    std_a = np.std([results['a'] for results in batch])
    avg_b = np.mean([results['b'] for results in batch])
    std_b = np.std([results['b'] for results in batch])
    avg_c = np.mean([results['c'] for results in batch])
    std_c = np.std([results['c'] for results in batch])
    avg_alpha = np.mean([results['alpha'] for results in batch])
    std_alpha = np.std([results['alpha'] for results in batch])
    avg_beta = np.mean([results['beta'] for results in batch])
    std_beta = np.std([results['beta'] for results in batch])
    avg_gamma = np.mean([results['gamma'] for results in batch])
    std_gamma = np.std([results['gamma'] for results in batch])

    if sg == cons_sg:
      tag = '===>'
    else:
      tag = ''

    uc_table.append("{:>4} {:<6} {:^4}:  {:<6.2f} ({:>4.2f}), {:<6.2f} ({:>4.2f}), "\
                     "{:<6.2f} ({:>4.2f}), {:<6.2f} ({:>4.2f}), "\
                     "{:<6.2f} ({:>4.2f}), {:<6.2f} ({:>4.2f})"
                     "".format(tag, '({})'.format(len(batch)), sg, avg_a, std_a,
                               avg_b, std_b, avg_c, std_c, avg_alpha, std_alpha,
                               avg_beta, std_beta, avg_gamma, std_gamma))
    uc_check = std_a >= avg_a * 0.1 or std_b >= avg_b * 0.1 or \
               std_c >= avg_c * 0.1 or std_alpha >= avg_alpha * 0.1 or \
               std_beta >= avg_beta * 0.1 or std_gamma >= avg_gamma * 0.1

    if uc_check and sg == cons_sg:
        morphology_warning = True
    else:
        morphology_warning = False

  return uc_table, morphology_warning


def print_results(clean_results, gs_range, logfile):
  """ Prints diagnostics from the final integration run.

      input: clean_results - list of integrated pickles w/ integration data
             gs_range - range of the grid search

  """
  images = [results['img'] for results in clean_results]
  spot_heights = [results['sph'] for results in clean_results]
  sig_heights = [results['sih'] for results in clean_results]
  spot_areas = [results['spa'] for results in clean_results]
  resolutions = [results['res'] for results in clean_results]
  num_spots = [results['strong'] for results in clean_results]
  mosaicities = [results['mos'] for results in clean_results]

  cons_s = Counter(spot_heights).most_common(1)[0][0]
  cons_h = Counter(spot_heights).most_common(1)[0][0]
  cons_a = Counter(spot_areas).most_common(1)[0][0]

  final_table = []
  final_table.append("\n\n{:-^80}\n".format('ANALYSIS OF RESULTS'))
  final_table.append("Total images:          {}".format(len(images)))
  final_table.append("Avg. signal height:    {:<8.3f}  std. dev:    {:<6.2f}"\
                     "  max: {:<3}  min: {:<3}  consensus: {:<3}"\
                     "".format(np.mean(sig_heights), np.std(sig_heights),
                               max(sig_heights), min(sig_heights), cons_s))
  final_table.append("Avg. spot height:      {:<8.3f}  std. dev:    {:<6.2f}"\
                     "  max: {:<3}  min: {:<3}  consensus: {:<3}"\
                     "".format(np.mean(spot_heights), np.std(spot_heights),
                               max(spot_heights), min(spot_heights), cons_h))
  final_table.append("Avg. spot areas:       {:<8.3f}  std. dev:    {:<6.2f}"\
                    "  max: {:<3}  min: {:<3}  consensus: {:<3}"\
                    "".format(np.mean(spot_areas), np.std(spot_areas),
                              max(spot_areas), min(spot_areas), cons_a))
  final_table.append("Avg. resolution:       {:<8.3f}  std. dev:    {:<6.2f}"\
                    "".format(np.mean(resolutions), np.std(resolutions)))
  final_table.append("Avg. number of spots:  {:<8.3f}  std. dev:    {:<6.2f}"\
                    "".format(np.mean(num_spots), np.std(num_spots)))

  final_table.append("Avg. mosaicity:        {:<8.3f}  std. dev:    {:<6.2f}"\
                    "".format(np.mean(mosaicities), np.std(mosaicities)))
  final_table.append("Found unit cells:\n")

  uc_table, morphology_warning = sort_uc(clean_results)

  for item in uc_table:
    final_table.append(item)

  if morphology_warning:
    final_table.append('\nWARNING: the consensus unit cell is too variable! '\
                       'This may be the result \nof multiple crystal forms '\
                       'present in the dataset. Please run cxi.cluster \nto '\
                       'test this more thoroughly.')

  for item in final_table:
      print item
      inp.main_log(logfile, item)

def print_summary(gs_params, n_img, logfile, iota_version, now):
  """ Prints summary by reading contents of files listing
      a) images not integrated
      b) images that failed unit cell filter
      c) total images input
      d) final images successfully processed

      Appends summary to general log file. Also outputs some of it on stdout.

      input: gs_params - parameters from *.param file in PHIL format
  """

  summary = []
  int_fail_count = 0
  bad_int_count = 0
  final_count = 0

  print "\n\n{:-^80}\n".format('SUMMARY')
  inp.main_log(logfile, "\n\n{:-^80}\n".format('SUMMARY'))

  summary.append('raw images processed:                {}'.format(n_img))

  if os.path.isfile('{0}/not_integrated.lst'.format(os.path.abspath(gs_params.output))):
    with open('{0}/not_integrated.lst'.format(os.path.abspath(gs_params.output)),
              'r') as int_fail_list:
      int_fail_list_contents = int_fail_list.read()
      int_fail_count = len(int_fail_list_contents.splitlines())
    summary.append('raw images not integrated:           {}'.format(int_fail_count))

  if os.path.isfile('{0}/prefilter_fail.lst'.format(os.path.abspath(gs_params.output))):
    with open('{0}/prefilter_fail.lst'.format(os.path.abspath(gs_params.output)),
            'r') as bad_int_list:
      bad_int_list_contents = bad_int_list.read()
      bad_int_count = len(bad_int_list_contents.splitlines())
    summary.append('images failed prefilter:             {}'.format(bad_int_count))

  if os.path.isfile('{0}/gs_selected.lst'.format(os.path.abspath(gs_params.output))):
    with open('{0}/gs_selected.lst'.format(os.path.abspath(gs_params.output)),
              'r') as sel_list:
      sel_list_contents = sel_list.read()
      sel_gs_count = len(sel_list_contents.splitlines())
    summary.append('images in grid search selection:     {}'.format(sel_gs_count))

  if os.path.isfile('{0}/integrated.lst'.format(os.path.abspath(gs_params.output))):
    with open('{0}/integrated.lst'.format(os.path.abspath(gs_params.output)),
              'r') as final_list:
      final_list_contents = final_list.read()
      final_count = len(final_list_contents.splitlines())
    summary.append('final integrated pickles:            {}'.format(final_count))

  summary.append('\n\nIOTA version {0}'.format(iota_version))

  for item in summary:
    print item
    inp.main_log(logfile, "{}".format(item))

  inp.main_log(logfile, now)
