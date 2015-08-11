from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/10/2014
Last Changed: 07/29/2015
Description : Runs cctbx.xfel integration module either in grid-search or final
              integration mode. Has options to output diagnostic visualizations.
              Includes selector class for best integration result selection
'''

import os
import sys
import numpy as np

import iota_vis_integration as viz
import prime.iota.iota_misc as misc


class Integrator(object):
  """ Class for image integration (w/ grid search params) """
  def __init__(self,
               source_image = None,
               output_image = None,
               min_sigma = 0,
               target = None,
               charts = False,
               viz = None,
               log = None,
               tag = 'grid search'):

    self.img = source_image
    self.out_img = output_image
    self.min_sigma = min_sigma
    self.target = os.path.abspath(target)
    self.viz = viz
    self.tag = tag
    self.int_log = log
    self.charts = charts


    self.args = ["target={}".format(self.target),
                 "indexing.data={}".format(self.img),
                 "beam_search_scope=0.5",
                 "lepage_max_delta = 3.0",
                 "spots_pickle = None",
                 "subgroups_pickle = None",
                 "refinements_pickle = None",
                 "rmsd_tolerance = 5.0",
                 "mosflm_rmsd_tolerance = 5.0",
                 "difflimit_sigma_cutoff=2.0",
                 "indexing.verbose_cv=True"]

  def integrate(self, grid_point):
    """ Runs the integration module in cctbx.xfel; used by either grid-search or
        final integration function. """
    from xfel.phil_preferences import load_cxi_phil
    from xfel.cxi.display_spots import run_one_index_core

    self.s = grid_point['sih']
    self.h = grid_point['sph']
    self.a = grid_point['spa']

    args = self.args

    # Generate advanced arguments (and PDF subfolder)
    if self.charts and self.tag == 'grid search':
      filename = os.path.basename(self.img).split('.')[0]
      pdf_folder = os.path.join(self.viz, 'pdf_{}/s{}_h{}_a{}'\
                   ''.format(filename, self.s, self.h, self.a))
      if not os.path.exists(pdf_folder):
        os.makedirs(pdf_folder)
      self.args.extend(["integration.enable_residual_map=True",
                        "integration.enable_residual_scatter=True",
                        "integration.mosaic.enable_AD14F7B=True",
                        "integration.graphics_backend=pdf",
                        "integration.pdf_output_dir={}".format(pdf_folder)])
    if self.tag == 'integrate':
      args.append("indexing.completeness_pickle={}".format(self.out_img))

    #Actually run integration (from run_one_index_core)
    error_message = ''
    with misc.Capturing() as index_log:
      try:
        arguments = ["distl.minimum_signal_height={}".format(str(self.s)),
                     "distl.minimum_spot_height={}".format(str(self.h)),
                     "distl.minimum_spot_area={}".format(str(self.a)),
                     "indexing.open_wx_viewer=False"] + list(args[1:])

        horizons_phil = load_cxi_phil(self.target, arguments)
        info = run_one_index_core(horizons_phil)
        int_final = info.organizer.info['best_integration']['integration']
        int_AD14 = int_final['AD14_parameters']

      except Exception, e:
        int_final = None
        if hasattr(e, "classname"):
          print e.classname, "for %s:"%file,
          error_message = "{}: {}".format(e.classname, e[0].replace('\n',' ')[:50])
        else:
          print "Integration error for %s:"%file,
          error_message = "{}".format(str(e).replace('\n', ' ')[:50])
        print e


    # Output results of integration (from the "info" object returned by
    # run_one_index_core)
    if int_final == None:
      if error_message != '':
        reason_for_failure = " - {}".format(error_message)
      else:
        reason_for_failure = ''
      int_status = 'not integrated' + reason_for_failure
      int_results = {'info': int_status}
    elif int_final['I_Observations'] == None:
      int_status = 'no data recorded'
      int_results = {'info': int_status}
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
        strong_spots = len([i for i in I_over_sigI if i >= self.min_sigma])

        # Mosaicity parameters
        dom_size = int_AD14['domain_sz_ang']
        mosaicity = round(int_AD14['fw_mos_deg'], 6)
        mos_quality = round(int_AD14['mosaic_model_area_under_green_curve_sampled'], 6)
        ewald_proximal_volume = round(int_AD14['ewald_proximal_volume'], 6)

        # Assemble output for log file and/or integration result file
        p_cell = "{:>6.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}"\
               "".format(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5])

        int_status = 'RES: {:<4.2f}  NSREF: {:<4}  SG: {:<5}  CELL: {}'\
                     ''.format(res, strong_spots, sg, p_cell)

        int_results = {'sg':sg, 'a':cell[0], 'b':cell[1], 'c':cell[2],
                        'alpha':cell[3], 'beta':cell[4], 'gamma':cell[5],
                        'strong':strong_spots, 'res':res, 'mos':mosaicity,
                        'epv':ewald_proximal_volume, 'info':int_status,
                        'ok':True}

      except ValueError:
        import traceback
        print
        print self.img
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))
        sys.exit()


    # write integration logfile
    if self.tag == 'integration':
        misc.main_log(self.int_log,
                      "{:-^100}\n{:-^100}\n{:-^100}\n"\
                      "".format("", " FINAL INTEGRATION: ", ""\
                      "S = {:>2}, H ={:>2}, A ={:>2} "\
                      "".format(self.s, self.h, self.a)))
    else:
      misc.main_log(self.int_log,
                    "{:-^100}\n".format(" INTEGRATION: "\
                    "S = {:>2}, H ={:>2}, A ={:>2} "\
                    "".format(self.s, self.h, self.a)))
    for item in index_log:
      misc.main_log(self.int_log, item)

    misc.main_log(self.int_log, "\n[ {:^100} ]\n\n".format(int_status))

    return int_results


class Selector(object):
  """ Class for selection of optimal spotfinding parameters from grid search """

  def __init__(self,
               grid,
               final,
               apply_prefilter = False,
               uc_tol = 0,
               pg = None,
               uc = None,
               min_ref = 0,
               min_res = None):

    self.grid = grid
    self.apply_prefilter = apply_prefilter
    self.uc = uc
    self.uc_tol = uc_tol
    self.pg = pg
    self.min_ref = min_ref
    self.min_res = min_res
    self.final = final
    self.best = final


  def prefilter(self):
    """ Unit cell pre-filter. Applies hard space-group constraint and stringent
        unit cell parameter restraints to filter out integration results that
        deviate. Optional step. Unit cell tolerance user-defined. """

    for i in self.grid:
      if self.uc != None:
        user_uc = [prm for prm in uc.parameters()]
        delta_a = abs(i['a'] - user_uc[0])
        delta_b = abs(i['b'] - user_uc[1])
        delta_c = abs(i['c'] - user_uc[2])
        delta_alpha = abs(i['alpha'] - user_uc[3])
        delta_beta = abs(i['beta'] - user_uc[4])
        delta_gamma = abs(i['gamma'] - user_uc[5])
        uc_check = (delta_a <= user_uc[0] * uc_tol and
                    delta_b <= user_uc[1] * uc_tol and
                    delta_c <= user_uc[2] * uc_tol and
                    delta_alpha <= user_uc[3] * uc_tol and
                    delta_beta <= user_uc[4] * uc_tol and
                    delta_gamma <= user_uc[5] * uc_tol)
      else:
        uc_check = True

      i_fail = i['strong'] <= self.min_ref or (self.min_res != None and\
               i['res'] >= self.min_res) or (self.pg != None and\
               self.pg.replace(" ","") != i['sg'].replace(" ","")) or\
               not uc_check

      if i_fail:
        i['ok'] = False
      else:
        i['ok'] = True

    return [j for j in self.grid if j['ok']]


  def select(self):
    """ First round of selection for results from the initial spotfinding grid
        search. Select the 25% with lowest mosaicities, then select for most
        spots. """
    log_entry = []
    if len(self.grid) == 0:
      log_entry = '\nNo integration results for {}\n'.format(self.final['img'])
      self.best['info'] = log_entry
      self.best['final'] = None
    else:
      if self.apply_prefilter:
        acceptable_results = self.prefilter()
      else:
        acceptable_results = self.grid

      if len(acceptable_results) == 0:
        log_entry = '\nAll {0} results from {1} failed prefilter' \
                    '\n'.format(len(self.grid), self.final['img'])
        self.best['info'] = log_entry
        self.best['final'] = None
        self.prefilter = False
      else:
        # Generate log summary of integration results
        log_entry.append('\nSelecting from {0} out '\
                        'of {1} integration results for ' \
                        '{2}:\n'.format(len(acceptable_results),
                        len(self.grid), self.final['img']))
        categories = '{:^4}{:^4}{:^4}{:^9}{:^8}{:^55}{:^12}{:^14}{:^14}'\
                     ''.format('S', 'H', 'A', 'RES', 'SG.',
                               'UNIT CELL', 'SPOTS', 'MOS', 'EPV')
        line = '{:-^4}{:-^4}{:-^4}{:-^9}{:-^8}{:-^55}{:-^16}{:-^14}{:-^14}'\
               ''.format('', '', '', '', '','', '', '', '')
        log_entry.append(categories)
        log_entry.append(line)

        for acc in acceptable_results:
          cell = '{:>8.2f}, {:>8.2f}, {:>8.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}'\
                 ''.format(acc['a'], acc['b'], acc['c'],
                           acc['alpha'], acc['beta'], acc['gamma'])
          info_line = '{:^4}{:^4}{:^4}{:^9.2f}{:^8}{:^55}{:^12}{:^14.8f}{:^14.8f}'\
                      ''.format(acc['sih'], acc['sph'], acc['spa'], acc['res'],
                                acc['sg'], cell, acc['strong'], acc['mos'],
                                acc['epv'])
          log_entry.append(info_line)

        # Perform selection
        sorted_entries = sorted(acceptable_results, key=lambda i: i['mos'])
        subset = [j[1] for j in enumerate(sorted_entries) \
                  if j[0] <= len(sorted_entries) * 0.25]
        sub_spots = [sp['strong'] for sp in subset]
        self.best.update(subset[np.argmax(sub_spots)])

        cell = '{:>8.2f}, {:>8.2f}, {:>8.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}'\
                 ''.format(self.best['a'], self.best['b'], self.best['c'],
                           self.best['alpha'], self.best['beta'], self.best['gamma'])
        info_line = '\n{:^4}{:^4}{:^4}{:^9.2f}{:^8}{:^55}{:^12}{:^14.8f}{:^14.8f}\n'\
                      ''.format(self.best['sih'], self.best['sph'], self.best['spa'],
                                self.best['res'], self.best['sg'], cell, self.best['strong'],
                                self.best['mos'], self.best['epv'])
        log_entry.append(info_line)
        self.best['info'] = "Selection results for {}: {}"\
                            "".format(self.best['img'], info_line)
        log_entry = "\n".join(log_entry)

    return self.best, log_entry
