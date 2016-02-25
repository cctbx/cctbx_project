from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/10/2014
Last Changed: 02/24/2015
Description : Runs cctbx.xfel integration module either in grid-search or final
              integration mode. Has options to output diagnostic visualizations.
              Includes selector class for best integration result selection
'''

import os
import sys
import uuid
import numpy as np

import iota_vis_integration as viz
import prime.iota.iota_misc as misc
from libtbx import easy_pickle, easy_run

class Empty: pass

class Triage(object):
  """ Currently only runs a single DISTL instance with default parameters and accepts or
      rejects an image based on number of spots found. In the works: a crude, wide, sparse
      grid search to establish starting spotfinding parameters.
  """

  def __init__(self,
               img,
               gain,  # Currently not used by DISTL, how awkward!
               params):
    self.img = img
    self.params = params

  def triage_image(self):
    """ Performs a quick DISTL spotfinding without grid search.
    """

    import spotfinder
    from spotfinder.command_line.signal_strength import master_params as sf_params
    from spotfinder.applications.wrappers import DistlOrganizer

    sf_params = sf_params.extract()
    sf_params.distl.image = self.img

    E = Empty()
    E.argv=['Empty']
    E.argv.append(sf_params.distl.image)

    selected_output = []
    total_output = []
    bragg_spots = []
    spotfinding_log = ['{}\n'.format(self.img)]

    # set spotfinding parameters for DISTL spotfinder
    sf_params.distl.minimum_spot_area = self.params.cctbx.grid_search.area_median
    sf_params.distl.minimum_spot_height = self.params.cctbx.grid_search.height_median
    sf_params.distl.minimum_signal_height = int(self.params.cctbx.grid_search.height_median / 2)

    # run DISTL spotfinder
    with misc.Capturing() as junk_output:
      Org = DistlOrganizer(verbose = False, argument_module=E,
                           phil_params=sf_params)

      Org.printSpots()

    # Extract relevant spotfinding info & make selection
    for frame in Org.S.images.keys():
      Bragg_spots = Org.S.images[frame]['N_spots_inlier']

    if Bragg_spots >= self.params.image_triage.min_Bragg_peaks:
      log_info = 'ACCEPTED! {} good Bragg peaks'.format(Bragg_spots)
      status = None
    else:
      log_info = 'REJECTED!'
      status = 'failed triage'

    return status, log_info


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
               tag = 'grid search',
               tmp_base = None,
               gain = 1,
               single_image = False):

    self.img = source_image
    self.out_img = output_image
    self.min_sigma = min_sigma
    self.target = os.path.abspath(target)
    self.viz = viz
    self.tag = tag
    self.int_log = log
    self.charts = charts
    self.tmp_base = tmp_base
    self.single_image = single_image

    self.args = ["target={}".format(self.target),
                 "indexing.data={}".format(self.img),
                 "beam_search_scope=0.5",
                 "lepage_max_delta=3.0",
                 "spots_pickle=None",
                 "subgroups_pickle=None",
                 "refinements_pickle=None",
                 "rmsd_tolerance=5.0",
                 "mosflm_rmsd_tolerance=5.0",
                 "difflimit_sigma_cutoff=2.0",
                 "integration.detector_gain={}".format(gain),
                 "indexing.verbose_cv=True"]

  def integrate(self, grid_point):
    """ Runs the integration module in cctbx.xfel; used by either grid-search or
        final integration function. """

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

    #Actually run integration using iota.bulletproof
    error_message = ''
    with misc.Capturing() as index_log:
      arguments = ["distl.minimum_signal_height={}".format(str(self.s)),
                   "distl.minimum_spot_height={}".format(str(self.h)),
                   "distl.minimum_spot_area={}".format(str(self.a)),
                   "indexing.open_wx_viewer=False"] + list(args[1:])

      tmppath = os.path.join(self.tmp_base, str(uuid.uuid4()) + ".pickle")
      assert not os.path.exists(tmppath)

      # invoke the indexer in a way that will protect iota from any crashes
      command = "iota.bulletproof %s %s %s"%(tmppath, self.target, " ".join(arguments))

      try:
        easy_run.fully_buffered(command,join_stdout_stderr=True).show_stdout()
        if not os.path.exists(tmppath):
          raise Exception("Indexing failed for an unknown reason")

        # iota.bulletproof saves the needed results from indexing in a tmp file
        result = easy_pickle.load(tmppath)
        os.remove(tmppath)
        if isinstance(result, str):
          raise Exception(result)
        else:
          int_final = result

      except Exception, e:
        int_final = None
        if hasattr(e, "classname"):
          print e.classname, "for %s:"%self.img,
          error_message = "{}: {}".format(e.classname, e[0].replace('\n',' ')[:50])
        else:
          print "Integration error for %s:"%self.img,
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
    elif int_final['observations'][0] == None:
      int_status = 'no data recorded'
      int_results = {'info': int_status}
    else:
      try:
        obs = int_final['observations'][0]
        cell = obs.unit_cell().parameters()
        sg = int_final['pointgroup']
        res = obs.d_min()

        # Calculate number of spots w/ high I / sigmaI
        Is = obs.data()
        sigmas = obs.sigmas()
        I_over_sigI = Is / sigmas
        spots = len(Is)
        strong_spots = len([i for i in I_over_sigI if i >= self.min_sigma])

        # Mosaicity parameters
        mosaicity = round((int_final.get('ML_half_mosaicity_deg', [0])[0]), 6)
        dom_size = int_final.get('ML_domain_size_ang', [0])[0]
        ewald_proximal_volume = int_final.get('ewald_proximal_volume', [0])[0]

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
    if self.tag == 'integrate':
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

    # In single-image mode, write a file with h, k, l, I, sigma
    if self.single_image == True and self.tag == 'integrate':
      hklI_filename = "{}.{}".format(os.path.basename(self.out_img).split('.')[0], 'hkli')
      hklI_file = os.path.join(os.path.dirname(self.out_img), hklI_filename)
      hklI = zip(obs.indices(), obs.data(), obs.sigmas())
      for i in hklI:
        with open(hklI_file, 'a') as f:
          entry = '{},{},{},{},{}'.format(i[0][0], i[0][1], i[0][2], i[1], i[2])
          f.write('{}\n'.format(entry))

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
               min_res = None,
               select_by = 'mosaicity'):

    self.grid = grid
    self.apply_prefilter = apply_prefilter
    self.uc = uc
    self.uc_tol = uc_tol
    self.pg = pg
    self.min_ref = min_ref
    self.min_res = min_res
    self.final = final
    self.best = final
    self.fail = None
    self.select_by = select_by


  def prefilter(self):
    """ Unit cell pre-filter. Applies hard space-group constraint and stringent
        unit cell parameter restraints to filter out integration results that
        deviate. Optional step. Unit cell tolerance user-defined. """

    for i in self.grid:
      if self.uc != None:
        user_uc = [prm for prm in self.uc.parameters()]
        delta_a = abs(i['a'] - user_uc[0])
        delta_b = abs(i['b'] - user_uc[1])
        delta_c = abs(i['c'] - user_uc[2])
        delta_alpha = abs(i['alpha'] - user_uc[3])
        delta_beta = abs(i['beta'] - user_uc[4])
        delta_gamma = abs(i['gamma'] - user_uc[5])
        uc_check = (delta_a <= user_uc[0] * self.uc_tol and
                    delta_b <= user_uc[1] * self.uc_tol and
                    delta_c <= user_uc[2] * self.uc_tol and
                    delta_alpha <= user_uc[3] * self.uc_tol and
                    delta_beta <= user_uc[4] * self.uc_tol and
                    delta_gamma <= user_uc[5] * self.uc_tol)
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
      self.fail = 'failed grid search'
    else:
      if self.apply_prefilter:
        acceptable_results = self.prefilter()
      else:
        acceptable_results = self.grid

      if len(acceptable_results) == 0:
        log_entry = '\nAll {0} results from {1} failed prefilter' \
                    '\n'.format(len(self.grid), self.final['img'])
        self.best['info'] = log_entry
        self.fail = 'failed prefilter'
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
        if self.select_by == 'mosaicity':
          sorted_entries = sorted(acceptable_results, key=lambda i: i['mos'])
        elif self.select_by == 'epv':
          sorted_entries = sorted(acceptable_results, key=lambda i: i['epv'])

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

    return self.fail, self.best, log_entry
