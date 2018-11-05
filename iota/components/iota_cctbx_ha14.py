from __future__ import division, print_function, absolute_import

'''
Author      : Lyubimov, A.Y.
Created     : 10/10/2014
Last Changed: 11/05/2018
Description : Runs cctbx.xfel integration module either in grid-search or final
              integration mode. Has options to output diagnostic visualizations.
              Includes selector class for best integration result selection
'''

import os
import uuid
import numpy as np

try:  # for Py3 compatibility
    import itertools.izip as zip
except ImportError:
    pass

from spotfinder.array_family import flex

import iota.components.iota_utils as util
from libtbx import easy_pickle, easy_run

class Empty:
  def __init__(self):
    pass

class Triage(object):
  """ Currently only runs a single DISTL instance with default parameters and accepts or
      rejects an image based on number of spots found. In the works: a crude, wide, sparse
      grid search to establish starting spotfinding parameters.
  """

  def __init__(self,
               img,
               params):

    self.img = img
    self.params = params

  def run_distl(self, params):
    """ Performs a quick DISTL spotfinding and returns Bragg spots information.
    """
    from spotfinder.applications import signal_strength
    # run DISTL spotfinder
    try:
      with util.Capturing() as distl_output:
        Org = signal_strength.run_signal_strength(params)
    except NotImplementedError as e:
      print ("NOT IMPLEMENTED ERROR FOR {}".format(self.img))

    # Extract relevant spotfinding info
    for frame in Org.S.images.keys():
      saturation = Org.Files.imageindex(frame).saturation
      Bragg_spots = [flex.sum(spot.wts) for spot in Org.S.images[frame]['inlier_spots']]

    return Bragg_spots

  def triage_image(self):
    """ Performs a quick DISTL spotfinding without grid search.
    """
    from spotfinder.command_line.signal_strength import master_params as sf_params

    sf_params = sf_params.extract()
    sf_params.distl.image = self.img

    E = Empty()
    E.argv=['Empty']
    E.argv.append(sf_params.distl.image)

    log_info = ['{}\n'.format(self.img)]
    img_filename = os.path.basename(self.img)

    # Perform spotfinding
    # ... using spotfinding grid search
    if self.params.cctbx_ha14.image_triage.type == 'grid_search':
      log_info.append('\n CCTBX TRIAGE grid search:')
      # Determine grid search extent
      a_min = self.params.cctbx_ha14.image_triage.grid_search.area_min
      a_max = self.params.cctbx_ha14.image_triage.grid_search.area_max
      a_step = (a_max - a_min) // self.params.cctbx_ha14.image_triage.grid_search.step_size
      h_min = self.params.cctbx_ha14.image_triage.grid_search.height_min
      h_max = self.params.cctbx_ha14.image_triage.grid_search.height_max
      h_step = (h_max - h_min) // self.params.cctbx_ha14.image_triage.grid_search.step_size

      # Cycle through grid points
      spotlist = []
      for spa in range(a_min, a_max + 1, a_step):
        for sph in range(h_min, h_max + 1, h_step):
          sf_params.distl.minimum_spot_area = spa
          sf_params.distl.minimum_spot_height = sph
          sf_params.distl.minimum_signal_height = sph

          # Perform spotfinding
          Bragg_spots = self.run_distl(sf_params)
          N_Bragg_spots = len(Bragg_spots)

          if N_Bragg_spots > 0:
            total_intensity = flex.sum(flex.double(Bragg_spots))
          else:
            total_intensity = 0

          spotlist.append({'bragg':N_Bragg_spots, 'ti':total_intensity,
                           'spa':spa, 'sph':sph})
          log_info.append('{:<{w}}: H = {:<2}, A = {:<2}, Bragg = {:<6.0f}  '\
                          'total intensity = {:<12.4f}'.format(img_filename, sph, spa,
                          N_Bragg_spots, total_intensity, w = len(img_filename)))

      # Pick best spotfinding result (highest total intensity seems to work for now
      pick = sorted(spotlist, key = lambda j: j['ti'])[-1]
      N_Bragg_spots = pick['bragg']
      start_sph = pick['sph']
      start_sih = pick['sph']
      start_spa = pick['spa']

    # ... using spotfinding without grid search
    else:
      # Set spotfinding params
      sf_params.distl.minimum_spot_area = self.params.cctbx_ha14.grid_search.area_median
      sf_params.distl.minimum_spot_height = self.params.cctbx_ha14.grid_search.height_median
      sf_params.distl.minimum_signal_height = self.params.cctbx_ha14.grid_search.height_median

      # Perform spotfinding
      Bragg_spots = self.run_distl(sf_params)

      # Extract spotfinding results
      N_Bragg_spots = len(Bragg_spots)
      start_sph = self.params.cctbx_ha14.grid_search.height_median
      start_sih = self.params.cctbx_ha14.grid_search.height_median
      start_spa = self.params.cctbx_ha14.grid_search.area_median

    # Determine triage success
    if N_Bragg_spots >= self.params.cctbx_ha14.image_triage.min_Bragg_peaks:
      log_info.append('ACCEPTED! Selected starting point:')
      log_info.append('{:<{w}}: S = {:<2}, H = {:<2}, A = {:<2}, Bragg = {:<6.0f}'\
                      ''.format(img_filename, start_sih, start_sph, start_spa,
                                N_Bragg_spots, w = len(img_filename)))
      status = None
    else:
      log_info.append('REJECTED! ({} Bragg peaks found)'.format(N_Bragg_spots))
      status = 'failed triage'

    log_entry = "\n".join(log_info)

    return status, log_entry, start_sph, start_spa


class Integrator():
  ''' Replaces img.process() function in old SingleImage object '''
  def __init__(self, init, verbose=True):
    self.init = init
    self.params = init.params
    self.img_object = None
    self.verbose = verbose

  def integrate_cctbx(self, tag, grid_point=0, single_image=False):
    """ Runs integration using the Integrator class """

    # Check to see if the image is suitable for grid search / integration
    if self.img_object.fail is not None:
      self.img_object.grid = []
      self.img_object.final['final'] = None
    else:
      integrator = Processor(params=self.init.params,
                             source_image=self.img_object.img_path,
                             output_image=self.img_object.int_file,
                             viz=self.img_object.viz_path,
                             log=self.img_object.int_log,
                             tag=tag,
                             tmp_base=self.init.tmp_base,
                             gain=self.img_object.gain,
                             single_image=single_image)
      if tag == 'grid search':
        self.img_object.log_info.append('\nCCTBX grid search:')
        for i in range(len(self.img_object.grid)):
          int_results = integrator.integrate(self.img_object.grid[i])
          self.img_object.grid[i].update(int_results)
          img_filename = os.path.basename(self.img_object.conv_img)
          log_entry ='{:<{width}}: S = {:<3} H = {:<3} ' \
                     'A = {:<3} ---> {}'.format(img_filename,
                      self.img_object.grid[i]['sih'],
                      self.img_object.grid[i]['sph'],
                      self.img_object.grid[i]['spa'],
                      self.img_object.grid[i]['info'],
                      width = len(img_filename) + 2)
          self.img_object.log_info.append(log_entry)
          self.img_object.gs_results.append(log_entry)

        # Throw out grid search results that yielded no integration
        self.img_object.grid = [i for i in self.img_object.grid if
                                "not integrated" not in i['info'] and
                                "no data recorded" not in i['info']]
        self.img_object.status = 'grid search'

      elif tag == 'split grid':
        self.img_object.log_info.append('\nCCTBX INTEGRATION grid search:')
        int_results = integrator.integrate(self.img_object.grid[grid_point])
        self.img_object.grid[grid_point].update(int_results)
        img_filename = os.path.basename(self.img_object.conv_img)
        log_entry = '{:<{width}}: S = {:<3} H = {:<3} A = {:<3} ---> {}' \
                    ''.format(img_filename,
                              self.img_object.grid[grid_point]['sih'],
                              self.img_object.grid[grid_point]['sph'],
                              self.img_object.grid[grid_point]['spa'],
                              self.img_object.grid[grid_point]['info'],
                              width=len(img_filename) + 2)
        self.img_object.log_info.append(log_entry)
        self.img_object.gs_results.append(log_entry)

      elif tag == 'integrate':
        self.img_object.log_info.append('\nCCTBX final integration:')
        final_results = integrator.integrate(self.img_object.final)
        self.img_object.final.update(final_results)
        self.img_object.status = 'final'
        img_filename = os.path.basename(self.img_object.conv_img)
        log_entry = '{:<{width}}: S = {:<3} H = {:<3} A = {:<3} ---> {}' \
                    ''.format(img_filename,
                              self.img_object.final['sih'],
                              self.img_object.final['sph'],
                              self.img_object.final['spa'],
                              self.img_object.final['info'],
                              width=len(img_filename) + 2)
        self.img_object.log_info.append(log_entry)

        import iota.components.iota_vis_integration as viz
        if self.params.analysis.viz == 'integration':
          viz.make_png(self.img_object.final['img'],
                       self.img_object.final['final'],
                       self.img_object.viz_file)
        elif self.params.analysis.viz == 'cv_vectors':
          viz.cv_png(self.img_object.final['img'],
                     self.img_object.final['final'],
                     self.img_object.viz_file)

  def select_cctbx(self):
    """ Selects best grid search result using the Selector class """
    if self.img_object.fail is None:
      selector = Selector(self.img_object.grid,
                          self.img_object.final,
                          self.params.cctbx_ha14.selection.prefilter.flag_on,
                          self.params.cctbx_ha14.selection.prefilter.target_uc_tolerance,
                          self.params.cctbx_ha14.selection.prefilter.target_pointgroup,
                          self.params.cctbx_ha14.selection.prefilter.target_unit_cell,
                          self.params.cctbx_ha14.selection.prefilter.min_reflections,
                          self.params.cctbx_ha14.selection.prefilter.min_resolution,
                          self.params.cctbx_ha14.selection.select_by)

      self.img_object.fail, self.img_object.final, log_entry = selector.select()
      self.img_object.status = 'selection'
      self.img_object.log_info.append(log_entry)

  def process(self, img_object, single_image=False):
    """ Image processing; selects method, runs requisite modules """
    self.img_object = img_object

    if self.img_object.status != 'bypass grid search':
      self.img_object.status = 'processing'

    #for CCTBX indexing / integration
    terminate = False
    prev_status = self.img_object.status
    prev_fail = 'first cycle'
    prev_final = self.img_object.final
    prev_epv = 9999

    while not terminate:
      # Run grid search if haven't already
      if (
              self.img_object.fail is None and
              'grid search' not in self.img_object.status
      ):
        self.integrate_cctbx('grid search', single_image=single_image)

      # Run selection if haven't already
      if (
              self.img_object.fail is None and
              self.img_object.status != 'selection'
      ):
        self.select_cctbx()

      # If smart grid search is active run multiple rounds until convergence
      if self.params.cctbx_ha14.grid_search.type == 'smart':
        if (
                self.img_object.fail is None and
                self.img_object.final['epv'] < prev_epv
        ):
          prev_epv = self.img_object.final['epv']
          prev_final = self.img_object.final
          prev_status = self.img_object.status
          prev_fail = self.img_object.fail
          self.hmed = self.img_object.final['sph']
          self.amed = self.img_object.final['spa']
          self.img_object.generate_grid()
          self.img_object.final['final'] = self.img_object.int_file
          if len(self.img_object.grid) == 0:
            self.img_object.final = prev_final
            self.img_object.status = prev_status
            self.img_object.fail = prev_fail
            terminate = True
            continue
          if self.verbose:
            log_entry = '\nNew starting point: H = {}, A = {}\n'\
                        ''.format(self.hmed, self.amed)
            self.img_object.log_info.append(log_entry)
        else:
          if prev_fail != 'first cycle':
            self.img_object.final = prev_final
            self.img_object.status = prev_status
            self.img_object.fail = prev_fail
            if self.verbose:
              log_entry = '\nFinal set of parameters: H = {}, A = {}'\
                          ''.format(self.img_object.final['sph'],
                                    self.img_object.final['spa'])
              self.img_object.log_info.append(log_entry)
          terminate = True

      # If brute force grid search is selected run one round
      else:
        terminate = True

    # Run final integration if haven't already
    if self.img_object.fail is None and self.img_object.status != 'final':
      self.integrate_cctbx('integrate', single_image=single_image)

    # If verbose output selected (default), write to main log
    if self.verbose:
       log_entry = "\n".join(self.img_object.log_info)
       util.main_log(self.init.logfile, log_entry)
       util.main_log(self.init.logfile, '\n\n')

    # Make a temporary process log into a final process log
    if os.path.isfile(self.img_object.int_log):
      l_fname = util.make_filename(self.img_object.int_log, new_ext='log')
      final_int_log = os.path.join(self.img_object.log_path, l_fname)
      os.rename(self.img_object.int_log, final_int_log)

    # Save results into a pickle file
    self.img_object.status = 'final'
    from libtbx import easy_pickle as ep
    ep.dump(self.img_object.obj_file, self.img_object)

    return self.img_object

  def run(self, img_object):
    return self.process(img_object=img_object, single_image=False)


class Processor(object):
  """ Class for image integration (w/ grid search params) """
  def __init__(self,
               params,
               source_image = None,
               output_image = None,
               viz = None,
               log = None,
               tag = 'grid search',
               tmp_base = None,
               gain = 1,
               single_image = False):

    self.params = params
    self.img = source_image
    self.out_img = output_image
    self.min_sigma = self.params.cctbx_ha14.selection.min_sigma
    self.target = os.path.abspath(self.params.cctbx_ha14.target)
    self.viz = viz
    self.tag = tag
    self.int_log = log
    self.charts = self.params.analysis.charts
    self.tmp_base = tmp_base
    self.single_image = single_image
    self.method = self.params.mp.method
    self.queue = self.params.mp.queue

    self.args = ["target={}".format(self.target),
                 "indexing.data={}".format(self.img),
                 "spots_pickle=None",
                 "subgroups_pickle=None",
                 "refinements_pickle=None",
                 "rmsd_tolerance=5.0",
                 "mosflm_rmsd_tolerance=5.0",
                 "integration.detector_gain={}".format(gain),
                 "indexing.verbose_cv=True"]

    # Add target unit cell if exists
    if self.params.cctbx_ha14.target_unit_cell is not None:
      t_uc = [str(i) for i in self.params.cctbx_ha14.target_unit_cell.parameters()]
      self.args.extend(['target_cell="{}"'.format(' '.join(t_uc))])

    # Translate / add target lattice if exists
    t_lat = util.makenone(self.params.cctbx_ha14.target_lattice_type)
    if t_lat:
      if t_lat == 'triclinic':
        known_setting = 1
      elif t_lat == 'monoclinic':
        known_setting = 2
      elif t_lat in ('orthorhombic', 'rhombohedral'):
        known_setting = 5
      elif t_lat == 'tetragonal':
        known_setting = 9
      elif t_lat == 'hexagonal':
        known_setting = 12
      elif t_lat == 'cubic':
        known_setting = 22
      else:
        known_setting = None

      if known_setting:
        self.args.extend(['known_setting={}'.format(known_setting)])

    # Centering type if exists
    t_ctype = self.params.cctbx_ha14.target_centering_type
    if t_ctype is not None:
      self.args.extend(['target_cell_centring_type={}'.format(t_ctype)])

    # Resolution, if exists
    hires = self.params.cctbx_ha14.resolution_limits.high
    lowres = self.params.cctbx_ha14.resolution_limits.low
    if hires is None:
      hires = 1.5
    if lowres is None:
      lowres = 99.9
    self.args.extend(['force_method2_resolution_limit={}'.format(hires),
                      'distl_lowres_limit={}'.format(lowres),
                      'distl_highres_limit={}'.format(hires),
                      'distl.res.inner={}'.format(lowres),
                      'distl.res.outer={}'.format(hires)])

    # Add gain (necessary now)
    self.args.extend(['integration.detector_gain={}'.format(gain)])

    with open(os.path.join(self.params.output, 'test.txt'), 'w') as tf:
      tf.write('\n'.join(self.args))

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
    with util.Capturing() as index_log:
      arguments = ["distl.minimum_signal_height={}".format(str(self.s)),
                   "distl.minimum_spot_height={}".format(str(self.h)),
                   "distl.minimum_spot_area={}".format(str(self.a)),
                   "indexing.open_wx_viewer=False"] + list(args[1:])

      tmppath = os.path.join(self.tmp_base, str(uuid.uuid4()) + ".pickle")
      assert not os.path.exists(tmppath)

      # invoke the indexer in a way that will protect iota from any crashes
      command = "iota.bulletproof {} {} {}" \
                "".format(tmppath, self.target, " ".join(arguments))
      try:
        easy_run.fully_buffered(command,join_stdout_stderr=True).show_stdout()
        if not os.path.exists(tmppath):
          print (tmppath)
          print (command)
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
          print (e.classname, "for %s:"%self.img,)
          error_message = "{}: {}".format(e.classname, e[0].replace('\n',' ')[:50])
        else:
          print ("Integration error for %s:"%self.img,)
          error_message = "{}".format(str(e).replace('\n', ' ')[:50])
        print (e)

    # Output results of integration (from the "info" object returned by
    # run_one_index_core)
    if int_final is None:
      if error_message != '':
        reason_for_failure = " - {}".format(error_message)
      else:
        reason_for_failure = ''
      int_status = 'not integrated' + reason_for_failure
      int_results = {'info': int_status}
    elif int_final['observations'][0] is None:
      int_status = 'no data recorded'
      int_results = {'info': int_status}
    else:
      try:
        obs = int_final['observations'][0]
        cell = obs.unit_cell().parameters()
        sg = int_final['pointgroup']
        lres, hres = obs.d_max_min()

        # Calculate number of spots w/ high I / sigmaI
        Is = obs.data()
        sigmas = obs.sigmas()
        I_over_sigI = Is / sigmas
        #spots = len(Is)
        strong_spots = len([i for i in I_over_sigI if i >= self.min_sigma])

        # Mosaicity parameters
        mosaicity = round((int_final.get('ML_half_mosaicity_deg', [0])[0]), 6)
        dom_size = int_final.get('ML_domain_size_ang', [0])[0]
        ewald_proximal_volume = int_final.get('ewald_proximal_volume', [0])[0]

        # Assemble output for log file and/or integration result file
        p_cell = "{:>6.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}"\
               "".format(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5])

        int_status = 'RES: {:<4.2f}  NSREF: {:<4}  SG: {:<5}  CELL: {}'\
                     ''.format(hres, strong_spots, sg, p_cell)

        int_results = {'sg':sg, 'a':cell[0], 'b':cell[1], 'c':cell[2],
                       'alpha':cell[3], 'beta':cell[4], 'gamma':cell[5],
                       'wavelength':int_final['wavelength'],
                       'distance':int_final['distance'],
                       'beamX': int_final['xbeam'], 'beamY': int_final['ybeam'],
                       'strong':strong_spots, 'res':hres, 'lres':lres,
                       'mos':mosaicity, 'epv':ewald_proximal_volume,
                       'info':int_status,'ok':True}
      except ValueError:
        print (self.img)


    # write integration logfile
    if self.tag == 'integrate':
      util.main_log(self.int_log,
                    "{:-^100}\n{:-^100}\n{:-^100}\n"\
                    "".format("", " FINAL INTEGRATION: ", ""\
                    "S = {:>2}, H ={:>2}, A ={:>2} "\
                    "".format(self.s, self.h, self.a)))
    else:
      util.main_log(self.int_log,
                    "{:-^100}\n".format(" INTEGRATION: "\
                    "S = {:>2}, H ={:>2}, A ={:>2} "\
                    "".format(self.s, self.h, self.a)))
    for item in index_log:
      util.main_log(self.int_log, item)

    util.main_log(self.int_log, "\n[ {:^100} ]\n\n".format(int_status))

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
      if self.uc is not None:
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

      i_fail = i['strong'] <= self.min_ref or (self.min_res is not None and
                                               i['res'] >= self.min_res) or (
                         self.pg is not None and
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
