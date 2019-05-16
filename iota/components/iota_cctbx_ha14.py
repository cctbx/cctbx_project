from __future__ import absolute_import, division, print_function
from past.builtins import range
from six.moves import range
from six.moves import zip

'''
Author      : Lyubimov, A.Y.
Created     : 10/10/2014
Last Changed: 03/06/2018
Description : Runs cctbx.xfel integration module either in grid-search or final
              integration mode. Has options to output diagnostic visualizations.
              Includes selector class for best integration result selection
'''

import os
import numpy as np
import math

try:  # for Py3 compatibility
    import itertools.izip as zip
except ImportError:
    pass

import dxtbx
from scitbx.array_family import flex
from xfel.cxi.cspad_ana.cspad_tbx import evt_timestamp, dpack
from xfel.cxi.display_spots import run_one_index_core
from xfel.phil_preferences import load_cxi_phil
from libtbx import easy_pickle as ep

import iota.components.iota_utils as util
from iota.components.iota_base import SingleImageBase, ImageImporterBase

ha14_str = '''
image_import
  .help = Parameters for raw image conversion to pickle format
  .alias = Image Import Options
{
  rename_pickle = None keep_file_structure *auto_filename custom_filename
    .type = choice
    .help = "keep_file_structure" retains the input filenames w/ folder tree
    .help = "auto_filename" is <your_login_name>_<conversion_run>_<#####>
  rename_pickle_prefix = None
    .type = str
    .help = Will only be used in conjunction with "custom_filename"
  convert_only = False
    .type = bool
    .help = Set to True (or use -c option) to convert and exit
  square_mode = None no_modification *pad crop
    .type = choice
    .help = Method to generate square image
  flip_beamXY = False
    .type = bool
    .help = flip beamX and beamY parameters when modifying image
}
cctbx_xfel
  .help = Settings for a 2014 version of cctbx.xfel
  .alias = Processing Options
{
  target = None
    .type = str
    .multiple = False
    .help = Target (.phil) file with integration parameters
  resolution_limits
    .help = Sets several resolution limits in Labelit settings
  {
    low = 50.0
      .type = float
    high = 1.5
      .type = float
  }
  target_lattice_type = *None triclinic monoclinic orthorhombic tetragonal rhombohedral hexagonal cubic
    .type = choice
    .help = Target Bravais lattice type if known
  target_centering_type = *None P C I R F
    .type = choice
    .help = Target lattice centering type if known
  target_unit_cell = None
    .type = unit_cell
    .help = Target unit cell parameters (if known)
  grid_search
    .help = "Parameters for the grid search."
  {
    type = None no_grid_search *brute_force smart
      .type = choice
      .help = Set to None to only use median spotfinding parameters
    area_median = 5
      .type = int
      .help = Median spot area.
    area_range = 2
      .type = int
      .help = Plus/minus range for spot area.
    height_median = 4
      .type = int
      .help = Median spot height.
    height_range = 2
      .type = int
      .help = Plus/minus range for spot height.
    sig_height_search = False
      .type = bool
      .help = Set to true to scan signal height in addition to spot height
  }
  selection
    .help = Parameters for integration result selection
  {
    select_by = *epv mosaicity
      .type = choice
      .help = Use mosaicity or Ewald proximal volume for optimal parameter selection
    prefilter
      .help = Used to throw out integration results that do not fit user-defined unit cell information
    {
      flag_on = False
        .type = bool
        .help = Set to True to activate prefilter
      target_pointgroup = None
        .type = str
        .help = Target point group, e.g. "P4"
      target_unit_cell = None
        .type = unit_cell
        .help = In format of "a, b, c, alpha, beta, gamma", e.g. 79.4, 79.4, 38.1, 90.0, 90.0, 90.0
      target_uc_tolerance = None
        .type = float
        .help = Maximum allowed unit cell deviation from target
      min_reflections = None
        .type = int
        .help = Minimum integrated reflections per image
      min_resolution = None
        .type = float
        .help = Minimum resolution for accepted images
    }
  }
}
analysis
  .alias = Analysis Options
  {
  charts = False
    .type = bool
    .help = Set to true to have cctbx.xfel HA14 output  analysis charts
  }
'''

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
    # Set spotfinding params
    sf_params.distl.minimum_spot_area = self.params.cctbx_xfel.grid_search.area_median
    sf_params.distl.minimum_spot_height = self.params.cctbx_xfel.grid_search.height_median
    sf_params.distl.minimum_signal_height = self.params.cctbx_xfel.grid_search.height_median

    # Perform spotfinding
    Bragg_spots = self.run_distl(sf_params)

    # Extract spotfinding results
    N_Bragg_spots = len(Bragg_spots)
    start_sph = self.params.cctbx_xfel.grid_search.height_median
    start_sih = self.params.cctbx_xfel.grid_search.height_median
    start_spa = self.params.cctbx_xfel.grid_search.area_median

    # Determine triage success
    if N_Bragg_spots >= self.params.image_import.minimum_Bragg_peaks:
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
                          self.params.cctbx_xfel.selection.prefilter.flag_on,
                          self.params.cctbx_xfel.selection.prefilter.target_uc_tolerance,
                          self.params.cctbx_xfel.selection.prefilter.target_pointgroup,
                          self.params.cctbx_xfel.selection.prefilter.target_unit_cell,
                          self.params.cctbx_xfel.selection.prefilter.min_reflections,
                          self.params.cctbx_xfel.selection.prefilter.min_resolution,
                          self.params.cctbx_xfel.selection.select_by)

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
      if self.params.cctbx_xfel.grid_search.type == 'smart':
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
    self.min_sigma = self.params.image_import.min_sigma
    self.target = os.path.abspath(self.params.cctbx_xfel.target)
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
    if self.params.cctbx_xfel.target_unit_cell is not None:
      t_uc = [str(i) for i in self.params.cctbx_xfel.target_unit_cell.parameters()]
      self.args.extend(['target_cell="{}"'.format(' '.join(t_uc))])

    # Translate / add target lattice if exists
    t_lat = util.makenone(self.params.cctbx_xfel.target_lattice_type)
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
    t_ctype = self.params.cctbx_xfel.target_centering_type
    if t_ctype is not None:
      self.args.extend(['target_cell_centring_type={}'.format(t_ctype)])

    # Resolution, if exists
    hires = self.params.cctbx_xfel.resolution_limits.high
    lowres = self.params.cctbx_xfel.resolution_limits.low
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

  def integrate(self, grid_point):
    """ Runs the integration module in cctbx.xfel; used by either grid-search or
        final integration function. """

    self.s = grid_point['sih']
    self.h = grid_point['sph']
    self.a = grid_point['spa']

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
      self.args.append("indexing.completeness_pickle={}".format(self.out_img))

    #Actually run integration
    error_message = ''
    with util.Capturing() as index_log:
      arguments = ["distl.minimum_signal_height={}".format(str(self.s)),
                   "distl.minimum_spot_height={}".format(str(self.h)),
                   "distl.minimum_spot_area={}".format(str(self.a)),
                   "indexing.open_wx_viewer=False"] + list(self.args[1:])
      try:
        horizons_phil = load_cxi_phil(self.target, ' '.join(arguments))
        info = run_one_index_core(horizons_phil)
        int_final = info.last_saved_best
      except Exception as e:
        int_final = None
        if hasattr(e, "classname"):
          print (e.classname, "for {}: ".format(self.img),)
          error_message = "{}: {}".format(e.classname, e[0].replace('\n',' ')[:50])
        else:
          print ("Integration error for {}:".format(self.img))
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


class SingleImage(SingleImageBase):
  ''' SingleImage object for use with the old HA14 processing (with
  conversion and all those goodies) '''
  def __init__(self, imgpath, idx=None):
    SingleImageBase.__init__(self, imgpath=imgpath, idx=idx)
    self.raw_img = imgpath
    self.conv_img = None
    self.gain = 1.0

    # Grid search parameters
    self.grid = []
    self.gs_results = []

    # Append DISTL spot-finding parameters to final info object
    self.final['sih'] = 0       # Signal height
    self.final['sph'] = 0       # Spot height
    self.final['spa'] = 0       # Spot area

  def generate_grid(self, gs_type='brute_force', sig_search=False,
                    hr=0, ar=0, hmed=0, amed=0):
    ''' HA14-specific function to run image triage and generate grid search '''

    if gs_type == 'smart':
      hrange = 1
      arange = 1
    elif gs_type in ('none', 'no_grid_search'):
      hrange = 0
      arange = 0
    else:
      hrange = hr
      arange = ar
    grid_points = []

    h_min = hmed - hrange
    h_max = hmed + hrange
    a_min = amed - arange
    a_max = amed + arange
    h_std = hr
    a_std = ar

    for spot_area in range(a_min, a_max + 1):
      for spot_height in range(h_min, h_max + 1):
        if sig_search:
          if spot_height >= 1 + h_std:
            sigs = range(spot_height - h_std, spot_height + 1)
          elif spot_height < 1 + h_std:
            sigs = range(1, spot_height + 1)
          elif spot_height == 1:
            sigs = [1]
        else:
          sigs = [spot_height]

        for sig_height in sigs:
          if (spot_area, spot_height, sig_height) not in grid_points:
            grid_points.append((spot_area, spot_height, sig_height))
            self.grid.append(
              {'sih': sig_height,
               'sph': spot_height,
               'spa': spot_area,
               'a': 0, 'b': 0, 'c': 0,
               'alpha': 0, 'beta': 0, 'gamma': 0,
               'sg': '', 'strong': 0, 'res': 0, 'lres': 0,
               'mos': 0, 'epv': 0, 'info': '', 'ok': True}
            )


class ImageImporter(ImageImporterBase):
  ''' Image importer for old-skool HA14 processing (outputs an image pickle) '''

  def __init__(self, init):
    ImageImporterBase.__init__(self, init=init)
    self.params = init.params
    self.modify = True
    self.conv_base = util.set_base_dir('converted_pickles',
                                       out_dir=self.params.output)

  def instantiate_image_object(self, filepath, idx=None):
    ''' Instantiate a SingleImage object for current backend
    :param filepath: path to image file
    '''
    self.img_object = SingleImage(imgpath=filepath, idx=idx)

  def image_triage(self):
    # Triage image (i.e. check for usable diffraction, using selected method)
    triage = Triage(self.img_object.conv_img, self.params)
    fail, log_entry, hmed, amed = triage.triage_image()

    if not fail:
      self.img_object.status = 'triaged'
    else:
      self.img_object.status = fail

    self.img_object.fail = fail
    self.img_object.log_info.append(log_entry)

    return hmed, amed

  def load_image_file(self, filepath):
    """ Reads a single image file (e.g. a CBF or MCCD) and returns raw data,
    experiment information, and an image type (e.g. pickle, raw image, or none)
    :param img: path to image file
    :return: data: raw image data and experiment info
             img_type: which type of image this was, or None if not loaded
    """

    error = None
    try:
      with util.Capturing() as junk_output:
        loaded_img = dxtbx.load(filepath)
    except Exception as e:
      error = 'IOTA IMPORTER ERROR: DXTBX failed! {}'.format(e)
      print (e)
      loaded_img = None

    # Extract image information
    if loaded_img is not None:
      raw_data   = loaded_img.get_raw_data()
      detector   = loaded_img.get_detector()[0]
      beam       = loaded_img.get_beam()
      scan       = loaded_img.get_scan()
      distance   = detector.get_distance()
      pixel_size = detector.get_pixel_size()[0]
      overload   = detector.get_trusted_range()[1]
      wavelength = beam.get_wavelength()
      beam_x     = detector.get_beam_centre(beam.get_s0())[0]
      beam_y     = detector.get_beam_centre(beam.get_s0())[1]

      if scan is None:
        timestamp = None
      else:
        msec, sec = math.modf(scan.get_epochs()[0])
        timestamp = evt_timestamp((sec,msec))

      # Assemble datapack
      data = dpack(data=raw_data,
                   distance=distance,
                   pixel_size=pixel_size,
                   wavelength=wavelength,
                   beam_center_x=beam_x,
                   beam_center_y=beam_y,
                   ccd_image_saturation=overload,
                   saturated_value=overload,
                   timestamp=timestamp
                   )

      if scan is not None:
        osc_start, osc_range = scan.get_oscillation()
        if osc_start != osc_range:
          data['OSC_START'] = 0 #osc_start
          data['OSC_RANGE'] = 0 #osc_start
          data['TIME'] = scan.get_exposure_times()[0]

      # Populate image object information
      self.img_object.final['pixel_size'] = pixel_size
      self.img_object.final['img_size'] = (data['SIZE1'], data['SIZE2'])
      self.img_object.final['beamX'] = beam_x
      self.img_object.final['beamY'] = beam_y
      self.img_object.final['gain'] = detector.get_gain()
      self.img_object.final['distance'] = distance
      self.img_object.final['wavelength'] = wavelength

    else:
      data = None

    return data, error

  def apply_mask_from_file(self, data, mask_file, invert=False):
    ''' Read a mask from file (should be a bool array) and set every pixel
    under mask to -2 (i.e. will be ignored by processing)

    :param data: raw image data
    :param mask_file: filepath to the mask
    :param invert: invert boolean array from True to False (need for HA14)
    :return: masked data '''

    img_raw_bytes = data['DATA']
    error = None

    try:
      full_mask = ep.load(mask_file)
      if type(full_mask) == tuple:
        full_mask = full_mask[0]

      if invert:
        mask = full_mask == False
      else:
        mask = full_mask
      if mask is not None:
        img_masked = img_raw_bytes.set_selected(mask, -2)
        data['DATA'] = img_masked

    except Exception as e:
      error =  'IOTA IMPORTER ERROR: Mask failed! {}'.format(e)
      print (error)

    return data, error

  def rename_converted_image(self, filepath):
    # Generate converted image pickle filename
    rename_choice = str(
      self.params.image_import.rename_pickle).lower()
    if rename_choice in ("keep_file_structure", "none"):
      img_dir = util.make_image_path(filepath, self.input_base, self.conv_base)
      img_filename = util.make_filename(filepath, new_ext='pickle')
    else:
      self.input_base = self.conv_base
      if rename_choice == "auto_filename":
        prefix = self.init.user_id
      elif rename_choice == "custom_filename":
        prefix = self.params.image_import.rename_pickle_prefix
      else:
        prefix = 'iota'
      img_dir = self.conv_base
      number = int(os.path.basename(self.conv_base))
      img_filename = "{}_{}_{:05d}.pickle" \
                 "".format(prefix, number, self.img_object.img_index)
    return os.path.abspath(os.path.join(img_dir, img_filename))

  def apply_threshold(self, data):
    """ Identifies beamstop shadow and sets pixels inside to -2 (to be ignored by
        processing software). Clunky now (merely finds all pixels with intensity
        less than 0.4 * image average intensity), to be refined later.
    """
    try:
      img_raw_bytes = data['DATA']
      beamstop = self.params.image_import.beamstop

      img_thresh = int((beamstop*flex.mean(img_raw_bytes.as_double())))
      beam_stop_sel = img_raw_bytes <= img_thresh
      img_masked = img_raw_bytes.set_selected(beam_stop_sel, -2)

      # mask extensive overloads, too
      #top_thresh = data['SATURATED_VALUE']
      #beam_stop_sel = img_raw_bytes >= top_thresh
      #img_masked_2 = img_masked.set_selected(beam_stop_sel, -2)

      data['DATA'] = img_masked
      return data, None

    except Exception as e:
      error = 'IOTA IMPORTER ERROR: Image thresholding failed! {}'.format(e)
      return None, error


  def square_pickle(self, data):
    """ A function to crop the image pickle to a square such that the beam center
        is in the center of image (mostly copied from cxi_image2pickle.py)
        CCTBX.XFEL ONLY """

    error = None

    # only one active area is allowed, and it should be the size of the image.
    try:
      test = flex.int([0, 0, data['SIZE1'], data['SIZE2']]) == data['ACTIVE_AREAS']
      assert test.all_eq(True)
    except AssertionError as e:
      error = 'IOTA IMPORTER ERROR: Image modification failed! {}'.format(e)
      return None, error

    # retrieve parameters from the dictionary
    pixel_size = data['PIXEL_SIZE']
    beam_x = int(round(data['BEAM_CENTER_X'] / pixel_size))
    beam_y = int(round(data['BEAM_CENTER_Y'] / pixel_size))
    width = data['SIZE1']
    height = data['SIZE2']
    pixels = data['DATA']
    left = beam_x
    right = width - left
    top = beam_y
    bottom = height - top
    new_pixels = pixels.deep_copy()
    new_size = None

    # the new image will be twice the size as the smallest distance between the
    # beam center and one of the edges of the images
    if self.params.image_import.square_mode == 'crop':
      try:
        new_half_size = min([right, left, top, bottom])
        new_size = new_half_size * 2
        min_x = int(beam_x - new_half_size)
        min_y = int(beam_y - new_half_size)
        max_x = int(beam_x + new_half_size)
        max_y = int(beam_y + new_half_size)

        if self.params.cctbx_xfel.flip_beamXY:
          new_beam_x = data['BEAM_CENTER_X'] - min_y * pixel_size
          new_beam_y = data['BEAM_CENTER_Y'] - min_x * pixel_size
          new_pixels = pixels[min_x:max_x, min_y:max_y]
        else:
          new_beam_x = data['BEAM_CENTER_X'] - min_x * pixel_size
          new_beam_y = data['BEAM_CENTER_Y'] - min_y * pixel_size
          new_pixels = pixels[min_y:max_y, min_x:max_x]

        assert new_pixels.focus()[0] == new_pixels.focus()[1]
      except (AssertionError, Exception) as e:
        error = 'IOTA IMPORTER ERROR: Image modification failed! {}'.format(e)
        return None, error

    # the new image will be twice the size as the LARGEST distance between the
    # beam center and one of the edges of the images
    elif self.params.image_import.square_mode == 'pad':
      new_half_size = max([right, left, top, bottom])
      new_size = new_half_size * 2
      delta_x = new_half_size - beam_x
      delta_y = new_half_size - beam_y
      new_pixels.resize(flex.grid(new_size, new_size))
      new_pixels.fill(-2)

      if self.params.image_import.flip_beamXY:
        new_beam_x = data['BEAM_CENTER_X'] + delta_y * pixel_size
        new_beam_y = data['BEAM_CENTER_Y'] + delta_x * pixel_size
        new_pixels.matrix_paste_block_in_place(pixels, delta_x, delta_y)
      else:
        new_beam_x = data['BEAM_CENTER_X'] + delta_x * pixel_size
        new_beam_y = data['BEAM_CENTER_Y'] + delta_y * pixel_size
        new_pixels.matrix_paste_block_in_place(pixels, delta_y, delta_x)

    else:
      new_beam_x = data['BEAM_CENTER_X']
      new_beam_y = data['BEAM_CENTER_Y']

    # save the results
    if new_size is not None:
      data['DATA'] = new_pixels
      data['SIZE1'] = new_size
      data['SIZE2'] = new_size
      data['ACTIVE_AREAS'] = flex.int([0, 0, new_size, new_size])
      data['BEAM_CENTER_X'] = new_beam_x
      data['BEAM_CENTER_Y'] = new_beam_y

    return data, error

  def modify_image(self, data=None):
    if not data:
      error = 'IOTA IMPORTER ERROR: Modification failed, no data!'
      return None, error

    # Deactivate image squaring if beam center is in the center of the image
    # CCTBX.XFEL ONLY (Need a better displacement cutoff)
    if (abs(data['BEAM_CENTER_X'] - data['BEAM_CENTER_Y']) < 0.1) :
      self.params.image_import.square_mode = 'no_modification'


    # Check if conversion/modification is required and carry them out
    if (
            self.params.image_import.square_mode != "no_modification" or
            self.params.image_import.beam_center.x != 0 or
            self.params.image_import.beam_center.y != 0 or
            self.params.image_import.beamstop != 0 or
            self.params.image_import.distance != 0
    ):

      # Convert raw image to image pickle
      beamstop = self.params.image_import.beamstop
      distance = self.params.image_import.distance
      beam_center = [self.params.image_import.beam_center.x,
                     self.params.image_import.beam_center.y]
      square = self.params.image_import.square_mode
      mask_file = self.params.image_import.mask

      # Apply mask
      if mask_file is not None:
        data, error = self.apply_mask_from_file(data, mask_file)
        if not data:
          return None, error

      # Override beam center
      if beam_center != [0,0]:
        pixel_size = data['PIXEL_SIZE']
        data['BEAM_CENTER_X'] = int(round(beam_center[0] * pixel_size))
        data['BEAM_CENTER_Y'] = int(round(beam_center[1] * pixel_size))
      if distance != 0:
        data['DISTANCE'] = distance
      if square in ('crop', 'pad'):
        data, error = self.square_pickle(data)
        if not data:
          return None, error
      if beamstop != 0:
        data, error = self.apply_threshold(data)
        if not data:
          return None, error

    # Generate new name for a converted image pickle
    new_path = self.rename_converted_image(self.img_object.img_path)
    self.img_object.raw_img = self.img_object.img_path
    self.img_object.conv_img = new_path
    self.img_object.img_path = new_path

    # Output converted image to file
    try:
      img_dir = os.path.dirname(new_path)
      if not os.path.isdir(img_dir):
        os.makedirs(img_dir)
    except OSError as e:
      error = 'IOTA IMPORTER ERROR: Cannot create folder! {}'.format(e)
      return error, None
    else:
      ep.dump(new_path, data)

    return data, None

  def make_image_object(self, input_entry):
    self.img_object, error = self.import_image(input_entry)

    # Triage image and generate grid search parameters
    if self.img_object.experiments:
      if self.params.image_import.image_triage:
        hmed, amed = self.image_triage()
      else:
        hmed = self.params.cctbx_xfel.grid_search.height_median
        amed = self.params.cctbx_xfel.grid_search.area_median
      hr = self.params.cctbx_xfel.grid_search.height_range
      ar = self.params.cctbx_xfel.grid_search.area_range
      sig_search = self.params.cctbx_xfel.grid_search.sig_height_search
      gs_type = str(self.params.cctbx_xfel.grid_search.type).lower()

      # Generate grid in image object (put that in there so that the HA14
      # integrator can work)
      if not self.img_object.fail:
        self.img_object.generate_grid(gs_type=gs_type, sig_search=sig_search,
                                      hr=hr, ar=ar, hmed=hmed, amed=amed)
        self.img_object.status = 'imported'
      else:
        self.img_object.status = 'final'

    return self.img_object
