from __future__ import absolute_import, division, print_function
from six.moves import range

'''
Author      : Lyubimov, A.Y.
Created     : 10/10/2014
Last Changed: 09/18/2019
Description : Runs spotfinding, indexing, refinement and integration using
              subclassed DIALS Stills Processor module. Selector class
              applies filters based on unit cell, space group, etc.
'''

import os
import numpy as np

from iotbx.phil import parse
from cctbx import sgtbx, crystal
import copy

from dials.util import log
from dials.array_family import flex
from dials.algorithms.indexing.symmetry import \
  refined_settings_factory_from_refined_triclinic
from dials.command_line.stills_process import phil_scope, Processor
from dials.command_line.refine_bravais_settings import phil_scope as sg_scope
from dials.command_line.refine_bravais_settings import \
  bravais_lattice_to_space_group_table

import iota.components.iota_utils as util

cctbx_str = '''
cctbx_xfel
  .help = Options for diffraction image processing with current cctbx.xfel
  .alias = Processing Options
{
  target = None
    .type = path
    .multiple = False
    .help = Target (.phil) file with integration parameters for DIALS
    .alias = Processing script
    .expert_level = 0
    .style = path:file
  target_space_group = None
    .type = space_group
    .help = Target space (or point) group (if known)
    .alias = Target Space Group
    .expert_level = 0
  target_unit_cell = None
    .type = unit_cell
    .help = Target unit cell parameters (if known)
    .alias = Target Unit Cell
    .expert_level = 0
  use_fft3d = True
    .type = bool
    .help = Set to True to use FFT3D in indexing
    .alias = Use FFT3D in indexing
    .expert_level = 2
  significance_filter
    .help = Set to True and add value to determine resolution based on I/sigI
    .alias = I/sig(I) cutoff
    .expert_level = 1
    .style = grid:auto
  {
    flag_on = True
      .type = bool
      .help = Set to true to activate significance filter
      .alias = Significance filter
      .style = scope_switch
    sigma = 1.0
      .type = float
      .help = Sigma level to determine resolution cutoff
  }
  determine_sg_and_reindex = True
    .type = bool
    .help = Will determine sg and reindex if no target space group supplied
    .alias = Determine space group and reindex
    .expert_level = 0
  auto_threshold = False
    .type = bool
    .help = Set to True to estimate global threshold for each image
    .alias = Estimate threshold for each image
    .expert_level = 2
  filter
      .help = Throw out results that do not fit user-defined parameters
      .alias = Filters
      .expert_level = 0
      .style = grid:auto
    {
      flag_on = False
        .type = bool
        .help = Set to True to activate filters
        .alias = Use filters
        .style = scope_switch
      crystal_system = None Triclinic Monoclinic Orthorhombic Tetragonal Trigonal Hexagonal Cubic
        .type = choice
        .help = Predicted crystal system for processed data
        .alias = Crystal System
      pointgroup = None
        .type = space_group
        .help = Target Bravais lattice, e.g. "P422"
        .alias = Bravais Lattice
        .expert_level = 1
      unit_cell = None
        .type = unit_cell
        .help = In format of "a, b, c, alpha, beta, gamma", e.g. 79.4, 79.4, 38.1, 90.0, 90.0, 90.0
        .alias = Unit cell
      uc_tolerance = None
        .type = float
        .help = Maximum allowed unit cell deviation from target
        .alias = Unit cell tolerance
        .expert_level = 1
      min_reflections = None
        .type = int
        .help = Minimum integrated reflections per image
        .alias = Min. no. reflections
        .expert_level = 1
      min_resolution = None
        .type = float
        .help = Minimum resolution for accepted images
        .alias = Resolution cutoff
    }
}
'''

class IOTAImageProcessor(Processor):
  """ Subclassed from dials.stills_process. Intended to only be used to
      process a single image; image import, pre-processing, triage,
      and process dispatching are handled by separate modules.

      Added features / overrides:
        - Streamlined writing of integration pickles.
        - Estimation of best-fit Bravais lattice and re-indexing """

  def __init__(self, iparams, write_pickle=True, write_logs=True,
               last_stage='integrate'):
    ''' Constructor
    :param iparams: IOTA params
    :param write_pickle: Set to True to write out an integration pickle
    '''

    self.iparams = iparams
    self.write_pickle = write_pickle
    self.write_logs = write_logs
    self.last_stage = last_stage
    self.dlog_bookmark = 0

    # Get Processor PHIL and initialize Processor
    if self.iparams.cctbx_xfel.target:
      with open(self.iparams.cctbx_xfel.target, 'r') as tf:
        tphil_string = tf.read()
      tparams = phil_scope.fetch(source=parse(tphil_string)).extract()
    else:
      tparams = phil_scope.extract()
    Processor.__init__(self, params=tparams)

    # IOTA-specific settings from here
    # Turn off all peripheral output
    self.params.output.experiments_filename = None
    self.params.output.indexed_filename = None
    self.params.output.strong_filename = None
    self.params.output.refined_experiments_filename = None
    self.params.output.integrated_experiments_filename = None
    self.params.output.integrated_filename = None
    self.params.output.profile_filename = None

    # Set customized parameters
    beamX = self.iparams.image_import.beam_center.x
    beamY = self.iparams.image_import.beam_center.y
    if beamX != 0 or beamY != 0:
      self.params.geometry.detector.slow_fast_beam_centre = '{} {}'.format(
        beamY, beamX)
    if self.iparams.image_import.distance != 0:
      self.params.geometry.detector.distance = self.iparams.image_import.distance
    if self.iparams.image_import.mask is not None:
      self.params.spotfinder.lookup.mask = self.iparams.image_import.mask
      self.params.integration.lookup.mask = self.iparams.image_import.mask
    if self.iparams.cctbx_xfel.target_space_group is not None:
      sg = self.iparams.cctbx_xfel.target_space_group
      self.params.indexing.known_symmetry.space_group = sg
    if self.iparams.cctbx_xfel.target_unit_cell is not None:
      uc = self.iparams.cctbx_xfel.target_unit_cell
      self.params.indexing.known_symmetry.unit_cell = uc
    if not self.params.indexing.stills.method_list:
      self.params.indexing.stills.method_list = ['fft1d',
                                                 'real_space_grid_search']
    if self.iparams.cctbx_xfel.use_fft3d:
      self.params.indexing.stills.method_list.insert(2, 'fft3d')
    if self.iparams.cctbx_xfel.significance_filter.flag_on:
      sigma = self.iparams.cctbx_xfel.significance_filter.sigma
      sigma = sigma if sigma else 1
      self.params.significance_filter.enable = True
      self.params.significance_filter.isigi_cutoff = sigma


    # Load reference geometry
    self.reference_detector = None
    if self.iparams.advanced.reference_geometry:

      from dxtbx.model.experiment_list import ExperimentListFactory
      try:
        ref_experiments = ExperimentListFactory.from_json_file(
          str(self.iparams.advanced.reference_geometry), check_format=False
        )
      except Exception as e:
        print ('DEBUG: Could not make exp. list because: ', e)
        try:
          import dxtbx
          img = dxtbx.load(str(self.iparams.advanced.reference_geometry))
        except Exception:
          print(
            "DEBUG: Couldn't load geometry file {}"
            "".format(self.iparams.advanced.reference_geometry)
          )
        else:
          self.reference_detector = img.get_detector()
      else:
        assert len(ref_experiments.detectors()) == 1
        self.reference_detector = ref_experiments.detectors()[0]


  def refine_bravais_settings(self, reflections, experiments):

    # configure DIALS logging
    log.config(verbosity=0, logfile=self.dials_log)

    proc_scope = phil_scope.format(python_object=self.params)
    sgparams = sg_scope.fetch(proc_scope).extract()
    sgparams.refinement.reflections.outlier.algorithm = 'tukey'

    crystal_P1 = copy.deepcopy(experiments[0].crystal)

    # Generate Bravais settings
    try:
      Lfat = refined_settings_factory_from_refined_triclinic(sgparams,
                                                             experiments,
                                                             reflections,
                                                             lepage_max_delta=5)
    except Exception as e:
      # If refinement fails, reset to P1 (experiments remain modified by Lfat
      # if there's a refinement failure, which causes issues down the line)
      for expt in experiments:
        expt.crystal = crystal_P1
      return None

    Lfat.labelit_printout()

    # Filter out not-recommended (i.e. too-high rmsd and too-high max angular
    # difference) solutions
    Lfat_recommended = [s for s in Lfat if s.recommended]

    # If none are recommended, return None (do not reindex)
    if len(Lfat_recommended) == 0:
      return None

    # Find the highest symmetry group
    possible_bravais_settings = set(solution['bravais'] for solution in
                                    Lfat_recommended)
    bravais_lattice_to_space_group_table(possible_bravais_settings)
    lattice_to_sg_number = {
      'aP': 1, 'mP': 3, 'mC': 5, 'oP': 16, 'oC': 20, 'oF': 22, 'oI': 23,
      'tP': 75, 'tI': 79, 'hP': 143, 'hR': 146, 'cP': 195, 'cF': 196, 'cI': 197
    }
    filtered_lattices = {}
    for key, value in lattice_to_sg_number.items():
      if key in possible_bravais_settings:
        filtered_lattices[key] = value

    highest_sym_lattice = max(filtered_lattices, key=filtered_lattices.get)
    highest_sym_solutions = [s for s in Lfat if s['bravais'] == highest_sym_lattice]
    if len(highest_sym_solutions) > 1:
      highest_sym_solution = sorted(highest_sym_solutions,
                                    key=lambda x: x['max_angular_difference'])[0]
    else:
      highest_sym_solution = highest_sym_solutions[0]

    return highest_sym_solution

  def reindex(self, reflections, experiments, solution):
    """ Reindex with newly-determined space group / unit cell """

    # Update space group / unit cell
    experiment = experiments[0]
    print ("Old crystal:")
    print (experiment.crystal, '\n')
    experiment.crystal.update(solution.refined_crystal)
    print ("New crystal:")
    print (experiment.crystal, '\n')

    # Change basis
    cb_op = solution['cb_op_inp_best'].as_abc()
    change_of_basis_op = sgtbx.change_of_basis_op(cb_op)
    miller_indices = reflections['miller_index']
    non_integral_indices = change_of_basis_op.apply_results_in_non_integral_indices(miller_indices)
    if non_integral_indices.size() > 0:
      print ("Removing {}/{} reflections (change of basis results in non-integral indices)" \
            "".format(non_integral_indices.size(), miller_indices.size()))
    sel = flex.bool(miller_indices.size(), True)
    sel.set_selected(non_integral_indices, False)
    miller_indices_reindexed = change_of_basis_op.apply(miller_indices.select(sel))
    reflections['miller_index'].set_selected(sel, miller_indices_reindexed)
    reflections['miller_index'].set_selected(~sel, (0, 0, 0))

    return experiments, reflections

  def write_integration_pickles(self, integrated, experiments, callback=None):
    """ This is streamlined vs. the code in stills_indexer, since the filename
        convention is set up upstream.
    """

    # Construct frame
    from xfel.command_line.frame_extractor import ConstructFrame
    self.frame = ConstructFrame(integrated, experiments[0]).make_frame()
    self.frame["pixel_size"] = experiments[0].detector[0].get_pixel_size()[0]

    if self.write_pickle:
      from libtbx import easy_pickle
      easy_pickle.dump(self.params.output.integration_pickle, self.frame)

  def pg_and_reindex(self, indexed, experiments):
    ''' Find highest-symmetry Bravais lattice '''
    solution = self.refine_bravais_settings(indexed, experiments)
    if solution is not None:
      experiments, indexed = self.reindex(indexed, experiments, solution)
    return experiments, indexed

  def error_handler(self, error, p_name, img_object, output=None):
    if not output:
      output = []

    if hasattr(error, "classname"):
      # print(error.classname, "for {}:".format(img_object.img_path), )
      error_message = "{}: {}".format(error.classname,
                                      error[0].replace('\n', ' ')[:50])
    else:
      p_name = p_name.lower().capitalize()
      # print("{} error for {}:".format(p_name, img_object.img_path), )
      error_message = "{}".format(str(error).replace('\n', ' ')[:50])
    # print(error_message)
    img_object.fail = 'failed {}'.format(p_name.lower())
    img_object.errors.append(error_message)

    # Generate log summary of integration results
    int_status = 'not integrated -- {}'.format(error)
    int_results = {'info': int_status}
    img_object.final['final'] = None
    img_object.final.update(int_results)

    log_entry = '\n{} - {}'.format(img_object.fail.upper(), error)
    img_object.log_info.append(log_entry)

    # Write log entry into log file
    output.append(error_message)
    if self.write_logs:
      self.write_int_log(path=img_object.int_log, output=output,
                         log_entry=log_entry)

    return img_object

  def process(self, img_object):

    # write out DIALS info (tied to self.write_pickle)
    if self.write_pickle:
      pfx = os.path.splitext(img_object.obj_file)[0]
      self.params.output.experiments_filename = pfx + '_imported.expt'
      self.params.output.indexed_filename = pfx + '_indexed.refl'
      self.params.output.strong_filename = pfx + '_strong.refl'
      self.params.output.refined_experiments_filename = pfx + '_refined.expt'
      self.params.output.integrated_experiments_filename = pfx + \
                                                           '_integrated.expt'
      self.params.output.integrated_filename = pfx + '_integrated.refl'

    # Set up integration pickle path and logfile
    self.params.output.integration_pickle = img_object.int_file
    self.int_log = img_object.int_log

    # configure DIALS logging
    dials_logging_dir = os.path.join(img_object.log_path, 'dials_logs')
    self.dials_log = os.path.join(dials_logging_dir,
                                  os.path.basename(self.int_log))
    log.config(verbosity=0, logfile=self.dials_log)

    # Create output folder if one does not exist
    if self.write_pickle:
      if not os.path.isdir(img_object.int_path):
        os.makedirs(img_object.int_path)

    if not img_object.experiments:
      from dxtbx.model.experiment_list import ExperimentListFactory as exp
      img_object.experiments = exp.from_filenames([img_object.img_path])[0]

    # Auto-set threshold and gain (not saved for target.phil)
    if self.iparams.cctbx_xfel.auto_threshold:
      center_int = img_object.center_int if img_object.center_int else 0
      threshold = int(center_int)
      self.params.spotfinder.threshold.dispersion.global_threshold = threshold
    if self.iparams.image_import.estimate_gain:
      self.params.spotfinder.threshold.dispersion.gain = img_object.gain

    # Update geometry if reference geometry was applied
    from dials.command_line.dials_import import ManualGeometryUpdater
    update_geometry = ManualGeometryUpdater(self.params)
    try:
      imagesets = img_object.experiments.imagesets()
      update_geometry(imagesets[0])
      experiment = img_object.experiments[0]
      experiment.beam = imagesets[0].get_beam()
      experiment.detector = imagesets[0].get_detector()
    except RuntimeError as e:
      print("DEBUG: Error updating geometry on {}, {}".format(
        img_object.img_path, e))

    # Set detector if reference geometry was applied
    if self.reference_detector is not None:
      try:
        from dxtbx.model import Detector
        imageset = img_object.experiments[0].imageset
        imageset.set_detector(
          Detector.from_dict(self.reference_detector.to_dict())
        )
        img_object.experiments[0].detector = imageset.get_detector()
      except Exception as e:
        print ('DEBUG: cannot set detector! ', e)

    # Write full params to file (DEBUG)
    param_string = phil_scope.format(python_object=self.params).as_str()
    full_param_dir = os.path.dirname(self.iparams.cctbx_xfel.target)
    full_param_fn = 'full_' + os.path.basename(self.iparams.cctbx_xfel.target)
    full_param_file = os.path.join(full_param_dir, full_param_fn)
    with open(full_param_file, 'w') as ftarg:
      ftarg.write(param_string)

    # **** SPOTFINDING **** #
    with util.Capturing() as output:
      try:
        # debug:
        print ("{:-^100}\n".format(" SPOTFINDING: "))
        print ('<--->')
        observed = self.find_spots(img_object.experiments)
        img_object.final['spots'] = len(observed)
      except Exception as e_spf:
        observed = None
      else:
        if (
                self.iparams.data_selection.image_triage and
                len(observed) >=
                self.iparams.data_selection.image_triage.minimum_Bragg_peaks
        ):
          msg = " FOUND {} SPOTS - IMAGE ACCEPTED!".format(len(observed))
          print("{:-^100}\n\n".format(msg))
        else:
          msg = " FOUND {} SPOTS - IMAGE REJECTED!".format(len(observed))
          print("{:-^100}\n\n".format(msg))
          e = 'Insufficient spots found ({})!'.format(len(observed))
          return self.error_handler(e, 'triage', img_object, output)
    if not observed:
      return self.error_handler(e_spf, 'spotfinding', img_object, output)

    if self.write_logs:
      self.write_int_log(path=img_object.int_log, output=output,
                         dials_log=self.dials_log)

    # Finish if spotfinding is the last processing stage
    if 'spotfind' in self.last_stage:
      try:
        detector = img_object.experiments.unique_detectors()[0]
        beam = img_object.experiments.unique_beams()[0]
      except AttributeError:
        detector = img_object.experiments.imagesets()[0].get_detector()
        beam = img_object.experiments.imagesets()[0].get_beam()

      s1 = flex.vec3_double()
      for i in range(len(observed)):
        s1.append(detector[observed['panel'][i]].get_pixel_lab_coord(
          observed['xyzobs.px.value'][i][0:2]))
      two_theta = s1.angle(beam.get_s0())
      d = beam.get_wavelength() / (2 * flex.asin(two_theta / 2))
      img_object.final['res'] = np.max(d)
      img_object.final['lres'] = np.min(d)
      return img_object

    # **** INDEXING **** #
    with util.Capturing() as output:
      try:
        print ("{:-^100}\n".format(" INDEXING"))
        print ('<--->')
        experiments, indexed = self.index(img_object.experiments, observed)
      except Exception as e_idx:
        indexed = None
      else:
        if indexed:
          img_object.final['indexed'] = len(indexed)
          print ("{:-^100}\n\n".format(" USED {} INDEXED REFLECTIONS "
                                     "".format(len(indexed))))
        else:
          img_object.fail = 'failed indexing'

    if indexed:
      if self.write_logs:
        self.write_int_log(path=img_object.int_log, output=output,
                           dials_log=self.dials_log)
    else:
      return self.error_handler(e_idx, 'indexing', img_object, output)

    with util.Capturing() as output:
      # Bravais lattice and reindex
      if self.iparams.cctbx_xfel.determine_sg_and_reindex:
        try:
          print ("{:-^100}\n".format(" DETERMINING SPACE GROUP"))
          print('<--->')
          experiments, indexed = self.pg_and_reindex(indexed, experiments)
          img_object.final['indexed'] = len(indexed)
          lat = experiments[0].crystal.get_space_group().info()
          sg = str(lat).replace(' ', '')
          if sg != 'P1':
            print ("{:-^100}\n".format(" REINDEXED TO SPACE GROUP {} ".format(sg)))
          else:
            print ("{:-^100}\n".format(" RETAINED TRICLINIC (P1) SYMMETRY "))
          reindex_success = True
        except Exception as e_ridx:
          reindex_success = False

    if reindex_success:
      if self.write_logs:
        self.write_int_log(path=img_object.int_log, output=output,
                           dials_log=self.dials_log)
    else:
      return self.error_handler(e_ridx, 'indexing', img_object, output)

    # **** INTEGRATION **** #
    with util.Capturing() as output:
      try:
        experiments, indexed = self.refine(experiments, indexed)
        print ("{:-^100}\n".format(" INTEGRATING "))
        print ('<--->')
        integrated = self.integrate(experiments, indexed)
      except Exception as e_int:

        integrated = None
      else:
        if integrated:
          img_object.final['integrated'] = len(integrated)
          print ("{:-^100}\n\n".format(" FINAL {} INTEGRATED REFLECTIONS "
                                       "".format(len(integrated))))
    if integrated:
      if self.write_logs:
        self.write_int_log(path=img_object.int_log, output=output,
                           dials_log=self.dials_log)
    else:
      return self.error_handler(e_int, 'integration', img_object, output)

    # Filter
    if self.iparams.cctbx_xfel.filter.flag_on:
      self.selector = Selector(frame=self.frame,
                               uc_tol=self.iparams.cctbx_xfel.filter.uc_tolerance,
                               xsys=self.iparams.cctbx_xfel.filter.crystal_system,
                               pg=self.iparams.cctbx_xfel.filter.pointgroup,
                               uc=self.iparams.cctbx_xfel.filter.unit_cell,
                               min_ref=self.iparams.cctbx_xfel.filter.min_reflections,
                               min_res=self.iparams.cctbx_xfel.filter.min_resolution)
      fail, e = self.selector.result_filter()
      if fail:
        return self.error_handler(e, 'filter', img_object, output)

    int_results, log_entry = self.collect_information(img_object=img_object)

    # Update final entry with integration results
    img_object.final.update(int_results)

    # Update image log
    log_entry = "\n".join(log_entry)
    img_object.log_info.append(log_entry)

    if self.write_logs:
      self.write_int_log(path=img_object.int_log, log_entry=log_entry)
    return img_object

  def write_int_log(self, path, output=None, log_entry=None, dials_log=None):
    if output:
      output = [i for i in output if 'cxi_version' not in i]
      output_string = '\n'.join(output)
    else:
      output_string = ''
    # insert DIALS logging output
    if dials_log is not None:
      with open(dials_log, 'r') as dlog:
        dlog.seek(self.dlog_bookmark)
        dials_log_lines = ['  {}'.format(l) for l in dlog.readlines()]
        dials_log_string = ''.join(dials_log_lines)
        self.dlog_bookmark = dlog.tell()
      output_string = output_string.replace('<--->', dials_log_string)

    # Add IOTA log entry if exists
    if log_entry:
      output_string += log_entry

    # Write log to file
    with open(path, 'a') as tf:
      tf.write(output_string)
      # for i in output:
      #   if 'cxi_version' not in i:
      #     tf.write('\n{}'.format(i))
      # if log_entry:
      #   tf.write('\n{}'.format(log_entry))


  def collect_information(self, img_object):
    # Collect information
    obs = self.frame['observations'][0]
    Bravais_lattice = self.frame['pointgroup']
    cell = obs.unit_cell().parameters()
    lres, hres = obs.d_max_min()

    # Calculate number of spots w/ high I / sigmaI
    Is = obs.data()
    sigmas = obs.sigmas()
    I_over_sigI = Is / sigmas
    strong_spots = len([i for i in I_over_sigI if
                        i >= self.iparams.data_selection.image_triage.strong_sigma])

    # Mosaicity parameters
    mosaicity = round((self.frame.get('ML_half_mosaicity_deg', [0])[0]), 6)
    ewald_proximal_volume = self.frame.get('ewald_proximal_volume', [0])[0]

    # Assemble output for log file and/or integration result file
    p_cell = "{:>6.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}"\
           "".format(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5])

    int_status = 'RES: {:<4.2f}  NSREF: {:<4}  SG: {:<5}  CELL: {}'\
                 ''.format(hres, strong_spots, Bravais_lattice, p_cell)

    int_results = {'sg':Bravais_lattice, 'a':cell[0], 'b':cell[1], 'c':cell[2],
                   'alpha':cell[3], 'beta':cell[4], 'gamma':cell[5],
                   'wavelength':self.frame['wavelength'],
                   'distance':self.frame['distance'],
                   'beamX':self.frame['xbeam'], 'beamY':self.frame['ybeam'],
                   'observations': obs, 'strong':strong_spots,
                   'res':hres, 'lres':lres,
                   'mos':mosaicity, 'epv':ewald_proximal_volume,
                   'info':int_status, 'ok':True}

    # Generate log summary of integration results
    img_filename = os.path.basename(img_object.img_path)
    log_entry = ['\nCCTBX.XFEL integration:']
    log_entry.append('{:<{width}} --->  {}'.format(img_filename, int_status,
                     width = len(img_filename) + 2))

    return int_results, log_entry

  def run(self, img_object):
    return self.process(img_object=img_object)


class Selector(object):
  """ Class for selection of optimal spotfinding parameters from grid search """

  def __init__(self,
               frame,
               uc_tol=0,
               xsys=None,
               pg=None,
               uc=None,
               min_ref=0,
               min_res=None):

    obs = frame['observations'][0]
    self.obs_pg = frame['pointgroup']
    self.obs_uc = [prm for prm in obs.unit_cell().parameters()]
    self.obs_res = obs.d_max_min()[1]
    self.obs_ref = len(obs.data())
    self.uc = uc
    self.uc_tol = uc_tol
    self.xsys = xsys
    self.pg = pg
    self.min_ref = min_ref
    self.min_res = min_res
    self.fail = False

  def result_filter(self):
    """ Unit cell pre-filter. Applies hard space-group constraint and stringent
        unit cell parameter restraints to filter out integration results that
        deviate. Optional step. Unit cell tolerance user-defined. """

    if self.uc is not None:
      user_uc = [prm for prm in self.uc.parameters()]
      delta_a = abs(self.obs_uc[0] - user_uc[0])
      delta_b = abs(self.obs_uc[1] - user_uc[1])
      delta_c = abs(self.obs_uc[2] - user_uc[2])
      delta_alpha = abs(self.obs_uc[3] - user_uc[3])
      delta_beta = abs(self.obs_uc[4] - user_uc[4])
      delta_gamma = abs(self.obs_uc[5] - user_uc[5])
      uc_check = (delta_a <= user_uc[0] * self.uc_tol and
                  delta_b <= user_uc[1] * self.uc_tol and
                  delta_c <= user_uc[2] * self.uc_tol and
                  delta_alpha <= user_uc[3] * self.uc_tol and
                  delta_beta <= user_uc[4] * self.uc_tol and
                  delta_gamma <= user_uc[5] * self.uc_tol)
    else:
      uc_check = True

    e = []
    if not uc_check:
      e.append('UC parameters outside {}% tolerance'.format(self.uc_tol*100))

    if self.obs_ref <= self.min_ref:
      e.append('Fewer than {} observed reflections'.format(self.min_ref))

    if self.min_res and self.obs_res >= self.min_res:
      e.append('Resolution above minimum of {}'.format(self.min_res))

    if self.pg or self.xsys:
      obs_sym = crystal.symmetry(space_group_symbol=self.obs_pg)
      obs_pg = obs_sym.space_group().conventional_centring_type_symbol() + \
               obs_sym.space_group().point_group_type()
      obs_cs = obs_sym.space_group().crystal_system()

      if self.pg:
        fil_sym = crystal.symmetry(space_group_symbol=self.pg)
        fil_pg = fil_sym.space_group().conventional_centring_type_symbol() + \
                 fil_sym.space_group().point_group_type()
        if fil_pg != obs_pg:
          e.append('Point group {} does not match expected ({})'
                   ''.format(self.obs_pg, self.pg))

      if self.xsys:
        if self.xsys != obs_cs:
          e.append('Crystal system {} does not match expected ({})'
                   ''.format(obs_cs, self.xsys))

    if bool(e):
      fail = 'failed filter'
      error = ', '.join(e)
    else:
      fail = None
      error = None

    return fail, error
