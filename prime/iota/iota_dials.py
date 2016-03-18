from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/10/2014
Last Changed: 02/24/2015
Description : Runs DIALS spotfinding, indexing, refinement and integration
              modules. So far, only spotfinding works. This is very much a work
              in progress
'''

import os
import sys
import prime.iota.iota_misc as misc
from iotbx.phil import parse
from dials.array_family import flex

# TEMPORARY: The settings in the phil_scope will eventually be constructed from
# settings in IOTA, reasonable defaults and an optional *.phil file
phil_scope = parse('''
  verbosity = 10
    .type = int(value_min=0)
    .help = "The verbosity level"

  input {
    template = None
      .type = str
      .help = "The image sweep template"
      .multiple = True
   }

  output {
    datablock_filename = datablock.json
      .type = str
      .help = "The filename for output datablock"

    strong_filename = strong.pickle
      .type = str
      .help = "The filename for strong reflections from spot finder output."

    indexed_filename = indexed.pickle
      .type = str
      .help = "The filename for indexed reflections."

    refined_experiments_filename = refined_experiments.json
      .type = str
      .help = "The filename for saving refined experimental models"

    integrated_filename = integrated.pickle
      .type = str
      .help = "The filename for final integrated reflections."

    profile_filename = profile.phil
      .type = str
      .help = "The filename for output reflection profile parameters"

    integration_pickle = int-%s_%d.pickle
      .type = str
      .help = Output integration results for each color data to separate cctbx.xfel-style pickle files
  }

  include scope dials.algorithms.spot_finding.factory.phil_scope
  include scope dials.algorithms.indexing.indexer.index_only_phil_scope
  include scope dials.algorithms.refinement.refiner.phil_scope
  include scope dials.algorithms.integration.integrator.phil_scope
  include scope dials.algorithms.profile_model.factory.phil_scope
  include scope dials.algorithms.spot_prediction.reflection_predictor.phil_scope

''', process_includes=True)


class Triage(object):
  """ Performs quick spotfinding (with mostly defaults) and determines if the number of
      found reflections is above the minimum, and thus if the image should be accepted
      for further processing.
  """

  def __init__(self, img, gain, params):
    """ Initialization and data read-in
    """
    from dxtbx.datablock import DataBlockFactory

    self.gain = gain
    self.params = params

    # Read settings from the DIALS target (.phil) file
    # If none is provided, use default settings (and may God have mercy)
    if self.params.dials.target != None:
      with open(self.params.dials.target, 'r') as settings_file:
        settings_file_contents = settings_file.read()
      settings = parse(settings_file_contents)
      current_phil = phil_scope.fetch(sources=[settings])
      self.phil = current_phil.extract()
    else:
      self.phil = phil_scope.extract()

    # Convert raw image into single-image datablock
    with misc.Capturing() as junk_output:
      self.datablock = DataBlockFactory.from_filenames([img])[0]

  def triage_image(self):
    """ Perform triage by running spotfinding and analyzing results
    """

    # Set spotfinding params
    self.phil.spotfinder.threshold.xds.global_threshold = self.params.dials.global_threshold
    self.phil.spotfinder.threshold.xds.gain = self.gain

    # Perform spotfinding
    observed = flex.reflection_table.from_observations(self.datablock, self.phil)

    # Triage the image
    if len(observed) >= self.params.image_triage.min_Bragg_peaks:
      log_info = 'ACCEPTED! {} observed reflections.'.format(len(observed))
      status = None
    else:
      log_info = 'REJECTED!'
      status = 'failed triage'

    return status, log_info


class Integrator(object):
  """ A class for indexing, integration, etc. using DIALS modules """

  def __init__(self,
               source_image,
               object_folder,
               final_folder,
               final_filename,
               final,
               gain = 0.32,
               params=None):
    '''Initialise the script.'''
    from dxtbx.datablock import DataBlockFactory

    self.params = params

    # Read settings from the DIALS target (.phil) file
    # If none is provided, use default settings (and may God have mercy)
    if self.params.dials.target != None:
      with open(self.params.dials.target, 'r') as settings_file:
        settings_file_contents = settings_file.read()
      settings = parse(settings_file_contents)
      current_phil = phil_scope.fetch(sources=[settings])
      self.phil = current_phil.extract()
    else:
      self.phil = phil_scope.extract()

    # Create the parser - probably don't need!
#     from dials.util.options import OptionParser
#     self.parser = OptionParser(
#       phil=current_phil,
#       read_datablocks=True,
#       read_datablocks_from_images=True)

    # Set general file-handling settings
    file_basename = os.path.basename(source_image).split('.')[0]
    self.phil.output.datablock_filename = "{}/{}.json".format(object_folder, file_basename)
    self.phil.output.indexed_filename = "{}/{}_indexed.pickle".format(object_folder, file_basename)
    self.phil.output.strong_filename = "{}/{}_strong.pickle".format(object_folder, file_basename)
    self.phil.output.refined_experiments_filename = "{}/{}_refined_experiments.json".format(object_folder, file_basename)
    self.phil.output.integrated_filename = "{}/{}_integrated.pickle".format(object_folder, file_basename)
    self.phil.output.profile_filename = "{}/{}_profile.phil".format(object_folder, file_basename)
    self.phil.output.integration_pickle = final_filename
    self.int_log = "{}/int_{}.log".format(final_folder, file_basename)

    self.img = [source_image]
    self.obj_base = object_folder
    self.gain = gain
    self.fail = None
    self.frame = None
    self.final = final
    self.final['final'] = final_filename
    with misc.Capturing() as junk_output:
      self.datablock = DataBlockFactory.from_filenames(self.img)[0]
    self.obj_filename = "int_{}".format(os.path.basename(self.img[0]))

  def find_spots(self):

    # Set spotfinding params (NEED TO SCREEN OR DETERMINE)
    self.phil.spotfinder.threshold.xds.global_threshold = self.params.dials.global_threshold
    self.phil.spotfinder.threshold.xds.gain = self.gain

    # Perform spotfinding
    self.observed = flex.reflection_table.from_observations(self.datablock, self.phil)

  def index(self):

    from dials.algorithms.indexing.indexer import indexer_base

    # Generate imagesets
    imagesets = self.datablock.extract_imagesets()

    # Necessary settings for stills processing
    self.phil.refinement.parameterisation.crystal.scan_varying = False
    self.phil.indexing.scan_range=[]

    # Check if unit cell / space group have been provided
    if self.phil.indexing.known_symmetry.space_group == None or \
       self.phil.indexing.known_symmetry.unit_cell == None:
      self.phil.indexing.method = 'fft1d'
    else:
      self.phil.indexing.method == 'real_space_grid_search'

    # Run indexing
    idxr = indexer_base.from_parameters(self.observed, imagesets, params=self.phil)
    self.indexed = idxr.refined_reflections
    self.experiments = idxr.refined_experiments

  def refine(self):
    # From Aaron Brewster: refinement step skipped as it's done in indexing
    # This writes out experiments to disc
    if self.phil.output.refined_experiments_filename:
      from dxtbx.model.experiment.experiment_list import ExperimentListDumper
      dump = ExperimentListDumper(self.experiments)
      dump.as_json(self.phil.output.refined_experiments_filename)

  def integrate(self):

    # Process reference reflections
    self.indexed,_ = self.process_reference(self.indexed)

    # Get integrator from input params
    from dials.algorithms.profile_model.factory import ProfileModelFactory
    from dials.algorithms.integration.integrator import IntegratorFactory

    # Compute the profile model
    self.experiments = ProfileModelFactory.create(self.phil, self.experiments, self.indexed)

    # Predict the reflections
    predicted = flex.reflection_table.from_predictions_multi(
      self.experiments,
      dmin=self.phil.prediction.d_min,
      dmax=self.phil.prediction.d_max,
      margin=self.phil.prediction.margin,
      force_static=self.phil.prediction.force_static)

    # Match the predictions with the reference
    predicted.match_with_reference(self.indexed)

    # Create the integrator
    integrator = IntegratorFactory.create(self.phil, self.experiments, predicted)

    # Integrate the reflections
    self.integrated = integrator.integrate()

    if self.integrated.has_key('intensity.prf.value'):
      method = 'prf' # integration by profile fitting
    elif self.integrated.has_key('intensity.sum.value'):
      method = 'sum' # integration by simple summation
    self.integrated = self.integrated.select(self.integrated['intensity.' + method + '.variance'] > 0) # keep only spots with sigmas above zero

    # Save the reflections if selected
    if self.phil.output.integrated_filename:
      self.save_reflections(self.integrated, self.phil.output.integrated_filename)

    self.write_integration_pickles()
    from dials.algorithms.indexing.stills_indexer import calc_2D_rmsd_and_displacements
    rmsd_indexed, _ = calc_2D_rmsd_and_displacements(self.indexed)
    rmsd_integrated, _ = calc_2D_rmsd_and_displacements(self.integrated)
    crystal_model = self.experiments.crystals()[0]

  def write_integration_pickles(self):
    ''' This is streamlined vs. the code in stills_indexer, since the filename
        convention is set up upstream.
    '''
    from libtbx import easy_pickle
    from xfel.command_line.frame_extractor import ConstructFrame

    self.frame = ConstructFrame(self.integrated, self.experiments[0]).make_frame()
    self.frame["pixel_size"] = self.experiments[0].detector[0].get_pixel_size()[0]
    easy_pickle.dump(self.phil.output.integration_pickle, self.frame)


  def process_reference(self, reference):
    ''' Load the reference spots. - from xfel_process.py'''
    from libtbx.utils import Sorry
    if reference is None:
      return None, None
    assert("miller_index" in reference)
    assert("id" in reference)
    mask = reference.get_flags(reference.flags.indexed)
    rubbish = reference.select(mask == False)
    if mask.count(False) > 0:
      reference.del_selected(mask == False)
    if len(reference) == 0:
      raise Sorry('''
        Invalid input for reference reflections.
        Expected > %d indexed spots, got %d
      ''' % (0, len(reference)))
    mask = reference['miller_index'] == (0, 0, 0)
    if mask.count(True) > 0:
      rubbish.extend(reference.select(mask))
      reference.del_selected(mask)
    mask = reference['id'] < 0
    if mask.count(True) > 0:
      raise Sorry('''
        Invalid input for reference reflections.
        %d reference spots have an invalid experiment id
      ''' % mask.count(True))
    return reference, rubbish

  def save_reflections(self, reflections, filename):
    ''' Save the reflections to file.  - from xfel_process.py (Aaron Brewster) '''
    reflections.as_pickle(filename)


  def run(self):

    log_entry = ['\n']
    with misc.Capturing() as output:
      e = None
      try:
        print "{:-^100}\n".format(" SPOTFINDING: ")
        self.find_spots()
        print "{:-^100}\n\n".format(" FOUND {} SPOTS: ".format(len(self.observed)))
      except Exception, e:
        if hasattr(e, "classname"):
          print e.classname, "for %s:"%self.img[0],
          error_message = "{}: {}".format(e.classname, e[0].replace('\n',' ')[:50])
        else:
          print "Spotfinding error for %s:"%self.img[0],
          error_message = "{}".format(str(e).replace('\n', ' ')[:50])
        print e
        self.fail = 'failed spotfinding'

      if self.fail == None:
        try:
          print "{:-^100}\n".format(" INDEXING: ")
          self.index()
          print "{:-^100}\n\n".format(" USED {} INDEXED REFLECTIONS: ".format(len(self.indexed)))
        except Exception, e:
          if hasattr(e, "classname"):
            print e.classname, "for %s:"%self.img[0],
            error_message = "{}: {}".format(e.classname, e[0].replace('\n',' ')[:50])
          else:
            print "Indexing error for %s:"%self.img[0],
            error_message = "{}".format(str(e).replace('\n', ' ')[:50])
          print e
          self.fail = 'failed indexing'

      if self.fail == None:
        try:
          self.refine()
          print "{:-^100}\n".format(" INTEGRATING: ")
          self.integrate()
          print "{:-^100}\n\n".format(" FINAL {} INTEGRATED REFLECTIONS: ".format(len(self.integrated)))
        except Exception, e:
          if hasattr(e, "classname"):
            print e.classname, "for %s:"%self.img[0],
            error_message = "{}: {}".format(e.classname, e[0].replace('\n',' ')[:50])
          else:
            print "Integration error for %s:"%self.img[0],
            error_message = "{}".format(str(e).replace('\n', ' ')[:50])
          print e
          self.fail = 'failed integration'

    with open(self.int_log, 'w') as tf:

      for i in output:
        if 'cxi_version' not in i:
          tf.write('\n{}'.format(i))

    if self.fail == None:
      # Collect information
      obs = self.frame['observations'][0]
      Bravais_lattice = self.frame['pointgroup']
      cell = obs.unit_cell().parameters()
      res = obs.d_min()

      # Calculate number of spots w/ high I / sigmaI
      Is = obs.data()
      sigmas = obs.sigmas()
      I_over_sigI = Is / sigmas
      spots = len(Is)
      strong_spots = len([i for i in I_over_sigI if i >= self.params.cctbx.selection.min_sigma])

      # Mosaicity parameters
      mosaicity = round((self.frame.get('ML_half_mosaicity_deg', [0])[0]), 6)
      dom_size = self.frame.get('ML_domain_size_ang', [0])[0]
      ewald_proximal_volume = self.frame.get('ewald_proximal_volume', [0])[0]

      # Assemble output for log file and/or integration result file
      p_cell = "{:>6.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}"\
             "".format(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5])

      int_status = 'RES: {:<4.2f}  NSREF: {:<4}  SG: {:<5}  CELL: {}'\
                   ''.format(res, strong_spots, Bravais_lattice, p_cell)

      int_results = {'sg':Bravais_lattice, 'a':cell[0], 'b':cell[1], 'c':cell[2],
                      'alpha':cell[3], 'beta':cell[4], 'gamma':cell[5],
                      'strong':strong_spots, 'res':res, 'mos':mosaicity,
                      'epv':ewald_proximal_volume, 'info':int_status,
                      'ok':True}

      # Update final entry with integration results
      self.final.update(int_results)

      # Generate log summary of integration results
      img_filename = os.path.basename(self.img[0])
      log_entry.append('DIALS integration:')
      log_entry.append('{:<{width}} --->  {}'.format(img_filename, int_status,
                       width = len(img_filename) + 2))

    else:
      # Generate log summary of integration results
      if 'spotfinding' in self.fail:
        step_id = 'SPOTFINDING'
      elif 'indexing' in self.fail:
        step_id = 'INDEXING'
      elif 'integration' in self.fail:
        step_id = 'INTEGRATION'
      log_entry.append('\n {} FAILED - {}'.format(step_id, e))
      int_status = 'not integrated -- {}'.format(e)
      int_results = {'info':int_status}
      self.final['final'] = None

    log_entry = "\n".join(log_entry)

    return self.fail, self.final, log_entry

# ============================================================================ #

if __name__ == "__main__":

  test = Integrator(sys.argv[1])
  test.find_spots()

  print len(test.observed)

  test.index()
  print len(test.indexed)
