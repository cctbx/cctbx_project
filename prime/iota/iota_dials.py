from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/10/2014
Last Changed: 02/05/2015
Description : Runs DIALS spotfinding, indexing, refinement and integration
              modules. So far, only spotfinding works. This is very much a work
              in progress
'''

import os
import sys
import prime.iota.iota_misc as misc
from iotbx.phil import parse

# TEMPORARY: The settings in the phil_scope will eventually be constructed from
# settings in IOTA, reasonable defaults and an optional *.phil file
phil_scope = parse('''
  verbosity = 1
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

  include scope dials.algorithms.peak_finding.spotfinder_factory.phil_scope
  include scope dials.algorithms.indexing.indexer.index_only_phil_scope
  include scope dials.algorithms.refinement.refiner.phil_scope
  include scope dials.algorithms.integration.integrator.phil_scope
  include scope dials.algorithms.profile_model.factory.phil_scope
  include scope dials.algorithms.spot_prediction.reflection_predictor.phil_scope

''', process_includes=True)

class Integrator(object):
  """ A class for indexing, integration, etc. using DIALS modules """

  def __init__(self,
               source_image=None,
               object_folder=None,
               gain = 0.32,
               params=None):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    from dxtbx.datablock import DataBlockFactory
    from dials.array_family import flex

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

    # Create the parser
    self.parser = OptionParser(
      phil=current_phil,
      read_datablocks=True,
      read_datablocks_from_images=True)


    self.img = [source_image]
    self.obj_base = object_folder
    self.gain = gain
    with misc.Capturing() as junk_output:
      self.datablock = DataBlockFactory.from_filenames(self.img)[0]
    self.obj_filename = "int_{}".format(os.path.basename(self.img[0]))

  def find_spots(self):
    from dials.array_family import flex
    self.observed = flex.reflection_table.from_observations(self.datablock, self.phil)

  def index(self):

    from dials.algorithms.indexing.indexer import indexer_base

    # Generate imagesets
    imagesets = self.datablock.extract_imagesets()

    # Necessary settings for stills processing
    self.phil.refinement.parameterisation.crystal.scan_varying = False
    self.phil.indexing.scan_range=[]

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
    from dials.array_family import flex

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

  def write_integration_pickles(self, callback=None):

    try:
      picklefilename = self.phil.output.integration_pickle
    except AttributeError:
      return

    if self.phil.output.integration_pickle is not None:

      from libtbx import easy_pickle
      from xfel.command_line.frame_extractor import ConstructFrame
      from dials.array_family import flex

      # Split everything into separate experiments for pickling
      # Is NOT NECESSARY as only a single experiment is always present. Need to
      # work on this later
      for e_number in xrange(len(self.experiments)):
        experiment = self.experiments[e_number]
        e_selection = flex.bool( [r['id']==e_number for r in self.integrated])
        reflections = self.integrated.select(e_selection)

        frame = ConstructFrame(reflections, experiment).make_frame()
        frame["pixel_size"] = experiment.detector[0].get_pixel_size()[0]

        try:
          # if the data was a file on disc, get the path
          event_timestamp = os.path.splitext(self.experiments[0].imageset.paths()[0])[0]
        except NotImplementedError:
          # if the data is in memory only, check if the reader set a timestamp on the format object
          event_timestamp = experiment.imageset.reader().get_format(0).timestamp
        event_timestamp = os.path.basename(event_timestamp)
        if event_timestamp.find("shot-")==0:
           event_timestamp = os.path.splitext(event_timestamp)[0] # micromanage the file name
        if hasattr(self.phil.output, "output_dir"):
          outfile = os.path.join(self.phil.output.output_dir, self.phil.output.integration_pickle%(event_timestamp,e_number))
        else:
          outfile = os.path.join(os.path.dirname(self.phil.output.integration_pickle), self.phil.output.integration_pickle%(event_timestamp,e_number))

        if callback is not None:
          callback(self.phil, outfile, frame)

        easy_pickle.dump(outfile, frame)


  def process_reference(self, reference):
    ''' Load the reference spots. - from xfel_process.py'''
    from dials.array_family import flex
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
    with misc.Capturing() as output:
      self.find_spots()
      self.index()
      self.refine()
      self.integrate()

    with open('test.txt', 'w') as tf:
      for i in output:
        tf.write('\n{}'.format(i))
      tf.write('\n\n{} {} {}'.format(len(self.observed), len(self.indexed), len(self.integrated)))

# ============================================================================ #

if __name__ == "__main__":

  test = Integrator(sys.argv[1])
  test.find_spots()

  print len(test.observed)

  test.index()
  print len(test.indexed)
