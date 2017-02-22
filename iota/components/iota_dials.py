from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/10/2014
Last Changed: 02/21/2017
Description : Runs DIALS spotfinding, indexing, refinement and integration
              modules. The entire thing works, but no optimization of parameters
              is currently available. This is very much a work in progress
'''

import os
import sys
import iota.components.iota_misc as misc
from iotbx.phil import parse
from dxtbx.datablock import DataBlockFactory
from dials.command_line.stills_process import phil_scope, Processor

class IOTADialsProcessor(Processor):
  ''' Subclassing the Processor module from dials.stills_process to introduce
  streamlined integration pickles output '''

  def __init__(self, params):
    self.phil = params
    Processor.__init__(self, params=params)

  def write_integration_pickles(self, integrated, experiments, callback=None):
    ''' This is streamlined vs. the code in stills_indexer, since the filename
        convention is set up upstream.
    '''
    from libtbx import easy_pickle
    from xfel.command_line.frame_extractor import ConstructFrame

    self.frame = ConstructFrame(integrated, experiments[0]).make_frame()
    self.frame["pixel_size"] = experiments[0].detector[0].get_pixel_size()[0]
    easy_pickle.dump(self.phil.output.integration_pickle, self.frame)

class Triage(object):
  """ Performs quick spotfinding (with mostly defaults) and determines if the number of
      found reflections is above the minimum, and thus if the image should be accepted
      for further processing.
  """

  def __init__(self, img, gain, params):
    """ Initialization and data read-in
    """
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

    self.processor = IOTADialsProcessor(params=self.phil)

    # Perform spotfinding
    observed = self.processor.find_spots(datablock=self.datablock)

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
               final_filename,
               final,
               logfile,
               gain = 0.32,
               params=None):
    '''Initialise the script.'''

    self.params = params

    # Read settings from the DIALS target (.phil) file
    # If none is provided, use default settings (and may God have mercy)
    if self.params.dials.target != None:
      with open(self.params.dials.target, 'r') as settings_file:
        settings_file_contents = settings_file.read()
      settings = parse(settings_file_contents)
      current_phil = phil_scope.fetch(source=settings)
      self.phil = current_phil.extract()
    else:
      self.phil = phil_scope.extract()

   # Set general file-handling settings
    file_basename = os.path.basename(source_image).split('.')[0]
    self.phil.output.datablock_filename = "{}/{}.json".format(object_folder, file_basename)
    self.phil.output.indexed_filename = "{}/{}_indexed.pickle".format(object_folder, file_basename)
    self.phil.output.strong_filename = "{}/{}_strong.pickle".format(object_folder, file_basename)
    self.phil.output.refined_experiments_filename = "{}/{}_refined_experiments.json".format(object_folder, file_basename)
    self.phil.output.integrated_filename = "{}/{}_integrated.pickle".format(object_folder, file_basename)
    self.phil.output.profile_filename = "{}/{}_profile.phil".format(object_folder, file_basename)
    self.phil.output.integration_pickle = final_filename
    self.int_log = logfile #"{}/int_{}.log".format(final_folder, file_basename)

    # Set customized parameters
    self.phil.spotfinder.threshold.xds.global_threshold = self.params.dials.global_threshold
    #self.phil.spotfinder.threshold.xds.gain = self.gain
    self.phil.spotfinder.filter.min_spot_size=self.params.dials.min_spot_size

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
    # Perform spotfinding
    self.observed = self.processor.find_spots(datablock=self.datablock)

  def index(self):
    # Run indexing
    self.experiments, self.indexed = self.processor.index(
      datablock=self.datablock, reflections=self.observed)

  def refine(self):
    # Run refinement
    self.experiments, self.indexed = self.processor.refine(
      experiments=self.experiments, centroids=self.indexed)

  def integrate(self):
    # Run integration
    self.integrated = self.processor.integrate(experiments=self.experiments,
                                               indexed=self.indexed)
    self.frame = self.processor.frame

  def run(self):

    self.processor = IOTADialsProcessor(params=self.phil)

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
      lres, hres = obs.d_max_min()

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
                   ''.format(hres, strong_spots, Bravais_lattice, p_cell)

      int_results = {'sg':Bravais_lattice, 'a':cell[0], 'b':cell[1], 'c':cell[2],
                     'alpha':cell[3], 'beta':cell[4], 'gamma':cell[5],
                     'strong':strong_spots, 'res':hres, 'lres':lres,
                     'mos':mosaicity, 'epv':ewald_proximal_volume,
                     'info':int_status, 'ok':True}

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
