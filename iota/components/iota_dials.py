from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/10/2014
Last Changed: 04/13/2017
Description : Runs DIALS spotfinding, indexing, refinement and integration
              modules. The entire thing works, but no optimization of parameters
              is currently available. This is very much a work in progress
'''

import os
import sys

from iotbx.phil import parse
from dxtbx.datablock import DataBlockFactory
from cctbx import sgtbx

from dials.array_family import flex
from dials.command_line.stills_process import phil_scope, Processor
from dials.command_line.refine_bravais_settings import phil_scope as sg_scope
from dials.command_line.refine_bravais_settings import \
  bravais_lattice_to_space_group_table

import iota.components.iota_misc as misc

class IOTADialsProcessor(Processor):
  ''' Subclassing the Processor module from dials.stills_process to introduce
  streamlined integration pickles output '''

  def __init__(self, params):
    self.phil = params
    Processor.__init__(self, params=params)

  def refine_bravais_settings(self, reflections, experiments):
    sgparams = sg_scope.extract()
    sgparams.refinement.reflections.outlier.algorithm = 'tukey'

    from dials.algorithms.indexing.symmetry \
      import refined_settings_factory_from_refined_triclinic

    # Generate Bravais settings
    Lfat = refined_settings_factory_from_refined_triclinic(
      sgparams, experiments, reflections,
      lepage_max_delta=5, nproc=1, refiner_verbosity=0)
    Lfat.labelit_printout()

    # Filter out not-recommended (i.e. too-high rmsd and too-high max angular
    #  difference) solutions
    Lfat_recommended = [s for s in Lfat if s.recommended]

    # If none are recommended, return P1
    if len(Lfat_recommended) == 0:
      return Lfat[-1]

    # Find the highest symmetry group
    possible_bravais_settings = set(solution['bravais'] for solution in
                                    Lfat_recommended)
    bravais_lattice_to_space_group_table(possible_bravais_settings)
    lattice_to_sg_number = {
      'aP': 1, 'mP': 3, 'mC': 5, 'oP': 16, 'oC': 20, 'oF': 22, 'oI': 23,
      'tP': 75, 'tI': 79, 'hP': 143, 'hR': 146, 'cP': 195, 'cF': 196, 'cI': 197
    }
    filtered_lattices = {}
    for key, value in lattice_to_sg_number.iteritems():
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
    ''' Reindex with newly-determined space group / unit cell '''

    # Update space group / unit cell
    experiment = experiments[0]
    print "Old crystal:"
    print experiment.crystal
    print
    experiment.crystal.update(solution.refined_crystal)
    print "New crystal:"
    print experiment.crystal
    print

    # Change basis
    cb_op = solution['cb_op_inp_best'].as_abc()
    change_of_basis_op = sgtbx.change_of_basis_op(cb_op)
    miller_indices = reflections['miller_index']
    non_integral_indices = change_of_basis_op.apply_results_in_non_integral_indices(miller_indices)
    if non_integral_indices.size() > 0:
      print "Removing {}/{} reflections (change of basis results in non-integral indices)" \
            "".format(non_integral_indices.size(), miller_indices.size())
    sel = flex.bool(miller_indices.size(), True)
    sel.set_selected(non_integral_indices, False)
    miller_indices_reindexed = change_of_basis_op.apply(
      miller_indices.select(sel))
    reflections['miller_index'].set_selected(sel, miller_indices_reindexed)
    reflections['miller_index'].set_selected(~sel, (0, 0, 0))

    return experiments, reflections


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

    # Modify settings
    self.phil.output.strong_filename = None

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
    self.int_log = logfile

    # Set customized parameters
    beamX = self.params.image_conversion.beam_center.x
    beamY = self.params.image_conversion.beam_center.y
    if beamX != 0 or beamY != 0:
      self.phil.geometry.detector.slow_fast_beam_centre = '{} {}'.format(
        beamY, beamX)
    if self.params.image_conversion.distance != 0:
      self.phil.geometry.detector.distance = self.params.image_conversion.distance
    if self.params.advanced.estimate_gain:
      self.phil.spotfinder.threshold.xds.gain = gain

    self.img = [source_image]
    self.obj_base = object_folder
    self.fail = None
    self.frame = None
    self.final = final
    self.final['final'] = final_filename
    self.datablock = DataBlockFactory.from_filenames(self.img)[0]
    self.obj_filename = "int_{}".format(os.path.basename(self.img[0]))

  def find_spots(self):
    # Perform spotfinding
    self.observed = self.processor.find_spots(datablock=self.datablock)

  def index(self):
    # Run indexing
    self.experiments, self.indexed = self.processor.index(
      datablock=self.datablock, reflections=self.observed)

  def refine_bravais_settings_and_reindex(self):
    solution = self.processor.refine_bravais_settings(
      reflections=self.indexed, experiments=self.experiments)
    self.experiments, self.indexed = self.processor.reindex(
      reflections=self.indexed,
      experiments=self.experiments,
      solution=solution)

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

      if (                                        self.fail is None and
              self.phil.indexing.known_symmetry.space_group is None and
                             self.params.dials.determine_sg_and_reindex
          ):
        try:
          print "{:-^100}\n".format(" DETERMINING SPACE GROUP : ")
          self.refine_bravais_settings_and_reindex()
          sg = self.experiments[0].crystal.get_space_group().info()
          print "{:-^100}\n".format(" REINDEXED TO SPACE GROUP {} ".format(sg))
        except Exception, e:
          print e

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
