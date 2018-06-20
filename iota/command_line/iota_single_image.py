from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME iota.single_image
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

'''
Author      : Lyubimov, A.Y.
Created     : 05/31/2018
Last Changed: 06/20/2018
Description : IOTA Single Image: can process single image using DIALS,
with an array of options (i.e. anything from only spotfinding, to indexing,
space group determination, refinement, integration)
'''

import os
import time

from scitbx.array_family import flex
import numpy as np

from iotbx import phil as ip
from dxtbx.datablock import DataBlockFactory

from threading import Thread

from iota.components.iota_dials import IOTADialsProcessor
# from iota.components.iota_dials import phil_scope
from dials.command_line.stills_process import phil_scope
from iota.components.iota_threads import IOTATermination
from iota.components.iota_misc import Capturing
from iota.components.iota_input import write_defaults


def parse_command_args():
  """ Parses command line arguments (only options for now) """
  parser = argparse.ArgumentParser(prog='iota.intercept')
  parser.add_argument('path', type=str, nargs = '?', default = None,
                      help = 'Path to data file')
  parser.add_argument('--backend', type=str, default='dials',
                      help='Backend for processing')
  parser.add_argument('--paramfile', type=str, default=None,
                      help='Parameter file for processing')
  parser.add_argument('--output_file', type=str, default=None,
                      help='Output filename')
  parser.add_argument('--output_dir', type=str, default=None,
                      help='Output directory (for BluIce)')
  parser.add_argument('--termfile', type=str, default='.stop',
                      help='Termination signal filename')
  parser.add_argument('--index', type=int, default=1,
                      help='Numerical index of the image')
  parser.add_argument('--min_bragg', type=int, default=10,
                      help='Minimum spots for successful spotfinding result')
  parser.add_argument('--nproc', type=int, default=None,
                      help='Number of processors')
  parser.add_argument('--action', type=str, default='spotfind',
                      help='Code for how far to go; available codes: '
                           'spotfind, index, integrate')
  parser.add_argument('--verbose', action = 'store_true',
                      help='Print information to stdout')

  return parser

class DIALSSpfIdx(Thread):
  def __init__(self,
               img,
               index=None,
               termfile=None,
               paramfile=None,
               output_file=None,
               output_dir=None,
               backend='dials',
               action_code='spotfind',
               min_bragg=10,
               n_processors=1,
               verbose=False
               ):

    self.img = img
    self.backend = backend
    self.paramfile = paramfile
    self.termfile = termfile
    self.n_processors = n_processors
    self.index = index
    self.verbose = verbose
    self.min_bragg = min_bragg

    if output_file is not None:
      if output_dir is not None:
        self.output = os.path.join(os.path.abspath(output_dir), output_file)
      else:
        self.output = os.path.abspath(output_file)
    else:
      self.output = None

    Thread.__init__(self)

    # Determine which processes will be included
    if action_code == 'spotfind':
      self.run_indexing = False
      self.run_integration = False
    elif action_code == 'index':
      self.run_indexing = True
      self.run_integration = False
    elif action_code == 'integrate':
      self.run_indexing = True
      self.run_integration = True

    # Initialize IOTA DIALS Processor
    if self.backend.lower() == 'dials':
      if self.paramfile is not None:
        with open(self.paramfile, 'r') as phil_file:
          phil_string = phil_file.read()
        user_phil = ip.parse(phil_string)
        self.dials_phil = phil_scope.fetch(source=user_phil)
      else:
        default_params, _ = write_defaults(method='dials',
                                           write_target_file=False,
                                           write_param_file=False)
        default_phil_string = '\n'.join(default_params)
        default_phil = ip.parse(default_phil_string)
        self.dials_phil = phil_scope.fetch(source=default_phil)

      self.params = self.dials_phil.extract()

    # Modify default DIALS parameters
    # These parameters will be set no matter what
    self.params.output.datablock_filename = None
    self.params.output.indexed_filename = None
    self.params.output.strong_filename = None
    self.params.output.refined_experiments_filename = None
    self.params.output.integrated_filename = None
    self.params.output.integrated_experiments_filename = None
    self.params.output.profile_filename = None
    self.params.output.integration_pickle = None

    # These parameters will be set only if there's no script
    if self.paramfile is None:
      self.params.indexing.stills.method_list = ['fft3d']
      self.params.spotfinder.threshold.dispersion.global_threshold = 75

    if self.backend == 'dials':
      self.processor = IOTADialsProcessor(params=self.params,
                                          write_pickle=False)

  def process_image(self):
    if os.path.isfile(self.termfile):
      raise IOTATermination('IOTA_TRACKER: Termination signal received!')
    else:
      with Capturing() as junk_output:
        err = []
        start = time.time()
        fail = False
        sg = None
        uc = None
        status = None
        score = 0
        try:
          datablock = DataBlockFactory.from_filenames([self.img])[0]
          observed = self.processor.find_spots(datablock=datablock)
          status = 'spots found'
        except Exception, e:
          fail = True
          observed = []
          err.append('SPOTFINDING ERROR: {}'.format(e))
          pass

        # TODO: Indexing / lattice determination very slow (how to speed up?)
        if self.run_indexing:
          if not fail:
            try:
              experiments, indexed = self.processor.index(
                datablock=datablock, reflections=observed)
              score = len(indexed)
            except Exception, e:
              fail = True
              err.append('INDEXING ERROR: {}'.format(e))
              pass

          if not fail:
            try:
              solution = self.processor.refine_bravais_settings(
                reflections=indexed, experiments=experiments)

              # Only reindex if higher-symmetry solution found
              if solution is not None:
                experiments, indexed = self.processor.reindex(
                  reflections=indexed,
                  experiments=experiments,
                  solution=solution)
              obs = experiments
              lat = experiments[0].crystal.get_space_group().info()
              sg = str(lat).replace(' ', '')
              status = 'indexed'

            except Exception:
              fail = True
              err.append('LATTICE ERROR: {}'.format(e))
              pass

          if not fail:
            unit_cell = experiments[0].crystal.get_unit_cell().parameters()
            uc = ' '.join(['{:.4f}'.format(i) for i in unit_cell])

          if self.run_integration:
            if not fail:
              try:
                # Run refinement
                experiments, indexed = self.processor.refine(
                  experiments=experiments,
                  centroids=indexed)
              except Exception, e:
                fail = True
                err.append('REFINEMENT ERROR: {}'.format(e))
                pass

            if not fail:
              try:
                print experiments
                print indexed
                integrated = self.processor.integrate(experiments=experiments,
                                                      indexed=indexed)
                frame = self.processor.frame
                status = 'integrated'
              except Exception, e:
                err.append('INTEGRATION ERROR: {}'.format(e))
                pass

      if status == 'integrated':
        res = frame['observations'][0].d_max_min()
      else:
        detector = datablock.unique_detectors()[0]
        beam = datablock.unique_beams()[0]

        s1 = flex.vec3_double()
        for i in xrange(len(observed)):
          s1.append(detector[observed['panel'][i]].get_pixel_lab_coord(
            observed['xyzobs.px.value'][i][0:2]))
        two_theta = s1.angle(beam.get_s0())
        d = beam.get_wavelength() / (2 * flex.asin(two_theta / 2))
        res = (np.max(d), np.min(d))

      if len(observed) < self.min_bragg:
        res = (99, 99)

      elapsed = time.time() - start
      info = [self.index, len(observed), self.img, sg, uc]
      return status, info, res, score, elapsed, err


  def run(self):
    errors = []
    n_spots = 0
    n_overloads = 0
    res = (99, 99)
    n_rings = 0
    avg_I = 0
    score = 0

    file_wait_start = time.time()
    while True:
      if time.time() - file_wait_start > 30:
        info = None
        elapsed = None
        errors.append('{} does not exist'.format(self.img))
        break
      if os.path.isfile(self.img):
        status, info, res, score, elapsed, err = self.process_image()
        # errors.extend(err)
        break

    if info is not None:
      idx, n_spots, img_path, sg, uc = info
      print 'IMAGE #{}: {}'.format(idx, img_path)
      print 'SPOTS FOUND: {}'.format(n_spots)
      print 'INDEXING: {} INDEXED SPOTS'.format(score)
      if res[0] != 99:
        print 'RESOLUTION: {:.2f} - {:.2f}'.format(res[0], res[1])
      if sg is not None and uc is not None:
        print 'BRAVAIS LATTICE: {}'.format(sg)
        print 'UNIT CELL: {}'.format(uc)
      print 'TOTAL PROCESSING TIME: {:.2f} SEC'.format(elapsed)

      if self.output is not None:
        with open(self.output, 'a') as outf:
          info_line = ' '.join([str(i) for i in info])
          outf.write('{}\n'.format(info_line))

    if self.verbose:
      if errors == []:
        err = ''
        print_errors = False
      else:
        err = errors[0]
        print_errors = True

      print '\n__RESULTS__'
      print '{} {} {} {:.2f} {} {} {} {} {{{}}}' .format(n_spots, n_overloads,
                                      score, res[1], n_rings, 0, avg_I, 0, err)

      if print_errors:
        print "__ERRORS__"
        for e in errors:
          print e




# ============================================================================ #
if __name__ == "__main__":
  import argparse
  args, unk_args = parse_command_args().parse_known_args()

  interceptor = DIALSSpfIdx(img=os.path.abspath(args.path),
                            index=args.index,
                            backend=args.backend,
                            paramfile=args.paramfile,
                            output_file=args.output_file,
                            output_dir=args.output_dir,
                            termfile=args.termfile,
                            action_code=args.action,
                            min_bragg=args.min_bragg,
                            verbose=args.verbose)
  if args.output_dir is not None:
    if not os.path.isdir(args.output_dir):
      os.makedirs(args.output_dir)

  interceptor.start()
