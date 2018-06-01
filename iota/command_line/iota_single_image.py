from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME iota.single_image
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

'''
Author      : Lyubimov, A.Y.
Created     : 05/31/2018
Last Changed: 05/31/2018
Description : IOTA Single Image: can process single image using DIALS, 
with an array of options (i.e. anything from only spotfinding, to indexing, 
space group determination, refinement, integration)
'''

import os
import time
import argparse

from iotbx import phil as ip
from libtbx.easy_mp import parallel_map
from dxtbx.datablock import DataBlockFactory

from threading import Thread

from iota.components.iota_dials import IOTADialsProcessor
from iota.components.iota_dials import phil_scope
from iota.components.iota_threads import IOTATermination
from iota.components.iota_misc import Capturing


def parse_command_args():
  """ Parses command line arguments (only options for now) """
  parser = argparse.ArgumentParser(prog='iota.intercept')
  parser.add_argument('path', type=str, nargs = '?', default = None,
                      help = 'Path to data file')
  parser.add_argument('--backend', type=str, default='dials',
                      help='Backend for processing')
  parser.add_argument('--paramfile', type=str, default=None,
                      help='Parameter file for processing')
  parser.add_argument('--output', type=str, default=None,
                      help='Output filename')
  parser.add_argument('--termfile', type=str, default='.stop',
                      help='Termination signal filename')
  parser.add_argument('--interval', type=int, default=1,
                      help='File-check interval')
  parser.add_argument('--nproc', type=int, default=None,
                      help='Number of processors')
  parser.add_argument('--action', type=str, default='spotfind',
                      help='Code for how far to go; available codes: '
                           'spotfind, index, integrate')

  return parser

class DIALSSpfIdx(Thread):
  def __init__(self,
               img,
               termfile=None,
               paramfile=None,
               output=None,
               backend='dials',
               action_code='spotfind',
               n_processors=1
               ):

    self.img = img
    self.backend = backend
    self.paramfile = paramfile
    self.termfile = termfile
    self.n_processors = n_processors
    self.output = output

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
        self.dials_phil = phil_scope
      self.params = self.dials_phil.extract()

    if self.backend == 'dials':
      self.processor = IOTADialsProcessor(params=self.params)


  def process_image(self):
    if os.path.isfile(self.termfile):
      raise IOTATermination('IOTA_TRACKER: Termination signal received!')
    else:
      with Capturing() as junk_output:
        start = time.time()
        fail = False
        sg = None
        uc = None
        try:
          datablock = DataBlockFactory.from_filenames([self.img])[0]
          observed = self.processor.find_spots(datablock=datablock)
        except Exception:
          fail = True
          pass

        # TODO: Indexing / lattice determination very slow (how to speed up?)
        if self.run_indexing:
          if not fail:
            try:
              experiments, indexed = self.processor.index(
                datablock=datablock, reflections=observed)
            except Exception:
              fail = True
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
              lat = experiments[0].crystal.get_space_group().info()
              sg = str(lat).replace(' ', '')
            except Exception:
              fail = True
              pass

        if self.run_integration:
          if not fail:
            try:
              # Run refinement
              experiments, indexed = self.processor.refine(
                experiments=experiments,
                centroids=indexed)

              integrated = self.processor.integrate(experiments=experiments,
                                                         indexed=indexed)
              frame = self.processor.frame
              unit_cell = frame['observations'][0].unit_cell().parameters()
              uc = ' '.join(['{:.1f}'.format(i) for i in unit_cell])

            except Exception:
              fail = True
              pass

        elapsed = time.time() - start
        return [self.img, len(observed), sg, uc, elapsed]


  def run(self):
    info = self.process_image()
    if info is not None:
      img_path = info[0]
      no_spots = info[1]
      sg_info = info[2]
      uc_info = info[3]
      elapsed = info[4]
      print 'RESULT: ', img_path, no_spots, sg_info, uc_info, '---> ', elapsed
    else:
      print 'RESULT: NONE'




# ============================================================================ #
if __name__ == "__main__":
  import argparse
  args, unk_args = parse_command_args().parse_known_args()

  interceptor = DIALSSpfIdx(img=os.path.abspath(args.path),
                            backend=args.backend,
                            paramfile=args.paramfile,
                            output=args.output,
                            termfile=args.termfile,
                            action_code=args.action)
  interceptor.run()