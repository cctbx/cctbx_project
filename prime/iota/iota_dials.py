from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/10/2014
Last Changed: 11/02/2015
Description : Runs DIALS spotfinding, indexing, refinement and integration
              modules. So far, only spotfinding works. This is very much a work
              in progress
'''

import os
import prime.iota.iota_misc as misc

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

    from iotbx.phil import parse
    from xfel.command_line.xfel_process import phil_scope

    phil_scope = parse('''
      include scope xfel.command_line.xtc_process.phil_scope
    ''', process_includes=True)

    sub_phil_scope = parse('''
      output {
        cxi_merge_picklefile = None
          .type = str
          .help = Output integration results for each color data to separate cctbx.xfel-style pickle files
      }
      indexing {
        stills {
          ewald_proximity_resolution_cutoff = 2.0
            .type = float
            .help = For calculating the area under the green curve, or the acceptable
            .help = volume of reciprocal space for spot prediction, use this high-resolution cutoff
        }
      }
      cxi_merge {
        include scope xfel.command_line.cxi_merge.master_phil
      }
    ''', process_includes=True)

    phil_scope.adopt_scope(sub_phil_scope)

    # Create the parser
    self.parser = OptionParser(
      phil=phil_scope,
      read_datablocks=True,
      read_datablocks_from_images=True)

    self.params = params
    self.img = [source_image]
    self.obj_base = object_folder
    self.phil = phil_scope.extract()
    with misc.Capturing() as junk_output:
      self.datablock = DataBlockFactory.from_filenames(self.img)[0]

    self.obj_filename = "int_{}".format(os.path.basename(self.img[0]))
    self.phil.output.cxi_merge_picklefile = os.path.join(self.obj_base, self.img[0])

  def find_spots(self):
    from dials.array_family import flex

#     self.phil.spotfinder.filter.min_spot_size = self.params.dials.min_spot_size
#     self.phil.spotfinder.threshold.xds.gain = 0.99
#     self.phil.spotfinder.threshold.xds.sigma_background = 6
    self.phil.spotfinder.threshold.xds.global_threshold=75

    self.observed = flex.reflection_table.from_observations(self.datablock, self.phil)
    print len(self.observed)

  def index(self):

    import copy
    from cxi_xdr_xes.two_color.stills_indexer import stills_indexer

    imagesets = self.datablock.extract_imagesets()
    params = copy.deepcopy(self.phil)
    params.refinement.parameterisation.crystal.scan_varying = False
    params.indexing.scan_range=[]



    idxr = stills_indexer.from_parameters(self.observed, imagesets, params=params)
    self.indexed = idxr.refined_reflections
    self.experiments = idxr.refined_experiments



  def refine(self):
    pass

  def integrate(self):
    pass

# ============================================================================ #

if __name__ == "__main__":

  test = Integrator(test_file)
  test.find_spots()

  print len(test.observed)

  test.index()
  print len(test.indexed)
