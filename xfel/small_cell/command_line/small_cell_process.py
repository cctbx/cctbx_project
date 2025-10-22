#!/usr/bin/env python
#
# LIBTBX_SET_DISPATCHER_NAME cctbx.small_cell_process

from __future__ import absolute_import, division, print_function

import logging

from dials.util import show_mail_on_error

logger = logging.getLogger('cctbx.small_cell_process')

help_message = '''
DIALS script for processing sparse images.
'''

from dials.command_line.stills_process import phil_scope, Processor as BaseProcessor
from xfel.small_cell.command_line.small_cell_index import small_cell_phil_str
from iotbx.phil import parse
phil_scope.adopt_scope(parse(small_cell_phil_str))

# Use the center of mass (com) for the centroid definition for small cell.
program_defaults_phil_str = """
dispatch {
  refine = True
  hit_finder {
    minimum_number_of_reflections = 3
    maximum_number_of_reflections = 120
  }
}
refinement.parameterisation.crystal.fix = cell
refinement.reflections.outlier.algorithm = null
profile {
  gaussian_rs {
    centroid_definition = *com s1
  }
}
"""
phil_scope = phil_scope.fetch(parse(program_defaults_phil_str))

class Processor(BaseProcessor):
  def index(self, experiments, reflections):
    from time import time
    import copy
    from xfel.small_cell.small_cell import small_cell_index_detail

    st = time()

    logger.info('*' * 80)
    logger.info('Indexing Strong Spots')
    logger.info('*' * 80)

    params = copy.deepcopy(self.params)

    max_clique_len, experiments, indexed = small_cell_index_detail(experiments, reflections, params, write_output=False)

    logger.info('')
    logger.info('Time Taken = %f seconds' % (time() - st))
    return experiments, indexed

if __name__ == '__main__':
  from dials.command_line import stills_process
  stills_process.Processor = Processor
  stills_process.phil_scope = phil_scope

  with show_mail_on_error():
    script = stills_process.Script()
    script.run()
