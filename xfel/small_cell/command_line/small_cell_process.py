#!/usr/bin/env python
#
# LIBTBX_SET_DISPATCHER_NAME cctbx.small_cell_process

from __future__ import absolute_import, division, print_function

import logging
logger = logging.getLogger('cctbx.small_cell_process')

help_message = '''
DIALS script for processing sparse images.
'''

from dials.command_line.stills_process import phil_scope, Processor as BaseProcessor
from xfel.small_cell.command_line.small_cell_index import small_cell_phil_str
from iotbx.phil import parse
phil_scope.adopt_scope(parse(small_cell_phil_str))

class Processor(BaseProcessor):
  def index(self, datablock, reflections):
    from time import time
    import copy
    from xfel.small_cell.small_cell import small_cell_index_detail

    st = time()

    logger.info('*' * 80)
    logger.info('Indexing Strong Spots')
    logger.info('*' * 80)

    params = copy.deepcopy(self.params)

    max_clique_len, experiments, indexed = small_cell_index_detail(datablock, reflections, params, write_output=False)

    logger.info('')
    logger.info('Time Taken = %f seconds' % (time() - st))
    return experiments, indexed

from dials.command_line import stills_process
stills_process.Processor = Processor

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = stills_process.Script()
    script.run()
  except Exception as e:
    halraiser(e)
