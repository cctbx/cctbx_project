# LIBTBX_SET_DISPATCHER_NAME phenix.HKLinfo
# LIBTBX_SET_DISPATCHER_NAME cctbx.HKLinfo
from __future__ import absolute_import, division, print_function

import sys, time
from iotbx.data_manager import DataManager
from iotbx.gui_tools.reflections import ArrayInfo
import textwrap

from crys3d.hklviewer import hklinfo

try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  from libtbx.program_template import ProgramTemplate
import os

master_phil_str = """
data_file_name = None
  .type = path
  .short_caption = Reflection file
  .multiple = False
  .help = Reflection file name
wrap_labels = 15
  .type = int
  .short_caption = Wrap width for labels
  .help = Number of letters for wrapping long miller array labels. If less than 1 no wrapping is done
delimiter = "|"
  .type = str
  .short_caption = column delimiter when printing table to standard output
  .help = column delimiter
"""


class Program(ProgramTemplate):
  datatypes = ['miller_array', 'phil' ]
  master_phil_str = master_phil_str + ArrayInfo.arrayinfo_phil_str

  def validate(self):
    pass

  def run(self):
    data_file = self.data_manager.get_miller_array_names()[0]
    arrays = self.data_manager.get_miller_arrays(filename = data_file)
    hklinfo.run(arrays, self.params)


from iotbx.cli_parser import CCTBXParser
from libtbx.utils import multi_out, show_total_time
from io import StringIO
#from iotbx.cli_parser import run_program

if __name__ == '__main__':
  time.sleep(10) # enough time to attach debugger
  #run_program(program_class=Program)
  #dmp = StringIO()
  dmp = sys.stdout
  # create parser
  #plogger = multi_out()
  #plogger.register('console_output', dmp)
  logger = multi_out()
  logger.register('console_output', sys.stdout)

  parser = CCTBXParser(program_class = Program, logger=logger)
  namespace = parser.parse_args(sys.argv[1:])
  #dmp.close()

  # start program
  print('Starting job', file=logger)
  print('='*79, file=logger)
  task = Program(parser.data_manager, parser.working_phil.extract(), logger=logger)

  # validate inputs
  task.validate()

  # run program
  fnames = task.run()

  # stop timer
  print('', file=logger)
  print('='*79, file=logger)
  print('Job complete', file=logger)
  show_total_time(out=logger)

