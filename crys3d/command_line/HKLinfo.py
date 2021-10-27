# LIBTBX_SET_DISPATCHER_NAME phenix.HKLinfo
# LIBTBX_SET_DISPATCHER_NAME cctbx.HKLinfo
from __future__ import absolute_import, division, print_function

import sys
from iotbx.cli_parser import CCTBXParser
from libtbx.utils import multi_out, show_total_time

try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  from libtbx.program_template import ProgramTemplate
import os

from iotbx.gui_tools.reflections import ArrayInfo
from crys3d.hklviewer import hklinfo


class Program(ProgramTemplate):
  prog = os.getenv('LIBTBX_DISPATCHER_NAME')
  description="""
  %(prog)s reflectionfile
or
  %(prog)s reflectionfile philinput.txt

where reflectionfile can be any of the conventional reflection file formats
such as .mtz, .cif, .sca or .hkl file. This will print a table to the screen
listing properties of the reflection data arrays present in the file. The
name of the properties to be listed can be shown by typing:
%(prog)s --show-defaults
and noting which PHIL parameters in the scope \"selected_info\" are set to True.
These can be changed either by specifying these on the command line or by
entering the assignments into a text file, say \"philinput.txt\" that is
submitted on the command line together with the name of the reflection file.
""" % locals()

  datatypes = ['miller_array', 'phil' ]
  master_phil_str ="""
merge_equivalents = False
  .type = bool
  .caption = "merging symmetry equivalent reflections into unique wedge in reciprocal space"
  .short_caption = "merging reflections into a symmetry unique wedge"
""" + ArrayInfo.arrayinfo_phil_str

  def validate(self):
    pass

  def run(self):
    data_file = self.data_manager.get_miller_array_names()[0]
    hklinfo.run(data_file, self.params)


if __name__ == '__main__':
  #import time
  #time.sleep(10) # enough time to attach debugger
  # create parser
  logger = multi_out()
  logger.register('console_output', sys.stdout)

  parser = CCTBXParser(program_class = Program, logger=logger)
  namespace = parser.parse_args(sys.argv[1:])

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
