#!/usr/bin/env python
#
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.small_cell_process

from __future__ import absolute_import, division, print_function

help_message = '''

Process small cell data with database logging.

'''

from iotbx.phil import parse

from xfel.command_line.xfel_process import control_phil_str, delete_shoeboxes_override_str, radial_average_phil_str
from xfel.small_cell.command_line.small_cell_index import small_cell_phil_str

from xfel.ui.db.frame_logging import DialsProcessorWithLogging
from xfel.small_cell.command_line.small_cell_process import Processor as SmallCellProcessor

class SmallCellProcessorWithLogging(DialsProcessorWithLogging, SmallCellProcessor):
  pass

from dials.util import show_mail_on_error
from dials.command_line.stills_process import dials_phil_str, program_defaults_phil_str, Script as DialsScript, control_phil_str as dials_control_phil_str
from xfel.ui import db_phil_str
from xfel.small_cell.command_line.small_cell_process import program_defaults_phil_str as small_cell_program_defaults_phil_str

phil_scope = parse(dials_control_phil_str + control_phil_str + dials_phil_str + db_phil_str + radial_average_phil_str + small_cell_phil_str,
                   process_includes=True).fetch(parse(program_defaults_phil_str + small_cell_program_defaults_phil_str))
phil_scope = phil_scope.fetch(parse(delete_shoeboxes_override_str))

class Script(DialsScript):
  '''A class for running the script.'''
  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import ArgumentParser
    import libtbx.load_env

    # The script usage
    usage = "usage: %s [options] [param.phil] filenames" % libtbx.env.dispatcher_name

    self.tag = None
    self.reference_detector = None
    self.debug_file_handle = None

    # Create the parser
    self.parser = ArgumentParser(
      usage=usage,
      phil=phil_scope,
      epilog=help_message
      )

if __name__ == '__main__':
  import dials.command_line.stills_process
  dials.command_line.stills_process.Processor = SmallCellProcessorWithLogging

  with show_mail_on_error():
    script = Script()
    script.run()
