from __future__ import absolute_import, division, print_function

from libtbx.program_template import ProgramTemplate
from qttbx.viewers.chimerax import ChimeraXViewer

class Program(ProgramTemplate):

  master_phil_str = '''
    command = None
      .type = str
      .help = The full path to the executable. By default, the standard \
              name will be searched in $PATH and some standard locations.
    port = None
      .type = int
      .help = The port for the ChimeraX REST server. By default, a random \
              port between 49152 and 65535 is chosen.
  '''

  datatypes = ['phil']

  def validate(self):
    pass

  def run(self):
    viewer = ChimeraXViewer()

    if self.params.command is not None:
      viewer.command = self.params.command
    if self.params.port is not None:
      viewer.port = self.params.port

    self.output = viewer.start_viewer()

  def get_results(self):
    return self.output
