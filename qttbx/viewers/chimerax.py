"""
Interface for ChimeraX using ISOLDE REST API client

https://www.cgl.ucsf.edu/chimerax/
https://isolde.cimr.cam.ac.uk/what-isolde/

"""

from __future__ import absolute_import, division, print_function

import os
import subprocess
import tempfile

from libtbx.utils import Sorry
from qttbx.viewers import ModelViewer

# =============================================================================
class ChimeraXViewer(ModelViewer):

  def start_viewer(self):
    # basic checks
    if self.command is None:
      raise Sorry('The command for ChimeraX is not set.')
    if not os.path.isfile(self.command):
      raise Sorry('The command for ChimeraX is not available.'
                  'Please check that {command} exists.'
                  .format(command=self.command))
    if self.port is None:
      raise Sorry('The port has not been set.')
    if not isinstance(self.port, int):
      raise Sorry('The port should be an integer.')
    if self.port > 65535:
      raise Sorry('The port number should be 65535 or lower')

    # write script
    self.script = tempfile.NamedTemporaryFile()
    self.script.write('isolde remote rest start port {port}'.format(port=self.port))

    # start viewer
    command = self.command + ' ' + self.script.name
    print(command)
    self.process = subprocess.Popen(args=command)

# =============================================================================
if __name__ == '__main__':
  viewer = ChimeraXViewer()
  viewer.command = '/Users/bkpoon/work/ChimeraX-0.91.app/Contents/bin/ChimeraX'
  viewer.port = 12546

  viewer.start_viewer()
# =============================================================================
