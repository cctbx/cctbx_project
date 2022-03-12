"""
Define basic interfaces for communicating with various viewers
"""

from __future__ import absolute_import, division, print_function

import os
import shutil
import sys

from libtbx.utils import Sorry

# =============================================================================
class ModelViewer(object):
  """
  Base class to define the interface for communicating with model
  viewers.

  Supported model viewers (future plans):
    * Coot - https://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/
    * ChimeraX - https://www.cgl.ucsf.edu/chimerax/
    * PyMOL - https://pymol.org/2/
    * MolStar - https://molstar.org
    * NGLViewer - http://nglviewer.org
  """

  viewer_name = None

  def __init__(self):
    self.command = None       # path to viewer executable
    self.search_paths = None  # paths to search for viewer executable
    self.flags = None         # flags to pass to executable
    self.script = None        # startup script for connecting with viewer
    self.port = None          # port for communicating with viewer
    self.process = None       # viewer process
    self.url = None           # URL for viewer

    self._connected = False

  def run_basic_checks(self):
    """
    Check that the minimum information is provided about the viewer.
    """
    if self.command is None:
      raise Sorry('The command for {viewer_name} is not set.'
                  .format(viewer_name=self.viewer_name))
    if not os.path.isfile(self.command):
      raise Sorry('The command for {viewer_name} is not available. '
                  'Please check that {command} exists.'
                  .format(viewer_name=self.viewer_name, command=self.command))
    if self.port is None:
      raise Sorry('The port has not been set.')
    if not isinstance(self.port, int):
      raise Sorry('The port should be an integer.')
    if self.port > 65535:
      raise Sorry('The port number should be 65535 or lower')

  # ---------------------------------------------------------------------------
  # Viewer functions
  def find_command(self, cmd=None, path=None):
    """
    Function for finding the program executable. This function can be
    used to set self.command

    Parameters
    ----------
      cmd: str
        The name of the command. If not set, the viewer_name is used.
      path: str or list
        The path or paths to search for the command. The shutil.which
        function is used to search the existing path if this is not set.

    Returns
    -------
      command: str
        The full path to the command
    """
    if cmd is None:
      cmd = self.viewer_name
    if isinstance(path, list):
      if sys.platform == 'win32':
        path = ';'.join(path)
      else:
        path = ':'.join(path)
    if sys.version_info.major == 3:
      return shutil.which(cmd=cmd, path=path)
    else:
      for p in path.split(':'):
        f = os.path.join(p, cmd)
        if os.path.isfile(f):
          return f
        else:
          return None

  def start_viewer(self):
    raise NotImplementedError

  def close_viewer(self):
    raise NotImplementedError

  # ---------------------------------------------------------------------------
  # Model functions
  def load_model(self, filename=None):
    raise NotImplementedError

  def close_model(self, model_id=None):
    raise NotImplementedError

  # ---------------------------------------------------------------------------
  # Map functions
  def load_map(self, filename=None):
    raise NotImplementedError

  def load_map_coefficients(self, filename=None, label=None):
    raise NotImplementedError

  def close_map(self, map_id=None):
    raise NotImplementedError

# =============================================================================
