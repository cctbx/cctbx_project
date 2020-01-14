"""
Interface for ChimeraX using ISOLDE REST API client

https://www.cgl.ucsf.edu/chimerax/
https://isolde.cimr.cam.ac.uk/what-isolde/

"""

from __future__ import absolute_import, division, print_function

import os
import subprocess
import sys
import tempfile

from libtbx.utils import to_str, Sorry
from qttbx.viewers import ModelViewer

# =============================================================================
class ChimeraXViewer(ModelViewer):

  viewer_name = 'ChimeraX'

  # ---------------------------------------------------------------------------
  def start_viewer(self):
    self.run_basic_checks()

    # write script
    with tempfile.NamedTemporaryFile(mode='w', suffix='.cxc', delete=False) as self.script:
      self.script.write('isolde remote rest start port {port}'.format(port=self.port))

    # start viewer
    command = [self.command, self.script.name]
    self.process = subprocess.Popen(args=command, stdout=subprocess.PIPE,
      stderr=subprocess.PIPE, shell=False)

    # connect to viewer
    #self._connect()

  # ---------------------------------------------------------------------------
  def close_viewer(self):
    if os.path.isfile(self.script.name):
      os.remove(self.script.name)

  # ---------------------------------------------------------------------------
  def load_model(self, filename=None):
    model_id = None
    if os.path.isfile(filename):
      model_id = self._client.load_model(filename)
    else:
      raise Sorry('Model file ({filename}) is not found.'.format(filename=filename))
    return model_id

  # ---------------------------------------------------------------------------
  def _connect(self):
    # find client.py in ISOLDE
    client_path = None
    with tempfile.NamedTemporaryFile(mode='w', suffix='.py', delete=False) as python_script:
      python_script.write("""\
from chimerax import isolde
print(isolde.__file__)
""")
    command = [self.command, '--nogui', '--exit', '--script', python_script.name]
    python_process = subprocess.Popen(args=command, stdout=subprocess.PIPE,
      stderr=subprocess.PIPE, shell=False)
    stdout, stderr = python_process.communicate()
    python_process.wait()
    if python_process.returncode != 0:
      raise Sorry('The ISOLDE installation location could not be found.')
    else:
      if os.path.isfile(python_script.name):
        try:
          os.remove(python_script.name)
        except IOError:
          pass
      stdout = to_str(stdout).split('\n')
      line = None
      for line in stdout:
        if 'isolde' in line:
          break
      client_path = os.path.abspath(
        os.path.join(line, '..', 'remote_control', 'rest_server', 'client.py'))
    if client_path is not None and os.path.isfile(client_path):
      # import from file
      if sys.version_info.major > 2:
        import importlib.util
        spec = importlib.util.spec_from_file_location('client', client_path)
        client = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(client)
      else:
        import imp
        client = imp.load_source('client', client_path)
      self._client = client.IsoldeRESTClient('localhost', self.port)
      self._client.connect()
      self._connected = True

# =============================================================================
