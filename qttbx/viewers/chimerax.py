"""
Interface for ChimeraX using REST API

https://www.cgl.ucsf.edu/chimerax/
https://www.cgl.ucsf.edu/chimerax/docs/user/commands/remotecontrol.html

"""

from __future__ import absolute_import, division, print_function

import glob
import random
import requests
import subprocess
import sys
import time

try:
  from urllib.parse import unquote
except ImportError:
  from urllib import unquote

from libtbx.utils import Sorry
from qttbx.viewers import ModelViewer

# =============================================================================
class ChimeraXViewer(ModelViewer):

  viewer_name = 'ChimeraX'

  # ---------------------------------------------------------------------------
  def start_viewer(self, timeout=60):
    '''
    Function for starting the ChimeraX REST server

    Parameters
    ----------
      timeout: int
        The number of seconds to wait for the REST server to become
        available before raising a Sorry

    Returns
    -------
      Nothing
    '''

    # set some standard search locations for each platform
    if self.command is None:
      search_paths = []
      if sys.platform == 'darwin':
        search_paths = glob.glob('/Applications/ChimeraX*.app/Contents/MacOS')
        search_paths.sort(reverse=True)

      if len(search_paths) > 0:
        self.command = self.find_command(cmd=self.viewer_name, path=search_paths)

    # randomly select port
    if self.port is None:
      self.port = random.randint(49152, 65535)

    # set REST server information
    self.flags = ['--cmd', 'remotecontrol rest start port {}'.format(self.port)]
    self.url = "http://127.0.0.1:{}/".format(self.port)

    self.run_basic_checks()

    # construct ChimeraX command
    # ChimeraX --cmd "remotecontrol rest start port <port>"
    cmd = [self.command] + self.flags

    self.process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # start ChimeraX server and wait until it is ready
    print()
    print('-'*79)
    print('Starting ChimeraX REST server')
    print(self.url)
    counter = 0
    while counter<timeout:
      output = self._check_status()
      if self._connected:
        break
      counter += 1
      time.sleep(1)
    if not self._connected:
      raise Sorry('The ChimeraX REST server is not reachable at {} after '
                  '{} seconds.'.format(self.url, counter))
    print('ChimeraX is ready')
    print('-'*79)
    print()

  # ---------------------------------------------------------------------------
  def _check_status(self):
    '''
    Check if the REST server is available
    '''
    output = None
    try:
      output = requests.get(url=self.url + 'cmdline.html')
      if output.status_code == 200:
        self._connected = True
    except requests.exceptions.ConnectionError:
      self._connected = False
    return output

  # ---------------------------------------------------------------------------
  def _run_command(self, params):
    '''
    Make requests call to REST server
    https://www.cgl.ucsf.edu/chimerax/docs/user/commands/remotecontrol.html
    '''
    output = None
    try:
      r = requests.Request('GET', self.url + 'run', params=params)
      p = r.prepare()
      p.url = unquote(p.url)
      s = requests.Session()
      output = s.send(p)
      # output = requests.get(url=self.url + 'run', params=params)
    except requests.exceptions.ConnectionError:
      pass
    return output

  # ---------------------------------------------------------------------------
  def close_viewer(self):
    print('='*79)
    if self._connected:
      print('Shutting down ChimeraX')
      params = {'command': 'quit'}
      self._run_command(params)
    else:
      print('ChimeraX already shut down')
    rc = self.process.returncode
    stdout, stderr = self.process.communicate()
    # print('-'*79)
    # print(stdout)
    # print('-'*79)
    # print(stderr)
    print('='*79)

  # ---------------------------------------------------------------------------
  def _load_file(self, filename=None):
    params = {'command': 'open+{}'.format(filename) }
    return self._run_command(params)

  # ---------------------------------------------------------------------------
  def load_model(self, filename=None):
    return self._load_file(filename)

  # ---------------------------------------------------------------------------
  def load_map(self, filename=None):
    return self._load_file(filename)

  # ---------------------------------------------------------------------------
  # functions for wxPython GUI
  def is_alive(self):
    self._check_status()
    return self._connected

  def quit(self):
    return self.close_viewer()

# =============================================================================
