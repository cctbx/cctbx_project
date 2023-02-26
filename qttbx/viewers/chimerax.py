"""
Interface for ChimeraX using REST API

https://www.cgl.ucsf.edu/chimerax/
https://www.cgl.ucsf.edu/chimerax/docs/user/commands/remotecontrol.html

"""

from __future__ import absolute_import, division, print_function

import glob
import os
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
  def is_available(self):
    '''
    Function for determining if ChimeraX is available

    Parameters
    ----------
      Nothing

    Returns
    -------
      True if ChimeraX is available
    '''
    self.find_viewer()
    if self.command:
      return True
    else:
      return False

  def find_viewer(self):
    '''
    Function for finding ChimeraX

    Parameters
    ----------
      Nothing

    Returns
    -------
      Command for running ChimeraX
    '''

    if self.command is not None:
      return self.command  # no work to do

    search_paths = os.getenv('PATH')
    search_command = self.viewer_name
    if search_paths is not None:
      if sys.platform == 'win32':
        search_paths = search_paths.split(';')
      else:
        search_paths = search_paths.split(':')

    # /Applications/ChimeraX<version>.app
    if sys.platform == 'darwin':
      known_paths = glob.glob('/Applications/ChimeraX*.app/Contents/MacOS')
      known_paths.sort(reverse=True)
      search_paths += known_paths
    # /usr/bin/chimerax
    elif sys.platform.startswith('linux'):
      search_paths += ['/usr/bin']
      search_command = 'chimerax'
    # C:\Program Files\ChimeraX <version>\bin
    elif sys.platform == 'win32':
      known_paths = glob.glob('C:\\Program Files\\ChimeraX*\\bin')
      known_paths.sort(reverse=True)
      search_paths += known_paths
      search_command = 'ChimeraX.exe'

    if len(search_paths) > 0:
      self.command = self.find_command(cmd=search_command, path=search_paths)

    return self.command

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

    # append some standard search locations for each platform
    # the existing $PATH is searched first
    if self.command is None:
      self.find_viewer()

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

    # clean environment for launching ChimeraX
    env = os.environ.copy()
    for v in ['PYTHONPATH', 'LD_LIBRARY_PATH', 'DYLD_LIBRARY_PATH', 'DYLD_FALLBACK_LIBRARY_PATH']:
      env.pop(v, None)

    # start ChimeraX server and wait until it is ready
    self.process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env)

    print()
    print('-'*79)
    print('Starting ChimeraX REST server')
    print(self.command)
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
  def send_command(self, cmds=None):
    params = {'command': "+".join(cmds)}
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
