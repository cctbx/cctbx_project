"""
Interface for ChimeraX using REST API

https://www.cgl.ucsf.edu/chimerax/
https://www.cgl.ucsf.edu/chimerax/docs/user/commands/remotecontrol.html

"""

import glob
import os
import random
import re
import requests
import subprocess
import sys
import tempfile
import threading
import time

try:
  from urllib.parse import unquote
except ImportError:
  from urllib import unquote

from libtbx import launchd_process
from libtbx.utils import Sorry
from qttbx.viewers import ModelViewer, residue_selection_terms

# =============================================================================
def chimerax_residue_spec(selection_string):
  '''
  Function for converting a Phenix atom selection into a ChimeraX spec
  for the residues it names

  ChimeraX cannot center on a point like Coot, so the Phenix GUI centers it
  on the residues of a validation result instead. Residues rather than
  atoms because Probe clash atoms are often hydrogens that the displayed
  model does not have.

  Parameters
  ----------
    selection_string: str
      A selection such as "chain A and resid 42A and name CB", or the forms
      stored by older versions ("chain 'A' and resseq 42", resseq ranges)

  Returns
  -------
    A spec such as "/A:42A | /B:10-20", or None if the selection does not
    name residues (see qttbx.viewers.residue_selection_terms)
  '''
  terms = residue_selection_terms(selection_string)
  if terms is None:
    return None
  return ' | '.join(spec for spec, single in _residue_specs(terms))

def _residue_specs(terms):
  # (spec, names one residue) for each residue_selection_terms term
  specs = []
  for term in terms:
    start, end = term['start'], term['end']
    icode = term.get('icode', '')
    if start != end and 'icode' in term:
      # resid 10:20 leaves out 20A, as 10-20 does
      number = '%d-%d' % (start, end)
    elif start != end:
      # resseq 10:20 takes 20A, which 10-20 alone leaves out
      number = '%d-%d,%d-%d' % (start, end, end, end)
    elif 'icode' in term:
      number = '%d%s' % (start, icode)
      if icode.islower():
        # read from PDB files as a capital, from mmCIF as it is
        number += ',%d%s' % (start, icode.upper())
    else:
      # resseq ignores the insertion code; so does 42-42, unlike 42
      number = '%d-%d' % (start, start)
    chain = term.get('chain')
    if chain is None:
      spec = ':' + number
    elif chain.isalnum():
      spec = '/%s:%s' % (chain, number)
    else:
      # :5 alone matches residue 5 of every chain, and /A-2:5 matches nothing
      spec = '(:%s & ::chain_id="%s")' % (number, chain or ' ')
    if spec not in [s for s, single in specs]:
      specs.append((spec, start == end))
  return specs

# =============================================================================
class ChimeraXViewer(ModelViewer):

  viewer_name = 'ChimeraX'

  # Angstroms per screen pixel after centering on a selection
  zoom_pixel_size = 0.04
  # a residue with a copy within this many Angstroms of a selection's
  # x,y,z is kept to that copy
  zoom_near = 6
  # model id and name of the temporary marker at that x,y,z
  zoom_marker_id = 9999
  zoom_marker_name = 'Phenix zoom point'

  def __init__(self):
    super(ChimeraXViewer, self).__init__()
    # held while Probe dots load: ChimeraX answers one request at a time
    self._busy = threading.Lock()

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

    # start ChimeraX server and wait until it is ready; on macOS 26 and later
    # as a launchd job, in a coalition of its own (libtbx.launchd_process)
    try:
      self.process = launchd_process.spawn(
        cmd, lambda: subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env),
        env=env, label_prefix='org.cctbx.qttbx.chimerax', warn=print,
        name='ChimeraX')
    except (launchd_process.LaunchError, OSError) as e:
      raise Sorry('Could not start ChimeraX ({}): {}'.format(self.command, e)) from e

    print()
    print('-'*79)
    print('Starting ChimeraX REST server')
    print(self.command)
    print(self.url)
    counter = 0
    while counter<timeout:
      output = self._check_status(timeout=1)
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
  def _check_status(self, timeout=3):
    '''
    Check if the REST server is available. A server that does not answer
    in time is busy (it answers one request at a time), not gone.
    '''
    output = None
    try:
      output = requests.get(url=self.url + 'cmdline.html', timeout=timeout)
      if output.status_code == 200:
        self._connected = True
    except requests.exceptions.Timeout:
      pass
    except requests.exceptions.ConnectionError:
      self._connected = False
    return output

  # ---------------------------------------------------------------------------
  def _run_command(self, params, timeout=30):
    '''
    Make requests call to REST server, waiting at most timeout seconds
    https://www.cgl.ucsf.edu/chimerax/docs/user/commands/remotecontrol.html
    '''
    output = None
    try:
      r = requests.Request('GET', self.url + 'run', params=params)
      p = r.prepare()
      p.url = unquote(p.url)
      s = requests.Session()
      output = s.send(p, timeout=timeout)
      # output = requests.get(url=self.url + 'run', params=params)
    except requests.exceptions.RequestException:
      pass
    return output

  # ---------------------------------------------------------------------------
  def close_viewer(self, timeout=10):
    print('='*79)
    if self._connected:
      print('Shutting down ChimeraX')
      params = {'command': 'quit'}
      self._run_command(params, timeout=timeout)
    else:
      print('ChimeraX already shut down')
    rc = self.process.returncode
    try:
      stdout, stderr = self.process.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
      # busy; it quits when it gets to the request
      print('ChimeraX has not quit yet')
    # print('-'*79)
    # print(stdout)
    # print('-'*79)
    # print(stderr)
    print('='*79)

  # ---------------------------------------------------------------------------
  def _load_file(self, filename=None, timeout=120):
    # quoted for spaces, and encoded so that '+' and '#' survive; ChimeraX
    # reads file names without escapes, so Windows paths need nothing more
    return self._run_encoded('open "{f}"'.format(f=filename), timeout=timeout)

  # ---------------------------------------------------------------------------
  def send_command(self, cmds=None, timeout=30):
    params = {'command': "+".join(cmds)}
    return self._run_command(params, timeout=timeout)

  # ---------------------------------------------------------------------------
  def zoom_to_selection(self, selection_string, xyz=None, timeout=3):
    '''
    Function for centering ChimeraX on the residues of a Phenix atom
    selection (see chimerax_residue_spec). The residues' atoms are shown
    first because view only frames displayed atoms. The ChimeraX selection
    is left alone.

    Given the selection's x,y,z, a residue that has a copy within zoom_near
    of it is kept to that copy: another open model with the same chain and
    residue IDs, or (for segid models) a chain somewhere else that a row
    names, is left out. Residues without a copy there (a model moved in
    ChimeraX, a pick of residues far apart) and ranges are viewed wherever
    they are. The point is a temporary marker.

    view resets or removes the near and far clip planes. The user's planes
    are put back without cutting into the residues: both planes the same
    distance apart (or wider apart if the residues need it), centered on
    the residues; a single plane just clear of them.

    Parameters
    ----------
      selection_string: str
        Phenix atom selection, e.g. from a validation result
      xyz: tuple of float or None
        Coordinates the selection is at, in the model file's frame
      timeout: float
        Seconds to wait for ChimeraX on each request

    Returns
    -------
      The REST response, or None if the selection names no residues,
      Probe dots are loading, or ChimeraX cannot be reached
    '''
    terms = residue_selection_terms(selection_string)
    if (terms is None) or self._busy.locked():
      return None
    specs = _residue_specs(terms)
    marker = self.zoom_marker_id
    state = self._run_encoded('info models #%d; clip list' % marker,
                              timeout=timeout)
    if state is None:
      return None
    state = self._reply_text(state)
    taken = re.search(r'model id #%d (.*)' % marker, state)
    mine = (taken is None) or (self.zoom_marker_name in taken.group(1))
    # a marker left by a click that timed out is closed
    leftover = ['close #%d' % marker] if (taken and mine) else []
    near = [False] * len(specs)
    zone = '#{m} :<{r:g}'.format(m=marker, r=self.zoom_near)
    use_marker = (xyz is not None) and mine and any(
      single for spec, single in specs)
    if use_marker:
      commands = leftover + [
        'marker #{m} position {x:.3f},{y:.3f},{z:.3f} radius 0.1'.format(
          m=marker, x=xyz[0], y=xyz[1], z=xyz[2]),
        'rename #{m} "{n}"'.format(m=marker, n=self.zoom_marker_name)]
      leftover = []
      # the marker's info line starts each residue's answer (echo would
      # write to the Log panel of the ChimeraX window, not to the reply)
      probed = [i for i, (spec, single) in enumerate(specs) if single]
      for i in probed:
        commands += ['info models #%d' % marker,
                     'info residues ({s}) & ({z})'.format(s=specs[i][0],
                                                          z=zone)]
      response = self._run_encoded('; '.join(commands), timeout=timeout)
      if response is None:
        return None  # the next click closes the marker
      answers = re.split(r'model id #%d ' % marker,
                         self._reply_text(response))[1:]
      for i, answer in zip(probed, answers):
        near[i] = 'residue id' in answer
    # each part in parentheses: a zone takes in everything after it
    target = ' | '.join(
      '(({s}) & ({z}))'.format(s=spec, z=zone) if near[i] else
      '({s})'.format(s=spec) for i, (spec, single) in enumerate(specs))
    planes = self._clip_planes(state)
    response = self._run_encoded(
      '; '.join(leftover + self._view_commands(target, planes)),
      timeout=timeout)
    if (response is not None) and isinstance(planes, float):
      # view brackets the residues; the user's slab is kept if thicker
      bracket = self._clip_planes(self._reply_text(response))
      if isinstance(bracket, float) and (planes > bracket):
        self._run_encoded(
          'clip near {n:.2f} far {f:.2f} position {t}'.format(
            n=-planes / 2, f=planes / 2, t=target), timeout=timeout)
    if use_marker and (response is not None):
      self._run_encoded('close #%d' % marker, timeout=timeout)
    return response

  def _view_commands(self, target, planes):
    '''
    Commands that show and view target and put the user's near and far
    planes back (planes from _clip_planes), never through target
    '''
    commands = ['show {t} atoms'.format(t=target)]
    if planes is None:
      commands.append('view {t} clip false'.format(t=target))
    else:
      # view puts both planes just clear of target
      commands.append('view {t}'.format(t=target))
      if planes in ('near', 'far'):
        other = 'far' if planes == 'near' else 'near'
        commands.append('clip {o} off'.format(o=other))
    commands.append('zoom pixelSize {p:g}'.format(p=self.zoom_pixel_size))
    if isinstance(planes, float):
      commands.append('clip list')  # how far apart view put them
    return commands

  @staticmethod
  def _clip_planes(state):
    '''
    The near and far planes listed by "clip list": None, 'near' or 'far'
    for one plane, or the distance between them
    '''
    number = r'([-+\d.eE]+)'
    plane = re.compile(r'clip (near|far) 0 axis {n},{n},{n} position\s+'
                       r'{n},{n},{n}'.format(n=number))
    # the axes point into the slab, so the planes' depths along them add
    # up to minus its thickness
    depths = {}
    for match in plane.finditer(state):
      values = [float(v) for v in match.groups()[1:]]
      axis, point = values[:3], values[3:]
      along = sum(a * p for a, p in zip(axis, point))
      depths.setdefault(match.group(1), along)
    if len(depths) < 2:
      return (list(depths) or [None])[0]
    return -depths['far'] - depths['near']

  # ---------------------------------------------------------------------------
  def _run_encoded(self, command, timeout=10):
    '''
    Run a command and return the REST response, or None if ChimeraX cannot
    be reached. requests encodes the command; _run_command's unquote would
    let a '#', '&' or '+' in it break the query.
    '''
    if self.url is None:
      return None
    try:
      return requests.get(self.url + 'run', params={'command': command},
                          timeout=timeout)
    except requests.exceptions.RequestException:
      return None

  @staticmethod
  def _reply_text(response):
    # ChimeraX sends text/plain with no charset, which requests reads as
    # ISO-8859-1; the log is UTF-8
    return response.content.decode('utf-8', 'replace')

  # ---------------------------------------------------------------------------
  def find_structure(self, model_file, timeout=10):
    '''
    Function for finding the ChimeraX model id of an open model file

    Parameters
    ----------
      model_file: str or list of str
        Path of the model file, or the paths of several files any of which
        may be open (e.g. the PDB and mmCIF files of one model)
      timeout: float
        Seconds to wait for ChimeraX

    Returns
    -------
      The id (e.g. "3") of the most recently opened structure read from
      model_file, or None if it is not open
    '''
    response = self._run_encoded('info models ##filename attribute filename',
                                 timeout=timeout)
    if response is None:
      return None
    def normalized(path):
      return os.path.normcase(os.path.realpath(path))
    model_files = [model_file] if isinstance(model_file, str) else model_file
    wanted = set(normalized(f) for f in model_files)
    def same_file(path):
      return normalized(path) in wanted
    found = None
    # lines like: model id #3 type AtomicStructure filename "/a b/model.pdb"
    for line in self._reply_text(response).splitlines():
      match = re.match(
        r'model id #([\d.]+) type AtomicStructure filename (.+)$', line.strip())
      if match is None:
        continue
      model_id, path = match.group(1), match.group(2).strip()
      if len(path) > 1 and path[0] == path[-1] == '"':
        path = path[1:-1]
      model_number = [int(n) for n in model_id.split('.')]
      if same_file(path) and (found is None or model_number > found[0]):
        found = (model_number, model_id)
    return found[1] if found is not None else None

  # ---------------------------------------------------------------------------
  def load_probe_dots(self, probe_file, model_file, timeout=600):
    '''
    Function for showing Probe clash dots on a model open in ChimeraX, with
    the PhenixUI bundle's reader for Probe's unformatted (-u) output. Only
    the bad overlaps (the clashes listed in the validation tables) are
    shown: the reader takes seconds per thousand dots. Dots loaded earlier
    onto the same model are closed first, so this also reloads. This can
    take a while; GUI callers run it in a thread.

    Parameters
    ----------
      probe_file: str
        Probe -u output, as written for Coot (e.g. .comm/probe.txt)
      model_file: str or list of str
        Path of the model file the dots belong to, or of several files of
        the model (see find_structure)
      timeout: float
        Seconds to wait for ChimeraX

    Returns
    -------
      The REST response, or None if the probe file is missing, the model is
      not open in ChimeraX, dots are already loading, or ChimeraX cannot be
      reached

    Raises
    ------
      Sorry if ChimeraX did not show the dots, e.g. without PhenixUI 1.7.0
      or newer
    '''
    if (not model_file) or not os.path.isfile(probe_file):
      return None
    if not self._busy.acquire(False):
      return None  # already loading
    try:
      return self._load_probe_dots(probe_file, model_file, timeout)
    finally:
      self._busy.release()

  def _load_probe_dots(self, probe_file, model_file, timeout):
    model_id = self.find_structure(model_file)
    if model_id is None:
      return None
    with open(probe_file) as f:
      clashes = [line for line in f if line.split(':')[2:3] == ['bo']]
    command = 'close #{id} & ##name="Probe markups"'.format(id=model_id)
    if len(clashes) == 0:
      return self._run_encoded(command, timeout=timeout)
    fd, dots_file = tempfile.mkstemp(suffix='.probe')
    with os.fdopen(fd, 'w') as f:
      f.writelines(clashes)
    command += ('; open "{f}" format probe structureModel #{id}'
                '; info models #{id}').format(f=dots_file, id=model_id)
    try:
      response = requests.get(self.url + 'run', params={'command': command},
                              timeout=timeout)
    except requests.exceptions.Timeout:
      return None  # ChimeraX may still be reading the file; leave it
    except requests.exceptions.RequestException:
      response = None
    os.remove(dots_file)
    if response is None:
      return None
    text = self._reply_text(response)
    dots = r'model id #%s\.\d+ type \S+ name "Probe markups"' % re.escape(
      model_id)
    if re.search(dots, text) is None:
      said = (text.strip().splitlines() or [''])[0]
      if len(said) > 200:  # e.g. the list of every format it can open
        said = said[:200] + '...'
      message = 'ChimeraX did not show the Probe clash dots.'
      if '"format"' in said:  # no reader for the probe format
        message += (' Showing them needs the PhenixUI bundle 1.7.0 or newer '
                    '(Tools > More Tools in ChimeraX).')
      raise Sorry('%s ChimeraX said:\n%s' % (message, said))
    return response

  # ---------------------------------------------------------------------------
  def load_model(self, filename=None):
    return self._load_file(filename)

  # ---------------------------------------------------------------------------
  def load_map(self, filename=None):
    return self._load_file(filename)

  # ---------------------------------------------------------------------------
  # functions for wxPython GUI
  def is_alive(self):
    if self._busy.locked():
      return True  # loading Probe dots; asking would wait for them
    self._check_status()
    return self._connected

  def is_running(self):
    '''
    Function for checking the ChimeraX process without asking its REST
    server, which may be busy for seconds

    Returns
    -------
      True if the ChimeraX started by start_viewer has not exited
    '''
    return self.process is not None and self.process.poll() is None

  def quit(self):
    return self.close_viewer()

# =============================================================================
