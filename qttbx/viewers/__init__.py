"""
Define basic interfaces for communicating with various viewers
"""

import os
import shutil
import sys

from libtbx.utils import Sorry

# =============================================================================
# Phenix atom selections for viewers that select by residue
def _resseq(text):
  from iotbx.pdb.hybrid_36 import hy36decode
  text = text.strip()
  try:
    return int(text)
  except ValueError:
    pass
  try:
    return hy36decode(4, text)
  except (ValueError, RuntimeError):
    raise NotImplementedError()

def _resid(value):
  # a number (decimal or hybrid-36) and any insertion code: 42, 42A, -12A,
  # A000, A000B
  for number, icode in ((value, ''), (value[:-1], value[-1:])):
    if icode.isalpha() or icode == '':
      try:
        return _resseq(number), icode
      except NotImplementedError:
        pass
  raise NotImplementedError()

def _residue_number(keyword, value):
  value = value.strip()
  if ':' in value:
    start, end = value.split(':', 1)
    term = {'start': _resseq(start), 'end': _resseq(end)}
    if keyword != 'resseq':
      # resid 10:20 stops at 20, without 20A; resseq 10:20 takes 20A
      term['icode'] = ''
    return term
  if keyword == 'resseq':
    # any insertion code, as in a Phenix selection
    number = _resseq(value)
    return {'start': number, 'end': number}
  number, icode = _resid(value)
  return {'start': number, 'end': number, 'icode': icode}

def _merge(a, b):
  merged = dict(a)
  for key, value in b.items():
    if merged.get(key, value) != value:
      raise NotImplementedError()
    merged[key] = value
  return merged

def _parse_residues(selection_string):
  '''
  Reads the subset of the Phenix selection syntax that names residues:
  chain, resseq, resid, icode, with name, altloc and resname accepted but not
  needed, joined by and/or with parentheses. Returns one dict per residue or
  residue range.
  '''
  from iotbx.pdb.atom_selection import selection_tokenizer
  words = selection_tokenizer(selection_string)

  def pop():
    word = words.try_pop()
    if word is None:
      raise NotImplementedError()
    return word

  def is_operator(word, operator):
    return (word is not None and word.quote_token is None
            and word.value.lower() == operator)

  def next_is(operator):
    word = words.try_pop()
    if word is not None:
      words.backup()
    return is_operator(word, operator)

  def expression():
    terms = conjunction()
    while next_is('or'):
      pop()
      terms = terms + conjunction()
    return terms

  def conjunction():
    terms = primary()
    while next_is('and'):
      pop()
      rhs = primary()
      terms = [_merge(a, b) for a in terms for b in rhs]
    return terms

  def primary():
    word = pop()
    if is_operator(word, '('):
      terms = expression()
      if not is_operator(pop(), ')'):
        raise NotImplementedError()
      return terms
    keyword = word.value.lower() if word.quote_token is None else None
    if keyword == 'chain':
      return [{'chain': pop().value.strip()}]
    if keyword in ('resseq', 'resid', 'resi'):
      return [_residue_number(keyword, pop().value)]
    if keyword == 'icode':
      return [{'icode': pop().value.strip()}]
    if keyword in ('name', 'altloc', 'altid', 'resname'):
      pop()
      return [{}]
    raise NotImplementedError()

  terms = expression()
  if words.try_pop() is not None:
    raise NotImplementedError()
  return terms

def residue_selection_terms(selection_string):
  '''
  Function for reading the residues named by a Phenix atom selection, for
  viewers that select by residue (ChimeraX, PyMOL)

  Reads the subset of the selection syntax that validation results and
  the Phenix GUI produce: chain, resseq, resid and icode, with name, altloc
  and resname accepted but not needed, joined by and/or with parentheses.

  Parameters
  ----------
    selection_string: str
      e.g. "(chain A and resid 42A and name CB) or (chain ' ' and resid 5)"

  Returns
  -------
    A list of dicts, one per residue or residue range, with start and end
    (residue numbers, hybrid-36 decoded), chain (absent: any chain, '': a
    blank chain ID) and icode (the insertion code of the residue, or of the
    end of a range; absent: any); or None if the selection names anything
    else
  '''
  if not selection_string:
    return None
  try:
    terms = _parse_residues(selection_string)
  except (NotImplementedError, RuntimeError):
    return None
  for term in terms:
    if 'start' not in term:
      return None
    if term.get('icode') and term['start'] != term['end']:
      return None
  return terms

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
