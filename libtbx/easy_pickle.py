"""
Utility functions for automatically opening and serializing a python object into
a file using the latest pickling protocol. Can optionally compress the output if
file_name ends with .gz. Also provides functions to read that object back.
"""

from __future__ import absolute_import, division, print_function

import os

from libtbx.str_utils import show_string
from libtbx import smart_open
import six
from six.moves import cPickle as pickle

def _open(file_name, mode):
  """
  Wraps libtbx.smart_open.

  Parameters
  ----------
  file_name : str
  mode : str

  Returns
  -------
  file
  """
  file_name = os.path.expanduser(file_name)
  try: return smart_open.file(file_name=file_name, mode=mode)
  except IOError as e:
    raise IOError("Cannot open pickle file %s (%s)" % (
      show_string(file_name), str(e)))

def dump(file_name, obj):
  """
  Wraps pickle.dump.

  Parameters
  ----------
  file_name : str
  obj : object

  Examples
  --------
  >>> from libtbx.easy_pickle import dump, load
  >>> dump("output.pkl.gz", [1, 2, 3])
  >>> print load("output.pkl.gz")
  [1, 2, 3]
  """
  with _open(file_name, "wb") as f:
    p = pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
  return p

def dumps(obj):
  """
  Wraps pickle.dumps.

  Parameters
  ----------
  obj : object

  Returns
  -------
  str
  """
  return pickle.dumps(obj, pickle.HIGHEST_PROTOCOL)

def load(file_name, faster_but_using_more_memory=True):
  """
  Wraps pickle.load.

  Parameters
  ----------
  file_name : str
  faster_but_using_more_memory : bool, optional
      Optionally read the entirety of a file into memory before converting it
      into a python object.

  Returns
  -------
  object
  """
  if (faster_but_using_more_memory):
    if six.PY2:
      with _open(file_name, "rb") as f:
        s = f.read()
      try:
        return pickle.loads(s)
      except ValueError as e:
        from libtbx.utils import Sorry
        import sys
        raise Sorry(
         "Please check that you are not using an old version "+
          "and reading a file from a new version.  "+
         "(The file %s " %(file_name) + "was saved with a version of Python "+
         "that is not supported in this version (%s): %s)" %(
          sys.version, str(e)))

    with _open(file_name, "rb") as f:
      s = f.read()
    return pickle.loads(s, encoding='bytes')
  with _open(file_name, "rb") as f:
    p = pickle.load(f)
  return p

def loads(string):
  """
  Wraps pickle.loads.

  Parameters
  ----------
  string : str

  Returns
  -------
  object

  Examples
  --------
  >>> from libtbx.easy_pickle import dumps, loads
  >>> print loads(dumps([1, 2, 3])
  [1, 2, 3]
  """
  return pickle.loads(string)

def dump_args(*args, **keyword_args):
  """
  Dumps args and keyword_args to args.pickle.
  """
  dump("args.pickle", (args, keyword_args))

def fix_py2_pickle_orig(p):
  '''
  Fix pickles from Python 2
  Original version

  Parameters
  ----------
  p: pickle

  Returns
  -------
  p: the fixed pickle
  '''
  from collections.abc import Mapping, MutableSequence
  from cctbx.crystal import symmetry
  from cctbx.sgtbx import empty, space_group
  from cctbx.xray.structure import structure
  from iotbx.pdb.hierarchy import root
  skip_types = (empty, root, space_group, structure, symmetry)
  if isinstance(p, Mapping):
    for key in list(p.keys()):
      if isinstance(key, bytes):
        str_key = key.decode('utf8')
        p[str_key] = p[key]
        del p[key]
        key = str_key
      if isinstance(p[key], skip_types) or callable(p[key]):
        pass
      else:
        # print(key, type(p[key]), p[key])
        p[key] = fix_py2_pickle_orig(p[key])
  if isinstance(p, MutableSequence):
    for i in range(len(p)):
      p[i] = fix_py2_pickle_orig(p[i])

  if hasattr(p, '__dict__'):
    p.__dict__ = fix_py2_pickle_orig(p.__dict__)
  # miller array object
  if hasattr(p, '_info') and hasattr(p._info, 'labels'):
    p._info.labels = fix_py2_pickle_orig(p._info.labels)

  if isinstance(p, bytes):
    p = p.decode('utf8')

  return p

def fix_py2_pickle(p):
  '''
  Fix pickles from Python 2
  Version 3

  Parameters
  ----------
  p: pickle

  Returns
  -------
  p: the fixed pickle

  Comments:
  ---------

  '''
  from mmtbx.model.model import get_hierarchy_and_run_hierarchy_method
  from collections.abc import Mapping, MutableSequence
  from libtbx import group_args
  from scitbx_array_family_flex_ext import std_string
  if isinstance(p, get_hierarchy_and_run_hierarchy_method):
    return p
  if isinstance(p, group_args):
    p = p() # now it is a dict
    for key in list(p.keys()):    # fix the key
      if isinstance(key, bytes):
        str_key = key.decode('utf8')
        p[str_key] = p[key]
        del p[key]
        key = str_key
      p[key] = fix_py2_pickle(p[key])
    # convert to dict, fix, convert back to group args
    p = group_args(**p)

  elif isinstance(p, bytes):
    p = p.decode('utf8')

  elif isinstance(p, MutableSequence):
    for i in range(len(p)):
      p[i] = fix_py2_pickle(p[i])

  elif isinstance(p, Mapping):
    for key in list(p.keys()):    # fix the key
      if isinstance(key, bytes):
        str_key = key.decode('utf8')
        p[str_key] = p[key]
        del p[key]
        key = str_key
      p[key] = fix_py2_pickle(p[key])


  elif isinstance(p,tuple):
    p = tuple(fix_py2_pickle(list(p)))

  elif isinstance(p, std_string):
    new_p = std_string()
    for x in p:
      new_p.append(fix_py2_pickle(x))
    p = new_p

  # Classes like mmtbx.monomer_library.cif_types.chem_mod_angle remain here
  elif hasattr(p, '__dict__'):
    for key in list(p.__dict__.keys()):    # fix the key
      if isinstance(key, bytes):
        str_key = key.decode('utf8')
        p.__dict__[str_key] = p.__dict__[key]
        del p.__dict__[key]
        key = str_key
      if not key.startswith("__"):
        p.__dict__[key] = fix_py2_pickle(p.__dict__[key])

  else:
    # We have no idea...skip conversion (should never be here)
    pass
  return p
