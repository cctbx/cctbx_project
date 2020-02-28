"""
Reads a libtbx_env file and dumps the entire contents.
Some structures e.g. module dictionaries refer to the same object in multiple
places. These will be dumped multiple times. Anything that refers back to a
previous level will show as "__selfreference__". Relocatable paths are just
shown as the regular, joined path (from inspection all base off of the build
path anyway).
Usage:
  libtbx.show_env <libtbx_env>
"""

from __future__ import absolute_import, division, print_function

import os
import pickle
import sys
from pprint import pprint
from types import ModuleType

import six


def _read_obj(obj, prev=None):
  if prev is None:
    prev = []
  if obj in prev:
    return "__selfreference__"
  prev = list(prev) + [obj]

  if isinstance(obj, prop_object):
    dic = {name: _read_obj(val, prev) for name, val in obj.__dict__.items()}
    dic["__type__"] = obj._pathed_type
    return dic
  elif isinstance(obj, list):
    p = []
    for i in obj:
      p.append(_read_obj(i, prev))
    return p
  elif isinstance(obj, dict):
    return {a: _read_obj(b, prev) for a, b in obj.items()}
  else:
    return obj


class prop_object(object):
  """Object that can convert itself to a dictionary"""

  def to_dict(self):
    return _read_obj(self)


def pathed_prop_object(path):
  "Create a class that knows the path it's supposed to represent"

  class _pathed_prop_object(prop_object):
    """Object that can convert itself to a dictionary"""

    _pathed_type = path

  return _pathed_prop_object


class relocatable_path(object):
  def __repr__(self):
    return os.path.normpath(os.path.join(self._anchor._path, self.relocatable))


class absolute_path(object):
  def __repr__(self):
    return self._path


def plainlify(thing):
  if (
    isinstance(thing, six.string_types)
    or isinstance(thing, six.integer_types)
    or isinstance(thing, (float, complex))
  ):
    return thing
  if thing in (None, True, False):
    return thing
  if isinstance(thing, tuple):
    return tuple(map(plainlify, thing))
  if isinstance(thing, list):
    return list(map(plainlify, thing))
  if isinstance(thing, dict):
    return {plainlify(key): plainlify(value) for key, value in thing.items()}
  if isinstance(thing, set):
    return {plainlify(item) for item in thing}
  return str(thing)


def new_module(name, doc=None):
  """Create a new module and inject it into sys.modules"""
  m = ModuleType(name, doc)
  m.__file__ = name + ".py"
  sys.modules[name] = m
  return m


# Create the fake libtbx environment
libtbx = new_module("libtbx")
libtbx.env_config = new_module("libtbx.env_config")
libtbx.path = new_module("libtbx.path")
libtbx.env_config.environment = pathed_prop_object("libtbx.env_config.environment")
libtbx.env_config.build_options = pathed_prop_object("libtbx.env_config.build_options")
libtbx.env_config.module = pathed_prop_object("libtbx.env_config.module")
libtbx.path.relocatable_path = relocatable_path
libtbx.path.absolute_path = absolute_path

# Look for arguments
if "--help" in sys.argv or "-h" in sys.argv:
  print(__doc__)
elif len(sys.argv) != 2:
  print("Usage: libtbx.show_env.py <libtbx_env>")
elif not os.path.isfile(sys.argv[1]):
  print("Error: {} is not a file.".format(sys.argv[1]))
else:
  # Load the environment dump and
  env = pickle.load(open(sys.argv[1], "rb"))
  d = plainlify(env.to_dict())
  pprint(d)
