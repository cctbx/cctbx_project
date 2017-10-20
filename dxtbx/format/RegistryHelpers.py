#!/usr/bin/env python
# RegistryHelpers.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Things to help the ImageFormat registry to work.

from __future__ import absolute_import, division

import os
import sys
import imp
import exceptions
import traceback

def LoadFormatClasses():
  '''Look for files named Format(something).py in the sensible
  places (i.e. in the xia2 distribution and in the users home area)
  and import the corresponding modules using their fully qualified
  names.'''

  import dxtbx.format

  # FIXME in here - dxtbx should already be in os.path - look for it there,
  # also wouldn't it be tidy to refer to a Phil parameter?

  format_dir = os.path.split(dxtbx.format.__file__)[0]

  home = os.curdir
  if os.name == 'nt' and \
     'HOMEDRIVE' in os.environ and \
     'HOMEPATH' in os.environ:
    home = os.path.join(os.environ['HOMEDRIVE'],
                        os.environ['HOMEPATH'])
  elif 'HOME' in os.environ:
    home = os.environ['HOME']

  for f in os.listdir(format_dir):
    if 'Format' in f[:6] and '.py' in f[-3:]:
      name = f[:-3]
      fqname = dxtbx.format.__name__ + '.' + name
      _LoadFormatModule(name, fqname, format_dir)

  format_dir = os.path.join(home, '.dxtbx')
  if os.path.exists(format_dir):
    if format_dir not in sys.path:
      sys.path.append(format_dir)
    for f in os.listdir(format_dir):
      if 'Format' in f[:6] and '.py' in f[-3:]:
        name = f[:-3]
        _LoadFormatModule(name, name, format_dir)

def _LoadFormatModule(name, fqname, path):
  '''Load a format class module, which will trigger the automated
  self-registration.  This module will therefore not need to publish
  anything as the module will self-publish.  The idea being that
  these format classes were found by the search procedure above.  On
  success, the function returns the loaded module, otherwise it
  returns None.'''

  # Early return if module already imported.
  try:
    return sys.modules[fqname]
  except KeyError:
    pass

  try:
    stream, pathname, description = imp.find_module(name, [path])
  except ImportError:
    return None

  try:
    module = imp.load_module(fqname, stream, pathname, description)
  except exceptions.Exception:
    traceback.print_exc(sys.stderr)
  finally:
    if stream:
      stream.close()

  return module

if __name__ == '__main__':

  LoadFormatClasses()
