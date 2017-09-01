#!/usr/bin/env python
#
# dxtbx.serialize.imageset.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import absolute_import, division
from contextlib import contextmanager

@contextmanager
def temp_chdir(path):
  ''' A context manager to temporarily change the current directory, perform a
  task and then change back to the previous working directory. '''
  from os import getcwd, chdir
  cwd = getcwd()
  try:
    chdir(path)
    yield
  finally:
    chdir(cwd)

def load_path(path):
  ''' Load a filename from a JSON file.

  First expand any environment and user variables. Then create the absolute path
  from the current directory (i.e. the directory in which the JSON file is
  located.

  Params:
    path The path to the file.

  '''
  from os.path import abspath, expanduser, expandvars
  if path is None or path == "":
    return ""
  return abspath(expanduser(expandvars(path)))
