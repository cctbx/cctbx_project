#
#   Copyright (C) 2017 Diamond Light Source, Markus Gerstel
#
#   This code is distributed under the BSD license, a copy of which is
#   attached at the end of this file.
#
# What does the code do?
# ----------------------
#
# This is a proof-of-concept test.
# In other words: ideally you should not notice that this code is run.
#
# The aim of the exercise is to fundamentally subvert the python import
# mechanism to speed up loading times of python commands over networked file
# systems. The python 'import' command used to load packages is rather
# inefficient, as it tries to locate each package in a long list of potential
# places under lots of different names using a trial-and-error strategy.
# The code below instead predicts exactly where the requested package lives,
# and immediately points python to the correct place.
#
#
# How do I use this code?
# -----------------------
#
# For now you have to set the environment variable
#   LIBTBX_IMPORTCACHE
# to enable the predictive caching mechanism.
# If you want to see a running commentary of what is happening, set
#   LIBTBX_IMPORTCACHE_DEBUG
# or
#   LIBTBX_IMPORTCACHE_TIME
# as well.
#
#
# Why is this file imported? Why is this file in this particular place?
# ---------------------------------------------------------------------
#
# This file overrides the python 'site' package. It will therefore always be
# imported by python (unless the python interpreter is invoked with '-S'), and
# before any other imports. When running python from a cctbx dispatcher the
# first item in sys.path is the cctbx_project directory. Therefore this file is
# pretty much the very first python file loaded when starting the interpreter.
#

# The following imports are all libraries that will be loaded by python in any
# case, so unconditionally importing them comes with little additional cost.
import imp
import os
import sys

_libtbx = {}

# Enable predictive import caching only if environment variable
# LIBTBX_IMPORTCACHE is set.
pkgutil = None
if hasattr(os, 'getenv') and os.getenv('LIBTBX_IMPORTCACHE'):
  try:
    import pkgutil
  except ImportError:
    pkgutil = None

if pkgutil:
  class _PredictiveImport(object):
    import_path = {}
    '''Dictionary of examined sys.path elements.
       Values are 'False' if the path element is not a valid path.
       Otherwise values are dictionaries of lists of found suffix types.'''
    loader_path = {}
    '''Dictionary of examined sys.path elements.
       Value is a special loader if the path element is handled by this.'''

    def __init__(self):
      if os.getenv('LIBTBX_IMPORTCACHE_DEBUG'):
        self.debug = self.debug_print
      if os.getenv('LIBTBX_IMPORTCACHE_TIME'):
        import time
        self.time = time.time
        self.start = time.time()
        self.debug = self.debug_print_time
      self.debug('libtbx.importcache enabled')

    def debug(self, message):
      pass

    def debug_print(self, message):
      print("predictive import: " + message)

    def debug_print_time(self, message):
      print("predictive import %8.3f: " % (self.time() - self.start) + message)

    def find_module(self, fullname, path=None):
      if '.' in fullname:
        # Caching deals only with top level packages
        rdot = fullname.rindex('.')
        head = fullname[:rdot] # 'logging'
        tail = fullname[rdot+1:] # 'os'
        #self.debug("ignoring call for " + fullname + " (" + str(path) + ")")
        self.debug("trying to predict location of " + fullname + " within " + str(path))
        for directory in path:
          solution = self.find_in_directory(tail, fullname, directory)
          if solution:
            return solution
        self.debug("no such luck")
        raise ImportError('No module named ' + fullname)
        return None

      self.debug("trying to predict location of " + fullname)
      for directory in sys.path:
        solution = self.find_in_directory(fullname, fullname, directory)
        if solution:
          return solution

      # No match found in paths. Check if this is a builtin module
      if imp.is_builtin(fullname):
        self.debug("identified " + fullname + " as builtin")
        solution = (None, fullname, ('', '', imp.C_BUILTIN))
        return pkgutil.ImpLoader(fullname, *solution)

      # TODO: Frozen?

      # Still no match. This might not exist
      self.debug("Import " + fullname + " failed: not found")
      raise ImportError('No module named ' + fullname)

    def find_in_directory(self, module, fullname, directory):
      if directory == '': directory = '.'

      if directory not in self.import_path:
        # If directory has not been seen before gather information about it
        self.examine_path(directory)

      if directory in self.loader_path:
        # Path element has a special handler. Ask for module.
        self.debug("checking " + fullname + " with " + str(self.loader_path[directory]))
        module_loader = self.loader_path[directory].find_module(module)
        if module_loader:
          self.debug("module found.")
          return module_loader

      if self.import_path[directory]:
        # Path element is a regular directory. Try finding the module in here.
        if self.import_path[directory].get(module):
          # Match found.
          if imp.PKG_DIRECTORY in self.import_path[directory][module]:
            full_pathname = os.path.join(directory, module)
            if os.path.isdir(full_pathname) and os.path.exists(os.path.join(full_pathname, '__init__.py')):
              solution = (None, full_pathname, ('', '', imp.PKG_DIRECTORY))
              self.debug("predicted solution: " + str(solution))
              return pkgutil.ImpLoader(fullname, *solution)
          for suffix in imp.get_suffixes():
            if suffix[0] in self.import_path[directory][module]:
              filename = os.path.join(directory, module + suffix[0])
              fh = open(filename, suffix[1])
              solution = (fh, filename, suffix)
              self.debug("predicted solution: " + str(solution))
              return pkgutil.ImpLoader(fullname, *solution)

    def examine_path(self, path):
      self.import_path[path] = False
      if not os.path.exists(path):
        self.debug("Ignoring new path " + path + ": Does not exist")
        return

      self.debug("Examining new path:" + path)
      if path not in sys.path_importer_cache:
        pkgutil.get_importer(path)
      external_importer = sys.path_importer_cache.get(path)
      if external_importer:
        self.debug(path + " is served by external importer " + str(external_importer))
        self.loader_path[path] = external_importer
        return

      if not os.path.isdir(path):
        self.debug("Disabling path " + path + ": Not a directory")
        return

      self.import_path[path] = {}
      suffixes = [suffix[0] for suffix in imp.get_suffixes()]
      self.debug("Found in " + path + ": " + str(os.listdir(path)))
      for entry in os.listdir(path):
        self.import_path[path][entry] = [ imp.PKG_DIRECTORY ] + self.import_path[path].get(entry, [])
        for suffix in suffixes:
          if entry.endswith(suffix):
            basename = entry[:-len(suffix)]
            self.import_path[path][basename] = [ suffix ] + self.import_path[path].get(basename, [])
            self.debug("Found in " + path + ": " + basename + " (by extension)")

  _libtbx_cache = _PredictiveImport()
  sys.meta_path.append(_libtbx_cache)


# Now hand over to the original python 'site' package
# Find the tail of sys.path not including the directory of this file
_libtbx['path'] = []
try:
  _libtbx['this_path'] = __file__
except NameError:
  _libtbx['this_path'] = '.' # May not be defined in all cases
_libtbx['this_path'] = os.path.abspath(os.path.dirname(_libtbx['this_path']))
if sys.hexversion >= 0x02020000:
  _libtbx['this_path'] = os.path.realpath(_libtbx['this_path'])
for _libtbx['path_candidate'] in sys.path:
  _libtbx['path_candidate'] = os.path.abspath(_libtbx['path_candidate'])
  if sys.hexversion >= 0x02020000:
    _libtbx['path_candidate'] = os.path.realpath(_libtbx['path_candidate'])
  if _libtbx['path_candidate'] == _libtbx['this_path']:
    _libtbx['path'] = []
  else:
    _libtbx['path'].append(_libtbx['path_candidate'])

# Attempt to find the original python 'site' package in that path
_libtbx['true_site'] = []
if _libtbx['path']:
  try:
    _libtbx['true_site'] = imp.find_module('site', _libtbx['path'])
  except ImportError:
    pass # Given that site should be in the python directory this
         # should not fail. May however only be true for cPython.
if _libtbx['true_site']:
  _libtbx['true_site_path'] = os.path.abspath(os.path.dirname(_libtbx['true_site'][1]))
  if sys.hexversion >= 0x02020000:
    _libtbx['true_site_path'] = os.path.realpath(_libtbx['true_site_path'])
  if _libtbx['true_site_path'] == _libtbx['this_path']:
    print("Error in site.py: Could not find original python site package")
  else:
    # Load the original python site package in place
    sys.modules['site'] = imp.load_module('site', *_libtbx['true_site'])
    __file__ = sys.modules['site'].__file__

# The end.

# Copyright (c) 2017 Diamond Light Source, Lawrence Berkeley National Laboratory
# and the Science and Technology Facilities Council, All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
#
# Redistributions in binary form must reproduce the above copyright notice, this
# list of conditions and the following disclaimer in the documentation and/or
# other materials provided with the distribution.
#
# Neither the name of the Diamond Light Source, Lawrence Berkeley National
# Laboratory or the Science and Technology Facilities Council, nor the names of
# its contributors may be used to endorse or promote products derived from this
# software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
