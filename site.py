#
# 2017-01-18
#
# This is a proof-of-concept test.
#
# This file overrides the python 'site' package. It will therefore always be
# imported by python (unless the python interpreter is invoked with '-S'), and
# before any other imports. For this test, all this code does is identify the
# location of the 'site' package that otherwise would have been loaded, and
# load this one instead.
#
# In other words: ideally you should not notice that this code is run.
#
# The aim of the exercise is to eventually include code here that can
# fundamentally subvert the python import mechanism in order to speed up
# loading times, particularly over networked file systems.
#


# The following imports are all libraries that will be loaded by python in any
# case, so they come with little additional cost.
import imp
import os
import sys

# Print a string if environment variable LIBTBX_IMPORT_CACHEDB is set
try:
  _libtbx_db = os.getenv('LIBTBX_IMPORT_CACHEDB')
  if _libtbx_db:
    print("Proof of concept")
except (ImportError, AttributeError):
  pass

# Now hand over to the original python site package
# Find a subset of sys.path that does not include the directory of this file
_path_remainder = []
for _path in sys.path:
  if __file__.startswith(_path):
    _path_remainder = []
  else:
    _path_remainder.append(_path)

# Attempt to find the original python site package in that path
_original_site_package = imp.find_module('site', _path_remainder)
if _original_site_package[1] == __file__:
  print("Error in site.py: Could not find original python site package")
else:
  # Load the original python site package in place
  sys.modules['site'] = imp.load_module('site', *_original_site_package)
  __file__ = sys.modules['site'].__file__

# Clean up
del(_original_site_package)
del(_path_remainder)
