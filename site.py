#
#   Copyright (C) 2017 Diamond Light Source, Markus Gerstel
#
#   This code is distributed under the BSD license, a copy of which is
#   attached at the end of this file.
#
# 2017-01-24
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

_libtbx = {}

# Print a string if environment variable LIBTBX_IMPORTCACHE is set
if hasattr(os, 'getenv'):
  if os.getenv('LIBTBX_IMPORTCACHE'):
    print("Proof of concept")

# Now hand over to the original python site package
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

# Attempt to find the original python site package in that path
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
