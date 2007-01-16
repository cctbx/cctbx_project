#!/usr/bin/env BOOST_ADAPTBX_DOCSTRING_OPTIONS=show_user_defined=True,show_signatures=False ../../../cctbx_build/bin/cctbx.python

"""
Synopsis

epydoc --conf=<path/to/configuration/file>

Description

Generate the epydoc documentation specified by the given configuration file.

Instructions for documentation writers

When a pure Python class A uses the Boost.Python wrapping of a C++ class B, 
the docstring of A should feature a link to the doxygen-generated
documentation of B. That link shall be written as e.g.
  U{least_squares_residual<CCTBX_DOXYGEN_DOC_ROOT/
  classcctbx_1_1xray_1_1targets_1_1least__squares__residual.html>}
(white spaces are irrelevant in the URL part between <> which can be taken
advantage of to break long lines nicely).

This relies on two requirements on the part of the C++ documentation writers:
(a) the doxygen configuration file is assumed to be at ../dox/Doxyfile w.r.t
this script;
(b) the C++ documentation root directory and the Python documentation root
directory will be in the same directory on the server.
"""

import sys, os.path
if os.path.exists(os.path.join(sys.path[0], 'epydoc.py')): 
  del sys.path[0]

import epydoc.cli
options, names = epydoc.cli.parse_arguments()
epydoc.cli.main(options, names)
output_dir = options.target
doxyfile = open('../dox/Doxyfile')
import re
html_output_pat = re.compile(r'^HTML_OUTPUT\s*=\s*(.*)$')
for li in doxyfile:
  m = html_output_pat.search(li)
  if m is not None:
    cpp_doc = "../" + m.group(1) + "/"
    cpp_doc = cpp_doc.replace('//', '/')
    break
import os, os.path, fnmatch
for root, dirs, files in os.walk(output_dir):
  for f in [ os.path.join(root,f) for f in fnmatch.filter(files, '*.html') ]:
    tmp = '%s.tmp' % f
    fin = open(f)
    fout = open(tmp, 'w')
    for li in fin:
      fout.write(li.replace('http://CCTBX_DOXYGEN_DOC_ROOT/', cpp_doc))
    os.rename(tmp, f)

    
