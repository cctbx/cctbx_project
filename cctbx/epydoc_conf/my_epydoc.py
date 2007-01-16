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
  U{DOXYCLASS:scitbx::lbfgs::minimizer}
and will give 
<a href="../c_plus_plus/classscitbx_1_1lbfgs_1_1drop__convergence__test.html">
class scitbx::lbfgs::minimizer</a>
The other magic keyword is DOXYSTRUCT for struct instead of class.

This relies on 3 requirements on the part of the C++ documentation writers:
(a) the doxygen configuration file is assumed to be at ../dox/Doxyfile w.r.t
this script;
(b) the C++ documentation root directory is given by the HTML_OUTPUT
configuration entry in the doxyfile (c_plus_plus  in our example above);
(c) the C++ documentation root directory and the Python documentation root
directory will be in the same directory on the server.
"""

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
    cpp_doc = "../" + m.group(1)
    cpp_doc = cpp_doc.replace('//', '/')
    break
import os, os.path, fnmatch
doxypat = re.compile(r'(href=")? DOXY ( CLASS | STRUCT ) : ([:\w]+)',
                   re.VERBOSE)
def repl(m):
  klass = m.group(2).lower()
  name = m.group(3)
  href = m.group(1)
  if href:
    return "%s%s/%s%s.html" % (
      href,
      cpp_doc, 
      klass,
      name.replace('_', '__').replace('::', '_1_1'),
    )
  else:
    return "%s %s" % (klass, name)
    
for root, dirs, files in os.walk(output_dir):
  for f in [ os.path.join(root,f) for f in fnmatch.filter(files, '*.html') ]:
    tmp = '%s.tmp' % f
    fin = open(f)
    fout = open(tmp, 'w')
    for li in fin: 
      li = doxypat.sub(repl, li)
      fout.write(li)
    os.rename(tmp, f)
