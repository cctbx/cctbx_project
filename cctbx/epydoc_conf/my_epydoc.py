#!/usr/bin/env BOOST_ADAPTBX_DOCSTRING_OPTIONS=show_user_defined=True,show_signatures=False ../../../cctbx_build/bin/cctbx.python

import sys, os.path
if os.path.exists(os.path.join(sys.path[0], 'epydoc.py')): 
  del sys.path[0]

import epydoc.cli
epydoc.cli.cli()
options, names = epydoc.cli.parse_arguments()
output_dir = options.target
import os, os.path, fnmatch
for root, dirs, files in os.walk(output_dir):
  for f in [ os.path.join(root,f) for f in fnmatch.filter(files, '*.html') ]:
    tmp = '%s.tmp' % f
    fin = open(f)
    fout = open(tmp, 'w')
    for li in fin:
      fout.write(li.replace('http://EPYDOC_RELATIVE_URL/', ''))
    os.rename(tmp, f)

    
