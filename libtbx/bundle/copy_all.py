from __future__ import absolute_import, division, print_function
from libtbx.bundle import copy_runtime_sources
from libtbx.bundle import copy_build_libtbx
import sys, os

def run(prefix, to_skip_src=[]):
  copy_runtime_sources.run(prefix+"_sources", to_skip=to_skip_src)
  copy_build_libtbx.run(prefix+"_build")

def options_to_list(opts):
  rv = []
  toks = opts.split(';')
  for t in toks:
    if not t: continue
    rv.append(os.sep.join([".", *t.split('/')]))
  return rv

if (__name__ == "__main__"):
  from optparse import OptionParser
  parser = OptionParser()
  parser.add_option('--skip_src',
        dest='skip_src',
        default="",
        help='Skip paths from the bundle, Usage: skip_src="langchain;maptbx/bcr"')
  options, argv = parser.parse_args()
  assert len(argv) == 1
  to_skip_src = options_to_list(options.skip_src)
  run(argv[0], to_skip_src=options_to_list(options.skip_src))
