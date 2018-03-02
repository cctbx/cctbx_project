from __future__ import absolute_import, division, print_function

import sys

from fable import simple_compilation

def run(args):
  if args not in (["--with-quotes"], ["--no-quotes"]):
    sys.exit("fable.fem_include_search_paths --with-quotes|--no-quotes")
  comp_env = simple_compilation.environment()
  print(comp_env.assemble_include_search_paths(
    no_quotes=(args[0]=="--no-quotes")))

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
