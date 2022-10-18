from __future__ import absolute_import, division, print_function
import sys

def parse_input():
  from iotbx.phil import parse
  master_phil="""
    context = kokkos_gpu *cuda
      .type = choice
      .optional = False
      .help = backend for parallel execution
  """
  phil_scope = parse(master_phil)
  # The script usage
  import libtbx.load_env # implicit import
  from dials.util.options import ArgumentParser
  # Create the parser
  parser = ArgumentParser(
        usage="\n libtbx.python %s context=[kokkos_gpu|cuda]"%(sys.argv[0]),
        phil=phil_scope,
        epilog="test exafel API, cuda vs. kokkos")
  # Parse the command line. quick_parse is required for MPI compatibility
  params, options = parser.parse_args(show_diff_phil=True,quick_parse=True)
  return params,options
