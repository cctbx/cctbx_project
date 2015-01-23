# LIBTBX_SET_DISPATCHER_NAME phenix.clashscore

from __future__ import division
from mmtbx.validation.clashscore import clashscore
from libtbx.utils import Usage
import iotbx.phil
import sys

master_phil_str = """
  model = None
    .type = path
    .optional = False
    .help = '''input PDB file'''

  verbose = True
    .type = bool

  keep_hydrogens = True
    .type = bool
    .help = '''Keep hydrogens in input file'''

  nuclear = False
    .type = bool
    .help = '''Use nuclear hydrogen positions'''

  time_limit = 120
    .type = int
    .help = '''Time limit (sec) for Reduce optimization'''

  b_factor_cutoff = None
    .type = int
    .help = '''B factor cutoff for use with MolProbity'''
"""

usage_string = """\
phenix.clashscore file.pdb [params.eff] [options ...]

Options:

  model=input_file          input PDB file
  keep_hydrogens=True       keep input hydrogen files (otherwise regenerate)
  nuclear=False             use nuclear x-H distances and vdW radii
  verbose=True              verbose text output
  b_factor_cutoff=40        B factor cutoff for clash analysis

Example:

  phenix.clashscore model=1ubq.pdb keep_hydrogens=True
"""

def run (args, out=sys.stdout, quiet=None) :
  """
  Calculates nonbonded clashscore using MolProbity (PROBE)

  Returns:
    When verbose=True the function print detailed results to log
    When verbose=False it will print clashscore
  """
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil_string=master_phil_str,
    pdb_file_def="model",
    usage_string=usage_string)
  params = cmdline.work.extract()
  if (params.model is None) :
    raise Usage(usage_string)

  pdb_in = cmdline.get_file(params.model, force_type="pdb")
  hierarchy = pdb_in.file_object.hierarchy
  result = clashscore(
    pdb_hierarchy=hierarchy,
    keep_hydrogens=params.keep_hydrogens,
    nuclear=params.nuclear,
    out=out,
    verbose=params.verbose and not quiet,
    b_factor_cutoff=params.b_factor_cutoff)
  if params.verbose:
    result.show_old_output(out=out)
  else:
    print >> out,round(result.get_clashscore(),2)

if (__name__ == "__main__") :
  run(sys.argv[1:])
