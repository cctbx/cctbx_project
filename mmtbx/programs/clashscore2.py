"""
All-atom contact analysis.
This is a rewrite of the original clashscore.  This version uses mmtbx.reduce and
mmtbx.probe to generate the contact information rather than stand-alone programs.
It take the same parameters as the original clashscore (except for time_limit)
and it also takes mmtbx.probe parameters.
"""

from __future__ import absolute_import, division, print_function

import os
from mmtbx.validation.clashscore2 import clashscore2
from libtbx.program_template import ProgramTemplate
from datetime import datetime
from mmtbx.probe import Helpers

try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  pass
#from libtbx.utils import Sorry

class Program(ProgramTemplate):
  prog = os.getenv('LIBTBX_DISPATCHER_NAME')
  description="""
%(prog)s file.pdb [params.eff] [options ...]

Options:

  model=input_file          input PDB file
  fast = False              Produce only clashscore number, without anything else
  condensed_probe=False     Run probe with -CON parameter
  keep_hydrogens=False      keep input hydrogen atoms if True, regenerate if False
  nuclear=False             use nuclear x-H distances and vdW radii
  json=False                Outputs results as JSON compatible dictionary
  verbose=True              verbose text output
  b_factor_cutoff=40        B factor cutoff for clash analysis
  do_flips=False            Do flips when adding Hs, overides keep_hydrogens

Example:

  %(prog)s model=1ubq.pdb keep_hydrogens=True
""" % locals()

  master_phil_str = """
    model = None
    .type = path
    .optional = False
    .help = '''input PDB file'''

  fast = False
    .type = bool
    .help = ''' Produce only clashscore number, without anything else'''

  condensed_probe = False
    .type = bool
    .help = ''' Run probe with -CON parameter '''

  json = False
    .type = bool
    .help = "Prints results as JSON format dictionary"

  verbose = True
    .type = bool

  keep_hydrogens = False
    .type = bool
    .help = '''Keep hydrogens in input file'''

  do_flips = False
    .type = bool
    .help = '''Do flips when adding Hsi, overides keep_hydrogens=True'''

  nuclear = False
    .type = bool
    .help = '''Use nuclear hydrogen positions'''

  b_factor_cutoff = None
    .type = int
    .help = '''B factor cutoff for use with MolProbity'''

  clash_cutoff = -0.4
    .type = float
    .help = '''dummy variable for MolProbity, will be removed after MP update'''
""" + Helpers.probe_phil_parameters

# Removed time_limit from the phil parameters
#  time_limit = 120
#    .type = int
#    .help = '''Time limit (sec) for Reduce optimization'''

  datatypes = ['model','phil']
  data_manager_options = ['model_skip_expand_with_mtrix']
  known_article_ids = ['molprobity']

  def validate(self):
    self.data_manager.has_models(raise_sorry=True)

  def run(self, quiet=None): #preserved how quiet was passed to the old run, not sure why
    """
  Calculates nonbonded clashscore using MolProbity (PROBE)

  Returns:
    When verbose=True the function print detailed results to log
    When verbose=False it will print clashscore
  """
    # if do_flips, make keep_hydrogens false
    if self.params.do_flips : self.params.keep_hydrogens = False
    self.info_json = {"model_name":self.data_manager.get_default_model_name(),
                      "time_analyzed": str(datetime.now())}
    self.results = clashscore2(
      self.params.probe,
      data_manager=self.data_manager,
      fast = self.params.fast,
      condensed_probe = self.params.condensed_probe,
      keep_hydrogens=self.params.keep_hydrogens,
      nuclear=self.params.nuclear,
      out=self.logger,
      verbose=self.params.verbose and not quiet,
      b_factor_cutoff=self.params.b_factor_cutoff,
      do_flips=self.params.do_flips)
    if self.params.json:
      print(self.results.as_JSON(self.info_json))
    elif self.params.verbose:
      self.results.show_old_output(out=self.logger)
    else:
      print(round(self.results.get_clashscore(),2), file=self.logger)

  def get_results(self):
    return self.results

  def get_results_as_JSON(self):
    return self.results.as_JSON(self.info_json)
