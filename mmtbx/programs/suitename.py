"""Get suitename for RNA"""
#        Copyright 2021  Richardson Lab at Duke University
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import absolute_import, division, print_function

import os
from mmtbx.suitename.suitealyze import suitealyze
from libtbx.program_template import ProgramTemplate
#from libtbx.utils import Sorry
from datetime import datetime

class Program(ProgramTemplate):
  prog = os.getenv('LIBTBX_DISPATCHER_NAME')
  description="""
  %(prog)s file.pdb [params.eff] [options ...]

Commandline interface for suitename validation
Supply a PDB or mmCIF file; get suitename validation

Options:

  model=input_file      input PDB file
  outliers_only=False   suppress non-outlier results
  report=True           standard text output
  json=False            JSON format output
  string=False          output in suitestring format, 3 characters per suite
  kinemage=False        kinemage-format visual markup

Example:

  %(prog)s model=1ubq.pdb
""" % locals()

  master_phil_str = """
  include scope mmtbx.validation.molprobity_cmdline_phil_str
    suitename {
    # input
      infile=""
        .type=str
        .help="the file to process"
      anglefields = 9
        .type=int
        .help="number of angle fields provided, for textual input only"
      pointidfields = 7
        .type=int
        .help="number of point id fields before the angle fields"
      ptid=0
        .type=int
        .help="number of point id fields before the angle fields"
      residuein=false
        .type=bool
        .help="expect dangle format giving residues"
      suitein=false
        .type=bool
        .help="expect kinemage format giving suites directly"
    # output
      string=False
        .type=bool
        .help="output in string format, 3 characters per suite"
      json=False
        .type=bool
        .help="output in JSON format, useful for machine parsing"
      kinemage=False
        .type=bool
        .help="output in kinemage format, useful for visualization"
      markup=False
        .type=bool
        .help="visual markup of suites and outliers for kinemage format"
      report=true
        .type=bool
        .help="output as a report, giving statistical details"
      chart=False
        .type=bool
        .help="modifier to standard report, output without statistical summary"
      nosequence = False
        .type=bool
        .help="modifier to string format, do not include base letters"
      causes=False
        .type=bool
        .help="output extra details concerning the causes of each assignment made"
      test=False
        .type=bool
        .help="display a lat of additional information about program internals"
    # compute
      satellites=False
        .type=bool
        .help="use the special satelliteWidths values for satellites"
      nowannabe=False
        .type=bool
        .help="do not consider 'wannabe' clusters"
      noinc=False
        .type=bool
        .help="do not display incomplete suites"
      etatheta=False
        .type=bool
      altid="A"
        .type=str
        .help="which alternate conformer to use (A, B, etc)"
      altidfield = 6
        .type=int
        .help="which field (1-based) gives the alternate conformer code"
      version=false
        .type=bool
        .help="give the version number of suite name"
    # deprecated and automatically true:
      oneline=false
        .type=bool
    }
"""
  datatypes = ['model','phil']
  data_manager_options = ['model_skip_expand_with_mtrix']
  known_article_ids = ['molprobity']

  def validate(self):
    self.data_manager.has_models(raise_sorry=True)

  def run(self):
    hierarchy = self.data_manager.get_model().get_hierarchy()
    self.info_json = {"model_name":self.data_manager.get_default_model_name(),
                      "time_analyzed": str(datetime.now())}
    self.suite_results = suitealyze(pdb_hierarchy=hierarchy, options=self.params)
    if self.params.suitename.string or self.params.suitename.oneline:# == "string":
      self.suite_results.display_suitestrings(blockform=True)
    elif self.params.suitename.markup:# == "markup":
      print(self.suite_results.as_kinemage_markup())
    elif self.params.suitename.json:
      print(self.suite_results.as_JSON())
    elif self.params.suitename.report:# == "report":
      self.suite_results.show_old_output(verbose=True)
    else:
      self.suite_results.show_old_output(verbose=True)
    #print(suite_results.as_kinemage())
    #suite_results.show_summary()

  def get_results(self):
    return self.suite_results

  def get_results_as_JSON(self):
    return self.suite_results.as_JSON(self.info_json)

    #hierarchy.atoms().reset_i_seq()
    #result = cbetadev(
    #  pdb_hierarchy=hierarchy,
    #  outliers_only=self.params.outliers_only,
    #  apply_phi_psi_correction=self.params.cbetadev.apply_phi_psi_correction,
    #  display_phi_psi_correction=self.params.cbetadev.display_phi_psi_correction,
    #  exclude_d_peptides=self.params.cbetadev.exclude_d_peptides,
    #  out=self.logger,
    #  quiet=False)
    #if self.params.cbetadev.output == "kin":
    #  self.logger.write(result.as_kinemage())
    #elif self.params.cbetadev.output == "bullseye":
    #  filebase = os.path.basename(self.data_manager.get_model_names()[0])
    #  self.logger.write(result.as_bullseye_kinemage(pdbid=filebase))
    #elif self.params.verbose:
    #  #pdb_file_str = os.path.basename(self.params.model)[:-4]
    #  #get input file name from data manager, strip file extension
    #  pdb_file_str = os.path.basename(self.data_manager.get_model_names()[0])[:-4]
    #  result.show_old_output(out=self.logger, prefix=pdb_file_str, verbose=True)
