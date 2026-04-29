"""Run CaBLAM on a model"""
from __future__ import absolute_import, division, print_function

import os
from mmtbx.validation.cablam import cablamalyze
from datetime import datetime
from libtbx.program_template import ProgramTemplate
try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  pass
#from libtbx.utils import Sorry

class Program(ProgramTemplate):
  prog = os.getenv('LIBTBX_DISPATCHER_NAME')
  description="""
  %(prog)s file.pdb [options ...]

Options:

  output=          text : default output.  Prints machine-readable
                          columnated and colon-separated validation text to
                          screen.
                        json : prints results as JSON compatible dictionary
                        kin : prints kinemage markup for validation to screen
                        full_kin : prints kinemage markup and struture kinamge
                          to screen
                        points_kin : prints point cloud of residues in cablam
                          space in kinemage format
                        records : prints pdb-style HELIX and SHEET records to
                          screen, based on CaBLAM's identification of secondary
                          structure
                        records_and_pdb : prints pdb-style HELIX and SHEET
                          records to screen, followed by PDB file formatted
                          coordinates for the input structure
                        oneline : prints single-line summary of CaBLAM
                          validation statistics

  outliers_only=False   compresses certain outputs (text) to show only outlier
                          residues



Example:

  %(prog)s model=1ubq.pdb outliers_only=True
""" % locals()

  master_phil_str = """
  include scope mmtbx.validation.molprobity_cmdline_phil_str
  json = False
    .type = bool
    .help = "Prints results as JSON format dictionary"
  pdb_infile = None
    .type = path
    .help = input PDB file
  output_type = *text kin full_kin points_kin records records_and_pdb oneline
    .type = choice
    .help = choose output type: \
     =text for default colon-separated residue-by-residue validation \
     =kin for outlier markup in kinemage format \
     =full_kin for outlier markup appended to structure - opens in KiNG \
     =points_kin for pointcloud in cablam space \
     =records for PDB-style HELIX/SHEET records \
     =records_and_pdb for PDB-style HELIX/SHEET records attached to a PDB file \
     =oneline for a one-line structure summary
"""
  datatypes = ['model','phil']
  data_manager_options = ['model_skip_expand_with_mtrix']
  known_article_ids = ['molprobity']

#TODO: get cablam.interpretation() to print again

  def validate(self):
    self.data_manager.has_models(raise_sorry=True)

  def run(self):
    hierarchy = self.data_manager.get_model().get_hierarchy()
    self.info_json = {"model_name":self.data_manager.get_default_model_name(),
                      "time_analyzed": str(datetime.now())}
    self.cablam = cablamalyze(
      pdb_hierarchy=hierarchy,
      outliers_only=self.params.outliers_only,
      out=self.logger,
      quiet=False)

    #output_type = *text kin full_kin points_kin records records_and_pdb oneline
    if self.params.json:
      print(self.get_results_as_JSON())
    elif self.params.output_type=='oneline':
      pdb_file_str = os.path.basename(self.data_manager.get_model_names()[0])
      self.cablam.as_oneline(pdbid=pdb_file_str)
    elif self.params.output_type=='kin':
      self.cablam.as_kinemage()
    elif self.params.output_type=='full_kin':
      self.cablam.as_full_kinemage(pdbid=pdbid)
    elif self.params.output_type=='points_kin':
      self.cablam.as_pointcloud_kinemage()
    elif self.params.output_type=='records':
      self.cablam.as_records()
    elif self.params.output_type=='records_and_pdb':
      self.cablam.as_records_and_pdb()
    else: #default text output
      self.cablam.as_text(outliers_only=self.params.outliers_only)

  def get_results(self):
    return self.cablam

  def get_results_as_JSON(self):
    return self.cablam.as_JSON(self.info_json)
