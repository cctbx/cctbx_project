"""Calculates ideal CB atom position based on mainchain geometry,
then displays deviation of modeled CB from that ideal positon.
"""
from __future__ import absolute_import, division, print_function

import os
from mmtbx.validation.cbetadev import cbetadev
from libtbx.program_template import ProgramTemplate
from libtbx.utils import Sorry
from datetime import datetime

class Program(ProgramTemplate):
  prog = os.getenv('LIBTBX_DISPATCHER_NAME')
  description="""
  %(prog)s file.pdb [params.eff] [options ...]

Calculates ideal CB atom position based on mainchain geometry,
then displays deviation of modeled CB from that ideal positon.
Deviations of >= 0.25A are considered outliers

Options:

  model=input_file      input PDB file
  output=text, kin, or bullseye    select type of output
  json=False            Outputs results as JSON compatible dictionary
  outliers_only=False   suppress non-outlier results

Example:

  %(prog)s model=1ubq.pdb
""" % locals()

  master_phil_str = """
  include scope mmtbx.validation.molprobity_cmdline_phil_str
    json = False
      .type = bool
      .help = Prints results as JSON format dictionary
    cbetadev {
      output = *text kin bullseye
        .type = choice
        .help = choose output type
      apply_positive_phi_correction = False
        .type = bool
        .help = Subtract the mean C-beta deviation for the residue (phi, psi) \
                bin, but only for residues with 0 < phi < 150. Cannot be \
                used with apply_phi_psi_correction.
      apply_phi_psi_correction = False
        .type = bool
        .help = Subtract the mean C-beta deviation for the residue (phi, psi) \
                bin from each deviation. Uses 10 degree bins and separate \
                tables for PRO, VAL/THR/ILE and other residues. Cannot be \
                used with apply_positive_phi_correction.
      use_parent_for_modified_residues = False
        .type = bool
        .help = Correct modified amino acids (e.g. SEP, HYP, MLY) using the \
                table for their parent amino acid, taken from \
                iotbx.pdb.modified_aa_names. Residues without a single \
                parent, including D-amino acids, are not corrected. \
                Requires apply_phi_psi_correction, \
                apply_positive_phi_correction or mahalanobis_z_score.
      fill_empty_phi_psi_bins = False
        .type = bool
        .help = For a residue with 0 < phi < 150 whose (phi, psi) bin has \
                no correction, use the mean correction of the populated \
                neighbouring bins, if at least 3 of the 8 are populated. \
                Requires apply_phi_psi_correction or \
                apply_positive_phi_correction.
      display_phi_psi_correction = False
        .type = bool
        .help = Plot (matplotlib) the original deviation, correction and \
                result for each residue that is an outlier before or after \
                correction. Requires apply_phi_psi_correction or \
                apply_positive_phi_correction.
      mahalanobis_z_score = False
        .type = bool
        .style = hidden
        .help = Also score each uncorrected deviation vector by its \
                Mahalanobis distance d from the reference distribution for \
                its (phi, psi) bin. An isotropic model error s is fitted \
                for each model so that the median d squared matches \
                chi-squared with two degrees of freedom. The 0.25 Angstrom \
                outliers are unchanged.
      mahalanobis_cutoff = 4.0
        .type = float
        .style = hidden
        .help = Residues with d at or above this value are counted as \
                Mahalanobis outliers. Used with mahalanobis_z_score.
      show_deviation_stats = False
        .type = bool
        .style = hidden
        .help = Print descriptive statistics of the C-beta deviations \
                projected onto the dihedralNABB = 45 degree axis.
      plot_deviation_histogram = False
        .type = bool
        .style = hidden
        .help = Plot (matplotlib) a histogram of the C-beta deviations \
                projected onto the dihedralNABB = 45 degree axis.
      exclude_d_peptides = False
        .type = bool
        .style = hidden
        .help = Attempts to exclude D-peptide using the large CBD
      }
"""
  datatypes = ['model','phil']
  data_manager_options = ['model_skip_expand_with_mtrix']
  known_article_ids = ['molprobity']

  def validate(self):
    self.data_manager.has_models(raise_sorry=True)
    if (self.params.cbetadev.apply_phi_psi_correction and
        self.params.cbetadev.apply_positive_phi_correction):
      raise Sorry('Please use either apply_phi_psi_correction or '
                  'apply_positive_phi_correction')
    if (self.params.cbetadev.display_phi_psi_correction and
        not (self.params.cbetadev.apply_phi_psi_correction or
             self.params.cbetadev.apply_positive_phi_correction)):
      raise Sorry('display_phi_psi_correction requires '
                  'apply_phi_psi_correction or apply_positive_phi_correction')
    if (self.params.cbetadev.use_parent_for_modified_residues and
        not (self.params.cbetadev.apply_phi_psi_correction or
             self.params.cbetadev.apply_positive_phi_correction or
             self.params.cbetadev.mahalanobis_z_score)):
      raise Sorry('use_parent_for_modified_residues requires '
                  'apply_phi_psi_correction, apply_positive_phi_correction '
                  'or mahalanobis_z_score')
    if (self.params.cbetadev.fill_empty_phi_psi_bins and
        not (self.params.cbetadev.apply_phi_psi_correction or
             self.params.cbetadev.apply_positive_phi_correction)):
      raise Sorry('fill_empty_phi_psi_bins requires '
                  'apply_phi_psi_correction or apply_positive_phi_correction')

  def get_results(self):
    return self.results

  def get_results_as_JSON(self):
    return self.results.as_JSON(self.info_json)

  def run(self):
    hierarchy = self.data_manager.get_model().get_hierarchy()
    self.info_json = {"model_name":self.data_manager.get_default_model_name(),
                      "time_analyzed": str(datetime.now())}
    self.results = cbetadev(
      pdb_hierarchy=hierarchy,
      outliers_only=self.params.outliers_only,
      apply_positive_phi_correction=self.params.cbetadev.apply_positive_phi_correction,
      apply_phi_psi_correction=self.params.cbetadev.apply_phi_psi_correction,
      display_phi_psi_correction=self.params.cbetadev.display_phi_psi_correction,
      use_parent_for_modified_residues=self.params.cbetadev.use_parent_for_modified_residues,
      fill_empty_phi_psi_bins=self.params.cbetadev.fill_empty_phi_psi_bins,
      exclude_d_peptides=self.params.cbetadev.exclude_d_peptides,
      mahalanobis_z_score=self.params.cbetadev.mahalanobis_z_score,
      mahalanobis_cutoff=self.params.cbetadev.mahalanobis_cutoff,
      show_deviation_stats=self.params.cbetadev.show_deviation_stats,
      plot_deviation_histogram=self.params.cbetadev.plot_deviation_histogram,
      out=self.logger,
      quiet=False)
    if self.params.cbetadev.output == "kin":
      self.logger.write(self.results.as_kinemage())
    elif self.params.cbetadev.output == "bullseye":
      filebase = os.path.basename(self.data_manager.get_model_names()[0])
      self.logger.write(self.results.as_bullseye_kinemage(pdbid=filebase))
    elif self.params.json:
      print(self.get_results_as_JSON())
    elif self.params.verbose:
      #pdb_file_str = os.path.basename(self.params.model)[:-4]
      #get input file name from data manager, strip file extension
      pdb_file_str = os.path.basename(self.data_manager.get_model_names()[0])[:-4]
      self.results.show_old_output(out=self.logger, prefix=pdb_file_str, verbose=True)
