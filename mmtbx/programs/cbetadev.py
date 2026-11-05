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
      apply_phi_psi_correction = False
        .type = bool
        .help = Subtract the mean C-beta deviation for the residue (phi, psi) \
                bin, scaled by a per-model strength, from each deviation. \
                Uses 10 degree bins and separate tables for PRO, \
                VAL/THR/ILE and other residues, with VAL/THR/ILE and others \
                split by whether a PRO follows.
      phi_psi_correction_strength = None
        .type = float(value_min=0)
        .help = Strength of the phi/psi correction, 1 = full. If not set, a \
                strength between 0 and 1 is fitted for each model, and no \
                correction is applied to a model whose mean deviation it \
                would raise. Requires apply_phi_psi_correction.
      phi_psi_correction_region = *all positive
        .type = choice
        .help = Residues the phi/psi correction applies to. all corrects every \
                residue with a table entry. positive corrects only residues \
                in the phi 0 to 150 bins, where the mean deviations share a \
                direction, and fits the strength and checks the mean \
                deviation on those residues alone. It adds far fewer new \
                outliers but also removes far fewer. Requires \
                apply_phi_psi_correction.
      use_parent_for_modified_residues = False
        .type = bool
        .help = Correct modified amino acids (e.g. SEP, HYP, MLY) using the \
                table for their parent amino acid, taken from \
                iotbx.pdb.modified_aa_names. Residues without a single \
                parent, including D-amino acids, are not corrected. \
                Requires apply_phi_psi_correction or mahalanobis_z_score.
      display_phi_psi_correction = False
        .type = bool
        .help = Plot (matplotlib) the original deviation, correction and \
                result for each residue that is an outlier before or after \
                correction. Requires apply_phi_psi_correction.
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
        .help = Plot (matplotlib) a histogram of the C-beta deviations, \
                after any phi/psi correction. With mahalanobis_z_score, \
                plot the Mahalanobis distances d instead.
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
    if (self.params.cbetadev.phi_psi_correction_strength is not None and
        not self.params.cbetadev.apply_phi_psi_correction):
      raise Sorry('phi_psi_correction_strength requires apply_phi_psi_correction')
    if (self.params.cbetadev.phi_psi_correction_region != 'all' and
        not self.params.cbetadev.apply_phi_psi_correction):
      raise Sorry('phi_psi_correction_region requires apply_phi_psi_correction')
    if (self.params.cbetadev.display_phi_psi_correction and
        not self.params.cbetadev.apply_phi_psi_correction):
      raise Sorry('display_phi_psi_correction requires apply_phi_psi_correction')
    if (self.params.cbetadev.use_parent_for_modified_residues and
        not (self.params.cbetadev.apply_phi_psi_correction or
             self.params.cbetadev.mahalanobis_z_score)):
      raise Sorry('use_parent_for_modified_residues requires '
                  'apply_phi_psi_correction or mahalanobis_z_score')

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
      apply_phi_psi_correction=self.params.cbetadev.apply_phi_psi_correction,
      phi_psi_correction_strength=self.params.cbetadev.phi_psi_correction_strength,
      phi_psi_correction_region=self.params.cbetadev.phi_psi_correction_region,
      display_phi_psi_correction=self.params.cbetadev.display_phi_psi_correction,
      use_parent_for_modified_residues=self.params.cbetadev.use_parent_for_modified_residues,
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
