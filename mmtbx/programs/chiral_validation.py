"""Counts of various chiral volume outlier classes"""
from __future__ import absolute_import, division, print_function

import os
from mmtbx.validation.restraints import chiralities
from mmtbx.model import manager
from libtbx.program_template import ProgramTemplate
from libtbx.utils import Sorry
from libtbx.utils import null_out
from datetime import datetime

class Program(ProgramTemplate):
  prog = os.getenv('LIBTBX_DISPATCHER_NAME')
  description="""
%(prog)s file.pdb [params.eff] [options ...]

Options:

  model=input_file      input PDB or mmCIF file
  outliers_only=True    Only return chiral outliers
  kinemage=False        Create kinemage markup (overrides text output)
  json=False            Outputs results as JSON compatible dictionary
  help=False          Prints this help message if true

  counts of various chiral volume outlier classes, including the following:
    tetrahedral geometry outliers (e.g. flattened geometry)
    chiral identity swaps (e.g. L vs D amino acids)
    pseudochiral naming issues (e.g. swapped chemically identical atoms with distinct names)

Example:

  %(prog)s model=1ubq.pdb kinemage=True
""" % locals()

  master_phil_str = """
    outliers_only = True
      .type = bool
      .help = "Only show outliers"
    kinemage = False
      .type = bool
      .help = "Prints kinemage markup for chiral volume outliers"
    json = False
      .type = bool
      .help = "Prints results as JSON format dictionary"
    result_file = None
      .type = path
      .help = "Path for output file"
"""
  datatypes = ['model','phil']
  data_manager_options = ['model_skip_expand_with_mtrix']

  def validate(self):
    self.data_manager.has_models(raise_sorry=True)

  def run(self):
    f = None
    if self.params.result_file is not None:
      try:
        f = open(self.params.result_file, 'w')
        self.logger.register('file', f)
      except IOError:
        raise Sorry("The output file could not be opened")
    model = self.data_manager.get_model()
    model.set_stop_for_unknowns(False)
    p = manager.get_default_pdb_interpretation_params()
    ##print(dir(p.pdb_interpretation))
    p.pdb_interpretation.allow_polymer_cross_special_position=True
    p.pdb_interpretation.flip_symmetric_amino_acids=False
    p.pdb_interpretation.clash_guard.nonbonded_distance_threshold = None
    model.set_log(log = null_out())
    model.process(make_restraints=True, pdb_interpretation_params=p)
    geometry_restraints_manager = model.get_restraints_manager().geometry
    pdb_hierarchy = model.get_hierarchy()
    self.info_json = {"model_name":self.data_manager.get_default_model_name(),
                      "time_analyzed": str(datetime.now())}
    pdb_hierarchy.atoms().reset_i_seq()
    xray_structure = model.get_xray_structure()
    from mmtbx import restraints
    restraints_manager = restraints.manager(
      geometry=geometry_restraints_manager)
    sites_cart = xray_structure.sites_cart()
    hd_selection = xray_structure.hd_selection()
    pdb_atoms = pdb_hierarchy.atoms()
    energies_sites = restraints_manager.energies_sites(
      sites_cart=sites_cart,
      compute_gradients=False).geometry
    restraint_proxies = getattr(restraints_manager.geometry, "chirality_proxies")
    self.results = chiralities(
        pdb_atoms=pdb_atoms,
        sites_cart=sites_cart,
        energies_sites=energies_sites,
        restraint_proxies=restraint_proxies,
        unit_cell=xray_structure.unit_cell(),
        ignore_hd=True,
        sigma_cutoff=4.0,
        outliers_only=self.params.outliers_only,
        use_segids_in_place_of_chainids=False)

    if self.params.kinemage:
      print(self.results.as_kinemage(), file=self.logger)
    elif self.params.json:
      print(self.results.as_JSON(self.info_json), file=self.logger)
    else:
      self.results.show(out=self.logger, verbose=True)
    if f:
      try:
        f.close()
      except Exception:
        raise Sorry("Could not close output file")

  def get_results(self):
    return self.results

  def get_results_as_JSON(self):
    return self.results.as_JSON(self.info_json)
