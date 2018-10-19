from __future__ import division, print_function

import iotbx.pdb
from libtbx import group_args

# =============================================================================

class validate_ligands(object):

  def __init__(self, model):
    self.model = model

    self.ligand_isels = None

  # ---------------------------------------------------------------------------

  def validate_inputs(self):
    if self.model is None:
      raise Sorry("No input model.")
      #return 0

  # ---------------------------------------------------------------------------

  def get_ligands(self, ph):
    # Store ligands as list of iselections --> better way?
    ligand_isels = []
    get_class = iotbx.pdb.common_residue_names_get_class
    exclude = ["common_amino_acid", "modified_amino_acid", "common_rna_dna",
               "modified_rna_dna", "ccp4_mon_lib_rna_dna", "common_water",
                "common_element"]
    # loop through models?
    for chain in ph.chains():
      for rg in chain.residue_groups():
        for resname in rg.unique_resnames():
          if (not get_class(name=resname) in exclude):
            iselection = rg.atoms().extract_i_seq()
            ligand_isels.append(iselection)
    self.ligand_isels = ligand_isels

  # ---------------------------------------------------------------------------

  def run(self):
    ph = self.model.get_hierarchy()

    # Get ligands
    self.get_ligands(ph = ph)

  # ---------------------------------------------------------------------------

  def get_results(self):
    return group_args(
      ligand_isels = self.ligand_isels)
