from __future__ import division, print_function

import iotbx.pdb

# =============================================================================

class validate_ligands(object):

  def __init__(self, model):
    self.model = model

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
    return ligand_isels

  # ---------------------------------------------------------------------------

  def run(self):
    ph = self.model.get_hierarchy()

    # Get ligands
    ligand_isels = self.get_ligands(ph = ph)
    print('\nThe following ligands were found in the input model:')
    for ligand_isel in ligand_isels:
      for rg in ph.select(ligand_isel).residue_groups():
        rn = ",".join(rg.unique_resnames())
        print(rn, rg.id_str())

  # ---------------------------------------------------------------------------

  def get_results(self):
    print('Returning results...')
