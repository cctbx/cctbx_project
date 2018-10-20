from __future__ import division, print_function

import iotbx.pdb
#from libtbx import group_args

# =============================================================================
# Manager class for ALL ligands

class manager(dict):

  def __init__(self, model):
    self.model = model

  # ---------------------------------------------------------------------------

  def run(self):
    ph = self.model.get_hierarchy()
    ligand_isels = self.get_ligands(ph = ph)
    for i, isel in zip(range(len(ligand_isels)), ligand_isels):
      lr = ligand_result(
        model = self.model,
        isel  = isel)
      lr.get_occupancies()
      self[i] = lr

  # ---------------------------------------------------------------------------

  @staticmethod
  def get_ligands(ph):
    # Store ligands as list of iselections --> better way? Careful if H will be
    # added at some point!
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

# =============================================================================
# Class storing info per ligand

class ligand_result(object):

  def __init__(self, model, isel):
    self.model = model
    self.isel = isel
    # results
    self._occupancies = None


    self.ph = self.model.get_hierarchy()

    rg_ligand = self.ph.select(self.isel).only_residue_group()
    self.resname = ",".join(rg_ligand.unique_resnames())
    self.id_str = rg_ligand.id_str()

  # ---------------------------------------------------------------------------

  def get_occupancies(self):
    if self._occupancies is None:
      self._occupancies = self.ph.select(self.isel).occupancy_counts()
    return self._occupancies
