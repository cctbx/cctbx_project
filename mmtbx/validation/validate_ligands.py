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
    ligand_isel_dict = self.get_ligands(ph = ph)
    for id_tuple, isel in ligand_isel_dict.items():
      for rg in ph.select(isel).residue_groups():
        ligand_dict = {}
        for conformer in rg.conformers():
          altloc = conformer.altloc
          conformer_isel = conformer.atoms().extract_i_seq()
          lr = ligand_result(
            model = self.model,
            isel = conformer_isel)
          lr.get_occupancies()
          ligand_dict[altloc] = lr
      self[id_tuple] = ligand_dict

  # ---------------------------------------------------------------------------

  @staticmethod
  def get_ligands(ph):
    # Store ligands as list of iselections --> better way? Careful if H will be
    # added at some point!
    ligand_isel_dict = {}
    get_class = iotbx.pdb.common_residue_names_get_class
    exclude = ["common_amino_acid", "modified_amino_acid", "common_rna_dna",
               "modified_rna_dna", "ccp4_mon_lib_rna_dna", "common_water",
                "common_element"]
    for model in ph.models():
      for chain in model.chains():
        for rg in chain.residue_groups():
          for resname in rg.unique_resnames():
            if (not get_class(name=resname) in exclude):
              iselection = rg.atoms().extract_i_seq()
              id_tuple = (model.id, chain.id, rg.resseq)
              ligand_isel_dict[id_tuple] = iselection
    return ligand_isel_dict

  # ---------------------------------------------------------------------------

  def print_ligand_counts(self):
    print('\nThe following ligands were found in the input model:')
    for id_tuple, ligand_dict in self.items():
      for altloc, lr in ligand_dict.items():
        print(lr.resname, lr.id_str, altloc)


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
