from .utils import make_reskey

class ResidueValidationData:
  """A simple class to store validation data for a single residue."""
  def __init__(self, rg, model_id, chain_id, resseq, icode):
    self.reskey = make_reskey(model_id, chain_id, resseq, icode)
    self.chain_id = chain_id
    self.resseq = resseq
    self.icode = icode
    self.resnames = self._get_resnames_from_hierarchy(rg)
    self.resname = self.resnames[sorted(self.resnames.keys())[0]]
    self.validations = {}
    self.has_outlier = False
    self.alternates = []
    self.modeled_alternates = self.get_alts_from_hierarchy(rg)

  def get_alternates(self):
    """
    Combines alternates from the model and from validation results.
    """
    # Start with a set of alternates found in the model structure
    all_alts = set(self.modeled_alternates)
    # Add alternates found in any of the validation results
    for validation_dict in self.validations.values():
        all_alts.update(validation_dict.keys())
    # Ensure there's at least one entry ('') for non-alternate residues
    if not all_alts:
        all_alts.add('')
    # Sort the final, unique list of alternates
    self.alternates = sorted(list(all_alts))

  def _get_resnames_from_hierarchy(self, rg):
    resnames = {}
    for ag in rg.atom_groups():
      resnames[ag.altloc] = ag.resname
    return resnames

  def get_alts_from_hierarchy(self, rg):
    #Help track whether an alt came from the model or from calculations
    modeled_alternates = []
    for ag in rg.atom_groups():
      modeled_alternates.append(ag.altloc.strip())
    return modeled_alternates

  def find_resname(self, alt=''):
    return self.resnames.get(alt, self.resname)

  def get_row_data(self, table_order):
    """
    Returns a dictionary of raw, unformatted data for this residue,
    ready to be passed to the HtmlBuilder.
    """
    data = {
      'chain_id': self.chain_id,
      'resseq': self.resseq,
      'icode': self.icode,
      'resnames': self.resnames,
      'alternates': self.alternates,
      'modeled_alternates': self.modeled_alternates,
      'reskey': self.reskey,
    }

    # Add the raw validation results for the requested columns
    for val_type in table_order:
      if val_type in self.validations:
        data[val_type] = self.validations[val_type]

    return data
