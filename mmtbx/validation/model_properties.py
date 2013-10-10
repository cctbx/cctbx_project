
"""
Analysis of model properties, independent of data.
"""

from __future__ import division
from mmtbx.validation import atom
from libtbx import slots_getstate_setstate

class model_statistics (slots_getstate_setstate) :
  """
  Atom statistics for the overall model, and various selections within.
  """
  __slots__ = [
    "all",
    "macromolecules",
    "ligands",
    "water",
    "n_models",
    "n_atoms",
    "n_waters",
    "n_polymer",
  ]
  def __init__ (self, pdb_hierarchy, xray_structure, ignore_hd=True) :
    self.n_atoms = xray_structure.scatterers().size()
    self.n_models = len(pdb_hierarchy.models())
    first_model = pdb_hierarchy.models()[0]
    counts = pdb_hierarchy.overall_counts()
    resnames = counts.resnames
    self.n_waters = resnames.get("HOH", 0) + resnames.get("WAT", 0)
    self.n_polymer = 0
    # FIXME not totally confident that this will always work...
    for chain in first_model.chains() :
      if (chain.is_protein()) or (chain.is_na()) :
        self.n_polymer += len(chain.residue_groups())
    self.all = xray_structure_statistics(
      pdb_hierarchy=pdb_hierarchy,
      xray_structure=xray_structure,
      ignore_hd=ignore_hd)
    self.macromolecules = None
    self.ligands = None
    self.water = None

class occupancy (atom) : # TODO
  pass

class atom_bfactor (atom) : # TODO
  pass

# FIXME redundant with model_vs_data, but what can I do???
# TODO add validation stuff
class xray_structure_statistics (slots_getstate_setstate) :
  """
  Occupancy and B-factor statistics.
  """
  __slots__ = [
    "n_atoms",
    "n_aniso",
    "n_zero_occ",
    "n_npd",
    "b_mean",
    "b_min",
    "b_max",
    "o_mean",
    "o_min",
    "o_max",
    "aniso_h",
    "zero_occ",
    "bad_adps",
    "b_histogram",
  ]
  def __init__ (self, pdb_hierarchy, xray_structure, ignore_hd=True) :
    assert len(xray_structure.scatterers()) != 0
    from cctbx import adptbx
    from scitbx.array_family import flex
    xrs = xray_structure
    pdb_atoms = pdb_hierarchy.atoms()
    hd_selection = xrs.hd_selection()
    if (ignore_hd) :
      xrs = xrs.select(~hd_selection)
      pdb_atoms = pdb_atoms.select(~hd_selection)
    u_isos = xrs.extract_u_iso_or_u_equiv()
    occ = xrs.scatterers().extract_occupancies()
    self.n_atoms = xrs.scatterers().size()
    self.n_aniso = xrs.use_u_aniso().count(True)
    self.n_npd = xrs.is_positive_definite_u().count(False)
    self.n_zero_occ = (occ == 0).count(True)
    self.b_mean = adptbx.u_as_b(flex.mean(u_isos))
    self.b_min = adptbx.u_as_b(flex.min(u_isos))
    self.b_max = adptbx.u_as_b(flex.max(u_isos))
    self.o_mean = flex.mean(occ)
    self.o_min = flex.min(occ)
    self.o_max = flex.max(occ)
    self.aniso_h = []
    self.bad_adps = []
    self.zero_occ = []
    self.b_histogram = None # TODO
    # these statistics cover all atoms!
    for i_seq, occ in enumerate(occ) :
      if hd_selection[i_seq] :
        continue
      if (occ <= 0) :
        atom = pdb_atoms[i_seq]
        labels = atom.fetch_labels()
        assert (atom.occ == occ), "%s: %s <--> %s" % (atom.id_str(), atom.occ,
          occ)
        outlier = occupancy(
          pdb_atom=atom,
          occupancy=occ,
          b_iso=adptbx.u_as_b(u_iso[i_seq]),
          outlier=True)
        self.zero_occ.append(outlier)
