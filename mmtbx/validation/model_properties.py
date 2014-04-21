
# TODO unit testing

"""
Analysis of model properties, independent of data.
"""

from __future__ import division
from mmtbx.validation import atom, residue, validation, dummy_validation
from libtbx import slots_getstate_setstate
import sys

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
    "n_hydrogens",
    "n_waters",
    "n_polymer",
    "n_protein",
    "n_nuc",
    "occupancy_outliers",
    "ignore_hd",
  ]
  def __init__ (self, pdb_hierarchy,
      xray_structure,
      all_chain_proxies=None,
      ignore_hd=True) :
    self.ignore_hd = ignore_hd
    self.n_atoms = xray_structure.scatterers().size()
    self.n_hydrogens = xray_structure.hd_selection().count(True)
    self.n_models = len(pdb_hierarchy.models())
    first_model = pdb_hierarchy.models()[0]
    counts = pdb_hierarchy.overall_counts()
    resnames = counts.resnames
    self.n_waters = resnames.get("HOH", 0) + resnames.get("WAT", 0)
    self.n_polymer = 0
    self.n_nuc = self.n_protein = 0
    self.occupancy_outliers = []
    # FIXME not totally confident that this will always work...
    for chain in first_model.chains() :
      is_protein = chain.is_protein()
      is_na = False
      chain_type = "other"
      if (not is_protein) :
        is_na = chain.is_na()
        chain_type = "nucleic acid"
      else :
        chain_type = "protein"
      residue_groups = chain.residue_groups()
      if (is_protein) or (is_na) :
        n_res = len(residue_groups)
        self.n_polymer += n_res
        if (is_protein) :
          self.n_protein += n_res
        else :
          self.n_nuc += n_res
      for residue_group in chain.residue_groups() :
        atom_groups = residue_group.atom_groups()
        if ((len(atom_groups) > 1) and (len(atom_groups[0].atoms()) != 1) and
            (len(set([ ag.resname for ag in atom_groups ])) == 1)) :
          total_occ = sum([ ag.atoms()[0].occ for ag in atom_groups ])
          if (total_occ != 1) and (total_occ != 0) :
            xyz = residue_group.atoms().extract_xyz().mean()
            outlier = residue_occupancy(
              chain_id=chain.id,
              resseq=residue_group.resseq,
              icode=residue_group.icode,
              resname=atom_groups[0].resname,
              occupancy=occupancy,
              chain_type=chain_type,
              outlier=True,
              xyz=xyz)
            self.occupancy_outliers.append(outlier)
    self.occupancy_outliers.sort(lambda a,b: cmp(b.occupancy, a.occupancy))
    self.all = xray_structure_statistics(
      pdb_hierarchy=pdb_hierarchy,
      xray_structure=xray_structure,
      ignore_hd=ignore_hd)
    self.macromolecules = dummy_validation()
    self.ligands = dummy_validation()
    self.water = dummy_validation()
    if (all_chain_proxies is not None) :
      macro_sel = all_chain_proxies.selection("protein or rna or dna")
      water_sel = all_chain_proxies.selection("water")
      ligand_sel = all_chain_proxies.selection(
        "not (protein or dna or rna or water)")
      if (macro_sel.count(True) > 0) :
        self.macromolecules = xray_structure_statistics(
          pdb_hierarchy=pdb_hierarchy.select(macro_sel),
          xray_structure=xray_structure.select(macro_sel),
          ignore_hd=ignore_hd)
      if (water_sel.count(True) > 0) :
        self.water = xray_structure_statistics(
          pdb_hierarchy=pdb_hierarchy.select(water_sel),
          xray_structure=xray_structure.select(water_sel),
          ignore_hd=ignore_hd)
      if (ligand_sel.count(True) > 0) :
        self.ligands = xray_structure_statistics(
          pdb_hierarchy=pdb_hierarchy.select(ligand_sel),
          xray_structure=xray_structure.select(ligand_sel),
          ignore_hd=ignore_hd)

  def show (self, out=sys.stdout, prefix="") :
    print >> out, prefix+"Overall:"
    self.all.show_summary(out=out, prefix=prefix+"  ")
    self.all.show_bad_occupancy(out=out, prefix=prefix+"  ")
    if (self.ligands) or (self.water) :
      for label, props in zip(["Macromolecules", "Ligands", "Waters"],
          [self.macromolecules, self.ligands, self.water]) :
        if (props) :
          print >> out, prefix+"%s:" % label
          props.show_summary(out=out, prefix=prefix+"  ")
    if (self.ignore_hd) and (self.n_hydrogens > 0) :
      print >> out, prefix+"(Hydrogen atoms not included in overall counts.)"

class residue_occupancy (residue) :
  __slots__ = residue.__slots__ + ["chain_type"]

  def as_string (self, prefix="") :
    return " %3s%2s%5s (all)  occ=%.2f" % (self.resname, self.chain_id,
      self.resid, self.occupancy)

class occupancy (atom) :
  """
  Container for single-atom occupancy outliers (usually atoms with zero
  occupancy).
  """
  def as_string (self, prefix="") :
    return "%s %4s   occ=%.2f" % (self.atom_group_id_str(), self.name,
      self.occupancy)

class atom_bfactor (atom) : # TODO
  pass

# FIXME redundant with model_vs_data, but what can I do???
# TODO add validation stuff
class xray_structure_statistics (validation) :
  """
  Occupancy and B-factor statistics.
  """
  __slots__ = validation.__slots__ + [
    "n_all",
    "n_atoms",
    "n_non_hd",
    "n_hd",
    "n_aniso",
    "n_aniso_h",
    "n_zero_occ",
    "n_npd",
    "b_mean",
    "b_min",
    "b_max",
    "o_mean",
    "o_min",
    "o_max",
    "zero_occ",
    "bad_adps",
    "b_histogram",
  ]
  def __init__ (self, pdb_hierarchy, xray_structure, ignore_hd=True,
      collect_outliers=True) :
    validation.__init__(self)
    assert len(xray_structure.scatterers()) != 0
    from cctbx import adptbx
    from scitbx.array_family import flex
    xrs = xray_structure
    self.n_total = xrs.scatterers().size() # always include H/D
    self.results = None
    pdb_atoms = pdb_hierarchy.atoms()
    hd_selection = xrs.hd_selection()
    subtract_hd = True
    if (ignore_hd) and (hd_selection.count(True) > 0) :
      xrs = xrs.select(~hd_selection)
      subtract_hd = False
    u_isos = xrs.extract_u_iso_or_u_equiv()
    occ = xrs.scatterers().extract_occupancies()
    self.n_all = hd_selection.size()
    self.n_atoms = xrs.scatterers().size()
    self.n_hd = hd_selection.count(True)
    self.n_non_hd = self.n_all - self.n_hd
    self.n_aniso = xrs.use_u_aniso().count(True)
    self.n_aniso_h = (xray_structure.use_u_aniso() & hd_selection).count(True)
    self.n_npd = xrs.is_positive_definite_u().count(False)
    self.n_zero_occ = (occ == 0).count(True)
    self.b_mean = adptbx.u_as_b(flex.mean(u_isos))
    self.b_min = adptbx.u_as_b(flex.min(u_isos))
    self.b_max = adptbx.u_as_b(flex.max(u_isos))
    self.o_mean = flex.mean(occ)
    self.o_min = flex.min(occ)
    self.o_max = flex.max(occ)
    self.n_outliers = self.n_aniso_h + self.n_zero_occ + self.n_npd
    self.zero_occ = []
    self.bad_adps = [] # TODO
    self.b_histogram = None # TODO
    # these statistics cover all atoms!
    occupancies = xray_structure.scatterers().extract_occupancies()
    u_isos = xray_structure.extract_u_iso_or_u_equiv()
    collected = flex.bool(occupancies.size(), False)
    if (collect_outliers) :
      for i_seq, occ in enumerate(occupancies) :
        if (hd_selection[i_seq] and ignore_hd) or collected[i_seq] :
          continue
        if (occ <= 0) :
          atom = pdb_atoms[i_seq]
          parent = atom.parent()
          group_atoms = parent.atoms()
          labels = atom.fetch_labels()
          if (len(group_atoms) > 1) and (group_atoms.extract_occ().all_eq(0)) :
            outlier = residue_occupancy(
              chain_id=labels.chain_id,
              resseq=labels.resseq,
              icode=labels.icode,
              altloc=labels.altloc,
              resname=labels.resname,
              occupancy=occ,
              outlier=True,
              xyz=group_atoms.extract_xyz().mean())
            self.zero_occ.append(outlier)
            collected.set_selected(group_atoms.extract_i_seq(), True)
          else :
            assert (atom.occ == occ), "%s: %s <--> %s" % (atom.id_str(),
              atom.occ, occ)
            outlier = occupancy(
              pdb_atom=atom,
              occupancy=occ,
              b_iso=adptbx.u_as_b(u_isos[i_seq]),
              outlier=True)
            self.zero_occ.append(outlier)

  def show_summary (self, out=sys.stdout, prefix="") :
    print >> out, prefix + "Number of atoms = %d  (anisotropic = %d)" % \
      (self.n_atoms, self.n_aniso)
    print >> out, prefix + "B_iso: mean = %5.1f  max = %5.1f  min = %5.1f" % \
      (self.b_mean, self.b_max, self.b_min)
    if (self.n_aniso_h > 0) :
      print >> out, prefix + "  warning: %d anisotropic hydrogen atoms" % \
        self.n_aniso_h
    if (self.o_min != 1.0) :
      print >> out, prefix + \
        "Occupancy: mean = %4.2f  max = %4.2f  min = %4.2f" % \
        (self.o_mean, self.o_max, self.o_min)
      if (self.n_zero_occ > 0) :
        print >> out, prefix + "  warning: %d atoms with zero occupancy" % \
          self.n_zero_occ
    if (self.n_outliers > 0) :
      print >> out, prefix + \
        "%d total B-factor or occupancy problems detected" % \
        self.n_outliers

  def show_bad_occupancy (self, out=sys.stdout, prefix="") :
    if (len(self.zero_occ) > 0) :
      print >> out, prefix + "Atoms or residues with zero occupancy:"
      for outlier in self.zero_occ :
        print >> out, prefix + str(outlier)

  def show_bfactors (self, out=sys.stdout, prefix="") :
    print >> out, prefix + "B_iso: mean = %5.1f  max = %5.1f  min = %5.1f" % \
      (self.b_mean, self.b_max, self.b_min)
