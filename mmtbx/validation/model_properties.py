
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
  This does not actually contain individual outliers, which are instead held
  in the xray_structure_statistics objects for subsets of the model.
  """
  __slots__ = [
    "all",
    "_macromolecules",
    "_ligands",
    "_water",
    "n_models",
    "n_atoms",
    "n_hydrogens",
    "n_waters",
    "n_polymer",
    "n_protein",
    "n_nuc",
    "ignore_hd",
  ]
  def __init__ (self, pdb_hierarchy,
      xray_structure,
      all_chain_proxies=None,
      ignore_hd=True) :
    for name in self.__slots__ :
      setattr(self, name, None)
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
    self.all = xray_structure_statistics(
      pdb_hierarchy=pdb_hierarchy,
      xray_structure=xray_structure,
      ignore_hd=ignore_hd)
    if (all_chain_proxies is not None) :
      macro_sel = all_chain_proxies.selection("protein or rna or dna")
      water_sel = all_chain_proxies.selection("water")
      ligand_sel = all_chain_proxies.selection(
        "not (protein or dna or rna or water)")
      if (macro_sel.count(True) > 0) :
        self._macromolecules = xray_structure_statistics(
          pdb_hierarchy=pdb_hierarchy.select(macro_sel),
          xray_structure=xray_structure.select(macro_sel),
          ignore_hd=ignore_hd)
      if (water_sel.count(True) > 0) :
        self._water = xray_structure_statistics(
          pdb_hierarchy=pdb_hierarchy.select(water_sel),
          xray_structure=xray_structure.select(water_sel),
          ignore_hd=ignore_hd)
      if (ligand_sel.count(True) > 0) :
        self._ligands = xray_structure_statistics(
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

  @property
  def macromolecules (self) :
    if (self._macromolecules is None) :
      return dummy_validation()
    return self._macromolecules

  @property
  def water (self) :
    if (self._water is None) :
      return dummy_validation()
    return self._water

  @property
  def ligands (self) :
    if (self._ligands is None) :
      return dummy_validation()
    return self._ligands

class residue_occupancy (residue) :
  __slots__ = residue.__slots__ + ["chain_type", "b_iso"]

  def as_string (self, prefix="") :
    return " %3s%2s%5s (all)  occ=%.2f" % (self.resname, self.chain_id,
      self.resid, self.occupancy)

  def as_table_row_phenix (self) :
    return [ self.id_str(), "residue", self.occupancy, self.b_iso ]

class atom_occupancy (atom) :
  """
  Container for single-atom occupancy outliers (usually atoms with zero
  occupancy).
  """
  def as_string (self, prefix="") :
    return "%s %4s   occ=%.2f" % (self.atom_group_id_str(), self.name,
      self.occupancy)

  def as_table_row_phenix (self) :
    return [ self.id_str(), "atom", self.occupancy, self.b_iso ]

class residue_bfactor (residue_occupancy) :

  def as_string (self, prefix="") :
    return "%s %4s   b_iso=%.2f" % (self.atom_group_id_str(), self.name,
      self.b_iso)

class atom_bfactor (atom_occupancy) :
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
    "n_zero_b",
    "n_zero_occ",
    "n_npd",
    "b_mean",
    "b_min",
    "b_max",
    "o_mean",
    "o_min",
    "o_max",
    "zero_occ",
    "partial_occ",
    "bad_adps",
    "b_histogram",
  ]
  gui_list_headers = ["Atom(s)", "Type", "Occupancy", "Isotropic B-factor"]
  gui_formats = ["%s","%s","%.2f", "%.2f"]
  wx_column_widths = [300,100,100,200]
  def __init__ (self, pdb_hierarchy, xray_structure, ignore_hd=True,
      collect_outliers=True) :
    for name in self.__slots__ :
      setattr(self, name, None)
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
    self.n_zero_b = (u_isos == 0).count(True)
    self.n_zero_occ = (occ == 0).count(True)
    mv = flex.mean_and_variance(u_isos.select(u_isos > 0))
    sigma = mv.unweighted_sample_standard_deviation()
    u_cutoff_high = mv.mean() + (4.0 * sigma)
    u_cutoff_low = mv.mean() - (4.0 * sigma)
    self.b_mean = adptbx.u_as_b(flex.mean(u_isos))
    self.b_min = adptbx.u_as_b(flex.min(u_isos))
    self.b_max = adptbx.u_as_b(flex.max(u_isos))
    self.o_mean = flex.mean(occ)
    self.o_min = flex.min(occ)
    self.o_max = flex.max(occ)
    self.n_outliers = self.n_aniso_h + self.n_npd
    self.zero_occ = []
    self.partial_occ = []
    self.bad_adps = [] # TODO
    self.b_histogram = None # TODO
    def is_u_iso_outlier (u) :
      return (u < u_cutoff_low) or (u > u_cutoff_high) or (u <= 0)
    # these statistics cover all atoms!
    occupancies = xray_structure.scatterers().extract_occupancies()
    u_isos = xray_structure.extract_u_iso_or_u_equiv()
    collected = flex.bool(occupancies.size(), False)
    if (collect_outliers) :
      for i_seq, occ in enumerate(occupancies) :
        if (hd_selection[i_seq] and ignore_hd) or collected[i_seq] :
          continue
        pdb_atom = pdb_atoms[i_seq]
        parent = pdb_atom.parent()
        if (occ <= 0) :
          group_atoms = parent.atoms()
          labels = pdb_atom.fetch_labels()
          if (len(group_atoms) > 1) and (group_atoms.extract_occ().all_eq(0)) :
            i_seqs = group_atoms.extract_i_seq()
            b_mean = adptbx.u_as_b(flex.mean(u_isos.select(i_seqs)))
            outlier = residue_occupancy(
              chain_id=labels.chain_id,
              resseq=labels.resseq,
              icode=labels.icode,
              altloc=labels.altloc,
              resname=labels.resname,
              occupancy=occ,
              outlier=True,
              xyz=group_atoms.extract_xyz().mean(),
              b_iso=b_mean)
            self.zero_occ.append(outlier)
            self.n_outliers += 1
            collected.set_selected(i_seqs, True)
          else :
            assert (pdb_atom.occ == occ), "%s: %s <--> %s" % (pdb_atom.id_str(),
              pdb_atom.occ, occ)
            outlier = atom_occupancy(
              pdb_atom=pdb_atom,
              occupancy=occ,
              b_iso=adptbx.u_as_b(u_isos[i_seq]),
              xyz=pdb_atom.xyz,
              outlier=True)
            self.zero_occ.append(outlier)
            self.n_outliers += 1
        elif is_u_iso_outlier(u_isos[i_seq]) :
          # zero displacements will always be recorded on a per-atom basis
          if (u_isos[i_seq] <= 0) :
            outlier = atom_bfactor(
              pdb_atom=pdb_atom,
              occupancy=occ,
              b_iso=adptbx.u_as_b(u_isos[i_seq]),
              xyz=pdb_atom.xyz,
              outlier=True)
            self.bad_adps.append(outlier)
            self.n_outliers += 1
          else :
            # if the average displacement for the entire residue falls outside
            # the cutoffs, save as a single residue outlier
            group_atoms = parent.atoms()
            i_seqs = group_atoms.extract_i_seq()
            u_mean = flex.mean(u_isos.select(i_seqs))
            if is_u_iso_outlier(u_mean) :
              labels = pdb_atom.fetch_labels()
              outlier = residue_bfactor(
                chain_id=labels.chain_id,
                resseq=labels.resseq,
                icode=labels.icode,
                altloc=labels.altloc,
                resname=labels.resname,
                occupancy=occ,
                outlier=True,
                xyz=group_atoms.extract_xyz().mean(),
                b_iso=adptbx.u_as_b(u_mean))
              self.bad_adps.append(outlier)
              self.n_outliers += 1
              collected.set_selected(i_seqs, True)
            # otherwise, just save this atom
            else :
              outlier = atom_bfactor(
                pdb_atom=pdb_atom,
                occupancy=occ,
                b_iso=adptbx.u_as_b(u_isos[i_seq]),
                xyz=pdb_atom.xyz,
                outlier=True)
              self.bad_adps.append(outlier)
              self.n_outliers += 1
      # FIXME not totally confident that this will always work...
      first_model = pdb_hierarchy.models()[0]
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
                occupancy=total_occ,
                chain_type=chain_type,
                outlier=True,
                xyz=xyz)
              self.partial_occ.append(outlier)
              self.n_outliers += 1

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

  def iter_results (self, property_type=None, outliers_only=True) :
    if (property_type is None) :
      outliers = self.zero_occ + self.partial_occ + self.bad_adps
    elif (property_type == "occupancy") :
      outliers = self.zero_occ + self.partial_occ
    elif (property_type == "b_factor") :
      outliers = self.bad_adps
    else :
      raise RuntimeError("Unknown property type '%s'" % property_type)
    for result in outliers :
      if (result.is_outlier()) or (not outliers_only) :
        yield result

  def as_gui_table_data (self, property_type=None, outliers_only=True,
      include_zoom=False) :
    table = []
    for result in self.iter_results(property_type=property_type,
        outliers_only=outliers_only) :
      extra = []
      if (include_zoom) :
        extra = result.zoom_info()
      row = result.as_table_row_phenix()
      assert (len(row) == len(self.gui_list_headers) == len(self.gui_formats))
      table.append(row + extra)
    return table
