from cctbx.array_family import flex
import mmtbx
from mmtbx.monomer_library.pdb_interpretation import monomer_mapping

from mmtbx.monomer_library.pdb_interpretation import (
geometry_restraints_proxy_registries,
ener_lib_as_nonbonded_params,
nonbonded_energy_type_registry,
master_params,
build_all_chain_proxies,
monomer_mapping,
is_same_model_as_before
)


class NullObject:
  # Returns None for all dot syntax access
  def __getattr__(self, name):
    return None


class AtomView:
  """
  A 'view' of a iotbx.pdb atom. Extreme effort should be made to not
  store data here.
  """
  def __init__(self,parent,atom):
    self._parent = parent # a ModelView object
    self._atom = atom
    self._residue_idx = None
    #assert self.parent.atoms[self.atom.i_seq] == self.atom, "AtomView i_seq in parent ModelView.atoms does not return AtomView's atom"


  @property
  def parent(self):
    return self._parent
  @property
  def i_seq(self):
    return self.atom.i_seq

  @property
  def conformer_idx(self):
    return self.parent.conformer_indices[self.i_seq]

  @property
  def model_h_idx(self):
    return self.parent.model_h_indices[self.i_seq]

  @property
  def model_h(self):
    return self.atom.parent().parent().parent().parent()

  @property
  def chain(self):
    return self.atom.parent().parent().parent()

  @property
  def ag(self):
    return self.atom.parent()

  @property
  def rg(self):
    return self.atom.parent().parent()

  @property
  def atom(self):
    return self._atom

  @property
  def mon_lib_srv(self):
    return self.parent.mon_lib_srv

  @property
  def mm_id(self):
    # monomer map id
    mm_id = (self.conformer_idx,self.residue_idx)
    return mm_id

  @property
  def residue(self):
    for conformer in self.chain.conformers():
      for residue in conformer.residues():
        if self.atom in residue.atoms():
          return residue

  @property
  def residue_idx(self):
    if self._residue_idx is None:
      for conformer in self.chain.conformers():
        for i,residue in enumerate(conformer.residues()):
          if self.atom in residue.atoms():
            self._residue_idx = i
            return self._residue_idx
    return self._residue_idx

  @property
  def prev_residue(self):
    for conformer in self.chain.conformers():
      for i,residue in enumerate(conformer.residues()):
        if self.atom in residue.atoms():
          if i>0:
            return conformer.residues()[i-1]

          else:
            return None
  @property
  def next_residue(self):
    for conformer in self.chain.conformers():
      for i,residue in enumerate(conformer.residues()):
        if self.atom in residue.atoms():
          if i<len(conformer.residues())-1:
            conformer.residues()[i+1]
          else:
            return None

  @property
  def mm(self):
    if self.mm_id in self.parent.monomer_mappings:
      return self.parent.monomer_mappings[self.mm_id]
    mm = monomer_mapping(
      pdb_atoms=self.residue.atoms,
      mon_lib_srv=self.mon_lib_srv,
      translate_cns_dna_rna_residue_names=False,
      rna_sugar_pucker_analysis_params=NullObject(),
      apply_cif_modifications={},
      apply_cif_links_mm_pdbres_dict={},
      i_model=self.model_h_idx,
      i_conformer=self.conformer_idx,
      is_first_conformer_in_chain=True, # TODO: ??
      conf_altloc=None,# TODO: ??
      pdb_residue=self.residue,
      next_pdb_residue=self.next_residue,
      chainid=self.chain.id,
      specific_residue_restraints=None)
    self.parent.monomer_mappings[self.mm_id] = mm
    return mm

  @property
  def charge(self):
    return self.atom.charge_as_int()
  @property
  def element(self):
    return self.atom.element.strip()

  @property
  def nb_energy_type(self):
    return self.parent.nb_energy_type_registry[self.i_seq]

class ModelView:
  """
  A 'view' of an mmtbx.model.manager. Extreme effort should be made to not
  store data here.

  This class also provides access to 'model-centric' data for pdb interpretation
  and geometry restraints.
  """
  def __init__(self,model):
    self._model = model
    self._atom_views = [AtomView(self,atom) for atom in self.atoms]
    self._conformer_indices = None
    self._model_indices = None
    self._mon_lib_srv = mmtbx.monomer_library.server.server() # slow
    self._residues = None
    self._nb_energy_type_registry = None
    self._monomer_mappings = {} # ag keys, mm values

  def __len__(self):
    return self.model.get_number_of_atoms()

  @property
  def model(self):
    return self._model
  @property
  def hierarchy(self):
    return self.model.get_hierarchy()
  @property
  def atoms(self):
    return self.model.get_atoms()
  @property
  def residues(self):
    if self._residues is None:
      self._residues = list(self.hierarchy.residue_groups())
    return self._residues
  @property
  def atom_views(self):
    return self._atom_views
  @property
  def mon_lib_srv(self):
    return self._mon_lib_srv
  @property
  def conformer_indices(self):
    if self._conformer_indices is None:
      altloc_i_conformer = {"":0}
      conformer_indices = flex.size_t(len(self), 0)
      altloc_indices = self.hierarchy.altloc_indices()
      if ("" in altloc_indices): p = 0
      else:                      p = 1
      altlocs = sorted(altloc_indices.keys())
      for i,altloc in enumerate(altlocs):
        if (altloc == ""): continue
        conformer_indices.set_selected(altloc_indices[altloc], i+p)
        altloc_i_conformer[altloc] = i+p
      self._conformer_indices = conformer_indices#.as_numpy_array()
    return self._conformer_indices

  @property
  def model_h_indices(self):
    if self._model_indices is None:
      n_seq = len(self)
      model_indices = flex.size_t(n_seq, n_seq)
      for i_model,m in enumerate(self.hierarchy.models()):
        model_indices.set_selected(m.atoms().extract_i_seq(), i_model)
      assert model_indices.count(n_seq) == 0
      self._model_indices = model_indices#.as_numpy_array()
    return self._model_indices

  @property
  def nb_energy_type_registry(self,strict_conflict_handling=False):
    if self._nb_energy_type_registry is None:
      nb_energy_type_registry = nonbonded_energy_type_registry(
        n_seq=self.model.get_atoms().size(),
        strict_conflict_handling=strict_conflict_handling)
      for atomv in self.atom_views:
        mm = atomv.mm
        nb_energy_type_registry.assign_from_monomer_mapping(mm.conf_altloc,mm)
        nb_energy_type_registry.assign_charge(atomv.i_seq,atomv.charge)
      self._nb_energy_type_registry = nb_energy_type_registry
    return self._nb_energy_type_registry

  @property
  def monomer_mappings(self):
    return self._monomer_mappings

  @property
  def all_monomer_mappings(self):
    for atomv in self.atom_views:
      mm = atomv.mm
    return list(self._monomer_mappings.values())

  @property
  def all_monomer_mapping_summaries(self):
    return [mm.summary() for mm in self.all_monomer_mappings]