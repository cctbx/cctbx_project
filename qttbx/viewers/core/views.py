

from cctbx.array_family import flex
import mmtbx
from mmtbx.monomer_library.pdb_interpretation import monomer_mapping
from mmtbx.monomer_library.pdb_interpretation import (
nonbonded_energy_type_registry,
monomer_mapping
)


class NullObject:
  # Returns None for all dot syntax access
  def __getattr__(self, name):
    return None


class AtomView:
  """
  A 'view' of a iotbx.pdb atom. Extreme effort should be made to not store data
  here.
  """
  def __init__(self,parent,conformer_idx,residue_idx,atom):
    self._parent = parent # a ModelView object
    self._conformer_idx = conformer_idx
    self._residue_idx = residue_idx # idx in conformer
    self._atom = atom
    self._mm = None # monomer mapping
    #assert self.parent.atoms[self.atom.i_seq] == self.atom, "AtomView i_seq in parent ModelView.atoms does not return AtomView's atom"

    # validate
   #assert self.conformer_idx == self.parent.conformer_indices[self.i_seq], f"Mismatch conformer indices: {self.conformer_idx} vs {self.parent.conformer_indices[self.i_seq]}"

  @property
  def parent(self): # a ModelView object
    return self._parent
  @property
  def i_seq(self):
    return self.atom.i_seq

  @property
  def chain(self):
    return self.rg.parent()

  @property
  def conformer_idx(self):
    return self._conformer_idx

  @property
  def conformer_index(self):
    return self.parent.conformer_indices[self.i_seq]

  @property
  def conformer(self):
    return self.chain.conformers()[self.conformer_idx]

  @property
  def model(self):
    return self.atom.parent().parent().parent().parent()


  @property
  def ag(self):
    return self.atom.parent()

  @property
  def rg(self):
    return self.ag.parent()

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
    return self.conformer.residues()[self.residue_idx]

  @property
  def residue_previous(self):
    if self.residue.link_to_previous:
      return self.conformer.residues()[self.residue_idx-1]

  @property
  def residue_next(self):
    next_i = self.residue_idx+1
    if next_i<len(self.conformer.residues()):
      residue_next = self.conformer.residues()[next_i]
      if residue_next.link_to_previous:
        return residue_next

  @property
  def residue_idx(self):
    return self._residue_idx

  # @property
  # def prev_residue(self):
  #   for conformer in self.chain.conformers():
  #     for i,residue in enumerate(conformer.residues()):
  #       if self.atom in residue.atoms():
  #         if i>0:
  #           return conformer.residues()[i-1]

  #         else:
  #           return None
  # @property
  # def next_residue(self):
  #   for conformer in self.chain.conformers():
  #     for i,residue in enumerate(conformer.residues()):
  #       if self.atom in residue.atoms():
  #         if i<len(conformer.residues())-1:
  #           conformer.residues()[i+1]
  #         else:
  #           return None

  @property
  def mm(self):
   return self._mm

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
  A 'view' of an mmtbx.model.manager. Extreme effort should be made to not store
  data here.

  This class also provides access to 'model-centric' data for pdb interpretation
  and geometry restraints.
  """
  def __init__(self,model,log=None):
    self._model = model
    self._atom_views = []
    self._conformer_indices = None
    self._model_indices = None
    self._mon_lib_srv = mmtbx.monomer_library.server.server() # slow
    self._residues = None
    self._nb_energy_type_registry = None
    self._monomer_mappings = {} # ag keys, mm values
    self._altloc_i_conformer = None
    self._mm_ids = None

    models = model.get_hierarchy().models()
    model_set = set()
    n_unique_models = 0
    model_type_indices = [-1] * len(models)
    self._all_monomer_mappings = []

    for i_model,model in enumerate(models):
      # NOTE: printing to null_out takes time. Faster to just check if log is None
      if (log is not None):
        print('  Model: "%s"' % model.id, file=log)
      is_unique_model = model not in model_set
      if (is_unique_model):
        n_unique_models += 1
        model_set.add(model)
      elif (log is not None):
        print("    Same as model", \
          models[model_type_indices[i_model]].id, file=log)
      if (is_unique_model and log is not None):
        print("    Number of chains:", model.chains_size(), file=log)
      for chain in model.chains():
        conformers = chain.conformers()
        if (is_unique_model and log is not None):
          print('    Chain: "%s"' % chain.id, file=log)
          print("      Number of atoms:", chain.atoms_size(), file=log)
          print("      Number of conformers:", len(conformers), file=log)
        for j_conformer,conformer in enumerate(conformers):
          if (is_unique_model and log is not None):
            print('      Conformer: "%s"' % conformer.altloc, file=log)
          self.altloc_i_conformer[conformer.altloc]
          # residues
          residues = conformer.residues()
          print("        Number of residues:", len(residues), file=log)
          for i_residue,residue in enumerate(residues):
            mm = None
            for i_atom,atom in enumerate(residue.atoms()):
              atomv = AtomView(self,j_conformer,i_residue,atom)

              # set up mm on the first atom view
              if mm is None:
                mm = monomer_mapping(
                  pdb_atoms=residue.atoms,
                  mon_lib_srv=self.mon_lib_srv,
                  translate_cns_dna_rna_residue_names=False,
                  rna_sugar_pucker_analysis_params=NullObject(),
                  apply_cif_modifications={},
                  apply_cif_links_mm_pdbres_dict={},
                  i_model=i_model,
                  i_conformer=j_conformer,
                  is_first_conformer_in_chain=True, # TODO: ??
                  conf_altloc=None,# TODO: ??
                  pdb_residue=residue,
                  next_pdb_residue=atomv.residue_next,
                  chainid=atomv.chain.id,
                  specific_residue_restraints=None)
                self._all_monomer_mappings.append(mm)
              atomv._mm = mm
              self._atom_views.append(atomv)

    print("Total monomer mappings:", len(self._all_monomer_mappings), file=log)
    print("Total atom views:",len(self._atom_views),file=log)

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
  # @property
  # def residues(self):
  #   if self._residues is None:
  #     self._residues = list(self.hierarchy.residue_groups())
  #   return self._residues

  @property
  def mm_ids(self):
    # monomer mapping ids (conformer_idx,residue_idx)
    if self._mm_ids is None:
      self._mm_ids = [atomv.mm_id for atomv in self.atom_views]
    return self._mm_ids
  @property
  def mm_id_set(self):
    return set(self.mm_ids)

  @property
  def atom_views(self):
    return self._atom_views
  @property
  def mon_lib_srv(self):
    return self._mon_lib_srv
  @property
  def conformer_indices(self):
    if self._conformer_indices is None:
      self._altloc_i_conformer = {"":0}
      conformer_indices = flex.size_t(len(self), 0)
      altloc_indices = self.hierarchy.altloc_indices()
      if ("" in altloc_indices): p = 0
      else:                      p = 1
      altlocs = sorted(altloc_indices.keys())
      for i,altloc in enumerate(altlocs):
        if (altloc == ""): continue
        conformer_indices.set_selected(altloc_indices[altloc], i+p)
        self._altloc_i_conformer[altloc] = i+p
      self._conformer_indices = conformer_indices
    return self._conformer_indices

  @property
  def altloc_i_conformer(self):
    if self._altloc_i_conformer is None:
      _ = self.conformer_indices # set it
    return self._altloc_i_conformer

  @property
  def model_indices(self):
    if self._model_indices is None:
      n_seq = len(self)
      model_indices = flex.size_t(n_seq, n_seq)
      for i_model,m in enumerate(self.hierarchy.models()):
        model_indices.set_selected(m.atoms().extract_i_seq(), i_model)
      assert model_indices.count(n_seq) == 0
      self._model_indices = model_indices
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
  def all_monomer_mappings(self):
    return self._all_monomer_mappings