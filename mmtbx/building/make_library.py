
# TODO tests, obviously

"""
Tools for assembling an ensemble of related structures for local rebuilding by
homology.
"""

from __future__ import absolute_import, division, print_function
from libtbx import slots_getstate_setstate_default_initializer
from libtbx import Auto, adopt_init_args
from libtbx.utils import null_out
from libtbx import easy_mp
import libtbx.phil
import os
import sys

master_phil = libtbx.phil.parse("""
min_identity = 0.95
  .type = float
xray_only = True
  .type = bool
single_chain = False
  .type = bool
min_resolution = 3.0
  .type = float
sort_by_resolution = True
  .type = bool
ligands = None
  .type = str
""")

def fetch_similar_pdb_ids(
    sequence,
    min_identity=0.95,
    min_resolution=3.0,
    xray_only=True,
    data_only=False,
    expect=0.0000001,
    sort_by_resolution=True,
    ligands=None,
    log=None):
  """
  Finds closely related structures in the RCSB database, optionally limited to
  structures with specific ligands.
  """
  if (log is None) : log = null_out()
  from mmtbx.wwpdb import rcsb_web_services
  allowed_set = None
  if (ligands is not None) and (len(ligands) > 0):
    for lig_str in ligands :
      resnames = lig_str.replace(" ", ",").split(",")
      for resname in resnames :
        ids = rcsb_web_services.chemical_id_search(
          resname=resname,
          d_max=min_resolution,
          xray_only=xray_only,
          data_only=data_only)
        if (allowed_set is None) : allowed_set = set([])
        allowed_set.update([ pdb_id.lower() for pdb_id in ids ])
  matching_ids = rcsb_web_services.sequence_search(
    sequence=sequence,
    min_identity=min_identity,
    expect=expect,
    xray_only=xray_only,
    data_only=data_only,
    d_max=min_resolution,
    sort_by_resolution=sort_by_resolution,
    log=log)
  filtered = []
  for pdb_id in matching_ids :
    pdb_id = pdb_id.lower()
    if (allowed_set is not None) and (not pdb_id in allowed_set):
      continue
    filtered.append(pdb_id)
  return filtered

def load_pdb_models(pdb_ids, fault_tolerant=False, log=null_out()):
  """
  Given a list of PDB IDs, load hierarchies for all of them (ignoring any
  unknown atoms).
  """
  import iotbx.pdb.fetch
  iotbx.pdb.fetch.validate_pdb_ids(pdb_ids)
  ids_and_hierarchies = []
  for pdb_id in pdb_ids :
    try :
      pdb_hierarchy, xray_structure = iotbx.pdb.fetch.load_pdb_structure(
        id=pdb_id,
        allow_unknowns=True,
        local_cache=os.getcwd())
    except Exception as e :
      if (not fault_tolerant):
        raise
      else :
        print("Error - skipping %s:" % pdb_id, file=log)
        print(str(e), file=log)
    else :
      if (len(pdb_hierarchy.models()) > 1):
        print("Warning: %s is multi-MODEL - skipping" % pdb_id, file=log)
      else :
        ids_and_hierarchies.append((pdb_id, pdb_hierarchy))
  return ids_and_hierarchies

def load_all_models_in_directory(dir_name,
    limit_extensions=True,
    recursive=False):
  """
  Load all models in the specified directory, returning a list of file names
  and iotbx.file_reader objects.
  """
  from iotbx.file_reader import any_file, guess_file_type
  assert os.path.isdir(dir_name)
  file_names_and_objects = []
  for file_name in os.listdir(dir_name):
    full_path = os.path.join(dir_name, file_name)
    if os.path.isdir(full_path) and recursive :
      file_names_and_objects.extend(
        load_all_models_in_directory(dir_name=full_path,
          limit_extensions=limit_extensions,
          recursive=True))
    elif os.path.isfile(full_path):
      if (limit_extensions) and (guess_file_type(full_path) != "pdb"):
        continue
      input_file = any_file(full_path,
        raise_sorry_if_not_expected_format=True)
      if (input_file.file_type == "pdb"):
        file_names_and_objects.append((full_path, input_file.file_object))
  return file_names_and_objects

class related_chain(slots_getstate_setstate_default_initializer):
  """
  Container for a protein chain (as hierarchy object) and metadata, extracted
  based on similarity to a target sequence.
  """
  __slots__ = [ "source_info", "chain_id", "pdb_hierarchy", "identity", ]

def find_similar_chains(pdb_hierarchy,
    sequence,
    source_info=None,
    first_chain_only=False,
    remove_alt_confs=True,
    reset_chain_id=None,
    min_identity=0.95,
    min_identity_epsilon=0.02,
    log=null_out()):
  """
  Find all chains with at least the specified fractional sequence identity to
  the target, and extract as related_chain objects with new pdb hierarchies.
  """
  import mmtbx.alignment
  import iotbx.pdb.hierarchy
  results = []
  for chain in pdb_hierarchy.only_model().chains():
    if (chain.is_protein()):
      chain_seq = chain.as_padded_sequence(pad=True,
        substitute_unknown='X',
        pad_at_start=False)
      alignment = mmtbx.alignment.align(
        seq_a=chain_seq,
        seq_b=sequence).extract_alignment()
      identity = alignment.calculate_sequence_identity(skip_chars=['X'])
      if (identity >= min_identity - min_identity_epsilon):
        root = iotbx.pdb.hierarchy.root()
        model = iotbx.pdb.hierarchy.model()
        root.append_model(model)
        chain_id = chain.id
        if (reset_chain_id is not None):
          chain_id = reset_chain_id
        new_chain = iotbx.pdb.hierarchy.chain(id=chain_id)
        model.append_chain(new_chain)
        for residue_group in chain.residue_groups():
          atom_groups = residue_group.atom_groups()
          new_rg = iotbx.pdb.hierarchy.residue_group(
            resseq=residue_group.resseq,
            icode=residue_group.icode)
          for k, atom_group in enumerate(atom_groups):
            if ((remove_alt_confs) and
                ((k > 0) or (not atom_group.altloc.strip() in ['','A']))):
              continue
            new_rg.append_atom_group(atom_group.detached_copy())
          new_chain.append_residue_group(new_rg)
        xrs = root.extract_xray_structure()
        xrs.convert_to_isotropic()
        root.atoms().set_adps_from_scatterers(xrs.scatterers(),
          xrs.unit_cell())
        results.append(
          related_chain(
            source_info=source_info,
            chain_id=chain.id,
            pdb_hierarchy=root,
            identity=identity))
        if (first_chain_only) : break
  return results

class extract_related_models(object):
  """
  Multiprocessing wrapper for pulling matching chains out of a collection of
  PDB hierarchies (either from local files or from a PDB mirror).
  """
  def __init__(self,
      sources_and_models,
      sequence,
      first_chain_only=False,
      reset_chain_id=None,
      min_identity=0.95,
      nproc=Auto,
      log=null_out()):
    adopt_init_args(self, locals())
    self._results = easy_mp.pool_map(
      fixed_func=self.examine_model,
      iterable=range(len(self.sources_and_models)),
      processes=self.nproc)
    for k, results in enumerate(self._results):
      if (len(results) == 0):
        print("  no matches for %s" % self.sources_and_models[k][0], file=log)

  def results(self):
    results = []
    for result_list in self._results : results.extend(result_list)
    return results

  def examine_model(self, i_model):
    source_info, pdb_hierarchy = self.sources_and_models[i_model]
    related = find_similar_chains(
      pdb_hierarchy=pdb_hierarchy,
      sequence=self.sequence,
      source_info=source_info,
      min_identity=self.min_identity,
      first_chain_only=self.first_chain_only,
      reset_chain_id=self.reset_chain_id,
      log=None) #log)
    return related

def load_related_models_in_directory(
      dir_name,
      sequence,
      limit_extensions=True,
      recursive=False,
      first_chain_only=False,
      reset_chain_id=None,
      min_identity=0.95,
      return_original_files=False,
      nproc=Auto,
      log=null_out()):
  """Combines load_all_models_in_directory and extract_related_models."""
  file_names_and_objects = load_all_models_in_directory(
    dir_name=dir_name,
    limit_extensions=limit_extensions,
    recursive=recursive)
  sources_and_models = [ (file_name, file_object.hierarchy)
      for file_name, file_object in file_names_and_objects ]
  related = extract_related_models(
    sources_and_models=sources_and_models,
    sequence=sequence,
    min_identity=min_identity,
    first_chain_only=first_chain_only,
    reset_chain_id=reset_chain_id,
    nproc=nproc,
    log=log).results()
  return related

def load_related_models_from_pdb(pdb_ids, **kwds):
  """
  Combines load_pdb_models and extract_related_models, with multiprocessing
  support.  Functionally equivalent to load_related_models_in_directory but
  with models from the PDB (local or remote).
  """
  sources_and_models = load_pdb_models(pdb_ids=pdb_ids,
    log=kwds.get('log', null_out()))
  return extract_related_models(sources_and_models, **kwds).results()

def selection_by_symmetry(crystal_symmetries,
    source_infos,
    exclude_space_group=None,
    only_space_group=None,
    log=null_out()):
  """
  Return a list of indices for those crystal symmetries which either match the
  only_space_group parameter, or do not match exclude_space_group.
  """
  from cctbx import sgtbx
  if (type(exclude_space_group).__name__ != 'info'):
    exclude_space_group = sgtbx.space_group_info(exclude_space_group)
  if (type(only_space_group).__name__ != 'info'):
    only_space_group = sgtbx.space_group_info(only_space_group)
  selection = []
  for k, symm in enumerate(symmetries):
    skip = False
    sg = None
    if (only_space_group is not None):
      if (symm is None):
        skip = True
      else :
        sg = str(symm.space_group_info())
        if (sg != str(only_space_group)):
          skip = True
    elif (exclude_space_group is not None):
      if (symm is not None):
        sg = str(symm.space_group_info())
        if (sg == str(exclude_space_group)):
          skip = True
    if (skip):
      print("  %s is in space group %s, skipping" %(source_infos[k],sg), file=log)
    else :
      selection.append(k)
  return selection

# FIXME needs to be much smarter about residue matching
class ensembler(object):
  """
  Superpose a collection of single-chain models on a reference chain.  This is
  intentionally very limited in function, but flexible with respect to the
  choice of atoms to superpose.
  """
  def __init__(self,
      reference_hierarchy,
      related_chains,
      atom_selection_string=None,
      backbone_only=Auto,
      calpha_only=None,
      nproc=Auto,
      sieve_fit=False,
      frac_discard=0.5,
      log=null_out()):
    adopt_init_args(self, locals())
    self.realign(atom_selection_string=atom_selection_string)

  def realign(self, atom_selection_string):
    self.atom_selection_string = atom_selection_string
    ref_atoms = self.reference_hierarchy.atoms()
    ref_atoms.reset_i_seq()
    if (self.atom_selection_string is not None):
      assert (not self.calpha_only) and (self.backbone_only != True)
    elif (self.calpha_only):
      self.atom_selection_string = "name CA"
    elif (self.backbone_only):
      self.atom_selection_string = \
        "name CA or name CB or name C or name N or name O"
    else :
      self.atom_selection_string = "all"
    from scitbx.array_family import flex
    sel_cache = self.reference_hierarchy.atom_selection_cache()
    self.atom_selection = sel_cache.selection(self.atom_selection_string)
    assert (self.atom_selection.count(True) > 0)
    self.atoms_ref = []
    ref_sel = flex.size_t()
    ref_chain = self.reference_hierarchy.only_model().only_chain()
    for residue_group in ref_chain.residue_groups():
      rg_atoms = residue_group.only_atom_group().atoms()
      for atom in rg_atoms :
        if (not self.atom_selection[atom.i_seq]):
          continue
        resid = residue_group.resid()
        self.atoms_ref.append("%s %s" % (resid, atom.name.strip()))
        ref_sel.append(atom.i_seq)
    assert (len(ref_sel) > 0)
    sites_ref = ref_atoms.extract_xyz()
    self.reference_sites = sites_ref.select(ref_sel)
    sites_moved = easy_mp.pool_map(
      iterable=range(len(self.related_chains)),
      fixed_func=self.align_model,
      processes=self.nproc)
    self.selection_moved = []
    for k, lsq_fit in enumerate(sites_moved):
      hierarchy = self.related_chains[k].pdb_hierarchy
      if (lsq_fit is None):
        print("No LSQ fit for model %s:%s" % \
          (self.related_chains[k].source_info, self.related_chains[k].chain_id), file=self.log)
        continue
      self.selection_moved.append(k)
      pdb_atoms = hierarchy.atoms()
      sites_cart = lsq_fit.r.elems * pdb_atoms.extract_xyz() + lsq_fit.t.elems
      pdb_atoms.set_xyz(sites_cart)

  # FIXME this entire function is very high on the list of things that need to
  # be replaced by a new, more flexible superposition module
  def align_model(self, i_model):
    from scitbx.array_family import flex
    from scitbx.math import superpose
    hierarchy_moving = self.related_chains[i_model].pdb_hierarchy
    mov_atoms = hierarchy_moving.atoms()
    mov_atoms.reset_i_seq()
    sel_cache = hierarchy_moving.atom_selection_cache()
    mov_atom_selection = sel_cache.selection(self.atom_selection_string)
    mov_chain = hierarchy_moving.only_model().only_chain()
    sel_ref = flex.size_t()
    sel_mov = flex.size_t()
    for residue_group in mov_chain.residue_groups():
      for atom in residue_group.only_atom_group().atoms():
        if (not mov_atom_selection[atom.i_seq]):
          continue
        resid = residue_group.resid()
        ref_name = "%s %s" % (resid, atom.name.strip())
        if (ref_name in self.atoms_ref):
          sel_mov.append(atom.i_seq)
          sel_ref.append(self.atoms_ref.index(ref_name))
    if (len(sel_ref) == 0):
      assert (self.atom_selection_string is not None)
      return None
    assert (len(sel_ref) > 0) and (len(sel_ref) == len(sel_mov))
    xyz_mov = mov_atoms.extract_xyz()
    sites_mov = xyz_mov.select(sel_mov)
    sites_ref = self.reference_sites.select(sel_ref)
    if (self.sieve_fit):
      return superpose.sieve_fit(
        sites_fixed=sites_ref,
        sites_moving=sites_mov,
        frac_discard=self.frac_discard)
    else :
      return superpose.least_squares_fit(
        reference_sites=sites_ref,
        other_sites=sites_mov)

  def as_multi_model_hierarchy(self):
    hierarchies = []
    for k in self.selection_moved :
      hierarchies.append(self.related_chains[k].pdb_hierarchy)
    return join_hierarchies(hierarchies)

def join_hierarchies(hierarchies, model_ids=None):
  import iotbx.pdb.hierarchy
  root = iotbx.pdb.hierarchy.root()
  for k, hierarchy in enumerate(hierarchies):
    model_id = k + 1
    if (model_ids is not None):
      model_id = model_ids[k]
    model = iotbx.pdb.hierarchy.model(id=str(model_id))
    root.append_model(model)
    chain = hierarchy.only_model().only_chain().detached_copy()
    model.append_chain(chain)
  return root

def extract_and_superpose(
    reference_hierarchy,
    search_directory,
    sequence,
    params,
    nproc=Auto,
    out=sys.stdout):
  if (search_directory is None):
    pdb_ids = fetch_similar_pdb_ids(
      sequence=sequence,
      min_identity=params.min_identity,
      min_resolution=params.min_resolution,
      ligands=params.ligands,
      xray_only=params.xray_only,
      log=out)
    if (params.exclude_ids is not None):
      exclude_ids = [ s.lower().split(",") for s in params.exclude_ids ]
      for other_id in exclude_ids :
        if (other_id in pdb_ids):
          pdb_ids.remove(other_id)
    if (len(pdb_ids) == 0):
      return None
    else :
      print("%d PDB IDs retrieved:" % len(pdb_ids), file=out)
      i = 0
      while (i < len(pdb_ids)):
        print("  %s" % " ".join(pdb_ids[i:i+16]), file=out)
        i += 16
    related_chains = load_related_models_from_pdb(
      pdb_ids=pdb_ids,
      sequence=sequence,
      first_chain_only=params.single_chain,
      min_identity=params.min_identity,
      reset_chain_id=reference_hierarchy.only_model().only_chain().id,
      nproc=nproc,
      log=out)
  else : # TODO search_directory support
    raise NotImplementedError()
  assert (related_chains is not None)
  return ensembler(
    reference_hierarchy=reference_hierarchy,
    related_chains=related_chains,
    nproc=params.nproc)

def extract_peptide_fragments_by_sequence(
    pdb_hierarchy,
    sequence,
    renumber_from=None,
    mainchain_only=False,
    remove_hydrogens=True,
    reset_chain_id='A',
    reset_icode=' '):
  from iotbx.pdb import amino_acid_codes
  import iotbx.pdb.hierarchy
  aa_3_as_1 = amino_acid_codes.one_letter_given_three_letter
  fragments = []
  def make_hierarchy():
    root = iotbx.pdb.hierarchy.root()
    model = iotbx.pdb.hierarchy.model()
    chain = iotbx.pdb.hierarchy.chain(id=reset_chain_id)
    root.append_model(model)
    model.append_chain(chain)
    return root
  def find_next_matching_atom_group(residue_groups, i_res, i_seq):
    for next_atom_group in residue_groups[i_res].atom_groups():
      seq_code = aa_3_as_1.get(next_atom_group.resname)
      if (seq_code == sequence[i_seq]) or (sequence[i_seq] in ['.','X']):
        return next_atom_group
    return None
  for chain in pdb_hierarchy.only_model().chains():
    if chain.is_protein():
      i_res = 0
      residue_groups = chain.residue_groups()
      while (i_res < len(residue_groups)):
        atom_group = find_next_matching_atom_group(residue_groups,i_res,0)
        if (atom_group is not None):
          fragment = [atom_group]
          for i_seq in range(1, len(sequence)):
            j_res = i_res + i_seq
            if (j_res >= len(residue_groups)):
              fragment = []
              break
            next_atom_group = find_next_matching_atom_group(residue_groups,
              i_res=j_res, i_seq=i_seq)
            if (next_atom_group is not None):
              fragment.append(next_atom_group)
            else :
              fragment = []
              break
          if (len(fragment) > 0):
            assert len(fragment) == len(sequence)
            root = make_hierarchy()
            new_chain = root.only_model().only_chain()
            for i_group, next_atom_group in enumerate(fragment):
              resseq = next_atom_group.parent().resseq
              icode = next_atom_group.parent().icode
              if (renumber_from is not None):
                resseq = "%4d" % (renumber_from + i_group)
              if (reset_icode is not None):
                icode = reset_icode
              new_residue_group = iotbx.pdb.hierarchy.residue_group(
                resseq=resseq,
                icode=icode)
              new_chain.append_residue_group(new_residue_group)
              new_residue_group.append_atom_group(
                next_atom_group.detached_copy())
            sel = None
            sel_cache = root.atom_selection_cache()
            if (mainchain_only):
              sel = sel_cache.selection(
                "name C or name O or name N or name CA or name CB")
            elif (remove_hydrogens):
              sel = sel_cache.selection("not (element H or element D)")
            if (sel is not None):
              assert (sel.count(True) > 0)
              root = root.select(sel)
            fragments.append(root)
        i_res += 1
  return fragments
