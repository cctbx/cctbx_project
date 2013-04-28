
from __future__ import division
from scitbx.matrix import col
from libtbx.utils import null_out
from libtbx import Auto
import sys

master_params = """
selection = None
  .type = atom_selection
build_hydrogens = Auto
  .type = bool
max_atoms_missing = None
  .type = int
use_rotamers = True
  .type = bool
anneal_residues = False
  .type = bool
skip_rsr = False
  .type = bool
"""

class conformation_scorer (object) :
  """
  Stand-in for the conformation scoring class in mmtbx.refinement.real_space;
  instead of calculating fit to the map, this simply uses the change in
  position of the first atom being moved at each rotation.  This allows us to
  superimpose the conformations for those atoms which are present in both the
  old and the new residues.
  """
  def __init__ (self, old_residue, new_residue) :
    from scitbx.array_family import flex
    old_residue_atoms = old_residue.atoms()
    self.new_residue_atoms = new_residue.atoms()
    n_atoms = self.new_residue_atoms.size()
    self.new_residue_selection = flex.bool(n_atoms, False)
    self.selection_mappings = flex.size_t(n_atoms, 0)
    for i_seq, old_atom in enumerate(old_residue_atoms) :
      for j_seq, new_atom in enumerate(self.new_residue_atoms) :
        if (old_atom.name == new_atom.name) :
          self.new_residue_selection[j_seq] = True
          self.selection_mappings[j_seq] = i_seq
    self.sites_old = old_residue_atoms.extract_xyz()
    self.sites_cart = self.new_residue_atoms.extract_xyz()
    self.dist_min = None

  def update (self, sites_cart, selection) :
    first_i_seq = selection[0]
    if (self.new_residue_selection[first_i_seq]) :
      dist = abs(col(sites_cart[first_i_seq]) -
        col(self.sites_old[self.selection_mappings[first_i_seq]]))
      if (dist < self.dist_min) :
        self.sites_cart = sites_cart
        self.dist_min = dist
        return True
    return False

  def reset_with (self, sites_cart, selection) :
    self.sites_cart = sites_cart
    first_i_seq = selection[0]
    if (self.new_residue_selection[first_i_seq]) :
      self.dist_min = abs(col(sites_cart[first_i_seq]) -
        col(self.sites_old[self.selection_mappings[first_i_seq]]))
    else :
      self.dist_min = sys.maxint

  def apply_final (self) :
    self.new_residue_atoms.set_xyz(self.sites_cart)

def extend_residue (residue,
    ideal_dict,
    mon_lib_srv=None,
    hydrogens=False,
    match_conformation=True) :
  """
  Rebuild a sidechain by substituting an ideal amino acid and (optionally)
  rotating the sidechain to match the old conformation as closely as possible.
  """
  from iotbx.pdb import hierarchy
  res_key = residue.resname.lower()
  if (hydrogens == True) : res_key += "_h"
  ideal = ideal_dict[res_key]
  tmp_residue = residue.detached_copy()
  new_residue = hierarchy.substitute_atom_group(
    current_group=tmp_residue,
    new_group=ideal.only_model().only_chain().only_residue_group().only_atom_group(),
    backbone_only=True,
    exclude_hydrogens=False)
  assert (new_residue.resname == residue.resname)
  if (match_conformation) :
    from mmtbx.building import generate_sidechain_clusters
    import mmtbx.refinement.real_space
    assert (mon_lib_srv is not None)
    clusters = generate_sidechain_clusters(residue=new_residue,
      mon_lib_srv=mon_lib_srv)
    scorer = conformation_scorer(
      old_residue=residue,
      new_residue=new_residue)
    mmtbx.refinement.real_space.torsion_search(
      scorer=scorer,
      clusters=clusters,
      sites_cart=new_residue.atoms().extract_xyz(),
      start=0,
      stop=355,
      step=5)
    scorer.apply_final()
  return new_residue

def extend_protein_model (pdb_hierarchy,
    selection=None,
    hydrogens=Auto,
    max_atoms_missing=None,
    log=None,
    modify_segids=True,
    prefilter_callback=None) :
  """
  Replace all sidechains with missing non-hydrogen atoms in a PDB hierarchy.
  """
  from mmtbx.monomer_library import idealized_aa
  from mmtbx.rotamer import rotamer_eval
  import mmtbx.monomer_library.server
  from iotbx.pdb import common_residue_names_get_class
  from scitbx.array_family import flex
  if (prefilter_callback is not None) :
    assert hasattr(prefilter_callback, "__call__")
  else :
    prefilter_callback = lambda r: True
  ideal_dict = idealized_aa.residue_dict()
  if (log is None) : log = null_out()
  mon_lib_srv = mmtbx.monomer_library.server.server()
  pdb_atoms = pdb_hierarchy.atoms()
  if (selection is None) :
    selection = flex.bool(pdb_atoms.size(), True)
  partial_sidechains = []
  for chain in pdb_hierarchy.only_model().chains() :
    if (not chain.is_protein()) :
      print >> log, "    skipping non-protein chain '%s'" % chain.id
      continue
    for residue_group in chain.residue_groups() :
      atom_groups = residue_group.atom_groups()
      if (len(atom_groups) > 1) :
        print >> log, "    %s %s has multiple conformations, skipping" % \
          (chain.id, residue_group.resid())
        continue
      residue = atom_groups[0]
      i_seqs = residue.atoms().extract_i_seq()
      residue_sel = selection.select(i_seqs)
      if (not residue_sel.all_eq(True)) :
        continue
      res_class = common_residue_names_get_class(residue.resname)
      if (res_class != "common_amino_acid") :
        print >> log, "    skipping non-standard residue %s" % residue.resname
        continue
      missing_atoms = rotamer_eval.eval_residue_completeness(
        residue=residue,
        mon_lib_srv=mon_lib_srv,
        ignore_hydrogens=True)
      if (len(missing_atoms) > 0) :
        print >> log, "    missing %d atoms in %s: %s" % (len(missing_atoms),
          residue.id_str(), ",".join(missing_atoms))
        if ((max_atoms_missing is None) or
            (len(missing_atoms) < max_atoms_missing)) :
          if (prefilter_callback(residue)) :
            partial_sidechains.append(residue)
  for residue in partial_sidechains :
    new_residue = extend_residue(residue=residue,
      ideal_dict=ideal_dict,
      hydrogens=hydrogens,
      mon_lib_srv=mon_lib_srv,
      match_conformation=True)
    if (modify_segids) :
      for atom in new_residue.atoms() :
        atom.segid = "XXXX"
    rg = residue.parent()
    rg.remove_atom_group(residue)
    rg.append_atom_group(new_residue.detached_copy())
  pdb_hierarchy.atoms().reset_i_seq()
  pdb_hierarchy.atoms().reset_serial()
  return len(partial_sidechains)

def refit_residues (
    pdb_hierarchy,
    cif_objects,
    fmodel,
    use_rotamers=True,
    anneal=False,
    verbose=True,
    out=sys.stdout) :
  from mmtbx.rotamer import rotamer_eval
  import mmtbx.monomer_library.server
  from mmtbx import building
  from scitbx.array_family import flex
  mon_lib_srv = mmtbx.monomer_library.server.server()
  rotamer_manager = rotamer_eval.RotamerEval()
  ppdb_out = box_out = out
  if (not verbose) :
    ppdb_out = null_out()
    box_out = null_out()
  processed_pdb = building.reprocess_pdb(
    pdb_hierarchy=pdb_hierarchy,
    crystal_symmetry=fmodel.f_obs().crystal_symmetry(),
    cif_objects=cif_objects,
    out=ppdb_out)
  target_map = fmodel.map_coefficients(
    map_type="2mFo-DFc",
    exclude_free_r_reflections=True).fft_map(
      resolution_factor=0.25).apply_sigma_scaling().real_map_unpadded()
  hierarchy = processed_pdb.all_chain_proxies.pdb_hierarchy
  xrs = processed_pdb.xray_structure()
  unit_cell = xrs.unit_cell()
  for chain in hierarchy.only_model().chains() :
    if (not chain.is_protein()) : continue
    residue_groups = chain.residue_groups()
    for i_rg, residue_group in enumerate(residue_groups) :
      atom_groups = residue_group.atom_groups()
      if (len(atom_groups) > 1) : continue
      residue = atom_groups[0]
      atoms = residue.atoms()
      segids = atoms.extract_segid()
      if (segids.all_eq("XXXX")) :
        sites_start = atoms.extract_xyz()
        iselection = atoms.extract_i_seq()
        assert (not iselection.all_eq(0))
        selection = flex.bool(xrs.scatterers().size(), False)
        selection.set_selected(iselection, True)
        # FIXME this doesn't really do what I want, unfortunately
        box = building.box_build_refine_base(
          pdb_hierarchy=hierarchy,
          xray_structure=xrs,
          processed_pdb_file=processed_pdb,
          target_map=target_map,
          selection=selection,
          d_min=fmodel.f_obs().d_min(),
          out=null_out())
        cc_start = box.cc_model_map()
        if (use_rotamers) :
          box.fit_residue_in_box(
            mon_lib_srv=mon_lib_srv,
            rotamer_manager=rotamer_manager)
        else :
          box.real_space_refine()
        if (anneal) :
          box.anneal()
          box.real_space_refine()
        cc_end = box.cc_model_map()
        flag = ""
        sites_backup = xrs.sites_cart().deep_copy()
        sites_cart = box.update_original_coordinates()
        hierarchy.atoms().set_xyz(sites_cart)
        sites_end = residue.atoms().extract_xyz()
        if (cc_end > cc_start) :
          flag = " <-- keep"
          xrs.set_sites_cart(sites_cart)
        else :
          hierarchy.atoms().set_xyz(sites_backup)
        print >> out, \
          "    residue '%s' : rmsd=%5.3f cc_start=%5.3f cc_end=%5.3f%s" % \
          (residue.id_str(), sites_end.rms_difference(sites_start), cc_start,
           cc_end, flag)
  return hierarchy

class prefilter (object) :
  """
  Optional filter for excluding residues with poor backbone density from being
  extended.  This is done as a separate callback to enable the main rebuilding
  routine to be independent of data/maps.
  """
  def __init__ (self, fmodel, out, backbone_min_sigma=1.0) :
    target_map = fmodel.map_coefficients(
      map_type="2mFo-DFc",
      exclude_free_r_reflections=True).fft_map(
        resolution_factor=0.25).apply_sigma_scaling().real_map_unpadded()
    self.unit_cell = fmodel.f_obs().unit_cell()
    self.map = target_map
    self.out = out
    self.backbone_min_sigma = backbone_min_sigma

  def __call__ (self, residue) :
    atoms = residue.atoms()
    sigma_mean = n_bb = 0
    for atom in atoms :
      if (atom.name.strip() in ["N","C","CA", "CB"]) :
        site_frac = self.unit_cell.fractionalize(site_cart=atom.xyz)
        sigma_mean += self.map.eight_point_interpolation(site_frac)
        n_bb += 1
    if (n_bb > 1) :
      sigma_mean /= n_bb
    if (sigma_mean < self.backbone_min_sigma) :
      print >> self.out, "      *** poor backbone density, skipping"
      return False
    return True
