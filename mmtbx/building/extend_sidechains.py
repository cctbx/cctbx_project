
from __future__ import division
from scitbx.matrix import col
from libtbx.str_utils import make_sub_header
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

def correct_sequence (pdb_hierarchy,
    sequences,
    truncate_to_cbeta=False,
    out=sys.stdout) :
  """
  Modify the sequence for the pdb hierarchy to match that of the aligned
  sequence.  This will remove incompatible atoms; the sidechains will still
  need to be extended separated.  For proteins only - mismatches in nucleic
  acids will only result in a warning.

  :param pdb_hierarchy: iotbx.pdb.hierarchy.root object
  :param sequences: list of iotbx.bioinformatics.sequence objects
  :param trucate_to_cbeta: chop off entire sidechain to C-beta (default: leave
                           common atoms in place)
  :param out: output filehandle (default = stdout)
  :returns: number of atom_group objects renamed
  """
  from mmtbx.monomer_library import idealized_aa
  import mmtbx.validation.sequence
  from iotbx.pdb.amino_acid_codes import three_letter_given_one_letter
  seq_validation = mmtbx.validation.sequence.validation(
    pdb_hierarchy=pdb_hierarchy,
    sequences=sequences,
    log=out)
  for chain_seq in seq_validation.chains :
    if (chain_seq.chain_type == mmtbx.validation.sequence.NUCLEIC_ACID) :
      if (len(chain_seq.mismatch) > 0) :
        print >> out, \
          "  WARNING: will skip %d mismatches in nucleic acid chain '%s'" % \
          chain_seq.chain_id
  res_dict = idealized_aa.residue_dict()
  expected_names = {}
  for resname in res_dict.keys() :
    if (not "_h" in resname) :
      ideal_res = res_dict[resname]
      expected_names[resname] = set([ a.name for a in ideal_res.atoms() ])
  n_changed = 0
  for chain in pdb_hierarchy.only_model().chains() :
    if (not chain.is_protein()) :
      continue
    for chain_seq in seq_validation.chains :
      if (chain.id == chain_seq.chain_id) and (len(chain_seq.mismatch) > 0) :
        for residue_group in chain.residue_groups() :
          resid = residue_group.resid()
          if (resid in chain_seq.mismatch) :
            idx = chain_seq.mismatch.index(resid)
            new_code = chain_seq.actual_code[idx]
            new_resname = three_letter_given_one_letter.get(new_code)
            if (new_resname is not None) :
              expected_atoms = expected_names[new_resname.lower()]
              if (truncate_to_cbeta) :
                expected_atoms = expected_names["ala"]
              for atom_group in residue_group.atom_groups() :
                n_changed += 1
                n_removed = 0
                atom_group.resname = new_resname
                for atom in atom_group.atoms() :
                  if (not atom.name in expected_atoms) :
                    atom_group.remove_atom(atom)
                    n_removed += 1
              print >> out, "  chain '%s' %s %s --> %s (%d atoms removed)" % \
                (chain.id, resid, residue_group.atom_groups()[0].resname,
                 new_resname, n_removed)
  pdb_hierarchy.atoms().reset_i_seq()
  return n_changed

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
    is_amino_acid=True, # XXX this could probably be smarter...
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
      sites_cart=new_residue.atoms().extract_xyz())
    scorer.apply_final()
  return new_residue

def extend_protein_model (pdb_hierarchy,
    selection=None,
    hydrogens=Auto,
    max_atoms_missing=None,
    log=None,
    modify_segids=True,
    prefilter_callback=None,
    idealized_residue_dict=None,
    skip_non_protein_chains=True) :
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
  ideal_dict = idealized_residue_dict
  if (ideal_dict is None) :
    ideal_dict = idealized_aa.residue_dict()
  if (log is None) : log = null_out()
  mon_lib_srv = mmtbx.monomer_library.server.server()
  pdb_atoms = pdb_hierarchy.atoms()
  if (selection is None) :
    selection = flex.bool(pdb_atoms.size(), True)
  partial_sidechains = []
  for chain in pdb_hierarchy.only_model().chains() :
    if (not chain.is_protein()) and (skip_non_protein_chains) :
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
      if (idealized_residue_dict is None) :
        res_class = common_residue_names_get_class(residue.resname)
        if (res_class != "common_amino_acid") :
          print >> log, "    skipping non-standard residue %s" % residue.resname
          continue
      else :
        key = residue.resname.lower()
        if (hydrogens == True) :
          key = key + "_h"
        if (not key in idealized_residue_dict.keys()) :
          pass
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
    allow_modified_residues=False,
    out=sys.stdout) :
  """
  Use real-space refinement tools to fit newly extended residues.
  """
  from mmtbx.refinement.real_space import fit_residue
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
  make_sub_header("Processing new model", out=ppdb_out)
  processed_pdb = building.reprocess_pdb(
    pdb_hierarchy=pdb_hierarchy,
    crystal_symmetry=fmodel.f_obs().crystal_symmetry(),
    cif_objects=cif_objects,
    out=ppdb_out)
  print >> ppdb_out, ""
  hierarchy = processed_pdb.all_chain_proxies.pdb_hierarchy
  xrs = processed_pdb.xray_structure()
  grm = building.get_restraints_manager(
    processed_pdb_file = processed_pdb,
    xray_structure     = xrs)
  make_sub_header("Fitting residues", out=out)
  target_map = fmodel.map_coefficients(
    map_type="2mFo-DFc",
    exclude_free_r_reflections=True).fft_map(
      resolution_factor=0.25).apply_sigma_scaling().real_map_unpadded()
  unit_cell = xrs.unit_cell()
  for chain in hierarchy.only_model().chains() :
    if (not chain.is_protein()) and (not allow_modified_residues) : continue
    residue_groups = chain.residue_groups()
    for i_rg, residue_group in enumerate(residue_groups) :
      prev_res = next_res = None
      atom_groups = residue_group.atom_groups()
      if (len(atom_groups) > 1) : continue
      residue = atom_groups[0]
      atoms = residue.atoms()
      atoms.reset_tmp()
      segids = atoms.extract_segid()
      if (segids.all_eq("XXXX")) :
        sites_start = atoms.extract_xyz()
        def get_two_fofc_mean (residue) :
          sum = 0
          n_atoms = 0
          for atom in residue.atoms() :
            if (not atom.element.strip() in ["H","D"]) :
              site_frac = unit_cell.fractionalize(site_cart=atom.xyz)
              sum += target_map.eight_point_interpolation(site_frac)
              n_atoms += 1
          assert (n_atoms > 0)
          return sum / n_atoms
        sites_start = atoms.extract_xyz().deep_copy()
        two_fofc_mean_start = get_two_fofc_mean(residue)
        refit = fit_residue.run(
          target_map=target_map,
          residue=residue,
          xray_structure=xrs,
          mon_lib_srv=mon_lib_srv,
          rotamer_manager=rotamer_manager,
          geometry_restraints_manager=grm,
          real_space_gradients_delta=fmodel.f_obs().d_min()*0.25,
          rms_bonds_limit=0.01,
          rms_angles_limit=1.0,
          backbone_sample_angle=20,
          allow_modified_residues=allow_modified_residues)
        two_fofc_mean_end = get_two_fofc_mean(residue)
        sites_end = atoms.extract_xyz()
        flag = ""
        if (two_fofc_mean_end > two_fofc_mean_start) :
          flag = " <-- keep"
          xrs = refit.xray_structure
        else :
          atoms.set_xyz(sites_start)
          for atom in atoms :
            atom.tmp = 1
        print >> out, \
          "    residue '%s' : rmsd=%5.3f 2fofc_start=%5.3f 2fofc_end=%5.3f%s" \
          % (residue.id_str(), sites_end.rms_difference(sites_start),
            two_fofc_mean_start, two_fofc_mean_end, flag)
  return hierarchy, xrs

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

class extend_and_refine (object) :
  """
  Run the combined sidechain extension and real-space fitting, and optionally
  write final model and map coefficients.
  """
  def __init__ (self,
      pdb_hierarchy,
      xray_structure,
      fmodel,
      params,
      cif_objects=(),
      out=sys.stdout,
      output_model=None,
      output_map_coeffs=None,
      prefix=None,
      write_files=True,
      reset_segid=True,
      verbose=True) :
    if (write_files) :
      assert ((prefix is not None) or
              (not None in [output_model,output_map_coeffs]))
    if (params.build_hydrogens is Auto) :
      params.build_hydrogens = xray_structure.hd_selection().count(True) > 0
    make_sub_header("Filling in partial sidechains", out=out)
    prefilter_callback = prefilter(
      fmodel=fmodel,
      out=out)
    n_atoms_start = xray_structure.sites_cart().size()
    self.n_extended = extend_protein_model(
      pdb_hierarchy=pdb_hierarchy,
      hydrogens=params.build_hydrogens,
      max_atoms_missing=params.max_atoms_missing,
      prefilter_callback=prefilter_callback,
      log=out)
    print >> out, "  %d sidechains extended." % self.n_extended
    if (self.n_extended > 0) and (not params.skip_rsr) :
      pdb_hierarchy, xray_structure = refit_residues(
        pdb_hierarchy=pdb_hierarchy,
        cif_objects=cif_objects,
        fmodel=fmodel,
        use_rotamers=params.use_rotamers,
        anneal=params.anneal_residues,
        out=out)
    else :
      xray_structure = pdb_hierarchy.extract_xray_structure(
        crystal_symmetry=xray_structure)
    fmodel.update_xray_structure(xray_structure, update_f_calc=True)
    n_atoms_end = xray_structure.sites_cart().size()
    self.r_work = fmodel.r_work()
    self.r_free = fmodel.r_free()
    self.n_new_atoms = n_atoms_end - n_atoms_start
    self.pdb_file = self.map_file = None
    if reset_segid :
      for atom in pdb_hierarchy.atoms() :
        atom.segid = ""
    if (write_files) :
      if (output_model is None) :
        output_model = prefix + "_extended.pdb"
      f = open(output_model, "w")
      f.write(pdb_hierarchy.as_pdb_string(fmodel.xray_structure))
      f.close()
      self.pdb_file = output_model
      print >> out, "  wrote new model to %s" % output_model
      if (output_map_coeffs is None) :
        output_map_coeffs = prefix + "_maps.mtz"
      from mmtbx.maps.utils import get_maps_from_fmodel
      import iotbx.map_tools
      two_fofc_map, fofc_map = get_maps_from_fmodel(fmodel)
      iotbx.map_tools.write_map_coeffs(
        fwt_coeffs=two_fofc_map,
        delfwt_coeffs=fofc_map,
        file_name=output_map_coeffs)
      print >> out, "  wrote map coefficients to %s" % output_map_coeffs
      self.map_file = output_map_coeffs
