
# TODO tests!

"""
Post-fitting cleanup of ligand positions to match NCS operations present in
protein model.  This can be used to recover cases where one copy is placed
successfully but another misses due to interfering protein atoms, weak density,
false positive in empty protein density, etc.  It can also accelerate ligand
placement for large structures if at least one copy is already placed with high
confidence.  Should only be used when the initial placement (e.g. LigandFit CC)
is sufficiently good.

Usually this will be called via the command mmtbx.apply_ncs_to_ligand.
"""

from __future__ import absolute_import, division, print_function
from libtbx.utils import Sorry, null_out
from libtbx.str_utils import make_sub_header
from libtbx import adopt_init_args, group_args
import copy
from math import sqrt
import sys
from six.moves import zip

debug = True

ncs_ligand_phil = """
d_min = 2.5
  .type = float
max_rmsd = 2.0
  .type = float
min_cc = 0.7
  .type = float
min_cc_reference = 0.85
  .type = float
min_2fofc = 1.0
  .type = float
min_dist_center = 5
  .type = float
remove_clashing_atoms = True
  .type = bool
clash_cutoff = 2.0
  .type = float
write_sampled_pdbs = False
  .type = bool
"""

def resid_str(atom):
  labels = atom.fetch_labels()
  return "%s %s" % (labels.resname, labels.resid())

def xyz_distance(xyz1, xyz2):
  x1,y1,z1 = xyz1
  x2,y2,z2 = xyz2
  return sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)

class group_operators(object):
  """
  Object for storing information about NCS relationships relative to an
  atom selection, i.e. given four identical chains A, B, C, and D, this might
  store information about chain B and the NCS operators for approximate
  superpositioning on chains A, C, and D.
  """
  def __init__(self, selection, sele_str, sites_cart):
    self.selection = selection
    self.selection_string = sele_str
    self.center_of_mass = sites_cart.select(selection).mean()
    self.operators = []
    self.op_selections = []

  def add_operator(self, ops, sele_str):
    """
    Save an NCS operator and corresponding atom selection.
    """
    self.operators.append(ops)
    self.op_selections.append(sele_str)

  def distance_from_center(self, sites):
    """
    Compute the distance between centers-of-mass of this selection and the
    given sites.  Used to determine ligand-protein chain relationships.
    """
    c_o_m = sites.mean()
    return xyz_distance(c_o_m, self.center_of_mass)

  def show_summary(self, out=None, prefix=""):
    """
    Print out selections and NCS operators.
    """
    if (out is None) : out = sys.stdout
    print(prefix+"Reference selection: %s" % self.selection_string, file=out)
    for op, op_sele in zip(self.operators, self.op_selections):
      print(prefix+"  Selection: %s" % op_sele, file=out)
      print(prefix+"  Rotation:", file=out)
      print(prefix+"    %6.4f  %6.4f  %6.4f" % op.r.elems[0:3], file=out)
      print(prefix+"    %6.4f  %6.4f  %6.4f" % op.r.elems[3:6], file=out)
      print(prefix+"    %6.4f  %6.4f  %6.4f" % op.r.elems[6:9], file=out)
      print(prefix+"  Translation:", file=out)
      print(prefix+"    %6.4f  %6.4f  %6.4f" % op.t.elems, file=out)
      print("", file=out)
    print("", file=out)

class sample_operators(object):
  """
  Determines an appropriate "reference" ligand, and samples the density around
  sites transformed by each operator, applying the best operator to the
  other ligand copies if meeting cutoff criteria.
  """
  def __init__(self,
      pdb_hierarchy,
      fmodel,
      ncs_operators,
      ligands,
      params,
      log=None):
    if (log is None) : log = sys.stdout
    if (params is None):
      params = master_phil().fetch().extract()
    adopt_init_args(self, locals())
    self.xray_structure = fmodel.xray_structure.deep_copy_scatterers()
    xrs_ncs = fmodel.xray_structure.deep_copy_scatterers()
    from iotbx.pdb import hierarchy
    self.setup_maps()
    best_cc = 0
    best_k = -1
    best_ligand = None
    other_ligands = []
    print("Identifying reference ligand...", file=log)
    def show_map_stats(prefix, stats):
      print("   %s: CC = %5.3f  mean = %6.2f" % (prefix, stats.cc,
        stats.map_mean), file=log)
    for k, ligand in enumerate(ligands):
      atoms = ligand.atoms()
      start = self.get_sites_cc(atoms)
      show_map_stats("Ligand %d" % (k+1), start)
      if (start.cc > best_cc) and (start.cc > params.min_cc_reference):
        best_ligand = ligand
        best_k = k
        best_cc = best_cc
    if (best_ligand is None):
      raise Sorry("No ligand with acceptable CC (>%.2f) found." %
        params.min_cc_reference)
    best_atoms = best_ligand.atoms()
    for k, ligand in enumerate(ligands):
      if (ligand is not best_ligand):
        other_ligands.append(ligand)
    print("Copy #%d was the best, using that as reference" % (best_k+1), file=log)
    print("", file=log)
    sites_ref = best_ligand.atoms().extract_xyz()
    min_dist = sys.maxsize
    best_group = None
    shifts = ncs_operators.get_ncs_groups_shifts(
      self.xray_structure.sites_cart(),
      sites_ref
      )
    for i, s in enumerate(shifts):
      dxyz = xyz_distance(s[0], (0,0,0))
      if (dxyz < min_dist):
        best_group = i
        min_dist = dxyz
    array_of_str_selections = ncs_operators.get_array_of_str_selections()[0]
    if (best_group is not None):
      print("This appears to be bound to the selection \"%s\"" % \
        array_of_str_selections[best_group], file=log)
    if best_group==0: pass
    else:
      print('best_group',best_group)
      assert 0
    # always have the first ligand in the master
    self.new_ligands = []
    for j, operator in enumerate(ncs_operators[0].copies):
      new_ligand = best_ligand.detached_copy()
      atoms = new_ligand.atoms()
      sites_new = operator.r.elems * sites_ref + operator.t.elems
      sites_mean = sites_new.mean()
      for other in other_ligands :
        sites_other_mean = other.atoms().extract_xyz().mean()
        dxyz = xyz_distance(sites_other_mean, sites_new.mean())
        if (dxyz < params.min_dist_center):
          print("  operator %d specifies an existing ligand" % (j+1), file=log)
          break
      else :
        atoms.set_xyz(sites_new)
        stats_new = self.get_sites_cc(best_atoms, sites_new)
        show_map_stats("NCS op. %2d" % (j+1), stats_new)
        if (params.write_sampled_pdbs):
          lig_rg = hierarchy.residue_group()
          lig_rg.resseq = j+1
          lig_rg.append_atom_group(new_ligand)
          f = open("ncs_ligand_%d.pdb" % (j+1), "w")
          for atom in new_ligand.atoms():
            f.write(atom.format_atom_record()+"\n")
          f.close()
        # XXX ideally, given multiple high-quality ligand placements, we should
        # probably try sampling NCS operations for all of these and pick the
        # best new CC, rather than assuming that the best starting ligand will
        # superpose best on the density.
        if (stats_new.cc > params.min_cc):
          print("  operator %d has acceptable CC (%.3f)" % (j+1,
            stats_new.cc), file=log)
          self.new_ligands.append(new_ligand)

  def setup_maps(self):
    """
    Create 2mFo-DFc, mFo-DFc, and Fc maps.
    """
    map_helper = self.fmodel.electron_density_map()
    map_coeffs = map_helper.map_coefficients("2mFo-DFc")
    diff_map_coeffs = map_helper.map_coefficients("mFo-DFc")
    fft_map = map_coeffs.fft_map(resolution_factor=1/3.)
    diff_fft_map = diff_map_coeffs.fft_map(resolution_factor=1/3.)
    fft_map.apply_sigma_scaling()
    diff_fft_map.apply_sigma_scaling()
    fcalc = map_helper.map_coefficients("Fc")
    fcalc_map = fcalc.fft_map(resolution_factor=1/3.)
    fcalc_map.apply_sigma_scaling()
    self.unit_cell = map_coeffs.unit_cell()
    self.real_map = fft_map.real_map()
    self.diff_map = diff_fft_map.real_map()
    self.n_real = self.real_map.focus()
    self.m_real = self.real_map.all()
    self.fcalc_real_map = fcalc_map.real_map()

  def get_new_fcalc_map(self, sites_new, i_seqs):
    xrs_new = self.xray_structure.deep_copy_scatterers()
    all_sites = xrs_new.sites_cart()
    all_sites.set_selected(i_seqs, sites_new)
    xrs_new.set_sites_cart(all_sites)
    self.fmodel.update_xray_structure(
      xray_structure=xrs_new,
      update_f_calc=True,
      update_f_mask=True)
    fcalc = self.fmodel.electron_density_map().map_coefficients("Fc")
    fcalc_map = fcalc.fft_map(resolution_factor=1/3.)
    fcalc_map.apply_sigma_scaling()
    # XXX now revert to original xray structure
    self.fmodel.update_xray_structure(
      xray_structure=self.xray_structure,
      update_f_calc=True,
      update_f_mask=True)
    return fcalc_map.real_map()

  def get_sites_cc(self, atoms, sites=None):
    from cctbx import maptbx
    from scitbx.array_family import flex
    radii = flex.double()
    for atom in atoms :
      if (atom.element.strip() in ["H", "D"]):
        radii.append(1.)
      else :
        radii.append(1.5)
    fcalc_map = self.fcalc_real_map
    if (sites is None):
      sites = atoms.extract_xyz()
    else :
      fcalc_map = self.get_new_fcalc_map(
        sites_new=sites,
        i_seqs=atoms.extract_i_seq())
    sel = maptbx.grid_indices_around_sites(
      unit_cell  = self.unit_cell,
      fft_n_real = self.n_real,
      fft_m_real = self.m_real,
      sites_cart = sites,
      site_radii = radii)
    m1 = self.real_map.select(sel)
    m2 = fcalc_map.select(sel)
    cc = flex.linear_correlation(x=m1, y=m2).coefficient()
    return group_args(
      cc=cc,
      map_mean=flex.mean(m1.as_1d()))

def _is_same_residue(atom1, atom2):
  labels1 = atom1.fetch_labels()
  labels2 = atom2.fetch_labels()
  return ((labels1.resid() == labels2.resid()) and
          (labels1.chain_id == labels2.chain_id) and
          (labels1.altloc == labels2.altloc))

def remove_clashing_atoms(
    xray_structure,
    pdb_hierarchy,
    ligands,
    params,
    log):
  """
  Since the transformed ligands will very frequently overlap with existing
  atoms, these need to be deleted if we are confident about the new positions.
  """
  from scitbx.array_family import flex
  pdb_atoms = pdb_hierarchy.atoms()
  selection = flex.bool(pdb_atoms.size(), True)
  pair_asu_table = xray_structure.pair_asu_table(
    distance_cutoff=params.clash_cutoff)
  #xray_structure.show_distances(params.clash_cutoff)
  # XXX for some unknown reason this doesn't work if I use the object returned
  # by pair_asu_table.extract_pair_sym_table() instead
  asu_table = pair_asu_table.table()
  remove_atoms = set([])
  for ligand in ligands :
    i_seqs = ligand.atoms().extract_i_seq()
    for i_seq in i_seqs :
      asu_dict = asu_table[i_seq]
      for j_seq, sym_ops in asu_dict.items():
        if (not _is_same_residue(pdb_atoms[i_seq], pdb_atoms[j_seq])):
          remove_atoms.add(j_seq)
  if (len(remove_atoms) > 0):
    def show_removed(atoms):
      for atom in atoms :
        print("  warning: deleting atom %s" % atom.id_str(), file=log)
    deleted = []
    bad_atoms = [ pdb_atoms[j_seq] for j_seq in remove_atoms ]
    for atom in bad_atoms :
      if (atom.i_seq in deleted) : continue
      labels = atom.fetch_labels()
      atom_group = atom.parent()
      residue_group = atom_group.parent()
      other_atoms = residue_group.atoms()
      chain = residue_group.parent()
      if (labels.resname == "HOH"):
        deleted.extend(other_atoms.extract_i_seq())
        chain.remove_residue_group(residue_group)
        show_removed(other_atoms)
      else :
        if (chain.is_protein()):
          if (atom.name.strip() in ["N","C","CA","CB","O","H","HA"]):
            show_removed(other_atoms)
            deleted.extend(other_atoms.extract_i_seq())
            chain.remove_residue_group(residue_group)
          else :
            for atom2 in other_atoms :
              if (atom2.name.strip() in ["N","C","CA","CB","O","H","HA"]):
                continue
              else :
                show_removed([atom2])
                deleted.append(atom2.i_seq)
                atom_group.remove_atom(atom2)
        else :
          deleted.append(atom.i_seq)
          show_removed([atom])
          atom_group.remove_atom(atom)
          if (len(other_atoms) == 1):
            chain.remove_residue_group(residue_group)
    n_removed = len(deleted)
    print("%d atoms removed due to clashes with ligand(s)." % n_removed, file=log)
    for i_seq in deleted :
      selection[i_seq] = False
    xrs_new = xray_structure.select(selection)
  else :
    xrs_new = xray_structure
  return xrs_new

class get_final_maps_and_cc(object):
  def __init__(self,
      fmodel,
      ligands,
      params,
      log):
    from cctbx import maptbx
    from scitbx.array_family import flex
    map_helper = fmodel.electron_density_map()
    self.two_fofc_map_coeffs = map_helper.map_coefficients("2mFo-DFc")
    self.fofc_map_coeffs = map_helper.map_coefficients("mFo-DFc")
    fft_map = self.two_fofc_map_coeffs.fft_map(resolution_factor=0.25)
    fft_map.apply_sigma_scaling()
    fcalc = map_helper.map_coefficients("Fc")
    fcalc_map = fcalc.fft_map(resolution_factor=0.25)
    fcalc_map.apply_sigma_scaling()
    real_map = fft_map.real_map()
    fcalc_real_map = fcalc_map.real_map()
    final_cc = []
    for k, ligand in enumerate(ligands):
      atoms = ligand.atoms()
      sites = flex.vec3_double()
      radii = flex.double()
      for atom in atoms :
        if (not atom.element.strip() in ["H","D"]):
          sites.append(atom.xyz)
          radii.append(1.5)
      sel = maptbx.grid_indices_around_sites(
        unit_cell  = self.two_fofc_map_coeffs.unit_cell(),
        fft_n_real = real_map.focus(),
        fft_m_real = real_map.all(),
        sites_cart = sites,
        site_radii = radii)
      m1 = real_map.select(sel)
      m2 = fcalc_real_map.select(sel)
      cc = flex.linear_correlation(x=m1, y=m2).coefficient()
      final_cc.append(cc)
      print("  Ligand %d: CC = %5.3f" % (k+1, cc), file=log)
    print("", file=log)
    self.final_cc = final_cc

  def write_maps(self, file_name):
    import iotbx.mtz
    dec = iotbx.mtz.label_decorator(phases_prefix="PH")
    mtz_dat = self.two_fofc_map_coeffs.as_mtz_dataset(
      column_root_label="2FOFCWT",
      label_decorator=dec)
    mtz_dat.add_miller_array(self.fofc_map_coeffs,
      column_root_label="FOFCWT",
      label_decorator=dec)
    mtz_dat.mtz_object().write(file_name)

def extract_ligand_residues(
    pdb_hierarchy,
    ligand_code,
    atom_selection=None,
    only_segid=None):
  """
  Extract the atom_group object(s) with given 3-character residue name.

  :param pdb_hierarchy: input model
  :param ligand_code: 3-letter residue ID
  :param atom_selection: optional flex.bool object specifying ligand selection
      to use
  :param only_segid: optional segid which ligand(s) must match
  :returns: list of atom_group objects
  """
  assert (len(pdb_hierarchy.models()) == 1)
  ligands = []
  for chain in pdb_hierarchy.models()[0].chains():
    for residue_group in chain.residue_groups():
      for atom_group in residue_group.atom_groups():
        if (atom_group.resname == ligand_code):
          use_ligand = True
          if (only_segid is not None):
            for atom in atom_group.atoms():
              if (atom.segid != only_segid):
                use_ligand = False
                break
          if (use_ligand) and (atom_selection is not None):
            for atom in atom_group.atoms():
              if (not atom_selection[atom.i_seq]):
                use_ligand = False
                break
          if (use_ligand):
            ligands.append(atom_group)
  return ligands

def combine_ligands_and_hierarchy(pdb_hierarchy, ligands, log=None):
  from iotbx.pdb import hierarchy
  if (log is None) : log = null_out()
  chain_id_counts = {}
  model = pdb_hierarchy.models()[0]
  for i_lig, ligand in enumerate(ligands):
    xyz_mean = ligand.atoms().extract_xyz().mean()
    best_chain = None
    min_dist = sys.maxsize
    for chain in model.chains():
      last_resseq = chain.residue_groups()[-1].resseq_as_int()
      if ((not chain.id in chain_id_counts) or
          (chain_id_counts[chain.id] < last_resseq)):
        chain_id_counts[chain.id] = last_resseq
      if (not chain.is_protein()) : continue
      chain_xyz_mean = chain.atoms().extract_xyz().mean()
      dist = xyz_distance(chain_xyz_mean, xyz_mean)
      if (dist < min_dist):
        min_dist = dist
        best_chain = chain
    best_chain_id = " "
    if (best_chain is not None):
      best_chain_id = best_chain.id
    new_chain = hierarchy.chain(id=best_chain_id)
    new_rg = hierarchy.residue_group()
    new_resseq = 1
    if (best_chain_id in chain_id_counts):
      new_resseq = chain_id_counts[best_chain_id] + 1
    print("  ligand %d: chain='%s' resseq=%s" % (i_lig+1,
      best_chain_id, new_resseq), file=log)
    new_rg.resseq = new_resseq
    new_rg.append_atom_group(ligand)
    new_chain.append_residue_group(new_rg)
    model.append_chain(new_chain)
    chain_id_counts[best_chain_id] = new_resseq

# Main function
class apply_ligand_ncs(object):
  """
  Wrapper class; this should be the primary entry point for external calling
  routines.

  :param pdb_hierarchy: initial model with one or more copies of the target
    ligand
  :param fmodel: mmtbx.f_model.manager object corresponding to the model
  :param ligand_code: three-letter residue name of ligand
  :param params: phil scope_extract object for ncs_ligand_phil
  :param atom_selection: selection for reference ligand (default: search for
    all copies and pick one with the best CC)
  :param add_new_ligands_to_pdb: generate combined PDB hierarchy and fmodel
    with new ligands incorporated
  :param only_segid:
  :param log: filehandle-like object
  """
  def __init__(self,
      pdb_hierarchy,
      fmodel,
      ligand_code,
      params,
      atom_selection=None,
      add_new_ligands_to_pdb=False,
      only_segid=None,
      log=None):
    if (log is None) : log = sys.stdout
    assert (ligand_code is not None) and (len(ligand_code) <= 3)
    make_sub_header("Determining NCS operators", log)
    import iotbx.ncs
    ncs_obj = iotbx.ncs.input(hierarchy=pdb_hierarchy)
    nrgl = ncs_obj.get_ncs_restraints_group_list()

    if 0:
      print("Summary of NCS operators:", file=log)
      for ncs_group in ncs_ops :
        for k, group in enumerate(ncs_group):
          group.show_summary(log, prefix="  ")
    print("Looking for ligands named %s..." % ligand_code, file=log)
    ligands = extract_ligand_residues(pdb_hierarchy, ligand_code,
      only_segid=only_segid)
    if (len(ligands) == 0):
      raise Sorry("No ligands found!")
    pdb_hierarchy.atoms().reset_i_seq()
    make_sub_header("Applying NCS to ligands", log)
    sampler = sample_operators(
      fmodel=fmodel,
      pdb_hierarchy=pdb_hierarchy,
      ncs_operators=nrgl, #ncs_ops_flat,
      params=params,
      ligands=ligands,
      log=log)
    if (len(sampler.new_ligands) > 0):
      if (params.remove_clashing_atoms):
        make_sub_header("Removing clashing atoms", log)
        xrs_new = remove_clashing_atoms(
          xray_structure=sampler.xray_structure,
          pdb_hierarchy=pdb_hierarchy,
          ligands=sampler.new_ligands,
          params=params,
          log=log)
        fmodel.update_xray_structure(
          xray_structure=xrs_new,
          update_f_calc=True,
          update_f_mask=True)
      if (add_new_ligands_to_pdb):
        combine_ligands_and_hierarchy(
          pdb_hierarchy=pdb_hierarchy,
          ligands=sampler.new_ligands,
          log=log)
        xrs_new = pdb_hierarchy.extract_xray_structure(
          crystal_symmetry=fmodel.xray_structure)
        fmodel.update_xray_structure(
          xray_structure=xrs_new,
          update_f_calc=True,
          update_f_mask=True)
    self.final_ligands = extract_ligand_residues(pdb_hierarchy, ligand_code,
      only_segid=only_segid)
    self.final_cc = get_final_maps_and_cc(
      fmodel=fmodel,
      ligands=self.final_ligands,
      params=params,
      log=log)
    self.n_ligands_start = len(ligands)
    self.n_ligands_new = len(sampler.new_ligands)
    self.pdb_hierarchy = pdb_hierarchy
    self.fmodel = fmodel

  @property
  def n_ligands(self):
    return self.n_ligands_start + self.n_ligands_new

  def write_pdb(self, file_name):
    f = open(file_name, "w")
    f.write(self.pdb_hierarchy.as_pdb_string(
      crystal_symmetry=self.fmodel.xray_structure))
    f.close()

  def write_maps(self, file_name):
    self.final_cc.write_maps(file_name)
