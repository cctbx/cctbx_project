
"""
Post-fitting cleanup of ligand positions to match NCS operations present in
protein model.  This can be used to recover cases where one copy is placed
successfully but another misses due to interfering protein atoms, weak density,
false positive in empty protein density, etc.  Should not be used when the
initial placement (e.g. LigandFit CC) is sufficiently good.
"""

from libtbx.utils import Sorry
from libtbx.str_utils import make_header
from libtbx.test_utils import approx_equal
from libtbx import adopt_init_args, group_args
import copy
from math import sqrt
import sys

debug = True

ncs_ligand_phil = """
max_rmsd = 2.0
  .type = float
min_cc = 0.7
  .type = float
min_cc_reference = 0.9
  .type = float
max_cc_to_replace = 0.9
  .type = float
min_2fofc = 1.0
  .type = float
remove_clashing_atoms = True
  .type = bool
clash_cutoff = 2.0
  .type = float
write_sampled_pdbs = False
  .type = bool
"""

master_phil_str = """
include scope mmtbx.utils.cmdline_input_phil_str
ligand_code = LIG
  .type = str
output_file = None
  .type = path
output_map = None
  .type = path
%s
""" % ncs_ligand_phil

def master_phil () :
  import libtbx.phil
  return libtbx.phil.parse(master_phil_str, process_includes=True)

def resid_str (atom) :
  labels = atom.fetch_labels()
  return "%s %s" % (labels.resname, labels.resid())

class group_operators (object) :
  def __init__ (self, selection, sele_str, sites_cart) :
    self.selection = selection
    self.selection_string = sele_str
    self.center_of_mass = sites_cart.select(selection).mean()
    self.operators = []
    self.op_selections = []

  def add_operator (self, ops, sele_str) :
    self.operators.append(ops)
    self.op_selections.append(sele_str)

  def distance_from_center (self, sites) :
    c_o_m = sites.mean()
    x1,y1,z1 = c_o_m
    x2,y2,z2 = self.center_of_mass
    return sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)

  def show_summary (self, out=None) :
    if (out is None) : out = sys.stdout
    print >> out, "Reference selection: %s" % self.selection_string
    for op, op_sele in zip(self.operators, self.op_selections) :
      print >> out, "  Selection: %s" % op_sele
      print >> out, "  Rotation:"
      print >> out, "    %6.4f  %6.4f  %6.4f" % op.r.elems[0:3]
      print >> out, "    %6.4f  %6.4f  %6.4f" % op.r.elems[3:6]
      print >> out, "    %6.4f  %6.4f  %6.4f" % op.r.elems[6:9]
      print >> out, "  Translation:"
      print >> out, "    %6.4f  %6.4f  %6.4f" % op.t.elems
      print >> out, ""
    print >> out, ""

def find_ncs_operators (pdb_hierarchy, max_rmsd=2.0, log=None) :
  """
  Determines all possible NCS transformation matrices for the input structure,
  based on sequence alignemnt and simple C-alpha superposition.  There may be
  multiple sets of operators but these will eventually become a flat list.
  """
  from mmtbx.torsion_restraints import torsion_ncs
  from scitbx.math import superpose
  from scitbx.array_family import flex
  ncs_groups = torsion_ncs.determine_ncs_groups(pdb_hierarchy, log=log)
  print ncs_groups
  if (len(ncs_groups) == 0) :
    raise Sorry("No NCS present in the input model.")
  for k, group in enumerate(ncs_groups) :
    print >> log, "Group %d:" % (k+1)
    for sele in group :
      print >> log, "  %s" % sele
  selection_cache = pdb_hierarchy.atom_selection_cache()
  pdb_atoms = pdb_hierarchy.atoms()
  sites_cart = pdb_atoms.extract_xyz()
  operators = []
  def get_selection (sele_str) :
    sele_str = "(%s) and name CA" % sele_str
    return selection_cache.selection(sele_str).iselection()
  for restraint_group in ncs_groups :
    group_ops = []
    assert (len(restraint_group) >= 2)
    for j, sele_str in enumerate(restraint_group) :
      sele_j = get_selection(sele_str)
      group = group_operators(sele_j, sele_str, sites_cart)
      assert (len(sele_j) > 0)
      calpha_ids = []
      for i_seq in sele_j :
        calpha_ids.append(resid_str(pdb_atoms[i_seq]))
      for k, sele_str_k in enumerate(restraint_group) :
        if (k == j) : continue
        sele_k = get_selection(sele_str_k)
        group_sele = flex.size_t()
        group_ids = []
        assert (len(sele_k) > 0)
        # poor man's sequence alignment
        for i_seq in sele_k :
          id_str = resid_str(pdb_atoms[i_seq])
          group_ids.append(id_str)
          if (id_str in calpha_ids) :
            group_sele.append(i_seq)
        first_sele_copy = flex.size_t() #first_sele.deep_copy()
        delete_indices = []
        for i_seq, id_str in zip(sele_j, calpha_ids) :
          if (id_str in group_ids) :
            first_sele_copy.append(i_seq)
        assert (len(first_sele_copy) == len(group_sele))
        assert (len(group_sele) > 0)
        sites_ref = sites_cart.select(first_sele_copy)
        sites_group = sites_cart.select(group_sele).deep_copy()
        lsq_fit = superpose.sieve_fit(
          sites_fixed=sites_ref,
          sites_moving=sites_group,
          frac_discard=0.25)
        sites_fit = lsq_fit.r.elems * sites_group + lsq_fit.t.elems
        rmsd = sites_ref.rms_difference(sites_fit)
        print >> log, "  %d versus %d RMSD = %.3f" % (j+1, k+1, rmsd)
        if (rmsd <= max_rmsd) :
          group.add_operator(lsq_fit.rt().inverse(), sele_str_k)
        else :
          print >> log, "  exceeds cutoff, will not use this operator"
      group_ops.append(group)
    operators.append(group_ops)
  return operators

class sample_operators (object) :
  """
  Determines an appropriate "reference" ligand, and samples the density around
  sites transformed by each operator, applying the best operator to the
  other ligand copies if meeting cutoff criteria.
  """
  def __init__ (self,
      pdb_hierarchy,
      fmodel,
      ncs_operators,
      ligands,
      params,
      log=None) :
    if (log is None) : log = sys.stdout
    if (params is None) :
      params = master_phil().fetch().extract()
    adopt_init_args(self, locals())
    self.xray_structure = fmodel.xray_structure.deep_copy_scatterers()
    xrs_ncs = fmodel.xray_structure.deep_copy_scatterers()
    from iotbx.pdb import hierarchy
    self.setup_maps()
    best_cc = 0
    best_k = -1
    best_ligand = None
    print >> log, "Identifying reference ligand..."
    def show_map_stats (prefix, stats) :
      print >> log, "   %s: CC = %5.3f  mean = %6.2f" % (prefix, stats.cc,
        stats.map_mean)
    for k, ligand in enumerate(ligands) :
      atoms = ligand.atoms()
      start = self.get_sites_cc(atoms)
      show_map_stats("Ligand %d" % (k+1), start)
      if (start.cc > best_cc) and (start.cc > params.min_cc_reference) :
        best_ligand = ligand
        best_k = k
        best_cc = best_cc
    if (best_ligand is None) :
      raise Sorry("No ligand with acceptable CC (>%.2f) found." %
        params.min_cc_reference)
    print >> log, "Copy #%d was the best, using that as reference" % (best_k+1)
    print >> log, ""
    sites_ref = best_ligand.atoms().extract_xyz()
    min_dist = sys.maxint
    best_group = None
    for op_group in ncs_operators :
      distance = op_group.distance_from_center(sites_ref)
      if (distance < min_dist) :
        best_group = op_group
        min_dist = distance
    if (best_group is not None) :
      print >> log, "This appears to be bound to the selection \"%s\"" % \
        best_group.selection_string
    self.n_moved = 0
    self.n_removed = 0
    used_ops = []
    for k, ligand in enumerate(ligands) :
      if (k == best_k) :
        continue
      print >> log, "Sampling new positions for ligand %d..." % (k+1)
      atoms = ligand.atoms()
      start = self.get_sites_cc(atoms)
      show_map_stats("     start", start)
      if ((start.cc >= params.max_cc_to_replace) and
          (start.map_mean > params.min_2fofc)) :
        print >> log, "  CC is above cutoff for moving (%g), will skip" % \
          params.max_cc_to_replace
        continue
      cc_best = start.cc
      map_best = start.map_mean
      sites_best = None
      op_best = None
      sites = atoms.extract_xyz()
      for j, operator in enumerate(best_group.operators) :
        if (operator in used_ops) :
          continue
        sites_new = operator.r.elems * sites_ref + operator.t.elems
        stats_new = self.get_sites_cc(atoms, sites_new)
        show_map_stats("NCS op. %2d" % (j+1), stats_new)
        if (params.write_sampled_pdbs) :
          lig_new = ligand.detached_copy()
          lig_rg = hierarchy.residue_group()
          lig_rg.resseq = k+1
          lig_rg.append_atom_group(lig_new)
          lig_new.atoms().set_xyz(sites_new)
          f = open("ncs_ligand_%d_%d.pdb" % (k+1, j+1), "w")
          for atom in lig_new.atoms() :
            f.write(atom.format_atom_record()+"\n")
          f.close()
        if ((stats_new.cc > params.min_cc) and
            (stats_new.cc > cc_best) and
            (stats_new.map_mean > map_best)) : # arams.min_2fofc) and
          cc_best = stats_new.cc
          sites_best = sites_new
          map_best = stats_new.map_mean
          op_best = operator
      if (op_best is not None) :
        used_ops.append(op_best)
      if (sites_best is None) :
        print >> log, "  no acceptable sites found."
      elif ((cc_best > start.cc) or (map_best > start.map_mean)) :
        print >> log, "  applying best sites"
        ligand.atoms().set_xyz(sites_best)
        all_sites = xrs_ncs.sites_cart()
        all_sites.set_selected(ligand.atoms().extract_i_seq(), sites_best)
        xrs_ncs.set_sites_cart(all_sites)
        self.n_moved += 1
    self.xray_structure = xrs_ncs
    if (debug) :
      sites_cart = xrs_ncs.sites_cart()
      for i_seq, atom in enumerate(pdb_hierarchy.atoms()) :
        assert approx_equal(atom.xyz, sites_cart[i_seq], eps=0.001)
    self.fmodel.update_xray_structure(
      xray_structure=xrs_ncs,
      update_f_calc=True,
      update_f_mask=True)
    if (self.n_moved == 0) :
      print >> log, "No NCS operators applied."

  def setup_maps (self) :
    map_helper = self.fmodel.electron_density_map()
    map_coeffs = map_helper.map_coefficients("2mFo-DFc")
    diff_map_coeffs = map_helper.map_coefficients("mFo-DFc")
    fft_map = map_coeffs.fft_map(resolution_factor=0.25)
    diff_fft_map = diff_map_coeffs.fft_map(resolution_factor=0.25)
    fft_map.apply_sigma_scaling()
    diff_fft_map.apply_sigma_scaling()
    fcalc = map_helper.map_coefficients("Fc")
    fcalc_map = fcalc.fft_map(resolution_factor=0.25)
    fcalc_map.apply_sigma_scaling()
    self.unit_cell = map_coeffs.unit_cell()
    self.real_map = fft_map.real_map()
    self.diff_map = diff_fft_map.real_map()
    self.n_real = self.real_map.focus()
    self.m_real = self.real_map.all()
    self.fcalc_real_map = fcalc_map.real_map()

  def get_new_fcalc_map (self, sites_new, i_seqs) :
    xrs_new = self.xray_structure.deep_copy_scatterers()
    all_sites = xrs_new.sites_cart()
    all_sites.set_selected(i_seqs, sites_new)
    xrs_new.set_sites_cart(all_sites)
    self.fmodel.update_xray_structure(
      xray_structure=xrs_new,
      update_f_calc=True,
      update_f_mask=True)
    fcalc = self.fmodel.electron_density_map().map_coefficients("Fc")
    fcalc_map = fcalc.fft_map(resolution_factor=0.25)
    fcalc_map.apply_sigma_scaling()
    # XXX now revert to original xray structure
    self.fmodel.update_xray_structure(
      xray_structure=self.xray_structure,
      update_f_calc=True,
      update_f_mask=True)
    return fcalc_map.real_map()

  def get_sites_cc (self, atoms, sites=None) :
    from cctbx import maptbx
    from scitbx.array_family import flex
    radii = flex.double()
    for atom in atoms :
      if (atom.element.strip() in ["H", "D"]) :
        radii.append(1.)
      else :
        radii.append(1.5)
    fcalc_map = self.fcalc_real_map
    if (sites is None) :
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

def is_same_residue (atom1, atom2) :
  labels1 = atom1.fetch_labels()
  labels2 = atom2.fetch_labels()
  return ((labels1.resid() == labels2.resid()) and
          (labels1.chain_id == labels2.chain_id) and
          (labels1.altloc == labels2.altloc))

def remove_clashing_atoms (
    xray_structure,
    pdb_hierarchy,
    ligands,
    params,
    log) :
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
      for j_seq, sym_ops in asu_dict.items() :
        if (not is_same_residue(pdb_atoms[i_seq], pdb_atoms[j_seq])) :
          remove_atoms.add(j_seq)
  if (len(remove_atoms) > 0) :
    def show_removed (atoms) :
      for atom in atoms :
        print >> log, "  warning: deleting atom %s" % atom.id_str()
    deleted = []
    bad_atoms = [ pdb_atoms[j_seq] for j_seq in remove_atoms ]
    for atom in bad_atoms :
      if (atom.i_seq in deleted) : continue
      labels = atom.fetch_labels()
      atom_group = atom.parent()
      residue_group = atom_group.parent()
      other_atoms = residue_group.atoms()
      chain = residue_group.parent()
      if (labels.resname == "HOH") :
        deleted.extend(other_atoms.extract_i_seq())
        chain.remove_residue_group(residue_group)
        show_removed(other_atoms)
      else :
        if (chain.conformers()[0].is_protein()) :
          if (atom.name.strip() in ["N","C","CA","CB","O","H","HA"]) :
            show_removed(other_atoms)
            deleted.extend(other_atoms.extract_i_seq())
            chain.remove_residue_group(residue_group)
          else :
            for atom2 in other_atoms :
              if (atom2.name.strip() in ["N","C","CA","CB","O","H","HA"]) :
                continue
              else :
                show_removed([atom2])
                deleted.append(atom2.i_seq)
                atom_group.remove_atom(atom2)
        else :
          deleted.append(atom.i_seq)
          show_removed([atom])
          atom_group.remove_atom(atom)
          if (len(other_atoms) == 1) :
            chain.remove_residue_group(residue_group)
    n_removed = len(deleted)
    print >> log, "%d atoms removed due to clashes with ligand(s)." % n_removed
    for i_seq in deleted :
      selection[i_seq] = False
    xrs_new = xray_structure.select(selection)
  else :
    xrs_new = xray_structure
  return xrs_new

def get_final_maps_and_cc (
    fmodel,
    ligands,
    params,
    log) :
  from cctbx import maptbx
  from scitbx.array_family import flex
  map_helper = fmodel.electron_density_map()
  map_coeffs = map_helper.map_coefficients("2mFo-DFc")
  fft_map = map_coeffs.fft_map(resolution_factor=0.25)
  fft_map.apply_sigma_scaling()
  fcalc = map_helper.map_coefficients("Fc")
  fcalc_map = fcalc.fft_map(resolution_factor=0.25)
  fcalc_map.apply_sigma_scaling()
  real_map = fft_map.real_map()
  fcalc_real_map = fcalc_map.real_map()
  final_cc = []
  for k, ligand in enumerate(ligands) :
    atoms = ligand.atoms()
    sites = atoms.extract_xyz()
    radii = flex.double()
    for atom in atoms :
      if (atom.element.strip() in ["H", "D"]) :
        radii.append(1.)
      else :
        radii.append(1.5)
    sel = maptbx.grid_indices_around_sites(
      unit_cell  = map_coeffs.unit_cell(),
      fft_n_real = real_map.focus(),
      fft_m_real = real_map.all(),
      sites_cart = sites,
      site_radii = radii)
    m1 = real_map.select(sel)
    m2 = fcalc_real_map.select(sel)
    cc = flex.linear_correlation(x=m1, y=m2).coefficient()
    final_cc.append(cc)
    print >> log, "  Ligand %d: CC = %5.3f" % (k+1, cc)
  print >> log, ""
  if (params is not None) and (params.output_map is not None) :
    import iotbx.mtz
    diff_map = map_helper.map_coefficients("mFo-DFc")
    dec = iotbx.mtz.label_decorator(phases_prefix="PH")
    mtz_dat = map_coeffs.as_mtz_dataset(
      column_root_label="2FOFCWT",
      label_decorator=dec)
    mtz_dat.add_miller_array(diff_map,
      column_root_label="FOFCWT",
      label_decorator=dec)
    mtz_dat.mtz_object().write(output_map)
    print >> log, "Wrote %s" % params.output_map
  return final_cc

def extract_ligand_residues (pdb_hierarchy, ligand_code) :
  assert (len(pdb_hierarchy.models()) == 1)
  ligands = []
  for chain in pdb_hierarchy.models()[0].chains() :
    for residue_group in chain.residue_groups() :
      for atom_group in residue_group.atom_groups() :
        if (atom_group.resname == ligand_code) :
          ligands.append(atom_group)
  return ligands

# Main function
def apply_ligand_ncs (
    pdb_hierarchy,
    fmodel,
    params=None,
    log=None) :
  if (log is None) : log = sys.stdout
  if (params is None) :
    params = master_phil().fetch().extract()
  make_header("Determining NCS operators", log)
  ncs_ops = find_ncs_operators(pdb_hierarchy,
    max_rmsd=params.max_rmsd,
    log=log)
  ncs_ops_flat = []
  for op_group in ncs_ops :
    ncs_ops_flat.extend(op_group)
  print "Summary of NCS operators:"
  for ncs_group in ncs_ops :
    for k, group in enumerate(ncs_group) :
      group.show_summary(log)
  ligands = extract_ligand_residues(pdb_hierarchy, params.ligand_code)
  if (len(ligands) == 0) :
    raise Sorry("No ligands found!")
  pdb_hierarchy.atoms().reset_i_seq()
  make_header("Applying NCS to ligands", log)
  sampler = sample_operators(
    fmodel=fmodel,
    pdb_hierarchy=pdb_hierarchy,
    ncs_operators=ncs_ops_flat,
    params=params,
    ligands=ligands,
    log=log)
  if (sampler.n_moved > 0) and (params.remove_clashing_atoms) :
    make_header("Removing clashing atoms", log)
    xrs_new = remove_clashing_atoms(
      xray_structure=sampler.xray_structure,
      pdb_hierarchy=pdb_hierarchy,
      ligands=ligands,
      params=params,
      log=log)
    fmodel.update_xray_structure(
      xray_structure=xrs_new,
      update_f_calc=True,
      update_f_mask=True)
  make_header("Statistics for final ligands", log)
  final_cc = get_final_maps_and_cc(
    fmodel=fmodel,
    ligands=ligands,
    params=params,
    log=log)
  if (params.output_file is not None) :
    f = open(params.output_file, "w")
    cryst1 = iotbx.pdb.format_cryst1_record(fmodel.xray_structure)
    f.write("%s\n" % cryst1)
    f.write(pdb_hierarchy.as_pdb_string())
    f.close()
    print >> out, ""
    print >> out, "Wrote %s" % params.output_file
  return final_cc

def run (args, out=None) :
  if (out is None) : out = sys.stdout
  from mmtbx.utils import cmdline_load_pdb_and_data
  import iotbx.pdb
  cmdline = cmdline_load_pdb_and_data(
    args=args,
    master_phil=master_phil(),
    out=out,
    process_pdb_file=False)
  pdb_hierarchy = cmdline.pdb_hierarchy
  fmodel = cmdline.fmodel
  params = cmdline.params
  if (params.output_file is None) :
    params.output_file = "ncs_ligands.pdb"
  if (params.output_map is None) :
    params.output_map = "ncs_ligands.mtz"
  final_cc = apply_ligand_ncs(
    pdb_hierarchy=pdb_hierarchy,
    fmodel=fmodel,
    params=params,
    log=out)
  return final_cc

if (__name__ == "__main__") :
  run(sys.argv[1:])
