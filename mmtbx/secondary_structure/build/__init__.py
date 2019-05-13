from __future__ import division
from __future__ import print_function
from scitbx.math import superpose

from libtbx.utils import Sorry, null_out
from libtbx import Auto
from cctbx.array_family import flex
from cctbx import crystal
import iotbx.pdb
from iotbx.pdb.amino_acid_codes import one_letter_given_three_letter as one_three
from iotbx.pdb.amino_acid_codes import three_letter_given_one_letter as three_one
from mmtbx import secondary_structure
from mmtbx.monomer_library import idealized_aa
from mmtbx.refinement.real_space.individual_sites import minimize_wrapper_with_map
from mmtbx.refinement.geometry_minimization import minimize_wrapper_for_ramachandran
from mmtbx.secondary_structure.ss_validation import gather_ss_stats
from mmtbx.rotamer.rotamer_eval import RotamerEval
from time import time
import scitbx.math
import mmtbx.idealized_aa_residues.rotamer_manager
import mmtbx.refinement.real_space.fit_residues

alpha_helix_str = """
ATOM      1  N   GLY A   1      -5.606  -2.251 -12.878  1.00  0.00           N
ATOM      2  CA  GLY A   1      -5.850  -1.194 -13.852  1.00  0.00           C
ATOM      3  C   GLY A   1      -5.186  -1.524 -15.184  1.00  0.00           C
ATOM      4  O   GLY A   1      -5.744  -1.260 -16.249  1.00  0.00           O
ATOM      6  N   GLY A   2      -3.991  -2.102 -15.115  1.00  0.00           N
ATOM      7  CA  GLY A   2      -3.262  -2.499 -16.313  1.00  0.00           C
ATOM      8  C   GLY A   2      -3.961  -3.660 -17.011  1.00  0.00           C
ATOM      9  O   GLY A   2      -4.016  -3.716 -18.240  1.00  0.00           O
"""

a310_helix_str = """\
ATOM      1  N   GLY A   1       8.836  -4.233 -14.408  1.00  0.00           N
ATOM      2  CA  GLY A   1      10.232  -4.071 -14.799  1.00  0.00           C
ATOM      3  C   GLY A   1      10.764  -5.331 -15.476  1.00  0.00           C
ATOM      4  O   GLY A   1      11.679  -5.262 -16.297  1.00  0.00           O
ATOM      6  N   GLY A   2      10.176  -6.478 -15.143  1.00  0.00           N
ATOM      7  CA  GLY A   2      10.582  -7.741 -15.750  1.00  0.00           C
ATOM      8  C   GLY A   2      10.381  -7.714 -17.262  1.00  0.00           C
ATOM      9  O   GLY A   2      11.080  -8.410 -17.999  1.00  0.00           O
"""

pi_helix_str = """\
ATOM      1  N   GLY A   1      -3.365  -3.446  -8.396  1.00  0.00           N
ATOM      2  CA  GLY A   1      -4.568  -4.249  -8.592  1.00  0.00           C
ATOM      3  C   GLY A   1      -5.809  -3.386  -8.805  1.00  0.00           C
ATOM      4  O   GLY A   1      -6.559  -3.591  -9.759  1.00  0.00           O
ATOM      6  N   GLY A   2      -6.025  -2.424  -7.914  1.00  0.00           N
ATOM      7  CA  GLY A   2      -7.221  -1.588  -7.976  1.00  0.00           C
ATOM      8  C   GLY A   2      -7.101  -0.486  -9.025  1.00  0.00           C
ATOM      9  O   GLY A   2      -8.089  -0.114  -9.659  1.00  0.00           O
"""

beta_pdb_str = """
ATOM      1  N   GLY A   1      27.961   0.504   1.988  1.00  0.00           N
ATOM      2  CA  GLY A   1      29.153   0.205   2.773  1.00  0.00           C
ATOM      3  C   GLY A   1      30.420   0.562   2.003  1.00  0.00           C
ATOM      4  O   GLY A   1      30.753  -0.077   1.005  1.00  0.00           O
ATOM      6  N   GLY A   2      31.123   1.587   2.474  1.00  0.00           N
ATOM      7  CA  GLY A   2      32.355   2.031   1.832  1.00  0.00           C
ATOM      8  C   GLY A   2      33.552   1.851   2.758  1.00  0.00           C
ATOM      9  O   GLY A   2      33.675   2.539   3.772  1.00  0.00           O
"""

helix_class_to_pdb_str = {'alpha':alpha_helix_str,
                          'pi':pi_helix_str,
                          '3_10': a310_helix_str}

ss_idealization_master_phil_str = """
ss_idealization
{
  enabled = False
    .type = bool
    .help = Enable secondary structure idealization
  file_name_before_regularization = None
    .type = path
  skip_empty_ss_elements = True
    .type = bool
    .help = skip SS element if there is no atoms in the structure
  skip_good_ss_elements = False
    .type = bool
    .help = Skip SS elements without wrong/outlier Ramachandran angles and \
      with reasonable hydrogen bonds length
  restrain_torsion_angles = False
    .type = bool
    .help = Restrain torsion angles
  fix_rotamer_outliers = True
    .type = bool
    .help = Fix rotamer outliers after replacing SS elements
  sigma_on_reference_non_ss = 1
    .type = float
    .help = Weight on original model coordinates restraints where no ss \
      present. Keeps loops close to initial model. \
      (bigger number gives lighter restraints)
  sigma_on_reference_helix = 1
    .type = float
    .help = Weight on original model coordinates restraints where helix \
      present. Bends helices a bit according to initial model. \
      (bigger number gives lighter restraints)
  sigma_on_reference_sheet = 0.5
    .type = float
    .help = Weight on original model coordinates restraints where sheet \
      present. Bends helices a bit according to initial model. \
      (bigger number gives lighter restraints)
  sigma_on_torsion_ss = 5
    .type = float
    .help = Weight on torsion angles restraints where ss present. \
      Keeps helices torsion angles close to ideal. \
      (bigger number gives lighter restraints)
  sigma_on_torsion_nonss = 5
    .type = float
  sigma_on_ramachandran = 1
    .type = float
  sigma_on_cbeta = 2.5
    .type = float
  n_macro = Auto
    .type = int
  n_iter=300
    .type = int
}
"""

master_phil = iotbx.phil.parse(ss_idealization_master_phil_str)

def print_hbond_proxies(geometry, hierarchy, pymol=False):
  """ Print hydrogen bonds in geometry restraints manager for debugging
  purposes"""
  assert 0, "need to rewrite due to reorganization of GRM"
  atoms = hierarchy.atoms()
  if pymol:
    dashes = open('dashes.pml', 'w')
  hbondlen=flex.double()
  for hb in geometry.generic_restraints_manager.hbonds_as_simple_bonds():
    hbondlen.append(atoms[hb[0]].distance(atoms[hb[1]]))
    print((atoms[hb[0]].id_str(), "<====>",atoms[hb[1]].id_str(),
        atoms[hb[0]].distance(atoms[hb[1]]), hb[0], hb[1]))
    if pymol:
      s1 = atoms[hb[0]].id_str()
      s2 = atoms[hb[1]].id_str()
      #print "pdbstr1:", s1
      #print "pdbstr1:",s2
      ps = "dist chain \"%s\" and resi %s and name %s, chain \"%s\" and resi %s and name %s\n" % (s1[14:15],
         s1[16:19], s1[5:7], s2[14:15], s2[16:19], s2[5:7])
      dashes.write(ps)
  print("min, max, mean, sd hbond lenghts", hbondlen.min_max_mean().as_tuple(),\
    hbondlen.standard_deviation_of_the_sample())
  if pymol:
    dashes.close()

def get_r_t_matrices_from_structure(pdb_str):
  """ Return rotation and translation matrices for the ideal structure.

  The function determines r and t matrices from alingment of 1st and 2nd
  residues of the structure passed in pdb_str.
  """
  pdb_hierarchy = iotbx.pdb.input(source_info=None, lines=pdb_str).\
    construct_hierarchy()
  conformer = pdb_hierarchy.models()[0].chains()[0].conformers()[0]
  residues = conformer.residues()
  fixed_sites = flex.vec3_double()
  moving_sites = flex.vec3_double()
  main_chain_atoms = ["N","CA","C","O"]
  if len(residues)>=2:
    for (r, arr) in [(residues[0], fixed_sites), (residues[1], moving_sites)]:
      for a in r.atoms():
        if a.name.strip() in main_chain_atoms:
          arr.append(a.xyz)
  else:
    raise Sorry('pdb_str should contain at least 2 residues')
  lsq_fit_obj = superpose.least_squares_fit(reference_sites = moving_sites,
                                            other_sites = fixed_sites)
  return lsq_fit_obj.r, lsq_fit_obj.t


def side_chain_placement(ag_to_place, current_reference_ag, rotamer_manager):
  """
  Works with poly_gly truncated hierarchy.
  Also used in fix_rama_outliers.
  """
  resname = current_reference_ag.resname.upper()
  c = one_three.get(resname, None)

  # seems to work with unusual residues...
  # if c is None:
  #   msg = "Only standard protein residues are currently supported.\n"
  #   msg += "The residue %s (chain %s, resid %s) chain is not standard." % (
  #       resname,
  #       current_reference_ag.parent().parent().id,
  #       current_reference_ag.parent().resid())
  #   raise Sorry(msg)
  ag_to_place.resname = three_one.get(c,resname)
  if c == 'G':
    return

  # align residue from ideal_res_dict to just placed ALA (ag_to_place)
  # or from pdb_hierarchy_template
  fixed_sites = flex.vec3_double()
  moving_sites = flex.vec3_double()
  reper_atoms = ["C","CA", "N"]
  for (ag, arr) in [(ag_to_place, fixed_sites),
                    (current_reference_ag, moving_sites)]:
    for a in ag.atoms():
      if a.name.strip() in reper_atoms:
        arr.append(a.xyz)
  assert len(fixed_sites) == 3, ag_to_place.id_str()
  if len(moving_sites) < 3:
    error_msg = "C, CA or N atoms are absent in secondary structure element." +\
        "\nPlease add them to the model and try again."
    raise Sorry(error_msg)
  assert len(moving_sites) == 3
  lsq_fit_obj = superpose.least_squares_fit(reference_sites = fixed_sites,
                                            other_sites = moving_sites)
  ideal_correct_ag = current_reference_ag.detached_copy()
  ideal_correct_ag.atoms().set_xyz(
      lsq_fit_obj.r.elems*ideal_correct_ag.atoms().extract_xyz()+\
      lsq_fit_obj.t.elems)
  ideal_correct_ag.atoms().set_xyz(
      rotamer_manager.nearest_rotamer_sites_cart(ideal_correct_ag))
  if len(ideal_correct_ag.atoms()) > 4:
    ag_to_place.pre_allocate_atoms(number_of_additional_atoms=\
                                                len(ideal_correct_ag.atoms())-4)
    for a in ideal_correct_ag.atoms():
      if a.name.strip() not in ["N","CA","C","O"]:
        at = a.detached_copy()
        at.uij_erase()
        ag_to_place.append_atom(atom=at)
  else:
    # This means something wrong with input model, e.g. only 3 atoms in
    # the residue and they happened to be N, CA, C
    pass


def secondary_structure_from_sequence(pdb_str,
      sequence=None,
      pdb_hierarchy_template=None,
      rotamer_manager=None):
  """ Return pdb.hierarchy with secondary structure according to sequence or
  reference hierarcy. If reference hierarchy provided, the resulting hierarchy
  will be rigid body aligned to it. Residue numbers will start from 1.

  pdb_str - "ideal" structure at least 2 residues long.
  sequence - string with sequence (one-letter codes)
  pdb_hierarchy_template - reference hierarchy.
  """
  if rotamer_manager is None:
    rotamer_manager = RotamerEval()
  pht = pdb_hierarchy_template
  assert [sequence, pht].count(None) == 1
  if pht is not None:
    lk = len(pht.altloc_indices().keys())
    if lk ==0:
      raise Sorry(
          "Hierarchy template in secondary_structure_from_sequence is empty")
    else:
      if len(pht.altloc_indices().keys()) != 1:
        raise Sorry("Alternative conformations are not supported")
  number_of_residues = len(sequence) if sequence!=None else \
    len(pht.models()[0].chains()[0].conformers()[0].residues())
  if number_of_residues<1:
    raise Sorry('sequence should contain at least one residue.')
  ideal_res_dict = idealized_aa.residue_dict()
  real_res_list = None
  if pht:
    real_res_list = pht.models()[0].chains()[0].residue_groups()
  pdb_hierarchy = iotbx.pdb.input(source_info=None, lines=pdb_str).\
      construct_hierarchy()
  pdb_hierarchy.truncate_to_poly_gly()
  chain = pdb_hierarchy.models()[0].chains()[0]
  current_gly_ag = chain.residue_groups()[0].atom_groups()[0]
  new_chain = iotbx.pdb.hierarchy.chain(id="A")
  new_chain.pre_allocate_residue_groups(number_of_additional_residue_groups=\
                                                            number_of_residues)
  r, t = get_r_t_matrices_from_structure(pdb_str)
  for j in range(number_of_residues):
    # put ALA
    rg = iotbx.pdb.hierarchy.residue_group(icode="")
    rg.resseq = j+1
    new_chain.append_residue_group(residue_group=rg)
    ag_to_place = current_gly_ag.detached_copy()
    rg.append_atom_group(atom_group=ag_to_place)
    current_gly_ag.atoms().set_xyz(
                          r.elems*current_gly_ag.atoms().extract_xyz()+t.elems)
    current_reference_ag = real_res_list[j].atom_groups()[0] if pht else \
        ideal_res_dict[three_one[sequence[j]].lower()].models()[0].chains()[0].\
        residue_groups()[0].atom_groups()[0]
    side_chain_placement(ag_to_place, current_reference_ag, rotamer_manager)
  new_pdb_h = iotbx.pdb.hierarchy.new_hierarchy_from_chain(new_chain)
  # align to real
  if pht != None:
    fixed_sites, moving_sites = get_matching_sites_cart_in_both_h(pht, new_pdb_h)
    assert len(fixed_sites) == len(moving_sites)
    lsq_fit_obj = superpose.least_squares_fit(reference_sites = fixed_sites,
                                              other_sites = moving_sites)
    new_pdb_h.atoms().set_xyz(
        lsq_fit_obj.r.elems*new_pdb_h.atoms().extract_xyz()+lsq_fit_obj.t.elems)
  return new_pdb_h

def get_helix(helix_class, rotamer_manager, sequence=None, pdb_hierarchy_template=None):
  if helix_class not in helix_class_to_pdb_str.keys():
    raise Sorry("Unsupported helix type.")
  return secondary_structure_from_sequence(
    pdb_str=helix_class_to_pdb_str[helix_class],
    sequence=sequence,
    rotamer_manager=rotamer_manager,
    pdb_hierarchy_template=pdb_hierarchy_template)

def calculate_rmsd_smart(h1, h2, backbone_only=False):
  # assert h1.atoms_size() == h2.atoms_size(), "%d!=%d" % (h1.atoms_size(),h2.atoms_size())
  rmsd = 0
  n = 0
  for atom in h1.atoms():
    if backbone_only and atom.name.strip() not in ['N', 'CA', 'C', 'O']:
      continue
    for c in h2.chains():
      if c.id != atom.parent().parent().parent().id:
        continue
      for rg in c.residue_groups():
        if (rg.resseq, rg.icode) != (atom.parent().parent().resseq, atom.parent().parent().icode):
          continue
        for ag in rg.atom_groups():
          if (ag.resname, ag.altloc) != (atom.parent().resname, atom.parent().altloc):
            continue
          a = ag.get_atom(atom.name.strip())
          if a is not None:
            rmsd += a.distance(atom)**2
            n += 1
  if n == 0:
    return rmsd ** 0.5
  return (rmsd/float(n)) ** 0.5

def set_xyz_smart(dest_h, source_h):
  """
  Even more careful setting of coordinates than set_xyz_carefully below
  Not clear why assert dest_h.atoms_size() =< source_h.atoms_size()
  """
  # try shortcut
  # print "SHORTCUT atoms:", dest_h.atoms_size(), source_h.atoms_size()
  # dest_h.write_pdb_file(file_name="dest_h.pdb")
  # source_h.write_pdb_file(file_name="source_h.pdb")
  if dest_h.atoms_size() == source_h.atoms_size():
    # print "TRYING SHORTCUT"
    good = True
    for a1, a2 in zip(dest_h.atoms(), source_h.atoms()):
      if a1.id_str() != a2.id_str():
        print(a1.id_str(), a2.id_str())
        good = False
        break
    if good:
      dest_h.atoms().set_xyz(source_h.atoms().extract_xyz())
      # print "SHORTCUT DONE"
      return
  # if dest_h.atoms_size() > source_h.atoms_size():
  #   print "Atom sizes:", dest_h.atoms_size(), source_h.atoms_size()
  #   dest_h.write_pdb_file("dest_h.pdb")
  #   source_h.write_pdb_file("source_h.pdb")
  #   assert 0
  for atom in source_h.atoms():
    chains = []
    if hasattr(dest_h, 'chains'):
      chains = dest_h.chains()
    else:
      chains = [dest_h]
    for c in chains:
      if c.id != atom.parent().parent().parent().id:
        continue
      for rg in c.residue_groups():
        if (rg.resseq, rg.icode) != (atom.parent().parent().resseq, atom.parent().parent().icode):
          continue
        for ag in rg.atom_groups():
          if (ag.resname, ag.altloc) != (atom.parent().resname, atom.parent().altloc):
            continue
          # print "atom name", atom.name
          a = ag.get_atom(atom.name.strip())
          if a is not None:
            # print "actually setting coordinates:", a.xyz, "->", atom.xyz
            a.set_xyz(atom.xyz)

def set_xyz_carefully(dest_h, source_h):
  # This assertion fails when in dest_h some main-chain atoms were missing
  # Nevertheless, the function should work anyway.
  # assert dest_h.atoms_size() >= source_h.atoms_size()
  for d_ag, s_ag in zip(dest_h.atom_groups(), source_h.atom_groups()):
    for s_atom in s_ag.atoms():
      d_atom = d_ag.get_atom(s_atom.name.strip())
      if d_atom is not None:
        d_atom.set_xyz(s_atom.xyz)

def get_matching_sites_cart_in_both_h(old_h, new_h):
  old_h.reset_atom_i_seqs()
  new_h.reset_atom_i_seqs()
  fixed_sites = flex.vec3_double()
  moving_sites = flex.vec3_double()
  isel_for_old = flex.size_t()
  isel_for_new = flex.size_t()
  if old_h.atoms_size() == new_h.atoms_size():
    good = True
    for a1, a2 in zip(old_h.atoms(), new_h.atoms()):
      if a1.id_str()[:-6] != a2.id_str()[:-6]:
        # print "No match: '%s', '%s'" % (a1.id_str()[:-6], a2.id_str()[:-6])
        good = False
        break
      else:
        fixed_sites.append(a1.xyz)
        moving_sites.append(a2.xyz)
    if good:
      # print "SHORTCUT DONE"
      assert fixed_sites.size() == moving_sites.size()
      return fixed_sites, moving_sites
  fixed_sites = flex.vec3_double()
  moving_sites = flex.vec3_double()
  for old_rg, new_rg in zip(old_h.only_chain().residue_groups(), new_h.only_chain().residue_groups()):
    for old_ag, new_ag in zip(old_rg.atom_groups(), new_rg.atom_groups()):
      for atom in old_ag.atoms():
        a = new_ag.get_atom(atom.name.strip())
        if a is not None:
          fixed_sites.append(atom.xyz)
          moving_sites.append(a.xyz)
  assert fixed_sites.size() == moving_sites.size()
  return fixed_sites, moving_sites

def process_params(params):
  min_sigma = 1e-5
  if params is None:
    params = master_phil.fetch().extract()
    params.ss_idealization.enabled = True
  if hasattr(params, "ss_idealization"):
    p_pars = params.ss_idealization
  else:
    assert hasattr(params, "enabled") and hasattr(params, "sigma_on_cbeta"), \
        "Something wrong with parameters passed to ss_idealization"
    p_pars = params
  assert isinstance(p_pars.enabled, bool)
  assert isinstance(p_pars.restrain_torsion_angles, bool)
  for par in ["sigma_on_reference_non_ss",
      "sigma_on_reference_helix", "sigma_on_reference_sheet",
      "sigma_on_torsion_ss", "sigma_on_torsion_nonss", "sigma_on_ramachandran",
      "sigma_on_cbeta"]:
    assert (isinstance(getattr(p_pars, par), float) and \
      getattr(p_pars, par) > min_sigma), "" + \
      "Bad %s parameter" % par
  for par in ["n_iter"]:
    assert (isinstance(getattr(p_pars, par), int) and \
      getattr(p_pars, par) >= 0), ""+ \
      "Bad %s parameter" % par
  return p_pars

def ss_element_is_good(ss_stats_obj, hsh_tuple):
  (n_hbonds, n_bad_hbonds, n_mediocre_hbonds, hb_lens, n_outliers,
      n_wrong_region) = ss_stats_obj(hsh_tuple)
  # check Rama angles. Any outlier is bad.
  if n_wrong_region+n_outliers > 0:
    return False
  # check hbonds. Theoretically, they could be corrected by refinement.
  if n_bad_hbonds > n_hbonds/2.:
    return False
  return True

def substitute_ss(
                    model, # changed in place
                    params = None,
                    log=null_out(),
                    reference_map=None,
                    verbose=False):
  """
  Substitute secondary structure elements in real_h hierarchy with ideal
  ones _in_place_.
  Returns reference torsion proxies - the only thing that cannot be restored
  with little effort outside the procedure.
  """

  ss_annotation = model.get_ss_annotation()

  t0 = time()
  if model.get_hierarchy().models_size() > 1:
    raise Sorry("Multi model files are not supported")
  for m in model.get_hierarchy().models():
    for chain in m.chains():
      if len(chain.conformers()) > 1:
        raise Sorry("Alternative conformations are not supported.")

  processed_params = process_params(params)
  if not processed_params.enabled:
    return None
  if ss_annotation is None:
    return None

  ann = ss_annotation
  if model.ncs_constraints_present():
    print("Using master NCS to reduce amount of work", file=log)

  expected_n_hbonds = 0
  for h in ann.helices:
    expected_n_hbonds += h.get_n_maximum_hbonds()
  edited_h = model.get_hierarchy().deep_copy()
  n_atoms_in_real_h = model.get_number_of_atoms()
  selection_cache = model.get_atom_selection_cache()

  # check the annotation for correctness (atoms are actually in hierarchy)
  error_msg = "The following secondary structure annotations result in \n"
  error_msg +="empty atom selections. They don't match the structre: \n"
  t1 = time()
  # Checking for SS selections
  deleted_annotations = ann.remove_empty_annotations(
      hierarchy=model.get_hierarchy(),
      asc=selection_cache)
  if not deleted_annotations.is_empty():
    if processed_params.skip_empty_ss_elements:
      if len(deleted_annotations.helices) > 0:
        print("Removing the following helices because there are", file=log)
        print("no corresponding atoms in the model:", file=log)
        for h in deleted_annotations.helices:
          print(h.as_pdb_str(), file=log)
          error_msg += "  %s\n" % h
      if len(deleted_annotations.sheets) > 0:
        print("Removing the following sheets because there are", file=log)
        print("no corresponding atoms in the model:", file=log)
        for sh in deleted_annotations.sheets:
          print(sh.as_pdb_str(), file=log)
          error_msg += "  %s\n" % sh.as_pdb_str(strand_id=st.strand_id)
    else:
      raise Sorry(error_msg)
  phil_str = ann.as_restraint_groups()

  # gathering initial special position atoms
  special_position_settings = crystal.special_position_settings(
      crystal_symmetry = model.crystal_symmetry())
  site_symmetry_table = \
      special_position_settings.site_symmetry_table(
        sites_cart = model.get_sites_cart(),
        unconditional_general_position_flags=(
          model.get_atoms().extract_occ() != 1))
  original_spi = site_symmetry_table.special_position_indices()

  t2 = time()
  # Actually idelizing SS elements
  fixed_ss_selection = flex.bool(n_atoms_in_real_h, False)
  log.write("Replacing ss-elements with ideal ones:\n")
  log.flush()
  ss_stats = gather_ss_stats(pdb_h=model.get_hierarchy())
  n_idealized_elements = 0
  master_bool_sel = model.get_master_selection()
  if master_bool_sel is None or master_bool_sel.size() == 0:
    master_bool_sel = flex.bool(model.get_number_of_atoms(), True)
  elif isinstance(master_bool_sel, flex.size_t):
    master_bool_sel = flex.bool(model.get_number_of_atoms(), master_bool_sel)
  assert master_bool_sel.size() == model.get_number_of_atoms()
  for h in ann.helices:
    log.write("  %s\n" % h.as_pdb_str())
    log.flush()
    if processed_params.skip_good_ss_elements and ss_element_is_good(ss_stats, ([h],[])):
      log.write("    skipping, good element.\n")
    else:
      selstring = h.as_atom_selections()
      sel = selection_cache.selection(selstring[0])
      isel = sel.iselection()
      if (master_bool_sel & sel).iselection().size() == 0:
        log.write("    skipping, not in NCS master.\n")
        continue
      n_idealized_elements += 1
      log.write("    substitute with idealized one.\n")
      fixed_ss_selection.set_selected(isel, True)
      all_bsel = flex.bool(n_atoms_in_real_h, False)
      all_bsel.set_selected(isel, True)
      sel_h = model.get_hierarchy().select(all_bsel, copy_atoms=True)
      ideal_h = get_helix(helix_class=h.helix_class,
                          pdb_hierarchy_template=sel_h,
                          rotamer_manager=model.get_rotamer_manager())
      # edited_h.select(all_bsel).atoms().set_xyz(ideal_h.atoms().extract_xyz())
      set_xyz_carefully(dest_h=edited_h.select(all_bsel), source_h=ideal_h)
      # set_xyz_smart(dest_h=edited_h.select(all_bsel), source_h=ideal_h) # does not work here
  for sh in ann.sheets:
    s = "  %s\n" % sh.as_pdb_str()
    ss = s.replace("\n", "\n  ")
    log.write(ss[:-2])
    log.flush()
    if processed_params.skip_good_ss_elements and ss_element_is_good(ss_stats, ([],[sh])):
      log.write("    skipping, good element.\n")
    else:
      full_sh_selection = flex.bool(n_atoms_in_real_h, False)
      for st in sh.strands:
        selstring = st.as_atom_selections()
        isel = selection_cache.iselection(selstring)
        full_sh_selection.set_selected(isel, True)
      if (master_bool_sel & full_sh_selection).iselection().size() == 0:
        log.write("    skipping, not in NCS master.\n")
        continue
      n_idealized_elements += 1
      log.write("    substitute with idealized one.\n")
      for st in sh.strands:
        selstring = st.as_atom_selections()
        isel = selection_cache.iselection(selstring)
        all_bsel = flex.bool(n_atoms_in_real_h, False)
        all_bsel.set_selected(isel, True)
        fixed_ss_selection.set_selected(isel, True)
        sel_h = model.get_hierarchy().select(all_bsel, copy_atoms=True)
        ideal_h = secondary_structure_from_sequence(
            pdb_str=beta_pdb_str,
            sequence=None,
            pdb_hierarchy_template=sel_h,
            rotamer_manager=model.get_rotamer_manager(),
            )
        set_xyz_carefully(edited_h.select(all_bsel), ideal_h)
        # edited_h.select(all_bsel).atoms().set_xyz(ideal_h.atoms().extract_xyz())
  if n_idealized_elements == 0:
    log.write("Nothing was idealized.\n")
    # Don't do geometry minimization and stuff if nothing was changed.
    return None

  # XXX here we want to adopt new coordinates
  model.set_sites_cart(sites_cart=edited_h.atoms().extract_xyz())
  if model.ncs_constraints_present():
    model.set_sites_cart_from_hierarchy(multiply_ncs=True)

  t3 = time()
  # pre_result_h = edited_h
  # pre_result_h.reset_i_seq_if_necessary()
  bsel = flex.bool(n_atoms_in_real_h, False)
  helix_selection = flex.bool(n_atoms_in_real_h, False)
  sheet_selection = flex.bool(n_atoms_in_real_h, False)
  other_selection = flex.bool(n_atoms_in_real_h, False)
  ss_for_tors_selection = flex.bool(n_atoms_in_real_h, False)
  nonss_for_tors_selection = flex.bool(n_atoms_in_real_h, False)
  # set all CA atoms to True for other_selection
  #isel = selection_cache.iselection("name ca")
  isel = selection_cache.iselection("name ca or name n or name o or name c")
  other_selection.set_selected(isel, True)
  n_main_chain_atoms = other_selection.count(True)
  isel = selection_cache.iselection("name ca or name n or name o or name c")
  nonss_for_tors_selection.set_selected(isel, True)
  main_chain_selection_prefix = "(name ca or name n or name o or name c) %s"

  t4 = time()
  print("Preparing selections...", file=log)
  log.flush()
  # Here we are just preparing selections
  for h in ann.helices:
    ss_sels = h.as_atom_selections()[0]
    selstring = main_chain_selection_prefix % ss_sels
    isel = selection_cache.iselection(selstring)
    helix_selection.set_selected(isel, True)
    other_selection.set_selected(isel, False)
    isel = selection_cache.iselection(selstring)
    ss_for_tors_selection.set_selected(isel, True)
    nonss_for_tors_selection.set_selected(isel, False)

  for sheet in ann.sheets:
    for ss_sels in sheet.as_atom_selections():
      selstring = main_chain_selection_prefix % ss_sels
      isel = selection_cache.iselection(selstring)
      sheet_selection.set_selected(isel, True)
      other_selection.set_selected(isel, False)
      isel = selection_cache.iselection(selstring)
      ss_for_tors_selection.set_selected(isel, True)
      nonss_for_tors_selection.set_selected(isel, False)
  t5 = time()
  isel = selection_cache.iselection(
      "not name ca and not name n and not name o and not name c")
  other_selection.set_selected(isel, False)
  helix_sheet_intersection = helix_selection & sheet_selection
  if helix_sheet_intersection.count(True) > 0:
    sheet_selection = sheet_selection & ~helix_sheet_intersection
  assert ((helix_selection | sheet_selection) & other_selection).count(True)==0

  from mmtbx.monomer_library.pdb_interpretation import grand_master_phil_str
  params_line = grand_master_phil_str
  params_line += "secondary_structure {%s}" % secondary_structure.sec_str_master_phil_str
  # print "params_line"
  # print params_line
  params = iotbx.phil.parse(input_string=params_line, process_includes=True)#.extract()
  # This does not work the same way for a strange reason. Need to investigate.
  # The number of resulting hbonds is different later.
  # w_params = params.extract()
  # w_params.pdb_interpretation.secondary_structure.protein.remove_outliers = False
  # w_params.pdb_interpretation.peptide_link.ramachandran_restraints = True
  # w_params.pdb_interpretation.c_beta_restraints = True
  # w_params.pdb_interpretation.secondary_structure.enabled = True
  # params.format(python_object=w_params)
  # params.show()
  # print "="*80
  # print "="*80
  # print "="*80
  grm = model.get_restraints_manager()
  ssm_log = null_out()
  if verbose:
    ssm_log = log
  ss_params = secondary_structure.sec_str_master_phil.fetch().extract()
  ss_params.secondary_structure.protein.remove_outliers=False
  ss_manager = secondary_structure.manager(
      pdb_hierarchy=model.get_hierarchy(),
      geometry_restraints_manager=grm.geometry,
      sec_str_from_pdb_file=ss_annotation,
      params=ss_params.secondary_structure,
      mon_lib_srv=None,
      verbose=-1,
      log=ssm_log)
  grm.geometry.set_secondary_structure_restraints(
      ss_manager=ss_manager,
      hierarchy=model.get_hierarchy(),
      log=ssm_log)
  model.get_hierarchy().reset_i_seq_if_necessary()
  from mmtbx.geometry_restraints import reference
  if reference_map is None:
    if verbose:
      print("Adding reference coordinate restraints...", file=log)
    grm.geometry.append_reference_coordinate_restraints_in_place(
        reference.add_coordinate_restraints(
            sites_cart = model.get_sites_cart().select(helix_selection),
            selection  = helix_selection,
            sigma      = processed_params.sigma_on_reference_helix))
    grm.geometry.append_reference_coordinate_restraints_in_place(
        reference.add_coordinate_restraints(
            sites_cart = model.get_sites_cart().select(sheet_selection),
            selection  = sheet_selection,
            sigma      = processed_params.sigma_on_reference_sheet))
    grm.geometry.append_reference_coordinate_restraints_in_place(
        reference.add_coordinate_restraints(
            sites_cart = model.get_sites_cart().select(other_selection),
            selection  = other_selection,
            sigma      = processed_params.sigma_on_reference_non_ss))

  # XXX Somewhere here we actually should check placed side-chains for
  # clashes because we used ones that were in original model and just moved
  # them to nearest allowed rotamer. The idealization may affect a lot
  # the orientation of side chain thus justifying changing rotamer on it
  # to avoid clashes.
  if processed_params.fix_rotamer_outliers:
    print("Fixing/checking rotamers...", file=log)
    # pre_result_h.write_pdb_file(file_name="before_rotamers.pdb")
    br_txt = model.model_as_pdb()
    with open("before_rotamers.pdb", 'w') as f:
      f.write(br_txt)
    if(reference_map is None):
      backbone_sample=False
    else:
      backbone_sample=True
    result = mmtbx.refinement.real_space.fit_residues.run(
        pdb_hierarchy     = model.get_hierarchy(),
        crystal_symmetry  = model.crystal_symmetry(),
        map_data          = reference_map,
        rotamer_manager   = mmtbx.idealized_aa_residues.rotamer_manager.load(),
        sin_cos_table     = scitbx.math.sin_cos_table(n=10000),
        backbone_sample   = backbone_sample,
        mon_lib_srv       = model.get_mon_lib_srv(),
        log               = log)
    model.set_sites_cart(
        sites_cart = result.pdb_hierarchy.atoms().extract_xyz(),
        update_grm = True)

  if verbose:
    print("Adding chi torsion restraints...", file=log)
  # only backbone
  grm.geometry.add_chi_torsion_restraints_in_place(
          pdb_hierarchy   = model.get_hierarchy(),
          sites_cart      = model.get_sites_cart().\
                                 select(ss_for_tors_selection),
          selection = ss_for_tors_selection,
          chi_angles_only = False,
          sigma           = processed_params.sigma_on_torsion_ss)
  grm.geometry.add_chi_torsion_restraints_in_place(
          pdb_hierarchy   = model.get_hierarchy(),
          sites_cart      = model.get_sites_cart().\
                                select(nonss_for_tors_selection),
          selection = nonss_for_tors_selection,
          chi_angles_only = False,
          sigma           = processed_params.sigma_on_torsion_nonss)

  # real_h.atoms().set_xyz(pre_result_h.atoms().extract_xyz())
  #
  # Check and correct for special positions
  #
  real_h = model.get_hierarchy() # just a shortcut here...
  special_position_settings = crystal.special_position_settings(
      crystal_symmetry = model.crystal_symmetry())
  site_symmetry_table = \
      special_position_settings.site_symmetry_table(
        sites_cart = model.get_sites_cart(),
        unconditional_general_position_flags=(
          model.get_atoms().extract_occ() != 1))
  spi = site_symmetry_table.special_position_indices()
  if spi.size() > 0:
    print("Moving atoms from special positions:", file=log)
    for spi_i in spi:
      if spi_i not in original_spi:
        new_coords = (
            real_h.atoms()[spi_i].xyz[0]+0.2,
            real_h.atoms()[spi_i].xyz[1]+0.2,
            real_h.atoms()[spi_i].xyz[2]+0.2)
        print("  ", real_h.atoms()[spi_i].id_str(), end=' ', file=log)
        print(tuple(real_h.atoms()[spi_i].xyz), "-->", new_coords, file=log)
        real_h.atoms()[spi_i].set_xyz(new_coords)
  model.set_sites_cart_from_hierarchy()


  t9 = time()
  if processed_params.file_name_before_regularization is not None:
    grm.geometry.pair_proxies(sites_cart=model.get_sites_cart())
    grm.geometry.update_ramachandran_restraints_phi_psi_targets(
        hierarchy=model.get_hierarchy())
    print("Outputting model before regularization %s" % processed_params.file_name_before_regularization, file=log)

    m_txt = model.model_as_pdb()
    g_txt = model.restraints_as_geo()
    with open(processed_params.file_name_before_regularization, 'w') as f:
      f.write(m_txt)

    geo_fname = processed_params.file_name_before_regularization[:-4]+'.geo'
    print("Outputting geo file for regularization %s" % geo_fname, file=log)
    with open(geo_fname, 'w') as f:
      f.write(g_txt)

  #testing number of restraints
  assert grm.geometry.get_n_den_proxies() == 0
  if reference_map is None:
    assert grm.geometry.get_n_reference_coordinate_proxies() == n_main_chain_atoms, "" +\
        "%d %d" % (grm.geometry.get_n_reference_coordinate_proxies(), n_main_chain_atoms)
  refinement_log = null_out()
  log.write(
      "Refining geometry of substituted secondary structure elements...")
  log.flush()
  if verbose:
    refinement_log = log
  t10 = time()
  if reference_map is None:
    minimize_wrapper_for_ramachandran(
        model = model,
        original_pdb_h = None,
        excl_string_selection = "",
        log = refinement_log,
        number_of_cycles = processed_params.n_iter)
  else:
    ref_xrs = model.crystal_symmetry()
    minimize_wrapper_with_map(
        model = model,
        target_map=reference_map,
        refine_ncs_operators=False,
        number_of_cycles=processed_params.n_macro,
        log=log)
  model.set_sites_cart_from_hierarchy()

  log.write(" Done\n")
  log.flush()
  t11 = time()
  # print >> log, "Initial checking, init   : %.4f" % (t1-t0)
  # print >> log, "Checking SS              : %.4f" % (t2-t1)
  # print >> log, "Initializing selections  : %.4f" % (t4-t3)
  # print >> log, "Looping for selections   : %.4f" % (t5-t4)
  # print >> log, "Finalizing selections    : %.4f" % (t6-t5)
  # print >> log, "PDB interpretation       : %.4f" % (t7-t6)
  # print >> log, "Get GRM                  : %.4f" % (t8-t7)
  # print >> log, "Adding restraints to GRM : %.4f" % (t9-t8)
  # print >> log, "Running GM               : %.4f" % (t11-t10)
  # print_hbond_proxies(grm.geometry,real_h)
  grm.geometry.remove_reference_coordinate_restraints_in_place()
  grm.geometry.remove_chi_torsion_restraints_in_place(nonss_for_tors_selection)
  return grm.geometry.get_chi_torsion_proxies()


def beta():
  pdb_hierarchy = secondary_structure_from_sequence(beta_pdb_str,
      "ACEDGFIHKMLNQPSRTWVY")
  pdb_hierarchy.write_pdb_file(file_name = "o_beta_seq.pdb")

def alpha_310():
  pdb_hierarchy = secondary_structure_from_sequence(a310_helix_str,
      "ACEDGFIHKMLNQPSRTWVY")
  pdb_hierarchy.write_pdb_file(file_name = "o_helix310_seq.pdb")

def alpha_pi():
  pdb_hierarchy = secondary_structure_from_sequence(pi_helix_str,
      "ACEDGFIHKMLNQPSRTWVY")
  pdb_hierarchy.write_pdb_file(file_name = "o_helix_pi_seq.pdb")

def alpha():
  pdb_hierarchy = secondary_structure_from_sequence(alpha_helix_str,
      "ACEDGFIHKMLNQPSRTWVY")
  pdb_hierarchy.write_pdb_file(file_name = "o_helix_seq.pdb")
