from __future__ import division
from scitbx.math import superpose
import iotbx.pdb
from cctbx.array_family import flex
from mmtbx.monomer_library import idealized_aa
from libtbx.utils import Sorry, null_out
from iotbx.pdb.amino_acid_codes import one_letter_given_three_letter as one_three
from iotbx.pdb.amino_acid_codes import three_letter_given_one_letter as three_one
from mmtbx.refinement.geometry_minimization import run2
from mmtbx.rotamer.rotamer_eval import RotamerEval
import mmtbx.utils
from iotbx.pdb import secondary_structure as ioss

global_counter = 1

alpha_pdb_str = """\
ATOM      1  N   ALA     2       1.643  -2.366  -1.408  1.00  0.00           N
ATOM      3  CA  ALA     2       1.280  -3.608  -2.069  1.00  0.00           C
ATOM     10  C   ALA     2      -0.114  -3.466  -2.684  1.00  0.00           C
ATOM     11  O   ALA     2      -0.327  -3.827  -3.840  1.00  0.00           O
ATOM      6  CB  ALA     2       1.361  -4.762  -1.068  1.00  0.00           C
ATOM     12  N   ALA     3      -1.028  -2.938  -1.882  1.00  0.00           N
ATOM     14  CA  ALA     3      -2.395  -2.743  -2.333  1.00  0.00           C
ATOM     21  C   ALA     3      -2.396  -1.855  -3.579  1.00  0.00           C
ATOM     22  O   ALA     3      -3.059  -2.167  -4.567  1.00  0.00           O
ATOM     17  CB  ALA     3      -3.228  -2.150  -1.194  1.00  0.00           C
"""

alpha310_pdb_str = """\
ATOM      1  N   ALA    1       -1.204  -0.514   0.643  1.00  0.00           N
ATOM      1  CA  ALA    1        0.000   0.000   0.000  1.00  0.00           C
ATOM      1  C   ALA    1        0.804  -1.124  -0.644  1.00  0.00           C
ATOM      1  O   ALA    1        1.628  -0.884  -1.526  1.00  0.00           O
ATOM      1  CB  ALA    1        0.870   0.757   1.006  1.00  0.00           C
ATOM      1  N   ALA    2        0.559  -2.352  -0.197  1.00  0.00           N
ATOM      1  CA  ALA    2        1.260  -3.515  -0.728  1.00  0.00           C
ATOM      1  C   ALA    2        1.116  -3.602  -2.244  1.00  0.00           C
ATOM      1  O   ALA    2        1.905  -4.266  -2.915  1.00  0.00           O
ATOM      1  CB  ALA    2        0.743  -4.801  -0.079  1.00  0.00           C
"""

alpha_pi_pdb_str = """\
ATOM      1  N   ALA     2       2.054  -2.383  -1.604  1.00  0.00           N
ATOM      3  CA  ALA     2       1.733  -3.620  -2.296  1.00  0.00           C
ATOM     10  C   ALA     2       0.216  -3.730  -2.460  1.00  0.00           C
ATOM     11  O   ALA     2      -0.301  -3.624  -3.570  1.00  0.00           O
ATOM      6  CB  ALA     2       2.324  -4.804  -1.527  1.00  0.00           C
ATOM     12  N   ALA     3      -0.454  -3.940  -1.336  1.00  0.00           N
ATOM     14  CA  ALA     3      -1.902  -4.065  -1.341  1.00  0.00           C
ATOM     21  C   ALA     3      -2.516  -2.807  -1.959  1.00  0.00           C
ATOM     22  O   ALA     3      -3.064  -2.855  -3.059  1.00  0.00           O
ATOM     17  CB  ALA     3      -2.398  -4.316   0.085  1.00  0.00           C
"""

beta_pdb_str = """\
ATOM      1  N   ALA    1       -1.204  -0.514   0.643  1.00  0.00           N
ATOM      1  CA  ALA    1        0.000   0.000   0.000  1.00  0.00           C
ATOM      1  C   ALA    1        1.258  -0.522   0.685  1.00  0.00           C
ATOM      1  O   ALA    1        1.316  -0.616   1.911  1.00  0.00           O
ATOM      1  CB  ALA    1       -0.000   1.530   0.000  1.00  0.00           C
ATOM      1  N   ALA    2        2.265  -0.860  -0.115  1.00  0.00           N
ATOM      1  CA  ALA    2        3.524  -1.373   0.412  1.00  0.00           C
ATOM      1  C   ALA    2        4.694  -0.478   0.018  1.00  0.00           C
ATOM      1  O   ALA    2        4.758   0.018  -1.107  1.00  0.00           O
ATOM      1  CB  ALA    2        3.770  -2.802  -0.076  1.00  0.00           C
"""

helix_class_to_pdb_str = {1:alpha_pdb_str, 3:alpha_pi_pdb_str, 5: alpha310_pdb_str}



def get_expected_n_hbonds_from_helix(helix):
  # assess helix length
  if helix.helix_class==1:
    return helix.length-4
  elif helix.helix_class==5:
    return helix.length-3
  elif htype==3:
    return helix.length-5
  else:
    raise Sorry("Unsupported helix type.")

def print_hbond_proxies(geometry, hierarchy):
  """ Print hydrogen bonds in geometry restraints manager for debugging
  purposes"""
  atoms = hierarchy.atoms()
  for hb in geometry.generic_restraints_manager.hbonds_as_simple_bonds():
    print (atoms[hb[0]].id_str(), "<====>",atoms[hb[1]].id_str(),
        atoms[hb[0]].distance(atoms[hb[1]]), hb[0], hb[1])

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



def make_ss_structure_from_sequence(pdb_str, sequence=None,
      pdb_hierarchy_template=None):
  """ Return pdb.hierarchy with secondary structure according to sequence or
  reference hierarcy. If reference hierarchy provided, the resulting hierarchy
  will be rigid body refined to it.

  pdb_str - "ideal" structure at least 2 residues long.
  sequence - string with sequence (one-letter codes)
  pdb_hierarchy_template - reference hierarchy.
  """
  pht = pdb_hierarchy_template
  assert [sequence, pht].count(None) == 1
  if pht:
    lk = len(pht.altloc_indices().keys())
    if lk ==0:
      raise Sorry(
          "Hierarchy template in make_ss_structure_from_sequence is empty")
    else:
      assert len(pht.altloc_indices().keys()) == 1, \
          "Alternative conformations are not supported"
  number_of_residues = len(sequence) if sequence!=None else \
    len(pht.models()[0].chains()[0].conformers()[0].residues())
  if number_of_residues<1:
    raise Sorry('sequence should contain at least one residue.')
  r, t = get_r_t_matrices_from_structure(pdb_str)
  ideal_res_dict = idealized_aa.residue_dict()
  real_res_list = None
  if pht:
    real_res_list = pht.models()[0].chains()[0].residue_groups()
  pdb_hierarchy = iotbx.pdb.input(source_info=None, lines=pdb_str).\
      construct_hierarchy()
  chain = pdb_hierarchy.models()[0].chains()[0]
  current_ala_ag = chain.residue_groups()[0].atom_groups()[0]
  new_chain = iotbx.pdb.hierarchy.chain(id=" ")
  new_chain.pre_allocate_residue_groups(number_of_additional_residue_groups=\
                                                            number_of_residues)
  rotamer_manager = RotamerEval()
  for j in range(number_of_residues):
    # put ALA
    c = sequence[j] if sequence!=None else \
        one_three[real_res_list[j].atom_groups()[0].resname.upper()]
    rg = iotbx.pdb.hierarchy.residue_group(icode="")
    rg.resseq = j+1
    new_chain.append_residue_group(residue_group=rg)
    ag_to_place = current_ala_ag.detached_copy()
    rg.append_atom_group(atom_group=ag_to_place)
    current_ala_ag.atoms().set_xyz(
                          r.elems*current_ala_ag.atoms().extract_xyz()+t.elems)
    if c == 'A':
      continue
    ag_to_place.resname = three_one[c]
    if c == 'G':
      for a in ag_to_place.atoms():
        if a.name.strip() == "CB":
          ag_to_place.remove_atom(atom=a)
          break
      continue
    # align residue from ideal_res_dict to just placed ALA (ag_to_place)
    # or from pdb_hierarchy_template
    fixed_sites = flex.vec3_double()
    moving_sites = flex.vec3_double()
    reper_atoms = ["CB","CA", "N"]
    current_reference_ag = real_res_list[j].atom_groups()[0] if pht else \
        ideal_res_dict[three_one[c].lower()].models()[0].chains()[0].\
        residue_groups()[0].atom_groups()[0]
    for (ag, arr) in [(ag_to_place, fixed_sites),
        (current_reference_ag, moving_sites)]:
      for a in ag.atoms():
        if a.name.strip() in reper_atoms:
          arr.append(a.xyz)
    lsq_fit_obj = superpose.least_squares_fit(reference_sites = fixed_sites,
                                              other_sites = moving_sites)
    ideal_correct_ag = current_reference_ag.detached_copy()
    ideal_correct_ag.atoms().set_xyz(
        lsq_fit_obj.r.elems*ideal_correct_ag.atoms().extract_xyz()+\
        lsq_fit_obj.t.elems)
    ideal_correct_ag.atoms().set_xyz(
        rotamer_manager.nearest_rotamer_sites_cart(ideal_correct_ag))
    ag_to_place.pre_allocate_atoms(number_of_additional_atoms=\
                                                len(ideal_correct_ag.atoms())-5)
    for a in ideal_correct_ag.atoms():
      if a.name.strip() not in ["N","CA","C","O", "CB"]:
        at = a.detached_copy()
        at.uij_erase()
        ag_to_place.append_atom(atom=at)
  new_pdb_h = iotbx.pdb.hierarchy.new_hierarchy_from_chain(new_chain)
  # align to real
  if pht != None:
    fixed_sites = pht.atoms().extract_xyz()
    moving_sites = new_pdb_h.atoms().extract_xyz()
    assert len(fixed_sites) == len(moving_sites)
    lsq_fit_obj = superpose.least_squares_fit(reference_sites = fixed_sites,
                                              other_sites = moving_sites)
    new_pdb_h.atoms().set_xyz(
        lsq_fit_obj.r.elems*new_pdb_h.atoms().extract_xyz()+lsq_fit_obj.t.elems)
  return new_pdb_h

def get_helix(helix_class, sequence=None, pdb_hierarchy_template=None):
  if helix_class not in helix_class_to_pdb_str.keys():
    raise Sorry("Unsupported helix type.")
  return make_ss_structure_from_sequence(
    pdb_str=helix_class_to_pdb_str[helix_class],
    sequence=sequence,
    pdb_hierarchy_template=pdb_hierarchy_template)

def substitute_ss(real_h,
                    crystal_symmetry,
                    helices,
                    sigma_on_reference_non_ss = 1,
                    sigma_on_reference_ss = 5,
                    sigma_on_torsion_ss = 5,
                    log=null_out()):
  """
  Substitute secondary structure elements in real_h hierarchy with ideal
  ones _in_place_. Currently only helices are supported.
  Returns geometry restraints manager for furhter refinements with all
  correct hydrogen bonds. It is not guaranteed that
  mmtbx.secondary_structure.manager.find_automatically()
  will be able to generate them again.
  real_h - hierarcy to substitute secondary structure elements.
  crystal_symmetry - symmetry object for the hierarchy provided.
  helices - list with HELIX records. Types supported:
      1:alpha_pdb_str, 3:alpha_pi_pdb_str, 5: alpha310_pdb_str

  Weights (bigger number gives lighter restraints):
  sigma_on_reference_non_ss - weight on original model coordinates restraints
      where no ss present. Keeps loops close to initial model.
  sigma_on_reference_ss - weight on original model coordinates restraints
      where ss present. Bends helices a bit according to initial model.
  sigma_on_torsion_ss - weight on torsion angles restraints where ss present.
      Keeps helices torsion angles close to ideal.
  """
  expected_n_hbonds = 0
  ann = ioss.annotation(helices=helices, sheets=[])
  phil_str = ann.as_restraint_groups()

  for h in helices:
    expected_n_hbonds += get_expected_n_hbonds_from_helix(h)
  edited_h = iotbx.pdb.input(source_info=None,
      lines=flex.split_lines(real_h.as_pdb_string())).construct_hierarchy()
  n_atoms_in_real_h = real_h.atoms().size()
  cumm_bsel = flex.bool(n_atoms_in_real_h, False)
  selection_cache = real_h.atom_selection_cache()
  for h in helices:
    selstring = h.as_atom_selections(params=None)
    isel = selection_cache.iselection(selstring[0])
    all_bsel = flex.bool(n_atoms_in_real_h, False)
    all_bsel.set_selected(isel, True)
    cumm_bsel.set_selected(isel, True)
    sel_h = real_h.select(all_bsel, copy_atoms=True)
    ideal_h = get_helix(helix_class=h.helix_class, pdb_hierarchy_template=sel_h)
    edited_h.select(all_bsel).atoms().set_xyz(ideal_h.atoms().extract_xyz())
  pre_result_h = edited_h
  # generate different atom selections.
  bsel = flex.bool(real_h.atoms().size(), False)
  ss_selection = flex.bool(real_h.atoms().size(), False)
  other_selection = flex.bool(real_h.atoms().size(), False)
  ss_for_tors_selection = flex.bool(real_h.atoms().size(), False)
  selection_cache = real_h.atom_selection_cache()
  # set all CA atoms to True for other_selection
  isel = selection_cache.iselection("name ca")
  other_selection.set_selected(isel, True)
  for h in helices:
    selstring = "name ca chain %s resseq %s:%s"\
                % (h.start_chain_id, h.start_resseq, h.end_resseq)
    isel = selection_cache.iselection(selstring)
    ss_selection.set_selected(isel, True)
    other_selection.set_selected(isel, False)
    selstring = "(name ca or name n or name o or name c) chain %s resseq %s:%s"\
                % (h.start_chain_id, h.start_resseq, h.end_resseq)
    isel = selection_cache.iselection(selstring)
    ss_for_tors_selection.set_selected(isel, True)
  processed_pdb_files_srv = mmtbx.utils.\
      process_pdb_file_srv(crystal_symmetry= crystal_symmetry, log=log)
  processed_pdb_file, junk = processed_pdb_files_srv.\
      process_pdb_files(raw_records=flex.split_lines(real_h.as_pdb_string()))
  defpars = mmtbx.secondary_structure.sec_str_master_phil.fetch()
  custom_pars = defpars.fetch(source = iotbx.phil.parse(
      "h_bond_restraints.remove_outliers=False\n%s" % phil_str))
  ss_manager = mmtbx.secondary_structure.manager(
      pdb_hierarchy=pre_result_h, params=custom_pars.extract())
  #ss_manager.find_automatically(log=log)
  #ss_manager.initialize(log=log)
  proxies_for_grm = ss_manager.create_hbond_proxies(
      log          = log,
      as_python_objects = False)
  custom_nb_excl = proxies_for_grm.exclude_nb_list
  grm = processed_pdb_file.geometry_restraints_manager(
      show_energies                = False,
      show_nonbonded_clashscore    = False,
      hydrogen_bond_proxies=proxies_for_grm.proxies,
      custom_nonbonded_exclusions = custom_nb_excl,
      assume_hydrogens_all_missing = True)
  grm.generic_restraints_manager.reference_manager.\
      add_coordinate_restraints(
          sites_cart = real_h.atoms().extract_xyz().select(ss_selection),
          selection  = ss_selection,
          sigma      = sigma_on_reference_ss)
  grm.generic_restraints_manager.reference_manager.\
      add_coordinate_restraints(
          sites_cart = real_h.atoms().extract_xyz().select(other_selection),
          selection  = other_selection,
          sigma      = sigma_on_reference_non_ss)
  grm.generic_restraints_manager.reference_manager.\
      add_coordinate_restraints(
          sites_cart = real_h.atoms().extract_xyz().select(other_selection),
          selection  = other_selection,
          sigma      = sigma_on_reference_non_ss)
  grm.generic_restraints_manager.reference_manager.\
      add_torsion_restraints(
          pdb_hierarchy   = pre_result_h,
          sites_cart      = pre_result_h.atoms().extract_xyz().\
                                select(ss_for_tors_selection),
          selection = ss_for_tors_selection,
          chi_angles_only = False,
          sigma           = sigma_on_torsion_ss)
  real_h.atoms().set_xyz(pre_result_h.atoms().extract_xyz())
  restraints_manager = mmtbx.restraints.manager(
      geometry=grm,
      normalization=True)
  assert restraints_manager.geometry.generic_restraints_manager.\
      get_n_hbonds() == expected_n_hbonds
  #real_h.write_pdb_file(file_name="before_geom_reg.pdb")
  obj = run2(
      restraints_manager       = restraints_manager,
      pdb_hierarchy            = real_h,
      max_number_of_iterations = 500,
      number_of_macro_cycles   = 5,
      bond                     = True,
      nonbonded                = True,
      angle                    = True,
      dihedral                 = True,
      chirality                = True,
      planarity                = True,
      log                      = log)
  # removing unnecessary restraints!!!
  restraints_manager.geometry.generic_restraints_manager.reference_manager.\
      remove_coordinate_restraints(selection=ss_selection)
  restraints_manager.geometry.generic_restraints_manager.reference_manager.\
      remove_coordinate_restraints(selection=other_selection)
  return restraints_manager


def beta():
  pdb_hierarchy = make_ss_structure_from_sequence(beta_pdb_str,
      "ACEDGFIHKMLNQPSRTWVY")
  pdb_hierarchy.write_pdb_file(file_name = "o_beta_seq.pdb")

def alpha_310():
  pdb_hierarchy = make_ss_structure_from_sequence(alpha310_pdb_str,
      "ACEDGFIHKMLNQPSRTWVY")
  pdb_hierarchy.write_pdb_file(file_name = "o_helix310_seq.pdb")

def alpha_pi():
  pdb_hierarchy = make_ss_structure_from_sequence(alpha_pi_pdb_str,
      "ACEDGFIHKMLNQPSRTWVY")
  pdb_hierarchy.write_pdb_file(file_name = "o_helix_pi_seq.pdb")

def alpha():
  pdb_hierarchy = make_ss_structure_from_sequence(alpha_pdb_str,
      "ACEDGFIHKMLNQPSRTWVY")
  pdb_hierarchy.write_pdb_file(file_name = "o_helix_seq.pdb")
