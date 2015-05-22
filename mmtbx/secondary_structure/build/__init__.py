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
#from mmtbx.command_line.geometry_minimization import master_params
from mmtbx import secondary_structure
from mmtbx.monomer_library.pdb_interpretation import grand_master_phil_str

alpha_helix_str = """
ATOM      1  N   ALA A   1      -5.606  -2.251 -12.878  1.00  0.00           N
ATOM      2  CA  ALA A   1      -5.850  -1.194 -13.852  1.00  0.00           C
ATOM      3  C   ALA A   1      -5.186  -1.524 -15.184  1.00  0.00           C
ATOM      4  O   ALA A   1      -5.744  -1.260 -16.249  1.00  0.00           O
ATOM      5  CB  ALA A   1      -5.350   0.142 -13.324  1.00  0.00           C
ATOM      6  N   ALA A   2      -3.991  -2.102 -15.115  1.00  0.00           N
ATOM      7  CA  ALA A   2      -3.262  -2.499 -16.313  1.00  0.00           C
ATOM      8  C   ALA A   2      -3.961  -3.660 -17.011  1.00  0.00           C
ATOM      9  O   ALA A   2      -4.016  -3.716 -18.240  1.00  0.00           O
ATOM     10  CB  ALA A   2      -1.829  -2.872 -15.965  1.00  0.00           C
"""

a310_helix_str = """\
ATOM      1  N   ALA A   1       8.836  -4.233 -14.408  1.00  0.00           N
ATOM      2  CA  ALA A   1      10.232  -4.071 -14.799  1.00  0.00           C
ATOM      3  C   ALA A   1      10.764  -5.331 -15.476  1.00  0.00           C
ATOM      4  O   ALA A   1      11.679  -5.262 -16.297  1.00  0.00           O
ATOM      5  CB  ALA A   1      11.088  -3.716 -13.593  1.00  0.00           C
ATOM      6  N   ALA A   2      10.176  -6.478 -15.143  1.00  0.00           N
ATOM      7  CA  ALA A   2      10.582  -7.741 -15.750  1.00  0.00           C
ATOM      8  C   ALA A   2      10.381  -7.714 -17.262  1.00  0.00           C
ATOM      9  O   ALA A   2      11.080  -8.410 -17.999  1.00  0.00           O
ATOM     10  CB  ALA A   2       9.815  -8.901 -15.134  1.00  0.00           C
"""

pi_helix_str = """\
ATOM      1  N   ALA A   1      -3.365  -3.446  -8.396  1.00  0.00           N
ATOM      2  CA  ALA A   1      -4.568  -4.249  -8.592  1.00  0.00           C
ATOM      3  C   ALA A   1      -5.809  -3.386  -8.805  1.00  0.00           C
ATOM      4  O   ALA A   1      -6.559  -3.591  -9.759  1.00  0.00           O
ATOM      5  CB  ALA A   1      -4.775  -5.185  -7.411  1.00  0.00           C
ATOM      6  N   ALA A   2      -6.025  -2.424  -7.914  1.00  0.00           N
ATOM      7  CA  ALA A   2      -7.221  -1.588  -7.976  1.00  0.00           C
ATOM      8  C   ALA A   2      -7.101  -0.486  -9.025  1.00  0.00           C
ATOM      9  O   ALA A   2      -8.089  -0.114  -9.659  1.00  0.00           O
ATOM     10  CB  ALA A   2      -7.511  -0.985  -6.610  1.00  0.00           C
"""

beta_pdb_str = """
ATOM      1  N   ALA A   1      27.961   0.504   1.988  1.00  0.00           N
ATOM      2  CA  ALA A   1      29.153   0.205   2.773  1.00  0.00           C
ATOM      3  C   ALA A   1      30.420   0.562   2.003  1.00  0.00           C
ATOM      4  O   ALA A   1      30.753  -0.077   1.005  1.00  0.00           O
ATOM      5  CB  ALA A   1      29.170  -1.262   3.172  1.00  0.00           C
ATOM      6  N   ALA A   2      31.123   1.587   2.474  1.00  0.00           N
ATOM      7  CA  ALA A   2      32.355   2.031   1.832  1.00  0.00           C
ATOM      8  C   ALA A   2      33.552   1.851   2.758  1.00  0.00           C
ATOM      9  O   ALA A   2      33.675   2.539   3.772  1.00  0.00           O
ATOM     10  CB  ALA A   2      32.232   3.483   1.399  1.00  0.00           C
"""

helix_class_to_pdb_str = {1:alpha_helix_str, 3:pi_helix_str, 5: a310_helix_str}

ss_idealization_master_phil_str = """
idealization
{
  enabled = False
    .type = bool
  restrain_torsion_angles = False
    .type = bool
  sigma_on_reference_non_ss = 1
    .type = float
  sigma_on_reference_helix = 1
    .type = float
  sigma_on_reference_sheet = 0.5
    .type = float
  sigma_on_torsion_ss = 5
    .type = float
  sigma_on_torsion_nonss = 5
    .type = float
  sigma_on_ramachandran = 1
    .type = float
  sigma_on_cbeta = 2.5
    .type = float
  n_macro=3
    .type = int
  n_iter=300
    .type = int
}
"""

def get_expected_n_hbonds_from_helix(helix):
  # assess helix length
  if helix.helix_class==1:
    return helix.length-4
  elif helix.helix_class==5:
    return helix.length-3
  elif helix.helix_class==3:
    return helix.length-5
  else:
    raise Sorry("Unsupported helix type.")

def print_hbond_proxies(geometry, hierarchy, pymol=False):
  """ Print hydrogen bonds in geometry restraints manager for debugging
  purposes"""
  atoms = hierarchy.atoms()
  if pymol:
    dashes = open('dashes.pml', 'w')
  hbondlen=flex.double()
  for hb in geometry.generic_restraints_manager.hbonds_as_simple_bonds():
    hbondlen.append(atoms[hb[0]].distance(atoms[hb[1]]))
    print (atoms[hb[0]].id_str(), "<====>",atoms[hb[1]].id_str(),
        atoms[hb[0]].distance(atoms[hb[1]]), hb[0], hb[1])
    if pymol:
      s1 = atoms[hb[0]].id_str()
      s2 = atoms[hb[1]].id_str()
      #print "pdbstr1:", s1
      #print "pdbstr1:",s2
      ps = "dist chain \"%s\" and resi %s and name %s, chain \"%s\" and resi %s and name %s\n" % (s1[14:15],
         s1[16:19], s1[5:7], s2[14:15], s2[16:19], s2[5:7])
      dashes.write(ps)
  print "min, max, mean, sd hbond lenghts", hbondlen.min_max_mean().as_tuple(),\
    hbondlen.standard_deviation_of_the_sample()
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
  c = one_three[current_reference_ag.resname.upper()]
  if c == 'A':
    return
  ag_to_place.resname = three_one[c]
  if c == 'G':
    for a in ag_to_place.atoms():
      if a.name.strip() == "CB":
        ag_to_place.remove_atom(atom=a)
        break
    return
  # align residue from ideal_res_dict to just placed ALA (ag_to_place)
  # or from pdb_hierarchy_template
  fixed_sites = flex.vec3_double()
  moving_sites = flex.vec3_double()
  reper_atoms = ["CB","CA", "N"]
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
      assert len(pht.altloc_indices().keys()) == 1, \
          "Alternative conformations are not supported"
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
  chain = pdb_hierarchy.models()[0].chains()[0]
  current_ala_ag = chain.residue_groups()[0].atom_groups()[0]
  new_chain = iotbx.pdb.hierarchy.chain(id="A")
  new_chain.pre_allocate_residue_groups(number_of_additional_residue_groups=\
                                                            number_of_residues)
  r, t = get_r_t_matrices_from_structure(pdb_str)
  for j in range(number_of_residues):
    # put ALA
    rg = iotbx.pdb.hierarchy.residue_group(icode="")
    rg.resseq = j+1
    new_chain.append_residue_group(residue_group=rg)
    ag_to_place = current_ala_ag.detached_copy()
    rg.append_atom_group(atom_group=ag_to_place)
    current_ala_ag.atoms().set_xyz(
                          r.elems*current_ala_ag.atoms().extract_xyz()+t.elems)
    current_reference_ag = real_res_list[j].atom_groups()[0] if pht else \
        ideal_res_dict[three_one[sequence[j]].lower()].models()[0].chains()[0].\
        residue_groups()[0].atom_groups()[0]
    side_chain_placement(ag_to_place, current_reference_ag, rotamer_manager)
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

def get_helix(helix_class, rotamer_manager, sequence=None, pdb_hierarchy_template=None):
  if helix_class not in helix_class_to_pdb_str.keys():
    raise Sorry("Unsupported helix type.")
  return secondary_structure_from_sequence(
    pdb_str=helix_class_to_pdb_str[helix_class],
    sequence=sequence,
    rotamer_manager=rotamer_manager,
    pdb_hierarchy_template=pdb_hierarchy_template)


def get_empty_ramachandran_proxies():
  import boost.python
  ext = boost.python.import_ext("mmtbx_ramachandran_restraints_ext")
  proxies = ext.shared_phi_psi_proxy()
  return proxies


def substitute_ss(real_h,
                    xray_structure,
                    ss_annotation,
                    sigma_on_reference_non_ss = 1,
                    sigma_on_reference_helix = 1,
                    sigma_on_reference_sheet = 0.5,
                    sigma_on_torsion_ss = 5,
                    sigma_on_torsion_nonss = 5,
                    sigma_on_ramachandran = 1,
                    sigma_on_cbeta = 2.5,
                    n_macro=3,
                    n_iter=300,
                    fname_before_regularization=None,
                    log=null_out(),
                    rotamer_manager=None,
                    verbose=False):
  """
  Substitute secondary structure elements in real_h hierarchy with ideal
  ones _in_place_.
  Returns reference torsion proxies - the only thing that cannot be restored
  with little effort outside the procedure.
  real_h - hierarcy to substitute secondary structure elements.
  xray_structure - xray_structure - needed to get crystal symmetry (to
      construct processed_pdb_file and xray_structure is needed to call
      get_geometry_restraints_manager for no obvious reason).
  ss_annotation - annotation records.

  Weights (bigger number gives lighter restraints):
  sigma_on_reference_non_ss - weight on original model coordinates restraints
      where no ss present. Keeps loops close to initial model.
  sigma_on_reference_ss - weight on original model coordinates restraints
      where ss present. Bends helices a bit according to initial model.
  sigma_on_torsion_ss - weight on torsion angles restraints where ss present.
      Keeps helices torsion angles close to ideal.
  """
  if rotamer_manager is None:
    rotamer_manager = RotamerEval()
  for model in real_h.models():
    for chain in model.chains():
      if len(chain.conformers()) > 1:
        raise Sorry("Secondary structure substitution does not support\n"+\
            "the presence of alternative conformations.")
  expected_n_hbonds = 0
  ann = ss_annotation
  phil_str = ann.as_restraint_groups()
  for h in ann.helices:
    expected_n_hbonds += get_expected_n_hbonds_from_helix(h)
  edited_h = real_h.deep_copy()
  n_atoms_in_real_h = real_h.atoms().size()
  cumm_bsel = flex.bool(n_atoms_in_real_h, False)
  selection_cache = real_h.atom_selection_cache()
  # check the annotation for correctness (atoms are actually in hierarchy)
  error_msg = "The following secondary structure annotations result in \n"
  error_msg +="empty atom selections. They don't match the structre: \n"
  error_flg = False
  for h in ann.helices:
    selstring = h.as_atom_selections()
    isel = selection_cache.iselection(selstring[0])
    if len(isel) == 0:
      error_flg = True
      error_msg += "  %s\n" % h
  for sh in ann.sheets:
    for st in sh.strands:
      selstring = st.as_atom_selections()
      isel = selection_cache.iselection(selstring)
      if len(isel) == 0:
        error_flg = True
        error_msg += "  %s\n" % sh.particular_strand_as_pdb_str(
                                      strand_id=st.strand_id)
  if error_flg:
    raise Sorry(error_msg)

  for h in ann.helices:
    selstring = h.as_atom_selections()
    isel = selection_cache.iselection(selstring[0])
    all_bsel = flex.bool(n_atoms_in_real_h, False)
    all_bsel.set_selected(isel, True)
    cumm_bsel.set_selected(isel, True)
    sel_h = real_h.select(all_bsel, copy_atoms=True)
    ideal_h = get_helix(helix_class=h.helix_class,
                        pdb_hierarchy_template=sel_h,
                        rotamer_manager=rotamer_manager)
    edited_h.select(all_bsel).atoms().set_xyz(ideal_h.atoms().extract_xyz())
  for sh in ann.sheets:
    for st in sh.strands:
      selstring = st.as_atom_selections()
      isel = selection_cache.iselection(selstring)
      all_bsel = flex.bool(n_atoms_in_real_h, False)
      all_bsel.set_selected(isel, True)
      cumm_bsel.set_selected(isel, True)
      sel_h = real_h.select(all_bsel, copy_atoms=True)
      ideal_h = secondary_structure_from_sequence(
          pdb_str=beta_pdb_str,
          sequence=None,
          pdb_hierarchy_template=sel_h,
          rotamer_manager=rotamer_manager,
          )
      edited_h.select(all_bsel).atoms().set_xyz(ideal_h.atoms().extract_xyz())

  pre_result_h = edited_h
  pre_result_h.reset_i_seq_if_necessary()
  n_atoms = real_h.atoms().size()
  bsel = flex.bool(n_atoms, False)
  helix_selection = flex.bool(n_atoms, False)
  sheet_selection = flex.bool(n_atoms, False)
  other_selection = flex.bool(n_atoms, False)
  ss_for_tors_selection = flex.bool(n_atoms, False)
  nonss_for_tors_selection = flex.bool(n_atoms, False)
  selection_cache = real_h.atom_selection_cache()
  # set all CA atoms to True for other_selection
  #isel = selection_cache.iselection("name ca")
  isel = selection_cache.iselection("name ca or name n or name o or name c")
  other_selection.set_selected(isel, True)
  n_main_chain_atoms = other_selection.count(True)
  isel = selection_cache.iselection("name ca or name n or name o or name c")
  nonss_for_tors_selection.set_selected(isel, True)
  main_chain_selection_prefix = "(name ca or name n or name o or name c) %s"
  if verbose:
    log.write("Replacing ss-elements with ideal ones:\n")
  for h in ann.helices:
    log.write("  %s\n" % h.as_pdb_str())
    ss_sels = h.as_atom_selections()[0]
    selstring = main_chain_selection_prefix % ss_sels
    isel = selection_cache.iselection(selstring)
    helix_selection.set_selected(isel, True)
    other_selection.set_selected(isel, False)
    isel = selection_cache.iselection(selstring)
    ss_for_tors_selection.set_selected(isel, True)
    nonss_for_tors_selection.set_selected(isel, False)

  for sheet in ann.sheets:
    log.write("  %s\n" % sheet.as_pdb_str())
    for ss_sels in sheet.as_atom_selections():
      selstring = main_chain_selection_prefix % ss_sels
      isel = selection_cache.iselection(selstring)
      sheet_selection.set_selected(isel, True)
      other_selection.set_selected(isel, False)
      isel = selection_cache.iselection(selstring)
      ss_for_tors_selection.set_selected(isel, True)
      nonss_for_tors_selection.set_selected(isel, False)

  isel = selection_cache.iselection(
      "not name ca and not name n and not name o and not name c")
  other_selection.set_selected(isel, False)
  helix_sheet_intersection = helix_selection & sheet_selection
  if helix_sheet_intersection.count(True) > 0:
    sheet_selection = sheet_selection & ~helix_sheet_intersection
  assert ((helix_selection | sheet_selection) & other_selection).count(True)==0

  params_line = grand_master_phil_str
  params_line += "secondary_structure {%s}" % secondary_structure.sec_str_master_phil_str
  params = iotbx.phil.parse(input_string=params_line, process_includes=True)
  custom_pars = params.fetch(source = iotbx.phil.parse("\n".join([
      "pdb_interpretation.secondary_structure {h_bond_restraints.remove_outliers = False\n%s}" \
          % phil_str,
      "pdb_interpretation.peptide_link.ramachandran_restraints = True",
      "c_beta_restraints = True",
      "pdb_interpretation.secondary_structure.enabled=True",
      "pdb_interpretation.secondary_structure.find_automatically=False"]))).extract()

  processed_pdb_files_srv = mmtbx.utils.\
      process_pdb_file_srv(
          crystal_symmetry= xray_structure.crystal_symmetry(),
          pdb_interpretation_params = custom_pars.pdb_interpretation,
          log=null_out())
  processed_pdb_file, junk = processed_pdb_files_srv.\
      process_pdb_files(raw_records=flex.split_lines(real_h.as_pdb_string()))
  has_hd = None
  if(xray_structure is not None):
    sctr_keys = xray_structure.scattering_type_registry().type_count_dict().keys()
    has_hd = "H" in sctr_keys or "D" in sctr_keys
  geometry = processed_pdb_file.geometry_restraints_manager(
    show_energies                = False,
    params_edits                 = custom_pars.geometry_restraints.edits,
    plain_pairs_radius           = 5,
    assume_hydrogens_all_missing = not has_hd)
  restraints_manager = mmtbx.restraints.manager(
    geometry      = geometry,
    normalization = True)
  if(xray_structure is not None):
    restraints_manager.crystal_symmetry = xray_structure.crystal_symmetry()
  grm = restraints_manager

  real_h.reset_i_seq_if_necessary()
  from mmtbx.geometry_restraints import reference
  grm.geometry.append_reference_coordinate_restraints_in_place(
      reference.add_coordinate_restraints(
          sites_cart = real_h.atoms().extract_xyz().select(helix_selection),
          selection  = helix_selection,
          sigma      = sigma_on_reference_helix))
  grm.geometry.append_reference_coordinate_restraints_in_place(
      reference.add_coordinate_restraints(
          sites_cart = real_h.atoms().extract_xyz().select(sheet_selection),
          selection  = sheet_selection,
          sigma      = sigma_on_reference_sheet))
  grm.geometry.append_reference_coordinate_restraints_in_place(
      reference.add_coordinate_restraints(
          sites_cart = real_h.atoms().extract_xyz().select(other_selection),
          selection  = other_selection,
          sigma      = sigma_on_reference_non_ss))
  grm.geometry.generic_restraints_manager.reference_manager.\
      add_torsion_restraints(
          pdb_hierarchy   = pre_result_h,
          sites_cart      = pre_result_h.atoms().extract_xyz().\
                                 select(ss_for_tors_selection),
          selection = ss_for_tors_selection,
          chi_angles_only = False,
          sigma           = sigma_on_torsion_ss)
  grm.geometry.generic_restraints_manager.reference_manager.\
      add_torsion_restraints(
          pdb_hierarchy   = pre_result_h,
          sites_cart      = real_h.atoms().extract_xyz().\
                                select(nonss_for_tors_selection),
          selection = nonss_for_tors_selection,
          chi_angles_only = False,
          sigma           = sigma_on_torsion_nonss)

  real_h.atoms().set_xyz(pre_result_h.atoms().extract_xyz())
  if fname_before_regularization is not None:
    real_h.write_pdb_file(file_name=fname_before_regularization)

  #testing number of restraints
  assert grm.geometry.get_n_den_proxies() == 0
  assert grm.geometry.get_n_reference_coordinate_proxies() == n_main_chain_atoms
  refinement_log = null_out()
  if verbose:
    refinement_log = log
    refinement_log.write(
      "Refining geometry of substituted secondary structure elements.\n")
  obj = run2(
      restraints_manager       = grm,
      pdb_hierarchy            = real_h,
      max_number_of_iterations = n_iter,
      number_of_macro_cycles   = n_macro,
      bond                     = True,
      nonbonded                = True,
      angle                    = True,
      dihedral                 = True,
      chirality                = True,
      planarity                = True,
      log                      = refinement_log)

  #print_hbond_proxies(grm.geometry,real_h)
  return grm.geometry.generic_restraints_manager.\
      reference_manager.reference_torsion_proxies


def beta():
  pdb_hierarchy = secondary_structure_from_sequence(beta_pdb_str,
      "ACEDGFIHKMLNQPSRTWVY")
  pdb_hierarchy.write_pdb_file(file_name = "o_beta_seq.pdb")

def alpha_310():
  pdb_hierarchy = secondary_structure_from_sequence(alpha310_pdb_str,
      "ACEDGFIHKMLNQPSRTWVY")
  pdb_hierarchy.write_pdb_file(file_name = "o_helix310_seq.pdb")

def alpha_pi():
  pdb_hierarchy = secondary_structure_from_sequence(alpha_pi_pdb_str,
      "ACEDGFIHKMLNQPSRTWVY")
  pdb_hierarchy.write_pdb_file(file_name = "o_helix_pi_seq.pdb")

def alpha():
  pdb_hierarchy = secondary_structure_from_sequence(alpha_pdb_str,
      "ACEDGFIHKMLNQPSRTWVY")
  pdb_hierarchy.write_pdb_file(file_name = "o_helix_seq.pdb")
