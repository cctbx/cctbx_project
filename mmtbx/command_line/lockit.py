from __future__ import division
import mmtbx.monomer_library.pdb_interpretation
import iotbx.pdb.atom_name_interpretation
import iotbx.mtz
import iotbx.phil, libtbx.phil
from cctbx import maptbx
import cctbx.maptbx.real_space_refinement_simple
import cctbx.geometry_restraints.flags
from cctbx.array_family import flex
import scitbx.rigid_body
import scitbx.graph.tardy_tree
import scitbx.lbfgs
import scitbx.math.superpose
from scitbx import matrix
from libtbx.str_utils import show_string
from libtbx.utils import Sorry, null_out
from libtbx import Auto, group_args
import libtbx
import sys, os
op = os.path

def real_space_rigid_body_gradients_simple(
      unit_cell,
      density_map,
      sites_cart_0,
      center_of_mass,
      q,
      unit_quaternion_delta=0.01,
      translation_delta=0.3):
  result = flex.double()
  q_delta = q.deep_copy()
  def get(i, delta):
    fs = []
    for signed_delta in [delta, -delta]:
      q_delta[i] = q[i] + signed_delta
      aja = matrix.rt(scitbx.rigid_body.joint_lib_six_dof_aja_simplified(
        center_of_mass=center_of_mass,
        q=q_delta))
      sites_cart_delta = aja * sites_cart_0
      rs_f = maptbx.real_space_target_simple(
        unit_cell=unit_cell,
        density_map=density_map,
        sites_cart=sites_cart_delta)
      fs.append(rs_f)
    result.append((fs[0]-fs[1])/(2*delta))
  for i in xrange(4): get(i=i, delta=unit_quaternion_delta)
  for i in xrange(3): get(i=i+4, delta=translation_delta)
  return result

class residue_refine_constrained(object):

  def __init__(O,
        pdb_hierarchy,
        residue,
        density_map,
        geometry_restraints_manager,
        real_space_target_weight,
        real_space_gradients_delta,
        lbfgs_termination_params):
    O.pdb_hierarchy = pdb_hierarchy
    O.residue = residue
    O.density_map = density_map
    O.geometry_restraints_manager = geometry_restraints_manager
    O.real_space_gradients_delta = real_space_gradients_delta
    O.real_space_target_weight = real_space_target_weight
    #
    O.unit_cell = geometry_restraints_manager.crystal_symmetry.unit_cell()
    O.sites_cart_all = pdb_hierarchy.atoms().extract_xyz()
    O.residue_i_seqs = residue.atoms().extract_i_seq()
    O.sites_cart_residue_0 = O.sites_cart_all.select(indices=O.residue_i_seqs)
    O.residue_center_of_mass = O.sites_cart_residue_0.mean()
    residue_tardy_tree = scitbx.graph.tardy_tree.construct(
      n_vertices=O.sites_cart_residue_0.size(),
      edge_list="all_in_one_rigid_body") \
        .build_tree() \
        .fix_near_singular_hinges(sites=None)
    O.residue_tardy_model = scitbx.rigid_body.tardy_model(
      labels=None,
      sites=O.sites_cart_residue_0,
      masses=flex.double(O.sites_cart_residue_0.size(), 1),
      tardy_tree=residue_tardy_tree,
      potential_obj=O)
    O.x = O.residue_tardy_model.pack_q()
    assert O.x.size() == 7 # other cases not implemented
    #
    O.number_of_function_evaluations = -1
    O.f_start, O.g_start = O.compute_functional_and_gradients()
    O.rs_f_start = O.rs_f
    O.minimizer = scitbx.lbfgs.run(
      target_evaluator=O,
      termination_params=lbfgs_termination_params)
    O.f_final, O.g_final = O.compute_functional_and_gradients()
    O.rs_f_final = O.rs_f
    del O.rs_f
    del O.x
    del O.residue_center_of_mass
    del O.sites_cart_residue_0
    del O.residue_i_seqs
    del O.sites_cart_all
    del O.unit_cell

  def compute_functional_and_gradients(O):
    if (O.number_of_function_evaluations == 0):
      O.number_of_function_evaluations += 1
      return O.f_start, O.g_start
    O.number_of_function_evaluations += 1
    O.residue_tardy_model.unpack_q(q_packed=O.x)
    O.sites_cart_residue = O.residue_tardy_model.sites_moved()
    rs_f = maptbx.real_space_target_simple(
      unit_cell=O.unit_cell,
      density_map=O.density_map,
      sites_cart=O.sites_cart_residue)
    rs_g = real_space_rigid_body_gradients_simple(
      unit_cell=O.unit_cell,
      density_map=O.density_map,
      sites_cart_0=O.sites_cart_residue_0,
      center_of_mass=O.residue_center_of_mass,
      q=O.x)
    O.rs_f = rs_f
    rs_f *= -O.real_space_target_weight
    rs_g *= -O.real_space_target_weight
    if (O.geometry_restraints_manager is None):
      f = rs_f
      g = rs_g
    else:
      O.sites_cart_all.set_selected(O.residue_i_seqs, O.sites_cart_residue)
      gr_e = O.geometry_restraints_manager.energies_sites(
        sites_cart=O.sites_cart_all, compute_gradients=True)
      O.__d_e_pot_d_sites = gr_e.gradients.select(indices=O.residue_i_seqs)
      f = rs_f + gr_e.target
      g = rs_g + O.residue_tardy_model.d_e_pot_d_q_packed()
    return f, g.as_double()

  def d_e_pot_d_sites(O, sites_moved):
    result = O.__d_e_pot_d_sites
    del O.__d_e_pot_d_sites
    return result

class residue_refine_restrained(object):

  def __init__(O,
        pdb_hierarchy,
        residue,
        density_map,
        geometry_restraints_manager,
        real_space_target_weight,
        real_space_gradients_delta,
        lbfgs_termination_params):
    O.pdb_hierarchy = pdb_hierarchy
    O.residue = residue
    O.density_map = density_map
    O.geometry_restraints_manager = geometry_restraints_manager
    O.real_space_gradients_delta = real_space_gradients_delta
    O.real_space_target_weight = real_space_target_weight
    #
    O.unit_cell = geometry_restraints_manager.crystal_symmetry.unit_cell()
    O.sites_cart_all = pdb_hierarchy.atoms().extract_xyz()
    O.residue_i_seqs = residue.atoms().extract_i_seq()
    O.x = O.sites_cart_all.select(indices=O.residue_i_seqs).as_double()
    #
    O.real_space_target = None
    O.number_of_function_evaluations = -1
    O.f_start, O.g_start = O.compute_functional_and_gradients()
    O.rs_f_start = O.rs_f
    O.minimizer = scitbx.lbfgs.run(
      target_evaluator=O,
      termination_params=lbfgs_termination_params)
    O.f_final, O.g_final = O.compute_functional_and_gradients()
    O.rs_f_final = O.rs_f
    del O.rs_f
    del O.x
    del O.residue_i_seqs
    del O.sites_cart_all
    del O.unit_cell

  def compute_functional_and_gradients(O):
    if (O.number_of_function_evaluations == 0):
      O.number_of_function_evaluations += 1
      return O.f_start, O.g_start
    O.number_of_function_evaluations += 1
    O.sites_cart_residue = flex.vec3_double(O.x)
    rs_f = maptbx.real_space_target_simple(
      unit_cell=O.unit_cell,
      density_map=O.density_map,
      sites_cart=O.sites_cart_residue)
    O.real_space_target = rs_f
    rs_g = maptbx.real_space_gradients_simple(
      unit_cell=O.unit_cell,
      density_map=O.density_map,
      sites_cart=O.sites_cart_residue,
      delta=O.real_space_gradients_delta)
    O.rs_f = rs_f
    rs_f *= -O.real_space_target_weight
    rs_g *= -O.real_space_target_weight
    if (O.geometry_restraints_manager is None):
      f = rs_f
      g = rs_g
    else:
      O.sites_cart_all.set_selected(O.residue_i_seqs, O.sites_cart_residue)
      gr_e = O.geometry_restraints_manager.energies_sites(
        sites_cart=O.sites_cart_all, compute_gradients=True)
      f = rs_f + gr_e.target
      g = rs_g + gr_e.gradients.select(indices=O.residue_i_seqs)
    return f, g.as_double()

def compute_functional_lite(pdb_hierarchy,
                            residue,
                            density_map,
                            geometry_restraints_manager,
                            real_space_target_weight):
  unit_cell = geometry_restraints_manager.crystal_symmetry.unit_cell()
  sites_cart_all = pdb_hierarchy.atoms().extract_xyz()
  residue_i_seqs = residue.atoms().extract_i_seq()
  x = sites_cart_all.select(indices=residue_i_seqs).as_double()
  sites_cart_residue = flex.vec3_double(x)
  rs_f = maptbx.real_space_target_simple(
    unit_cell=unit_cell,
    density_map=density_map,
    sites_cart=sites_cart_residue)
  return rs_f

def get_rotamer_iterator(mon_lib_srv, residue, atom_selection_bool):
  atoms = residue.atoms()
  if (atom_selection_bool is not None):
    if (atom_selection_bool.select(
          indices=residue.atoms().extract_i_seq()).all_eq(False)):
      return None
  rotamer_iterator = mon_lib_srv.rotamer_iterator(
    comp_id=residue.resname,
    atom_names=residue.atoms().extract_name(),
    sites_cart=residue.atoms().extract_xyz())
  if (rotamer_iterator.problem_message is not None):
    return None
  if (rotamer_iterator.rotamer_info is None):
    return None
  return rotamer_iterator

def rotamer_score_and_choose_best(
      mon_lib_srv,
      density_map,
      pdb_hierarchy,
      geometry_restraints_manager,
      atom_selection_bool,
      real_space_target_weight,
      real_space_gradients_delta,
      lbfgs_termination_params,
      force_rotamer=False):
  n_other_residues = 0
  n_amino_acids_ignored = 0
  n_amino_acids_scored = 0
  get_class = iotbx.pdb.common_residue_names_get_class
  def refine_constrained():
    refined = residue_refine_constrained(
      pdb_hierarchy=pdb_hierarchy,
      residue=residue,
      density_map=density_map,
      geometry_restraints_manager=geometry_restraints_manager,
      real_space_target_weight=real_space_target_weight,
      real_space_gradients_delta=real_space_gradients_delta,
      lbfgs_termination_params=lbfgs_termination_params)
    print residue.id_str(), "constr. refined(%s): %.6g -> %.6g" % (
      rotamer_id, refined.rs_f_start, refined.rs_f_final)
    return refined
  def refine_restrained():
    refined = residue_refine_restrained(
      pdb_hierarchy=pdb_hierarchy,
      residue=residue,
      density_map=density_map,
      geometry_restraints_manager=geometry_restraints_manager,
      real_space_target_weight=real_space_target_weight,
      real_space_gradients_delta=real_space_gradients_delta,
      lbfgs_termination_params=lbfgs_termination_params)
    print residue.id_str(), "restr. refined(%s): %.6g -> %.6g" % (
      rotamer_id, refined.rs_f_start, refined.rs_f_final)
    return refined
  def refine():
    residue.atoms().set_xyz(new_xyz=refine_constrained().sites_cart_residue)
    return refine_restrained()
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for residue in chain.only_conformer().residues():
        if (get_class(residue.resname) != "common_amino_acid"):
          n_other_residues += 1
        else:
          rotamer_iterator = get_rotamer_iterator(
            mon_lib_srv=mon_lib_srv,
            residue=residue,
            atom_selection_bool=atom_selection_bool)
          if (rotamer_iterator is None):
            n_amino_acids_ignored += 1
          else:
            best = None
            rotamer_id = "as_given"
            if not force_rotamer:
              best = group_args(rotamer_id=rotamer_id, refined=refine())
            n_amino_acids_scored += 1
            for rotamer,rotamer_sites_cart in rotamer_iterator:
              residue.atoms().set_xyz(new_xyz=rotamer_sites_cart)
              rotamer_id=rotamer.id
              trial = group_args(rotamer_id=rotamer.id, refined=refine())
              if best is None:
                best = trial
              if (trial.refined.rs_f_final > best.refined.rs_f_final):
                best = trial
            print residue.id_str(), "best rotamer:", best.rotamer_id
            residue.atoms().set_xyz(new_xyz=best.refined.sites_cart_residue)
            print

  print "number of amino acid residues scored:", n_amino_acids_scored
  print "number of amino acid residues ignored:", n_amino_acids_ignored
  print "number of other residues:", n_other_residues
  print
  sys.stdout.flush()

def get_best_rotamer(processed_pdb_file,
                     fft_map,
                     d_max,
                     d_min):
  master_phil = get_master_phil()
  work_params = master_phil.extract()
  get_class = iotbx.pdb.common_residue_names_get_class
  mon_lib_srv = mmtbx.monomer_library.server.server()
  grm = processed_pdb_file.geometry_restraints_manager(
    params_edits=work_params.geometry_restraints.edits,
    params_remove=work_params.geometry_restraints.remove)
  density_map = fft_map.real_map()
  atom_selection_bool=get_atom_selection_bool(
        scope_extract=work_params.rotamer_score_and_choose_best,
        attr="atom_selection",
        processed_pdb_file=processed_pdb_file,
        work_params=work_params)
  n_other_residues = 0
  n_amino_acids_ignored = 0
  n_amino_acids_scored = 0
  pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for residue in chain.only_conformer().residues():
        if (get_class(residue.resname) != "common_amino_acid"):
          n_other_residues += 1
        else:
          rotamer_iterator = get_rotamer_iterator(
            mon_lib_srv=mon_lib_srv,
            residue=residue,
            atom_selection_bool=atom_selection_bool)
          if (rotamer_iterator is None):
            n_amino_acids_ignored += 1
          else:
            best = None
            best_sites_cart = None
            n_amino_acids_scored += 1
            for rotamer,rotamer_sites_cart in rotamer_iterator:
              residue.atoms().set_xyz(new_xyz=rotamer_sites_cart)
              trial = group_args(rotamer_id=rotamer.id,
                                 rs_f=compute_functional_lite(
                pdb_hierarchy=pdb_hierarchy,
                residue=residue,
                density_map=density_map,
                geometry_restraints_manager=grm,
                real_space_target_weight=work_params.real_space_target_weight))
              #print residue.id_str(),trial.rotamer_id, trial.rs_f
              if best is None:
                best = trial
                best_sites_cart = rotamer_sites_cart
              if (trial.rs_f > best.rs_f):
                best = trial
                best_sites_cart = rotamer_sites_cart
            print residue.id_str(), "best rotamer:", best.rotamer_id
            residue.atoms().set_xyz(new_xyz=best_sites_cart)
            #print

real_space_params_phil_str = """\
real_space_target_weight = 10
  .type = float
real_space_gradients_delta_resolution_factor = 1/3
  .type = float
"""

coordinate_refinement_params_phil_str = """\
coordinate_refinement {
  run = False
    .type = bool
  atom_selection = Auto
    .type = str
  compute_final_correlation = True
    .type = bool
  finishing_geometry_minimization
  {
    cycles_max = 100
      .type = int
    first_weight_scale = 0.1
      .type = float
    cycle_weight_multiplier = 1.0
      .type = float
    superpose_cycle_end_with_cycle_start = False
      .type = bool
    dihedral_restraints = False
      .type = bool
    output_file = Auto
      .type = str
  }
  home_restraints
    .multiple = True
  {
    selection = None
      .type = str
    sigma = 0.05
      .type = float
    slack = 0
      .type = float
  }
  real_space_target_weights {
    first_sample = 10
      .type = float
    sampling_step = 30
      .type = float
    number_of_samples = 10
      .type = int
    bond_rmsd_target = 0.03
      .type = float
    worst_acceptable_bond_rmsd {
      pool_size = 10
        .type = int
      max_pool_average = 0.1
        .type = float
    }
  }
  lbfgs_max_iterations = 500
    .type = int
}
"""

def get_master_phil():
  return iotbx.phil.parse(
    input_string="""\
atom_selection = None
  .type = str
strict_processing = True
  .type = bool

map {
  coeff_labels {
    f = 2FOFCWT,PH2FOFCWT
      .type = str
    phases = None
      .type = str
    weights = None
      .type = str
  }
  low_resolution = None
    .type = float
  high_resolution = None
    .type = float
  grid_resolution_factor = 1/3
    .type = float
}

output_file = Auto
  .type = str
final_geo_file = Auto
  .type = str

%s

%s

rotamer_score_and_choose_best {
  run = False
    .type = bool
  atom_selection = Auto
    .type = str
  lbfgs_max_iterations = 50
    .type = int
  output_file = Auto
    .type = str
}

include scope mmtbx.monomer_library.pdb_interpretation.grand_master_phil_str
""" % (real_space_params_phil_str, coordinate_refinement_params_phil_str),
    process_includes=True)

def extract_map_coeffs(params, miller_arrays):
  def find(labels):
    for miller_array in miller_arrays:
      if (",".join(miller_array.info().labels) == labels):
        return miller_array
    matching_array = None
    for miller_array in miller_arrays:
      if (",".join(miller_array.info().labels).lower() == labels.lower()):
        if (matching_array is not None):
          return None
        matching_array = miller_array
    return matching_array
  def raise_sorry(msg_intro, name):
    msg = [
      msg_intro,
      "  %s = %s" % params.__phil_path_and_value__(object_name=name),
      "  List of available labels:"]
    for miller_array in miller_arrays:
      msg.append("    %s" % ",".join(miller_array.info().labels))
    raise Sorry("\n".join(msg))
  if (params.f is None):
    raise_sorry(msg_intro="Missing assignment:", name="f")
  f = find(labels=params.f)
  if (f is None):
    raise_sorry(msg_intro="Cannot find map coefficients:", name="f")
  if (not f.is_complex_array() or params.phases is not None):
    if (params.phases is None):
      raise_sorry(msg_intro="Missing assignment:", name="phases")
    phases = find(labels=params.phases)
    if (phases is None):
      raise_sorry(
        msg_intro="Cannot find map coefficient phases:", name="phases")
    if (phases.is_hendrickson_lattman_array()
          and phases.data().all_eq((0,0,0,0))):
      raise Sorry(
        "All Hendrickson-Lattman coefficient zero:\n"
        "  %s = %s" % params.__phil_path_and_value__(object_name="phases"))
    cf, cp = f.common_sets(other=phases)
    if (cf.indices().size() != f.indices().size()):
      raise Sorry(
        "Number of missing map coefficient phases: %d" % (
          f.indices().size() - cf.indices().size()))
    f = cf.phase_transfer(phase_source=cp, deg=True)
  if (params.weights is not None):
    weights = find(labels=params.weights)
    if (weights is None):
      raise_sorry(
        msg_intro="Cannot find map coefficient weights:",
        name="weights")
    cf, cw = f.common_sets(other=weights)
    if (cf.indices().size() != f.indices().size()):
      raise Sorry(
        "Number of missing map coefficient weights: %d" % (
          f.indices().size() - cf.indices().size()))
    f = cf.customized_copy(data=cw.data()*cf.data())
  return f

class home_restraints(object):

  __slots__ = ["iselection", "weight", "slack"]

  def __init__(O, iselection, weight, slack):
    O.iselection = iselection
    O.weight = weight
    O.slack = slack

def process_home_restraints_params(work_params, processed_pdb_file):
  result = []
  for params in work_params:
    sigma = params.sigma
    if (sigma is None or sigma <= 0):
      continue
    bsel = processed_pdb_file.all_chain_proxies.phil_atom_selection(
      cache=None,
      scope_extract=params,
      attr="selection",
      allow_none=True,
      allow_auto=False)
    if (bsel is not None):
      slack = params.slack
      if (slack is None or slack <= 0):
        slack = 0
      result.append(home_restraints(
        iselection=bsel.iselection(), weight=1/sigma**2, slack=slack))
  return result

class geometry_restraints_manager_plus(object):

  __slots__ = [
    "manager",
    "home_sites_cart",
    "home_restraints_list",
    "home_dihedral_proxies",
    "crystal_symmetry",
    "pair_proxies",
    "angle_proxies",
    "dihedral_proxies",
    "site_symmetry_table",
    "plain_pair_sym_table",
    "plain_pairs_radius"]

  def __init__(O,
               manager,
               home_sites_cart=None,
               home_restraints_list=None,
               home_dihedral_proxies=None):
    O.manager = manager
    #O.sites_cart = sites_cart
    O.home_sites_cart = home_sites_cart
    O.home_restraints_list = home_restraints_list
    O.home_dihedral_proxies = home_dihedral_proxies
    O.crystal_symmetry = manager.crystal_symmetry
    O.pair_proxies = manager.pair_proxies
    O.angle_proxies = manager.angle_proxies
    O.dihedral_proxies = manager.dihedral_proxies
    O.site_symmetry_table = manager.site_symmetry_table
    O.plain_pair_sym_table = manager.plain_pair_sym_table
    O.plain_pairs_radius = manager.plain_pairs_radius

  def energies_add(O, energies_obj):
    if (O.manager.site_symmetry_table is not None):
      site_symmetry_table_indices = O.manager.site_symmetry_table.indices()
    else:
      site_symmetry_table_indices = None
    from cctbx.geometry_restraints import bond
    n_restraints = 0
    r_sum = 0
    hr_summation = cctbx.geometry_restraints \
      .home_restraints_summation_skip_special_positions
    for hr in O.home_restraints_list:
      r_sum += hr_summation(
        sites_cart=energies_obj.sites_cart,
        gradients=energies_obj.gradients,
        site_symmetry_table_indices=site_symmetry_table_indices,
        home_sites_cart=O.home_sites_cart,
        iselection=hr.iselection,
        weight=hr.weight,
        slack=hr.slack)
      n_restraints += hr.iselection.size()
    energies_obj.n_home_restraints = n_restraints
    energies_obj.home_restraints_residual_sum = r_sum
    energies_obj.number_of_restraints += n_restraints
    energies_obj.residual_sum += r_sum
    r_sum=cctbx.geometry_restraints.dihedral_residual_sum(
      sites_cart=energies_obj.sites_cart,
      proxies=O.home_dihedral_proxies,
      gradient_array=energies_obj.gradients)
    if O.home_dihedral_proxies is not None:
      n_restraints = len(O.home_dihedral_proxies)
    else:
      n_restraints = 0
    energies_obj.n_dihedral_restraints = n_restraints
    energies_obj.dihedral_restraints_residual_sum = r_sum
    energies_obj.residual_sum += r_sum
    energies_obj.number_of_restraints += n_restraints

  def energies_show(O, energies_obj, f, prefix):
    print >> f, prefix+"  home_restraints_residual_sum (n=%d): %.6g" % (
      energies_obj.n_home_restraints,
      energies_obj.home_restraints_residual_sum)

  def dihedral_energies_show(O, energies_obj, f, prefix):
    print >> f, prefix+"  dihedral_restraints_residual_sum (n=%d): %.6g" % (
      energies_obj.n_dihedral_restraints,
      energies_obj.dihedral_restraints_residual_sum)

  def energies_sites(O,
        sites_cart,
        flags=None,
        custom_nonbonded_function=None,
        compute_gradients=False,
        gradients=None,
        disable_asu_cache=False,
        normalization=False):
    return O.manager.energies_sites(
      sites_cart=sites_cart,
      flags=flags,
      custom_nonbonded_function=custom_nonbonded_function,
      compute_gradients=compute_gradients,
      gradients=gradients,
      disable_asu_cache=disable_asu_cache,
      normalization=normalization,
      extension_objects=[O])

  def show_sorted(O,
        flags=None,
        sites_cart=None,
        site_labels=None,
        f=None):
    return O.manager.show_sorted(
      flags=flags,
      sites_cart=sites_cart,
      site_labels=site_labels,
      f=f)

  def select(O, selection=None, iselection=None):
    return O.manager.select(selection=selection, iselection=iselection)

def fit_rotamers_simple(processed_pdb_file,
                        fft_map,
                        d_max,
                        d_min):
  master_phil = get_master_phil()
  work_params = master_phil.extract()
  mon_lib_srv = mmtbx.monomer_library.server.server()
  #ener_lib = mmtbx.monomer_library.server.ener_lib()
  grm = processed_pdb_file.geometry_restraints_manager(
    params_edits=work_params.geometry_restraints.edits,
    params_remove=work_params.geometry_restraints.remove)
  density_map = fft_map.real_map()
  real_space_gradients_delta = \
    d_min * work_params.real_space_gradients_delta_resolution_factor
  rotamer_score_and_choose_best(
      mon_lib_srv=mon_lib_srv,
      density_map=density_map,
      pdb_hierarchy=processed_pdb_file.all_chain_proxies.pdb_hierarchy,
      geometry_restraints_manager=grm,
      atom_selection_bool=get_atom_selection_bool(
        scope_extract=work_params.rotamer_score_and_choose_best,
        attr="atom_selection",
        processed_pdb_file=processed_pdb_file,
        work_params=work_params),
      real_space_target_weight=work_params.real_space_target_weight,
      real_space_gradients_delta=real_space_gradients_delta,
      lbfgs_termination_params=scitbx.lbfgs.termination_parameters(
        max_iterations=work_params
          .rotamer_score_and_choose_best.lbfgs_max_iterations),
      force_rotamer=True)

def get_atom_selection_bool(scope_extract,
                            attr,
                            processed_pdb_file,
                            work_params):
  common_atom_selection_bool_cache=[]
  result = processed_pdb_file.all_chain_proxies \
    .phil_atom_selection(
      cache=None,
      scope_extract=scope_extract,
      attr="atom_selection",
      allow_none=True,
      allow_auto=True)
  if (result is None or result is not Auto):
    return result
  if (len(common_atom_selection_bool_cache) == 0):
    common_atom_selection_bool_cache.append(
      processed_pdb_file.all_chain_proxies
        .phil_atom_selection(
          cache=None,
          scope_extract=work_params,
          attr="atom_selection",
          allow_none=True))
  return common_atom_selection_bool_cache[0]

def compose_output_file_name(input_pdb_file_name, new_suffix):
  result = op.basename(input_pdb_file_name)
  if (   result.endswith(".pdb")
      or result.endswith(".ent")):
    result = result[:-4]
  result += new_suffix
  return result

def write_pdb(
      file_name,
      input_pdb_file_name,
      processed_pdb_file,
      grm,
      new_suffix):
  if (file_name is not None):
    if (file_name is Auto):
      file_name = compose_output_file_name(
        input_pdb_file_name=input_pdb_file_name,
        new_suffix=new_suffix)
    pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
    print "Writing file: %s" % show_string(file_name)
    sys.stdout.flush()
    pdb_hierarchy.write_pdb_file(
      file_name=file_name,
      crystal_symmetry=grm.crystal_symmetry)
    print

def run_coordinate_refinement_driver(
      processed_pdb_file,
      geometry_restraints_manager,
      fft_map,
      real_space_gradients_delta,
      work_params,
      write_pdb_callback=None,
      log=None):
  if (log is None): log = null_out()
  atom_selection_bool = get_atom_selection_bool(
    scope_extract=work_params.coordinate_refinement,
    attr="atom_selection",
    processed_pdb_file=processed_pdb_file,
    work_params=work_params)
  home_restraints_list = process_home_restraints_params(
    work_params=work_params.coordinate_refinement.home_restraints,
    processed_pdb_file=processed_pdb_file)
  if (not work_params.coordinate_refinement.compute_final_correlation):
    work_scatterers = None
  else:
    work_scatterers = processed_pdb_file.xray_structure(
      show_summary=False).scatterers()
    if (atom_selection_bool is None):
      work_scatterers = work_scatterers.deep_copy()
    else:
      work_scatterers = work_scatterers.select(atom_selection_bool)
  return run_coordinate_refinement(
    pdb_atoms=processed_pdb_file.all_chain_proxies.pdb_atoms,
    geometry_restraints_manager=geometry_restraints_manager,
    selection_variable=atom_selection_bool,
    density_map=fft_map.real_map(),
    real_space_gradients_delta=real_space_gradients_delta,
    work_params=work_params,
    home_restraints_list=home_restraints_list,
    work_scatterers=work_scatterers,
    unit_cell=fft_map.unit_cell(),
    d_min=fft_map.d_min(),
    write_pdb_callback=write_pdb_callback,
    log=log)

def run_coordinate_refinement(
      pdb_atoms,
      geometry_restraints_manager,
      selection_variable,
      density_map,
      real_space_gradients_delta,
      work_params,
      home_restraints_list=[],
      work_scatterers=None,
      unit_cell=None,
      d_min=None,
      write_pdb_callback=None,
      log=None):
  if (work_scatterers is not None):
    assert unit_cell is not None
    assert d_min is not None
  best_info = None
  sites_cart_start = pdb_atoms.extract_xyz()
  site_labels = [atom.id_str() for atom in pdb_atoms]
  grmp = geometry_restraints_manager_plus(
    manager=geometry_restraints_manager,
    home_sites_cart=sites_cart_start,
    home_restraints_list=home_restraints_list)
  print >> log, "Before coordinate refinement:"
  grmp.energies_sites(sites_cart=sites_cart_start).show(f=log)
  print >> log
  log.flush()
  rstw_params = work_params.coordinate_refinement.real_space_target_weights
  if (rstw_params.number_of_samples is None):
    rstw_list = [work_params.real_space_target_weight]
  else:
    rstw_list = [rstw_params.first_sample + i * rstw_params.sampling_step
      for i in xrange(rstw_params.number_of_samples)]
  lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
    max_iterations=work_params.coordinate_refinement
      .lbfgs_max_iterations)
  lbfgs_exception_handling_params = \
    scitbx.lbfgs.exception_handling_parameters(
      ignore_line_search_failed_step_at_lower_bound=True)
  bond_rmsd_list = []
  for rstw in rstw_list:
    refined = maptbx.real_space_refinement_simple.lbfgs(
      sites_cart=sites_cart_start,
      density_map=density_map,
      selection_variable=selection_variable,
      geometry_restraints_manager=grmp,
      real_space_target_weight=rstw,
      real_space_gradients_delta=real_space_gradients_delta,
      lbfgs_termination_params=lbfgs_termination_params,
      lbfgs_exception_handling_params=lbfgs_exception_handling_params)
    print >> log, "After coordinate refinement" \
      " with real-space target weight %.1f:" % rstw
    grmp.energies_sites(sites_cart=refined.sites_cart).show(f=log)
    bond_proxies = grmp.pair_proxies().bond_proxies
    bond_proxies.show_sorted(
      by_value="residual",
      sites_cart=refined.sites_cart,
      site_labels=site_labels,
      f=log,
      prefix="  ",
      max_items=3)
    print >> log
    print >> log, "  number_of_function_evaluations:", \
      refined.number_of_function_evaluations
    print >> log, "  real+geo target start: %.6g" % refined.f_start
    print >> log, "  real+geo target final: %.6g" % refined.f_final
    deltas = bond_proxies.deltas(sites_cart=refined.sites_cart)
    deltas_abs = flex.abs(deltas)
    deltas_abs_sorted = deltas_abs.select(
      flex.sort_permutation(deltas_abs), reverse=True)
    wabr_params = rstw_params.worst_acceptable_bond_rmsd
    acceptable = (
      flex.mean(deltas_abs_sorted[:wabr_params.pool_size])
        < wabr_params.max_pool_average)
    bond_rmsd = flex.mean_sq(deltas)**0.5
    bond_rmsd_list.append(bond_rmsd)
    print >> log, "  Bond RMSD: %.3f" % bond_rmsd
    #
    if (work_scatterers is None):
      region_cc = None
    else:
      region_cc = maptbx.region_density_correlation(
        large_unit_cell=unit_cell,
        large_d_min=d_min,
        large_density_map=density_map,
        sites_cart=refined.sites_cart_variable,
        site_radii=flex.double(refined.sites_cart_variable.size(), 1),
        work_scatterers=work_scatterers)
    #
    rmsd_diff = abs(rstw_params.bond_rmsd_target - bond_rmsd)
    if (   best_info is None
        or best_info.rmsd_diff > rmsd_diff
          and (acceptable or not best_info.acceptable)):
      best_info = group_args(
        rstw=rstw,
        refined=refined,
        acceptable=acceptable,
        rmsd_diff=rmsd_diff,
        region_cc=region_cc)
    print >> log
  if (best_info is not None):
    print >> log, "Table of real-space target weights vs. bond RMSD:"
    print >> log, "  weight   RMSD"
    for w,d in zip(rstw_list, bond_rmsd_list):
      print >> log, "  %6.1f  %5.3f" % (w,d)
    print >> log, "Best real-space target weight: %.1f" % best_info.rstw
    print >> log, \
      "Associated refined final value: %.6g" % best_info.refined.f_final
    if (best_info.region_cc is not None):
      print >> log, "Associated region correlation: %.4f" % best_info.region_cc
    print >> log
    pdb_atoms.set_xyz(new_xyz=refined.sites_cart)
    if (write_pdb_callback is not None):
      write_pdb_callback(situation="after_best_weight_determination")
    #
    fgm_params = work_params.coordinate_refinement \
      .finishing_geometry_minimization
    if (fgm_params.cycles_max is not None and fgm_params.cycles_max > 0):
      print >> log, "As previously obtained with target weight %.1f:" \
        % best_info.rstw
      grmp.energies_sites(sites_cart=refined.sites_cart).show(
        f=log, prefix="  ")
      print >> log, "Finishing refinement to idealize geometry:"
      print >> log, "            number of function"
      print >> log, "    weight     evaluations      cycle RMSD"
      number_of_fgm_cycles = 0
      rstw = best_info.rstw * fgm_params.first_weight_scale
      sites_cart_start = best_info.refined.sites_cart.deep_copy()
      for i_cycle in xrange(fgm_params.cycles_max):
        fgm_refined = maptbx.real_space_refinement_simple.lbfgs(
          sites_cart=sites_cart_start,
          density_map=density_map,
          selection_variable=selection_variable,
          geometry_restraints_manager=grmp,
          energies_sites_flags=cctbx.geometry_restraints.flags.flags(
            default=True, dihedral=fgm_params.dihedral_restraints),
          real_space_target_weight=rstw,
          real_space_gradients_delta=real_space_gradients_delta,
          lbfgs_termination_params=lbfgs_termination_params,
          lbfgs_exception_handling_params=lbfgs_exception_handling_params)
        cycle_rmsd = sites_cart_start.rms_difference(fgm_refined.sites_cart)
        print >> log, "   %6.1f     %10d          %6.3f" % (
          rstw, fgm_refined.number_of_function_evaluations, cycle_rmsd)
        number_of_fgm_cycles += 1
        rstw *= fgm_params.cycle_weight_multiplier
        if (cycle_rmsd < 1.e-4):
          break
        if (fgm_params.superpose_cycle_end_with_cycle_start):
          fit = scitbx.math.superpose.least_squares_fit(
            reference_sites=best_info.refined.sites_cart,
            other_sites=fgm_refined.sites_cart)
          fgm_refined.sites_cart = fit.other_sites_best_fit()
        sites_cart_start = fgm_refined.sites_cart.deep_copy()
      print >> log, "After %d refinements to idealize geometry:" % (
        number_of_fgm_cycles)
      grmp.energies_sites(sites_cart=fgm_refined.sites_cart).show(
        f=log, prefix="  ")
      if (work_scatterers is None):
        fgm_region_cc = None
      else:
        fgm_region_cc = maptbx.region_density_correlation(
          large_unit_cell=unit_cell,
          large_d_min=d_min,
          large_density_map=density_map,
          sites_cart=fgm_refined.sites_cart_variable,
          site_radii=flex.double(fgm_refined.sites_cart_variable.size(), 1),
          work_scatterers=work_scatterers)
        print >> log, "  Associated region correlation: %.4f" % fgm_region_cc
      print >> log
      best_info.fgm_refined = fgm_refined
      best_info.fgm_region_cc = fgm_region_cc
      pdb_atoms.set_xyz(new_xyz=fgm_refined.sites_cart)
      if (write_pdb_callback is not None):
        write_pdb_callback(situation="after_finishing_geo_min")
    else:
      best_info.fgm_refined = None
      best_info.fgm_region_cc = None
  return best_info

def run(args):
  show_times = libtbx.utils.show_times(time_start="now")
  master_phil = get_master_phil()
  import iotbx.utils
  input_objects = iotbx.utils.process_command_line_inputs(
    args=args,
    master_phil=master_phil,
    input_types=("mtz", "pdb", "cif"))
  work_phil = master_phil.fetch(sources=input_objects["phil"])
  work_phil.show()
  print
  print "#phil __OFF__"
  print
  work_params = work_phil.extract()
  #
  assert len(input_objects["mtz"]) == 1
  map_coeffs = extract_map_coeffs(
    miller_arrays=input_objects["mtz"][0].file_content.as_miller_arrays(),
    params=work_params.map.coeff_labels)
  #
  mon_lib_srv = mmtbx.monomer_library.server.server()
  ener_lib = mmtbx.monomer_library.server.ener_lib()
  for file_obj in input_objects["cif"]:
    print "Processing CIF file: %s" % show_string(file_obj.file_name)
    for srv in [mon_lib_srv, ener_lib]:
      srv.process_cif_object(
        cif_object=file_obj.file_content,
        file_name=file_obj.file_name)
  if (len(input_objects["cif"]) != 0):
    print
  #
  assert len(input_objects["pdb"]) == 1 # TODO not implemented
  file_obj = input_objects["pdb"][0]
  input_pdb_file_name = file_obj.file_name
  print "Crystal symmetry (from map file):"
  map_coeffs.crystal_symmetry().show_summary(prefix="  ")
  pdb_crystal_symmetry = file_obj.file_content.crystal_symmetry()
  if (    pdb_crystal_symmetry is not None
      and not pdb_crystal_symmetry.is_similar_symmetry(
            map_coeffs.crystal_symmetry())):
    print "  NOTE: Crystal symmetry from PDB file is different:"
    pdb_crystal_symmetry.show_summary(prefix="    Ignored: ")
  print
  processed_pdb_file = mmtbx.monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    params=work_params.pdb_interpretation,
    file_name=file_obj.file_name,
    pdb_inp=file_obj.file_content,
    strict_conflict_handling=work_params.strict_processing,
    substitute_non_crystallographic_unit_cell_if_necessary=True,
    crystal_symmetry=map_coeffs.crystal_symmetry(),
    force_symmetry=True,
    log=sys.stdout)
  if (work_params.strict_processing):
    msg = processed_pdb_file.all_chain_proxies.fatal_problems_message()
    if (msg is not None):
      raise Sorry(msg)
  #
  grm = processed_pdb_file.geometry_restraints_manager(
    params_edits=work_params.geometry_restraints.edits,
    params_remove=work_params.geometry_restraints.remove)
  print
  sys.stdout.flush()
  #
  def show_completeness(annotation):
    print "Completeness of %s map coefficients:" % annotation
    map_coeffs.setup_binner(auto_binning=True)
    if (map_coeffs.binner().n_bins_used() > 12):
      map_coeffs.setup_binner(n_bins=12)
    map_coeffs.completeness(use_binning=True).show(prefix="  ")
    print
    sys.stdout.flush()
  show_completeness("input")
  map_coeffs_input = map_coeffs
  #
  low_res = work_params.map.low_resolution
  high_res = work_params.map.high_resolution
  d_max, d_min = map_coeffs.d_max_min()
  d_max_apply, d_min_apply = None, None
  if (low_res is not None and low_res < d_max):
    d_max_apply = low_res
    print "Applying low resolution cutoff to map coefficients:" \
      " d_max=%.6g" % d_max_apply
  if (high_res is not None and high_res > d_min):
    d_min_apply = high_res
    print "Applying high resolution cutoff to map coefficients:" \
      " d_min=%.6g" % d_min_apply
  if (d_max_apply is not None or d_min_apply is not None):
    map_coeffs = map_coeffs.resolution_filter(
      d_max=d_max_apply, d_min=d_min_apply)
    if (d_min_apply is not None):
      d_min = d_min_apply
    print
    sys.stdout.flush()
  #
  if (map_coeffs is not map_coeffs_input):
    show_completeness("final")
  #
  fft_map = map_coeffs.fft_map(
    d_min=d_min,
    resolution_factor=work_params.map.grid_resolution_factor)
  fft_map.apply_sigma_scaling()
  #
  real_space_gradients_delta = \
    d_min * work_params.real_space_gradients_delta_resolution_factor
  print "real_space_gradients_delta: %.6g" % real_space_gradients_delta
  print
  sys.stdout.flush()
  #
  if (work_params.rotamer_score_and_choose_best.run):
    rotamer_score_and_choose_best(
      mon_lib_srv=mon_lib_srv,
      density_map=fft_map.real_map(),
      pdb_hierarchy=processed_pdb_file.all_chain_proxies.pdb_hierarchy,
      geometry_restraints_manager=grm,
      atom_selection_bool=get_atom_selection_bool(
        scope_extract=work_params.rotamer_score_and_choose_best,
        attr="atom_selection",
        processed_pdb_file=processed_pdb_file,
        work_params=work_params),
      real_space_target_weight=work_params.real_space_target_weight,
      real_space_gradients_delta=real_space_gradients_delta,
      lbfgs_termination_params=scitbx.lbfgs.termination_parameters(
        max_iterations=work_params
          .rotamer_score_and_choose_best.lbfgs_max_iterations))
    write_pdb(
      file_name=work_params.rotamer_score_and_choose_best.output_file,
      input_pdb_file_name=input_pdb_file_name,
      processed_pdb_file=processed_pdb_file,
      grm=grm,
      new_suffix="_lockit_rotamer_score_and_choose_best.pdb")
  #
  best_info = None
  if (work_params.coordinate_refinement.run):
    def write_pdb_callback(situation):
      if (situation == "after_best_weight_determination"):
        write_pdb(
          file_name=work_params.output_file,
          input_pdb_file_name=input_pdb_file_name,
          processed_pdb_file=processed_pdb_file,
          grm=grm,
          new_suffix="_lockit_best_weight.pdb")
      elif (situation == "after_finishing_geo_min"):
        write_pdb(
          file_name=work_params.coordinate_refinement \
            .finishing_geometry_minimization.output_file,
          input_pdb_file_name=input_pdb_file_name,
          processed_pdb_file=processed_pdb_file,
          grm=grm,
          new_suffix="_lockit_finishing_geo_min.pdb")
      else:
        raise AssertionError
    best_info = run_coordinate_refinement_driver(
      processed_pdb_file=processed_pdb_file,
      geometry_restraints_manager=grm,
      fft_map=fft_map,
      real_space_gradients_delta=real_space_gradients_delta,
      work_params=work_params,
      write_pdb_callback=write_pdb_callback,
      log=sys.stdout)
  else:
    write_pdb(
      file_name=work_params.output_file,
      input_pdb_file_name=input_pdb_file_name,
      processed_pdb_file=processed_pdb_file,
      grm=grm,
      new_suffix="_lockit_best_weight.pdb")
  #
  file_name = work_params.final_geo_file
  if (file_name is not None):
    if (file_name is Auto):
      file_name = compose_output_file_name(
        input_pdb_file_name=input_pdb_file_name,
        new_suffix="_lockit_final.geo")
    pdb_atoms = processed_pdb_file.all_chain_proxies.pdb_atoms
    sites_cart = pdb_atoms.extract_xyz()
    site_labels = [atom.id_str() for atom in pdb_atoms]
    print "Writing file: %s" % show_string(file_name)
    sys.stdout.flush()
    f = open(file_name, "w")
    grm.show_sorted(sites_cart=sites_cart, site_labels=site_labels, f=f)
    del f
    print
  #
  show_times()
  sys.stdout.flush()
  return best_info

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
