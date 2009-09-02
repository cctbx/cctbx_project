from __future__ import division
import mmtbx.monomer_library.pdb_interpretation
import iotbx.pdb.atom_name_interpretation
import iotbx.mtz
import iotbx.phil, libtbx.phil
from cctbx import maptbx
import cctbx.maptbx.real_space_refinement_simple
import cctbx.geometry_restraints
from cctbx.array_family import flex
import scitbx.rigid_body
import scitbx.graph.tardy_tree
import scitbx.lbfgs
from scitbx import matrix
from libtbx.str_utils import show_string
from libtbx.utils import Sorry
from libtbx import group_args
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
      lbfgs_termination_params):
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
            rotamer_id = "as_given"
            best = group_args(rotamer_id=rotamer_id, refined=refine())
            n_amino_acids_scored += 1
            for rotamer,rotamer_sites_cart in rotamer_iterator:
              residue.atoms().set_xyz(new_xyz=rotamer_sites_cart)
              trial = group_args(rotamer_id=rotamer.id, refined=refine())
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

class try_read_file(object):

  def __init__(O, file_name):
    def set(file_type, file_content):
      O.file_name = file_name
      O.file_type = file_type
      O.file_content = file_content
    lead = open(file_name, "rb").read(3)
    if (lead == "MTZ"):
      mtz_obj = iotbx.mtz.object(file_name=file_name)
      set(file_type="mtz", file_content=mtz_obj)
      return
    try:
      pdb_inp = iotbx.pdb.input(file_name=file_name)
    except KeyboardInterrupt: raise
    except:
      if (iotbx.pdb.is_pdb_file(file_name=file_name)):
        raise
      pdb_inp = None
    else:
      if (pdb_inp.atoms().size() != 0):
        set(file_type="pdb", file_content=pdb_inp)
        return
    try:
      cif_obj = mmtbx.monomer_library.server.read_cif(file_name=file_name)
    except KeyboardInterrupt: raise
    except: pass
    else:
      if (len(cif_obj) != 0):
        set(file_type="cif", file_content=cif_obj)
        return
    try:
      phil_obj = iotbx.phil.parse(file_name=file_name)
    except KeyboardInterrupt: raise
    except: pass
    else:
      set(file_type="phil", file_content=phil_obj)
      return
    if (pdb_inp is not None):
      if (pdb_inp.unknown_section().size() != 0):
        set(file_type=None, file_content=None)
        return
      if (pdb_inp.header_section().size() != 0):
        set(file_type="pdb", file_content=pdb_inp)
    set(file_type=None, file_content=None)
    return

def get_master_phil():
  return iotbx.phil.parse(
    input_string="""\
atom_selection = None
  .type = str

symmetry_from_file = None
  .type = path
  .multiple = True
unit_cell = None
  .type=unit_cell
space_group = None
  .type=space_group

map_coeff_labels = 2FOFCWT PH2FOFCWT
  .type = strings
map_resolution_factor = 1/3
  .type = float

real_space_target_weight = 1
  .type = float
real_space_gradients_delta_resolution_factor = 1/3
  .type = float

all_coordinate_refinement {
  run = False
    .type = bool
  lbfgs_max_iterations = 500
    .type = int
}

rotamer_score_and_choose_best {
  run = False
    .type = bool
  lbfgs_max_iterations = 50
    .type = int
}

pdb_interpretation {
  include scope mmtbx.monomer_library.pdb_interpretation.master_params
}

geometry_restraints.edits {
  include scope \
    mmtbx.monomer_library.pdb_interpretation.geometry_restraints_edits_str
}

geometry_restraints.remove {
  include scope \
    mmtbx.monomer_library.pdb_interpretation.geometry_restraints_remove_str
}
""", process_includes=True)

def run(args):
  show_times = libtbx.utils.show_times(time_start="now")
  master_phil = get_master_phil()
  argument_interpreter = libtbx.phil.command_line.argument_interpreter(
    master_phil=master_phil)
  phil_objects = []
  file_objects = {
    "mtz": [],
    "pdb": [],
    "cif": []}
  for arg in args:
    if (len(arg) == 0): continue
    def try_as_file():
      if (not op.isfile(arg)): return False
      obj = try_read_file(file_name=arg)
      if (obj.file_type is None): return False
      if (obj.file_type == "phil"):
        phil_objects.append(obj.file_content)
      else:
        file_objects[obj.file_type].append(obj)
      return True
    def try_as_command_line_params():
      try: command_line_params = argument_interpreter.process(arg=arg)
      except KeyboardInterrupt: raise
      except:
        if (op.isfile(arg)):
          raise Sorry(
            "Error processing file: %s" % show_string(arg))
        raise Sorry(
          "Command-line argument not recognized: %s" % show_string(arg))
      phil_objects.append(command_line_params)
    if (not try_as_file()):
      try_as_command_line_params()
  work_phil = master_phil.fetch(sources=phil_objects)
  work_phil.show()
  print
  work_params = work_phil.extract()
  #
  assert len(work_params.symmetry_from_file) == 0 # TODO not implemented
  assert work_params.unit_cell is None # TODO not implemented
  assert work_params.space_group is None # TODO not implemented
  #
  assert len(file_objects["mtz"]) == 1
  miller_arrays = file_objects["mtz"][0].file_content.as_miller_arrays()
  map_coeffs = None
  for miller_array in miller_arrays:
    if (miller_array.info().labels == work_params.map_coeff_labels):
      map_coeffs = miller_array
      break
  assert map_coeffs is not None
  #
  mon_lib_srv = mmtbx.monomer_library.server.server()
  ener_lib = mmtbx.monomer_library.server.ener_lib()
  for file_obj in file_objects["cif"]:
    print "Processing CIF file: %s" % show_string(file_obj.file_name)
    for srv in [mon_lib_srv, ener_lib]:
      srv.process_cif_object(
        cif_object=file_obj.file_content,
        file_name=file_obj.file_name)
  #
  assert len(file_objects["pdb"]) == 1 # TODO not implemented
  file_obj = file_objects["pdb"][0]
  input_pdb_file_name = file_obj.file_name
  processed_pdb_file = mmtbx.monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    params=work_params.pdb_interpretation,
    file_name=file_obj.file_name,
    pdb_inp=file_obj.file_content,
    strict_conflict_handling=True,
    substitute_non_crystallographic_unit_cell_if_necessary=True,
    log=sys.stdout)
  #
  grm = processed_pdb_file.geometry_restraints_manager(
    params_edits=work_params.geometry_restraints.edits,
    params_remove=work_params.geometry_restraints.remove)
  print
  sys.stdout.flush()
  #
  d_min = map_coeffs.d_min()
  fft_map = map_coeffs.fft_map(
    d_min=d_min,
    resolution_factor=work_params.map_resolution_factor)
  fft_map.apply_sigma_scaling()
  density_map = fft_map.real_map()
  real_space_gradients_delta = \
    d_min * work_params.real_space_gradients_delta_resolution_factor
  print "real_space_gradients_delta: %.6g" % real_space_gradients_delta
  print
  sys.stdout.flush()
  #
  if (work_params.rotamer_score_and_choose_best.run):
    if (work_params.atom_selection is None):
      atom_selection_bool = None
    else:
      atom_selection_bool = processed_pdb_file.all_chain_proxies \
        .phil_atom_selection(
          cache=None,
          scope_extract=work_params,
          attr="atom_selection")
    rotamer_score_and_choose_best(
      mon_lib_srv=mon_lib_srv,
      density_map=density_map,
      pdb_hierarchy=processed_pdb_file.all_chain_proxies.pdb_hierarchy,
      geometry_restraints_manager=grm,
      atom_selection_bool=atom_selection_bool,
      real_space_target_weight=work_params.real_space_target_weight,
      real_space_gradients_delta=real_space_gradients_delta,
      lbfgs_termination_params=scitbx.lbfgs.termination_parameters(
        max_iterations=work_params
          .rotamer_score_and_choose_best.lbfgs_max_iterations))
  #
  if (work_params.all_coordinate_refinement.run != 0):
    pdb_atoms = processed_pdb_file.all_chain_proxies.pdb_atoms
    sites_cart = pdb_atoms.extract_xyz()
    print "Before all coordinate refinement:"
    grm.energies_sites(sites_cart=sites_cart).show()
    print
    sys.stdout.flush()
    refined = maptbx.real_space_refinement_simple.lbfgs(
      sites_cart=sites_cart,
      density_map=density_map,
      geometry_restraints_manager=grm,
      real_space_target_weight=work_params.real_space_target_weight,
      real_space_gradients_delta=real_space_gradients_delta,
      lbfgs_termination_params=scitbx.lbfgs.termination_parameters(
        max_iterations=work_params.all_coordinate_refinement
          .lbfgs_max_iterations))
    print "After all coordinate refinement:"
    grm.energies_sites(sites_cart=refined.sites_cart).show()
    pdb_atoms.set_xyz(new_xyz=refined.sites_cart)
    print
    print "number_of_function_evaluations:", \
      refined.number_of_function_evaluations
    print "real+geo target start: %.6g" % refined.f_start
    print "real+geo target final: %.6g" % refined.f_final
    print
  #
  file_name = op.basename(input_pdb_file_name)
  if (   file_name.endswith(".pdb")
      or file_name.endswith(".ent")):
    file_name = file_name[:-4]
  file_name += "_lockit.pdb"
  pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  print "Writing file: %s" % show_string(file_name)
  sys.stdout.flush()
  pdb_hierarchy.write_pdb_file(
    file_name=file_name,
    crystal_symmetry=grm.crystal_symmetry)
  print
  #
  show_times()
  sys.stdout.flush()

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
