"""Developers tool to simulate Fobs data"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.fake_f_obs

from cctbx import adptbx
from cctbx.array_family import flex
import random, math, sys, os
import iotbx.pdb
import mmtbx.utils
from libtbx import easy_run
import mmtbx.dynamics.cartesian_dynamics as cartesian_dynamics
from mmtbx import monomer_library
import mmtbx.monomer_library.pdb_interpretation
import mmtbx.monomer_library.server
from mmtbx.tls import ladp
from mmtbx.utils import run_reduce_with_timeout
import mmtbx.tls.tools
import mmtbx.f_model
import iotbx.phil
import mmtbx.masks
from libtbx.utils import Sorry
from six.moves import range
import mmtbx.model

if(1):
  random.seed(0)
  flex.set_random_seed(0)

master_params_str="""\
f_obs {
  high_resolution = 2.0
    .type = float
  low_resolution = 15.0
    .type = float
  scattering_table = wk1995 it1992 *n_gaussian neutron
  f_calc {
    atomic_model {
      ensemble_size = 20
        .type = int
      add_hydrogens = False
        .type = bool
      tls {
        max_tl = 2
          .type = float
        min_tl = 0
          .type = float
      }
      apply_cartesian_dynamics = True
        .type = bool
      regularize_geometry {
        rmsd_bonds_target = 0.025
          .type = float
        rmsd_angles_target = 2.5
          .type = float
      }
      ladp_angle = 3.0
        .type = float
      switch_rotamers = True
        .type = bool
      shake_sites_rmsd = 0.01
        .type = float
      rigid_body_shift {
        rotation_angle = 1.0
          .type = float
        translation_length = 0.1
          .type = float
      }
      stop_cartesian_dynamics_at_diff = 0.5
        .type = float
      use_ramachandran_plot_restraints = True
        .type = bool
      output_file_name = fake_model.pdb
        .type = str
    }
    accuracy {
      include scope mmtbx.f_model.sf_and_grads_accuracy_master_params
    }
  }
  f_bulk {
    k_sol = 0.35
      .type = float
    b_sol = 50.0
      .type = float
    mask {
      include scope mmtbx.masks.mask_master_params
    }
  }
  overall_scale = 1.0
  overall_anisotropic_scale_matrix_b_cart {
    max = 10
      .type = float
    min = 0
      .type = float
  }
  experimental_noise {
    add_random_error_to_amplitudes_percent = 5
      .type = float
  }
  output_file_name = fake_f_obs.mtz
    .type = str
}

"""

class show(object):
  def __init__(self,
               xrs,
               xrs_start,
               grm,
               prefix=""):
    esg = grm.energies_sites(
      sites_cart = xrs.sites_cart(), compute_gradients = False).geometry
    self.bond_rmsd = esg.bond_deviations()[2]
    self.angle_rmsd = esg.angle_deviations()[2]
    self.error = flex.mean(xrs.distances(other = xrs_start))
    print("  %s err=%8.3f rmsd: bonds=%6.3f angles=%6.3f"%(prefix, self.error,
      self.bond_rmsd, self.angle_rmsd))

def switch_rotamers(xray_structure, pdb_hierarchy):
  x = xray_structure.deep_copy_scatterers()
  p = pdb_hierarchy.deep_copy()
  p.atoms().reset_i_seq()
  p = mmtbx.utils.switch_rotamers(
    pdb_hierarchy = p,
    mode = "min_distant")
  x.set_sites_cart(sites_cart = p.atoms().extract_xyz())
  return x, p

def set_ladp(xray_structure, pdb_hierarchy, angle):
  axes_and_atoms_i_seqs = ladp.get_axes_and_atoms_i_seqs(
    pdb_hierarchy = pdb_hierarchy,
    mon_lib_srv   = monomer_library.server.server())
  xray_structure = xray_structure.set_b_iso(value=random.randrange(5,10))
  xray_structure.convert_to_isotropic()
  xray_structure = ladp.set_ladp(
    xray_structure        = xray_structure,
    axes_and_atoms_i_seqs = axes_and_atoms_i_seqs,
    value                 = angle,
    enable_recursion      = True,
    depth                 = 0)
  return xray_structure

def random_aniso_adp(space_group, unit_cell, u_scale=2, u_min=0):
  return adptbx.u_star_as_u_cart(unit_cell, space_group.average_u_star(
    u_star = adptbx.u_cart_as_u_star(unit_cell, adptbx.random_u_cart(
      u_scale=u_scale, u_min=u_min))))

def apply_tls(xray_structure, params):
  uc = xray_structure.unit_cell()
  sg = xray_structure.space_group()
  selections_1d = flex.bool(xray_structure.scatterers().size(),True)
  selections = [selections_1d.iselection()]
  T=random_aniso_adp(space_group=sg, unit_cell=uc, u_scale=params.max_tl,
    u_min=params.min_tl)
  L=random_aniso_adp(space_group=sg, unit_cell=uc, u_scale=params.max_tl,
    u_min=params.min_tl)
  print("  T: %s"%",".join([("%7.3f"%i).strip() for i in T]))
  print("  L: %s"%",".join([("%7.3f"%i).strip() for i in L]))
  tlsos = mmtbx.tls.tools.generate_tlsos(
    selections     = selections,
    xray_structure = xray_structure,
    T=[T],
    L=[L],
    S=[[0,0,0,0,0,0,0,0,0]])
  u_cart_from_tls = mmtbx.tls.tools.u_cart_from_tls(
    sites_cart = xray_structure.sites_cart(),
    selections = selections,
    tlsos      = tlsos)
  xray_structure.convert_to_anisotropic()
  u_cart = xray_structure.scatterers().extract_u_cart(uc)
  utot = u_cart_from_tls+u_cart
  xray_structure.set_u_cart(u_cart=utot, selection = selections_1d.iselection())
  xray_structure.tidy_us()
  return xray_structure

def apply_rigid_body_shift(xray_structure, params):
  import scitbx.matrix
  mt = flex#.mersenne_twister(seed=0)
  rot_axis = scitbx.matrix.col(mt.random_double_point_on_sphere())
  rot_matrix = scitbx.math.r3_rotation_axis_and_angle_as_matrix(
    axis=rot_axis, angle=params.rotation_angle, deg=True)
  run_away_counter = 0
  while True:
    transl = mt.random_double_point_on_sphere()
    transl_no_cont_sh = scitbx.matrix.col(xray_structure.crystal_symmetry()
      .subtract_continuous_allowed_origin_shifts(translation_cart=transl))
    l = abs(transl_no_cont_sh)
    if(l > 0.1):
      break
    run_away_counter += 1
    print("run_away_counter", run_away_counter, l)

  assert run_away_counter < 100
  transl = transl_no_cont_sh * (params.translation_length/l)
  sites_cart = xray_structure.sites_cart()
  cm = xray_structure.center_of_mass()
  ns = rot_matrix * (sites_cart-cm) + transl + cm
  xray_structure.set_sites_cart(sites_cart =
    rot_matrix * (sites_cart-cm) + transl + cm)
  return xray_structure

def simulate_f_obs(root, crystal_symmetry, params):
  f_calc_data = None
  f_masks_data = []
  for i_m, m in enumerate(root.models()):
    raw_records = flex.std_string()
    raw_records.append(
      iotbx.pdb.format_cryst1_record(crystal_symmetry = crystal_symmetry))
    for atom in m.atoms():
      ra = atom.format_atom_record()
      ru = atom.format_anisou_record()
      raw_records.append(ra[:])
      raw_records.append(ru[:])
    xrs = iotbx.pdb.input(lines = raw_records,
      source_info=None).xray_structure_simple()
    if(i_m==0):
      dummy = abs(xrs.structure_factors(
        d_min=params.f_obs.high_resolution).f_calc())
      dummy = dummy.resolution_filter(d_max = params.f_obs.low_resolution)
    fmodel = mmtbx.f_model.manager(
      f_obs          = dummy,
      xray_structure = xrs,
      mask_params    = params.f_obs.f_bulk.mask,
      sf_and_grads_accuracy_params = params.f_obs.f_calc.accuracy)
    fcd = fmodel.f_calc().data()
    fms = fmodel.f_masks()
    if(i_m==0):
      f_calc_data = fcd
      f_masks_data = []
      for f in fms:
        f_masks_data.append(f.data())
    else:
      f_calc_data += fcd
      fmsks = fms
      assert len(f_masks_data) == len(fmsks)
      for ifmd in range(len(f_masks_data)):
        f_masks_data[ifmd] += fmsks[ifmd].data()
  fcalc_average = fmodel.f_obs().array(data = f_calc_data)
  f_masks_data_average = []
  for f in f_masks_data:
    f_masks_data_average.append(fmodel.f_obs().array(data = f/len(root.models())))
  b_cart = None
  if([params.f_obs.overall_anisotropic_scale_matrix_b_cart.max,
      params.f_obs.overall_anisotropic_scale_matrix_b_cart.min].count(None)==0):
    b_cart = random_aniso_adp(
      space_group=crystal_symmetry.space_group(),
      unit_cell=crystal_symmetry.unit_cell(),
      u_scale=params.f_obs.overall_anisotropic_scale_matrix_b_cart.max,
      u_min=params.f_obs.overall_anisotropic_scale_matrix_b_cart.min)
    print("\noverall_anisotropic_scale_matrix_b_cart: %s"%",".join(
      [("%7.3f"%i).strip() for i in b_cart]))
  fmodel = mmtbx.f_model.manager(
    f_obs  = dummy,
    f_calc = fcalc_average,
    f_mask = f_masks_data_average,
    k_sol  = params.f_obs.f_bulk.k_sol,
    b_sol  = params.f_obs.f_bulk.b_sol,
    b_cart = b_cart)
  #
  f_obs = abs(fmodel.f_model())
  f_obs.set_observation_type_xray_amplitude()
  mtz_dataset = f_obs.as_mtz_dataset(column_root_label="F(ake)obs")
  r_free_flags = f_obs.generate_r_free_flags()
  mtz_dataset.add_miller_array(
    miller_array=r_free_flags, column_root_label="R-free-flags")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name=params.f_obs.output_file_name)

def regularize_geometry(xray_structure, restraints_manager, params):
  from mmtbx.refinement import geometry_minimization as gm
  import scitbx.lbfgs
  sites_cart = xray_structure.sites_cart()
  minimized = gm.lbfgs(
    sites_cart = sites_cart,
    correct_special_position_tolerance = 1.0,
    geometry_restraints_manager = restraints_manager.geometry,
    geometry_restraints_flags = gm.geometry_restraints.flags.flags(default=True),
    rmsd_bonds_termination_cutoff=params.rmsd_bonds_target,
    rmsd_angles_termination_cutoff=params.rmsd_angles_target,
    lbfgs_termination_params=scitbx.lbfgs.termination_parameters(
      max_iterations=500))
  xray_structure = xray_structure.replace_sites_cart(new_sites = sites_cart)
  return xray_structure

def cd(xray_structure, restraints_manager, params):
  gradients_calculator=cartesian_dynamics.gradients_calculator_reciprocal_space(
    restraints_manager = restraints_manager,
    sites_cart         = xray_structure.sites_cart(),
    wc                 = 1)
  cartesian_dynamics.run(
    gradients_calculator             = gradients_calculator,
    xray_structure                   = xray_structure,
    temperature                      = 3000,
    n_steps                          = 500000,
    time_step                        = 0.0005,
    initial_velocities_zero_fraction = 0,
    n_print                          = 100,
    stop_cm_motion                   = True,
    log                              = None,
    stop_at_diff                     = params.stop_cartesian_dynamics_at_diff,
    verbose                          = -1)

def loop_2(params, xray_structure, pdb_hierarchy, restraints_manager, root):
  print("model:")
  amp = params.f_obs.f_calc.atomic_model
  grm = restraints_manager
  xrs = xray_structure.deep_copy_scatterers()
  show(xrs = xrs, xrs_start = xrs, grm = grm, prefix = "start:")
  xrs_sh = xrs.deep_copy_scatterers()
  if(amp.shake_sites_rmsd is not None):
    xrs_sh.shake_sites_in_place(rms_difference = amp.shake_sites_rmsd)
  if(amp.apply_cartesian_dynamics):
    cd(xray_structure = xrs_sh, restraints_manager = grm, params = amp)
    show(xrs = xrs_sh, xrs_start = xrs, grm = grm, prefix = "cd:   ")
  if([amp.regularize_geometry.rmsd_bonds_target,
      amp.regularize_geometry.rmsd_angles_target].count(None)==0):
    xrs_sh = regularize_geometry(xray_structure = xrs_sh,
      restraints_manager = grm, params = amp.regularize_geometry)
    show(xrs = xrs_sh, xrs_start = xrs, grm = grm, prefix = "min:  ")
  if(amp.ladp_angle is not None):
    xrs_sh = set_ladp(xray_structure = xrs_sh, pdb_hierarchy = pdb_hierarchy,
      angle = amp.ladp_angle)
  if([amp.tls.max_tl, amp.tls.min_tl].count(None)==0):
    xrs_sh = apply_tls(xray_structure = xrs_sh, params = amp.tls)
  if([amp.rigid_body_shift.rotation_angle,
      amp.rigid_body_shift.translation_length].count(None)==0):
    if xrs.crystal_symmetry().space_group().type().number()!=1:
      xrs_sh = apply_rigid_body_shift(xray_structure = xrs_sh,
        params = amp.rigid_body_shift)
      show(xrs = xrs_sh, xrs_start = xrs, grm = grm, prefix = "rb:   ")
  #
  h = pdb_hierarchy.deep_copy()
  h.atoms().reset_i_seq() # XXX
  h.atoms().set_xyz(xrs_sh.sites_cart().deep_copy())
  h.atoms().set_uij(xrs_sh.scatterers().extract_u_cart(xrs_sh.unit_cell()))
  h.atoms().set_b(xrs_sh.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1.))
  m = h.models()[0].detached_copy()
  m.id = str(None)
  root.append_model(m)

def loop_1(params, root, xray_structure, pdb_hierarchy, restraints_manager):
  xh = [(xray_structure,pdb_hierarchy)]
  if(params.f_obs.f_calc.atomic_model.switch_rotamers):
    xh.append(switch_rotamers(
      xray_structure = xray_structure.deep_copy_scatterers(),
      pdb_hierarchy = pdb_hierarchy.deep_copy()))
  counter = 0
  size = int(math.ceil(params.f_obs.f_calc.atomic_model.ensemble_size/len(xh)))
  for xh_ in xh:
    x_, h_ = xh_
    for mc in range(size):
      loop_2(
        params         = params,
        xray_structure = x_,
        pdb_hierarchy  = h_,
        restraints_manager = restraints_manager,
        root               = root)
  for i_model, model in enumerate(root.models()):
    model.id = str(i_model)
  root.atoms().set_occ(root.atoms().extract_occ()/len(root.models()))

def defaults(log):
  print("Default params::\n", file=log)
  parsed = iotbx.phil.parse(master_params_str, process_includes=True)
  print(file=log)
  return parsed

def run(args, log = sys.stdout):
  if(len(args)==0):
    parsed = defaults(log=log)
    parsed.show(prefix="  ", out=log)
    return
  parsed = defaults(log=log)
  processed_args = mmtbx.utils.process_command_line_args(args = args,
    log = sys.stdout, master_params = parsed)
  processed_args.params.show()
  params = processed_args.params.extract()
  if(len(processed_args.pdb_file_names)==0):
    raise Sorry("No PDB file found.")
  if(len(processed_args.pdb_file_names)>1):
    raise Sorry("More than one PDB file found.")
  pdb_file_name = processed_args.pdb_file_names[0]
  if(params.f_obs.f_calc.atomic_model.add_hydrogens):
    pdb_file_name_r = os.path.basename(pdb_file_name)+"_reduce"
    # easy_run.go("phenix.reduce %s > %s"% (pdb_file_name, pdb_file_name_r))
    run_reduce_with_timeout(file_name=pdb_file_name, parameters=" > %s" % pdb_file_name_r)
    pdb_file_name = pdb_file_name_r
  pdbi_params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  if(params.f_obs.f_calc.atomic_model.use_ramachandran_plot_restraints):
    pdbi_params.pdb_interpretation.ramachandran_plot_restraints.enabled=True
  model = mmtbx.model.manager(
    model_input = iotbx.pdb.input(file_name = pdb_file_name))
  model.process(make_restraints=True,
    pdb_interpretation_params = pdbi_params)
  root = iotbx.pdb.hierarchy.root()
  loop_1(
    params             = params,
    root               = root,
    xray_structure     = model.get_xray_structure(),
    pdb_hierarchy      = model.get_hierarchy(),
    restraints_manager = model.get_restraints_manager())
  root.write_pdb_file(
    file_name = params.f_obs.f_calc.atomic_model.output_file_name,
    crystal_symmetry = model.crystal_symmetry())
  simulate_f_obs(root=root, crystal_symmetry=model.crystal_symmetry(),
    params = params)

if (__name__ == "__main__"):
  run(sys.argv[1:])

