from cctbx.array_family import flex
from iotbx import pdb
import iotbx.pdb.interpretation
import mmtbx.f_model
from libtbx import introspection
import libtbx.load_env
from iotbx import pdb
import sys, os, time
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
import mmtbx.refinement.rigid_body
import mmtbx.model
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times
import random, math
from cctbx import xray

random.seed(0)
flex.set_random_seed(0)

def test_matrices():
  rot_obj = mmtbx.refinement.rigid_body.rb_mat(phi = 10,
                                               psi = 20,
                                               the = 30)

  rot = [0.92541657839832336, -0.16317591116653482,  0.34202014332566871,
         0.31879577759716782,  0.82317294464550095, -0.46984631039295416,
        -0.20487412870286215,  0.54383814248232554,  0.8137976813493738]
  assert approx_equal(rot, rot_obj.rot_mat())

  r_phi = [-0.16317591116653482, -0.92541657839832336, 0,
            0.82317294464550095, -0.31879577759716782, 0,
            0.54383814248232554,  0.20487412870286215, 0]
  assert approx_equal(r_phi, rot_obj.r_phi())

  r_psi = [-0.33682408883346515,  0.059391174613884698,  0.93969262078590843,
            0.46270828919916163, -0.081587955583267396,  0.17101007166283433,
           -0.80143426597622169,  0.14131448435589203 , -0.29619813272602386]
  assert approx_equal(r_psi, rot_obj.r_psi())

  r_the = [0                  ,                    0,                    0,
           0.20487412870286215, -0.54383814248232554, -0.8137976813493738 ,
           0.31879577759716782,  0.82317294464550095, -0.46984631039295416]
  assert approx_equal(r_the, rot_obj.r_the())

def test_1(fmodel, model):
  rot_obj = mmtbx.refinement.rigid_body.rb_mat(phi = 0.0,
                                               psi = 0.0,
                                               the = 0.0)
  dim = fmodel.xray_structure.scatterers().size()
  selections = [flex.bool(dim, True)]
  ss_ = []
  for s in selections:
    ss_.append(s.iselection())
  selections = ss_
  new_sites_frac, new_sites_cart, centers_of_mass = \
                             mmtbx.refinement.rigid_body.apply_transformation_(
                  xray_structure      = fmodel.xray_structure,
                  sites_cart          = fmodel.xray_structure.sites_cart(),
                  sites_frac          = fmodel.xray_structure.sites_frac(),
                  rotation_matrices   = [rot_obj.rot_mat()],
                  translation_vectors = [(2.5,2.5,2.5)],
                  selections          = selections,
                  atomic_weights      = fmodel.xray_structure.atomic_weights())
  fmodel.xray_structure.set_sites_frac(new_sites_frac)
  new_xray_structure = fmodel.xray_structure

  fmodel.update_xray_structure(xray_structure = new_xray_structure,
                               update_f_calc  = True)
  assert approx_equal(fmodel.r_work(), 0.52, 1.e-2)
  assert approx_equal(fmodel.r_free(), 0.52, 1.e-2)
  rb = mmtbx.refinement.rigid_body.manager(fmodel           = fmodel,
                                           selections       = selections,
                                           refine_r         = False,
                                           refine_t         = True,
                                           convergence_test = True,
                                           protocol         = "multiple_zones")
  assert approx_equal(rb.translation()[0], [-2.5,-2.5,-2.5], 1.e-4)
  assert approx_equal(rb.rotation()[0], [0.0,0.0,0.0])
  assert approx_equal(fmodel.r_work(), 0.0, 1.e-3)
  assert approx_equal(fmodel.r_free(), 0.0, 1.e-3)

def test_2(fmodel, model):
  rot_obj = mmtbx.refinement.rigid_body.rb_mat(phi = 0.0,
                                               psi = 0.0,
                                               the = 0.0)
  dim = fmodel.xray_structure.scatterers().size()
  selections = [flex.bool(dim, True)]
  ss_ = []
  for s in selections:
    ss_.append(s.iselection())
  selections = ss_
  new_sites_frac, new_sites_cart, centers_of_mass = \
                             mmtbx.refinement.rigid_body.apply_transformation_(
                  xray_structure      = fmodel.xray_structure,
                  sites_cart          = fmodel.xray_structure.sites_cart(),
                  sites_frac          = fmodel.xray_structure.sites_frac(),
                  rotation_matrices   = [rot_obj.rot_mat()],
                  translation_vectors = [(1.5,1.5,1.5)],
                  selections          = selections,
                  atomic_weights      = fmodel.xray_structure.atomic_weights())
  fmodel.xray_structure.set_sites_frac(new_sites_frac)
  new_xray_structure = fmodel.xray_structure

  fmodel.update_xray_structure(xray_structure = new_xray_structure,
                               update_f_calc  = True)
  assert approx_equal(fmodel.r_work(), 0.53, 1.e-2)
  assert approx_equal(fmodel.r_free(), 0.51, 1.e-2)
  rb = mmtbx.refinement.rigid_body.manager(fmodel           = fmodel,
                                           selections       = selections,
                                           refine_r         = True,
                                           refine_t         = True,
                                           convergence_test = False,
                                           protocol         = "multiple_zones")
  assert approx_equal(rb.translation()[0], [-1.5,-1.5,-1.5], 1.e-4)
  assert approx_equal(rb.rotation()[0], [0.0,0.0,0.0], 1.e-2)
  assert approx_equal(fmodel.r_work(), 0.0, 1.e-3)
  assert approx_equal(fmodel.r_free(), 0.0, 1.e-3)

def test_3(fmodel, model):
  rot_obj = mmtbx.refinement.rigid_body.rb_mat(phi = 0.0,
                                               psi = 0.0,
                                               the = 0.0)
  dim = fmodel.xray_structure.scatterers().size()
  selections = [flex.bool(dim, True)]
  ss_ = []
  for s in selections:
    ss_.append(s.iselection())
  selections = ss_
  new_sites_frac, new_sites_cart, centers_of_mass = \
                             mmtbx.refinement.rigid_body.apply_transformation_(
                  xray_structure      = fmodel.xray_structure,
                  sites_cart          = fmodel.xray_structure.sites_cart(),
                  sites_frac          = fmodel.xray_structure.sites_frac(),
                  rotation_matrices   = [rot_obj.rot_mat()],
                  translation_vectors = [(0,0,0)],
                  selections          = selections,
                  atomic_weights      = fmodel.xray_structure.atomic_weights())
  fmodel.xray_structure.set_sites_frac(new_sites_frac)
  new_xray_structure = fmodel.xray_structure
  fmodel.update_xray_structure(xray_structure = new_xray_structure,
                               update_f_calc  = True)
  assert approx_equal(fmodel.r_work(), 0.0, 1.e-3)
  assert approx_equal(fmodel.r_free(), 0.0, 1.e-3)
  rb = mmtbx.refinement.rigid_body.manager(fmodel           = fmodel,
                                           selections       = selections,
                                           refine_r         = True,
                                           refine_t         = True,
                                           convergence_test = True,
                                           protocol         = "multiple_zones")
  assert approx_equal(rb.translation()[0], [0.0,0.0,0.0], 1.e-4)
  assert approx_equal(rb.rotation()[0], [0.0,0.0,0.0], 1.e-3)
  assert approx_equal(fmodel.r_work(), 0.0, 1.e-3)
  assert approx_equal(fmodel.r_free(), 0.0, 1.e-3)

def test_4(fmodel, model):
  rot_obj = mmtbx.refinement.rigid_body.rb_mat(phi = 0.5,
                                               psi = 0.5,
                                               the = 0.5)
  dim = fmodel.xray_structure.scatterers().size()
  selections = [flex.bool(dim, True)]
  ss_ = []
  for s in selections:
    ss_.append(s.iselection())
  selections = ss_
  new_sites_frac, new_sites_cart, centers_of_mass = \
                             mmtbx.refinement.rigid_body.apply_transformation_(
                  xray_structure      = fmodel.xray_structure,
                  sites_cart          = fmodel.xray_structure.sites_cart(),
                  sites_frac          = fmodel.xray_structure.sites_frac(),
                  rotation_matrices   = [rot_obj.rot_mat()],
                  translation_vectors = [(0,0,0)],
                  selections          = selections,
                  atomic_weights      = fmodel.xray_structure.atomic_weights())
  fmodel.xray_structure.set_sites_frac(new_sites_frac)
  new_xray_structure = fmodel.xray_structure
  fmodel.update_xray_structure(xray_structure = new_xray_structure,
                               update_f_calc  = True)
  assert approx_equal(fmodel.r_work(), 0.12, 1.e-2)
  assert approx_equal(fmodel.r_free(), 0.13, 1.e-2)
  rb = mmtbx.refinement.rigid_body.manager(fmodel           = fmodel,
                                           selections       = selections,
                                           refine_r         = True,
                                           refine_t         = False,
                                           convergence_test = True,
                                           protocol         = "multiple_zones")
  assert approx_equal(rb.translation()[0], [0.0,0.0,0.0], 1.e-4)
  assert approx_equal(rb.rotation()[0], [-0.5,-0.5,-0.5], 1.e-1)
  assert approx_equal(fmodel.r_work(), 0.0, 1.e-3)
  assert approx_equal(fmodel.r_free(), 0.0, 1.e-3)

def test_5(fmodel, model):
  rot_obj = mmtbx.refinement.rigid_body.rb_mat(phi = 0.5,
                                               psi = 0.5,
                                               the = 0.5)
  dim = fmodel.xray_structure.scatterers().size()
  selections = [flex.bool(dim, True)]
  ss_ = []
  for s in selections:
    ss_.append(s.iselection())
  selections = ss_
  new_sites_frac, new_sites_cart, centers_of_mass = \
                             mmtbx.refinement.rigid_body.apply_transformation_(
                  xray_structure      = fmodel.xray_structure,
                  sites_cart          = fmodel.xray_structure.sites_cart(),
                  sites_frac          = fmodel.xray_structure.sites_frac(),
                  rotation_matrices   = [rot_obj.rot_mat()],
                  translation_vectors = [(1,1,1)],
                  selections          = selections,
                  atomic_weights      = fmodel.xray_structure.atomic_weights())
  fmodel.xray_structure.set_sites_frac(new_sites_frac)
  new_xray_structure = fmodel.xray_structure
  fmodel.update_xray_structure(xray_structure = new_xray_structure,
                               update_f_calc  = True)
  assert approx_equal(fmodel.r_work(), 0.51, 1.e-2)
  assert approx_equal(fmodel.r_free(), 0.52, 1.e-2)
  fmodel.show_comprehensive(reflections_per_bin = 250,
                            max_number_of_bins  = 30)

  rb = mmtbx.refinement.rigid_body.manager(fmodel           = fmodel,
                                           selections       = selections,
                                           refine_r         = True,
                                           refine_t         = True,
                                           convergence_test = True,
                                           protocol         = "multiple_zones")
  fmodel.show_comprehensive(reflections_per_bin = 250,
                            max_number_of_bins  = 30)
  assert approx_equal(rb.translation()[0], [-1.0,-1.0,-1.0], 1.e-4)
  assert approx_equal(rb.rotation()[0], [-0.5,-0.5,-0.5], 1.e-1)
  assert approx_equal(fmodel.r_work(), 0.0, 1.e-3)
  assert approx_equal(fmodel.r_free(), 0.0, 1.e-3)

def test_6(fmodel, model):
  elections = []
  dim = fmodel.xray_structure.scatterers().size()
  sel1 = flex.bool()
  sel2 = flex.bool()
  for i in xrange(dim):
    if(i<500):
       sel1.append(True)
       sel2.append(False)
    else:
       sel1.append(False)
       sel2.append(True)
  selections = [sel1, sel2]
  ss_ = []
  for s in selections:
    ss_.append(s.iselection())
  selections = ss_

  rot_obj_1 = mmtbx.refinement.rigid_body.rb_mat(phi = 0.5,
                                                 psi = 0.5,
                                                 the = 0.5)
  rot_obj_2 = mmtbx.refinement.rigid_body.rb_mat(phi = 0.3,
                                                 psi = 0.3,
                                                 the = 0.3)
  new_sites_frac, new_sites_cart, centers_of_mass = \
                             mmtbx.refinement.rigid_body.apply_transformation_(
                  xray_structure      = model.xray_structure,
                  sites_cart          = model.xray_structure.sites_cart(),
                  sites_frac          = model.xray_structure.sites_frac(),
                  rotation_matrices   = [rot_obj_1.rot_mat(),rot_obj_2.rot_mat()],
                  translation_vectors = [(1,1,1),(1.2,1.2,1.2)],
                  selections          = selections,
                  atomic_weights      = fmodel.xray_structure.atomic_weights())
  fmodel.xray_structure.set_sites_frac(new_sites_frac)
  new_xray_structure = fmodel.xray_structure
  fmodel.update_xray_structure(xray_structure = new_xray_structure,
                               update_f_calc  = True)
  assert approx_equal(fmodel.r_work(), 0.52, 1.e-1)
  assert approx_equal(fmodel.r_free(), 0.53, 1.e-1)
  rb = mmtbx.refinement.rigid_body.manager(fmodel           = fmodel,
                                           selections       = selections,
                                           refine_r         = True,
                                           refine_t         = True,
                                           convergence_test = True,
                                           protocol         = "multiple_zones")
  assert approx_equal(rb.translation()[0], [-1.0,-1.0,-1.0], 1.e-4)
  assert approx_equal(rb.translation()[1], [-1.2,-1.2,-1.2], 1.e-4)
  assert approx_equal(rb.rotation()[0], [-0.5,-0.5,-0.5], 1.e-1)
  assert approx_equal(rb.rotation()[1], [-0.3,-0.3,-0.3], 1.e-1)
  assert approx_equal(fmodel.r_work(), 0.0, 1.e-3)
  assert approx_equal(fmodel.r_free(), 0.0, 1.e-3)

def run_tests(sf_algorithm = "fft"):
  pdb_file = libtbx.env.find_in_repositories(
        relative_path="regression/pdb/lys_rigid_body.pdb", test=os.path.isfile)
  processed_pdb_file = monomer_library.pdb_interpretation.process(
                            mon_lib_srv    = monomer_library.server.server(),
                            ener_lib       = monomer_library.server.ener_lib(),
                            file_name      = pdb_file,
                            raw_records    = None,
                            force_symmetry = True)
################
  geometry = processed_pdb_file.geometry_restraints_manager(
                                                    show_energies      = False,
                                                    plain_pairs_radius = 5.0)
  restraints_manager = mmtbx.restraints.manager(geometry      = geometry,
                                                normalization = False)
  xray_structure = processed_pdb_file.xray_structure()
  selection = flex.bool(xray_structure.scatterers().size(), True)
  restraints_manager_ini = mmtbx.restraints.manager(
                                  geometry      = geometry.select(selection),
                                  normalization = False)
  aal= processed_pdb_file.all_chain_proxies.stage_1.atom_attributes_list
  model = mmtbx.model.manager(
             restraints_manager     = restraints_manager,
             restraints_manager_ini = restraints_manager_ini,
             xray_structure         = xray_structure,
             atom_attributes_list   = aal)
  model.xray_structure.scattering_type_registry(table = "wk1995")
################
  dummy = model.xray_structure.structure_factors(algorithm = sf_algorithm,
                                                 d_min     = 2.0).f_calc()
  f_obs = abs(dummy.structure_factors_from_scatterers(
                                         xray_structure = model.xray_structure,
                                         algorithm      = sf_algorithm,
                                         cos_sin_table  = False).f_calc())
  flags =f_obs.array(data=flex.size_t(xrange(1,f_obs.data().size()+1))%10 == 0)
  fmodel = mmtbx.f_model.manager(xray_structure    = model.xray_structure,
                                 f_obs             = f_obs,
                                 r_free_flags      = flags,
                                 target_name       = "ls_wunit_k1",
                                 sf_cos_sin_table  = False,
                                 #target_name       = "ml",
                                 sf_algorithm      = sf_algorithm)
  fmodel.show_comprehensive(reflections_per_bin = 250,
                            max_number_of_bins  = 30)

  print "test 1: "
  test_1(fmodel = fmodel.deep_copy(), model  = model.deep_copy())
  if (not "--comprehensive" in sys.argv[1:]): return
  print "test 2: "
  test_2(fmodel = fmodel.deep_copy(), model  = model.deep_copy())
  print "test 3: "
  test_3(fmodel = fmodel.deep_copy(), model  = model.deep_copy())
  print "test 4: "
  test_4(fmodel = fmodel.deep_copy(), model  = model.deep_copy())
  print "test 5: "
  test_5(fmodel = fmodel.deep_copy(), model  = model.deep_copy())
  print "test 6: "
  test_6(fmodel = fmodel.deep_copy(), model  = model.deep_copy())


def finite_differences_test(sf_algorithm = "direct"):
  print "finite_differences_test: "
  pdb_file = libtbx.env.find_in_repositories(
        relative_path="regression/pdb/enk_rbr.pdb", test=os.path.isfile)
  processed_pdb_file = monomer_library.pdb_interpretation.process(
                            mon_lib_srv    = monomer_library.server.server(),
                            ener_lib       = monomer_library.server.ener_lib(),
                            file_name      = pdb_file,
                            raw_records    = None,
                            force_symmetry = True)
################
  geometry = processed_pdb_file.geometry_restraints_manager(
                                                    show_energies      = False,
                                                    plain_pairs_radius = 5.0)
  restraints_manager = mmtbx.restraints.manager(geometry      = geometry,
                                                normalization = False)
  xray_structure = processed_pdb_file.xray_structure()
  selection = flex.bool(xray_structure.scatterers().size(), True)
  restraints_manager_ini = mmtbx.restraints.manager(
                                  geometry      = geometry.select(selection),
                                  normalization = False)
  aal= processed_pdb_file.all_chain_proxies.stage_1.atom_attributes_list
  model = mmtbx.model.manager(
             restraints_manager     = restraints_manager,
             restraints_manager_ini = restraints_manager_ini,
             xray_structure         = xray_structure,
             atom_attributes_list   = aal)
  model.xray_structure.scattering_type_registry(table = "wk1995")
################
  dummy = model.xray_structure.structure_factors(algorithm = sf_algorithm,
                                                 d_min     = 1.0).f_calc()
  f_obs = abs(dummy.structure_factors_from_scatterers(
                                         xray_structure = model.xray_structure,
                                         algorithm      = sf_algorithm,
                                         cos_sin_table  = False).f_calc())
  flags = f_obs.array(data = flex.bool(f_obs.data().size(), False))
  model.xray_structure.shake_sites(mean_error = 0.1)
  fmodel = mmtbx.f_model.manager(xray_structure    = model.xray_structure,
                                 f_obs             = f_obs,
                                 r_free_flags      = flags,
                                 target_name       = "ls_wunit_k1",
                                 sf_cos_sin_table  = False,
                                 sf_algorithm      = sf_algorithm)
  fmodel.show_essential()
  xray.set_scatterer_grad_flags(scatterers = fmodel.xray_structure.scatterers(),
                                site       = True)
  for convention in ["zyz","xyz"]:
      rot_obj = mmtbx.refinement.rigid_body.euler(
                            phi = 0, psi = 0, the = 0, convention = convention)
      dim = fmodel.xray_structure.scatterers().size()
      selections = [flex.bool(dim, True)]
      selections = [flex.bool(dim, True)]
      ss_ = []
      for s in selections:
        ss_.append(s.iselection())
      selections = ss_
      new_xray_structure = mmtbx.refinement.rigid_body.apply_transformation(
                       xray_structure      = model.xray_structure,
                       rotation_matrices   = [rot_obj.rot_mat()],
                       translation_vectors = [(1.0,2.0,3.0)],
                       selections          = selections)
      fmodel_copy = fmodel.deep_copy()
      fmodel_copy.update_xray_structure(xray_structure = new_xray_structure,
                                        update_f_calc  = True)
      centers_of_mass = []
      for s in selections:
        xrs = fmodel_copy.xray_structure.select(s)
        centers_of_mass.append(xrs.center_of_mass())
      tg_obj = mmtbx.refinement.rigid_body.target_and_grads(
                                             centers_of_mass = centers_of_mass,
                                             sites_cart      = fmodel_copy.xray_structure.sites_cart(),
                                             fmodel          = fmodel_copy,
                                             alpha           = None,
                                             beta            = None,
                                             rot_objs        = [rot_obj],
                                             selections      = selections,
                                             suppress_gradients = False)
      assert approx_equal(tg_obj.target(),fmodel_copy.target_w())
      g_rot, g_transl = tg_obj.gradients_wrt_r(), tg_obj.gradients_wrt_t()
      fd_transl = fd_translation(fmodel_copy, e = 0.00001)
      assert approx_equal(list(g_transl[0]), fd_transl)
      fd_rot = fd_rotation(fmodel     = fmodel_copy,
                           e          = 0.00001,
                           convention = convention)
      assert approx_equal(list(g_rot[0]), fd_rot)


def fd_translation(fmodel, e):
  grads = []
  for shift in [[e,0,0],[0,e,0],[0,0,e]]:
      xrs1 = fmodel.xray_structure.deep_copy_scatterers()
      xrs2 = fmodel.xray_structure.deep_copy_scatterers()
      fm = fmodel.deep_copy()
      xrs_g1_1 = xrs1.translate(x = shift[0], y = shift[1], z = shift[2])
      xrs_g1_2 = xrs2.translate(x =-shift[0], y =-shift[1], z =-shift[2])
      fm.update_xray_structure(xray_structure = xrs_g1_1,
                               update_f_calc  = True)
      t1 = fm.target_w()
      fm.update_xray_structure(xray_structure = xrs_g1_2,
                               update_f_calc  = True)
      t2 = fm.target_w()
      grads.append( (t1-t2)/(2*e) )
  return grads

def fd_rotation(fmodel, e, convention):
  grads = []
  for shift in [[e,0,0],[0,e,0],[0,0,e]]:
      xrs1 = fmodel.xray_structure.deep_copy_scatterers()
      xrs2 = fmodel.xray_structure.deep_copy_scatterers()
      fm = fmodel.deep_copy()
      #
      rot_obj = mmtbx.refinement.rigid_body.euler(phi = shift[0],
                                                  psi = shift[1],
                                                  the = shift[2],
                                                  convention = convention)
      selections = [flex.bool(xrs1.scatterers().size(), True)]
      xrs_g1_1 = mmtbx.refinement.rigid_body.apply_transformation(
                         xray_structure      = xrs1,
                         rotation_matrices   = [rot_obj.rot_mat()],
                         translation_vectors = [(0.0,0.0,0.0)],
                         selections          = selections)
      #
      rot_obj = mmtbx.refinement.rigid_body.euler(phi = -shift[0],
                                                  psi = -shift[1],
                                                  the = -shift[2],
                                                  convention = convention)
      selections = [flex.bool(xrs1.scatterers().size(), True)]
      xrs_g1_2 = mmtbx.refinement.rigid_body.apply_transformation(
                         xray_structure      = xrs2,
                         rotation_matrices   = [rot_obj.rot_mat()],
                         translation_vectors = [(0.0,0.0,0.0)],
                         selections          = selections)
      fm.update_xray_structure(xray_structure = xrs_g1_1,
                               update_f_calc  = True)
      t1 = fm.target_w()
      fm.update_xray_structure(xray_structure = xrs_g1_2,
                               update_f_calc  = True)
      t2 = fm.target_w()
      grads.append( (t1-t2)/(2*e*math.pi/180)  )
  return grads


def exercise():
  test_matrices()
  run_tests()
  finite_differences_test()
  print format_cpu_times()

if (__name__ == "__main__"):
  exercise()
