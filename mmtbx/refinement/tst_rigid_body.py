from cctbx.array_family import flex
from iotbx import pdb
import libtbx.load_env
from iotbx import pdb
import os, random
import mmtbx.refinement.rigid_body
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times
import random, math
from cctbx import xray
import mmtbx.utils
import iotbx.pdb
from scitbx import matrix
import mmtbx.command_line.fmodel

random.seed(0)
flex.set_random_seed(0)

def test_matrices_zyz():
  for i in xrange(10000):
    aa,bb,cc = random.randrange(-361,361),\
               random.randrange(-361,361),\
               random.randrange(-361,361)
    a,b,c = aa * math.pi/180, bb * math.pi/180, cc * math.pi/180
    r1 = matrix.sqr((math.cos(a), -math.sin(a), 0,
                     math.sin(a),  math.cos(a), 0,
                               0,            0, 1))
    r2 = matrix.sqr((math.cos(b), 0, math.sin(b),
                               0, 1,           0,
                    -math.sin(b), 0, math.cos(b)))
    r3 = matrix.sqr((math.cos(c), -math.sin(c), 0,
                     math.sin(c),  math.cos(c), 0,
                               0,            0, 1))
    r_zyz_1 = (r1*r2*r3)
    r_zyz_2 = mmtbx.refinement.rigid_body.euler(phi = cc, psi = bb, the = aa,
      convention = "zyz").rot_mat()
    r_zyz_3 = matrix.sqr((
      math.cos(a)*math.cos(b)*math.cos(c)-math.sin(a)*math.sin(c),-math.cos(a)*math.cos(b)*math.sin(c)-math.sin(a)*math.cos(c), math.cos(a)*math.sin(b),
      math.sin(a)*math.cos(b)*math.cos(c)+math.cos(a)*math.sin(c),-math.sin(a)*math.cos(b)*math.sin(c)+math.cos(a)*math.cos(c), math.sin(a)*math.sin(b),
     -math.sin(b)*math.cos(c)                                    , math.sin(b)*math.sin(c)                                    , math.cos(b)))
    assert approx_equal(r_zyz_1, r_zyz_2)
    assert approx_equal(r_zyz_1, r_zyz_3)

def test_1(fmodel, convention, phi = 0.0, psi = 0.0, the = 0.0,
           trans = [0.5,0.7,0.9]):
  rot_obj = mmtbx.refinement.rigid_body.euler(
    phi = phi, psi = psi, the = the, convention = convention)
  size = fmodel.xray_structure.scatterers().size()
  fmodel.xray_structure.apply_rigid_body_shift(
    rot = rot_obj.rot_mat(),
    trans = trans)
  fmodel.update_xray_structure(update_f_calc = True)
  assert fmodel.r_work() > 0.5
  params = mmtbx.refinement.rigid_body.master_params.extract()
  params.refine_rotation = False
  params.refine_translation = True
  params.target="ls_wunit_k1"
  rb = mmtbx.refinement.rigid_body.manager(
    fmodel     = fmodel,
    selections = [flex.bool(size,True).iselection()],
    params     = params)
  assert approx_equal(rb.translation()[0], [-0.5,-0.7,-0.9], eps=1.e-5)
  assert approx_equal(rb.rotation()[0], [0.0,0.0,0.0])
  assert approx_equal(fmodel.r_work(), 0.0, eps=1.e-5)

def test_2(fmodel, convention, phi = 0.0, psi = 0.0, the = 0.0,
           trans = [0.0,0.0,0.0]):
  rot_obj = mmtbx.refinement.rigid_body.euler(
    phi = phi, psi = psi, the = the, convention = convention)
  size = fmodel.xray_structure.scatterers().size()
  fmodel.xray_structure.apply_rigid_body_shift(
    rot = rot_obj.rot_mat(),
    trans = trans)
  fmodel.update_xray_structure(update_f_calc = True)
  assert approx_equal(fmodel.r_work(), 0.0)
  params = mmtbx.refinement.rigid_body.master_params.extract()
  params.refine_rotation = True
  params.refine_translation = True
  params.target="ls_wunit_k1"
  rb = mmtbx.refinement.rigid_body.manager(
    fmodel     = fmodel,
    selections = [flex.bool(size,True).iselection()],
    params     = params)
  assert approx_equal(rb.translation()[0], [0.0,0.0,0.0])
  assert approx_equal(rb.rotation()[0], [0.0,0.0,0.0])
  assert approx_equal(fmodel.r_work(), 0.0)

def test_3(fmodel, convention, phi = 1, psi = 2, the = 3, trans = [0,0,0]):
  rot_obj = mmtbx.refinement.rigid_body.euler(
    phi = phi, psi = psi, the = the, convention = convention)
  size = fmodel.xray_structure.scatterers().size()
  fmodel.xray_structure.apply_rigid_body_shift(
    rot = rot_obj.rot_mat(),
    trans = trans)
  fmodel.update_xray_structure(update_f_calc = True)
  assert fmodel.r_work() > 0.15
  params = mmtbx.refinement.rigid_body.master_params.extract()
  params.refine_rotation = True
  params.refine_translation = False
  params.high_resolution = 1.0
  params.max_iterations = 50
  params.lbfgs_line_search_max_function_evaluations = 50
  params.target="ls_wunit_k1"
  rb = mmtbx.refinement.rigid_body.manager(
    fmodel     = fmodel,
    selections = [flex.bool(size,True).iselection()],
    params     = params)
  if(convention == "xyz"):
    assert approx_equal(rb.rotation()[0], [-1.0,-2.0,-3.0], 0.2) # XXX
  assert approx_equal(rb.translation()[0], [0.0,0.0,0.0])
  assert approx_equal(fmodel.r_work(), 0.0)

def test_4(fmodel, convention, phi =1, psi =2, the =3, trans =[0.5,1.0,1.5]):
  rot_obj = mmtbx.refinement.rigid_body.euler(
    phi = phi, psi = psi, the = the, convention = convention)
  size = fmodel.xray_structure.scatterers().size()
  fmodel.xray_structure.apply_rigid_body_shift(
    rot = rot_obj.rot_mat(),
    trans = trans)
  fmodel.update_xray_structure(update_f_calc = True)
  assert fmodel.r_work() > 0.15
  params = mmtbx.refinement.rigid_body.master_params.extract()
  params.refine_rotation = True
  params.refine_translation = True
  params.high_resolution = 1.0
  params.target="ls_wunit_k1"
  rb = mmtbx.refinement.rigid_body.manager(
    fmodel     = fmodel,
    selections = [flex.bool(size,True).iselection()],
    params     = params)
  if(convention == "xyz"):
    assert approx_equal(rb.rotation()[0], [-1.0,-2.0,-3.0], 0.2) # XXX
  assert approx_equal(rb.translation()[0], [-0.5,-1.0,-1.5], 1.e-3)
  assert approx_equal(fmodel.r_work(), 0.0, 1.e-3)

def test_5(fmodel, convention):
  size = fmodel.xray_structure.scatterers().size()
  sel1 = flex.bool()
  sel2 = flex.bool()
  for i in xrange(size):
    if(i<500):
       sel1.append(True)
       sel2.append(False)
    else:
       sel1.append(False)
       sel2.append(True)
  selections = [sel1.iselection(), sel2.iselection()]
  rot_obj_1 = mmtbx.refinement.rigid_body.euler(
    phi = 1, psi = 2, the = 3, convention = convention)
  rot_obj_2 = mmtbx.refinement.rigid_body.euler(
    phi = 3, psi = 2, the = 1, convention = convention)
  fmodel.xray_structure.apply_rigid_body_shift(
    rot = rot_obj_1.rot_mat(),
    trans = [0.5,1.0,1.5],
    selection = sel1.iselection())
  fmodel.xray_structure.apply_rigid_body_shift(
    rot = rot_obj_2.rot_mat(),
    trans = [1.5,0.5,1.0],
    selection = sel2.iselection())
  fmodel.update_xray_structure(update_f_calc = True)
  assert fmodel.r_work() > 0.35
  params = mmtbx.refinement.rigid_body.master_params.extract()
  params.refine_rotation = True
  params.refine_translation = True
  params.high_resolution = 1.0
  params.target="ls_wunit_k1"
  rb = mmtbx.refinement.rigid_body.manager(fmodel     = fmodel,
                                           selections = selections,
                                           params     = params)
  if(convention == "xyz"):
    assert approx_equal(rb.rotation()[0], [-1,-2,-3], 0.2)
    assert approx_equal(rb.rotation()[1], [-3,-2,-1], 0.2)
  assert approx_equal(rb.translation()[0], [-0.5,-1.0,-1.5], 1.e-4)
  assert approx_equal(rb.translation()[1], [-1.5,-0.5,-1.0], 1.e-4)
  assert approx_equal(fmodel.r_work(), 0.0, 0.0005)
  assert approx_equal(fmodel.r_free(), 0.0, 0.0005)

def get_fmodel_from_pdb(pdb_file_name,
                        algorithm = "direct",
                        d_min = 2.0,
                        target = "ls_wunit_k1"):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/%s"%pdb_file_name,test=os.path.isfile)
  xray_structure = iotbx.pdb.input(file_name =pdb_file).xray_structure_simple()
  params = mmtbx.command_line.fmodel.fmodel_from_xray_structure_master_params.extract()
  assert params.high_resolution is None
  assert params.scattering_table == 'n_gaussian'
  assert params.structure_factors_accuracy.algorithm == 'fft'
  params.high_resolution = d_min
  params.scattering_table = "wk1995"
  params.structure_factors_accuracy.algorithm = algorithm
  fmodel = mmtbx.utils.fmodel_from_xray_structure(
    xray_structure        = xray_structure,
    params                = params,
    target                = target,
    r_free_flags_fraction = 0.01).fmodel
  assert approx_equal(fmodel.r_work(), 0)
  return fmodel

def run_tests():
  fmodel_small = get_fmodel_from_pdb(pdb_file_name = "enk_gbr.pdb")
  fmodel_big = get_fmodel_from_pdb(pdb_file_name = "lys_rigid_body.pdb",
    algorithm = "fft")
  for convention in ["xyz", "zyz"]:
    print "test 1: ", convention
    test_1(fmodel = fmodel_small.deep_copy(), convention = convention)
    print "test 2: ", convention
    test_2(fmodel = fmodel_small.deep_copy(), convention = convention)
    print "test 3: ", convention
    test_3(fmodel = fmodel_small.deep_copy(), convention = convention)
    print "test 4: ", convention
    test_4(fmodel = fmodel_big.deep_copy(), convention = convention)
    print "test 5: ", convention
    test_5(fmodel = fmodel_big.deep_copy(), convention = convention)


def finite_differences_test():
  print "finite_differences_test: "
  fmodel = get_fmodel_from_pdb(pdb_file_name = "enk_rbr.pdb",
                               algorithm = "direct",
                               d_min = 2.0,
                               target = "ls_wunit_k1")
  xray.set_scatterer_grad_flags(scatterers = fmodel.xray_structure.scatterers(),
                                site       = True)
  for convention in ["zyz","xyz"]:
      rot_obj = mmtbx.refinement.rigid_body.euler(
        phi = 0, psi = 0, the = 0, convention = convention)
      size = fmodel.xray_structure.scatterers().size()
      selections = [flex.bool(size, True).iselection()]
      fmodel.xray_structure.apply_rigid_body_shift(
        rot = rot_obj.rot_mat(), trans = [1,2,3])
      fmodel.update_xray_structure(update_f_calc = True)
      fmodel_copy = fmodel.deep_copy()
      centers_of_mass = []
      for s in selections:
        xrs = fmodel_copy.xray_structure.select(s)
        centers_of_mass.append(xrs.center_of_mass())
      tg_obj = mmtbx.refinement.rigid_body.target_and_grads(
        centers_of_mass = centers_of_mass,
        sites_cart      = fmodel_copy.xray_structure.sites_cart(),
        target_functor  = fmodel_copy.target_functor(),
        rot_objs        = [rot_obj],
        selections      = selections,
        suppress_gradients = False)
      assert approx_equal(tg_obj.target(),fmodel_copy.target_w())
      g_rot, g_transl = tg_obj.gradients_wrt_r(), tg_obj.gradients_wrt_t()
      fd_transl = fd_translation(fmodel_copy, e = 0.0001)
      assert approx_equal(list(g_transl[0]), fd_transl)
      fd_rot = fd_rotation(fmodel     = fmodel_copy,
                           e          = 0.0001,
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
      t1 = fm.target_functor()().target_work()
      fm.update_xray_structure(xray_structure = xrs_g1_2,
                               update_f_calc  = True)
      t2 = fm.target_functor()().target_work()
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
      t1 = fm.target_functor()().target_work()
      fm.update_xray_structure(xray_structure = xrs_g1_2,
                               update_f_calc  = True)
      t2 = fm.target_functor()().target_work()
      grads.append( (t1-t2)/(2*e*math.pi/180)  )
  return grads


def exercise():
  test_matrices_zyz()
  run_tests()
  finite_differences_test()
  print format_cpu_times()

if (__name__ == "__main__"):
  exercise()
