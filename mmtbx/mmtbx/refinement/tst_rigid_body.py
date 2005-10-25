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
  new_xray_structure = mmtbx.refinement.rigid_body.apply_transformation(
                         xray_structure      = model.xray_structure,
                         rotation_matrices   = [rot_obj.rot_mat()],
                         translation_vectors = [(2.5,2.5,2.5)],
                         selections          = selections)
  fmodel.update_xray_structure(xray_structure = new_xray_structure,
                               update_f_calc  = True)
  assert approx_equal(fmodel.r_work(), 0.52, 1.e-2)
  assert approx_equal(fmodel.r_free(), 0.52, 1.e-2)
  rb = mmtbx.refinement.rigid_body.manager(fmodel           = fmodel,
                                           selections       = selections,
                                           refine_r         = False,
                                           refine_t         = True,
                                           convergence_test = True)
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
  new_xray_structure = mmtbx.refinement.rigid_body.apply_transformation(
                         xray_structure      = model.xray_structure,
                         rotation_matrices   = [rot_obj.rot_mat()],
                         translation_vectors = [(1.5,1.5,1.5)],
                         selections          = selections)
  fmodel.update_xray_structure(xray_structure = new_xray_structure,
                               update_f_calc  = True)
  assert approx_equal(fmodel.r_work(), 0.53, 1.e-2)
  assert approx_equal(fmodel.r_free(), 0.51, 1.e-2)
  rb = mmtbx.refinement.rigid_body.manager(fmodel           = fmodel,
                                           selections       = selections,
                                           refine_r         = True,
                                           refine_t         = True,
                                           convergence_test = False)
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
  new_xray_structure = mmtbx.refinement.rigid_body.apply_transformation(
                         xray_structure      = model.xray_structure,
                         rotation_matrices   = [rot_obj.rot_mat()],
                         translation_vectors = [(0,0,0)],
                         selections          = selections)
  fmodel.update_xray_structure(xray_structure = new_xray_structure,
                               update_f_calc  = True)
  assert approx_equal(fmodel.r_work(), 0.0, 1.e-3)
  assert approx_equal(fmodel.r_free(), 0.0, 1.e-3)
  rb = mmtbx.refinement.rigid_body.manager(fmodel           = fmodel,
                                           selections       = selections,
                                           refine_r         = True,
                                           refine_t         = True,
                                           convergence_test = True)
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
  new_xray_structure = mmtbx.refinement.rigid_body.apply_transformation(
                         xray_structure      = model.xray_structure,
                         rotation_matrices   = [rot_obj.rot_mat()],
                         translation_vectors = [(0,0,0)],
                         selections          = selections)
  fmodel.update_xray_structure(xray_structure = new_xray_structure,
                               update_f_calc  = True)
  assert approx_equal(fmodel.r_work(), 0.12, 1.e-2)
  assert approx_equal(fmodel.r_free(), 0.13, 1.e-2)
  rb = mmtbx.refinement.rigid_body.manager(fmodel           = fmodel,
                                           selections       = selections,
                                           refine_r         = True,
                                           refine_t         = False,
                                           convergence_test = True)
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
  new_xray_structure = mmtbx.refinement.rigid_body.apply_transformation(
                         xray_structure      = model.xray_structure,
                         rotation_matrices   = [rot_obj.rot_mat()],
                         translation_vectors = [(1,1,1)],
                         selections          = selections)
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
                                           convergence_test = True)
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

  rot_obj_1 = mmtbx.refinement.rigid_body.rb_mat(phi = 0.5,
                                                 psi = 0.5,
                                                 the = 0.5)
  rot_obj_2 = mmtbx.refinement.rigid_body.rb_mat(phi = 0.3,
                                                 psi = 0.3,
                                                 the = 0.3)
  new_xray_structure = mmtbx.refinement.rigid_body.apply_transformation(
                         xray_structure      = model.xray_structure,
                         rotation_matrices   = [rot_obj_1.rot_mat(),rot_obj_2.rot_mat()],
                         translation_vectors = [(1,1,1),(1.2,1.2,1.2)],
                         selections          = selections)
  fmodel.update_xray_structure(xray_structure = new_xray_structure,
                               update_f_calc  = True)
  assert approx_equal(fmodel.r_work(), 0.52, 1.e-2)
  assert approx_equal(fmodel.r_free(), 0.53, 1.e-2)
  rb = mmtbx.refinement.rigid_body.manager(fmodel           = fmodel,
                                           selections       = selections,
                                           refine_r         = True,
                                           refine_t         = True,
                                           convergence_test = True)
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
  model = mmtbx.model.manager(processed_pdb_file = processed_pdb_file)
  model.xray_structure.scattering_dict(table = "wk1995")
  model.setup_restraints_manager()
  dummy = model.xray_structure.structure_factors(algorithm = sf_algorithm,
                                                 d_min     = 2.0).f_calc()
  f_obs = abs(dummy.structure_factors_from_scatterers(
                                         xray_structure = model.xray_structure,
                                         algorithm      = sf_algorithm,
                                         cos_sin_table  = True).f_calc())
  flags =f_obs.array(data=flex.size_t(xrange(1,f_obs.data().size()+1))%10 == 0)
  #fout = open("test.hkl","w")
  #for i,d,fl in zip(f_obs.indices(),f_obs.data(), flags.data()):
  #  print >> fout, (" INDE %5d%5d%5d FOBS=%10.4f SIGMA=%10.1f TEST=%10d") %\
  #           (i[0],i[1],i[2],d,1.0, fl)
  fmodel = mmtbx.f_model.manager(xray_structure    = model.xray_structure,
                                 f_obs             = f_obs,
                                 r_free_flags      = flags,
                                 target_name       = "ls_wunit_k1",
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
  ##############################################################################
  ##rot_obj = mmtbx.refinement.rigid_body.rb_mat(phi = 0,
  ##                                             psi = 0,
  ##                                             the = 0)
  ##fmodel_copy = fmodel.deep_copy()
  ##fmodel_copy = fmodel_copy.resolution_filter(d_min=10.0, d_max= 9999.0)
  ##print fmodel_copy.f_obs.data().size()
  ##dim = fmodel_copy.xray_structure.scatterers().size()
  ##tg_obj = mmtbx.refinement.rigid_body.target_and_grads(
  ##                                         fmodel     = fmodel_copy,
  ##                                         rot_objs   = [rot_obj],
  ##                                         selections = [flex.bool(dim, True)])
  ##f = tg_obj.target()
  ##g1, g2 = tg_obj.gradients_wrt_r(), tg_obj.gradients_wrt_t()
  ##print list(g2[0])
  ##rg = [i/1000000. for i in range(1,10500)] #[i/100000. for i in range(1,1050)]
  ##sum = 99999.999
  ##for i in rg:
  ##    qq = fd(fmodel_copy, e = i)
  ##    tmp = flex.sum(flex.abs(flex.abs(g2[0]) - flex.abs(flex.double(qq))))
  ##    if( tmp < sum ):
  ##        sum = tmp
  ##        print qq, i, sum
  ##assert 0
  ##############################################################################


def fd(fmodel, e):
  xrs1 = fmodel.xray_structure.deep_copy_scatterers()
  xrs2 = fmodel.xray_structure.deep_copy_scatterers()
  fm = fmodel.deep_copy()
  xrs_g1_1 = xrs1.translate(x = e)
  xrs_g1_2 = xrs2.translate(x =-e)
  fm.update_xray_structure(xray_structure = xrs_g1_1,
                           update_f_calc  = True)
  t1 = fm.target_w()
  fm.update_xray_structure(xray_structure = xrs_g1_2,
                           update_f_calc  = True)
  t2 = fm.target_w()
  gtx = (t1-t2)/(2*e)

  xrs = fmodel.xray_structure.deep_copy_scatterers()
  fm = fmodel.deep_copy()
  xrs_g1_1 = xrs.translate(y = e)
  xrs_g1_2 = xrs.translate(y =-e)
  fm.update_xray_structure(xray_structure = xrs_g1_1,
                           update_f_calc  = True)
  t1 = fm.target_w()
  fm.update_xray_structure(xray_structure = xrs_g1_2,
                           update_f_calc  = True)
  t2 = fm.target_w()
  gty = (t1-t2)/(2*e)

  xrs = fmodel.xray_structure.deep_copy_scatterers()
  fm = fmodel.deep_copy()
  xrs_g1_1 = xrs.translate(z = e)
  xrs_g1_2 = xrs.translate(z =-e)
  fm.update_xray_structure(xray_structure = xrs_g1_1,
                           update_f_calc  = True)
  t1 = fm.target_w()
  fm.update_xray_structure(xray_structure = xrs_g1_2,
                           update_f_calc  = True)
  t2 = fm.target_w()
  gtz = (t1-t2)/(2*e)
  return gtx,gty,gtz

def exercise():
  test_matrices()
  run_tests()
  print format_cpu_times()

if (__name__ == "__main__"):
  exercise()
