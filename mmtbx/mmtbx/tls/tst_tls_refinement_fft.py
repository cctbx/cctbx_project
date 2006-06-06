from iotbx import pdb
import iotbx.pdb.interpretation
import iotbx.pdb.remark_3_interpretation
from cctbx.array_family import flex
from libtbx.utils import format_cpu_times
import sys, math, time, os
from libtbx.test_utils import approx_equal
import libtbx.load_env
from cctbx import adptbx
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
import mmtbx.f_model
import mmtbx.model


from mmtbx.tls import tools
from mmtbx_tls_ext import *


def exercise_2(eps = 1.e-6):
###> Get started from PDB
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="regression/pdb/phe_abc_tlsanl_out_geometry_minimized.pdb",
    test=os.path.isfile)
  processed_pdb_file = monomer_library.pdb_interpretation.process(
                                       mon_lib_srv               = mon_lib_srv,
                                       ener_lib                  = ener_lib,
                                       file_name                 = pdb_file,
                                       raw_records               = None,
                                       force_symmetry            = True)
  xray_structure = processed_pdb_file.xray_structure()
  xray_structure.scattering_type_registry(table = "wk1995")
  xray_structure.convert_to_isotropic()
  u_iso_start = xray_structure.extract_u_iso_or_u_equiv()
  xray_structure.convert_to_anisotropic()
  stage_1 = processed_pdb_file.all_chain_proxies.stage_1
  selections = []
  for string in ["chain A", "chain B", "chain C"]:
      selections.append(processed_pdb_file.all_chain_proxies.selection(
                                                              string = string))
###> Get TLS <-> Ucart
  T_initial = []
  L_initial = []
  S_initial = []
  T_initial.append([0.11,0.22,0.33,0.12,0.13,0.23])
  L_initial.append([1.11,1.22,1.33,1.12,1.13,1.23])
  S_initial.append([0.11,0.12,0.13,0.21,0.22,0.23,0.31,0.32,-0.33])

  T_initial.append([0.22,0.44,0.66,0.24,0.26,0.46])
  L_initial.append([2.22,2.44,2.66,2.24,2.26,2.46])
  S_initial.append([0.22,0.24,0.26,0.42,0.44,0.46,0.62,0.64,-0.66])

  T_initial.append([0.33,0.66,0.99,0.36,0.39,0.69])
  L_initial.append([2.33,2.66,2.99,2.36,2.39,2.69])
  S_initial.append([0.22,0.24,0.26,0.42,0.44,0.46,0.62,0.64,-0.66])

  tlsosA = tools.generate_tlsos(selections     = selections,
                                xray_structure = xray_structure,
                                T              = T_initial,
                                L              = L_initial,
                                S              = S_initial)

  tlsos = tools.generate_tlsos(selections     = selections,
                               xray_structure = xray_structure,
                               T              = T_initial,
                               L              = L_initial,
                               S              = S_initial)
  tlsos = tools.make_tlso_compatible_with_u_positive_definite(
                  tlsos                                       = tlsos,
                  xray_structure                              = xray_structure.deep_copy_scatterers(),
                  selections                                  = selections,
                  max_iterations                              = 50,
                  number_of_u_nonpositive_definite            = 0,
                  eps                                         = eps,
                  number_of_macro_cycles_for_tls_from_uanisos = 30)

  u_cart_answer = tools.uanisos_from_tls(sites_cart = xray_structure.sites_cart(),
                                         selections = selections,
                                         tlsos      = tlsos)
  xray_structure.scatterers().set_u_cart(xray_structure.unit_cell(), u_cart_answer)

  assert approx_equal(u_cart_answer,
        xray_structure.scatterers().extract_u_cart(xray_structure.unit_cell()))


  tools.show_tls(tlsos = tlsos, text = "ANSWER")

###> Set up fmodel
  sf_algorithm = "direct"
  dummy = xray_structure.structure_factors(algorithm = sf_algorithm,
                                           d_min     = 1.5).f_calc()
  f_obs = abs(dummy.structure_factors_from_scatterers(
                                         xray_structure = xray_structure,
                                         algorithm      = sf_algorithm,
                                         cos_sin_table  = False).f_calc())
  flags = f_obs.generate_r_free_flags(fraction=0.05, max_free=2000)
  fmodel = mmtbx.f_model.manager(xray_structure    = xray_structure,
                                 f_obs             = f_obs,
                                 r_free_flags      = flags,
                                 target_name       = "ls_wunit_k1",
                                 sf_algorithm      = sf_algorithm,
                                 sf_cos_sin_table  = False)
  fmodel.show_comprehensive(reflections_per_bin = 250,
                            max_number_of_bins  = 30)
  xray_structure.convert_to_isotropic()
  xray_structure.set_b_iso(value = 25.0)
  fmodel.update_xray_structure(xray_structure = xray_structure,
                               update_f_calc  = True)
  fmodel.show_comprehensive(reflections_per_bin = 250,
                            max_number_of_bins  = 30)
  print "*"*80
###> TLS refinement against xray data
  if (not "--comprehensive" in sys.argv[1:]):
          number_of_macro_cycles   = 1
          max_number_of_iterations = 3
  else:
          number_of_macro_cycles   = 100
          max_number_of_iterations = 50

  for start_tls_value in [0.0, tlsosA, None]:
      print " \n "+str(start_tls_value) + " \n "
      fmodel_cp = fmodel.deep_copy()
      #for sc in fmodel_cp.xray_structure.scatterers():
      #  sc.flags.set_use_u_aniso(True)
      fmodel_cp.xray_structure.convert_to_anisotropic()

      if(start_tls_value is None):
         run_finite_differences_test = True
      else: run_finite_differences_test = False
      tls_refinement_manager = tools.tls_refinement(
                     fmodel                      = fmodel_cp,
                     selections                  = selections,
                     refine_T                    = 1,
                     refine_L                    = 1,
                     refine_S                    = 1,
                     number_of_macro_cycles      = number_of_macro_cycles,
                     max_number_of_iterations    = max_number_of_iterations,
                     start_tls_value             = start_tls_value,
                     run_finite_differences_test = run_finite_differences_test,
                     eps                         = eps)
      u_cart = tls_refinement_manager.fmodel.xray_structure.scatterers().extract_u_cart(
                                                        xray_structure.unit_cell())
      if("--comprehensive" in sys.argv[1:]):
         format   = "%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f"
         counter = 0
         if(start_tls_value == tlsosA): tolerance = 1.e-6
         else: tolerance = 0.02
         for m1,m2 in zip(u_cart_answer, u_cart):
             counter += 1
             if(counter < 10):
                print "1=" + format % (m1[0],m1[1],m1[2],m1[3],m1[4],m1[5])
                print "2=" + format % (m2[0],m2[1],m2[2],m2[3],m2[4],m2[5])
             assert approx_equal(m1,m2, tolerance)



if (__name__ == "__main__"):
  exercise_2()
  print format_cpu_times()
