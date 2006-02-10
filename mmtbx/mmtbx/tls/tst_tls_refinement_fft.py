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

from mmtbx.tls import tools
from mmtbx_tls_ext import *


def run():
###> Get start from PDB
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="regression/pdb/phe_abc_tlsanl_out.pdb", test=os.path.isfile)
  #pdb_file = libtbx.env.find_in_repositories(
  #            relative_path="regression/pdb/1OC2_tst.pdb", test=os.path.isfile)
  processed_pdb_file = monomer_library.pdb_interpretation.process(
                                       mon_lib_srv               = mon_lib_srv,
                                       ener_lib                  = ener_lib,
                                       file_name                 = pdb_file,
                                       raw_records               = None,
                                       force_symmetry            = True)
  xray_structure = processed_pdb_file.xray_structure()
  stage_1 = processed_pdb_file.all_chain_proxies.stage_1
  selections = []
  for string in ["chain A", "chain B", "chain C"]:
      selections.append(processed_pdb_file.all_chain_proxies.selection(
                                                              string = string))
  input_tls_data = iotbx.pdb.remark_3_interpretation.extract_tls_parameters(
                                                      stage_1.remark_3_records)
  tls_params = []
  origins = []
  for item in input_tls_data:
      origins.append(item.origin)
      tls_params.append(tools.tlso(t      = item.T,
                                   l      = item.L,
                                   s      = item.S,
                                   origin = item.origin))
  tools.show_tls(tlsos = tls_params)
###> Set up fmodel
  dummy = xray_structure.structure_factors(algorithm = "direct",
                                           d_min     = 1.5).f_calc()
  f_obs = abs(dummy.structure_factors_from_scatterers(
                                         xray_structure = xray_structure,
                                         algorithm      = "direct",
                                         cos_sin_table  = True).f_calc())
  flags =f_obs.array(data=flex.size_t(xrange(1,f_obs.data().size()+1))%3 == 0)
  fmodel = mmtbx.f_model.manager(xray_structure    = xray_structure,
                                 f_obs             = f_obs,
                                 r_free_flags      = flags,
                                 target_name       = "ls_wunit_k1",
                                 sf_algorithm      = "direct")
  fmodel.show_comprehensive(reflections_per_bin = 250,
                            max_number_of_bins  = 30)

  print "="*80
# refine only T
  T_initial = []
  L_initial = []
  S_initial = []
  T_initial.append([2.0,2.0,2.0,2.0,2.0,2.0])
  L_initial.append([1.11,1.22,1.33,1.12,1.13,1.23])
  S_initial.append([0.11,0.12,0.13,0.21,0.22,0.23,0.31,0.32,-0.33])

  T_initial.append([1.0,1.0,1.0,1.0,1.0,1.0])
  L_initial.append([2.22,2.44,2.66,2.24,2.26,2.46])
  S_initial.append([0.22,0.24,0.26,0.42,0.44,0.46,0.62,0.64,-0.66])

  T_initial.append([1.0,1.0,1.0,1.0,1.0,1.0])
  L_initial.append([2.33,2.66,2.99,2.36,2.39,2.69])
  S_initial.append([0.22,0.24,0.26,0.42,0.44,0.46,0.62,0.64,-0.66])

  tlsosT = tools.generate_tlsos(selections     = selections,
                                xray_structure = fmodel.xray_structure,
                                T              = T_initial,
                                L              = L_initial,
                                S              = S_initial)
# refine only L
  T_initial = []
  L_initial = []
  S_initial = []
  T_initial.append([0.11,0.22,0.33,0.12,0.13,0.23])
  L_initial.append([100.0,100.0,100.0,100.0,100.0,100.0])
  S_initial.append([0.11,0.12,0.13,0.21,0.22,0.23,0.31,0.32,-0.33])

  T_initial.append([0.22,0.44,0.66,0.24,0.26,0.46])
  L_initial.append([-100.0,100.0,-100.0,100.0,-100.0,100.0])
  S_initial.append([0.22,0.24,0.26,0.42,0.44,0.46,0.62,0.64,-0.66])

  T_initial.append([0.33,0.66,0.99,0.36,0.39,0.69])
  L_initial.append([100.0,-100.0,100.0,-100.0,100.0,-100.0])
  S_initial.append([0.22,0.24,0.26,0.42,0.44,0.46,0.62,0.64,-0.66])

  tlsosL = tools.generate_tlsos(selections     = selections,
                                xray_structure = fmodel.xray_structure,
                                T              = T_initial,
                                L              = L_initial,
                                S              = S_initial)
# refine only S
  T_initial = []
  L_initial = []
  S_initial = []
  T_initial.append([0.11,0.22,0.33,0.12,0.13,0.23])
  L_initial.append([1.11,1.22,1.33,1.12,1.13,1.23])
  S_initial.append([9.0,9.0,9.0,9.0,9.0,9.0,9.0,9.0,9.0])

  T_initial.append([0.22,0.44,0.66,0.24,0.26,0.46])
  L_initial.append([2.22,2.44,2.66,2.24,2.26,2.46])
  S_initial.append([9.0,9.0,9.0,9.0,9.0,9.0,9.0,9.0,9.0])

  T_initial.append([0.33,0.66,0.99,0.36,0.39,0.69])
  L_initial.append([2.33,2.66,2.99,2.36,2.39,2.69])
  S_initial.append([10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0])

  tlsosS = tools.generate_tlsos(selections     = selections,
                                xray_structure = fmodel.xray_structure,
                                T              = T_initial,
                                L              = L_initial,
                                S              = S_initial)
# Answer A
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
                                xray_structure = fmodel.xray_structure,
                                T              = T_initial,
                                L              = L_initial,
                                S              = S_initial)


  for set in ([1,0,0,tlsosT,"T"],[0,1,0,tlsosL,"L"],[0,0,1,tlsosS,"S"]):
      tls_refinement_manager = tools.tls_refinement(
                           fmodel                   = fmodel,
                           selections               = selections,
                           refine_T                 = set[0],
                           refine_L                 = set[1],
                           refine_S                 = set[2],
                           number_of_macro_cycles   = 10,
                           max_number_of_iterations = 50,
                           start_tls_value          = set[3])
      tlsos = tls_refinement_manager.tlsos
      if(set[4] == "T"):
         for i1,i2 in zip(tlsos, tlsosA):
             approx_equal(i1.t, i2.t, 0.01)
             approx_equal(i1.l, i2.l)
             approx_equal(i1.s, i2.s)
      if(set[4] == "L"):
         for i1,i2 in zip(tlsos, tlsosA):
             approx_equal(i1.t, i2.t)
             approx_equal(i1.l, i2.l, 0.45)
             approx_equal(i1.s, i2.s)
      if(set[4] == "S"):
         for i1,i2 in zip(tlsos, tlsosA):
             approx_equal(i1.t, i2.t)
             approx_equal(i1.l, i2.l)
             approx_equal(i1.s, i2.s, 0.02)

if (__name__ == "__main__"):
  run()
  print format_cpu_times()
