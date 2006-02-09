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


def run(n_macro_cycles):
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
  print "TLS REFINEMENT : "
  T_initial = []
  L_initial = []
  S_initial = []
  if(0):
    for i in tls_params:
      T_initial.append([1.0,1.0,1.0,1.0,1.0,1.0])
      L_initial.append([1.0,1.0,1.0,1.0,1.0,1.0])
      S_initial.append([1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0])
  if(1):
    for i in tls_params:
      T_initial.append([0.0,0.0,0.0,0.0,0.0,0.0])
      L_initial.append([0.0,0.0,0.0,0.0,0.0,0.0])
      S_initial.append([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])

# refine only T
  if(0):
        T_initial.append([2.0,2.0,2.0,2.0,2.0,2.0])
        L_initial.append([1.11,1.22,1.33,1.12,1.13,1.23])
        S_initial.append([0.11,0.12,0.13,0.21,0.22,0.23,0.31,0.32,-0.33])

        T_initial.append([1.0,1.0,1.0,1.0,1.0,1.0])
        L_initial.append([2.22,2.44,2.66,2.24,2.26,2.46])
        S_initial.append([0.22,0.24,0.26,0.42,0.44,0.46,0.62,0.64,-0.66])

        T_initial.append([1.0,1.0,1.0,1.0,1.0,1.0])
        L_initial.append([2.33,2.66,2.99,2.36,2.39,2.69])
        S_initial.append([0.22,0.24,0.26,0.42,0.44,0.46,0.62,0.64,-0.66])
# refine only L
  if(0):
        T_initial.append([0.11,0.22,0.33,0.12,0.13,0.23])
        L_initial.append([130.0,130.0,130.0,130.0,130.0,130.0])
        S_initial.append([0.11,0.12,0.13,0.21,0.22,0.23,0.31,0.32,-0.33])

        T_initial.append([0.22,0.44,0.66,0.24,0.26,0.46])
        L_initial.append([-130.0,130.0,-130.0,130.0,-130.0,130.0])
        S_initial.append([0.22,0.24,0.26,0.42,0.44,0.46,0.62,0.64,-0.66])

        T_initial.append([0.33,0.66,0.99,0.36,0.39,0.69])
        L_initial.append([130.0,-130.0,130.0,-130.0,130.0,-130.0])
        S_initial.append([0.22,0.24,0.26,0.42,0.44,0.46,0.62,0.64,-0.66])
# refine only S
  if(0):
        T_initial.append([0.11,0.22,0.33,0.12,0.13,0.23])
        L_initial.append([1.11,1.22,1.33,1.12,1.13,1.23])
        S_initial.append([9.0,9.0,9.0,9.0,9.0,9.0,9.0,9.0,9.0])

        T_initial.append([0.22,0.44,0.66,0.24,0.26,0.46])
        L_initial.append([2.22,2.44,2.66,2.24,2.26,2.46])
        S_initial.append([9.0,9.0,9.0,9.0,9.0,9.0,9.0,9.0,9.0])

        T_initial.append([0.33,0.66,0.99,0.36,0.39,0.69])
        L_initial.append([2.33,2.66,2.99,2.36,2.39,2.69])
        S_initial.append([10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0])


  ipd = fmodel.xray_structure.is_positive_definite_u()
  print "fmodel.xray_structure.is_positive_definite_u() = ", ipd.count(True), ipd.count(False)


  tools.tls_refinement(fmodel                   = fmodel,
                       selections               = selections,
                       refine_T                 = 1,
                       refine_L                 = 1,
                       refine_S                 = 1,
                       number_of_macro_cycles   = 5,
                       max_number_of_iterations = 10)

  #for i in xrange(n_macro_cycles):
  for (i,j,k) in [(1,1,1),]*5:
  #for (i,j,k) in [(1,0,0),(0,0,1),(1,0,1),(1,1,1)]*3 + [(0,1,0),(1,1,1)]*3:
  #for (i,j,k) in [(1,0,0),(0,0,1),(1,0,1),(1,1,1), (1,0,0)]*3+[(0,1,0),(1,1,0),(1,1,1)]*3+[(0,0,1),(0,1,0),(0,1,1),(1,1,1)]*3:
  #for (i,j,k) in [(1,0,0),(0,0,1),(1,0,1),(1,1,1)]*10:
    flag = False
    print "========================================================"
    #print "MACRO_CYCLE = ", i
    print "MACRO_CYCLE = ",i,j,k
    fmodel.xray_structure.tidy_us(u_min = 1.e-6)
    minimized = tools.tls_xray_target_minimizer(
                                          fmodel    = fmodel,
                                          T_initial = T_initial,
                                          L_initial = L_initial,
                                          S_initial = S_initial,
                                          refine_T  = 1,
                                          refine_L  = 1,
                                          refine_S  = 1,
                                          selections= selections,
                                          origins   = origins)
    fmodel.xray_structure.tidy_us(u_min = 1.e-6)
    T_initial = minimized.T_min
    L_initial = minimized.L_min
    S_initial = minimized.S_min
    T_ = minimized.T_min
    L_ = minimized.L_min
    S_ = minimized.S_min
    print "target = ", minimized.f
    for T,L,S in zip(T_,L_,S_):
      print
      print "T11=%8.4f T22=%8.4f T33=%8.4f T12=%8.4f T13=%8.4f T23=%8.4f"%\
             (T[0],T[1],T[2],T[3],T[4],T[5])
      print "L11=%8.4f L22=%8.4f L33=%8.4f L12=%8.4f L13=%8.4f L23=%8.4f"%\
             (L[0],L[1],L[2],L[3],L[4],L[5])
      format1 = "S11=%8.4f S22=%8.4f S33=%8.4f S12=%8.4f S13=%8.4f "
      format2 = "S23=%8.4f S21=%8.4f S31=%8.4f S32=%8.4f"
      print (format1+format2)% (S[0],S[4],S[8],S[1],S[2],S[5],S[3],S[6],S[7])

if (__name__ == "__main__"):
  if (len(sys.argv) == 1):
    n_macro_cycles = 10
  else:
    n_macro_cycles = int(sys.argv[1])
  run(n_macro_cycles=n_macro_cycles)
  print format_cpu_times()
