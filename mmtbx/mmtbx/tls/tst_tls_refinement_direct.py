from iotbx import pdb
import iotbx.pdb.interpretation
from cctbx.array_family import flex
import sys, math, os
from mmtbx.tls.tls import *
from libtbx.test_utils import not_approx_equal
from libtbx.utils import format_cpu_times
import libtbx.load_env
from cctbx import adptbx
from cctbx import miller
from mmtbx.tls.tls import *
from mmtbx_tls_ext import *
import iotbx.pdb.remark_3_interpretation
from libtbx.test_utils import approx_equal

def run(n_macro_cycles):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/phe_abc_tlsanl_out.pdb", test=os.path.isfile)
  if (pdb_file is None):
    print "Skipping run(): input file not available"
    return
  xray_structure = pdb.input(file_name=pdb_file).xray_structure_simple()
  stage_1 = pdb.interpretation.stage_1(file_name   = pdb_file)
  xray_structure.show_summary()

  if(1):
    print "="*80
    print "INPUT TLS :"
    tls_data = iotbx.pdb.remark_3_interpretation.extract_tls_parameters(stage_1.remark_3_records)

    tls_data.print_as_matrices()
    tls_obj = tls(tls_parameters = tls_data,
                  xray_structure = xray_structure,
                  stage_1        = stage_1)
    u_from_tls = tls_obj.u_from_tls()
    u_from_pdb = xray_structure.scatterers().extract_u_cart(xray_structure.unit_cell())

    print "="*80
    print "INPUT U_ANISO MODIFIED BY CCTBX:"
    for i in xrange(len(u_from_tls)):
      if( not_approx_equal(u_from_tls[i],u_from_pdb[i], 1e-3) ):
        print "tls=", u_from_tls[i]
        print "pdb=", u_from_pdb[i]

    print "="*80
    print "Fcalc calculation:"
    xray_structure.scatterers().set_u_cart(xray_structure.unit_cell(),
                                           flex.sym_mat3_double(u_from_tls))
    xray_structure.apply_symmetry_u_stars()
    fcalc_cctbx = xray_structure.structure_factors(d_min = 2.0,
                                                   algorithm = "direct").f_calc()
    fcalc_cctbx.show_summary()
    print "n_ref = ", fcalc_cctbx.data().size()

    print "="*80
    print "TLS REFINEMENT : "
    T_initial = []
    L_initial = []
    S_initial = []
    if(0):
      for i in xrange(tls_data.n_groups):
        T_initial.append([1.0,1.0,1.0,1.0,1.0,1.0])
        L_initial.append([1.0,1.0,1.0,1.0,1.0,1.0])
        S_initial.append([1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0])
    if(0):
      for i in xrange(tls_data.n_groups):
        T_initial.append([0.0,0.0,0.0,0.0,0.0,0.0])
        L_initial.append([0.0,0.0,0.0,0.0,0.0,0.0])
        S_initial.append([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    if(0):
        T_initial.append([2.0,2.0,2.0,2.0,2.0,2.0])
        #T_initial.append([0.11,0.22,0.33,0.12,0.13,0.23])
        L_initial.append([1.11,1.22,1.33,1.12,1.13,1.23])
        S_initial.append([0.11,0.12,0.13,0.21,0.22,0.23,0.31,0.32,-0.33])

        T_initial.append([1.0,1.0,1.0,1.0,1.0,1.0])
        #T_initial.append([0.22,0.44,0.66,0.24,0.26,0.46])
        L_initial.append([2.22,2.44,2.66,2.24,2.26,2.46])
        S_initial.append([0.22,0.24,0.26,0.42,0.44,0.46,0.62,0.64,-0.66])

        T_initial.append([1.0,1.0,1.0,1.0,1.0,1.0])
        #T_initial.append([0.33,0.66,0.99,0.36,0.39,0.69])
        L_initial.append([2.33,2.66,2.99,2.36,2.39,2.69])
        S_initial.append([0.22,0.24,0.26,0.42,0.44,0.46,0.62,0.64,-0.66])
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
# refine TLS:
    if(1):
        T_initial.append([0.10,0.20,0.30,0.10,0.10,0.20])
        #L_initial.append([1.10,1.20,1.30,1.10,1.10,1.20])
        L_initial.append([0.0,0.0,0.0,0.0,0.0,0.0])
        S_initial.append([0.10,0.10,0.10,0.20,0.20,0.20,0.30,0.30,-0.30])

        T_initial.append([0.20,0.40,0.60,0.20,0.20,0.40])
        #L_initial.append([2.20,2.40,2.60,2.20,2.20,2.40])
        L_initial.append([0.0,0.0,0.0,0.0,0.0,0.0])
        S_initial.append([0.20,0.20,0.20,0.40,0.40,0.40,0.60,0.60,-0.60])

        T_initial.append([0.30,0.60,0.90,0.30,0.30,0.60])
        #L_initial.append([2.30,2.60,2.90,2.30,2.30,2.60])
        L_initial.append([0.0,0.0,0.0,0.0,0.0,0.0])
        S_initial.append([0.20,0.20,0.20,0.40,0.40,0.40,0.60,0.60,-0.60])


    assert tls_data.n_groups == len(tls_obj.xray_structures)
    for i in xrange(n_macro_cycles):
      flag = False
      i=0
      #while not flag:
      i += 1
      print "========================================================"
      print "MACRO_CYCLE = ", i
      f1 = approx_equal(T_initial, tls_data.T, 1.e-2, out=None)
      f2 = approx_equal(L_initial, tls_data.L, 1.e-2, out=None)
      f3 = approx_equal(S_initial, tls_data.S, 1.e-2, out=None)
      flag = f1 and f2 and f3
      print f1, f2, f3
      minimized = tls_xray_target_minimizer(
           fcalc_cctbx,
           T_initial = T_initial,
           L_initial = L_initial,
           S_initial = S_initial,
           refine_T = 1,
           refine_L = 1,
           refine_S = 1,
           origins = tls_obj.origin,
           xray_structures = tls_obj.xray_structures)
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
        print (format1+format2)% (S[0],S[1],S[2],S[3],S[4],S[5],S[6],S[7],S[8])

if (__name__ == "__main__"):
  if (len(sys.argv) == 1):
    n_macro_cycles = 2
  else:
    n_macro_cycles = int(sys.argv[1])
  run(n_macro_cycles=n_macro_cycles)
  print format_cpu_times()
