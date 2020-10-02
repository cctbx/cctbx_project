from __future__ import absolute_import, division, print_function
from iotbx import pdb
from cctbx.array_family import flex
import sys, time
from scitbx import lbfgs
from mmtbx_tls_ext import *
from libtbx import adopt_init_args
from libtbx.test_utils import approx_equal
from cctbx import adptbx
from cctbx import xray
from libtbx.utils import user_plus_sys_time, Sorry
from libtbx.str_utils import line_breaker
from libtbx import group_args
import mmtbx.utils
import iotbx
import math
import iotbx.pdb.remark_3_interpretation
from scitbx.linalg import eigensystem
from six.moves import zip
from six.moves import range
from scitbx_array_family_flex_ext import reindexing_array

def combine_tls_and_u_local(xray_structure, tls_selections, tls_groups):
  assert len(tls_selections) == len(tls_groups)
  for sel in tls_selections:
    xray_structure.convert_to_anisotropic(selection = sel)
  tlsos = []
  for tls_group in tls_groups:
    tlsos.append(tlso(t = tls_group.t,
                      l = tls_group.l,
                      s = tls_group.s,
                      origin = tls_group.origin))
  u_cart_tls = u_cart_from_tls(
    sites_cart = xray_structure.sites_cart(),
    selections = tls_selections,
    tlsos      = tlsos)
  unit_cell = xray_structure.unit_cell()
  for i_seq, sc in enumerate(xray_structure.scatterers()):
    if(u_cart_tls[i_seq] != (0,0,0,0,0,0)):
      assert sc.flags.use_u_aniso()
      u_star_tls = adptbx.u_cart_as_u_star(unit_cell,
        tuple(u_cart_tls[i_seq]))
      sc.u_star = tuple(flex.double(sc.u_star) + flex.double(u_star_tls))

def tls_from_pdb_inp(remark_3_records, pdb_hierarchy):
  chain_ids = []
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      chain_ids.append(chain.id)
  return iotbx.pdb.remark_3_interpretation.extract_tls_parameters(
    remark_3_records = remark_3_records,
    pdb_hierarchy    = pdb_hierarchy,
    chain_ids        = chain_ids)

def extract_tls_from_pdb(pdb_inp_tls, model):
  if(len(pdb_inp_tls.tls_params)>0):
    tls_selection_strings = []
    for i_seq, tls_group in enumerate(pdb_inp_tls.tls_params):
      tls_selection_strings.append(tls_group.selection_string)
    try:
      selections = mmtbx.utils.get_atom_selections(
        model = model,
        selection_strings = tls_selection_strings)
      return group_args(pdb_inp_tls           = pdb_inp_tls,
                        tls_selections        = selections,
                        tls_selection_strings = tls_selection_strings)
    except Sorry:
      return group_args(pdb_inp_tls           = pdb_inp_tls,
                        tls_selections        = [],
                        tls_selection_strings = [])
  else:
    return group_args(pdb_inp_tls           = pdb_inp_tls,
                      tls_selections        = [],
                      tls_selection_strings = [])

class tls_group(object):
  def __init__(self, tlso, selection_string = None, selection_array = None):
    self.tlso = tlso
    self.selection_string = selection_string
    self.selection_array = selection_array

time_u_cart_from_tls                              = 0.0
time_tls_from_uanisos                              = 0.0
time_update_xray_structure_with_tls                = 0.0
time_split_u                                       = 0.0
time_tls_from_u_cart                               = 0.0
time_make_tlso_compatible_with_u_positive_definite = 0.0
time_generate_tlsos                                = 0.0
time_tls_total                                     = 0.0

def show_times(out = None):
  if(out is None): out = sys.stdout
  total = time_u_cart_from_tls                            +\
        time_tls_from_uanisos                              +\
        time_update_xray_structure_with_tls                +\
        time_split_u                                       +\
        time_tls_from_u_cart                               +\
        time_make_tlso_compatible_with_u_positive_definite +\
        time_generate_tlsos
  if(total > 0.01):
    print("TLS refinement:", file=out)
    print("  time_u_cart_from_tls                              = %-7.2f" % time_u_cart_from_tls, file=out)
    print("  time_tls_from_uanisos                              = %-7.2f" % time_tls_from_uanisos, file=out)
    print("  time_update_xray_structure_with_tls                = %-7.2f" % time_update_xray_structure_with_tls, file=out)
    print("  time_split_u                                       = %-7.2f" % time_split_u, file=out)
    print("  time_tls_from_u_cart                               = %-7.2f" % time_tls_from_u_cart, file=out)
    print("  time_make_tlso_compatible_with_u_positive_definite = %-7.2f" % time_make_tlso_compatible_with_u_positive_definite, file=out)
    print("  time_generate_tlsos                                = %-7.2f" % time_generate_tlsos, file=out)
    print("  sum_of_partial_contributions                       = %-7.2f" % total, file=out)
    print("  time_tls_total                                     = %-7.2f" % time_tls_total, file=out)
  return total

class tls_groups(object):
  def __init__(self, tlsos = None, selection_strings = None, iselections=None):
    self.tlsos = tlsos
    self.selection_strings = selection_strings
    # These selections are used to produce mmCIF file. They are not
    # used in refinement
    self.iselections = iselections

  def select(self, selection, n_atoms):
    # only select self.iselections
    selected_isels = []
    if self.iselections is not None:
      for isel in self.iselections:
        ra = reindexing_array(n_atoms, selection.iselection().as_int())
        bsel = flex.bool(n_atoms, isel)
        selected_isel = (bsel&selection).iselection()
        new_isel = flex.size_t([ra[i] for i in selected_isel ])
        selected_isels.append(new_isel)
    else:
      selected_isels = None
    return tls_groups(
        tlsos=self.tlsos,
        selection_strings=self.selection_strings,
        iselections=selected_isels)

  def as_cif_block(self, hierarchy, cif_block=None, scattering_type=None):
    import iotbx.cif.model
    if cif_block is None:
      cif_block = iotbx.cif.model.block()

    if (len(self.selection_strings) == 0):
      assert len(self.tlsos) == 1
      self.selection_strings = ["all"]
    else:
      assert len(self.tlsos) == len(self.selection_strings), "%d, %d" % (
            len(self.tlsos), len(self.selection_strings))
    if self.iselections is None or len(self.selection_strings) != len(self.iselections):
      # unfortunate event we have to regenerate iselections
      asc = hierarchy.atom_selection_cache()
      self.iselections = []
      for str_sel in self.selection_strings:
        self.iselections.append(asc.iselection(str_sel))

    if scattering_type is None:
      scattering_type = '?'

    tls_loop = iotbx.cif.model.loop(header=(
      "_pdbx_refine_tls.id",
      "_pdbx_refine_tls.details",
      "_pdbx_refine_tls.pdbx_refine_id",
      "_pdbx_refine_tls.method",
      "_pdbx_refine_tls.origin_x",
      "_pdbx_refine_tls.origin_y",
      "_pdbx_refine_tls.origin_z",
      "_pdbx_refine_tls.T[1][1]",
      "_pdbx_refine_tls.T[2][2]",
      "_pdbx_refine_tls.T[3][3]",
      "_pdbx_refine_tls.T[1][2]",
      "_pdbx_refine_tls.T[1][3]",
      "_pdbx_refine_tls.T[2][3]",
      "_pdbx_refine_tls.L[1][1]",
      "_pdbx_refine_tls.L[2][2]",
      "_pdbx_refine_tls.L[3][3]",
      "_pdbx_refine_tls.L[1][2]",
      "_pdbx_refine_tls.L[1][3]",
      "_pdbx_refine_tls.L[2][3]",
      "_pdbx_refine_tls.S[1][1]",
      "_pdbx_refine_tls.S[1][2]",
      "_pdbx_refine_tls.S[1][3]",
      "_pdbx_refine_tls.S[2][1]",
      "_pdbx_refine_tls.S[2][2]",
      "_pdbx_refine_tls.S[2][3]",
      "_pdbx_refine_tls.S[3][1]",
      "_pdbx_refine_tls.S[3][2]",
      "_pdbx_refine_tls.S[3][3]",
    ))

    tls_group_loop = iotbx.cif.model.loop(header=(
        "_pdbx_refine_tls_group.id",
        "_pdbx_refine_tls_group.refine_tls_id",
        "_pdbx_refine_tls_group.pdbx_refine_id",
        "_pdbx_refine_tls_group.beg_auth_asym_id",
        "_pdbx_refine_tls_group.beg_auth_seq_id",
        "_pdbx_refine_tls_group.beg_label_asym_id",
        "_pdbx_refine_tls_group.beg_label_seq_id",
        "_pdbx_refine_tls_group.end_auth_asym_id",
        "_pdbx_refine_tls_group.end_auth_seq_id",
        "_pdbx_refine_tls_group.end_label_asym_id",
        "_pdbx_refine_tls_group.end_label_seq_id",
        "_pdbx_refine_tls_group.selection",
        "_pdbx_refine_tls_group.selection_details",
    ))
    for i_tls, (tlso, selection_string, isel) in enumerate(
        zip(self.tlsos, self.selection_strings, self.iselections)):
      tls_row = [i_tls+1, "?", scattering_type, "refined"]
      tls_row.extend(list(tlso.origin))
      tls_row.extend(list(tlso.t))
      tls_row.extend(list(tlso.l))
      tls_row.extend(list(tlso.s))
      tls_loop.add_row(tls_row)
      tls_group_loop.add_row({
          "_pdbx_refine_tls_group.id" : i_tls+1,
          "_pdbx_refine_tls_group.refine_tls_id" : i_tls+1,
          "_pdbx_refine_tls_group.pdbx_refine_id" : scattering_type,
          "_pdbx_refine_tls_group.beg_auth_asym_id" : hierarchy.get_auth_asym_id_iseq(isel[0]),
          "_pdbx_refine_tls_group.beg_auth_seq_id" : hierarchy.get_auth_seq_id_iseq(isel[0]),
          "_pdbx_refine_tls_group.beg_label_asym_id" : hierarchy.get_label_asym_id_iseq(isel[0]),
          "_pdbx_refine_tls_group.beg_label_seq_id" : hierarchy.get_label_seq_id_iseq(isel[0]),
          "_pdbx_refine_tls_group.end_auth_asym_id" : hierarchy.get_auth_asym_id_iseq(isel[-1]),
          "_pdbx_refine_tls_group.end_auth_seq_id" : hierarchy.get_auth_seq_id_iseq(isel[-1]),
          "_pdbx_refine_tls_group.end_label_asym_id" : hierarchy.get_label_asym_id_iseq(isel[-1]),
          "_pdbx_refine_tls_group.end_label_seq_id" : hierarchy.get_label_seq_id_iseq(isel[-1]),
          "_pdbx_refine_tls_group.selection" : '.', # Some ccp4 stuff
          "_pdbx_refine_tls_group.selection_details" : selection_string,
      })
    cif_block.add_loop(tls_loop)
    cif_block.add_loop(tls_group_loop)
    return cif_block

def remark_3_tls(tlsos, selection_strings, out = None):
  if(out is None): out = sys.stdout
  if (len(selection_strings) == 0):
    assert len(tlsos) == 1
    selection_strings = [None]
  else:
    assert len(tlsos) == len(selection_strings)
  print("REMARK   3  TLS DETAILS.", file=out)
  print("REMARK   3   NUMBER OF TLS GROUPS: %-6d"%len(tlsos), file=out)
  print("REMARK   3   ORIGIN: CENTER OF MASS", file=out)
  r3 = "REMARK   3   "
  counter = 0
  for tlso, selection_string in zip(tlsos, selection_strings):
    if(selection_string is None):
      selection_string = "all"
    counter += 1
    t = tlso.t
    l = tlso.l
    s = tlso.s
    o = tlso.origin
    print(r3+"TLS GROUP : %-6d"%(counter), file=out)
    lines = line_breaker(selection_string, width=45)
    for i_line, line in enumerate(lines):
      if(i_line == 0):
        print(r3+" SELECTION: %s"%line, file=out)
      else:
        print(r3+"          : %s"%line, file=out)
    print(r3+" ORIGIN FOR THE GROUP (A):%9.4f%9.4f%9.4f"%(o[0],o[1],o[2]), file=out)
    print(r3+" T TENSOR                                            ", file=out)
    print(r3+"   T11:%9.4f T22:%9.4f                       "%(t[0], t[1]), file=out)
    print(r3+"   T33:%9.4f T12:%9.4f                       "%(t[2], t[3]), file=out)
    print(r3+"   T13:%9.4f T23:%9.4f                       "%(t[4], t[5]), file=out)
    print(r3+" L TENSOR                                            ", file=out)
    print(r3+"   L11:%9.4f L22:%9.4f                       "%(l[0], l[1]), file=out)
    print(r3+"   L33:%9.4f L12:%9.4f                       "%(l[2], l[3]), file=out)
    print(r3+"   L13:%9.4f L23:%9.4f                       "%(l[4], l[5]), file=out)
    print(r3+" S TENSOR                                            ", file=out)
    print(r3+"   S11:%9.4f S12:%9.4f S13:%9.4f         "%(s[0],s[1],s[2]), file=out)
    print(r3+"   S21:%9.4f S22:%9.4f S23:%9.4f         "%(s[3],s[4],s[5]), file=out)
    print(r3+"   S31:%9.4f S32:%9.4f S33:%9.4f         "%(s[6],s[7],s[8]), file=out)

class show_tls(object):
  def __init__(self, tlsos, text="", out=None):
    if(out is None): out = sys.stdout
    counter = 0
    formatT ="|T11=%8.4f T22=%8.4f T33=%8.4f T12=%8.4f T13=%8.4f T23=%8.4f|"
    formatL ="|L11=%8.4f L22=%8.4f L33=%8.4f L12=%8.4f L13=%8.4f L23=%8.4f|"
    formatS1="|S11=%8.4f S22=%8.4f S33=%8.4f S12=%8.4f S13=%8.4f S21=%8.4f|\n"
    formatS2="|S23=%8.4f S31=%8.4f S32=%8.4f"+" "*39+"|"
    formatS = formatS1 + formatS2
    part1 = "|-"+text
    part2 = "-|"
    n = 79 - len(part1+part2)
    print(part1 + "-"*n + part2, file=out)
    for item in tlsos:
      counter += 1
      T = item.t
      L = item.l
      S = item.s
      origin = item.origin
      line = "|TLS group number %d: " % counter
      n = 79 - len(line+"|")
      print(line+" "*n+"|", file=out)
      print("|"+" "*15+"Origin (x,y,z) = %9.4f %9.4f %9.4f"%\
            (origin[0],origin[1],origin[2])+" "*16+"|", file=out)
      print(formatT % (T[0],T[1],T[2],T[3],T[4],T[5]), file=out)
      print(formatL % (L[0],L[1],L[2],L[3],L[4],L[5]), file=out)
      print(formatS % (S[0],S[4],S[8],S[1],S[2],S[3],S[5],S[6],S[7]), file=out)
    print("|" +"-"*77+"|", file=out)

class show_tls_one_group(object):
  def __init__(self, tlso, text="", out=None):
    if(out is None): out = sys.stdout
    counter = 0
    formatT ="|T11=%8.4f T22=%8.4f T33=%8.4f T12=%8.4f T13=%8.4f T23=%8.4f|"
    formatL ="|L11=%8.4f L22=%8.4f L33=%8.4f L12=%8.4f L13=%8.4f L23=%8.4f|"
    formatS1="|S11=%8.4f S22=%8.4f S33=%8.4f S12=%8.4f S13=%8.4f S21=%8.4f|\n"
    formatS2="|S23=%8.4f S31=%8.4f S32=%8.4f"+" "*39+"|"
    formatS = formatS1 + formatS2
    part1 = "|-"+text
    part2 = "-|"
    n = 79 - len(part1+part2)
    print(part1 + "-"*n + part2, file=out)
    counter += 1
    T = tlso.t
    L = tlso.l
    S = tlso.s
    origin = tlso.origin
    print("|"+" "*15+"Origin (x,y,z) = %9.4f %9.4f %9.4f"%\
          (origin[0],origin[1],origin[2])+" "*16+"|", file=out)
    print(formatT % (T[0],T[1],T[2],T[3],T[4],T[5]), file=out)
    print(formatL % (L[0],L[1],L[2],L[3],L[4],L[5]), file=out)
    print(formatS % (S[0],S[4],S[8],S[1],S[2],S[3],S[5],S[6],S[7]), file=out)
    print("|" +"-"*77+"|", file=out)

def u_cart_from_tls(sites_cart, selections, tlsos):
  global time_u_cart_from_tls
  t1 = time.time()
  uanisos = flex.sym_mat3_double(sites_cart.size(), [0,0,0,0,0,0])
  for selection, tlso in zip(selections, tlsos):
    u = uaniso_from_tls_one_group(tlso       = tlso,
                                  sites_cart = sites_cart.select(selection),
                                  zeroize_trace=True)
    uanisos.set_selected(selection, u)
  t2 = time.time()
  time_u_cart_from_tls += (t2 - t1)
  return uanisos

def tls_from_uanisos(xray_structure,
                     selections,
                     tlsos_initial,
                     number_of_macro_cycles       = 3000,
                     max_iterations               = 1000,
                     refine_T                     = True,
                     refine_L                     = True,
                     refine_S                     = True,
                     verbose                      = -1,
                     enforce_positive_definite_TL = True,
                     out                          = None):
  global time_tls_from_uanisos
  t1 = time.time()
  if(out is None): out = sys.stdout
  if(verbose > 0):
    show_tls(tlsos = tlsos_initial,
             text  = "TLS from ADP: start TLS values", out = out)
  u_cart=xray_structure.scatterers().extract_u_cart(xray_structure.unit_cell())
  T_min = []
  L_min = []
  S_min = []
  group_counter = 0
  for tlso_initial, selection in zip(tlsos_initial, selections):
    group_counter += 1
    T_initial = tlso_initial.t
    L_initial = tlso_initial.l
    S_initial = tlso_initial.s
    if(enforce_positive_definite_TL):
      T_initial = adptbx.eigenvalue_filtering(T_initial)
      L_initial = adptbx.eigenvalue_filtering(L_initial)
    stop_flag = 0
    target_stop = -1.0
    sites_cart_selected = xray_structure.sites_cart().select(selection)
    u_cart_selected = u_cart.select(selection)
    for i in range(1, number_of_macro_cycles+1):
      target_start = target_stop
      minimized = tls_from_uaniso_minimizer(uaniso    = u_cart_selected,
                                            T_initial = T_initial,
                                            L_initial = L_initial,
                                            S_initial = S_initial,
                                            refine_T  = refine_T,
                                            refine_L  = refine_L,
                                            refine_S  = refine_S,
                                            max_iterations = max_iterations,
                                            origin    = tlso_initial.origin,
                                            sites     = sites_cart_selected)
      if(refine_T):  T_initial = minimized.T_min
      else:          assert approx_equal(T_initial, minimized.T_min)
      if(refine_L):  L_initial = minimized.L_min
      else:          assert approx_equal(L_initial, minimized.L_min)
      if(refine_S):  S_initial = minimized.S_min
      else:          assert approx_equal(S_initial, minimized.S_min)
    if(verbose > 0):
      print("TLS group %d: minimized target = " %(group_counter),minimized.f, file=out)
    T_min_ = minimized.T_min
    L_min_ = minimized.L_min
    if(enforce_positive_definite_TL):
      T_min_ = adptbx.eigenvalue_filtering(T_min_)
      L_min_ = adptbx.eigenvalue_filtering(L_min_)
    T_min.append(T_min_)
    L_min.append(L_min_)
    S_min.append(minimized.S_min)
    if(enforce_positive_definite_TL):
      assert adptbx.is_positive_definite(T_min_, 1.e-6)
      assert adptbx.is_positive_definite(L_min_, 1.e-6)
  tlsos_result = generate_tlsos(selections     = selections,
                                xray_structure = xray_structure,
                                T              = T_min,
                                L              = L_min,
                                S              = S_min)
  if(verbose > 0):
    show_tls(tlsos = tlsos_result,
             text  = "TLS from ADP: final TLS values", out = out)
  t2 = time.time()
  time_tls_from_uanisos += (t2 - t1)
  return tlsos_result

class tls_from_uaniso_minimizer(object):
  def __init__(self,
               uaniso,
               T_initial,
               L_initial,
               S_initial,
               refine_T,
               refine_L,
               refine_S,
               origin,
               sites,
               max_iterations):
    adopt_init_args(self, locals())
    assert uaniso.size() == sites.size()
    self.dim_T = len(self.T_initial)
    self.dim_L = len(self.L_initial)
    self.dim_S = len(self.S_initial)
    assert self.dim_T == 6 and self.dim_S == 9
    self.T_min = self.T_initial
    self.L_min = self.L_initial
    self.S_min = self.S_initial
    self.x = self.pack(self.T_min, self.L_min, self.S_min)
    self.n = self.x.size()
    self.minimizer = lbfgs.run(
      target_evaluator = self,
      termination_params = lbfgs.termination_parameters(
        max_iterations = max_iterations,
        max_calls      = int(max_iterations*1.5)),
        exception_handling_params =
        lbfgs.exception_handling_parameters(
          ignore_line_search_failed_step_at_lower_bound = True,
          ignore_line_search_failed_step_at_upper_bound = True,
          ignore_line_search_failed_maxfev              = True)
        )
    self.compute_functional_and_gradients()
    del self.x

  def pack(self, T, L, S):
    v = []
    if (self.refine_T): v += list(T)
    if (self.refine_L): v += list(L)
    if (self.refine_S): v += list(S)
    return flex.double(tuple(v))

  def unpack_x(self):
    i = 0
    if (self.refine_T):
      self.T_min = tuple(self.x)[i:self.dim_T]
      i = self.dim_T
    if (self.refine_L):
      self.L_min = tuple(self.x)[i:i+self.dim_L]
      i += self.dim_L
    if (self.refine_S):
      self.S_min = tuple(self.x)[i:i+self.dim_S]

  def compute_functional_and_gradients(self):
    self.unpack_x()
    manager = tls_from_uaniso_target_and_grads(self.T_min,
                                               self.L_min,
                                               self.S_min,
                                               self.origin,
                                               self.sites,
                                               self.uaniso)
    self.f = manager.target()
    self.g = self.pack(manager.grad_T(), manager.grad_L(), manager.grad_S())
    return self.f, self.g

#######
class tls_from_uiso_minimizer(object):
  def __init__(self,
               uiso,
               T_initial,
               L_initial,
               S_initial,
               refine_T,
               refine_L,
               refine_S,
               origin,
               sites,
               max_iterations):
    adopt_init_args(self, locals())
    assert uiso.size() == sites.size()
    self.dim_T = len(self.T_initial)
    self.dim_L = len(self.L_initial)
    self.dim_S = len(self.S_initial)
    assert self.dim_T == 1 and self.dim_S == 3 and self.dim_L == 6
    self.T_min = self.T_initial
    self.L_min = self.L_initial
    self.S_min = self.S_initial
    self.x = self.pack(self.T_min, self.L_min, self.S_min)
    self.n = self.x.size()
    self.minimizer = lbfgs.run(
      target_evaluator = self,
      termination_params = lbfgs.termination_parameters(
        max_iterations = max_iterations,
        max_calls      = int(max_iterations*1.5)),
        exception_handling_params =
        lbfgs.exception_handling_parameters(
          ignore_line_search_failed_step_at_lower_bound = True,
          ignore_line_search_failed_step_at_upper_bound = True,
          ignore_line_search_failed_maxfev              = True)
        )
    self.compute_functional_and_gradients()
    del self.x

  def pack(self, T, L, S):
    v = []
    if (self.refine_T): v += list(flex.double([T]))
    if (self.refine_L): v += list(L)
    if (self.refine_S): v += list(S)
    return flex.double(tuple(v))

  def unpack_x(self):
    i = 0
    if (self.refine_T):
      self.T_min = tuple(self.x)[i:self.dim_T]
      i = self.dim_T
    if (self.refine_L):
      self.L_min = tuple(self.x)[i:i+self.dim_L]
      i += self.dim_L
    if (self.refine_S):
      self.S_min = tuple(self.x)[i:i+self.dim_S]

  def compute_functional_and_gradients(self):
    self.unpack_x()
    manager = tls_from_uiso_target_and_grads(self.T_min[0],
                                             self.L_min,
                                             self.S_min,
                                             self.origin,
                                             self.sites,
                                             self.uiso)
    self.f = manager.target()
    self.g = self.pack(manager.grad_T(), manager.grad_L(), manager.grad_S())
    return self.f, self.g
#######

class tls_xray_target_minimizer(object):
  def __init__(self,
               fmodel,
               tlsos_initial,
               refine_T,
               refine_L,
               refine_S,
               selections,
               selections_1d,
               max_iterations,
               run_finite_differences_test = False,
               correct_adp = True):
    adopt_init_args(self, locals())
    fmodel.xray_structure.scatterers().flags_set_grads(state=False)
    xray.set_scatterer_grad_flags(scatterers = fmodel.xray_structure.scatterers(),
                                  u_aniso     = True)
    if(self.run_finite_differences_test): self.correct_adp = False
    self.fmodel_copy = self.fmodel.deep_copy()
    self.target_functor = self.fmodel_copy.target_functor()
    self.run_finite_differences_test_counter = 0
    self.T_initial = []
    self.L_initial = []
    self.S_initial = []
    self.origins   = []
    for tlso_ in tlsos_initial:
      self.T_initial.append(tlso_.t)
      self.L_initial.append(tlso_.l)
      self.S_initial.append(tlso_.s)
      self.origins.append(tlso_.origin)
    self.counter = 0
    self.n_groups = len(self.T_initial)
    self.dim_T = len(self.T_initial[0])
    self.dim_L = len(self.L_initial[0])
    self.dim_S = len(self.S_initial[0])
    self.T_min = self.T_initial
    self.L_min = self.L_initial
    self.S_min = self.S_initial
    self.x = self.pack(self.T_min, self.L_min, self.S_min)
    self.minimizer = lbfgs.run(
      target_evaluator = self,
      core_params = lbfgs.core_parameters(maxfev = 10),
      termination_params = lbfgs.termination_parameters(
        min_iterations = max_iterations,
        max_calls = int(max_iterations*1.5)),
        exception_handling_params = lbfgs.exception_handling_parameters(
          ignore_line_search_failed_step_at_lower_bound = True,
          ignore_line_search_failed_step_at_upper_bound = True,
          ignore_line_search_failed_maxfev              = True))
    self.compute_functional_and_gradients()
    del self.x
    self.tlsos_result = generate_tlsos(
      selections     = self.selections,
      xray_structure = self.fmodel.xray_structure,
      T              = self.T_min,
      L              = self.L_min,
      S              = self.S_min)

  def pack(self, T, L, S):
    v = []
    for Ti,Li,Si in zip(T,L,S):
      if (self.refine_T): v += list(Ti)
      if (self.refine_L): v += list(Li)
      if (self.refine_S): v += list(Si)
    return flex.double(tuple(v))

  def unpack_x(self):
    i = 0
    T_min = []
    L_min = []
    S_min = []
    for j in range(self.n_groups):
      if (self.refine_T):
        self.T_min[j] = tuple(self.x)[i:i+self.dim_T]
        i += self.dim_T
      if (self.refine_L):
        self.L_min[j] = tuple(self.x)[i:i+self.dim_L]
        i += self.dim_L
      if (self.refine_S):
        self.S_min[j] = tuple(self.x)[i:i+self.dim_S]
        i += self.dim_S

  def compute_functional_and_gradients(self):
    self.counter += 1
    self.unpack_x()
    tlsos = generate_tlsos(selections     = self.selections,
                           xray_structure = self.fmodel_copy.xray_structure,
                           T              = self.T_min,
                           L              = self.L_min,
                           S              = self.S_min)
    update_xray_structure_with_tls(
      xray_structure = self.fmodel_copy.xray_structure,
      selections     = self.selections,
      selections_1d  = self.selections_1d,
      tlsos          = tlsos,
      correct_adp    = self.correct_adp)
    self.fmodel_copy.update_xray_structure(update_f_calc=True)
    t_r = self.target_functor(compute_gradients=True)
    self.f = t_r.target_work()
    grad_manager = tls_xray_grads(
      target_result=t_r,
      selections=self.selections,
      tlsos=tlsos)
    self.g = self.pack(grad_manager.grad_T,
                       grad_manager.grad_L,
                       grad_manager.grad_S)
    if(self.run_finite_differences_test and
       self.run_finite_differences_test_counter < 2):
      tolerance = 1.e-3
      self.run_finite_differences_test_counter += 1
      GT,GL,GS = finite_differences_grads_of_xray_target_wrt_tls(
        target_functor=self.target_functor,
        T=self.T_min,
        L=self.L_min,
        S=self.S_min,
        origins=self.origins,
        selections=self.selections,
        delta=0.00001)
      format   = "%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f"
      formats="%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f"
      for m1,m2 in zip(grad_manager.grad_T, GT):
        if(0):
          print("T1=" + format % (m1[0],m1[1],m1[2],m1[3],m1[4],m1[5]))
          print("T2=" + format % (m2[0],m2[1],m2[2],m2[3],m2[4],m2[5]))
        assert approx_equal(m1,m2,tolerance)
      for m1,m2 in zip(grad_manager.grad_L, GL):
        if(0):
          print("L1=" + format % (m1[0],m1[1],m1[2],m1[3],m1[4],m1[5]))
          print("L2=" + format % (m2[0],m2[1],m2[2],m2[3],m2[4],m2[5]))
        assert approx_equal(m1,m2,tolerance)
      for m1,m2 in zip(grad_manager.grad_S, GS):
        if(0):
          print("S1=" + formats %\
                (m1[0],m1[1],m1[2],m1[3],m1[4],m1[5],m1[6],m1[7],m1[8]))
          print("S2=" + formats %\
                (m2[0],m2[1],m2[2],m2[3],m2[4],m2[5],m2[6],m2[7],m2[8]))
        assert approx_equal(m1,m2,tolerance)
    return self.f, self.g

class tls_xray_grads(object):

  def __init__(self, target_result, selections, tlsos):
    self.grad_T = []
    self.grad_L = []
    self.grad_S = []
    d_target_d_uaniso = target_result.gradients_wrt_atomic_parameters(
      u_aniso=True)
    for sel, tlso in zip(selections, tlsos):
      d_target_d_tls_manager = d_target_d_tls(
        sites=target_result.manager.xray_structure.sites_cart().select(sel),
        origin            = tlso.origin,
        d_target_d_uaniso = d_target_d_uaniso.select(sel),
        scale_l_and_s     = True,# False will brake f.d. test
        use_trace_s_zero_constraint = True)
      self.grad_T.append(list(d_target_d_tls_manager.grad_T()))
      self.grad_L.append(list(d_target_d_tls_manager.grad_L()))
      self.grad_S.append(list(d_target_d_tls_manager.grad_S()))

def update_xray_structure_with_tls(xray_structure,
                                   selections,
                                   tlsos,
                                   selections_1d = None,
                                   correct_adp = True):
  global time_update_xray_structure_with_tls
  timer = user_plus_sys_time()
  u_cart_from_tls_ = u_cart_from_tls(sites_cart = xray_structure.sites_cart(),
                                     selections = selections,
                                     tlsos      = tlsos)
  xray_structure.set_u_cart(u_cart=u_cart_from_tls_, selection = selections_1d)
  if(correct_adp): xray_structure.tidy_us(u_min = 1.e-6)
  time_update_xray_structure_with_tls += timer.elapsed()

def split_u(xray_structure, tls_selections, offset):
  global time_split_u
  timer = user_plus_sys_time()
  uc = xray_structure.unit_cell()
  u_iso = xray_structure.scatterers().extract_u_iso()
  u_eq_1  = xray_structure.extract_u_iso_or_u_equiv()
  for tls_selection in tls_selections:
    u_iso_sel = u_iso.select(tls_selection)
    u_iso_min = flex.min(u_iso_sel)
    if(offset):
      offset_ = adptbx.b_as_u(5.0)
    else: offset_ = 0.0
    if u_iso_min >= offset_:
      u_iso_min = u_iso_min - offset_
    t = adptbx.u_iso_as_u_star(uc, u_iso_min)
    for i_seq in tls_selection:
      sc = xray_structure.scatterers()[i_seq]
      assert sc.u_iso == u_iso[i_seq]
      u_iso_new = sc.u_iso - u_iso_min
      assert u_iso_new >= 0.0
      sc.u_iso = u_iso_new
      assert sc.flags.use_u_aniso()
      assert sc.flags.use_u_iso()
      if(sc.u_star == (-1.0,-1.0,-1.0,-1.0,-1.0,-1.0)):
        sc.u_star = t
      else:
        x = flex.double(sc.u_star)
        y = flex.double(t)
        z = list(x + y)
        sc.u_star = z
  u_iso = xray_structure.scatterers().extract_u_iso().select(
    xray_structure.use_u_iso())
  assert (u_iso < 0.0).count(True) == 0
  u_eq_2  = xray_structure.extract_u_iso_or_u_equiv()
  assert approx_equal(u_eq_1, u_eq_2)
  time_split_u += timer.elapsed()

def tls_from_u_cart(xray_structure,
                    tlsos_initial,
                    tls_selections,
                    number_of_macro_cycles = 100,
                    max_iterations         = 100):
  global time_tls_from_u_cart
  timer = user_plus_sys_time()
  uc = xray_structure.unit_cell()
  xray_structure.tidy_us(u_min = 1.e-6)
  ueq = xray_structure.extract_u_iso_or_u_equiv()
  assert (ueq < 0.0).count(True) == 0
  u_cart = xray_structure.scatterers().extract_u_cart(uc)
  for tls_selection in tls_selections:
    u_cart_selected = u_cart.select(tls_selection)
    assert adptbx.is_positive_definite(u_cart_selected,1.e-6).count(False)==0
  xray_structure.tidy_us(u_min = 1.e-6)
  t = []
  l = []
  s = []
  lim_l = 10.0
  for tls_selection, tlso in zip(tls_selections, tlsos_initial):
    t.append( tlso.t )
    #if(abs(tlso.t[0]) < eps and abs(tlso.t[1]) < eps and abs(tlso.t[2]) < eps and
    #   abs(tlso.t[3]) < eps and abs(tlso.t[4]) < eps and abs(tlso.t[5]) < eps):
    #   t.append( t_from_u_cart(u_cart.select(tls_selection), 1.e-6) )
    #else:
    #   t.append( tlso.t )
    #l.append( [0,0,0,0,0,0] )
    if abs(tlso.l[0])>lim_l: l1 = tlso.l[0]/5.
    else:                    l1 = tlso.l[0]
    if abs(tlso.l[1])>lim_l: l2 = tlso.l[1]/5.
    else:                    l2 = tlso.l[1]
    if abs(tlso.l[2])>lim_l: l3 = tlso.l[2]/5.
    else:                    l3 = tlso.l[2]
    if abs(tlso.l[3])>lim_l: l4 = tlso.l[3]/5.
    else:                    l4 = tlso.l[3]
    if abs(tlso.l[4])>lim_l: l5 = tlso.l[4]/5.
    else:                    l5 = tlso.l[4]
    if abs(tlso.l[5])>lim_l: l6 = tlso.l[5]/5.
    else:                    l6 = tlso.l[5]
    l.append( [l1,l2,l3,l4,l5,l6] )
    s.append( tlso.s )
  tlsos = generate_tlsos(selections     = tls_selections,
                         xray_structure = xray_structure,
                         T = t, L = l, S = s)
  #for rt,rl,rs in [[0,1,1],[1,0,0]]*3:
  for rt,rl,rs in [[0,1,1],[1,0,0],[1,1,1]]:
    tlsos_ = tls_from_uanisos(xray_structure        = xray_structure,
                              selections             = tls_selections,
                              tlsos_initial          = tlsos,
                              number_of_macro_cycles = number_of_macro_cycles,
                              max_iterations         = max_iterations,
                              refine_T               = rt,
                              refine_L               = rl,
                              refine_S               = rs,
                              verbose                = -1,
                              out                    = None)
    tlsos = tlsos_
  t = []
  l = []
  s = []
  for tlso in tlsos:
    t.append( tlso.t )
    if abs(tlso.l[0])>lim_l: l1 = lim_l/5.
    else:                    l1 = tlso.l[0]
    if abs(tlso.l[1])>lim_l: l2 = lim_l/5.
    else:                    l2 = tlso.l[1]
    if abs(tlso.l[2])>lim_l: l3 = lim_l/5.
    else:                    l3 = tlso.l[2]
    if abs(tlso.l[3])>lim_l: l4 = lim_l/5.
    else:                    l4 = tlso.l[3]
    if abs(tlso.l[4])>lim_l: l5 = lim_l/5.
    else:                    l5 = tlso.l[4]
    if abs(tlso.l[5])>lim_l: l6 = lim_l/5.
    else:                    l6 = tlso.l[5]
    l.append( [l1,l2,l3,l4,l5,l6] )
    s.append( tlso.s )
  tlsos = generate_tlsos(selections     = tls_selections,
                         xray_structure = xray_structure,
                         T = t, L = l, S = s)
  time_tls_from_u_cart += timer.elapsed()
  return tlsos

class tls_refinement(object):
  def __init__(self,
               fmodel,
               model,
               selections,
               selections_1d,
               refine_T,
               refine_L,
               refine_S,
               number_of_macro_cycles,
               max_number_of_iterations,
               start_tls_value = None,
               run_finite_differences_test = False,
               eps = 1.e-6,
               out = None,
               macro_cycle = None,
               verbose = True):
    global time_tls_total
    timer = user_plus_sys_time()
    if(out is None): out = sys.stdout
    prefix = "TLS refinement:"
    fmodel.info().show_targets(text = prefix+" start model", out = out)
    fmodel.xray_structure.show_u_statistics(text = prefix+" start model",
                                            out  = out)
    xrs = fmodel.xray_structure
    xrs.tidy_us(u_min = 1.e-6)
    if(start_tls_value is not None):
      try:
        crash_or_not = abs(start_tls_value + 0)
        tlsos = generate_tlsos(value          = start_tls_value,
                               selections     = selections,
                               xray_structure = xrs)
      except Exception:
        tlsos = start_tls_value
    else:
      tlsos = tls_from_u_cart(xray_structure = xrs,
                              tlsos_initial  = model.tls_groups.tlsos,
                              tls_selections = selections,
                              number_of_macro_cycles = 100,
                              max_iterations         = 100)
    if (verbose):
      show_tls(tlsos = tlsos, text = prefix+" start parameters",out = out)
    for macro_cycle in range(1, number_of_macro_cycles+1):
      print(file=out)
      prefix = "TLS refinement: after macrocycle "+str(macro_cycle)
      minimized = tls_xray_target_minimizer(
        fmodel                      = fmodel,
        tlsos_initial               = tlsos,
        refine_T                    = refine_T,
        refine_L                    = refine_L,
        refine_S                    = refine_S,
        selections                  = selections,
        selections_1d               = selections_1d,
        max_iterations              = max_number_of_iterations,
        run_finite_differences_test = run_finite_differences_test)
      xrs = minimized.fmodel_copy.xray_structure
      xrs.show_u_statistics(text = prefix, out  = out)
      if(verbose):
        show_tls(tlsos = minimized.tlsos_result, text = prefix, out = out)
      fmodel.update_xray_structure(xray_structure = xrs,
                                   update_f_calc  = True)
      fmodel.info().show_targets(text = prefix, out = out)
      if(xrs.is_positive_definite_u().count(False) > 0):
        xrs.tidy_us(u_min = 1.e-6)
        xrs.show_u_statistics(
          text = prefix+": after making positive definite",
          out  = out)
        fmodel.update_xray_structure(xray_structure = xrs,
                                     update_f_calc  = True)
        fmodel.info().show_targets(text=prefix+": after making positive definite",
                                   out = out)
        tlsos = make_tlso_compatible_with_u_positive_definite(
          tlsos                            = minimized.tlsos_result,
          xray_structure                   = xrs.deep_copy_scatterers(),
          selections                       = selections,
          max_iterations                   = 10,
          number_of_u_nonpositive_definite = 0,
          eps                              = eps,
          refine_T                         = refine_T,
          refine_L                         = refine_L,
          refine_S                         = refine_S,
          out                              = out,
          number_of_macro_cycles_for_tls_from_uanisos = 10)
      else: tlsos = minimized.tlsos_result
    if (verbose):
      show_tls(tlsos = tlsos,
               text = "TLS refinement: final values", out = out)
    self.tlsos = tlsos
    model.tls_groups.tlsos = tlsos
    self.fmodel = fmodel
    time_tls_total += timer.elapsed()


def make_tlso_compatible_with_u_positive_definite(
  xray_structure,
  selections,
  max_iterations,
  number_of_u_nonpositive_definite,
  eps,
  number_of_macro_cycles_for_tls_from_uanisos,
  tlsos=None,
  refine_T=True,
  refine_L=True,
  refine_S=True,
  out=None):
  global time_make_tlso_compatible_with_u_positive_definite
  t1 = time.time()
  if(out is None): out = sys.stdout
  for i in range(1, max_iterations+1):
    update_xray_structure_with_tls(
      xray_structure = xray_structure,
      selections     = selections,
      tlsos          = tlsos)
    ipd_1 = xray_structure.is_positive_definite_u()
    if(i == 1 or i == max_iterations):
      xray_structure.show_u_statistics(out = out)
    xray_structure.tidy_us(u_min = 1.e-6)
    tlsos = tls_from_uanisos(
      xray_structure         = xray_structure,
      selections             = selections,
      tlsos_initial          = tlsos,
      refine_T               = refine_T,
      refine_L               = refine_L,
      refine_S               = refine_S,
      max_iterations         = 100,
      number_of_macro_cycles = number_of_macro_cycles_for_tls_from_uanisos)
    if(i == max_iterations): xray_structure.show_u_statistics(out = out)
    if(ipd_1.count(False) == number_of_u_nonpositive_definite):
      break
  assert xray_structure.is_positive_definite_u().count(False) == 0
  t2 = time.time()
  time_make_tlso_compatible_with_u_positive_definite += (t2 - t1)
  return tlsos

def generate_tlsos(selections,
                   xray_structure,
                   value=None, T=None, L=None, S=None):
  global time_generate_tlsos
  t1 = time.time()
  if(value is None): assert [T,L,S].count(None) == 0
  else:              assert [T,L,S].count(None) == 3
  if(value is not None):
    T       = []
    L       = []
    S       = []
    v       = value
    vtl     = [v,v,v,v,v,v]
    vs      = [v,v,v,v,v,v,v,v,v]
    for selection in selections:
      T.append(vtl)
      L.append(vtl)
      S.append(vs)
  origins = []
  for selection in selections:
    xrs = xray_structure.select(selection)
    origins.append(xrs.center_of_mass())
  tlsos = []
  for T_,L_,S_,origin_ in zip(T, L, S, origins):
    tlsos.append(tlso(t = T_, l = L_, s = S_, origin = origin_))
  t2 = time.time()
  time_generate_tlsos += (t2 - t1)
  return tlsos

def finite_differences_grads_of_xray_target_wrt_tls(target_functor,
                                                    T,
                                                    L,
                                                    S,
                                                    origins,
                                                    selections,
                                                    delta=0.00001):
  fmodel = target_functor.manager
  derivative_T = []
  for j in range(len(T)):
    dT = []
    for i in range(6):
      target_values = []
      for d_sign in (-1, 1):
        T_ = []
        for item in T:
          T_.append(list(item))
        d = d_sign*delta
        T_[j][i] += d
        #
        tlsos = []
        for Ti,Li,Si,origini in zip(T_, L, S, origins):
          tlsos.append(tlso(t      = Ti,
                            l      = Li,
                            s      = Si,
                            origin = origini))
        update_xray_structure_with_tls(
          xray_structure = fmodel.xray_structure,
          selections     = selections,
          tlsos          = tlsos,
          correct_adp    = False)
        fmodel.update_xray_structure(update_f_calc=True)
        t_w = target_functor(compute_gradients=False).target_work()
        #
        target_values.append(t_w)
      derivative = (target_values[1] - target_values[0]) / (2 * delta)
      dT.append(derivative)
    derivative_T.append(dT)

  derivative_L = []
  for j in range(len(L)):
    dL = []
    for i in range(6):
      target_values = []
      for d_sign in (-1, 1):
        L_ = []
        for item in L:
          L_.append(list(item))
        d = d_sign*delta
        L_[j][i] += d
        #
        tlsos = []
        for Ti,Li,Si,origini in zip(T, L_, S, origins):
          tlsos.append(tlso(t      = Ti,
                            l      = Li,
                            s      = Si,
                            origin = origini))
        update_xray_structure_with_tls(
          xray_structure = fmodel.xray_structure,
          selections     = selections,
          tlsos          = tlsos,
          correct_adp    = False)
        fmodel.update_xray_structure(update_f_calc=True)
        t_w = target_functor(compute_gradients=False).target_work()
        #
        target_values.append(t_w)
      derivative = (target_values[1] - target_values[0]) / (2 * delta)
      dL.append(derivative)
    derivative_L.append(dL)

  derivative_S = []
  for j in range(len(L)):
    dS = []
    for i in range(9):
      target_values = []
      for d_sign in (-1, 1):
        S_ = []
        for item in S:
          S_.append(list(item))
        d = d_sign*delta
        S_[j][i] += d
        #
        tlsos = []
        for Ti,Li,Si,origini in zip(T, L, S_, origins):
          tlsos.append(tlso(t      = Ti,
                            l      = Li,
                            s      = Si,
                            origin = origini))
        update_xray_structure_with_tls(
          xray_structure = fmodel.xray_structure,
          selections     = selections,
          tlsos          = tlsos,
          correct_adp    = False)
        fmodel.update_xray_structure(update_f_calc=True)
        t_w = target_functor(compute_gradients=False).target_work()
        #
        target_values.append(t_w)
      derivative = (target_values[1] - target_values[0]) / (2 * delta)
      dS.append(derivative)
    derivative_S.append(dS)
  return derivative_T,derivative_L,derivative_S

def water_in_tls_selections(
    tls_selections,
    pdb_hierarchy):
  cache = pdb_hierarchy.atom_selection_cache()
  water_sel = cache.selection("resname HOH")
  tls_sel = flex.bool(water_sel.size(), False)
  for sele_str in tls_selections :
    group_sel = cache.selection(sele_str)
    tls_sel |= group_sel
  if ((tls_sel & water_sel).count(True) > 0):
    return True
  return False


def u_cart_from_ensemble(models):
  xyz_all = []
  for m in models:
    xyz_all.append(m.atoms().extract_xyz())
  n_atoms = xyz_all[0].size()
  xyz_atoms_all = []
  for i in range(n_atoms):
    xyz_atoms = flex.vec3_double()
    for xyzs in xyz_all:
      xyz_atoms.append(xyzs[i])
    xyz_atoms_all.append(xyz_atoms)
  result = flex.sym_mat3_double()
  for i in range(n_atoms):
    result.append(u_cart_from_xyz(sites_cart=xyz_atoms_all[i]))
  return result

# TEST utils

def get_u_cart(o_tfm, origin, sites_cart):
  import math
  scale = 180./math.pi # trick to work-around c++ implementation
  L = [v*scale**2 for v in o_tfm.L_M.as_sym_mat3()]
  S = [v*scale for v in o_tfm.S_M.as_mat3()]
  tlso_ = tlso(
    t      = o_tfm.T_M.as_sym_mat3(),
    l      = L,
    s      = S,
    origin = origin)
  return uaniso_from_tls_one_group(
    tlso          = tlso_,
    sites_cart    = sites_cart,
    zeroize_trace = False)

class u_tls_vs_u_ens(object):
  def __init__(self,
               pdb_str,
               dx=0,dy=0,dz=0,
               sx=0,sy=0,sz=0,
               lx=[1,0,0],ly=[0,1,0],lz=[0,0,1],
               tx=0,ty=0,tz=0,
               vx=[1,0,0],vy=[0,1,0],vz=[0,0,1],
               w_M_lx=[0,0,0], w_M_ly=[0,0,0], w_M_lz=[0,0,0],
               origin=None,
               n_models=10000,
               assert_similarity=True,
               show=False,
               log = sys.stdout,
               write_pdb_files=False,
               smear_eps = 0):
    from mmtbx.tls import analysis, tls_as_xyz
    from scitbx import matrix
    from libtbx.utils import null_out
    if(show):
      print("INPUTS:","-"*73, file=log)
      print("dx    :", dx, file=log)
      print("dy    :", dy, file=log)
      print("dz    :", dz, file=log)
      print("sx    :", sx, file=log)
      print("sy    :", sy, file=log)
      print("sz    :", sz, file=log)
      print("lx    :", [i for i in lx], file=log)
      print("ly    :", [i for i in ly], file=log)
      print("lz    :", [i for i in lz], file=log)
      print("tx    :", tx, file=log)
      print("ty    :", ty, file=log)
      print("tz    :", tz, file=log)
      print("vx    :", [i for i in vx], file=log)
      print("vy    :", [i for i in vy], file=log)
      print("vz    :", [i for i in vz], file=log)
      print("w_M_lx:", [i for i in w_M_lx], file=log)
      print("w_M_ly:", [i for i in w_M_ly], file=log)
      print("w_M_lz:", [i for i in w_M_lz], file=log)
      print("origin:", origin, file=log)
      print("-"*79, file=log)
    #
    pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
    ph = pdb_inp.construct_hierarchy()
    xrs = ph.extract_xray_structure(
      crystal_symmetry = pdb_inp.crystal_symmetry())
    sites_cart = xrs.sites_cart()
    ph.atoms().set_xyz(sites_cart)
    if(origin is None):
      origin = sites_cart.mean()
    #
    o_tfm = analysis.tls_from_motions(
      dx=dx,dy=dy,dz=dz,
      l_x=matrix.col(lx),l_y=matrix.col(ly),l_z=matrix.col(lz),
      sx=sx,sy=sy,sz=sz,
      tx=tx,ty=ty,tz=tz,
      v_x=matrix.col(vx),v_y=matrix.col(vy),v_z=matrix.col(vz),
      w_M_lx=matrix.col(w_M_lx),
      w_M_ly=matrix.col(w_M_ly),
      w_M_lz=matrix.col(w_M_lz))
    #
    self.u_cart_tls = get_u_cart(
      o_tfm=o_tfm, origin=origin, sites_cart=sites_cart)
    tlso_ = tlso(
      t      = o_tfm.T_M.as_sym_mat3(),
      l      = o_tfm.L_M.as_sym_mat3(),
      s      = o_tfm.S_M.as_mat3(),
      origin = origin)
    if(assert_similarity):
      T = matrix.sym(sym_mat3=tlso_.t)
      L = matrix.sym(sym_mat3=tlso_.l)
      S = matrix.sqr(tlso_.s)
      o_tfm = analysis.run(T=T, L=L, S=S, log=null_out()).self_check()
    #
    r = tls_as_xyz.ensemble_generator(
      tls_from_motions_object = o_tfm,
      pdb_hierarchy        = ph,
      xray_structure       = xrs,
      n_models             = n_models,
      origin               = origin,
      use_states           = write_pdb_files,
      log                  = null_out())
    if(write_pdb_files):
      r.write_pdb_file(file_name="ensemble_%s.pdb"%str(n_models))
    #
    xyz_all = r.sites_cart_ens
    n_atoms = xyz_all[0].size()
    ###
    xyz_atoms_all = all_vs_all(xyz_all = xyz_all)
    ###
    self.u_cart_ens = flex.sym_mat3_double()
    for i in range(n_atoms):
      self.u_cart_ens.append(u_cart_from_xyz(sites_cart=xyz_atoms_all[i]))
    u1 = self.u_cart_tls.as_double()
    u2 = self.u_cart_ens.as_double()
    #self.r = flex.sum(flex.abs(u1-u2))/\
    #         flex.sum(flex.abs(flex.abs(u1)+flex.abs(u2)))*2
    # LS
    #diff = u1-u2
    #self.rLS = math.sqrt(flex.sum(diff*diff)/(9.*diff.size()))
    #
    # Merritt / Murshudov
    e = smear_eps
    eps = matrix.sqr(
      [e, 0, 0,
       0, e, 0,
       0, 0, e])
    I = matrix.sqr(
      [2, 0, 0,
       0, 2, 0,
       0, 0, 2])

    def add_const(u):
      es = eigensystem.real_symmetric(u)
      vecs = es.vectors()
      l_z = matrix.col((vecs[0], vecs[1], vecs[2]))
      l_y = matrix.col((vecs[3], vecs[4], vecs[5]))
      l_x = matrix.col((vecs[6], vecs[7], vecs[8]))
      #l_x = l_y.cross(l_z)
      u = matrix.sym(sym_mat3=u)
      R = matrix.sqr(
        [l_x[0], l_y[0], l_z[0],
         l_x[1], l_y[1], l_z[1],
         l_x[2], l_y[2], l_z[2]])
      uD = R.transpose()*u*R
      result = R*(uD+eps)*R.transpose()
      tmp = R*uD*R.transpose()
      for i in range(6):
        assert approx_equal(tmp[i], u[i])
      return R*(uD+eps)*R.transpose()

    self.KL=0
    self.CC=0
    n1,n2,d1,d2=0,0,0,0
    for u1, u2 in zip(self.u_cart_tls, self.u_cart_ens):
      for i in range(6):
        n1 += abs(u1[i]-u2[i])
        d1 += (abs(u1[i])+abs(u2[i]))
      u1 = add_const(u=u1)
      u2 = add_const(u=u2)
      for i in range(6):
        n2 += abs(u1[i]-u2[i])
        d2 += (abs(u1[i])+abs(u2[i]))
      iu1 = u1.inverse()
      iu2 = u2.inverse()
      self.KL += (u1*iu2+u2*iu1-I).trace()
      diu1 = iu1.determinant()
      diu2 = iu2.determinant()
      den = (iu1+iu2).determinant()
      self.CC += (diu1*diu2)**0.25/(den/8)**0.5
    self.KL = self.KL/self.u_cart_ens.size()
    self.CC = self.CC/self.u_cart_ens.size()
    self.R1 = n1/d1*2
    self.R2 = n2/d2*2
    self.r = self.R1
    #
    ###
    for i in range(n_atoms):
      ut=["%8.5f"%u for u in self.u_cart_tls[i]]
      ue=["%8.5f"%u for u in self.u_cart_ens[i]]
      if(assert_similarity):
        for j in range(6):
          assert approx_equal(abs(float(ut[j])), abs(float(ue[j])), 1.e-3)
    #
    if(write_pdb_files):
      ph.atoms().set_uij(self.u_cart_tls)
      ph.write_pdb_file(
        file_name = "u_from_tls.pdb",
        crystal_symmetry = xrs.crystal_symmetry())
      ph.atoms().set_uij(self.u_cart_ens)
      ph.write_pdb_file(
        file_name = "u_from_ens.pdb",
        crystal_symmetry = xrs.crystal_symmetry())
