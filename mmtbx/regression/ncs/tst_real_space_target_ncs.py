from __future__ import absolute_import, division, print_function
import time
import mmtbx.refinement.minimization_ncs_constraints as mnc
from mmtbx.regression.ncs import tst_minimization_ncs_constraints_real_space as t0
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
import iotbx.ncs
import iotbx.pdb
import mmtbx.utils
import mmtbx.restraints

def get_inputs(d_min):
  pdb_inp_answer = iotbx.pdb.input(source_info=None, lines=t0.pdb_str_answer)
  ph_answer = pdb_inp_answer.construct_hierarchy()
  ph_answer.atoms().reset_i_seq()
  xrs_answer = pdb_inp_answer.xray_structure_simple()
  xrs_poor = iotbx.pdb.input(
    source_info=None, lines=t0.pdb_str_poor).xray_structure_simple()
  ppf = mmtbx.utils.process_pdb_file_srv().process_pdb_files(
    raw_records=t0.pdb_str_poor.splitlines())[0]
  restraints_manager = mmtbx.restraints.manager(
    geometry      = ppf.geometry_restraints_manager(show_energies = False),
    normalization = True)
  fft_map = xrs_answer.structure_factors(d_min=d_min).f_calc().fft_map(
    resolution_factor = 0.25)
  fft_map.apply_sigma_scaling()
  ncs_obj = iotbx.ncs.input(hierarchy=ph_answer)
  nrgl = ncs_obj.get_ncs_restraints_group_list()
  refine_selection = flex.size_t(range(xrs_poor.scatterers().size()))
  return xrs_answer, xrs_poor, restraints_manager, \
    fft_map.real_map_unpadded(), nrgl, refine_selection

def get_tfg(map_data, xrs, nrgl, refine_selection, restraints_manager, d_min,
            **kw):
  return mnc.target_function_and_grads_real_space(
    map_data                   = map_data,
    xray_structure             = xrs.deep_copy_scatterers(),
    ncs_restraints_group_list  = nrgl,
    refine_selection           = refine_selection,
    real_space_gradients_delta = d_min/4,
    restraints_manager         = restraints_manager,
    data_weight                = 1,
    refine_sites               = True,
    **kw)

def exercise_gradients(d_min=1.0):
  """
  Map term of the NCS real-space target: for each gradients_method the
  gradients agree with finite differences of the target of the same method.
  """
  xrs_answer, xrs_poor, rm, map_data, nrgl, rsel = get_inputs(d_min=d_min)
  for method, h, eps in [("fd",        d_min/4, 1.e-6),
                         ("linear",    1.e-5,   1.e-4),
                         ("tricubic",  1.e-5,   1.e-4)]:
    tfg = get_tfg(map_data, xrs_poor, nrgl, rsel, rm, d_min,
      gradients_method = method)
    t, g = tfg.data_target_and_grads(compute_gradients=True, x=None)
    sites = xrs_poor.sites_cart()
    for i_site in [0, 7, 13]:
      for i_axis in [0,1,2]:
        ts = []
        for sign in [1,-1]:
          s = list(sites[i_site])
          s[i_axis] += sign*h
          sites_ = sites.deep_copy()
          sites_[i_site] = s
          tfg.xray_structure.set_sites_cart(sites_)
          ts.append(tfg.data_target_and_grads(
            compute_gradients=False, x=None)[0])
        g_fd = (ts[0]-ts[1])/(2*h)
        assert approx_equal(g[i_site][i_axis], g_fd, eps), [
          method, i_site, i_axis, g[i_site][i_axis], g_fd]

def exercise_selection_real_space(d_min=1.0):
  """
  Atoms excluded with selection_real_space do not contribute to the map term.
  """
  xrs_answer, xrs_poor, rm, map_data, nrgl, rsel = get_inputs(d_min=d_min)
  n = xrs_poor.scatterers().size()
  sel = flex.bool(n, True)
  excluded = [0, 1, 2, n-1]
  for i in excluded: sel[i] = False
  for method in ["fd", "linear", "tricubic"]:
    t_all, g_all = get_tfg(map_data, xrs_poor, nrgl, rsel, rm, d_min,
      gradients_method = method).data_target_and_grads(True, None)
    t_sel, g_sel = get_tfg(map_data, xrs_poor, nrgl, rsel, rm, d_min,
      gradients_method     = method,
      selection_real_space = sel).data_target_and_grads(True, None)
    assert abs(t_all-t_sel) > 1.e-3
    for i in range(n):
      if(i in excluded):
        assert approx_equal(g_sel[i], [0,0,0], 1.e-12)
        assert flex.max(flex.abs(flex.double(g_all[i]))) > 1.e-3
      else:
        assert approx_equal(g_sel[i], g_all[i], 1.e-9)

def exercise_minimization(d_min=1.0):
  """
  NCS constrained minimization gets to the answer with each gradients_method.
  """
  for method in ["fd", "linear", "tricubic"]:
    xrs_answer, xrs, rm, map_data, nrgl, rsel = get_inputs(d_min=d_min)
    d_start = flex.mean(xrs_answer.distances(xrs))
    for i in range(5):
      tfg = get_tfg(map_data, xrs, nrgl, rsel, rm, d_min,
        gradients_method = method)
      mnc.lbfgs(
        target_and_grads_object      = tfg,
        xray_structure               = xrs,
        ncs_restraints_group_list    = nrgl,
        refine_selection             = rsel,
        finite_grad_differences_test = False,
        max_iterations               = 60,
        refine_sites                 = True)
      xrs = tfg.xray_structure
    d_final = flex.mean(xrs_answer.distances(xrs))
    print(method, "start: %6.3f final: %6.3f"%(d_start, d_final))
    assert d_start > 0.5
    assert d_final < 0.05, [method, d_final]

if (__name__ == "__main__"):
  t_start=time.time()
  exercise_gradients()
  exercise_selection_real_space()
  exercise_minimization()
  print("Time: %6.4f"%(time.time()-t_start))
  print("OK")
