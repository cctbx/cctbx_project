from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from libtbx import group_args
from libtbx.test_utils import approx_equal
import mmtbx.refinement.llgi_curvature as llgi_curvature
import mmtbx.refinement.llgi_hybrid as llgi_hybrid

""" mmtbx.refinement.llgi_curvature: the Fisher information of the llgi
target per reflection, with and without the tNCS frame, and the inverse
diagonal Hessian it gives for coordinate refinement. """

def exercise():
  from mmtbx.regression.llgi_test_utils import build_llgi_fmodel
  fmodel = build_llgi_fmodel(40, 2.2, seed=3)    # TEPS = 1
  n = fmodel.f_obs().size()
  free = fmodel.r_free_flags().data()
  info = llgi_curvature.fisher_information(fmodel)
  assert info.size() == n
  assert flex.min(info) >= 0
  assert info.select(free).all_eq(0)
  assert (info.select(~free) > 0).count(True) > 0.9 * (~free).count(True)
  # tNCS arrays with u = t = 1 give the same; a reflection with t = 0 (an
  # invalid frame) gets no information, without an exception
  ld = fmodel.llgi_data()
  one = fmodel.f_obs().array(data=flex.double(n, 1.))
  fmodel.set_llgi_data(llgi_hybrid.replace_llgi_data(ld, tncs_u=one,
    tncs_t=one))
  assert approx_equal(llgi_curvature.fisher_information(fmodel), info,
    eps=1e-12)
  t = flex.double(n, 1.)
  i_zero = (~free).iselection()[0]
  t[i_zero] = 0.
  fmodel.set_llgi_data(llgi_hybrid.replace_llgi_data(ld, tncs_u=one,
    tncs_t=fmodel.f_obs().array(data=t)))
  info_t = llgi_curvature.fisher_information(fmodel)
  assert info_t[i_zero] == 0
  assert approx_equal(info_t.select(t > 0), info.select(t > 0), eps=1e-12)
  fmodel.set_llgi_data(ld)
  # inverse diagonal for the coordinates of every scatterer
  xrs = fmodel.xray_structure
  xrs.scatterers().flags_set_grads(state=False)
  xrs.scatterers().flags_set_grad_site(
    iselection=flex.size_t(range(xrs.scatterers().size())))
  diag, h_data, h_restr = llgi_curvature.xyz_inverse_diagonal(fmodel=fmodel,
    model=None, weights=group_args(wx=1., wx_scale=1.), xray_structure=xrs,
    restraints_weight=0.)
  assert diag.size() == 3 * xrs.scatterers().size()
  assert flex.min(diag) > 0 and flex.min(h_data) > 0 and h_restr == 0

def exercise_symmetry_sum():
  """ the data curvature against sum_h I_h |dF_h/dx|^2 / N_work with dF_h/dx
  by finite differences of direct summation over all symmetry operators,
  centring included """
  from mmtbx.regression.llgi_test_utils import build_fmodel, \
    synthetic_llgi_data
  for sg in ("P 1 21 1", "C 1 2 1", "F 2 2 2"):
    fmodel = build_fmodel(n_atoms=20, d_min=2.0, seed=1, space_group=sg,
      update_scales=True)
    fmodel.set_llgi_data(synthetic_llgi_data(fmodel=fmodel, seed=101,
      feff_scale=1.0))
    fmodel.set_target_name("llgi")
    fmodel.update_llgi_sigmaa_scatfrac()
    xrs = fmodel.xray_structure
    info = llgi_curvature.fisher_information(fmodel)
    n_work = (~fmodel.r_free_flags().data()).count(True)
    h = llgi_curvature.data_curvatures(fmodel, xrs)
    ratios = []
    for j in range(4):
      x1 = xrs.select(flex.size_t([j]))
      s0 = x1.sites_cart()[0]
      total = 0.
      for c in range(3):
        f = []
        for step in (1.e-4, -1.e-4):
          s = list(s0); s[c] += step
          x1.set_sites_cart(flex.vec3_double([s]))
          f.append(fmodel.f_obs().structure_factors_from_scatterers(
            xray_structure=x1, algorithm="direct").f_calc().data())
        total += flex.sum(info * flex.norm((f[0] - f[1]) / 2.e-4)) / 3
      ratios.append(total / n_work / h[j])
    assert approx_equal(sum(ratios) / len(ratios), 1, eps=0.05), (sg, ratios)

if(__name__ == "__main__"):
  exercise()
  exercise_symmetry_sum()
  print("OK")
