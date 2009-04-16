import os
import libtbx.load_env
from iotbx import shelx
from iotbx.reflection_file_reader import any_reflection_file
from cctbx import xray
from smtbx import refinement
from scitbx.matrix import *
import cctbx.euclidean_model_matching as emma

def exercise():
  smtbx_regression = libtbx.env.find_in_repositories('smtbx_regression')
  if smtbx_regression is None:
    print "Skipped (check out smtbx_regression if you need to run this)"
    return
  xs_name = os.path.join(smtbx_regression, 'Co110')
  xs = xray.structure.from_shelx(filename=xs_name+'.res',
                                 set_grad_flags=True)
  xs = xs.select(xs.element_selection('H'), negate=True)
  uc = xs.unit_cell()
  scatt_named = dict([ (s.label, s) for s in xs.scatterers() ])
  for s in xs.scatterers():
    f = s.flags
    assert f.grad_site() == f.grad_u_aniso() == f.use_u_aniso() == True
    assert f.grad_occupancy() == f.grad_fp() == f.grad_fdp() == False
    assert f.grad_u_iso() == f.use_u_iso() == False
  bragg = any_reflection_file(xs_name + ".hkl=hklf+ins/res")
  f_o_sq = bragg.as_miller_arrays()[0].merge_equivalents().array()
  f_o = f_o_sq.f_sq_as_f()
  xs0 = xs.deep_copy_scatterers()
  c28 = scatt_named['C28']
  c27 = scatt_named['C27']
  c26 = scatt_named['C26']
  n22 = scatt_named['N22']
  u = (row(c26.site) - row(c27.site)).cross(row(n22.site) - row(c27.site))
  u /= uc.length(u)
  c28.site = row(c28.site) + 1.*u
  print "Before refinement:"
  print "------------------"
  matches = emma.model_matches(
    xs0.as_emma_model(), xs.as_emma_model(),
    break_if_match_with_no_singles=False).refined_matches
  matches[0].show()
  u0 = c28.u_star
  fit = refinement.minimization.lbfgs(
    #target_functor=xray.unified_least_squares_residual(
      #f_o_sq,
      #weighting=xray.weighting_schemes.shelx_weighting()),
    #target_functor=xray.unified_least_squares_residual(f_o_sq),
    target_functor=xray.unified_least_squares_residual(f_o),
    #target_functor=xray.unified_least_squares_residual(
      #f_o, weighting=xray.weighting_schemes.pure_statistical_weighting()),
    xray_structure=xs
  )
  print "After refinement:"
  print "-----------------"
  matches = emma.model_matches(
    xs0.as_emma_model(), xs.as_emma_model(),
    break_if_match_with_no_singles=False).refined_matches
  matches[0].show()
  print "delta u(C28)=", row(c28.u_star) - row(u0)

def run():
  exercise()
  print 'OK'

if __name__ == '__main__':
  run()
