from __future__ import absolute_import, division, print_function

from cctbx import crystal, xray
from cctbx.array_family import flex
from smtbx.refinement import constraints, model
from scitbx import matrix
import math
from libtbx.test_utils import approx_equal
import os
from smtbx.regression.test_data import fnames

def exercise_basics():
  # construct a simple structure whose sites and u_iso's are to be refined
  xs = xray.structure(
    crystal_symmetry=crystal.symmetry(
      unit_cell=(10, 10, 10, 90, 90, 90),
      space_group_symbol='hall: P 1'),
    scatterers=flex.xray_scatterer((
      xray.scatterer('C0', site=(0, -1/2, 0), u=0.1),
      xray.scatterer('C1', site=(0,  1/2, 0), u=0.2),
      xray.scatterer('C0a', site=( 1/2, 0, 0), u=0.1),
      xray.scatterer('C1a', site=(-1/2, 0, 0), u=0.2),
    )))
  for sc in xs.scatterers():
    sc.flags.set_grad_site(True)
    sc.flags.set_grad_u_iso(True)

  # copy the original structure as a reference to test against later
  xs_ref = xs.deep_copy_scatterers()

  # mess up the symmetries that the forthcoming constraints shall impose
  c0, c1, c0a, c1a = xs.scatterers()
  c0a.site = (0, 0, 0)
  c1a.site = (1/2, 1/2, 1/2)
  c0a.u_iso = 0.5
  c1a.u_iso = 0.6

  # construct a reparametrisation for the following constraints:
  # (C0a, C1a) is the image of (C0, C1) through a rotation of 90 degrees
  # about the z-axis
  r = constraints.ext.reparametrisation(xs.unit_cell())
  sc_params = constraints.shared_scatterer_parameters(xs.scatterers())
  c0_site_param  = r.add(constraints.independent_site_parameter, c0)
  sc_params[0].site = c0_site_param
  c1_site_param = r.add(constraints.independent_site_parameter, c1)
  sc_params[1].site = c1_site_param
  move_param = r.add(
    constraints.independent_small_6_vector_parameter,
    (0,0,0, 0, 0, math.pi/2))
  c0a_c1a_site_param = r.add(
    constraints.same_group_xyz,
    scatterers=(c0a, c1a),
    sites=(c0_site_param, c1_site_param),
    alignment_matrix=matrix.identity(3),
    shifts_and_angles=move_param)
  sc_params[2].site = sc_params[3].site = c0a_c1a_site_param
  c0_u_iso_param = r.add(constraints.independent_u_iso_parameter, c0)
  sc_params[0].u = c0_u_iso_param
  c1_u_iso_param = r.add(constraints.independent_u_iso_parameter, c1)
  sc_params[1].u = c1_u_iso_param
  c0a_c1a_u_iso_param = r.add(
    constraints.same_group_u_iso,
    scatterers=(c0a, c1a),
    u_isos=(c0_u_iso_param, c1_u_iso_param))
  sc_params[2].u = sc_params[3].u = c0a_c1a_u_iso_param
  r.finalise()

  # put the reparametrisation to work
  r.linearise()
  r.store()
  c0_ref, c1_ref, c0a_ref, c1a_ref = xs_ref.scatterers()
  assert approx_equal(c0a.site, c0a_ref.site, eps=1e-12)
  assert approx_equal(c1a.site, c1a_ref.site, eps=1e-12)
  assert approx_equal(c0a.u_iso, c0a_ref.u_iso, eps=1e-12)
  assert approx_equal(c1a.u_iso, c1a_ref.u_iso, eps=1e-12)

  # check that origin fixing restraints work in the presence of
  # that reparametrisation
  # this is a regression test as they used not to (bug reported by Oleg)
  from smtbx.refinement.restraints import origin_fixing_restraints
  from scitbx import lstbx
  orig_fixing = origin_fixing_restraints.homogeneous_weighting(xs.space_group())
  normal_eqn = lstbx.normal_eqns.ext.linear_ls(n_parameters=(3 + 1)*2 + 6)
  jacobian_transpose_matching_grad_fc = r.jacobian_transpose_matching(
    sc_params.mapping_to_grad_fc())
  orig_fixing.add_to(normal_eqn,
                     jacobian_transpose_matching_grad_fc,
                     sc_params)

def exercise_real_life_structure():
  working_dir = os.path.dirname(__file__)
  res = fnames.sucrose_p1_res
  xs = xray.structure.from_shelx(filename=res)
  fo_sq = xs.structure_factors(
    d_min=0.5, algorithm='direct').f_calc().intensities()
  fo_sq.set_sigmas(fo_sq.data()*0.05)
  m = model.from_shelx(res, fo_sq=fo_sq)
  ls = m.least_squares()
  ls.build_up()
  # wip

if __name__ == '__main__':
  exercise_basics()
  exercise_real_life_structure()
  print('OK')
