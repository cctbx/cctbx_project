from __future__ import absolute_import, division, print_function
from iotbx.kriber import strudat
from iotbx.option_parser import iotbx_option_parser
from cctbx.array_family import flex
from scitbx import matrix
from libtbx.utils import format_cpu_times
from libtbx.test_utils import approx_equal
import libtbx.load_env
from six.moves import cStringIO as StringIO
import math
import sys, os
from six.moves import range

flex.set_random_seed(0)

class cos_alpha:

  def __init__(self, h, site_constraints, independent_params):
    self.h = matrix.row(h)
    self.site_constraints = site_constraints
    self.independent_params = independent_params
    self.site = matrix.col(self.site_constraints.all_params(
      independent_params=self.independent_params))

  def f(self):
    return math.cos(self.h.dot(self.site))

  def df_d_site(self):
    all_gradients = -math.sin(self.h.dot(self.site)) * self.h
    return flex.double(self.site_constraints.independent_gradients(
      all_gradients=flex.double(all_gradients)))

  def d2f_d_site(self):
    all_curvatures = -math.cos(self.h.dot(self.site)) * self.h.outer_product()
    return self.site_constraints.independent_curvatures(
      all_curvatures=flex.double(
        all_curvatures.as_list_of_lists()).matrix_symmetric_as_packed_u())

def df_d_site_finite(h, site_constraints, independent_params, eps=1.e-8):
  result = flex.double()
  independent_params_eps = list(independent_params)
  for ip in range(len(independent_params)):
    vs = []
    for signed_eps in [eps, -eps]:
      independent_params_eps[ip] = independent_params[ip] + signed_eps
      ca = cos_alpha(
        h=h,
        site_constraints=site_constraints,
        independent_params=independent_params_eps)
      vs.append(ca.f())
    result.append((vs[0]-vs[1])/(2*eps))
    independent_params_eps[ip] = independent_params[ip]
  return result

def d2f_d_site_finite(h, site_constraints, independent_params, eps=1.e-8):
  result = flex.double()
  independent_params_eps = list(independent_params)
  for ip in range(len(independent_params)):
    vs = []
    for signed_eps in [eps, -eps]:
      independent_params_eps[ip] = independent_params[ip] + signed_eps
      ca = cos_alpha(
        h=h,
        site_constraints=site_constraints,
        independent_params=independent_params_eps)
      vs.append(ca.df_d_site())
    result.extend((vs[0]-vs[1])/(2*eps))
    independent_params_eps[ip] = independent_params[ip]
  np = len(independent_params)
  result.reshape(flex.grid(np,np))
  return result.matrix_symmetric_as_packed_u(relative_epsilon=1.e-5)

def exercise(structure, out):
  structure.show_summary(f=out)
  unit_cell = structure.unit_cell()
  for i_scatterer,scatterer in enumerate(structure.scatterers()):
    site = scatterer.site
    site_symmetry_ops = structure.site_symmetry_table().get(i_scatterer)
    assert approx_equal(site_symmetry_ops.special_op()*site, site)
    site_constraints = site_symmetry_ops.site_constraints()
    assert site_constraints.row_echelon_form().focus()[0] \
        == site_constraints.n_dependent_params()
    nip = site_constraints.n_independent_params()
    if (site_symmetry_ops.n_matrices() == 1):
      assert nip == 3
      assert site_constraints.n_dependent_params() == 0
    else:
      assert nip < 3
      assert site_constraints.n_dependent_params() > 0
    independent_params = site_constraints.independent_params(all_params=site)
    all_params = site_constraints.all_params(
      independent_params=independent_params)
    assert approx_equal(all_params, site)
    for i_trial in range(10):
      shifted_params = flex.double(independent_params) \
                     + (flex.random_double(size=nip)-0.5)
      shifted_site = site_constraints.all_params(
        independent_params=list(shifted_params))
      assert approx_equal(
        site_symmetry_ops.special_op()*shifted_site, shifted_site)
      if (nip == 0 or unit_cell.distance(shifted_site, site) > 0.1):
        break
    else:
      raise AssertionError
    independent_gradients = site_constraints.independent_gradients(
      all_gradients=flex.double(3, 0))
    assert len(independent_gradients) == nip
    assert approx_equal(independent_gradients, [0]*nip)
    independent_curvatures = site_constraints.independent_curvatures(
      all_curvatures=flex.double(6, 0))
    assert len(independent_curvatures) == nip*(nip+1)//2
    assert approx_equal(independent_curvatures, [0]*(nip*(nip+1)//2))
    #
    h = math.pi * flex.random_double(size=3)*2-1
    grads_fin = df_d_site_finite(
      h=h,
      site_constraints=site_constraints,
      independent_params=independent_params)
    print("grads_fin:", list(grads_fin), file=out)
    ca = cos_alpha(
      h=h,
      site_constraints=site_constraints,
      independent_params=independent_params)
    grads_ana = ca.df_d_site()
    print("grads_ana:", list(grads_ana), file=out)
    assert approx_equal(grads_ana, grads_fin)
    curvs_fin = d2f_d_site_finite(
      h=h,
      site_constraints=site_constraints,
      independent_params=independent_params)
    print("curvs_fin:", list(curvs_fin), file=out)
    curvs_ana = ca.d2f_d_site()
    print("curvs_ana:", list(curvs_ana), file=out)
    assert approx_equal(curvs_ana, curvs_fin)
  print(file=out)

def exercise_shake_sites_in_place(structure):
  for target_difference in [0, 0.7, 1.0, 1.3]:
    for selection in [None,
                      flex.random_bool(
                        size=structure.scatterers().size(), threshold=0.5)]:
      if (target_difference == 0): # the structure is not changed in this case
        n_variable = 1 # any value except 0 (the flex.sum() below must be 0)
      else:
        if (selection is None):
          n_variable = structure.scatterers().size()
        else:
          n_variable = selection.count(True)
        n_variable -= structure.coordinate_degrees_of_freedom_counts(
          selection=selection)[0]
      shaken = structure.deep_copy_scatterers()
      try:
        shaken.shake_sites_in_place(
          rms_difference=target_difference,
          selection=selection)
      except RuntimeError as e:
        if (selection is not None):
          if (selection.count(True) == 0):
            assert str(e) == "No scatterers selected."
          else:
            assert str(e) \
                == "All selected scatterers are fixed on special positions."
        else:
          assert str(e) \
              == "All scatterers are fixed on special positions."
      else:
        assert approx_equal(
          (flex.sum((  structure.sites_cart()
                     - shaken.sites_cart()).dot()) / n_variable) ** 0.5,
          target_difference)
        #
        shaken = structure.deep_copy_scatterers()
        shaken.shake_sites_in_place(
          mean_distance=target_difference,
          selection=selection)
        assert approx_equal(
          flex.sum(flex.sqrt((  structure.sites_cart()
                              - shaken.sites_cart()).dot())) / n_variable,
          target_difference)

def process_zeolite_atlas(args):
  command_line = (iotbx_option_parser()
    .option(None, "--tag",
      action="store",
      type="string")
    .option(None, "--verbose",
      action="store_true")
  ).process(args=args)
  atlas_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/misc/strudat_zeolite_atlas",
    test=os.path.isfile)
  if (atlas_file is None):
    print("Skipping process_zeolite_atlas(): input file not available")
    return
  if (command_line.options.verbose):
    out = sys.stdout
  else:
    out = StringIO()
  with open(atlas_file) as f:
    all_entries = strudat.read_all_entries(f)
  for i,entry in enumerate(all_entries.entries):
    structure = entry.as_xray_structure()
    if (command_line.options.tag is not None):
      if (command_line.options.tag != entry.tag):
        continue
    print("strudat tag:", entry.tag, file=out)
    exercise(structure=structure, out=out)
    exercise_shake_sites_in_place(structure=structure)

def run():
  process_zeolite_atlas(sys.argv[1:])
  print(format_cpu_times())

if (__name__ == "__main__"):
  run()
