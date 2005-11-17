from iotbx.kriber import strudat
from iotbx.option_parser import iotbx_option_parser
from cctbx.array_family import flex
from scitbx import matrix
from libtbx.test_utils import approx_equal
import libtbx.load_env
from cStringIO import StringIO
import math
import sys, os

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
    return self.site_constraints.independent_gradients(
      all_gradients=all_gradients)

  def d2f_d_site(self):
    all_curvatures = -math.cos(self.h.dot(self.site)) * self.h.outer_product()
    return self.site_constraints.independent_curvatures(
      all_curvatures=flex.double(all_curvatures.as_list_of_lists()))

def df_d_site_finite(h, site_constraints, independent_params, eps=1.e-8):
  result = flex.double()
  independent_params_eps = independent_params.deep_copy()
  for ip in xrange(len(independent_params)):
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
  independent_params_eps = independent_params.deep_copy()
  for ip in xrange(len(independent_params)):
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
  return result

def exercise(structure, out):
  structure.show_summary(f=out)
  unit_cell = structure.unit_cell()
  for i_scatterer,scatterer in enumerate(structure.scatterers()):
    site = scatterer.site
    site_symmetry_ops = structure.site_symmetry_table().get(i_scatterer)
    assert approx_equal(site_symmetry_ops.special_op()*site, site)
    site_constraints = site_symmetry_ops.site_constraints(
      initialize_gradient_handling=True)
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
    for i_trial in xrange(10):
      shifted_params = independent_params + (flex.random_double(size=nip)-0.5)
      shifted_site = list(site_constraints.all_params(
        independent_params=shifted_params))
      assert approx_equal(
        site_symmetry_ops.special_op()*shifted_site, shifted_site)
      if (nip == 0 or unit_cell.distance(shifted_site, site) > 0.1):
        break
    else:
      raise AssertionError
    independent_gradients = site_constraints.independent_gradients(
      all_gradients=[0,0,0])
    assert len(independent_gradients) == nip
    assert approx_equal(independent_gradients, [0]*nip)
    independent_curvatures = site_constraints.independent_curvatures(
      all_curvatures=flex.double([[0,0,0],[0,0,0],[0,0,0]]))
    assert independent_curvatures.size() == nip**2
    assert approx_equal(independent_curvatures, [0]*(nip**2))
    #
    h = math.pi * flex.random_double(size=3)*2-1
    grads_fin = df_d_site_finite(
      h=h,
      site_constraints=site_constraints,
      independent_params=independent_params)
    print >> out, "grads_fin:", list(grads_fin)
    ca = cos_alpha(
      h=h,
      site_constraints=site_constraints,
      independent_params=independent_params)
    grads_ana = ca.df_d_site()
    print >> out, "grads_ana:", list(grads_ana)
    assert approx_equal(grads_ana, grads_fin)
    curvs_fin = d2f_d_site_finite(
      h=h,
      site_constraints=site_constraints,
      independent_params=independent_params)
    print >> out, "curvs_fin:", list(curvs_fin)
    curvs_ana = ca.d2f_d_site()
    print >> out, "curvs_ana:", list(curvs_ana)
    assert approx_equal(curvs_ana, curvs_fin)
  print >> out

def process_zeolite_atlas(args):
  command_line = (iotbx_option_parser()
    .option(None, "--tag",
      action="store",
      type="string",
      dest="tag")
    .option(None, "--verbose",
      action="store_true",
      dest="verbose")
  ).process(args=args)
  atlas_file = libtbx.env.find_in_repositories(
    relative_path="regression/misc/strudat_zeolite_atlas",
    test=os.path.isfile)
  if (atlas_file is None):
    print "Skipping process_zeolite_atlas(): input file not available"
    return
  if (command_line.options.verbose):
    out = sys.stdout
  else:
    out = StringIO()
  all_entries = strudat.read_all_entries(open(atlas_file))
  for i,entry in enumerate(all_entries.entries):
    structure = entry.as_xray_structure()
    if (command_line.options.tag is not None):
      if (command_line.options.tag != entry.tag):
        continue
    print >> out, "strudat tag:", entry.tag
    exercise(structure=structure, out=out)

def run():
  process_zeolite_atlas(sys.argv[1:])
  print "OK"

if (__name__ == "__main__"):
  run()
