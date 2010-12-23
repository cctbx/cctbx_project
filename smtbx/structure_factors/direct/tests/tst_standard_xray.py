from __future__ import division

from cctbx.array_family import flex
from cctbx import sgtbx, xray, crystal, miller, eltbx
import cctbx.eltbx.wavelengths
import smtbx.development
from cctbx.development.space_group_option_parser\
     import space_group_option_parser
import smtbx.structure_factors.direct as structure_factors
from libtbx.test_utils import approx_equal
from scitbx.math import approx_equal_relatively
import libtbx.utils
import random
from itertools import islice, izip
from scitbx import matrix
from scitbx.math import median_statistics


class test_case(object):

  def __init__(self, inelastic_scattering):
    self.inelastic_scattering = inelastic_scattering

  def random_structure(self, space_group_info, set_grads):
    xs = smtbx.development.random_xray_structure(
      space_group_info,
      elements="random",
      n_scatterers=5,
      u_iso_xor_u_aniso=False,
      use_u_aniso=True,
      use_u_iso=True,
      random_u_iso=True,
    )
    if self.inelastic_scattering:
      xs.set_inelastic_form_factors(
        photon=eltbx.wavelengths.characteristic('Mo'),
        table="sasaki")
    for sc in xs.scatterers():
      (
        sc.flags.set_grad_site(set_grads)
                .set_grad_u_iso(set_grads).set_grad_u_aniso(set_grads)
                .set_grad_occupancy(set_grads)
      )
      assert sc.flags.use_u_iso()
      assert sc.flags.use_u_aniso()
    return xs

  def miller_indices(self, space_group_info):
    space_group = space_group_info.group()
    return flex.miller_index(miller.index_generator(
      space_group.type(),
      anomalous_flag=True,
      max_index=(10, 10, 10)))

  def exercise(self, xray_structure=None, space_group_info=None,
               verbose=False, fixed_random_seed=True, **kwds):
    assert [xray_structure, space_group_info].count(None) == 1
    if xray_structure is None:
      self.xs = self.random_structure(space_group_info, set_grads=True)
    else:
      self.xs = xray_structure

    if fixed_random_seed:
      random.seed(1)
      flex.set_random_seed(1)

    self.do_exercise(verbose)


class consistency_test_cases(test_case):

  def __init__(self, n_directions, inelastic_scattering):
    super(consistency_test_cases, self).__init__(inelastic_scattering)
    self.n_directions = n_directions

  def structures_forward(self, xs, xs_forward, eta_norm):
    while True:
      direction = flex.random_double(xs.n_parameters())
      direction /= direction.norm()
      eta = eta_norm * direction

      i = 0
      for sc_forward, sc in izip(xs_forward.scatterers(), xs.scatterers()):
        eta_site = matrix.col(eta[i:i+3])
        eta_iso = eta[i+3]
        eta_aniso = matrix.col(eta[i+4:i+10])
        eta_occ = eta[i+10]
        i += 11

        sc_forward.site = matrix.col(sc.site) + eta_site
        sc_forward.u_iso = sc.u_iso + eta_iso
        sc_forward.u_star = matrix.col(sc.u_star) + eta_aniso
        sc_forward.occupancy = sc.occupancy + eta_occ
      yield direction
  structures_forward = classmethod(structures_forward)

  def do_exercise(self, verbose=False):
    xs = self.xs
    sg = xs.space_group_info().group()
    origin_centric_case = sg.is_origin_centric()

    indices = self.miller_indices(xs.space_group_info())
    f = structure_factors.f_calc_modulus_squared(xs)
    f1 = structure_factors.f_calc_modulus_squared(xs)

    for h in indices:
      f.linearise(h)
      fl = f.f_calc
      f1.evaluate(h)
      fe = f1.f_calc
      assert f1.grad_f_calc is None
      assert approx_equal_relatively(fe, fl, relative_error=1e-12), (fe, fl)

    if (xs.space_group().is_origin_centric() and not self.inelastic_scattering):
      for h in indices:
        f.linearise(h)
        assert f.f_calc.imag == 0
        assert flex.imag(f.grad_f_calc).all_eq(0)

    eta = 1e-8
    xs_forward = xs.deep_copy_scatterers()
    f_forward = structure_factors.f_calc_modulus_squared(xs_forward)

    deltas = flex.double()
    for direction in islice(self.structures_forward(xs, xs_forward, eta),
                            self.n_directions):
      for h in indices:
        f.linearise(h)
        assert approx_equal(abs(f.f_calc)**2, f.observable)
        f_forward.linearise(h)
        diff_num = (f_forward.observable - f.observable) / eta
        diff = f.grad_observable.dot(direction)
        delta = abs(1 - diff/diff_num)
        deltas.append(delta)
    stats = median_statistics(deltas)
    tol = 1e-3
    assert stats.median < tol, (str(space_group_info), stats.median)
    assert stats.median_absolute_deviation < tol, (str(space_group_info),
      stats.median_absolute_deviation)

class smtbx_against_cctbx_test_case(test_case):

  def do_exercise(self, verbose=False):
    xs = self.xs
    indices = self.miller_indices(xs.space_group_info())
    cctbx_structure_factors = xray.structure_factors.from_scatterers_direct(
      xray_structure=xs,
      miller_set=miller.set(
        crystal.symmetry(unit_cell=xs.unit_cell(),
                         space_group_info=xs.space_group_info()),
        indices))
    f = structure_factors.f_calc_modulus_squared(xs)
    for h, fc in cctbx_structure_factors.f_calc():
      f.linearise(h)
      if fc == 0:
        assert f.f_calc == 0
      else:
        delta = abs((f.f_calc - fc)/fc)
        assert delta < 1e-6

class custom_vs_std_test_case(test_case):

  def __init__(self, *args, **kwds):
    self.has_printed_header = False
    test_case.__init__(self, *args, **kwds)

  def do_exercise(self, verbose=False):
    xs = self.xs
    indices = self.miller_indices(xs.space_group_info())
    exp_i_2pi_functor = cctbx.math_module.cos_sin_table(1024)
    custom_fc_sq = (
      structure_factors.f_calc_modulus_squared(
        xs, exp_i_2pi_functor))
    std_fc_sq = (
      structure_factors.f_calc_modulus_squared(xs))
    deltas = flex.double()
    for h in indices:
      custom_fc_sq.linearise(h)
      std_fc_sq.linearise(h)
      deltas.append(abs(custom_fc_sq.f_calc - std_fc_sq.f_calc)
                    /abs(std_fc_sq.f_calc))
    stats = median_statistics(deltas)
    if verbose:
      if not self.has_printed_header:
        print "f_calc and sin/cos: |tabulated - std|/|std|"
        print "median & median absolute deviation"
        self.has_printed_header = True
      print "%s: %.12g +/- %.12g" % (xs.space_group_info().type().hall_symbol(),
                                     stats.median,
                                     stats.median_absolute_deviation)
    assert stats.median < 0.01, (str(xs.space_group_info()), stats.median)
    assert stats.median_absolute_deviation < 0.005, (
      str(xs.space_group_info()), stats.median_absolute_deviation)


class f_vs_f_sq_test_case(test_case):

  skip_warning = True

  def do_exercise(self, verbose=False):
    import boost.python
    vers = boost.python.gcc_version()
    if (vers is not None and str(vers).startswith("4.1.")):
      if self.__class__.skip_warning:
        print "Skip F vs F^2 test because platform with GCC %s" % vers
        self.__class__.skip_warning = False
      return
    xs = self.xs
    indices = self.miller_indices(xs.space_group_info())
    f = structure_factors.f_calc_modulus(xs)
    f_sq = structure_factors.f_calc_modulus_squared(xs)
    for h in indices:
      f.linearise(h)
      f_sq.linearise(h)
      assert approx_equal_relatively(f.observable**2,
                                     f_sq.observable,
                                     relative_error=1e-12)
      grad_f_sq = f_sq.grad_observable
      two_f_grad_f = (2*f.observable*f.grad_observable)
      assert two_f_grad_f.all_approx_equal_relatively(grad_f_sq,
                                                      relative_error=1e-12)


def exercise_trigonometric_ff():
  from math import cos, sin,pi
  sgi = sgtbx.space_group_info("P1")
  cs = sgi.any_compatible_crystal_symmetry(volume=1000)
  miller_set = miller.build_set(cs, anomalous_flag=False, d_min=1)
  miller_set = miller_set.select(flex.random_double(miller_set.size()) < 0.2)
  for i in xrange(5):
    sites = flex.random_double(9)
    x1, x2, x3 = (matrix.col(sites[:3]),
                  matrix.col(sites[3:6]),
                  matrix.col(sites[6:]))
    xs = xray.structure(crystal.special_position_settings(cs))
    for x in (x1, x2, x3):
      sc = xray.scatterer(site=x, scattering_type="const")
      sc.flags.set_grad_site(True)
      xs.add_scatterer(sc)
    f_sq = structure_factors.f_calc_modulus_squared(xs)
    for h in miller_set.indices():
      h = matrix.col(h)
      phi1, phi2, phi3 = 2*pi*h.dot(x1), 2*pi*h.dot(x2), 2*pi*h.dot(x3)
      fc_mod_sq = 3 + 2*(cos(phi1 - phi2) + cos(phi2 - phi3) + cos(phi3 - phi1))
      g = []
      g.extend( -2*(sin(phi1 - phi2) - sin(phi3 - phi1))*2*pi*h )
      g.extend( -2*(sin(phi2 - phi3) - sin(phi1 - phi2))*2*pi*h )
      g.extend( -2*(sin(phi3 - phi1) - sin(phi2 - phi3))*2*pi*h )
      grad_fc_mod_sq = g

      f_sq.linearise(h)
      assert approx_equal(f_sq.observable, fc_mod_sq)
      assert approx_equal(f_sq.grad_observable, grad_fc_mod_sq)

def run(args):
  libtbx.utils.show_times_at_exit()
  parser = space_group_option_parser()
  parser.option('-x', '--xray_structure_pickle', default=None)
  parser.option(None, '--fixed_random_seed', default=True)
  commands = parser.process(args)

  n_directions = 2

  if hasattr(commands.options, 'xray_structure_pickle'):
    from libtbx import easy_pickle
    xs = easy_pickle.load(commands.options.xray_structure_pickle)

    t = consistency_test_cases(n_directions, inelastic_scattering=False)
    t.exercise(xray_structure=xs)

  else:
    exercise_trigonometric_ff()

    t = f_vs_f_sq_test_case(inelastic_scattering=False)
    commands.loop_over_space_groups(t.exercise)

    t = f_vs_f_sq_test_case(inelastic_scattering=True)
    commands.loop_over_space_groups(t.exercise)

    t = custom_vs_std_test_case(inelastic_scattering=None)
    commands.loop_over_space_groups(t.exercise)

    t = smtbx_against_cctbx_test_case(inelastic_scattering=False)
    commands.loop_over_space_groups(t.exercise)

    t = consistency_test_cases(n_directions, inelastic_scattering=False)
    commands.loop_over_space_groups(t.exercise)

    t = consistency_test_cases(n_directions, inelastic_scattering=True)
    commands.loop_over_space_groups(t.exercise)

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
