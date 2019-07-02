from __future__ import absolute_import, division, print_function
from six.moves import zip
from six.moves import cStringIO as StringIO
from libtbx.test_utils import approx_equal, Exception_expected
from scitbx.random import variate, bernoulli_distribution
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx import adptbx
from iotbx import shelx

def exercise(flags, space_group_info):
  # Prepare a structure compatible with the ShelX model
  xs = random_structure.xray_structure(
    space_group_info,
    elements="random",
    n_scatterers=10,
    use_u_iso=True, random_u_iso=True,
    use_u_aniso=True)
  xs.apply_symmetry_sites()
  xs.apply_symmetry_u_stars()
  for isotropic, sc in zip(variate(bernoulli_distribution(0.4)),
                                      xs.scatterers()):
    sc.flags.set_grad_site(True)
    if isotropic:
      sc.flags.set_use_u_iso(True)
      sc.flags.set_use_u_aniso(False)
      sc.flags.set_grad_u_iso(True)
    else:
      sc.flags.set_use_u_iso(False)
      sc.flags.set_use_u_aniso(True)
      sc.flags.set_grad_u_aniso(True)

  not_origin_centric = (
            xs.space_group().is_centric()
    and not xs.space_group().is_origin_centric())

  try:
    ins = list(shelx.writer.generator(
      xs,
      full_matrix_least_squares_cycles=4,
      weighting_scheme_params=(0,0),
      sort_scatterers=False))
  except AssertionError:
    if (not_origin_centric):
      print ("Omitted %s\n  because it is centric but not origin centric"
             % xs.space_group().type().hall_symbol())
      return
    raise
  else:
    if (not_origin_centric):
      raise Exception_expected

  ins = StringIO("".join(ins))
  xs1 = xs.from_shelx(file=ins)
  xs.crystal_symmetry().is_similar_symmetry(
    xs1.crystal_symmetry(),
    relative_length_tolerance=1e-3,
    absolute_angle_tolerance=1e-3)
  uc = xs.unit_cell()
  uc1 = xs1.unit_cell()
  for sc, sc1 in zip(xs.scatterers(), xs1.scatterers()):
    assert sc.label.upper() == sc1.label.upper()
    assert sc.scattering_type == sc1.scattering_type
    assert sc.flags.bits == sc1.flags.bits
    assert approx_equal(sc.site, sc1.site, eps=1e-6)
    if sc.flags.use_u_iso():
      assert approx_equal(sc.u_iso, sc1.u_iso, eps=1e-5)
    else:
      assert approx_equal(adptbx.u_star_as_u_cif(uc, sc.u_star),
                          adptbx.u_star_as_u_cif(uc1, sc1.u_star),
                          eps=1e-5)

def run():
  import sys
  debug_utils.parse_options_loop_space_groups(
    sys.argv[1:], exercise)

if __name__ == '__main__':
  run()
