from __future__ import absolute_import, division, print_function
import sys
from libtbx import group_args
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
from six.moves import cStringIO as StringIO
from iotbx.map_model_manager import map_model_manager

# Averaging of local anisotropic scaling information over locations
#  (map_model_manager._average_scale_factor_info_over_xyz) and the
#  "most anisotropy" section of _summarize_scale_factor_info.
# Each case can be run alone:
#  libtbx.python tst_map_model_manager_aniso_averaging.py complete_unchanged

ANISO_NOTE = "Scaling information for location with most anisotropy"

TSF = (1., 0.9, 0.8, 0.7)    # target_scale_factors
CC = (0.9, 0.8, 0.7, 0.6)    # cc_list
RMS_FO = (100., 80., 50., 30.)
U_CART = (1., 2., 3., 0.1, 0.2, 0.3)

def get_mmm():
  mmm = map_model_manager()
  mmm.generate_map(d_min = 3, wrapping = False)
  log = StringIO()
  mmm.set_log(log)
  return mmm, log

def u_cart(scale = 1.):
  return tuple([scale * x for x in U_CART])

def scaled(values, scale):
  return flex.double([scale * x for x in values])

def make_si(scale = 1., rms_fo_list = True):
  sthol2 = flex.double([0.25 / (d * d) for d in (8., 6., 4., 3.)])
  return group_args(
    group_args_type = 'synthetic scaling_info',
    target_sthol2 = sthol2,
    cc_list = scaled(CC, scale),
    rms_fo_list = scaled(RMS_FO, scale) if rms_fo_list else None,
    rms_fc_list = flex.double([110., 85., 55., 32.]),
    target_scale_factors = scaled(TSF, scale),
   )

def make_location(scale = 1., n_directions = 2, complete = True,
     empty = False, fo_u_cart = None, overall_si = True):
  '''Synthetic location; all values are proportional to scale'''
  if empty:
    scaling_info_list = []
    direction_vectors = []
  else:
    scaling_info_list = [make_si(scale = scale, rms_fo_list = complete)
       for i in range(n_directions)]
    direction_vectors = [(1., 0., 0.), (0., 1., 0.), (0., 0., 1.)][
       :n_directions]
  return group_args(
    group_args_type = 'synthetic scaling_group_info',
    direction_vectors = direction_vectors,
    scaling_info_list = scaling_info_list,
    overall_si = make_si(scale = scale) if overall_si else None,
    overall_scale = scaled(TSF, scale),
    aa_b_cart_as_u_cart = u_cart(1. * scale),
    bb_b_cart_as_u_cart = u_cart(2. * scale),
    ss_b_cart_as_u_cart = u_cart(3. * scale),
    uu_b_cart_as_u_cart = u_cart(4. * scale),
    fo_b_cart_as_u_cart = tuple(fo_u_cart) if fo_u_cart is not None
       else u_cart(5. * scale),
   )

XYZ = [(10., 11., 12.), (20., 21., 22.), (30., 31., 32.)]

def make_scale_factor_info(location_list):
  return group_args(
    group_args_type = 'synthetic scale_factor_info',
    value_list = location_list,
    xyz_list = flex.vec3_double(XYZ[:len(location_list)]),
   )

def check_si_values(si, scale, keys = ('target_scale_factors', 'cc_list',
     'rms_fo_list')):
  for key, values in (('target_scale_factors', TSF), ('cc_list', CC),
      ('rms_fo_list', RMS_FO)):
    if key in keys:
      assert approx_equal(list(getattr(si, key)), list(scaled(values, scale))
        ), (key, list(getattr(si, key)), list(scaled(values, scale)))

def check_b_cart(average, scale):
  for key, factor in (('aa_b_cart_as_u_cart', 1.),
     ('bb_b_cart_as_u_cart', 2.), ('ss_b_cart_as_u_cart', 3.),
     ('uu_b_cart_as_u_cart', 4.), ('fo_b_cart_as_u_cart', 5.)):
    assert approx_equal(list(average.get(key)), list(u_cart(factor * scale))
      ), (key, list(average.get(key)), list(u_cart(factor * scale)))

def tst_overall_si_partial_skip():
  '''Middle location lacks rms_fo_list and is skipped: overall_si average
     is over the two used locations (scales 1 and 4 -> 2.5; dividing the
     sum by all three locations would give 5/3)'''
  mmm, log = get_mmm()
  sfi = make_scale_factor_info([make_location(scale = 1.),
     make_location(scale = 2., complete = False),
     make_location(scale = 4.)])
  avg = mmm._average_scale_factor_info_over_xyz(sfi).value_list[0]
  for si in avg.scaling_info_list:
    check_si_values(si, 2.5)
  check_si_values(avg.overall_si, 2.5, keys = ('target_scale_factors',
    'cc_list', 'rms_fo_list'))

def tst_empty_location_skipped():
  '''Middle location has an empty scaling_info_list (other values set):
     it is not counted (scales 1 and 4 -> 2.5)'''
  mmm, log = get_mmm()
  sfi = make_scale_factor_info([make_location(scale = 1.),
     make_location(scale = 2., empty = True),
     make_location(scale = 4.)])
  avg = mmm._average_scale_factor_info_over_xyz(sfi).value_list[0]
  assert len(avg.scaling_info_list) == 2, avg.scaling_info_list
  for si in avg.scaling_info_list:
    check_si_values(si, 2.5, keys = ('target_scale_factors',))
  check_b_cart(avg, 2.5)

def tst_complete_unchanged():
  '''All locations complete: averages are plain means
     (scales 1, 2, 6 -> 3)'''
  mmm, log = get_mmm()
  sfi = make_scale_factor_info([make_location(scale = 1.),
     make_location(scale = 2.),
     make_location(scale = 6.)])
  avg = mmm._average_scale_factor_info_over_xyz(sfi).value_list[0]
  assert len(avg.scaling_info_list) == 2, avg.scaling_info_list
  for si in avg.scaling_info_list:
    check_si_values(si, 3.)
  check_si_values(avg.overall_si, 3.)
  check_b_cart(avg, 3.)

# fo_b_cart_as_u_cart (U11, U22, U33, U12, U13, U23) for most_anisotropy
FO_ISO = (2., 2., 2., 0., 0., 0.)      # A: isotropic, spread 0, mean |U| 1.000
FO_ANISO = (4., 0.2, 0.2, 0., 0., 0.)  # B: spread 3.8, mean |U| 0.733
FO_ROTATED = (1.25, 1.25, 0., 1.25, 0., 0.) # C: diag(2.5,0,0) rotated 45
                                       #  deg about z; spread 2.5, mean 0.625

def tst_most_anisotropy_section():
  '''Location with largest eigenvalue spread of fo_b_cart_as_u_cart (B,
     the middle one) is reported, with that spread and its xyz. The
     isotropic location A has the largest mean |component| and must not
     be picked.'''
  mmm, log = get_mmm()
  sfi = make_scale_factor_info([make_location(scale = 1., fo_u_cart = fo)
     for fo in (FO_ISO, FO_ANISO, FO_ROTATED)])
  mmm._summarize_scale_factor_info(sfi)
  text = log.getvalue()
  start = text.find(ANISO_NOTE)
  assert start > -1, text
  heading = text[start:].splitlines()[0]
  # XYZ printed by the per-location summary that follows the heading
  rest = text[start + len(heading):]
  xyz_start = rest.find("XYZ:")
  assert xyz_start > -1, text
  xyz_line = rest[xyz_start:].splitlines()[0]
  assert xyz_line == "XYZ: (%.3f, %.3f, %.3f)" %(XYZ[1]), (xyz_line, XYZ[1])
  assert heading.find("(Anisotropy (eigenvalue spread of U): 3.80 )") > -1,\
     heading

  # Spread is zero for isotropic U and does not change under rotation
  assert approx_equal(mmm._get_anisotropy_spread(FO_ISO), 0., eps = 1.e-6)
  assert approx_equal(mmm._get_anisotropy_spread((1., -1., 0., 0., 0., 0.)),
     mmm._get_anisotropy_spread((0., 0., 0., 1., 0., 0.)), eps = 1.e-6)
  assert approx_equal(mmm._get_anisotropy_spread(FO_ANISO), 3.8, eps = 1.e-6)
  assert approx_equal(mmm._get_anisotropy_spread(FO_ROTATED), 2.5,
     eps = 1.e-6)

def tst_overall_si_none():
  '''Single-location summary with overall_si None completes'''
  mmm, log = get_mmm()
  sfi = make_scale_factor_info([make_location(scale = 1.,
     overall_si = False)])
  mmm._summarize_scale_factor_info(sfi)
  text = log.getvalue()
  assert text.find("Scale factors by direction vector") > -1, text

def run(args = None):
  all_tests = [
     ('overall_si_partial_skip', tst_overall_si_partial_skip),
     ('empty_location_skipped', tst_empty_location_skipped),
     ('complete_unchanged', tst_complete_unchanged),
     ('most_anisotropy_section', tst_most_anisotropy_section),
     ('overall_si_none', tst_overall_si_none),
    ]
  if args:
    tests = [(name, f) for name, f in all_tests if name in args]
  else:
    tests = all_tests
  for name, f in tests:
    print("Running %s" %(name))
    f()
  print("OK")

if __name__ == "__main__":
  run(sys.argv[1:])
