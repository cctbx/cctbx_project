from __future__ import absolute_import, division, print_function
import sys
from libtbx import group_args
from scitbx.array_family import flex
from six.moves import cStringIO as StringIO
from iotbx.map_model_manager import map_model_manager

# Summary of local anisotropic scaling information
#  (map_model_manager._summarize_scale_factor_info) when some or all
#  locations have incomplete scaling information.
# Each case can be run alone:
#  libtbx.python tst_map_model_manager_aniso_summary.py complete all_skipped

SKIP_NOTE = "No locations with complete scaling information"
NO_CC_NOTE = "No CC* values available for location with"
ANISO_NOTE = "Scaling information for location with most anisotropy"
CC_HEADING = "Estimated CC* by direction vector "

def get_mmm():
  mmm = map_model_manager()
  mmm.generate_map(d_min = 3, wrapping = False)
  log = StringIO()
  mmm.set_log(log)
  return mmm, log

def u_cart(scale = 1.):
  return tuple([scale * x for x in (1., 2., 3., 0.1, 0.2, 0.3)])

def make_si(shift = 0., cc_list = True, rms_fo_list = True,
    target_scale_factors = True):
  sthol2 = flex.double([0.25 / (d * d) for d in (8., 6., 4., 3.)])
  return group_args(
    group_args_type = 'synthetic scaling_info',
    target_sthol2 = sthol2,
    cc_list = flex.double([0.9, 0.8, 0.7, 0.6]) - shift
       if cc_list else None,
    rms_fo_list = flex.double([100., 80., 50., 30.])
       if rms_fo_list else None,
    rms_fc_list = flex.double([110., 85., 55., 32.]),
    target_scale_factors = flex.double([1., 0.9, 0.8, 0.7])
       if target_scale_factors else None,
   )

def make_location(shift = 0., n_directions = 2, complete = True,
     empty = False, cc_list = True, fo_scale = 5.):
  if empty:
    scaling_info_list = []
    direction_vectors = []
  else:
    scaling_info_list = [make_si(shift = shift, rms_fo_list = complete,
       cc_list = cc_list)
       for i in range(n_directions)]
    direction_vectors = [(1., 0., 0.), (0., 1., 0.), (0., 0., 1.)][
       :n_directions]
  return group_args(
    group_args_type = 'synthetic scaling_group_info',
    direction_vectors = direction_vectors,
    scaling_info_list = scaling_info_list,
    overall_si = make_si(shift = shift),
    overall_scale = flex.double([1., 0.9, 0.8, 0.7]),
    aa_b_cart_as_u_cart = u_cart(1.),
    bb_b_cart_as_u_cart = u_cart(2.),
    ss_b_cart_as_u_cart = u_cart(3.),
    uu_b_cart_as_u_cart = u_cart(4.),
    # Read by _get_scale_factor_info_most_anisotropy: the summary shows the
    #  location with largest eigenvalue spread of fo_b_cart_as_u_cart; on a
    #  tie the first location is used
    fo_b_cart_as_u_cart = u_cart(fo_scale),
   )

def make_scale_factor_info(location_list):
  return group_args(
    group_args_type = 'synthetic scale_factor_info',
    value_list = location_list,
    xyz_list = flex.vec3_double([(10., 10., 10.), (20., 20., 20.),
      (30., 30., 30.)][:len(location_list)]),
   )

def cc_heading_lines(text):
  return [line for line in text.splitlines() if line.find(CC_HEADING) > -1]

def tst_complete():
  '''Positive control: complete data at two locations gives full summary'''
  mmm, log = get_mmm()
  sfi = make_scale_factor_info([make_location(shift = 0.),
     make_location(shift = 0.2)])
  mmm._summarize_scale_factor_info(sfi)
  text = log.getvalue()
  assert text.find("Scaling information for location with highest " +
     "resolution (Mean CC*: 0.75 )") > -1, text
  assert text.find("Scaling information for location with lowest " +
     "resolution (Mean CC*: 0.55 )") > -1, text
  assert text.find("Summary of overall average scaling information") > -1, \
     text
  assert text.find(SKIP_NOTE) < 0, text
  assert text.find(NO_CC_NOTE) < 0, text
  # Overall summary printed: averaged CC* (0.65) appears in heading
  assert text.find("Estimated CC* by direction vector  (Mean CC*: 0.65)") \
     > -1, text
  assert text.find(ANISO_NOTE) > -1, text
  # Every CC* heading carries a mean when CC* values are present
  for line in cc_heading_lines(text):
    assert line.find("(Mean CC*:") > -1, text

def tst_all_skipped():
  '''Every location lacks rms_fo_list, so all are skipped in the average'''
  mmm, log = get_mmm()
  sfi = make_scale_factor_info([make_location(shift = 0., complete = False),
     make_location(shift = 0.2, complete = False)])
  mmm._summarize_scale_factor_info(sfi)
  text = log.getvalue()
  assert text.find("Summary of overall average scaling information") > -1, \
     text
  assert text.find(SKIP_NOTE) > -1, text
  # Per-location summaries still printed (cc_list is present)
  assert text.find("(Mean CC*: 0.75 )") > -1, text

def tst_no_cc_star():
  '''No location has any CC* values (empty scaling_info_list)'''
  mmm, log = get_mmm()
  sfi = make_scale_factor_info([make_location(empty = True),
     make_location(empty = True)])
  mmm._summarize_scale_factor_info(sfi)
  text = log.getvalue()
  assert text.find(NO_CC_NOTE + " highest resolution") > -1, text
  assert text.find(NO_CC_NOTE + " lowest resolution") > -1, text
  assert text.find(SKIP_NOTE) > -1, text
  # Worst-anisotropy location is summarized; with no directions nothing
  #  is tabulated, so no CC* heading is printed
  assert text.find(ANISO_NOTE) > -1, text
  assert not cc_heading_lines(text), text

def tst_no_cc_list():
  '''Every si has cc_list None: no CC* anywhere, but data otherwise present'''
  mmm, log = get_mmm()
  sfi = make_scale_factor_info([make_location(shift = 0., cc_list = False),
     make_location(shift = 0.2, cc_list = False)])
  mmm._summarize_scale_factor_info(sfi)
  text = log.getvalue()
  assert text.find(NO_CC_NOTE + " highest resolution") > -1, text
  assert text.find(NO_CC_NOTE + " lowest resolution") > -1, text
  assert text.find(SKIP_NOTE) > -1, text
  # Worst-anisotropy location is summarized and its CC* heading has no mean
  assert text.find(ANISO_NOTE) > -1, text
  lines = cc_heading_lines(text)
  assert lines, text
  for line in lines:
    assert line.find("(Mean CC*:") < 0, text

def overall_text(text):
  marker = "Summary of overall average scaling information"
  assert text.find(marker) > -1, text
  return text[text.find(marker):]

def scale_factor_rows(text):
  '''Table rows following the first Scale factors heading in text'''
  lines = text.splitlines()
  for i, line in enumerate(lines):
    if line.find("Scale factors by direction vector") > -1:
      # Skip the ALL header line and blank lines; keep the 4 bin rows
      rows = [l.split() for l in lines[i+2:] if l.strip()]
      return rows[:4]
  return []

def tst_mixed_locations():
  '''One location has CC* values, the other has none'''
  # Variant A: second location has cc_list None in every direction
  # Variant B: second location has an empty scaling_info_list
  for missing in (make_location(shift = 0.2, cc_list = False),
      make_location(shift = 0.2, empty = True)):
    mmm, log = get_mmm()
    sfi = make_scale_factor_info([make_location(shift = 0.), missing])
    mmm._summarize_scale_factor_info(sfi)
    text = log.getvalue()
    # Highest and lowest both use the only location with CC* values
    assert text.find("Scaling information for location with highest " +
       "resolution (Mean CC*: 0.75 )") > -1, text
    assert text.find("Scaling information for location with lowest " +
       "resolution (Mean CC*: 0.75 )") > -1, text
    assert text.find(NO_CC_NOTE) < 0, text
    assert text.find(SKIP_NOTE) < 0, text
    assert text.find(ANISO_NOTE) > -1, text
    overall = overall_text(text)
    assert cc_heading_lines(overall), text
    # The location without CC* (cc_list None, or no scaling information) is
    #  skipped, so the average is the complete location (0.38 if it were
    #  counted without contributing values)
    assert overall.find(
       "Estimated CC* by direction vector  (Mean CC*: 0.75)") > -1, text
    # Overall scale factors (ALL column) and directional values are those
    #  of the complete location (ALL would be halved if divided by the
    #  number of locations rather than the number used)
    rows = scale_factor_rows(overall)
    assert [r[1:] for r in rows] == [["1.00", "1.00", "1.00"],
       ["0.90", "0.90", "0.90"], ["0.80", "0.80", "0.80"],
       ["0.70", "0.70", "0.70"]], text

def tst_one_direction_missing():
  '''One direction at one location has cc_list None'''
  # Single location: summary completes and the CC* heading has no mean
  partial = make_location(shift = 0.2)
  partial.scaling_info_list[0].cc_list = None
  mmm, log = get_mmm()
  mmm._summarize_scale_factor_info(make_scale_factor_info([partial]))
  text = log.getvalue()
  assert mmm._get_mean_cc_star(partial) is None
  lines = cc_heading_lines(text)
  assert len(lines) == 1, text
  assert lines[0].find("(Mean CC*:") < 0, text

  # Two locations: the partial one is treated as having no CC* values
  # Make it the location with most anisotropy so it is summarized too
  partial = make_location(shift = 0.2, fo_scale = 6.)
  partial.scaling_info_list[0].cc_list = None
  mmm, log = get_mmm()
  sfi = make_scale_factor_info([make_location(shift = 0.), partial])
  mmm._summarize_scale_factor_info(sfi)
  text = log.getvalue()
  # Highest and lowest both use the complete location (0.75), not 0.55
  assert text.find("Scaling information for location with highest " +
     "resolution (Mean CC*: 0.75 )") > -1, text
  assert text.find("Scaling information for location with lowest " +
     "resolution (Mean CC*: 0.75 )") > -1, text
  assert text.find(NO_CC_NOTE) < 0, text
  assert text.find(SKIP_NOTE) < 0, text
  # Most-anisotropy summary of the partial location has no mean CC*
  aniso = text[text.find(ANISO_NOTE):text.find("Anisotropy vs position")]
  assert text.find(ANISO_NOTE) > -1, text
  lines = cc_heading_lines(aniso)
  assert len(lines) == 1 and lines[0].find("(Mean CC*:") < 0, text
  # Overall average skips the partial location (0.75; 0.65 if included)
  assert overall_text(text).find(
     "Estimated CC* by direction vector  (Mean CC*: 0.75)") > -1, text

def run(args = None):
  all_tests = [('complete', tst_complete), ('no_cc_star', tst_no_cc_star),
     ('no_cc_list', tst_no_cc_list),
     ('all_skipped', tst_all_skipped),
     ('mixed_locations', tst_mixed_locations),
     ('one_direction_missing', tst_one_direction_missing)]
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
