"""
Image reading tests against the dials_regression suite
"""
from __future__ import division, print_function

import os
import re

import pytest
import py.path

from dxtbx.format.Registry import Registry
import libtbx.load_env
from iotbx.command_line.detector_image_as_png import convert_image
# from iotbx.command_line.detector_image_as_png import run as pngrun
from rstbx.slip_viewer.slip_viewer_image_factory import SlipViewerImageFactory
import scitbx.matrix


def _dials_regression():
  """
  Duplicate of dials_regression search logic, outside of fixture.

  It is currently not possible to use a fixture value to generate a
  parametrization for tests. Using the dials_regression fixture directly
  has issues with a) safely importing and b) skipping outside of tests
  is prohibited. The downside of this is that the logic is duplicated
  """
  try:
    import dials_regression as dr
  except ImportError:
    raise RuntimeError("dials_regression required for this test")
  return os.path.dirname(dr.__file__)


def _generate_all_test_images():
  """Internal convenience function to generates a list of test images.

  Along with a custom LBL-only file, iterates over all contents of
  dials_regression/image_examples and yields data for each file, for
  converting to a fixture.

  Do not use this function directly for tests - instead use the
  test_image fixture.

  Generates (path, name) where:
    path could be a pytest-marked value - and could be a skipped dummy
    name is intended as a short test-parameter name
  """

  # Handle the special berkeley-only h5 file
  special_h5 = "/net/viper/raid1/dectris/eiger16MNov2015/2015_11_10/insu6_1_master.h5"
  try:
    import h5py
    assert h5py.version, "Supress unused import warnings"
  except ImportError:
    yield pytest.mark.skip(reason="h5py not found")(special_h5), special_h5
  else:
    yield pytest.mark.skipif(
        not os.path.isfile(special_h5),
        reason="LBL-only file not present")(special_h5), special_h5

  # Filename patterns to ignore completely
  ignore_files = [
      r"\..*",  # Hidden files
      r"\.svn",  # SVN directories
      r"ED_From_TIFF",  # removed [2017:jmp@r1618] for unknown reason
      r"LCLS_cspad_nexus/cxi78513_bslz4_r0014_subset4_\d+\.h5$",  # Non-master files
      r"LCLS_cspad_nexus/[^/]*\.h5$",  # Incompatible with current nexus definition [2018:asmit@r1902]
      r"DLS_eBIC/image_0001\.tif$",  # The format is in testing_dxtbx_format_classes/, not dxtbx, so skip the example image [2017:upintheair@r1693]
      r"README$",
      r".*\.(?:pyc?|log|json)$",  # Extensions to ignore
  ]
  # Pattern to match item at start of string or after directory indicator
  DIRSTART = r"^(?:.*/)?"
  # Handle windows paths, and ignore case when matching
  ignore_files = [
      re.compile((DIRSTART + x).replace(r"/", re.escape(os.sep)),
                 re.IGNORECASE) for x in ignore_files
  ]

  # Try to find dials_regression
  try:
    dials_regression = _dials_regression()
  except RuntimeError:
    # Have one, skipped placeholder test, if there is no dials_regression
    yield (
        pytest.mark.skip(reason="dials_regression required for image format tests")
        ("image_examples"), "image_examples")
    return

  image_dir = os.path.join(dials_regression, "image_examples")
  # Now look at everything in image_examples
  for root, dirs, files in os.walk(image_dir):
    for file in files:
      full_path = os.path.join(root, file)
      rel_path = os.path.relpath(full_path, image_dir)

      # Test the filename patterns
      if any(pattern.search(full_path) for pattern in ignore_files):
        print("Ignoring {}".format(full_path))
        continue

      # Give this file back
      yield (full_path, rel_path)


# Generate this list once to use as a fixture
_all_images = list(_generate_all_test_images())


@pytest.fixture(
    params=[x[0] for x in _all_images], ids=[x[1] for x in _all_images])
def test_image(request, dials_regression):
  """Fixture to allow tests to be parametrized for every test image"""

  # Validate that our custom dials_regression is the same as the fixture
  assert _dials_regression() == dials_regression, "dials_regression mismatch"
  return request.param


@pytest.fixture
def test_image_for_reading(test_image):
  """Test images, plus xfail marking for things we expect to not read"""

  # Explicit list of directories that we expect to fail.
  # References are to original dials_regression commit that added the exclusion
  expected_fail_dirs = {
    'ADSC_CBF': "something wrong with ADSC CBF format reader         [2013:sauter@r12]",
    'APS_19BM': "APS beamline 19; appears to be an uncorrected image [2013:sauter@r16]",
    'SACLA_MPCCD_Cheetah': "MPCCD not supported by iotbx             [2016:nakane@r1540]",
    'putative_imgCIF_HDF5_mapping': "test hdf5 image set not yet supported [2013:aaron@r66]",
  }
  # Check that the name is an actual folder in the file path
  path_parts = {x.basename for x in py.path.local(test_image).parts()}
  is_expected_fail = set(expected_fail_dirs.keys()) & path_parts
  if is_expected_fail:
    reason = expected_fail_dirs[is_expected_fail.pop()]
    pytest.xfail(reason=reason)

  return test_image


@pytest.mark.regression
def test_read_image(test_image_for_reading):
  """Test that all the regression images can be read"""
  if "LCLS" in test_image_for_reading and not libtbx.env.has_module('xfel'):
    pytest.skip("Ignoring LCLS because xfel missing")
    return

  format_instance = Registry.find(test_image_for_reading)
  instance = format_instance(test_image_for_reading)

  print("Reading", test_image_for_reading)
  print("Format:", format_instance)

  # Test metadata reading
  instance.get_goniometer()
  instance.get_beam()
  instance.get_scan()
  detector = instance.get_detector()
  # From old test_dxtbx; get the image size
  if detector is not None:
    detector[0].get_image_size()

  for panel in detector:
    d_mat = scitbx.matrix.sqr(panel.get_d_matrix())
    if d_mat.determinant() == 0:
      print("  d matrix with zero determinant")
    if d_mat.determinant() < 0:
      print("  d matrix with negative determinant")

  # Test reading of the raw data
  # XDS, HKL we expect to fail for this  - so skip this part for those
  if not test_image_for_reading[-3:].upper() in {"XDS", "HKL"}:
    try:
      R_raw_data = instance.get_raw_data()
    except TypeError:
      R_raw_data = instance.get_raw_data(0)

    if not isinstance(R_raw_data, tuple):
      R_raw_data = (R_raw_data,)

    print("%-40s" % format_instance.__name__, R_raw_data[0].focus())

    # Set instance.detectorbase, if available. There used to be a blacklist
    # of files not to call this with, but relying on the attribute test
    # seems to be just as effective
    instance.get_detectorbase()
    print("  Have detectorbase? ", hasattr(instance, "detectorbase"))

    # test the older detectorbase interface if available
    if hasattr(instance, "detectorbase"):
      imgfactory = SlipViewerImageFactory(test_image_for_reading)
      imgfactory.read()
      print("  Detectorbase:", instance.detectorbase.__class__.__name__)

      try:
        print(imgfactory.rawdata.focus())
      except AttributeError:
        # Not all instances have this attribute
        print("  multireadout")

      I_raw_data = imgfactory.get_raw_data()
      if not isinstance(I_raw_data, tuple):
        I_raw_data = (I_raw_data,)

      # NOTE dxtbx and image factory arrays are compared here for identical values.
      for Ip, Rp in zip(I_raw_data, R_raw_data):
        assert (Ip == Rp).all_eq(True)

      convert_image(test_image_for_reading, graphics_bin=2).output().getvalue()

    print

    # optional test;actually write out the images to PNG files
    # pngrun([fk,])


def test_no_multiple_format_understanding(test_image):
  """ for a given image file, walks the whole tree of Format objects.
  If the file can be understood by multiple Format objects at the same
  inheiritance level, it will return False, otherwise True."""

  Registry.setup()

  global any_understood, highest_level, found_a_repeat

  any_understood = False

  def recurse(format, image_file, level):
    global any_understood, highest_level, found_a_repeat
    for child in format._children:
      understood = child.understand(image_file)
      if understood:
        any_understood = True
        print ("level: %d" % (level), child.__name__)
        found_a_repeat = level == highest_level
        highest_level = level
        recurse(child, image_file, level + 1)

  level = 1
  highest_level = 0
  found_a_repeat = False
  for format in Registry._formats:
    understood = format.understand(test_image)
    if understood:
      any_understood = True
      print ("level: %d" % (level), format.__name__)
      found_a_repeat = level == highest_level
      highest_level = 1
      recurse(format, test_image, level + 1)

  assert not found_a_repeat, "image file understood by multiple Format objects"
  # It's a failure if nothing could understand this file
  assert any_understood, "No formatter could be found"
