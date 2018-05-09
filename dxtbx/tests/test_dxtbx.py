from __future__ import absolute_import, division, print_function

import os
import pytest

def pytest_generate_tests(metafunc):
  if 'dxtbx_test_image' not in metafunc.fixturenames:
    return
  try:
    from dials_regression.image_examples.get_all_working_images import get_all_working_images
  except ImportError:
    metafunc.parametrize('dxtbx_test_image',
        [pytest.param(None, marks=pytest.mark.skip('dials_regression is required for this group of tests'))],
        ids=["noregression"])
    return
  imagelist = get_all_working_images()
  metafunc.parametrize('dxtbx_test_image', imagelist, ids=[i[1] for i in imagelist])

def test_dxtbx(dials_regression, dxtbx_test_image):
  #from boost.python import streambuf
  #from dxtbx import read_uint16
  from dxtbx.format.Registry import Registry

  directory, image = dxtbx_test_image
  file_path = os.path.join(dials_regression, 'image_examples', directory, image)

  format = Registry.find(file_path)
  i = format(file_path)
  det = i.get_detector()
  if det is not None:
    size = det[0].get_image_size()
  b = i.get_beam()
  g = i.get_goniometer()
  s = i.get_scan()
  try:
    d = i.get_raw_data()
  except IOError:
    pass

@pytest.mark.skip("Test disabled in https://github.com/cctbx/cctbx_project/commit/e6ae0e03afe87500ad573b2fde42dca9334ea2e5")
def test_dxtbx_compressed(dials_regression, tmpdir):
  from dxtbx.format.Registry import Registry

  from dials_regression.image_examples.get_all_working_images import \
      get_all_working_images

  # test that reading gz or bz2 compressed files works: it doesn't!

  from libtbx import smart_open
  import shutil

  tmpdir.chdir()

  for directory, image in get_all_working_images():
    file_path = os.path.join(dials_regression, 'image_examples',
                             directory, image)
    for ext in ('.gz', '.bz2')[:]:
      compressed_path = os.path.basename(file_path) + ext
      with open(file_path, 'rb') as f_in, smart_open.for_writing(compressed_path) as f_out:
        shutil.copyfileobj(f_in, f_out)
      print(file_path, compressed_path)
      format = Registry.find(compressed_path)
      try:
        i = format(compressed_path)
      except Exception:
        print('Error reading compressed file: %s' % compressed_path)
        import traceback
        traceback.print_exc()
      else:
        print('Successfully read compressed file: %s' % compressed_path)
        det = i.get_detector()
        if det is not None:
          size = det[0].get_image_size()
        b = i.get_beam()
        g = i.get_goniometer()
        s = i.get_scan()
        try:
          d = i.get_raw_data()
        except IOError:
          pass

def test_sweep():
  from dxtbx.sweep_filenames import template_regex

  questions_answers = {
      'foo_bar_001.img':'foo_bar_###.img',
      'foo_bar001.img':'foo_bar###.img',
      'foo_bar_1.8A_001.img':'foo_bar_1.8A_###.img',
      'foo_bar.001':'foo_bar.###',
      'foo_bar_001.img1000':'foo_bar_###.img1000',
      'foo_bar_00001.img':'foo_bar_#####.img'
      }

  for filename in questions_answers:
    answer = template_regex(filename)
    assert answer[0] == questions_answers[filename]
