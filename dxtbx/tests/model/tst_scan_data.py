from __future__ import absolute_import, division

from dxtbx.model import Scan

def tst_is_batch_valid(scan):
  """Check that the is_angle_valid function behaves properly."""
  br1, br2 = scan.get_batch_range()
  for i in range(0, br1):
    assert(scan.is_batch_valid(i) == False)
  for i in range(br1, br2+1):
    assert(scan.is_batch_valid(i) == True)
  for i in range(br2+1, br2+100):
    assert(scan.is_batch_valid(i) == False)
  print "OK"

def tst_is_angle_valid(scan):
  """Check that the is_angle_valid function behaves properly."""
  oscillation_range = scan.get_oscillation_range()
  os1 = int(oscillation_range[0])
  os2 = int(oscillation_range[1])
  for i in range(0, os1):
    assert(scan.is_angle_valid(i) == False)
  for i in range(os1, os2):
    assert(scan.is_angle_valid(i) == True)
  for i in range(os2+1, 360):
    assert(scan.is_angle_valid(i) == False)
  print "OK"

def tst_is_frame_valid(scan):
  """Check that the is_frame_valid function behaves properly."""
  image_range = scan.get_image_range()
  for i  in range(image_range[0] - 100, image_range[0]):
    assert(scan.is_image_index_valid(i) == False)
  for i in range(image_range[0], image_range[1]+1):
    assert(scan.is_image_index_valid(i) == True)
  for i in range(image_range[1]+1, image_range[1] + 100):
    assert(scan.is_image_index_valid(i) == False)
  print "OK"

def tst_get_angle_from_frame(scan):
  pass

def tst_get_frame_from_angle(scan):
  pass

def tst_get_frames_with_angle(scan):
  pass

def tst_scan_oscillation_recycle(scan):
  for deg in (True, False):
    oscillation = scan.get_oscillation(deg=deg)
    scan.set_oscillation(oscillation, deg=deg)
    assert scan.get_oscillation(deg=deg) == oscillation
  print 'OK'

def tst_scan_360_append():

  scan1 = Scan((1, 360), (0.0, 1.0))
  scan2 = Scan((361, 720), (0.0, 1.0))

  scan = scan1 + scan2
  eps = 1e-7
  assert(scan.get_num_images() == 720)
  assert(abs(scan.get_oscillation()[0] - 0.0) < eps)
  assert(abs(scan.get_oscillation()[1] - 1.0) < eps)
  assert(scan.get_image_range() == (1, 720))
  assert(scan.get_batch_range() == (1, 720))
  print 'OK'

  scan1 = Scan((1, 360), (0.0, 1.0))
  scan2 = Scan((361, 720), (360.0, 1.0))

  scan = scan1 + scan2
  eps = 1e-7
  assert(scan.get_num_images() == 720)
  assert(abs(scan.get_oscillation()[0] - 0.0) < eps)
  assert(abs(scan.get_oscillation()[1] - 1.0) < eps)
  assert(scan.get_image_range() == (1, 720))
  assert(scan.get_batch_range() == (1, 720))

  from libtbx.test_utils import Exception_expected
  scan2.set_batch_offset(10)
  try: scan = scan1 + scan2
  except Exception: pass
  else: raise Exception_expected

  print 'OK'

def tst_swap():

  scan1 = Scan((1, 20), (0.0, 1.0))
  scan2 = Scan((40, 60), (10.0, 2.0))
  scan1.swap(scan2)
  assert(scan2.get_image_range() == (1, 20))
  assert(scan1.get_image_range() == (40, 60))
  assert(scan2.get_oscillation() == (0.0, 1.0))
  assert(scan1.get_oscillation() == (10.0, 2.0))
  print 'OK'

def tst_from_phil():
  from dxtbx.model.scan import ScanFactory, scan_phil_scope
  from libtbx.phil import parse

  params = scan_phil_scope.fetch(parse('''
    scan {
      image_range = 1, 10
      oscillation = (-4, 0.1)
    }
  ''')).extract()

  s1 = ScanFactory.from_phil(params)

  assert s1.get_num_images() == 10
  assert s1.get_image_range() == (1, 10)
  assert s1.get_oscillation() == (-4, 0.1)
  assert s1.get_batch_offset() == 0
  assert s1.get_batch_range() == s1.get_image_range()
  for i in range(s1.get_image_range()[0], s1.get_image_range()[1]+1):
    assert s1.get_batch_for_image_index(i) == i
    assert s1.get_batch_for_array_index(i-1) == i

  params = scan_phil_scope.fetch(parse('''
    scan {
      image_range = 1, 20
      extrapolate_scan = True
      oscillation = (20, 0.01)
      batch_offset = 10
    }
  ''')).extract()

  s2 = ScanFactory.from_phil(params, s1)
  assert s2.get_num_images() == 20
  assert s2.get_image_range() == (1, 20)
  assert s2.get_oscillation() == (20, 0.01)
  assert s2.get_batch_offset() == 10
  assert s2.get_batch_range() == (11, 30)
  ir1, ir2 = s2.get_image_range()
  for i in range(ir1, ir2):
    assert s2.get_batch_for_image_index(i) == i + s2.get_batch_offset()
    assert s2.is_batch_valid(s2.get_batch_for_image_index(i))

  print 'OK'

def run():
  image_range = (0, 1000)
  oscillation = (0, 0.1)
  scan = Scan(image_range, oscillation)
  scan.set_batch_offset(100)
  tst_scan_oscillation_recycle(scan)
  tst_is_angle_valid(scan)
  tst_is_batch_valid(scan)
  tst_is_frame_valid(scan)
  tst_get_angle_from_frame(scan)
  tst_get_frame_from_angle(scan)
  tst_get_frames_with_angle(scan)
  tst_scan_360_append()
  tst_swap()
  tst_from_phil()

if __name__ == '__main__':
  run()
