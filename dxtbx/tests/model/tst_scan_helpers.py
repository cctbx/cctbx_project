from __future__ import absolute_import, division
from dxtbx.model import is_angle_in_range
from dxtbx.model import get_range_of_mod2pi_angles
from dxtbx.model import get_mod2pi_angles_in_range

def tst_is_angle_in_range():
  """Test that for a range of angles and angular ranges, the
  is_angle_in_range function correctly calculates if the angle
  is in the range.

  """
  from random import random

  # Some helper lambda functions
  mod_360 = lambda x: x - 360.0 * floor(x / 360.0)
  mod_360_range = lambda x: (mod_360(x[0]), mod_360(x[1]))
  random_0_360 = lambda: int(random() * 360.0)

  # Create a number of ranges between 0 and 360 and check 360 degrees worth
  # of angles to see if they are within the range
  num_range = 100
  for n in range(num_range):
    angular_range = (random_0_360(), random_0_360())

    # If A < B or A > B
    if angular_range[0] < angular_range[1]:

      # Check that the following are true
      #   angle in range 0 -> A = False
      #   angle in range A -> B = True
      #   angle in range B -> 360 = False
      for angle in range(0, angular_range[0]):
        assert(is_angle_in_range(angular_range, angle, True) == False)
      for angle in range(angular_range[0], angular_range[1]+1):
        assert(is_angle_in_range(angular_range, angle, True) == True)
      for angle in range(angular_range[1]+1, 360):
        assert(is_angle_in_range(angular_range, angle, True) == False)
    else:

      # Check that the following are true
      #   angle in range 0 -> B = True
      #   angle in range B -> A = False
      #   angle in range A -> 360 = True
      for angle in range(0, angular_range[1]+1):
        assert(is_angle_in_range(angular_range, angle, True) == True)
      for angle in range(angular_range[1]+1, angular_range[0]):
        assert(is_angle_in_range(angular_range, angle, True) == False)
      for angle in range(angular_range[0], 360):
        assert(is_angle_in_range(angular_range, angle, True) == True)

  # Create a range over 360 and make sure all angles are valid
  angular_range = (-10, 370)
  for angle in range(0, 360):
    assert(is_angle_in_range(angular_range, angle, True) == True)

def tst_get_range_of_mod2pi_angles():
  """Get the range of equivalent within a given angular range.

  """
  # In a 360 deg range, have only 1 angle
  a0, a1 = get_range_of_mod2pi_angles((0, 360), 180, deg=True)
  assert(a0 == 180 and a1 == 180)

  # If no angles within range, a0 > a1
  a0, a1 = get_range_of_mod2pi_angles((181, 360), 180, deg=True)
  assert(a0 > a1)

  # With 720 deg range, have 2 angles
  a0, a1 = get_range_of_mod2pi_angles((0, 720), 180, deg=True)
  assert(a0 == 180 and a1 == 540)

  # With 1080 deg range, have 3 angles
  a0, a1 = get_range_of_mod2pi_angles((-360, 720), 180, deg=True)
  assert(a0 == -180 and a1 == 540)

def tst_get_mod2pi_angles_in_range():
  """Get the list of equivalent within a given angular range.

  """
  # In a 360 deg range, have only 1 angle
  a = get_mod2pi_angles_in_range((0, 360), 180, deg=True)
  assert(len(a) == 1 and a[0] == 180)

  # If no angles within range have no angles
  a = get_mod2pi_angles_in_range((181, 360), 180, deg=True)
  assert(len(a) == 0)

  # With 720 deg range, have 2 angles
  a = get_mod2pi_angles_in_range((0, 720), 180, deg=True)
  assert(len(a) == 2 and a[0] == 180 and a[1] == 540)

  # With 1080 deg range, have 3 angles
  a = get_mod2pi_angles_in_range((-360, 720), 180, deg=True)
  assert(len(a) == 3 and a[0] == -180 and a[1] == 180 and a[2] == 540)

def run():
  """Run tests for the scan_helpers.h module."""
  tst_is_angle_in_range()
  tst_get_range_of_mod2pi_angles()
  tst_get_mod2pi_angles_in_range()

if __name__ == '__main__':
  run()
