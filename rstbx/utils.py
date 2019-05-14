from __future__ import absolute_import, division, print_function

# XXX please keep all shared-library imports inline
import math

def get_scattering_angle(x,
                          y,
                          center_x,
                          center_y,
                          distance,
                          detector_two_theta, # in radians, please
                          distance_is_corrected=False):
  if (not distance_is_corrected) and (detector_two_theta != 0):
    distance = distance / math.cos(detector_two_theta)
  r_x = center_x - x
  r_y = center_y - y
  r = math.sqrt(r_x**2 + r_y**2)
  if (detector_two_theta == 0):
    return math.atan(r / distance)
  else :
    if (r_y <= 0):
      gamma = (math.pi / 2) - detector_two_theta
    else :
      gamma = (math.pi / 2) + detector_two_theta
    # distance from crystal to (0,y)
    d2 = math.sqrt(distance**2 + r_y**2 - (2*distance*abs(r_y)*math.cos(gamma)))
    # distance from crystal to (x,y)
    d3 = math.sqrt(d2**2 + r_x**2)
    angle = math.acos((distance**2 + d3**2 - r**2) / (2*distance*d3))
    return angle

def angle_between_points(x1,
                          y1,
                          x2,
                          y2,
                          center_x,
                          center_y,
                          distance,
                          detector_two_theta,
                          distance_is_corrected=False):
  if (not distance_is_corrected) and (detector_two_theta != 0):
    distance = distance / math.cos(detector_two_theta)
  r_x_1 = center_x - x1
  r_x_2 = center_x - x2
  r_y_1 = center_y - y1
  r_y_2 = center_y - y2
  r = math.sqrt((x2-x1)**2 + (y2-y1)**2)
  if (r_y_1 <= 0):
    gamma_1 = (math.pi / 2) - detector_two_theta
  else :
    gamma_1 = (math.pi / 2) + detector_two_theta
  # distance from crystal to (0,y1)
  d2_1 = math.sqrt(distance**2 + r_y_1**2 - (2*distance*abs(r_y_1) *
    math.cos(gamma_1)))
  # distance from crystal to (x1,y1)
  d3_1 = math.sqrt(d2_1**2 + r_x_1**2)
  if (r_y_2 <= 0):
    gamma_2 = (math.pi / 2) - detector_two_theta
  else :
    gamma_2 = (math.pi / 2) + detector_two_theta
  # distance from crystal to (0,y2)
  d2_2 = math.sqrt(distance**2 + r_y_2**2 - (2*distance*abs(r_y_2) *
    math.cos(gamma_2)))
  # distance from crystal to (x2,y2)
  d3_2 = math.sqrt(d2_2**2 + r_x_2**2)
  angle = math.acos((d3_1**2 + d3_2**2 - r**2) / (2 * d3_1 * d3_2))
  return angle

def reciprocal_space_distance(x1,
                               y1,
                               x2,
                               y2,
                               wavelength,
                               **kwds):
  assert (wavelength > 0)
  angle = angle_between_points(x1, y1, x2, y2, **kwds)
  ewald_radius = 1 / wavelength
  return math.sqrt(2*(ewald_radius**2) - 2*(ewald_radius**2)*math.cos(angle))
