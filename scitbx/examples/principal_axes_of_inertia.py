from __future__ import absolute_import, division, print_function
from scitbx.math import principal_axes_of_inertia
from scitbx.array_family import flex

def run():
  points = flex.vec3_double([
    ( 8.292,  1.817,  6.147),
    ( 9.159,  2.144,  7.299),
    (10.603,  2.331,  6.885),
    (11.041,  1.811,  5.855),
    ( 9.061,  1.065,  8.369),
    ( 7.665,  0.929,  8.902),
    ( 6.771,  0.021,  8.327),
    ( 7.210,  1.756,  9.920),
    ( 5.480, -0.094,  8.796),
    ( 5.904,  1.649, 10.416),
    ( 5.047,  0.729,  9.831),
    ( 3.766,  0.589, 10.291),
    (11.358,  2.999,  7.612)])
  pai = principal_axes_of_inertia(points=points)
  print(pai.center_of_mass())
  print(pai.inertia_tensor())
  es = pai.eigensystem()
  print(list(es.values()))
  print(list(es.vectors()))

if (__name__ == "__main__"):
  run()
