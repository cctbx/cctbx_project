from libtbx.test_utils import show_diff
from libtbx import easy_run
import sys

def run_and_compare(commands, expected_output):
  for command in commands:
    lines = easy_run.fully_buffered(
      command=command).raise_if_errors().stdout_lines
    assert not show_diff("\n".join(lines)+"\n", expected_output)

def exercise(args):
  assert len(args) == 0
  run_and_compare(
    ['iotbx.lattice_symmetry'
       ' --unit_cell="12.7923,12.7923,29.4356,102.846,102.846,22.7475"'],
    """\

Input
=====

Unit cell: (12.7923, 12.7923, 29.4356, 102.846, 102.846, 22.7475)
Space group: P 1 (No. 1)

Angular tolerance: 3.000 degrees

Similar symmetries
==================

Symmetry in minimum-lengths cell: C 1 2/m 1 (x+y,z,2*x-z) (No. 12)
      Input minimum-lengths cell: (5.04549, 12.7923, 29.3711, 100.302, 94.9273, 101.374)
           Symmetry-adapted cell: (5.04549, 12.7923, 29.3711, 100.302, 94.9273, 101.374)
            Conventional setting: C 1 2/m 1 (No. 12)
                       Unit cell: (25.0822, 5.04549, 29.4356, 90, 103.108, 90)
                 Change of basis: 1/2*x+1/2*y,-1/2*x+1/2*y,z
                         Inverse: x-y,x+y,z
      Maximal angular difference: 0.000 degrees

Symmetry in minimum-lengths cell: P -1 (No. 2)
      Input minimum-lengths cell: (5.04549, 12.7923, 29.3711, 100.302, 94.9273, 101.374)
           Symmetry-adapted cell: (5.04549, 12.7923, 29.3711, 100.302, 94.9273, 101.374)
            Conventional setting: P -1 (No. 2)
                       Unit cell: (5.04549, 12.7923, 29.3711, 100.302, 94.9273, 101.374)
                 Change of basis: y-z,x+y-z,-z
                         Inverse: -x+y,x-z,-z
      Maximal angular difference: 0.000 degrees

""")
  #
  run_and_compare(
    ['iotbx.lattice_symmetry --unit_cell=12,12,12.1,89,90,92 F',
     'cctbx.lattice_symmetry'],
    """\

Input
=====

Unit cell: (12, 12, 12.1, 89, 90, 92)
Space group: P 1 (-a+b+c,a-b+c,a+b-c) (No. 1)

Angular tolerance: 3.000 degrees

Similar symmetries
==================

Symmetry in minimum-lengths cell: F m -3 m (-x+y+z,x-y+z,x+y-z) (No. 225)
      Input minimum-lengths cell: (8.3359, 8.44603, 8.52071, 61.1613, 61.2992, 61.0215)
           Symmetry-adapted cell: (8.50892, 8.50892, 8.50892, 60, 60, 60)
            Conventional setting: F m -3 m (No. 225)
                       Unit cell: (12.0334, 12.0334, 12.0334, 90, 90, 90)
                 Change of basis: z,-x,-y
                         Inverse: -y,-z,x
      Maximal angular difference: 2.236 degrees

Symmetry in minimum-lengths cell: I 4/m m m (y-z,-x+z,x+z) (No. 139)
      Input minimum-lengths cell: (8.3359, 8.44603, 8.52071, 61.1613, 61.2992, 61.0215)
           Symmetry-adapted cell: (8.48528, 8.52071, 8.52071, 59.7251, 60.1374, 60.1374)
            Conventional setting: I 4/m m m (No. 139)
                       Unit cell: (8.48528, 8.48528, 12.1, 90, 90, 90)
                 Change of basis: -x+y,-x-y,z
                         Inverse: -1/2*x-1/2*y,1/2*x-1/2*y,z
      Maximal angular difference: 2.236 degrees

Symmetry in minimum-lengths cell: I 4/m m m (-y+z,x-z,y+z) (No. 139)
      Input minimum-lengths cell: (8.3359, 8.44603, 8.52071, 61.1613, 61.2992, 61.0215)
           Symmetry-adapted cell: (8.50301, 8.52071, 8.50301, 59.9311, 60.1377, 59.9311)
            Conventional setting: I 4/m m m (No. 139)
                       Unit cell: (8.52071, 8.52071, 12, 90, 90, 90)
                 Change of basis: -y+z,y+z,-x
                         Inverse: -z,-1/2*x+1/2*y,1/2*x+1/2*y
      Maximal angular difference: 2.236 degrees

Symmetry in minimum-lengths cell: I 4/m m m (y+z,-y+z,x-z) (No. 139)
      Input minimum-lengths cell: (8.3359, 8.44603, 8.52071, 61.1613, 61.2992, 61.0215)
           Symmetry-adapted cell: (8.50301, 8.50301, 8.52071, 59.9311, 59.9311, 60.1377)
            Conventional setting: I 4/m m m (No. 139)
                       Unit cell: (8.52071, 8.52071, 12, 90, 90, 90)
                 Change of basis: -x+z,-x-z,-y
                         Inverse: -1/2*x-1/2*y,-z,1/2*x-1/2*y
      Maximal angular difference: 2.236 degrees

Symmetry in minimum-lengths cell: R -3 m :R (No. 166)
      Input minimum-lengths cell: (8.3359, 8.44603, 8.52071, 61.1613, 61.2992, 61.0215)
           Symmetry-adapted cell: (8.43456, 8.43456, 8.43456, 61.1649, 61.1649, 61.1649)
            Conventional setting: R -3 m :H (No. 166)
                       Unit cell: (8.58263, 8.58263, 20.4766, 90, 90, 120)
                 Change of basis: 4/3*x-2/3*y+2/3*z,2/3*x+2/3*y+4/3*z,-1/3*x-1/3*y+1/3*z
                         Inverse: 1/2*x-z,-1/2*x+1/2*y-z,1/2*y+z
      Maximal angular difference: 1.486 degrees

Symmetry in minimum-lengths cell: R -3 m :H (x+z,-y+z,-3*z) (No. 166)
      Input minimum-lengths cell: (8.3359, 8.44603, 8.52071, 61.1613, 61.2992, 61.0215)
           Symmetry-adapted cell: (8.43456, 8.43456, 8.58263, 60.5691, 60.5691, 60)
            Conventional setting: R -3 m :H (No. 166)
                       Unit cell: (8.43456, 8.43456, 21.2021, 90, 90, 120)
                 Change of basis: -4/3*x-2/3*y-2/3*z,-2/3*x+2/3*y-4/3*z,1/3*x-1/3*y-1/3*z
                         Inverse: -1/2*x+z,-1/2*x+1/2*y-z,-1/2*y-z
      Maximal angular difference: 1.510 degrees

Symmetry in minimum-lengths cell: R -3 m :H (-y+z,-3*z,x+z) (No. 166)
      Input minimum-lengths cell: (8.3359, 8.44603, 8.52071, 61.1613, 61.2992, 61.0215)
           Symmetry-adapted cell: (8.48448, 8.53328, 8.48448, 60.189, 60, 60.189)
            Conventional setting: R -3 m :H (No. 166)
                       Unit cell: (8.48448, 8.48448, 20.9617, 90, 90, 120)
                 Change of basis: -2/3*x+2/3*y+4/3*z,2/3*x+4/3*y+2/3*z,-1/3*x+1/3*y-1/3*z
                         Inverse: -1/2*x+1/2*y-z,1/2*y+z,1/2*x-z
      Maximal angular difference: 2.177 degrees

Symmetry in minimum-lengths cell: R -3 m :H (-3*z,x+z,-y+z) (No. 166)
      Input minimum-lengths cell: (8.3359, 8.44603, 8.52071, 61.1613, 61.2992, 61.0215)
           Symmetry-adapted cell: (8.48448, 8.53328, 8.53328, 60, 59.8096, 59.8096)
            Conventional setting: R -3 m :H (No. 166)
                       Unit cell: (8.53328, 8.53328, 20.7226, 90, 90, 120)
                 Change of basis: 2/3*x-4/3*y+2/3*z,4/3*x-2/3*y-2/3*z,1/3*x+1/3*y+1/3*z
                         Inverse: 1/2*y+z,-1/2*x+z,1/2*x-1/2*y+z
      Maximal angular difference: 2.177 degrees

Symmetry in minimum-lengths cell: I m m m (y-z,-x+z,x+z) (No. 71)
      Input minimum-lengths cell: (8.3359, 8.44603, 8.52071, 61.1613, 61.2992, 61.0215)
           Symmetry-adapted cell: (8.3359, 8.52071, 8.52071, 60.8666, 60.7149, 60.7149)
            Conventional setting: I m m m (No. 71)
                       Unit cell: (8.3359, 8.63208, 12.1, 90, 90, 90)
                 Change of basis: x+y,-x+y,z
                         Inverse: 1/2*x-1/2*y,1/2*x+1/2*y,z
      Maximal angular difference: 1.001 degrees

Symmetry in minimum-lengths cell: I m m m (-y+z,x-z,y+z) (No. 71)
      Input minimum-lengths cell: (8.3359, 8.44603, 8.52071, 61.1613, 61.2992, 61.0215)
           Symmetry-adapted cell: (8.50301, 8.44603, 8.50301, 60.2214, 60.715, 60.2214)
            Conventional setting: I m m m (No. 71)
                       Unit cell: (8.44603, 8.59474, 12, 90, 90, 90)
                 Change of basis: -y+z,y+z,-x
                         Inverse: -z,-1/2*x+1/2*y,1/2*x+1/2*y
      Maximal angular difference: 2.000 degrees

Symmetry in minimum-lengths cell: I m m m (y+z,-y+z,x-z) (No. 71)
      Input minimum-lengths cell: (8.3359, 8.44603, 8.52071, 61.1613, 61.2992, 61.0215)
           Symmetry-adapted cell: (8.50301, 8.50301, 8.52071, 59.9311, 59.9311, 60.1377)
            Conventional setting: I m m m (No. 71)
                       Unit cell: (8.52071, 8.52071, 12, 90, 90, 90)
                 Change of basis: -x+z,-x-z,-y
                         Inverse: -1/2*x-1/2*y,-z,1/2*x-1/2*y
      Maximal angular difference: 2.236 degrees

Symmetry in minimum-lengths cell: F m m m (-x+y+z,x-y+z,x+y-z) (No. 69)
      Input minimum-lengths cell: (8.3359, 8.44603, 8.52071, 61.1613, 61.2992, 61.0215)
           Symmetry-adapted cell: (8.48528, 8.52071, 8.52071, 59.7251, 60.1374, 60.1374)
            Conventional setting: F m m m (No. 69)
                       Unit cell: (12, 12, 12.1, 90, 90, 90)
                 Change of basis: -x,-y,z
                         Inverse: -x,-y,z
      Maximal angular difference: 2.236 degrees

Symmetry in minimum-lengths cell: C 1 2/m 1 (x+y,-x+y,x-y+z) (No. 12)
      Input minimum-lengths cell: (8.3359, 8.44603, 8.52071, 61.1613, 61.2992, 61.0215)
           Symmetry-adapted cell: (8.48528, 8.52071, 8.52071, 59.7251, 60.1374, 60.1374)
            Conventional setting: I 1 2/m 1 (No. 12)
                       Unit cell: (8.52071, 12, 8.52071, 90, 90.4755, 90)
                 Change of basis: x+z,-y,x-z
                         Inverse: 1/2*x+1/2*z,-y,1/2*x-1/2*z
      Maximal angular difference: 2.236 degrees

Symmetry in minimum-lengths cell: C 1 2/m 1 (x-y+z,x+y,-x+y) (No. 12)
      Input minimum-lengths cell: (8.3359, 8.44603, 8.52071, 61.1613, 61.2992, 61.0215)
           Symmetry-adapted cell: (8.3359, 8.52071, 8.52071, 60.8666, 60.7149, 60.7149)
            Conventional setting: I 1 2/m 1 (No. 12)
                       Unit cell: (8.63208, 12.1, 8.3359, 90, 90, 90)
                 Change of basis: -x+y,-z,-x-y
                         Inverse: -1/2*x-1/2*z,1/2*x-1/2*z,-y
      Maximal angular difference: 1.001 degrees

Symmetry in minimum-lengths cell: C 1 2/m 1 (-x+y,x-y+z,x+y) (No. 12)
      Input minimum-lengths cell: (8.3359, 8.44603, 8.52071, 61.1613, 61.2992, 61.0215)
           Symmetry-adapted cell: (8.48528, 8.44603, 8.52071, 60.013, 60.7147, 60.4302)
            Conventional setting: I 1 2/m 1 (No. 12)
                       Unit cell: (8.44603, 12, 8.59474, 90, 90.4756, 90)
                 Change of basis: y-z,-x,y+z
                         Inverse: -y,1/2*x+1/2*z,-1/2*x+1/2*z
      Maximal angular difference: 2.000 degrees

Symmetry in minimum-lengths cell: C 1 2/m 1 (z,x-y,x+y) (No. 12)
      Input minimum-lengths cell: (8.3359, 8.44603, 8.52071, 61.1613, 61.2992, 61.0215)
           Symmetry-adapted cell: (8.3359, 8.48345, 8.48345, 61.1625, 61.1613, 61.1613)
            Conventional setting: I 1 2/m 1 (No. 12)
                       Unit cell: (8.3359, 8.63208, 12.1, 90, 90.7198, 90)
                 Change of basis: x+y,x-y,-z
                         Inverse: 1/2*x+1/2*y,1/2*x-1/2*y,-z
      Maximal angular difference: 0.695 degrees

Symmetry in minimum-lengths cell: C 1 2/m 1 (x+y,z,x-y) (No. 12)
      Input minimum-lengths cell: (8.3359, 8.44603, 8.52071, 61.1613, 61.2992, 61.0215)
           Symmetry-adapted cell: (8.42881, 8.44603, 8.42881, 61.0941, 61.3067, 61.0941)
            Conventional setting: I 1 2/m 1 (No. 12)
                       Unit cell: (8.44603, 8.59474, 12, 90, 91.4206, 90)
                 Change of basis: y-z,y+z,x
                         Inverse: z,1/2*x+1/2*y,-1/2*x+1/2*y
      Maximal angular difference: 1.486 degrees

Symmetry in minimum-lengths cell: C 1 2/m 1 (x-y,x+y,z) (No. 12)
      Input minimum-lengths cell: (8.3359, 8.44603, 8.52071, 61.1613, 61.2992, 61.0215)
           Symmetry-adapted cell: (8.39115, 8.39115, 8.52071, 61.2305, 61.2305, 61.0242)
            Conventional setting: I 1 2/m 1 (No. 12)
                       Unit cell: (8.52071, 8.52071, 12, 90, 92.1185, 90)
                 Change of basis: x-z,-x-z,y
                         Inverse: 1/2*x-1/2*y,z,-1/2*x-1/2*y
      Maximal angular difference: 0.860 degrees

Symmetry in minimum-lengths cell: C 1 2/m 1 (-z,-2*x+z,x+y) (No. 12)
      Input minimum-lengths cell: (8.3359, 8.44603, 8.52071, 61.1613, 61.2992, 61.0215)
           Symmetry-adapted cell: (8.46631, 8.53956, 8.52071, 60.073, 59.7872, 60.1365)
            Conventional setting: I 1 2/m 1 (No. 12)
                       Unit cell: (8.52071, 8.52071, 12, 90, 90.6981, 90)
                 Change of basis: x+z,x-z,y
                         Inverse: 1/2*x+1/2*y,z,1/2*x-1/2*y
      Maximal angular difference: 2.177 degrees

Symmetry in minimum-lengths cell: C 1 2/m 1 (-z,-x+y,2*x+z) (No. 12)
      Input minimum-lengths cell: (8.3359, 8.44603, 8.52071, 61.1613, 61.2992, 61.0215)
           Symmetry-adapted cell: (8.42881, 8.44603, 8.57657, 60.5022, 60.7101, 59.9324)
            Conventional setting: I 1 2/m 1 (No. 12)
                       Unit cell: (8.59474, 8.44603, 12, 90, 91.3961, 90)
                 Change of basis: -y-z,-y+z,-x
                         Inverse: -z,-1/2*x-1/2*y,-1/2*x+1/2*y
      Maximal angular difference: 1.510 degrees

Symmetry in minimum-lengths cell: C 1 2/m 1 (-x+y,z,2*x-z) (No. 12)
      Input minimum-lengths cell: (8.3359, 8.44603, 8.52071, 61.1613, 61.2992, 61.0215)
           Symmetry-adapted cell: (8.3359, 8.48345, 8.5578, 60.8654, 60.8541, 60.5737)
            Conventional setting: I 1 2/m 1 (No. 12)
                       Unit cell: (8.63208, 8.3359, 12.1, 90, 90.6951, 90)
                 Change of basis: x-y,x+y,z
                         Inverse: 1/2*x+1/2*y,-1/2*x+1/2*y,z
      Maximal angular difference: 0.720 degrees

Symmetry in minimum-lengths cell: P -1 (No. 2)
      Input minimum-lengths cell: (8.3359, 8.44603, 8.52071, 61.1613, 61.2992, 61.0215)
           Symmetry-adapted cell: (8.3359, 8.44603, 8.52071, 61.1613, 61.2992, 61.0215)
            Conventional setting: P -1 (No. 2)
                       Unit cell: (8.3359, 8.44603, 8.52071, 61.1613, 61.2992, 61.0215)
                 Change of basis: -x-y-z,x-y+z,-x+y+z
                         Inverse: -1/2*x-1/2*z,-1/2*x-1/2*y,1/2*y+1/2*z
      Maximal angular difference: 0.000 degrees

""")
  #
  run_and_compare(
    ['iotbx.lattice_symmetry --unit-cel="22.54 22.54 6.35 90 90 90" I',
     'cctbx.lattice_symmetry 1'],
    """\

Input
=====

Unit cell: (22.54, 22.54, 6.35, 90, 90, 90)
Space group: P 1 (b+c,a+c,a+b) (No. 1)

Angular tolerance: 3.000 degrees

Similar symmetries
==================

Symmetry in minimum-lengths cell: I 4/m m m (-y+z,x+y,-x+y) (No. 139)
      Input minimum-lengths cell: (6.35, 16.2514, 16.2514, 87.8126, 78.7338, 78.7338)
           Symmetry-adapted cell: (6.35, 16.2514, 16.2514, 87.8126, 78.7338, 78.7338)
            Conventional setting: I 4/m m m (No. 139)
                       Unit cell: (22.54, 22.54, 6.35, 90, 90, 90)
                 Change of basis: x,y,z
                         Inverse: x,y,z
      Maximal angular difference: 0.000 degrees

Symmetry in minimum-lengths cell: I m m m (y-z,-x+z,x+z) (No. 71)
      Input minimum-lengths cell: (6.35, 16.2514, 16.2514, 87.8126, 78.7338, 78.7338)
           Symmetry-adapted cell: (6.35, 16.2514, 16.2514, 87.8126, 78.7338, 78.7338)
            Conventional setting: I m m m (No. 71)
                       Unit cell: (6.35, 22.54, 22.54, 90, 90, 90)
                 Change of basis: -z,-y,-x
                         Inverse: -z,-y,-x
      Maximal angular difference: 0.000 degrees

Symmetry in minimum-lengths cell: F m m m (x-y+z,-2*z,2*y) (No. 69)
      Input minimum-lengths cell: (6.35, 16.2514, 16.2514, 87.8126, 78.7338, 78.7338)
           Symmetry-adapted cell: (6.35, 16.2514, 16.2514, 87.8126, 78.7338, 78.7338)
            Conventional setting: F m m m (No. 69)
                       Unit cell: (6.35, 31.8764, 31.8764, 90, 90, 90)
                 Change of basis: z,-1/2*x-1/2*y,1/2*x-1/2*y
                         Inverse: -y+z,-y-z,x
      Maximal angular difference: 0.000 degrees

Symmetry in minimum-lengths cell: C 1 2/m 1 (-x+y,z,2*x-z) (No. 12)
      Input minimum-lengths cell: (6.35, 16.2514, 16.2514, 87.8126, 78.7338, 78.7338)
           Symmetry-adapted cell: (6.35, 16.2514, 16.2514, 87.8126, 78.7338, 78.7338)
            Conventional setting: I 1 2/m 1 (No. 12)
                       Unit cell: (22.54, 6.35, 22.54, 90, 90, 90)
                 Change of basis: y,-z,-x
                         Inverse: -z,x,-y
      Maximal angular difference: 0.000 degrees

Symmetry in minimum-lengths cell: C 1 2/m 1 (x-y,2*y,z) (No. 12)
      Input minimum-lengths cell: (6.35, 16.2514, 16.2514, 87.8126, 78.7338, 78.7338)
           Symmetry-adapted cell: (6.35, 16.2514, 16.2514, 87.8126, 78.7338, 78.7338)
            Conventional setting: C 1 2/m 1 (No. 12)
                       Unit cell: (6.35, 31.8764, 16.2514, 90, 101.266, 90)
                 Change of basis: 1/2*x+1/2*y-z,-1/2*x+1/2*y,x+y
                         Inverse: -y+1/2*z,y+1/2*z,-x+1/2*z
      Maximal angular difference: 0.000 degrees

Symmetry in minimum-lengths cell: C 1 2/m 1 (x+y,z,-2*y) (No. 12)
      Input minimum-lengths cell: (6.35, 16.2514, 16.2514, 87.8126, 78.7338, 78.7338)
           Symmetry-adapted cell: (6.35, 16.2514, 16.2514, 87.8126, 78.7338, 78.7338)
            Conventional setting: C 1 2/m 1 (No. 12)
                       Unit cell: (6.35, 31.8764, 16.2514, 90, 101.266, 90)
                 Change of basis: 1/2*x-1/2*y-z,1/2*x+1/2*y,x-y
                         Inverse: y+1/2*z,y-1/2*z,-x+1/2*z
      Maximal angular difference: 0.000 degrees

Symmetry in minimum-lengths cell: C 1 2/m 1 (z,x-y,x+y) (No. 12)
      Input minimum-lengths cell: (6.35, 16.2514, 16.2514, 87.8126, 78.7338, 78.7338)
           Symmetry-adapted cell: (6.35, 16.2514, 16.2514, 87.8126, 78.7338, 78.7338)
            Conventional setting: I 1 2/m 1 (No. 12)
                       Unit cell: (6.35, 22.54, 22.54, 90, 90, 90)
                 Change of basis: -z,y,x
                         Inverse: z,y,-x
      Maximal angular difference: 0.000 degrees

Symmetry in minimum-lengths cell: C 1 2/m 1 (x-y+z,x+y,-x+y) (No. 12)
      Input minimum-lengths cell: (6.35, 16.2514, 16.2514, 87.8126, 78.7338, 78.7338)
           Symmetry-adapted cell: (6.35, 16.2514, 16.2514, 87.8126, 78.7338, 78.7338)
            Conventional setting: I 1 2/m 1 (No. 12)
                       Unit cell: (6.35, 22.54, 22.54, 90, 90, 90)
                 Change of basis: -z,x,-y
                         Inverse: y,-z,-x
      Maximal angular difference: 0.000 degrees

Symmetry in minimum-lengths cell: P -1 (No. 2)
      Input minimum-lengths cell: (6.35, 16.2514, 16.2514, 87.8126, 78.7338, 78.7338)
           Symmetry-adapted cell: (6.35, 16.2514, 16.2514, 87.8126, 78.7338, 78.7338)
            Conventional setting: P -1 (No. 2)
                       Unit cell: (6.35, 16.2514, 16.2514, 87.8126, 78.7338, 78.7338)
                 Change of basis: x+z,-x+y,-x-y
                         Inverse: -1/2*y-1/2*z,1/2*y-1/2*z,x+1/2*y+1/2*z
      Maximal angular difference: 0.000 degrees

""")
  #
  run_and_compare(
    ['iotbx.lattice_symmetry'
       ' --unit_cell="78.9 82.3 57.0 90 93.4 90" C --delta=1.4'],
    """\

Input
=====

Unit cell: (78.9, 82.3, 57, 90, 93.4, 90)
Space group: P 1 (a+b,a-b,-c) (No. 1)

Angular tolerance: 1.400 degrees

Similar symmetries
==================

Symmetry in minimum-lengths cell: R -3 m :R (No. 166)
      Input minimum-lengths cell: (57, 57.0055, 57.0055, 92.4166, 92.3522, 92.3522)
           Symmetry-adapted cell: (57.0037, 57.0037, 57.0037, 92.3737, 92.3737, 92.3737)
            Conventional setting: R -3 m :H (No. 166)
                       Unit cell: (82.2678, 82.2678, 94.5557, 90, 90, 120)
                 Change of basis: -1/3*x-y+1/3*z,-2/3*x+2/3*z,-2/3*x-1/3*z
                         Inverse: -1/2*y-z,-x+1/2*y,y-z
      Maximal angular difference: 0.045 degrees

Symmetry in minimum-lengths cell: C 1 2/m 1 (x-y,x+y,z) (No. 12)
      Input minimum-lengths cell: (57, 57.0055, 57.0055, 92.4166, 92.3522, 92.3522)
           Symmetry-adapted cell: (57.0027, 57.0027, 57.0055, 92.3844, 92.3844, 92.3522)
            Conventional setting: C 1 2/m 1 (No. 12)
                       Unit cell: (78.9424, 82.2517, 57.0055, 90, 93.4445, 90)
                 Change of basis: -1/2*x-1/2*y-1/2*z,-1/2*x-1/2*y+1/2*z,-x+y
                         Inverse: -1/2*x-1/2*y-1/2*z,-1/2*x-1/2*y+1/2*z,-x+y
      Maximal angular difference: 0.045 degrees

Symmetry in minimum-lengths cell: C 1 2/m 1 (x+y,z,x-y) (No. 12)
      Input minimum-lengths cell: (57, 57.0055, 57.0055, 92.4166, 92.3522, 92.3522)
           Symmetry-adapted cell: (57.0027, 57.0055, 57.0027, 92.3844, 92.3522, 92.3844)
            Conventional setting: C 1 2/m 1 (No. 12)
                       Unit cell: (78.9424, 82.2517, 57.0055, 90, 93.4445, 90)
                 Change of basis: -1/2*x+1/2*y-1/2*z,1/2*x-1/2*y-1/2*z,-x-y
                         Inverse: -1/2*x+1/2*y-1/2*z,1/2*x-1/2*y-1/2*z,-x-y
      Maximal angular difference: 0.045 degrees

Symmetry in minimum-lengths cell: C 1 2/m 1 (z,x-y,x+y) (No. 12)
      Input minimum-lengths cell: (57, 57.0055, 57.0055, 92.4166, 92.3522, 92.3522)
           Symmetry-adapted cell: (57, 57.0055, 57.0055, 92.4166, 92.3522, 92.3522)
            Conventional setting: C 1 2/m 1 (No. 12)
                       Unit cell: (78.9, 82.3, 57, 90, 93.4, 90)
                 Change of basis: x,y,z
                         Inverse: x,y,z
      Maximal angular difference: 0.000 degrees

Symmetry in minimum-lengths cell: P -1 (No. 2)
      Input minimum-lengths cell: (57, 57.0055, 57.0055, 92.4166, 92.3522, 92.3522)
           Symmetry-adapted cell: (57, 57.0055, 57.0055, 92.4166, 92.3522, 92.3522)
            Conventional setting: P -1 (No. 2)
                       Unit cell: (57, 57.0055, 57.0055, 92.4166, 92.3522, 92.3522)
                 Change of basis: -z,-x-y,-x+y
                         Inverse: -1/2*y-1/2*z,-1/2*y+1/2*z,-x
      Maximal angular difference: 0.000 degrees

""")
  print "OK"

if (__name__ == "__main__"):
  exercise(sys.argv[1:])
