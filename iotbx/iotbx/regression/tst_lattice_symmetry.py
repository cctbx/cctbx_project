from libtbx.test_utils import show_diff
from libtbx import easy_run
import sys

def run_and_compare(commands, expected_output):
  def digest(lines):
    result = []
    for line in lines:
      if (line.strip().startswith("Change of basis: ")): continue
      if (line.strip().startswith("Inverse: ")): continue
      result.append(line)
    result.sort()
    return result
  for command in commands:
    lines = easy_run.fully_buffered(
      command=command).raise_if_errors().stdout_lines
    assert not show_diff(digest(lines), digest(expected_output.splitlines()))

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
    ['iotbx.lattice_symmetry --unit_cell=12,12.2,12.1,89,90,92 F',
     'cctbx.lattice_symmetry'],
    """\

Input
=====

Unit cell: (12, 12.2, 12.1, 89, 90, 92)
Space group: P 1 (-a+b+c,a-b+c,a+b-c) (No. 1)

Angular tolerance: 3.000 degrees

Similar symmetries
==================

Symmetry in minimum-lengths cell: F m -3 m (-x+y+z,x-y+z,x+y-z) (No. 225)
      Input minimum-lengths cell: (8.40567, 8.51612, 8.52071, 61.4489, 61.5879, 60.4641)
           Symmetry-adapted cell: (8.55619, 8.55619, 8.55619, 60, 60, 60)
            Conventional setting: F m -3 m (No. 225)
                       Unit cell: (12.1003, 12.1003, 12.1003, 90, 90, 90)
                 Change of basis: z,-x,-y
                         Inverse: -y,-z,x
      Maximal angular difference: 2.236 degrees

Symmetry in minimum-lengths cell: I 4/m m m (y-z,-x+z,x+z) (No. 139)
      Input minimum-lengths cell: (8.40567, 8.51612, 8.52071, 61.4489, 61.5879, 60.4641)
           Symmetry-adapted cell: (8.55628, 8.55614, 8.55614, 60.0011, 59.9994, 59.9994)
            Conventional setting: I 4/m m m (No. 139)
                       Unit cell: (8.55628, 8.55628, 12.1, 90, 90, 90)
                 Change of basis: -x+y,-x-y,z
                         Inverse: -1/2*x-1/2*y,1/2*x-1/2*y,z
      Maximal angular difference: 2.236 degrees

Symmetry in minimum-lengths cell: I 4/m m m (-y+z,x-z,y+z) (No. 139)
      Input minimum-lengths cell: (8.40567, 8.51612, 8.52071, 61.4489, 61.5879, 60.4641)
           Symmetry-adapted cell: (8.53852, 8.59142, 8.53852, 59.7948, 60.4103, 59.7948)
            Conventional setting: I 4/m m m (No. 139)
                       Unit cell: (8.59142, 8.59142, 12, 90, 90, 90)
                 Change of basis: -y+z,y+z,-x
                         Inverse: -z,-1/2*x+1/2*y,1/2*x+1/2*y
      Maximal angular difference: 2.236 degrees

Symmetry in minimum-lengths cell: I 4/m m m (y+z,-y+z,x-z) (No. 139)
      Input minimum-lengths cell: (8.40567, 8.51612, 8.52071, 61.4489, 61.5879, 60.4641)
           Symmetry-adapted cell: (8.57387, 8.57387, 8.52071, 60.2049, 60.2049, 59.5902)
            Conventional setting: I 4/m m m (No. 139)
                       Unit cell: (8.52071, 8.52071, 12.2, 90, 90, 90)
                 Change of basis: -x+z,-x-z,-y
                         Inverse: -1/2*x-1/2*y,-z,1/2*x-1/2*y
      Maximal angular difference: 2.236 degrees

Symmetry in minimum-lengths cell: R -3 m :R (No. 166)
      Input minimum-lengths cell: (8.40567, 8.51612, 8.52071, 61.4489, 61.5879, 60.4641)
           Symmetry-adapted cell: (8.481, 8.481, 8.481, 61.1714, 61.1714, 61.1714)
            Conventional setting: R -3 m :H (No. 166)
                       Unit cell: (8.63072, 8.63072, 20.5883, 90, 90, 120)
                 Change of basis: 4/3*x-2/3*y+2/3*z,2/3*x+2/3*y+4/3*z,-1/3*x-1/3*y+1/3*z
                         Inverse: 1/2*x-z,-1/2*x+1/2*y-z,1/2*y+z
      Maximal angular difference: 1.474 degrees

Symmetry in minimum-lengths cell: R -3 m :H (x+z,-y+z,-3*z) (No. 166)
      Input minimum-lengths cell: (8.40567, 8.51612, 8.52071, 61.4489, 61.5879, 60.4641)
           Symmetry-adapted cell: (8.481, 8.481, 8.63072, 60.5722, 60.5722, 60)
            Conventional setting: R -3 m :H (No. 166)
                       Unit cell: (8.481, 8.481, 21.3218, 90, 90, 120)
                 Change of basis: -4/3*x-2/3*y-2/3*z,-2/3*x+2/3*y-4/3*z,1/3*x-1/3*y-1/3*z
                         Inverse: -1/2*x+z,-1/2*x+1/2*y-z,-1/2*y-z
      Maximal angular difference: 1.498 degrees

Symmetry in minimum-lengths cell: R -3 m :H (-y+z,-3*z,x+z) (No. 166)
      Input minimum-lengths cell: (8.40567, 8.51612, 8.52071, 61.4489, 61.5879, 60.4641)
           Symmetry-adapted cell: (8.53148, 8.58082, 8.53148, 60.19, 60, 60.19)
            Conventional setting: R -3 m :H (No. 166)
                       Unit cell: (8.53148, 8.53148, 21.0788, 90, 90, 120)
                 Change of basis: -2/3*x+2/3*y+4/3*z,2/3*x+4/3*y+2/3*z,-1/3*x+1/3*y-1/3*z
                         Inverse: -1/2*x+1/2*y-z,1/2*y+z,1/2*x-z
      Maximal angular difference: 2.177 degrees

Symmetry in minimum-lengths cell: R -3 m :H (-3*z,x+z,-y+z) (No. 166)
      Input minimum-lengths cell: (8.40567, 8.51612, 8.52071, 61.4489, 61.5879, 60.4641)
           Symmetry-adapted cell: (8.53148, 8.58082, 8.58082, 60, 59.8085, 59.8085)
            Conventional setting: R -3 m :H (No. 166)
                       Unit cell: (8.58082, 8.58082, 20.8371, 90, 90, 120)
                 Change of basis: 2/3*x-4/3*y+2/3*z,4/3*x-2/3*y-2/3*z,1/3*x+1/3*y+1/3*z
                         Inverse: 1/2*y+z,-1/2*x+z,1/2*x-1/2*y+z
      Maximal angular difference: 2.177 degrees

Symmetry in minimum-lengths cell: I m m m (y-z,-x+z,x+z) (No. 71)
      Input minimum-lengths cell: (8.40567, 8.51612, 8.52071, 61.4489, 61.5879, 60.4641)
           Symmetry-adapted cell: (8.40567, 8.55614, 8.55614, 61.1489, 60.58, 60.58)
            Conventional setting: I m m m (No. 71)
                       Unit cell: (8.40567, 8.70429, 12.1, 90, 90, 90)
                 Change of basis: x+y,-x+y,z
                         Inverse: 1/2*x-1/2*y,1/2*x+1/2*y,z
      Maximal angular difference: 1.187 degrees

Symmetry in minimum-lengths cell: I m m m (-y+z,x-z,y+z) (No. 71)
      Input minimum-lengths cell: (8.40567, 8.51612, 8.52071, 61.4489, 61.5879, 60.4641)
           Symmetry-adapted cell: (8.53852, 8.51612, 8.53852, 60.0867, 60.9908, 60.0867)
            Conventional setting: I m m m (No. 71)
                       Unit cell: (8.51612, 8.66606, 12, 90, 90, 90)
                 Change of basis: -y+z,y+z,-x
                         Inverse: -z,-1/2*x+1/2*y,1/2*x+1/2*y
      Maximal angular difference: 2.000 degrees

Symmetry in minimum-lengths cell: I m m m (y+z,-y+z,x-z) (No. 71)
      Input minimum-lengths cell: (8.40567, 8.51612, 8.52071, 61.4489, 61.5879, 60.4641)
           Symmetry-adapted cell: (8.57387, 8.57387, 8.52071, 60.2049, 60.2049, 59.5902)
            Conventional setting: I m m m (No. 71)
                       Unit cell: (8.52071, 8.52071, 12.2, 90, 90, 90)
                 Change of basis: -x+z,-x-z,-y
                         Inverse: -1/2*x-1/2*y,-z,1/2*x-1/2*y
      Maximal angular difference: 2.236 degrees

Symmetry in minimum-lengths cell: F m m m (-x+y+z,x-y+z,x+y-z) (No. 69)
      Input minimum-lengths cell: (8.40567, 8.51612, 8.52071, 61.4489, 61.5879, 60.4641)
           Symmetry-adapted cell: (8.55628, 8.59142, 8.52071, 60, 60.4101, 59.5899)
            Conventional setting: F m m m (No. 69)
                       Unit cell: (12, 12.1, 12.2, 90, 90, 90)
                 Change of basis: x,z,-y
                         Inverse: x,-z,y
      Maximal angular difference: 2.236 degrees

Symmetry in minimum-lengths cell: C 1 2/m 1 (x+y,-x+y,x-y+z) (No. 12)
      Input minimum-lengths cell: (8.40567, 8.51612, 8.52071, 61.4489, 61.5879, 60.4641)
           Symmetry-adapted cell: (8.55628, 8.59142, 8.52071, 60, 60.4101, 59.5899)
            Conventional setting: I 1 2/m 1 (No. 12)
                       Unit cell: (8.52071, 12.2, 8.52071, 90, 90.4755, 90)
                 Change of basis: x-z,y,x+z
                         Inverse: 1/2*x+1/2*z,y,-1/2*x+1/2*z
      Maximal angular difference: 2.236 degrees

Symmetry in minimum-lengths cell: C 1 2/m 1 (x-y+z,x+y,-x+y) (No. 12)
      Input minimum-lengths cell: (8.40567, 8.51612, 8.52071, 61.4489, 61.5879, 60.4641)
           Symmetry-adapted cell: (8.40567, 8.59142, 8.52071, 61.1478, 61.0005, 60.1608)
            Conventional setting: I 1 2/m 1 (No. 12)
                       Unit cell: (8.40567, 12.1, 8.70429, 90, 90.9476, 90)
                 Change of basis: x+y,z,x-y
                         Inverse: 1/2*x+1/2*z,1/2*x-1/2*z,y
      Maximal angular difference: 1.001 degrees

Symmetry in minimum-lengths cell: C 1 2/m 1 (-x+y,x-y+z,x+y) (No. 12)
      Input minimum-lengths cell: (8.40567, 8.51612, 8.52071, 61.4489, 61.5879, 60.4641)
           Symmetry-adapted cell: (8.55628, 8.51612, 8.52071, 60.2943, 60.9905, 59.8794)
            Conventional setting: I 1 2/m 1 (No. 12)
                       Unit cell: (8.51612, 12, 8.66606, 90, 90.4716, 90)
                 Change of basis: y-z,x,-y-z
                         Inverse: y,1/2*x-1/2*z,-1/2*x-1/2*z
      Maximal angular difference: 2.000 degrees

Symmetry in minimum-lengths cell: C 1 2/m 1 (z,x-y,x+y) (No. 12)
      Input minimum-lengths cell: (8.40567, 8.51612, 8.52071, 61.4489, 61.5879, 60.4641)
           Symmetry-adapted cell: (8.40567, 8.51842, 8.51842, 61.4489, 61.0277, 61.0277)
            Conventional setting: I 1 2/m 1 (No. 12)
                       Unit cell: (8.40567, 8.70429, 12.1, 90, 90.7257, 90)
                 Change of basis: x+y,x-y,-z
                         Inverse: 1/2*x+1/2*y,1/2*x-1/2*y,-z
      Maximal angular difference: 1.172 degrees

Symmetry in minimum-lengths cell: C 1 2/m 1 (x+y,z,x-y) (No. 12)
      Input minimum-lengths cell: (8.40567, 8.51612, 8.52071, 61.4489, 61.5879, 60.4641)
           Symmetry-adapted cell: (8.46339, 8.51612, 8.46339, 60.9617, 61.5908, 60.9617)
            Conventional setting: I 1 2/m 1 (No. 12)
                       Unit cell: (8.51612, 8.66606, 12, 90, 91.4324, 90)
                 Change of basis: y-z,y+z,x
                         Inverse: z,1/2*x+1/2*y,-1/2*x+1/2*y
      Maximal angular difference: 1.474 degrees

Symmetry in minimum-lengths cell: C 1 2/m 1 (x-y,x+y,z) (No. 12)
      Input minimum-lengths cell: (8.40567, 8.51612, 8.52071, 61.4489, 61.5879, 60.4641)
           Symmetry-adapted cell: (8.46108, 8.46108, 8.52071, 61.5187, 61.5187, 60.4668)
            Conventional setting: I 1 2/m 1 (No. 12)
                       Unit cell: (8.52071, 8.52071, 12.2, 90, 92.1185, 90)
                 Change of basis: x-z,-x-z,y
                         Inverse: 1/2*x-1/2*y,z,-1/2*x-1/2*y
      Maximal angular difference: 0.860 degrees

Symmetry in minimum-lengths cell: C 1 2/m 1 (-z,-2*x+z,x+y) (No. 12)
      Input minimum-lengths cell: (8.40567, 8.51612, 8.52071, 61.4489, 61.5879, 60.4641)
           Symmetry-adapted cell: (8.53686, 8.61072, 8.52071, 60.3452, 60.0626, 59.589)
            Conventional setting: I 1 2/m 1 (No. 12)
                       Unit cell: (8.52071, 8.52071, 12.2, 90, 90.6981, 90)
                 Change of basis: x+z,x-z,y
                         Inverse: 1/2*x+1/2*y,z,1/2*x-1/2*y
      Maximal angular difference: 2.177 degrees

Symmetry in minimum-lengths cell: C 1 2/m 1 (-z,-x+y,2*x+z) (No. 12)
      Input minimum-lengths cell: (8.40567, 8.51612, 8.52071, 61.4489, 61.5879, 60.4641)
           Symmetry-adapted cell: (8.46339, 8.51612, 8.61299, 60.3713, 60.9859, 59.7937)
            Conventional setting: I 1 2/m 1 (No. 12)
                       Unit cell: (8.66606, 8.51612, 12, 90, 91.4076, 90)
                 Change of basis: -y-z,-y+z,-x
                         Inverse: -z,-1/2*x-1/2*y,-1/2*x+1/2*y
      Maximal angular difference: 1.498 degrees

Symmetry in minimum-lengths cell: C 1 2/m 1 (-x+y,z,2*x-z) (No. 12)
      Input minimum-lengths cell: (8.40567, 8.51612, 8.52071, 61.4489, 61.5879, 60.4641)
           Symmetry-adapted cell: (8.40567, 8.51842, 8.59369, 61.1477, 60.7211, 60.4369)
            Conventional setting: I 1 2/m 1 (No. 12)
                       Unit cell: (8.70429, 8.40567, 12.1, 90, 90.7008, 90)
                 Change of basis: x-y,x+y,z
                         Inverse: 1/2*x+1/2*y,-1/2*x+1/2*y,z
      Maximal angular difference: 1.187 degrees

Symmetry in minimum-lengths cell: P -1 (No. 2)
      Input minimum-lengths cell: (8.40567, 8.51612, 8.52071, 61.4489, 61.5879, 60.4641)
           Symmetry-adapted cell: (8.40567, 8.51612, 8.52071, 61.4489, 61.5879, 60.4641)
            Conventional setting: P -1 (No. 2)
                       Unit cell: (8.40567, 8.51612, 8.52071, 61.4489, 61.5879, 60.4641)
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
