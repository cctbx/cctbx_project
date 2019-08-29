from __future__ import absolute_import, division, print_function

from scitbx.array_family import flex, shared
import gltbx.viewer_utils
from six.moves import range

# this is the pdb file used to create the bonds table for testing
ethanol_pdb = """
COMPND      METHANOL
AUTHOR      DAVE WOODCOCK  96 01 03
CRYST1  100.000  100.000  100.000  90.00  90.00  90.00  P 1
SCALE1      0.010000  0.000000  0.000000       0.000000
SCALE2      0.000000  0.010000  0.000000       0.000000
SCALE3      0.000000  0.000000  0.010000       0.000000
ATOM      1  C1  EOH A   1      -0.426  -0.115  -0.147  1.00  0.00
ATOM      2  O   EOH A   1      -0.599   1.244  -0.481  1.00  0.00
ATOM      3  H11 EOH A   1      -0.750  -0.738  -0.981  1.00  0.00
ATOM      4  H12 EOH A   1      -1.022  -0.351   0.735  1.00  0.00
ATOM      5  HO  EOH A   1      -1.642   1.434  -0.689  1.00  0.00
ATOM      6  C2  EOH A   1       1.047  -0.383   0.147  1.00  0.00
ATOM      7  H21 EOH A   1       1.370   0.240   0.981  1.00  0.00
ATOM      8  H22 EOH A   1       1.642  -0.147  -0.735  1.00  0.00
ATOM      9  H23 EOH A   1       1.180  -1.434   0.405  1.00  0.00
HETATM   10  O   HOH W   2      29.478  23.354  61.364  1.00  8.67
END
"""

#---
def exercise():
  bonds = shared.stl_set_unsigned([
    [1, 2, 3, 5],
    [0, 4],
    [0],
    [0],
    [1],
    [0, 6, 7, 8],
    [5],
    [5],
    [5],
    []
  ])

  atoms_drawable = flex.bool([True for i in range(0, 10) ])
  atoms_drawable_non_h = flex.bool(
    [True, True, False, False, False, True, False, False, False, True]
  )

  # all atoms drawable
  visibility = gltbx.viewer_utils.atom_visibility(
    bonds            = bonds,
    atoms_drawable   = atoms_drawable,
    flag_show_points = True
  )
  assert (list(visibility.bonds_visible) ==
    [True, True, True, True, True, True, True, True, True, False]
  )
  assert (list(visibility.points_visible) ==
    [False, False, False, False, False, False, False, False, False, True]
  )

  # hydrogens off, points on
  visibility = gltbx.viewer_utils.atom_visibility(
    bonds            = bonds,
    atoms_drawable   = atoms_drawable_non_h,
    flag_show_points = True
  )
  assert (list(visibility.bonds_visible) ==
    [True, True, False, False, False, True, False, False, False, False]
  )
  assert (list(visibility.points_visible) ==
    [False, False, False, False, False, False, False, False, False, True]
  )

  # equivalent to "not element H"
  atoms_selected_non_h = flex.bool(
    [True, True, False, False, False, True, False, False, False, True]
  )

  visibility.get_selection_visibility(
    bonds          = bonds,
    atoms_selected = atoms_selected_non_h
  )
  assert (list(visibility.selected_points_visible) ==
    [False, False, False, False, False, False, False, False, False, True]
  )

  # no hydrogens, no points
  visibility = gltbx.viewer_utils.atom_visibility(
    bonds            = bonds,
    atoms_drawable   = atoms_drawable_non_h,
    flag_show_points = False
  )
  assert (list(visibility.bonds_visible) ==
    [True, True, False, False, False, True, False, False, False, False]
  )
  assert (not True in list(visibility.points_visible))

  visibility.get_selection_visibility(
    bonds          = bonds,
    atoms_selected = atoms_selected_non_h
  )
  assert (not True in list(visibility.selected_points_visible))

if (__name__ == "__main__"):
  exercise()
