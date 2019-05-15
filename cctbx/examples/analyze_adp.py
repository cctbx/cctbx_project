from __future__ import absolute_import, division, print_function
# Simple example for the use of the adptbx.

from cctbx import crystal
from cctbx import adptbx # anisotropic displacement parameter toolbox
from six.moves import range

def run():
  symmetry = crystal.symmetry(
    unit_cell=(10.67, 10.67, 4.68, 90, 90, 120),
    space_group_symbol="P 3")

  special_position_settings = crystal.special_position_settings(
    crystal_symmetry=symmetry,
    min_distance_sym_equiv=0.5)

  site = (0, 0, 0.236)
  u_cif = ((0.17, 0.17, 0.19, 0.09, 0, 0))

  site_symmetry = special_position_settings.site_symmetry(site)

  print("Input Ucif:", u_cif)
  u_star = adptbx.u_cif_as_u_star(symmetry.unit_cell(), u_cif)
  if (not site_symmetry.is_compatible_u_star(u_star)):
    print("Warning: ADP tensor is incompatible with site symmetry.")
  u_star = site_symmetry.average_u_star(u_star)
  u_cif = adptbx.u_star_as_u_cif(symmetry.unit_cell(), u_star)
  print("Averaged Ucif:", u_cif)

  u_cart = adptbx.u_star_as_u_cart(symmetry.unit_cell(), u_star)
  eigenvalues = adptbx.eigenvalues(u_cart)
  if (not adptbx.is_positive_definite(eigenvalues)):
    print("ADP tensor is not positive definite.")

  print("Eigenvectors and values:")
  eigensystem = adptbx.eigensystem(u_cart)
  for i in range(3):
    print("  v=(%.5f %.5f %.5f) " % eigensystem.vectors(i), end=' ')
    print("lambda=%.4f" % (eigensystem.values()[i],))

if (__name__ == "__main__"):
  run()
