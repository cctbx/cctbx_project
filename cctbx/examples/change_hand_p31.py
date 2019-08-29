"""
This example shows that a change of hand ("flipping coordinates") involves
changing the space group if the space group is enantimorphic.
Note that the interatomic distances do not change if the space group
symmetry is transformed correctly, but do change if the original space
group symmetry is simply retained.
"""
from __future__ import absolute_import, division, print_function

from cctbx import xray
from cctbx import crystal
from cctbx.array_family import flex

def run():
  print(__doc__)
  crystal_symmetry = crystal.symmetry(
    unit_cell=(5, 5, 6, 90, 90, 120),
    space_group_symbol="P 31")
  distance_cutoff = 2.5
  scatterers = flex.xray_scatterer()
  for i,site in enumerate([(0.7624, 0.5887, 0.3937),
                           (0.2813, 0.9896, 0.9449),
                           (0.4853, 0.8980, 0.4707)]):
    scatterers.append(xray.scatterer(
      label="Se%d"%(i+1), site=site))
  given_structure = xray.structure(
    crystal_symmetry=crystal_symmetry,
    scatterers=scatterers)
  print("==================")
  print("Original structure")
  print("==================")
  given_structure.show_summary().show_scatterers()
  print("Interatomic distances:")
  given_structure.show_distances(distance_cutoff=distance_cutoff)
  print()
  print()
  print("=====================================================")
  print("Other hand with sites flipped and space group changed")
  print("=====================================================")
  other_hand = given_structure.change_hand()
  other_hand.show_summary().show_scatterers()
  print("Interatomic distances:")
  other_hand.show_distances(distance_cutoff=distance_cutoff)
  print()
  print("==================================")
  print("Other hand with sites flipped only")
  print("==================================")
  other_sites_orig_symmetry = xray.structure(
    crystal_symmetry=crystal_symmetry,
    scatterers=other_hand.scatterers())
  other_sites_orig_symmetry.show_summary().show_scatterers()
  print("Interatomic distances:")
  other_sites_orig_symmetry.show_distances(distance_cutoff=distance_cutoff)
  print()

if (__name__ == "__main__"):
  run()
