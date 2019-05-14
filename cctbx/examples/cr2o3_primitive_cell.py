from __future__ import absolute_import, division, print_function
from cctbx import xray
from cctbx import crystal
from cctbx.array_family import flex

def demo():
  """
  Result of ICSD query:
    N * -Cr2O3-[R3-CH] Baster, M.;Bouree, F.;Kowalska, A.;Latacz, Z(2000)
    C 4.961950 4.961950 13.597400 90.000000 90.000000 120.000000
    S GRUP R -3 C
    A Cr1    0.000000 0.000000 0.347570 0.000000
    A O1    0.305830 0.000000 0.250000
  """
  crystal_symmetry = crystal.symmetry(
    unit_cell="4.961950 4.961950 13.597400 90.000000 90.000000 120.000000",
    space_group_symbol="R -3 C")
  scatterers = flex.xray_scatterer()
  scatterers.append(xray.scatterer(
    label="Cr1", site=(0.000000,0.000000,0.347570)))
  scatterers.append(xray.scatterer(
    label="O1", site=(0.305830,0.000000,0.250000)))
  icsd_structure = xray.structure(
    crystal_symmetry=crystal_symmetry,
    scatterers=scatterers)
  icsd_structure.show_summary().show_scatterers()
  print()
  icsd_structure.show_distances(distance_cutoff=2.5)
  print()
  primitive_structure = icsd_structure.primitive_setting()
  primitive_structure.show_summary().show_scatterers()
  print()
  p1_structure = primitive_structure.expand_to_p1()
  p1_structure.show_summary().show_scatterers()
  print()
  print("OK")

if (__name__ == "__main__"):
  demo()
