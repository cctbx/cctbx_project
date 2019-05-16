from __future__ import absolute_import, division, print_function
from cctbx import xray
from cctbx import crystal
import cctbx.crystal.coordination_sequences
from cctbx.array_family import flex
from six.moves import zip

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
  icsd_pairs = icsd_structure.show_distances(
    distance_cutoff=2.5, keep_pair_asu_table=True)
  print()
  primitive_structure = icsd_structure.primitive_setting()
  primitive_structure.show_summary().show_scatterers()
  print()
  p1_structure = primitive_structure.expand_to_p1()
  p1_structure.show_summary().show_scatterers()
  print()
  p1_pairs = p1_structure.show_distances(
    distance_cutoff=2.5, keep_pair_asu_table=True)
  print()
  for label,structure,pairs in [("ICSD", icsd_structure,icsd_pairs),
                                ("P1", p1_structure,p1_pairs)]:
    print("Coordination sequences for", label, "structure")
    term_table = crystal.coordination_sequences.simple(
      pair_asu_table=pairs.pair_asu_table,
      max_shell=10)
    crystal.coordination_sequences.show_terms(
      structure=structure,
      term_table=term_table)
    print()
  icsd_f_calc = icsd_structure.structure_factors(
    d_min=1, algorithm="direct").f_calc()
  icsd_f_calc_in_p1 = icsd_f_calc.primitive_setting().expand_to_p1()
  p1_f_calc = icsd_f_calc_in_p1.structure_factors_from_scatterers(
    xray_structure=p1_structure, algorithm="direct").f_calc()
  for h,i,p in zip(icsd_f_calc_in_p1.indices(),
                   icsd_f_calc_in_p1.data(),
                   p1_f_calc.data()):
    print(h, abs(i), abs(p)*3)
  print("OK")

if (__name__ == "__main__"):
  demo()
